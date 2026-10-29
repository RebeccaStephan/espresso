#
# Copyright (C) 2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Solvation coupling of particles near the box edge and near MPI subdomain
boundaries in two-component color gradient LB.

A particle within half a cell of the box edge or of a subdomain boundary
has several periodic images in the LB halo, and is coupled on two ranks
when the box is split. Tests:
  - Total momentum (particle + fluid) is conserved for particles at the
    droplet interface, inside the box, at the box edge, at a box corner and
    at the subdomain boundary (x = 6 when the box is split in two along x),
    with solvation and with friction coupling only
  - The force on a particle does not change when the particle and the
    droplet are shifted together by whole cells from the interior to the
    box edge or to the subdomain boundary

The tests run on the CPU and on the GPU, in double precision, on 1 and on
2 MPI ranks.
"""

import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb


DOMAIN_SIZE = 12
AGRID = 1.0
TAU = 1.0
RHO_0 = 1.0
EPSILON = 1e-6
VISCOSITY = 1.0 / 6.0
RADIUS = 4.0
SMOOTHING_WIDTH = 2.0
GAMMA = 0.5
SIGMA = 0.005
BETA = 0.8
DELTA_MU = 2.0
# the droplet center lies at this offset from the particle, so that the
# interface passes through the particle
CENTER_OFFSET = np.array([-RADIUS, -0.1, 0.1])

# particle positions and the interior positions they are shifted from by
# whole cells
POSITIONS = {
    "interior": ([8.3, 6.1, 5.9], [8.3, 6.1, 5.9]),
    "box edge lower": ([0.1, 6.1, 5.9], [8.1, 6.1, 5.9]),
    "box edge upper": ([11.9, 6.1, 5.9], [7.9, 6.1, 5.9]),
    "box corner": ([0.1, 0.1, 5.9], [8.1, 6.1, 5.9]),
    "subdomain boundary lower": ([5.9, 6.1, 5.9], [8.9, 6.1, 5.9]),
    "subdomain boundary upper": ([6.1, 6.1, 5.9], [9.1, 6.1, 5.9]),
}


def droplet_densities(center):
    """Densities of a droplet of component b in a matrix of component a,
    centered at ``center`` (in MD units) with periodic distances."""
    grid_size = int(DOMAIN_SIZE / AGRID)
    x = (np.arange(grid_size) + 0.5) * AGRID
    grid = np.meshgrid(x, x, x, indexing="ij")
    distances = [g - c for g, c in zip(grid, center)]
    distances = [d - DOMAIN_SIZE * np.round(d / DOMAIN_SIZE)
                 for d in distances]
    distance = np.sqrt(sum(d**2 for d in distances)) / AGRID
    rho_a = (RHO_0 - EPSILON * RHO_0) / 2.0 * (
        np.tanh((distance - RADIUS / AGRID) / (SMOOTHING_WIDTH / AGRID)
                * 2.64665) - 1.0) + RHO_0
    rho_b = RHO_0 + EPSILON * RHO_0 - rho_a
    return np.stack([rho_a, rho_b], axis=-1)


class ColorGradientSolvationMomentumTest:
    """Solvation coupling across the box edge and subdomain boundaries."""

    system = espressomd.System(box_l=[DOMAIN_SIZE] * 3)
    system.time_step = TAU
    system.cell_system.skin = 0.4
    system.periodicity = [True, True, True]

    lb_params = {}

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.lb = None

    def setup_particle(self, pos, delta_mu):
        """Droplet with its interface through a resting particle."""
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=TAU,
            kinematic_viscosity=[VISCOSITY, VISCOSITY], sigma=SIGMA,
            beta=BETA, **self.lb_params)
        self.system.lb = lbf
        lbf[:, :, :].density = droplet_densities(
            np.array(pos) + CENTER_OFFSET)
        lbf.init_two_component()
        self.system.integrator.run(1)
        p = self.system.part.add(pos=pos, v=[0., 0., 0.],
                                 solvation_delta_mu=delta_mu)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)
        return lbf, p

    @staticmethod
    def total_momentum(lbf, p):
        density = np.copy(lbf[:, :, :].density).sum(axis=-1)
        velocity = np.copy(lbf[:, :, :].velocity)
        fluid = np.sum(density[..., np.newaxis] * velocity,
                       axis=(0, 1, 2)) * AGRID**3
        return fluid + p.mass * np.copy(p.v)

    def check_momentum_conservation(self, delta_mu):
        for label, (pos, _) in POSITIONS.items():
            with self.subTest(position=label):
                lbf, p = self.setup_particle(pos, delta_mu)
                momentum_0 = self.total_momentum(lbf, p)
                for n_steps in (1, 9):
                    self.system.integrator.run(n_steps)
                    np.testing.assert_allclose(
                        self.total_momentum(lbf, p), momentum_0,
                        rtol=0., atol=1e-10)
                self.tearDown()

    def test_momentum_conservation_solvation(self):
        self.check_momentum_conservation(DELTA_MU)

    def test_momentum_conservation_friction(self):
        self.check_momentum_conservation(0.)

    def test_translation_invariance(self):
        def force(pos):
            _, p = self.setup_particle(pos, DELTA_MU)
            self.system.integrator.run(2)
            f = np.copy(p.f)
            self.tearDown()
            return f

        for label, (pos, interior_pos) in POSITIONS.items():
            with self.subTest(position=label):
                f = force(pos)
                f_ref = force(interior_pos)
                self.assertGreater(np.linalg.norm(f_ref), 0.1)
                np.testing.assert_allclose(f, f_ref, rtol=0., atol=1e-10)


@utx.skipIfMissingFeatures(["WALBERLA"])
class ColorGradientSolvationMomentumTestCPU(
        ColorGradientSolvationMomentumTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class ColorGradientSolvationMomentumTestGPU(
        ColorGradientSolvationMomentumTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}


if __name__ == "__main__":
    ut.main()
