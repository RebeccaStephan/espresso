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
Velocity handed to the MD coupling by the two-component color gradient LB,
compared to the single-component LB.

A color gradient fluid whose second component is negligible, without
surface tension, evolves like a single-component fluid. Tests:
  - Under a constant external force, the fluid velocity after an LB step is
    the post-collision velocity of the single-component LB, i.e. it
    contains the force of that step once and a half
  - A particle coupled to the fluid by friction follows the same trajectory
    as with the single-component LB

The tests run on the CPU and on the GPU, in double and in single precision.
"""

import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb


GRID_SIZE = 8
AGRID = 1.0
TAU = 1.0
RHO_0 = 1.0
EPSILON = 1e-6
VISCOSITY = 1.0 / 6.0
GAMMA = 1.0
EXT_FORCE_DENSITY = np.array([1e-3, -2e-3, 5e-4])
N_STEPS = 40


class ColorGradientSingleComponentParityTest:
    """Color gradient LB with a negligible second component against the
    single-component LB."""

    system = espressomd.System(box_l=[GRID_SIZE * AGRID] * 3)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    lb_params = {}

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.lb = None
        self.system.time = 0.

    def tol(self, double, single):
        return single if self.lb_params["single_precision"] else double

    def single_component_fluid(self, **kwargs):
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, tau=TAU, density=RHO_0,
            kinematic_viscosity=VISCOSITY, **kwargs, **self.lb_params)
        self.system.lb = lbf
        return lbf

    def color_gradient_fluid(self, **kwargs):
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, tau=TAU, density=RHO_0,
            kinematic_viscosity=[VISCOSITY, VISCOSITY], sigma=0., beta=0.8,
            **kwargs, **self.lb_params)
        self.system.lb = lbf
        densities = np.zeros(3 * (GRID_SIZE,) + (2,))
        densities[..., 0] = RHO_0 * (1. - EPSILON)
        densities[..., 1] = RHO_0 * EPSILON
        lbf[:, :, :].density = densities
        lbf.init_two_component()
        return lbf

    def test_constant_force_velocity(self):
        velocities = {}
        for label, make_fluid in [("sc", self.single_component_fluid),
                                  ("cg", self.color_gradient_fluid)]:
            lbf = make_fluid(ext_force_density=EXT_FORCE_DENSITY)
            velocities[label] = []
            for _ in range(3):
                self.system.integrator.run(1)
                velocities[label].append(np.copy(lbf[:, :, :].velocity))
            self.tearDown()

        for n, velocity in enumerate(velocities["cg"], start=1):
            np.testing.assert_allclose(
                velocity, np.broadcast_to(
                    (n + 0.5) * TAU * EXT_FORCE_DENSITY / RHO_0,
                    velocity.shape),
                rtol=self.tol(1e-10, 1e-5), atol=0.)
            np.testing.assert_allclose(
                velocity, velocities["sc"][n - 1],
                rtol=self.tol(1e-10, 1e-5), atol=0.)

    def test_friction_coupling(self):
        trajectories = {}
        for label, make_fluid in [("sc", self.single_component_fluid),
                                  ("cg", self.color_gradient_fluid)]:
            lbf = make_fluid()
            p = self.system.part.add(pos=[3.3, 4.1, 3.9],
                                     v=[0.02, 0.01, -0.005])
            self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=1)
            trajectories[label] = []
            for _ in range(N_STEPS):
                self.system.integrator.run(1)
                trajectories[label].append(np.copy(p.v))
            self.tearDown()

        np.testing.assert_allclose(
            trajectories["cg"], trajectories["sc"],
            rtol=0., atol=self.tol(1e-7, 1e-6))


@utx.skipIfMissingFeatures(["WALBERLA"])
class ColorGradientSingleComponentParityTestCPU(
        ColorGradientSingleComponentParityTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class ColorGradientSingleComponentParityTestCPUSinglePrecision(
        ColorGradientSingleComponentParityTest, ut.TestCase):
    lb_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class ColorGradientSingleComponentParityTestGPU(
        ColorGradientSingleComponentParityTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class ColorGradientSingleComponentParityTestGPUSinglePrecision(
        ColorGradientSingleComponentParityTest, ut.TestCase):
    lb_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
