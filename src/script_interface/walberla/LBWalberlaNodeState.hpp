/*
 * Copyright (C) 2021-2026 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <utils/Vector.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>

/** Checkpoint data for a LB node. */
struct LBWalberlaNodeState {
  /** One stencil per component (SC: 19 values, CG: 2*19 values). */
  std::vector<double> populations;
  /** One force per component. */
  std::vector<Utils::Vector3d> last_applied_force;
  /** Two-component only: per-component densities. */
  std::vector<double> density;
  /** Two-component only: order parameter. */
  double phasefield = 0.;
  /** Two-component only: stored barycentric velocity. */
  Utils::Vector3d velocity = {};
  Utils::Vector3d slip_velocity = {};
  bool is_boundary = false;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, unsigned int /* version */) {
    ar & populations & last_applied_force & density & phasefield & velocity &
        slip_velocity & is_boundary;
  }
};
