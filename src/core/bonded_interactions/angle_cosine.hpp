/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

/** \file
 *  Routines to calculate the angle energy or/and and force
 *  for a particle triple using the potential described in
 *  @ref bondedIA_angle_cosine.
 */

#include "BoxGeometry.hpp"
#include "angle_common.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <tuple>

/** Parameters for three-body angular potential (cosine). */
struct AngleCosineBond {
  /** bending constant */
  double bend;
  /** equilibrium angle (default is 180 degrees) */
  double phi0;
  /** cosine of @p phi0 (internal parameter) */
  double cos_phi0;
  /** sine of @p phi0 (internal parameter) */
  double sin_phi0;

  double cutoff() const { return 0.; }

  static constexpr int num = 2;

  AngleCosineBond(double bend, double phi0) {
    this->bend = bend;
    this->phi0 = phi0;
    this->cos_phi0 = cos(phi0);
    this->sin_phi0 = sin(phi0);
  }

  std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
  forces(Utils::Vector3d const &vec1, Utils::Vector3d const &vec2) const;
  double energy(Utils::Vector3d const &vec1, Utils::Vector3d const &vec2) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar & bend;
    ar & phi0;
    ar & cos_phi0;
    ar & sin_phi0;
  }
};

/** Compute the three-body angle interaction force.
 *  @param[in]  vec1  Vector from central particle to left particle.
 *  @param[in]  vec2  Vector from central particle to right particle.
 *  @return Forces on the second, first and third particles, in that order.
 */
inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
AngleCosineBond::forces(Utils::Vector3d const &vec1,
                        Utils::Vector3d const &vec2) const {

  auto forceFactor = [this](double const cos_phi) {
    auto const sin_phi = sqrt(1 - Utils::sqr(cos_phi));
    // angle force term: K(phi) = -k * sin(phi - phi0) / sin(phi)
    // trig identity: sin(a - b) = sin(a)cos(b) - cos(a)sin(b)
    return -bend * (sin_phi * cos_phi0 - cos_phi * sin_phi0) / sin_phi;
  };

  return angle_generic_force(vec1, vec2, forceFactor, false);
}

/** Computes the three-body angle interaction energy.
 *  @param[in]  vec1  Vector from central particle to left particle.
 *  @param[in]  vec2  Vector from central particle to right particle.
 */
inline double AngleCosineBond::energy(Utils::Vector3d const &vec1,
                                      Utils::Vector3d const &vec2) const {
  auto const cos_phi = calc_cosine(vec1, vec2, true);
  auto const sin_phi = sqrt(1 - Utils::sqr(cos_phi));
  // potential: U(phi) = k * [1 - cos(phi - phi0)]
  // trig identity: cos(phi - phi0) = cos(phi)cos(phi0) + sin(phi)sin(phi0)
  return bend * (1 - (cos_phi * cos_phi0 + sin_phi * sin_phi0));
}
