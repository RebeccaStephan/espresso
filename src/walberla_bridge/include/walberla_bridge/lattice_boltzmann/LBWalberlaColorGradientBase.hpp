/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

/**
 * @file
 * Color-gradient (two-component) LB extension interface.
 * Inherits the common fluid interface and adds model-specific
 * configuration, initialization, and per-component accessors.
 */

#include "LBWalberlaBase.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <optional>
#include <vector>

/**
 * @brief Negative component densities found on this MPI rank in one LB step,
 * after streaming and before collision.
 */
struct NegativeDensityReport {
  /** Component of the most negative density (0 = a, 1 = b). */
  int component;
  /** Global node of the most negative density. */
  Utils::Vector3i node;
  /** Most negative density, in lattice units. */
  double density;
  /** Number of nodes on this rank with a negative density of either component. */
  int n_nodes;
  /** LB step in which it was found, counted from 1. */
  unsigned int time_step;
};

class LBWalberlaColorGradientBase : public virtual LBWalberlaBase {
public:
  ~LBWalberlaColorGradientBase() override = default;

  // CG-specific configuration
  virtual void set_collision_model_color_gradient(double sigma,
                                                  double beta) = 0;
  virtual void init_pdfs_from_components() = 0;

  // CG-specific observables / forces
  virtual void
  add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &positions,
                              std::vector<double> const &delta_mus) = 0;
  virtual std::vector<Utils::Vector3d>
  get_solvation_particle_forces_at_pos(
      std::vector<Utils::Vector3d> const &pos,
      std::vector<double> const &delta_mus) = 0;

  // per-component viscosity
  virtual void set_component_viscosities(std::array<double, 2> const &nu) = 0;
  [[nodiscard]] virtual std::array<double, 2>
  get_component_viscosities() const = 0;

  /** @brief Order parameter phi = (rho_a - rho_b) / (rho_a + rho_b). */
  virtual std::optional<double>
  get_node_phasefield(Utils::Vector3i const &node,
                      bool consider_ghosts = false) const = 0;

  /** @brief Overwrite the stored order parameter (checkpointing only). */
  virtual bool set_node_phasefield(Utils::Vector3i const &node,
                                   double phasefield) = 0;

  /**
   * @brief Overwrite the stored barycentric velocity (checkpointing only).
   * Unlike @ref LBWalberlaBase::set_node_velocity, the populations are not
   * modified; this only restores the velocity field written by the stream
   * sweep, which the particle coupling reads before the next LB step.
   */
  virtual bool set_node_velocity_raw(Utils::Vector3i const &node,
                                     Utils::Vector3d const &v) = 0;

  /**
   * @brief Read the stored barycentric velocity (checkpointing only).
   * Unlike @ref LBWalberlaBase::get_node_velocity, boundary nodes return the
   * velocity field value instead of the boundary slip velocity.
   */
  virtual std::optional<Utils::Vector3d>
  get_node_velocity_raw(Utils::Vector3i const &node,
                        bool consider_ghosts = false) const = 0;

  /**
   * @brief Negative component densities found in the most recent LB step on
   * this MPI rank, or @c std::nullopt if there were none.
   */
  [[nodiscard]] virtual std::optional<NegativeDensityReport>
  get_negative_density_report() const = 0;
};
