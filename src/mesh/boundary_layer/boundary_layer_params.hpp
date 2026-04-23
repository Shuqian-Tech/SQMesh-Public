// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace sqmesh::mesh {

struct BoundaryLayerParams {
  enum class FirstHeightMode { Absolute, Aspect };

  // -- Size Control ---------------------------------------------------
  int num_layers = 5;
  double first_height = 0.01;
  double first_height_aspect = 0.5;
  double growth_rate = 1.2;
  std::vector<double> variable_growth_rates; // Per-layer growth rates
  int additional_outer_layers = 0;
  FirstHeightMode first_height_mode = FirstHeightMode::Absolute;

  // Orthogonal layers (strict normal, no smoothing) before transitioning
  int orthogonal_layers = 0;

  // -- Normal/Vector Control ------------------------------------------
  bool smooth_normals = true;
  double max_normal_deviation_deg = 60.0;
  int smooth_iterations = 5;
  bool separate_vectors_at_sharp_angles = false;
  double separation_angle_deg = 165.0;

  // -- Quality Control ------------------------------------------------
  double max_skewness = 0.98;
  double max_warping_deg = 30.0;
  double min_first_height = 0.005;
  double min_layer_aspect = 0.01;
  double max_layer_aspect = 100.0;

  // -- Problematic Area Treatment -------------------------------------
  enum class ProblematicAreaTreatment { Collapse, Exclude, Stop };
  ProblematicAreaTreatment treatment = ProblematicAreaTreatment::Collapse;
  bool allow_squeeze = false;
  int max_consecutive_collapsed_layers = 3;

  // -- Proximity Control ----------------------------------------------
  bool allow_proximity = true;
  double proximity_factor = 0.4;

  // -- Side Connection ------------------------------------------------
  bool adjust_to_neighbors = true;
  double adjust_connect_angle_limit_deg = 85.0;

  // -- y+ Calculation -------------------------------------------------
  double target_y_plus = 0.0; // >0 enables automatic first height
  double freestream_velocity = 0.0;
  double density = 1.205;          // kg/m³ (air at sea level)
  double dynamic_viscosity = 1.82e-5; // kg/(m·s)
  double reference_length = 1.0;   // m (characteristic length)

  // -- Quality gate thresholds ----------------------------------------
  bool retain_valid_layers = true; // Keep good layers even if later ones fail

  [[nodiscard]] bool use_aspect_first_height() const noexcept
  {
    return first_height_mode == FirstHeightMode::Aspect;
  }

  [[nodiscard]] double layer_growth_factor(int layer_index) const noexcept
  {
    if(layer_index <= 0) return 1.0;
    if(!variable_growth_rates.empty()) {
      double factor = 1.0;
      for(int i = 0; i < layer_index; ++i) {
        const auto gr_idx = static_cast<std::size_t>(i);
        const double gr = gr_idx < variable_growth_rates.size()
          ? variable_growth_rates[gr_idx]
          : variable_growth_rates.back();
        factor *= gr;
      }
      return factor;
    }
    return std::pow(growth_rate, layer_index);
  }

  [[nodiscard]] double layer_height_from_base(
    double base_first_height,
    int layer_index
  ) const noexcept
  {
    return base_first_height * layer_growth_factor(layer_index);
  }

  // -- Computed layer height for a given layer index ------------------
  [[nodiscard]] double layer_height(int layer_index) const noexcept
  {
    if(layer_index < 0) return first_height;
    return layer_height_from_base(first_height, layer_index);
  }

  [[nodiscard]] double total_height(int num) const noexcept
  {
    double total = 0.0;
    for(int i = 0; i < num; ++i) {
      total += layer_height(i);
    }
    return total;
  }

  // -- y+ -> first cell height (Schlichting flat-plate correlation) ------
  static double compute_first_height_from_yplus(
    double y_plus,
    double velocity,
    double rho,
    double mu,
    double ref_length
  ) noexcept
  {
    if(velocity <= 0.0 || rho <= 0.0 || mu <= 0.0 || ref_length <= 0.0) {
      return 0.01;
    }
    const double Re = (rho * velocity * ref_length) / mu;
    const double Cf = std::pow(2.0 * std::log10(Re) - 0.65, -2.3);
    const double tau_w = Cf * 0.5 * rho * velocity * velocity;
    const double u_tau = std::sqrt(tau_w / rho);
    return (y_plus * mu) / (rho * u_tau);
  }

  // Apply y+ auto-calculation if target_y_plus > 0
  void apply_yplus_auto_height() noexcept
  {
    if(target_y_plus <= 0.0 || use_aspect_first_height()) return;
    double h = compute_first_height_from_yplus(
      target_y_plus, freestream_velocity, density, dynamic_viscosity, reference_length
    );
    if(h > 0.0 && std::isfinite(h)) {
      first_height = h;
      min_first_height = h * 0.5;
    }
  }
};

} // namespace sqmesh::mesh
