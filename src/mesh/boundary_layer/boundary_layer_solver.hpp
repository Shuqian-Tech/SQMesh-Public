// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include <array>
#include "boundary_layer_params.hpp"
#include "boundary_layer_work_data.hpp"

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace sqmesh::mesh {

// BoundaryLayerSolver: core boundary layer mesh generation algorithm.
//
// Pipeline:
//   1. Initialize from surface mesh + non-BL faces
//   2. Compute face normals
//   3. Classify node types (flat/concave/convex)
//   4. Build the node march field (normals, distances, smoothing, vectors)
//   5. Per-layer loop:
//      a. Build the march field for the current front with collision/proximity limits
//      b. Propagate the active front by one layer
//      c. Quality check (skewness, warping)
//      d. Handle problematic areas (collapse/exclude/stop)
//      e. Update active front
//   6. Extract top-cap faces for subsequent tet fill
class BoundaryLayerSolver {
public:
  // Set parameters
  void set_params(const BoundaryLayerParams &params);

  // Set force flip: when true, flip face winding before generation.
  // When false, use auto-detection (centroid method).
  // Call before generate().
  void set_force_flip(bool flip) { force_flip_set_ = true; force_flip_ = flip; }

  // Set target point: BL pushes toward this point's region.
  // Per-face flip: each face normal is oriented to point toward target_point.
  void set_target_point(double x, double y, double z) {
    target_point_ = {x, y, z};
    has_target_point_ = true;
  }

  // Load BL surface faces from the working data (must be populated already)
  void initialize(BLWorkData &data);

  // Run the full layer generation loop
  void generate();

  // Access results
  [[nodiscard]] const BLWorkData &work_data() const noexcept { return *data_; }
  [[nodiscard]] BLWorkData &work_data() noexcept { return *data_; }

private:
  // -- Stage functions --

  // Compute face normals for active faces
  void update_face_normals();

  // Detect boundary nodes, classify node types, and update node angles
  void classify_nodes_and_update_angles();

  // Build the node march field for the current front.
  void build_node_march_field();

  // Build the per-node base first-layer heights.
  void initialize_base_first_heights();

  // Initialize thickness per node from layer height
  void initialize_thickness_of_nodes();

  // Smooth thickness values (Laplacian)
  void smooth_thickness();

  // Squeeze thickness by shortest distance (proximity/collision)
  void squeeze_thickness_by_shortest_distance();

  // Compute march vector = direction * thickness
  void compute_march_vectors();

  // Propagate the current front by one layer, validate the candidate layer,
  // apply problematic-area treatment, and update active state.
  // Returns false when generation should stop after retaining valid layers.
  [[nodiscard]] bool propagate_front_one_layer(int layer);

  // Deactivate faces/nodes whose valid layer budget has already been exhausted.
  void deactivate_reached_front_entities(int layer);

  // Node-corridor intersection screen before accepting a pushed layer.
  [[nodiscard]] bool check_node_corridor_intersections(
    const std::vector<Vec3> &current_positions,
    const std::vector<Vec3> &next_positions,
    std::unordered_set<std::uint32_t> &modify_nodes,
    std::unordered_set<std::uint32_t> &modify_faces
  ) const;

  // Supplemental geometric intersection screen on the pushed top front.
  [[nodiscard]] bool check_top_face_intersections(
    const std::vector<Vec3> &top_positions,
    std::unordered_set<std::uint32_t> &modify_nodes,
    std::unordered_set<std::uint32_t> &modify_faces
  ) const;

  // Geometric intersection screen on the newly swept side-wall quads.
  [[nodiscard]] bool check_side_face_intersections(
    const std::vector<Vec3> &bottom_positions,
    const std::vector<Vec3> &top_positions,
    std::unordered_set<std::uint32_t> &modify_nodes,
    std::unordered_set<std::uint32_t> &modify_faces
  ) const;

  // Quality check on top faces and newly created prism/hexa cells.
  [[nodiscard]] bool check_top_face_quality(
    const std::vector<Vec3> &bottom_positions,
    const std::vector<Vec3> &top_positions,
    std::unordered_set<std::uint32_t> &modify_nodes,
    std::unordered_set<std::uint32_t> &modify_faces
  ) const;

  // Apply Collapse / Exclude handling to nodes and faces marked as problematic.
  void apply_problematic_area_treatment(
    int layer,
    std::unordered_set<std::uint32_t> &modify_nodes,
    std::unordered_set<std::uint32_t> &modify_faces
  );

  void enforce_max_consecutive_collapsed_layers(
    const std::unordered_set<std::uint32_t> &seed_faces
  );

  void recompute_node_thickness();
  [[nodiscard]] double compute_node_visibility_angle(std::uint32_t node) const;

  void clear_problem_marks();
  void mark_problem_face(std::uint32_t face_id);
  void build_active_face_bvh(
    const std::vector<Vec3> &positions,
    FaceBVH &bvh,
    std::vector<std::array<std::uint32_t, 4>> &face_verts,
    std::vector<int> &face_nv,
    std::vector<std::uint32_t> &local_to_global,
    std::vector<int> &global_to_local
  ) const;
  void export_first_problem_debug(
    const char *reason,
    std::int32_t face_id,
    const std::vector<std::array<Vec3, 4>> &faces,
    const std::vector<int> &face_num_vertices,
    const std::string &details_markdown,
    const std::vector<Vec3> &polyline
  ) const;

  // -- Normal computation helpers --
  Vec3 compute_face_normal_tri(std::uint32_t fi) const;
  Vec3 compute_face_normal_quad(std::uint32_t fi) const;

  // Compute the average node normal from adjacent active face normals.
  Vec3 compute_average_node_normal(std::uint32_t node) const;

  // Iterative normal direction.
  Vec3 iterate_by_angle(std::uint32_t node) const;

  // Smooth normals using Laplacian + angle-weighted approach
  void smooth_normals_laplacian(std::vector<Vec3> &normals);

  // Distance-based growth coefficient (alpha for concave/convex)
  std::vector<double> compute_distance_coefficient() const;

  [[nodiscard]] double layer_height_for_node(
    std::uint32_t node,
    int layer_index
  ) const noexcept;

  // -- Quality metrics --
  static double compute_skewness_tri(
    const Vec3 &a, const Vec3 &b, const Vec3 &c);
  static double compute_skewness_quad(
    const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d);
  static double compute_warping_quad(
    const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d);
  static double signed_corner_jacobian(
    const Vec3 &e0, const Vec3 &e1, const Vec3 &e2);
  static bool prism_has_negative_jacobian(
    const Vec3 &b0, const Vec3 &b1, const Vec3 &b2,
    const Vec3 &t0, const Vec3 &t1, const Vec3 &t2);
  static bool hexa_has_negative_jacobian(
    const Vec3 &b0, const Vec3 &b1, const Vec3 &b2, const Vec3 &b3,
    const Vec3 &t0, const Vec3 &t1, const Vec3 &t2, const Vec3 &t3);
  static double compute_prism_skewness(
    const Vec3 &b0, const Vec3 &b1, const Vec3 &b2,
    const Vec3 &t0, const Vec3 &t1, const Vec3 &t2);
  static double compute_hexa_skewness(
    const Vec3 &b0, const Vec3 &b1, const Vec3 &b2, const Vec3 &b3,
    const Vec3 &t0, const Vec3 &t1, const Vec3 &t2, const Vec3 &t3);

  BoundaryLayerParams params_;
  BLWorkData *data_ = nullptr;
  int current_layer_ = 0;
  bool normals_flipped_ = false; // true if face normals were flipped to point inward
  bool force_flip_set_ = false;  // true if set_force_flip was called
  bool force_flip_ = false;      // the value passed by set_force_flip
  Vec3 target_point_ = {0, 0, 0};
  bool has_target_point_ = false;

  // Temporary storage for normals during smoothing
  std::vector<Vec3> initial_normals_;
  std::vector<Vec3> final_normals_;
  std::vector<double> base_first_heights_;
  std::vector<double> initial_distances_;
  std::vector<double> final_distances_;
  std::vector<bool> problem_nodes_;
  std::vector<bool> problem_faces_;
};

} // namespace sqmesh::mesh
