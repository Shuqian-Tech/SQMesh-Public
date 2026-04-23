// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "boundary_layer_solver.hpp"
#include "core/log.hpp"

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <map>
#include <numeric>
#include <sstream>
#include <string_view>
#include <stdexcept>
#include <unordered_set>

namespace sqmesh::mesh {

	namespace {

		constexpr double kPI = 3.14159265358979323846;
		constexpr double kHALF_PI = kPI / 2.0;
		constexpr double kEPS = 1e-12;
		constexpr double kVISIANGLE_EPS = 1e-10;

		constexpr double BL3d_CONCAVE_ANGLE = kPI * 0.85;  // ~153°
		constexpr double BL3d_CONVEX_ANGLE = kPI * 1.15;  // ~207°
		constexpr double BL_ANGLE_EPS = 1e-6;
		constexpr double kIntersectionEps = 1e-9;

		// Default layer ratio for squeeze
		constexpr double BL_DEFAULT_LAYER_RATIO = 0.5;

		void bl_log_vprint(const char* level, const char* fmt, std::va_list args)
		{
			std::fprintf(stderr, "[%s] ", level);
			std::vfprintf(stderr, fmt, args);
			std::fputc('\n', stderr);
		}

		void bl_log_info(const char* fmt, ...)
		{
			std::va_list args;
			va_start(args, fmt);
			bl_log_vprint("info", fmt, args);
			va_end(args);
		}

		void bl_log_warn(const char* fmt, ...)
		{
			std::va_list args;
			va_start(args, fmt);
			bl_log_vprint("warn", fmt, args);
			va_end(args);
		}

		double clamp_acos(double val) noexcept
		{
			return std::acos(std::clamp(val, -1.0, 1.0));
		}

		Vec3 calc_face_normal_from_three(const Vec3& a, const Vec3& b, const Vec3& c) noexcept
		{
			auto e1 = vec3_sub(b, a);
			auto e2 = vec3_sub(c, a);
			return vec3_normalized(vec3_cross(e1, e2));
		}

		double get_vec3d_angle(const Vec3& a, const Vec3& b) noexcept
		{
			double d = vec3_dot(a, b);
			d = std::clamp(d, -1.0, 1.0);
			return std::acos(d);
		}

		// Compute angle between two face normals, considering orientation
		double get_face_dihedral_angle(
			const Vec3& n1, const Vec3& n2,
			const Vec3& edge_dir) noexcept
		{
			// Dihedral angle = the angle between two faces measured from the inside.
			// For outward-pointing normals:
			//   - coplanar faces: dot(n1,n2) = 1 -> angle_between = 0 -> dihedral = π
			//   - 90° convex edge: dot(n1,n2) = 0 -> angle_between = π/2 -> dihedral = π/2
			//   - 90° concave edge: dot(n1,n2) = 0 -> angle_between = π/2 -> dihedral = 3π/2
			//
			// Formula: dihedral = π + atan2(dot(cross(n1,n2), edge_dir), dot(n1,n2))
			// This gives [0, 2π): < π = convex, = π = flat, > π = concave

			double dot = vec3_dot(n1, n2);
			auto cross_n = vec3_cross(n1, n2);
			double cross_dot = vec3_dot(cross_n, edge_dir);

			// atan2 gives signed angle in [-π, π]
			double signed_angle = std::atan2(cross_dot, dot);

			// Dihedral = π - signed_angle -> maps to [0, 2π]
			double dihedral = kPI - signed_angle;

			// Clamp to [0, 2π]
			if (dihedral < 0) dihedral += 2.0 * kPI;
			if (dihedral > 2.0 * kPI) dihedral -= 2.0 * kPI;

			return dihedral;
		}

		struct Triangle3 {
			Vec3 a;
			Vec3 b;
			Vec3 c;
		};

		struct ExplicitFace {
			std::array<Vec3, 4> points{};
			int num_vertices = 0;
		};

		[[nodiscard]] std::array<Vec3, 4> make_debug_face(
			const Vec3& p0,
			const Vec3& p1,
			const Vec3& p2
		) noexcept
		{
			std::array<Vec3, 4> face{};
			face[0] = p0;
			face[1] = p1;
			face[2] = p2;
			return face;
		}

		[[nodiscard]] std::array<Vec3, 4> make_debug_face(
			const Vec3& p0,
			const Vec3& p1,
			const Vec3& p2,
			const Vec3& p3
		) noexcept
		{
			std::array<Vec3, 4> face{};
			face[0] = p0;
			face[1] = p1;
			face[2] = p2;
			face[3] = p3;
			return face;
		}

		[[nodiscard]] std::array<Vec3, 4> copy_face_points(
			const std::vector<Vec3>& positions,
			const SurfaceFace& face
		) noexcept
		{
			std::array<Vec3, 4> points{};
			for (int i = 0; i < face.num_vertices; ++i) {
				points[static_cast<std::size_t>(i)] = positions[face.vertices[i]];
			}
			return points;
		}

		[[nodiscard]] std::array<Vec3, 4> copy_face_points(const ExplicitFace& face) noexcept
		{
			return face.points;
		}

		[[nodiscard]] std::string format_vec3(const Vec3& value)
		{
			std::ostringstream out;
			out << std::fixed << std::setprecision(10)
				<< value[0] << ", " << value[1] << ", " << value[2];
			return out.str();
		}

		[[nodiscard]] double triangle_area(const Triangle3& tri) noexcept
		{
			return 0.5 * vec3_length(vec3_cross(
				vec3_sub(tri.b, tri.a),
				vec3_sub(tri.c, tri.a)));
		}

		[[nodiscard]] AABB make_face_aabb(
			const std::vector<Vec3>& positions,
			const SurfaceFace& face
		) noexcept
		{
			AABB box;
			for (int i = 0; i < face.num_vertices; ++i) {
				box.expand(positions[face.vertices[i]]);
			}
			box.pad(kIntersectionEps);
			return box;
		}

		[[nodiscard]] AABB make_explicit_face_aabb(
			const ExplicitFace& face
		) noexcept
		{
			AABB box;
			for (int i = 0; i < face.num_vertices; ++i) {
				box.expand(face.points[static_cast<std::size_t>(i)]);
			}
			box.pad(kIntersectionEps);
			return box;
		}

		[[nodiscard]] bool faces_share_vertex(
			const SurfaceFace& lhs,
			const SurfaceFace& rhs
		) noexcept
		{
			for (int li = 0; li < lhs.num_vertices; ++li) {
				for (int ri = 0; ri < rhs.num_vertices; ++ri) {
					if (lhs.vertices[li] == rhs.vertices[ri]) {
						return true;
					}
				}
			}
			return false;
		}

		[[nodiscard]] bool face_contains_node(
			const SurfaceFace& face,
			std::uint32_t node
		) noexcept
		{
			for (int vi = 0; vi < face.num_vertices; ++vi) {
				if (face.vertices[vi] == node) {
					return true;
				}
			}
			return false;
		}

		[[nodiscard]] int face_triangle_count(const SurfaceFace& face) noexcept
		{
			return face.num_vertices == 4 ? 2 : 1;
		}

		[[nodiscard]] int explicit_face_triangle_count(const ExplicitFace& face) noexcept
		{
			return face.num_vertices == 4 ? 2 : 1;
		}

		[[nodiscard]] Triangle3 face_triangle(
			const std::vector<Vec3>& positions,
			const SurfaceFace& face,
			int tri_index
		) noexcept
		{
			if (face.num_vertices == 4 && tri_index != 0) {
				return {
				  positions[face.vertices[0]],
				  positions[face.vertices[2]],
				  positions[face.vertices[3]]
				};
			}

			return {
			  positions[face.vertices[0]],
			  positions[face.vertices[1]],
			  positions[face.vertices[2]]
			};
		}

		[[nodiscard]] Triangle3 explicit_face_triangle(
			const ExplicitFace& face,
			int tri_index
		) noexcept
		{
			if (face.num_vertices == 4 && tri_index != 0) {
				return {
				  face.points[0],
				  face.points[2],
				  face.points[3]
				};
			}

			return {
			  face.points[0],
			  face.points[1],
			  face.points[2]
			};
		}

		[[nodiscard]] double explicit_face_area(const ExplicitFace& face) noexcept
		{
			double area = 0.0;
			for (int tri_index = 0; tri_index < explicit_face_triangle_count(face); ++tri_index) {
				area += triangle_area(explicit_face_triangle(face, tri_index));
			}
			return area;
		}

		[[nodiscard]] std::array<double, 2> project_point_2d(
			const Vec3& point,
			int drop_axis
		) noexcept
		{
			switch (drop_axis) {
			case 0: return { point[1], point[2] };
			case 1: return { point[0], point[2] };
			default: return { point[0], point[1] };
			}
		}

		[[nodiscard]] double orient_2d(
			const std::array<double, 2>& a,
			const std::array<double, 2>& b,
			const std::array<double, 2>& c
		) noexcept
		{
			return (b[0] - a[0]) * (c[1] - a[1]) -
				(b[1] - a[1]) * (c[0] - a[0]);
		}

		[[nodiscard]] bool point_in_triangle_2d_strict(
			const std::array<double, 2>& p,
			const std::array<double, 2>& a,
			const std::array<double, 2>& b,
			const std::array<double, 2>& c
		) noexcept
		{
			const double ab = orient_2d(a, b, p);
			const double bc = orient_2d(b, c, p);
			const double ca = orient_2d(c, a, p);

			const bool positive =
				ab > kIntersectionEps && bc > kIntersectionEps && ca > kIntersectionEps;
			const bool negative =
				ab < -kIntersectionEps && bc < -kIntersectionEps && ca < -kIntersectionEps;
			return positive || negative;
		}

		[[nodiscard]] bool segments_intersect_2d_strict(
			const std::array<double, 2>& a0,
			const std::array<double, 2>& a1,
			const std::array<double, 2>& b0,
			const std::array<double, 2>& b1
		) noexcept
		{
			const double o1 = orient_2d(a0, a1, b0);
			const double o2 = orient_2d(a0, a1, b1);
			const double o3 = orient_2d(b0, b1, a0);
			const double o4 = orient_2d(b0, b1, a1);

			const bool separated =
				(o1 > kIntersectionEps && o2 > kIntersectionEps) ||
				(o1 < -kIntersectionEps && o2 < -kIntersectionEps) ||
				(o3 > kIntersectionEps && o4 > kIntersectionEps) ||
				(o3 < -kIntersectionEps && o4 < -kIntersectionEps);
			if (separated) {
				return false;
			}

			const bool straddles_a =
				(o1 > kIntersectionEps && o2 < -kIntersectionEps) ||
				(o1 < -kIntersectionEps && o2 > kIntersectionEps);
			const bool straddles_b =
				(o3 > kIntersectionEps && o4 < -kIntersectionEps) ||
				(o3 < -kIntersectionEps && o4 > kIntersectionEps);
			return straddles_a && straddles_b;
		}

		[[nodiscard]] bool coplanar_triangles_intersect_strict(
			const Triangle3& lhs,
			const Triangle3& rhs,
			const Vec3& plane_normal
		) noexcept
		{
			int drop_axis = 0;
			double max_abs = std::abs(plane_normal[0]);
			for (int axis = 1; axis < 3; ++axis) {
				const double cur_abs = std::abs(plane_normal[axis]);
				if (cur_abs > max_abs) {
					max_abs = cur_abs;
					drop_axis = axis;
				}
			}

			const auto la = project_point_2d(lhs.a, drop_axis);
			const auto lb = project_point_2d(lhs.b, drop_axis);
			const auto lc = project_point_2d(lhs.c, drop_axis);
			const auto ra = project_point_2d(rhs.a, drop_axis);
			const auto rb = project_point_2d(rhs.b, drop_axis);
			const auto rc = project_point_2d(rhs.c, drop_axis);

			const std::array<std::pair<std::array<double, 2>, std::array<double, 2>>, 3> lhs_edges = {
			  std::pair{la, lb},
			  std::pair{lb, lc},
			  std::pair{lc, la}
			};
			const std::array<std::pair<std::array<double, 2>, std::array<double, 2>>, 3> rhs_edges = {
			  std::pair{ra, rb},
			  std::pair{rb, rc},
			  std::pair{rc, ra}
			};

			for (const auto& lhs_edge : lhs_edges) {
				for (const auto& rhs_edge : rhs_edges) {
					if (segments_intersect_2d_strict(
						lhs_edge.first, lhs_edge.second,
						rhs_edge.first, rhs_edge.second)) {
						return true;
					}
				}
			}

			return point_in_triangle_2d_strict(la, ra, rb, rc) ||
				point_in_triangle_2d_strict(ra, la, lb, lc);
		}

		[[nodiscard]] bool segment_triangle_intersect_strict(
			const Vec3& p0,
			const Vec3& p1,
			const Triangle3& tri
		) noexcept
		{
			const Vec3 dir = vec3_sub(p1, p0);
			const Vec3 e1 = vec3_sub(tri.b, tri.a);
			const Vec3 e2 = vec3_sub(tri.c, tri.a);
			const Vec3 h = vec3_cross(dir, e2);
			const double det = vec3_dot(e1, h);
			if (std::abs(det) < kIntersectionEps) {
				return false;
			}

			const double inv_det = 1.0 / det;
			const Vec3 s = vec3_sub(p0, tri.a);
			const double u = inv_det * vec3_dot(s, h);
			if (u <= kIntersectionEps || u >= 1.0 - kIntersectionEps) {
				return false;
			}

			const Vec3 q = vec3_cross(s, e1);
			const double v = inv_det * vec3_dot(dir, q);
			if (v <= kIntersectionEps || u + v >= 1.0 - kIntersectionEps) {
				return false;
			}

			const double t = inv_det * vec3_dot(e2, q);
			return t > kIntersectionEps && t < 1.0 - kIntersectionEps;
		}

		[[nodiscard]] bool triangles_intersect_strict(
			const Triangle3& lhs,
			const Triangle3& rhs
		) noexcept
		{
			const Vec3 lhs_normal = vec3_cross(vec3_sub(lhs.b, lhs.a), vec3_sub(lhs.c, lhs.a));
			const Vec3 rhs_normal = vec3_cross(vec3_sub(rhs.b, rhs.a), vec3_sub(rhs.c, rhs.a));

			if (vec3_length(lhs_normal) < kIntersectionEps ||
				vec3_length(rhs_normal) < kIntersectionEps) {
				return false;
			}

			const std::array<double, 3> rhs_to_lhs = {
			  vec3_dot(lhs_normal, vec3_sub(rhs.a, lhs.a)),
			  vec3_dot(lhs_normal, vec3_sub(rhs.b, lhs.a)),
			  vec3_dot(lhs_normal, vec3_sub(rhs.c, lhs.a))
			};
			const std::array<double, 3> lhs_to_rhs = {
			  vec3_dot(rhs_normal, vec3_sub(lhs.a, rhs.a)),
			  vec3_dot(rhs_normal, vec3_sub(lhs.b, rhs.a)),
			  vec3_dot(rhs_normal, vec3_sub(lhs.c, rhs.a))
			};

			const auto same_side = [](const std::array<double, 3>& distances) noexcept {
				bool has_pos = false;
				bool has_neg = false;
				for (double d : distances) {
					if (d > kIntersectionEps) has_pos = true;
					if (d < -kIntersectionEps) has_neg = true;
				}
				return !(has_pos && has_neg);
				};

			const bool coplanar =
				std::abs(rhs_to_lhs[0]) <= kIntersectionEps &&
				std::abs(rhs_to_lhs[1]) <= kIntersectionEps &&
				std::abs(rhs_to_lhs[2]) <= kIntersectionEps &&
				std::abs(lhs_to_rhs[0]) <= kIntersectionEps &&
				std::abs(lhs_to_rhs[1]) <= kIntersectionEps &&
				std::abs(lhs_to_rhs[2]) <= kIntersectionEps;
			if (coplanar) {
				return coplanar_triangles_intersect_strict(lhs, rhs, lhs_normal);
			}

			if (same_side(rhs_to_lhs) || same_side(lhs_to_rhs)) {
				return false;
			}

			const std::array<std::pair<Vec3, Vec3>, 3> lhs_edges = {
			  std::pair{lhs.a, lhs.b},
			  std::pair{lhs.b, lhs.c},
			  std::pair{lhs.c, lhs.a}
			};
			const std::array<std::pair<Vec3, Vec3>, 3> rhs_edges = {
			  std::pair{rhs.a, rhs.b},
			  std::pair{rhs.b, rhs.c},
			  std::pair{rhs.c, rhs.a}
			};

			for (const auto& edge : lhs_edges) {
				if (segment_triangle_intersect_strict(edge.first, edge.second, rhs)) {
					return true;
				}
			}
			for (const auto& edge : rhs_edges) {
				if (segment_triangle_intersect_strict(edge.first, edge.second, lhs)) {
					return true;
				}
			}

			return false;
		}

		[[nodiscard]] bool faces_intersect_strict(
			const std::vector<Vec3>& lhs_positions,
			const SurfaceFace& lhs_face,
			const std::vector<Vec3>& rhs_positions,
			const SurfaceFace& rhs_face
		) noexcept
		{
			for (int li = 0; li < face_triangle_count(lhs_face); ++li) {
				const Triangle3 lhs_tri = face_triangle(lhs_positions, lhs_face, li);
				for (int ri = 0; ri < face_triangle_count(rhs_face); ++ri) {
					const Triangle3 rhs_tri = face_triangle(rhs_positions, rhs_face, ri);
					if (triangles_intersect_strict(lhs_tri, rhs_tri)) {
						return true;
					}
				}
			}

			return false;
		}

		[[nodiscard]] bool explicit_face_intersects_strict(
			const ExplicitFace& lhs_face,
			const std::vector<Vec3>& rhs_positions,
			const SurfaceFace& rhs_face
		) noexcept
		{
			for (int li = 0; li < explicit_face_triangle_count(lhs_face); ++li) {
				const Triangle3 lhs_tri = explicit_face_triangle(lhs_face, li);
				for (int ri = 0; ri < face_triangle_count(rhs_face); ++ri) {
					const Triangle3 rhs_tri = face_triangle(rhs_positions, rhs_face, ri);
					if (triangles_intersect_strict(lhs_tri, rhs_tri)) {
						return true;
					}
				}
			}

			return false;
		}

		[[nodiscard]] bool explicit_faces_intersect_strict(
			const ExplicitFace& lhs_face,
			const ExplicitFace& rhs_face
		) noexcept
		{
			for (int li = 0; li < explicit_face_triangle_count(lhs_face); ++li) {
				const Triangle3 lhs_tri = explicit_face_triangle(lhs_face, li);
				for (int ri = 0; ri < explicit_face_triangle_count(rhs_face); ++ri) {
					const Triangle3 rhs_tri = explicit_face_triangle(rhs_face, ri);
					if (triangles_intersect_strict(lhs_tri, rhs_tri)) {
						return true;
					}
				}
			}

			return false;
		}

	} // namespace

	// ===========================================================================
	// Public interface
	// ===========================================================================

	void BoundaryLayerSolver::set_params(const BoundaryLayerParams& params)
	{
		params_ = params;
		params_.apply_yplus_auto_height();
	}

	void BoundaryLayerSolver::initialize(BLWorkData& data)
	{
		data_ = &data;

		const auto nn = data_->surface_nodes.size();
		const auto nf = data_->surface_faces.size();

		// Build connectivity
		data_->build_connectivity();

		// Initialize state arrays
		data_->node_active.assign(nn, true);
		data_->is_boundary_node.assign(nn, false);
		data_->node_type.assign(nn, NodeTypeBL::FLAT);
		data_->march_direction.assign(nn, { 0, 0, 0 });
		data_->max_node_level.assign(nn, params_.num_layers);
		data_->node_thickness.assign(nn, 0.0);
		data_->node_dihedral_angles.assign(nn, {});
		data_->node_manifold_angle.assign(nn, kPI);

		data_->face_active.assign(nf, true);
		data_->face_normal.assign(nf, { 0, 0, 0 });
		data_->max_face_level.assign(nf, params_.num_layers);

		// Layer 0 = original surface
		data_->layer_positions.clear();
		data_->layer_positions.push_back(data_->surface_nodes);

		initial_normals_.resize(nn, { 0, 0, 0 });
		final_normals_.resize(nn, { 0, 0, 0 });
		base_first_heights_.resize(nn, params_.first_height);
		initial_distances_.resize(nn, 0.0);
		final_distances_.resize(nn, 0.0);
		problem_nodes_.assign(nn, false);
		problem_faces_.assign(nf, false);

		bl_log_info(
			"BL solver initialized: %zu nodes, %zu faces, %d layers",
			nn,
			nf,
			params_.num_layers);
	}

	void BoundaryLayerSolver::clear_problem_marks()
	{
		std::fill(problem_nodes_.begin(), problem_nodes_.end(), false);
		std::fill(problem_faces_.begin(), problem_faces_.end(), false);
	}

	void BoundaryLayerSolver::mark_problem_face(std::uint32_t face_id)
	{
		if (face_id >= problem_faces_.size()) {
			return;
		}

		problem_faces_[face_id] = true;
		const auto& face = data_->surface_faces[face_id];
		for (int i = 0; i < face.num_vertices; ++i) {
			problem_nodes_[face.vertices[i]] = true;
		}
	}

	void BoundaryLayerSolver::export_first_problem_debug(
		const char* reason,
		std::int32_t face_id,
		const std::vector<std::array<Vec3, 4>>& faces,
		const std::vector<int>& face_num_vertices,
		const std::string& details_markdown,
		const std::vector<Vec3>& polyline
	) const
	{
		static_cast<void>(reason);
		static_cast<void>(face_id);
		static_cast<void>(faces);
		static_cast<void>(face_num_vertices);
		static_cast<void>(details_markdown);
		static_cast<void>(polyline);
	}

	void BoundaryLayerSolver::build_active_face_bvh(
		const std::vector<Vec3>& positions,
		FaceBVH& bvh,
		std::vector<std::array<std::uint32_t, 4>>& face_verts,
		std::vector<int>& face_nv,
		std::vector<std::uint32_t>& local_to_global,
		std::vector<int>& global_to_local
	) const
	{
		const auto nf = data_->surface_faces.size();
		face_verts.clear();
		face_nv.clear();
		local_to_global.clear();
		global_to_local.assign(nf, -1);

		face_verts.reserve(nf);
		face_nv.reserve(nf);
		local_to_global.reserve(nf);

		for (std::uint32_t fi = 0; fi < nf; ++fi) {
			if (!data_->face_active[fi]) {
				continue;
			}

			const auto local_id = static_cast<int>(face_verts.size());
			global_to_local[fi] = local_id;
			local_to_global.push_back(fi);
			face_verts.push_back(data_->surface_faces[fi].vertices);
			face_nv.push_back(data_->surface_faces[fi].num_vertices);
		}

		if (!face_verts.empty()) {
			bvh.build(positions, face_verts, face_nv);
		}
	}

	void BoundaryLayerSolver::generate()
	{
		auto t_start = std::chrono::high_resolution_clock::now();
		bl_log_info("Boundary layer mesh generation starting...");

		// -- Ensure face normals point toward the BL target region --
		{
			const auto nf = data_->surface_faces.size();
			int flipped_count = 0;

			if (has_target_point_) {
				// Connected-component approach:
				// 1. Find connected components (faces sharing edges)
				// 2. For each component, compute signed volume → current normal orientation
				// 3. Check if target_point is inside/outside using winding number
				// 4. target inside → normals inward; target outside → normals outward
				// 5. Flip entire component if needed

				// BFS to find connected components
				std::vector<int> comp_id(nf, -1);
				int num_comps = 0;
				for (std::uint32_t seed = 0; seed < nf; ++seed) {
					if (comp_id[seed] >= 0) continue;
					int cid = num_comps++;
					std::vector<std::uint32_t> queue = { seed };
					comp_id[seed] = cid;
					for (std::size_t qi = 0; qi < queue.size(); ++qi) {
						auto fi = queue[qi];
						for (auto nfi : data_->face_to_faces[fi]) {
							if (comp_id[nfi] < 0) {
								comp_id[nfi] = cid;
								queue.push_back(nfi);
							}
						}
					}
				}

				// For each component: signed volume + winding number test
				for (int cid = 0; cid < num_comps; ++cid) {
					// Compute signed volume: V = (1/6) Σ v0·(v1×v2)
					// V > 0 → normals outward, V < 0 → normals inward
					double signed_vol = 0;
					for (std::uint32_t fi = 0; fi < nf; ++fi) {
						if (comp_id[fi] != cid) continue;
						const auto& face = data_->surface_faces[fi];
						const auto& v0 = data_->surface_nodes[face.vertices[0]];
						const auto& v1 = data_->surface_nodes[face.vertices[1]];
						const auto& v2 = data_->surface_nodes[face.vertices[2]];
						// v0 · (v1 × v2)
						double cx = v1[1] * v2[2] - v1[2] * v2[1];
						double cy = v1[2] * v2[0] - v1[0] * v2[2];
						double cz = v1[0] * v2[1] - v1[1] * v2[0];
						signed_vol += v0[0] * cx + v0[1] * cy + v0[2] * cz;
					}
					signed_vol /= 6.0;
					bool normals_outward = (signed_vol > 0);

					// Winding number: sum of signed solid angles / 4π
					// If ≈ 1 → target inside, if ≈ 0 → target outside
					double winding = 0;
					for (std::uint32_t fi = 0; fi < nf; ++fi) {
						if (comp_id[fi] != cid) continue;
						const auto& face = data_->surface_faces[fi];
						Vec3 a = vec3_sub(data_->surface_nodes[face.vertices[0]], target_point_);
						Vec3 b = vec3_sub(data_->surface_nodes[face.vertices[1]], target_point_);
						Vec3 c = vec3_sub(data_->surface_nodes[face.vertices[2]], target_point_);
						double la = vec3_length(a), lb = vec3_length(b), lc = vec3_length(c);
						if (la < 1e-15 || lb < 1e-15 || lc < 1e-15) continue;
						Vec3 bxc = vec3_cross(b, c);
						double num = vec3_dot(a, bxc);
						double den = la * lb * lc + vec3_dot(a, b) * lc + vec3_dot(a, c) * lb + vec3_dot(b, c) * la;
						winding += 2.0 * std::atan2(num, den);
					}
					winding /= (4.0 * kPI);
					bool target_inside = (std::abs(winding) > 0.5);

					// Decide: target inside → normals should point inward (toward target)
					//         target outside → normals should point outward (away from component, toward target region)
					bool want_outward = !target_inside;
					bool need_flip = (normals_outward != want_outward);

					int comp_face_count = 0;
					if (need_flip) {
						for (std::uint32_t fi = 0; fi < nf; ++fi) {
							if (comp_id[fi] != cid) continue;
							auto& face = data_->surface_faces[fi];
							if (face.num_vertices == 3) std::swap(face.vertices[1], face.vertices[2]);
							else if (face.num_vertices == 4) std::swap(face.vertices[1], face.vertices[3]);
							flipped_count++;
							comp_face_count++;
						}
					}
					else {
						for (std::uint32_t fi = 0; fi < nf; ++fi)
							if (comp_id[fi] == cid) comp_face_count++;
					}

					bl_log_info(
						"BL component %d: %d faces, vol=%.1f (%s), winding=%.2f (target %s) -> %s",
						cid,
						comp_face_count,
						signed_vol,
						normals_outward ? "outward" : "inward",
						winding,
						target_inside ? "inside" : "outside",
						need_flip ? "FLIP" : "keep");
#if 0
					std::fprintf(stderr, "[BL] Component %d: %d faces, vol=%.1f (%s), winding=%.2f (target %s) → %s\n",
						cid, comp_face_count, signed_vol,
						normals_outward ? "outward" : "inward",
						winding, target_inside ? "inside" : "outside",
						need_flip ? "FLIP" : "keep");
#endif
				}

				if (flipped_count > 0) data_->build_connectivity();
				normals_flipped_ = false; // per-component, not global
			}
			else if (force_flip_set_) {
				// Global flip from region-based detection
				if (force_flip_) {
					for (std::uint32_t fi = 0; fi < nf; ++fi) {
						auto& face = data_->surface_faces[fi];
						if (face.num_vertices == 3) std::swap(face.vertices[1], face.vertices[2]);
						else if (face.num_vertices == 4) std::swap(face.vertices[1], face.vertices[3]);
					}
					data_->build_connectivity();
					normals_flipped_ = true;
				}
			}
			else {
				// Auto-detect: centroid method (fallback)
				const auto nn = data_->surface_nodes.size();
				Vec3 centroid = { 0, 0, 0 };
				for (std::uint32_t i = 0; i < nn; ++i) {
					centroid[0] += data_->surface_nodes[i][0];
					centroid[1] += data_->surface_nodes[i][1];
					centroid[2] += data_->surface_nodes[i][2];
				}
				centroid[0] /= nn; centroid[1] /= nn; centroid[2] /= nn;

				int outward = 0, inward = 0;
				const int step = std::max(1, static_cast<int>(nf) / 100);
				for (int i = 0; i < static_cast<int>(nf); i += step) {
					const auto& face = data_->surface_faces[i];
					auto n = calc_face_normal_from_three(
						data_->surface_nodes[face.vertices[0]],
						data_->surface_nodes[face.vertices[1]],
						data_->surface_nodes[face.vertices[2]]);
					Vec3 fc = {
					  (data_->surface_nodes[face.vertices[0]][0] + data_->surface_nodes[face.vertices[1]][0] + data_->surface_nodes[face.vertices[2]][0]) / 3.0,
					  (data_->surface_nodes[face.vertices[0]][1] + data_->surface_nodes[face.vertices[1]][1] + data_->surface_nodes[face.vertices[2]][1]) / 3.0,
					  (data_->surface_nodes[face.vertices[0]][2] + data_->surface_nodes[face.vertices[1]][2] + data_->surface_nodes[face.vertices[2]][2]) / 3.0
					};
					if (vec3_dot(n, vec3_sub(centroid, fc)) > 0) ++inward; else ++outward;
				}
				if (outward > inward) {
					for (std::uint32_t fi = 0; fi < nf; ++fi) {
						auto& face = data_->surface_faces[fi];
						if (face.num_vertices == 3) std::swap(face.vertices[1], face.vertices[2]);
						else if (face.num_vertices == 4) std::swap(face.vertices[1], face.vertices[3]);
					}
					data_->build_connectivity();
					normals_flipped_ = true;
				}
			}
		}

		// Build BVH on non-BL faces for collision detection
		if (!data_->non_bl_face_verts.empty()) {
			data_->non_bl_bvh.build(
				data_->non_bl_nodes, data_->non_bl_face_verts, data_->non_bl_face_nv);
		}

		// Initialize thickness
		initialize_thickness_of_nodes();

		current_layer_ = 0;

		for (int layer = 0; layer < params_.num_layers; ++layer) {
			current_layer_ = layer;

			auto layer_start = std::chrono::high_resolution_clock::now();

			// 1. Update face normals
			update_face_normals();
			// 2. Classify nodes and update boundary/feature angles
			classify_nodes_and_update_angles();

			// 3. Build the node march field for the current front
			build_node_march_field();

			// 4. Propagate the active front by one layer
			const bool continue_generation = propagate_front_one_layer(layer);

			auto layer_end = std::chrono::high_resolution_clock::now();
			double layer_ms = std::chrono::duration<double, std::milli>(layer_end - layer_start).count();

			// Count active nodes (nodes still marching forward)
			std::size_t active = 0;
			for (std::size_t i = 0; i < data_->node_active.size(); ++i) {
				if (data_->node_active[i]) ++active;
			}
			bl_log_info(
				"BL layer %d: %zu / %zu nodes still active (%.1f%%) in %.1f ms",
				layer + 1,
				active,
				data_->node_active.size(),
				100.0 * static_cast<double>(active) /
				std::max<std::size_t>(1, data_->node_active.size()),
				layer_ms);


			if (!continue_generation) {
				bl_log_info(
					"Boundary layer generation stopped after layer %d while retaining valid layers",
					layer + 1);
				break;
			}

			if (active == 0) {
				break;
			}
		}

		auto t_end = std::chrono::high_resolution_clock::now();
		double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
		bl_log_info(
			"Boundary layer generation complete: %zu layers in %.1f ms",
			data_->layer_positions.size() - 1,
			total_ms);
	}

	// ===========================================================================
	// Face normals
	// ===========================================================================

	void BoundaryLayerSolver::update_face_normals()
	{
		const auto nf = data_->surface_faces.size();
		data_->face_normal.resize(nf);

		const auto& positions = data_->layer_positions.back(); // Current front

		for (std::uint32_t fi = 0; fi < nf; ++fi) {
			if (!data_->face_active[fi]) {
				data_->face_normal[fi] = { 0, 0, 0 };
				continue;
			}

			const auto& face = data_->surface_faces[fi];
			if (face.num_vertices == 3) {
				data_->face_normal[fi] = calc_face_normal_from_three(
					positions[face.vertices[0]],
					positions[face.vertices[1]],
					positions[face.vertices[2]]);
			}
			else if (face.num_vertices == 4) {
				// Quad: average of 4 sub-triangle normals
				auto n0 = calc_face_normal_from_three(
					positions[face.vertices[0]], positions[face.vertices[1]], positions[face.vertices[2]]);
				auto n1 = calc_face_normal_from_three(
					positions[face.vertices[0]], positions[face.vertices[2]], positions[face.vertices[3]]);
				auto n2 = calc_face_normal_from_three(
					positions[face.vertices[1]], positions[face.vertices[2]], positions[face.vertices[3]]);
				auto n3 = calc_face_normal_from_three(
					positions[face.vertices[0]], positions[face.vertices[1]], positions[face.vertices[3]]);
				Vec3 avg = {
				  (n0[0] + n1[0] + n2[0] + n3[0]) * 0.25,
				  (n0[1] + n1[1] + n2[1] + n3[1]) * 0.25,
				  (n0[2] + n1[2] + n2[2] + n3[2]) * 0.25
				};
				data_->face_normal[fi] = vec3_normalized(avg);
			}
		}
	}

	// ===========================================================================
	// Node classification and angle updates
	// ===========================================================================

	void BoundaryLayerSolver::classify_nodes_and_update_angles()
	{
		const auto nn = data_->surface_nodes.size();

		// 1. Detect boundary nodes: if #VFConn != #VVConn, it's a boundary node
		if (data_->is_boundary_node.empty()) {
			data_->is_boundary_node.resize(nn, false);
			for (std::uint32_t node = 0; node < nn; ++node) {
				const auto nface = data_->node_to_faces[node].size();
				const auto nneigh = data_->node_to_nodes[node].size();
				data_->is_boundary_node[node] = (nface != nneigh);
			}
		}


		data_->node_dihedral_angles.resize(nn);
		data_->node_type.resize(nn);
		data_->node_manifold_angle.resize(nn, kPI);
		const auto& pos = data_->layer_positions.back();

		constexpr double kGroupAngleCos = 0.9; // cos(~25°) grouping threshold

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) {
				data_->node_type[node] = NodeTypeBL::OTHER;
				continue;
			}

			auto& angles = data_->node_dihedral_angles[node];
			angles.clear();

			const auto& node_faces = data_->node_to_faces[node];

			if (data_->is_boundary_node[node]) {
				data_->node_type[node] = NodeTypeBL::FLAT;
				data_->node_manifold_angle[node] = kPI;
				continue;
			}

			// -- Step A: Group faces by normal similarity --
			struct NormalGroup {
				Vec3 avg_normal;
				std::vector<std::uint32_t> face_ids;
			};
			std::vector<NormalGroup> groups;

			for (auto fi : node_faces) {
				if (!data_->face_active[fi]) continue;
				const auto& fn = data_->face_normal[fi];

				bool found = false;
				for (auto& g : groups) {
					if (vec3_dot(g.avg_normal, fn) > kGroupAngleCos) {
						// Merge into this group, update averaged normal
						auto scaled = vec3_scale(g.avg_normal,
							static_cast<double>(g.face_ids.size()));
						g.avg_normal = vec3_normalized(vec3_add(scaled, fn));
						g.face_ids.push_back(fi);
						found = true;
						break;
					}
				}
				if (!found) {
					groups.push_back({ fn, {fi} });
				}
			}

			// Only one normal group -> flat node
			if (groups.size() <= 1) {
				data_->node_type[node] = NodeTypeBL::FLAT;
				data_->node_manifold_angle[node] = kPI;
				continue;
			}

			// -- Step B: Compute dihedral angles between adjacent groups --
			for (std::size_t gi = 0; gi < groups.size(); ++gi) {
				for (std::size_t gj = gi + 1; gj < groups.size(); ++gj) {
					// Find a pair of faces (one from each group) sharing an edge
					std::uint32_t fa_id = 0, fb_id = 0, shared_v = 0;
					bool adjacent = false;

					for (auto fa : groups[gi].face_ids) {
						if (adjacent) break;
						for (auto fb : groups[gj].face_ids) {
							if (adjacent) break;
							const auto& faceA = data_->surface_faces[fa];
							const auto& faceB = data_->surface_faces[fb];
							for (int vi = 0; vi < faceA.num_vertices; ++vi) {
								auto v = faceA.vertices[vi];
								if (v == node) continue;
								for (int vj = 0; vj < faceB.num_vertices; ++vj) {
									if (faceB.vertices[vj] == v) {
										fa_id = fa; fb_id = fb; shared_v = v;
										adjacent = true;
										goto found_adj;
									}
								}
							}
						}
					}
				found_adj:
					if (!adjacent) continue;

					// Edge direction follows face A's winding
					const auto& faceA = data_->surface_faces[fa_id];
					bool node_then_v = false;
					for (int k = 0; k < faceA.num_vertices; ++k) {
						if (faceA.vertices[k] == node &&
							faceA.vertices[(k + 1) % faceA.num_vertices] == shared_v) {
							node_then_v = true;
							break;
						}
					}
					Vec3 edge_dir;
					if (node_then_v) {
						edge_dir = vec3_normalized(vec3_sub(pos[shared_v], pos[node]));
					}
					else {
						edge_dir = vec3_normalized(vec3_sub(pos[node], pos[shared_v]));
					}

					// Compute dihedral using GROUP normals (not individual face normals)
					double beta = get_face_dihedral_angle(
						groups[gi].avg_normal,
						groups[gj].avg_normal,
						edge_dir);
					if (normals_flipped_) {
						beta = 2.0 * kPI - beta;
					}
					angles.push_back(beta);
				}
			}

			if (angles.empty()) {
				angles.push_back(kPI);
			}


			double max_angle = 0.0, min_angle = 2.0 * kPI;
			for (double a : angles) {
				max_angle = std::max(max_angle, a);
				min_angle = std::min(min_angle, a);
			}

			bool has_convex = (max_angle > BL3d_CONVEX_ANGLE);
			bool has_concave = (min_angle < BL3d_CONCAVE_ANGLE);

			if (has_convex && has_concave) {
				data_->node_type[node] = NodeTypeBL::CONCAVE_CONVEX;
				double sum = 0;
				for (double a : angles) sum += a;
				data_->node_manifold_angle[node] = sum / static_cast<double>(angles.size());
			}
			else if (has_convex) {
				data_->node_type[node] = NodeTypeBL::CONVEX;
				data_->node_manifold_angle[node] = max_angle;
			}
			else if (has_concave) {
				data_->node_type[node] = NodeTypeBL::CONCAVE;
				data_->node_manifold_angle[node] = min_angle;
			}
			else {
				data_->node_type[node] = NodeTypeBL::FLAT;
				data_->node_manifold_angle[node] = kPI;
			}
		}

		// Debug: print node classification stats
	}

	// ===========================================================================
	// Node march field construction
	// ===========================================================================

	Vec3 BoundaryLayerSolver::compute_face_normal_tri(std::uint32_t fi) const
	{
		const auto& f = data_->surface_faces[fi];
		const auto& pos = data_->layer_positions.back();
		return calc_face_normal_from_three(
			pos[f.vertices[0]], pos[f.vertices[1]], pos[f.vertices[2]]);
	}

	Vec3 BoundaryLayerSolver::compute_face_normal_quad(std::uint32_t fi) const
	{
		return data_->face_normal[fi]; // Already computed in update_face_normals
	}

	Vec3 BoundaryLayerSolver::compute_average_node_normal(std::uint32_t node) const
	{
		const auto& faces = data_->node_to_faces[node];
		if (faces.empty()) {
			return { 0, 0, 0 };
		}

		std::vector<Vec3> active_normals;
		active_normals.reserve(faces.size());
		for (auto fi : faces) {
			if (!data_->face_active[fi]) {
				continue;
			}
			active_normals.push_back(data_->face_normal[fi]);
		}

		if (active_normals.empty()) {
			return { 0, 0, 0 };
		}
		if (active_normals.size() == 1U) {
			return active_normals.front();
		}

		struct NormalGroup {
			Vec3 sum = { 0.0, 0.0, 0.0 };
			std::size_t count = 0;
		};

		const double kGroupAngleCos = std::cos(5.0 * kPI / 180.0);
		std::vector<NormalGroup> groups;
		groups.reserve(active_normals.size());

		for (const auto& fn : active_normals) {
			bool assigned = false;
			for (auto& group : groups) {
				const Vec3 group_dir = vec3_normalized(group.sum);
				if (vec3_dot(group_dir, fn) < kGroupAngleCos) {
					continue;
				}
				group.sum = vec3_add(group.sum, fn);
				++group.count;
				assigned = true;
				break;
			}

			if (!assigned) {
				groups.push_back({ fn, 1U });
			}
		}

		Vec3 normal = { 0.0, 0.0, 0.0 };
		for (const auto& group : groups) {
			if (group.count == 0U) {
				continue;
			}
			const Vec3 group_average = vec3_scale(group.sum, 1.0 / static_cast<double>(group.count));
			normal = vec3_add(normal, vec3_normalized(group_average));
		}

		if (vec3_length(normal) < kEPS) {
			return { 0, 0, 0 };
		}
		return vec3_normalized(normal);
	}

	Vec3 BoundaryLayerSolver::iterate_by_angle(std::uint32_t node) const
	{
		const auto& faces = data_->node_to_faces[node];

		// Collect active face normals
		std::vector<Vec3> normals;
		for (auto fi : faces) {
			if (!data_->face_active[fi]) continue;
			normals.push_back(data_->face_normal[fi]);
		}

		if (normals.empty()) return { 0, 0, 0 };
		if (normals.size() == 1) return normals[0];

		// Start with area-weighted average
		Vec3 best = compute_average_node_normal(node);

		constexpr int max_iter = 20;
		constexpr double convergence = 0.00001;

		const auto compute_min_angle = [&](const Vec3& direction) {
			double min_angle = kPI;
			for (const auto& normal : normals) {
				min_angle = std::min(min_angle, get_vec3d_angle(direction, normal));
			}
			return min_angle;
			};

		Vec3 candidate = best;
		for (int iter = 0; iter < max_iter; ++iter) {
			// Find the normal with the smallest angle to candidate
			double min_angle = kPI;
			int min_idx = -1;
			for (int i = 0; i < static_cast<int>(normals.size()); ++i) {
				double a = get_vec3d_angle(candidate, normals[i]);
				if (a < min_angle) {
					min_angle = a;
					min_idx = i;
				}
			}

			if (min_idx < 0) break;

			Vec3 adjustment = normals[min_idx];
			bool accepted = false;
			double delta = 0.0;

			for (const double step : {0.1, 0.01}) {
				Vec3 new_dir = vec3_add(candidate, vec3_scale(adjustment, step));
				new_dir = vec3_normalized(new_dir);

				const double new_min_angle = compute_min_angle(new_dir);
				if (new_min_angle < min_angle) {
					continue;
				}

				delta = get_vec3d_angle(new_dir, candidate);
				candidate = new_dir;
				accepted = true;
				break;
			}

			if (!accepted) {
				break;
			}

			if (delta < convergence) break;
		}

		return candidate;
	}

	void BoundaryLayerSolver::build_node_march_field()
	{
		const auto nn = data_->surface_nodes.size();
		const auto& pos = data_->layer_positions.back();
		const bool strict_aspect_height = params_.use_aspect_first_height();

		// ==================================================================
		// 1. Compute initial normals
		// ==================================================================
		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) continue;

			switch (data_->node_type[node]) {
			case NodeTypeBL::FLAT:
				initial_normals_[node] = compute_average_node_normal(node);
				break;
			case NodeTypeBL::CONCAVE:
			case NodeTypeBL::CONVEX:
			case NodeTypeBL::CONCAVE_CONVEX:
				initial_normals_[node] = iterate_by_angle(node);
				break;
			default:
				initial_normals_[node] = compute_average_node_normal(node);
				break;
			}
		}
		final_normals_ = initial_normals_;

		// ==================================================================
		// 2. Compute initial distances.
		// Keep the historical alpha-amplification path in comments for now.
		// ==================================================================
		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) continue;
			initial_distances_[node] = layer_height_for_node(node, current_layer_);
		}
		final_distances_ = initial_distances_;

		// auto alphas = compute_distance_coefficient();
		//
		// for(std::uint32_t node = 0; node < nn; ++node) {
		//   if(!data_->node_active[node]) continue;
		//   initial_distances_[node] *= alphas[node];
		// }
		// final_distances_ = initial_distances_;


		constexpr double omega_dis = 0.5;
		constexpr int dist_smooth_iters = 5;
		for (int iter = 0; iter < dist_smooth_iters; ++iter) {
			auto dist_copy = final_distances_;
			for (std::uint32_t node = 0; node < nn; ++node) {
				if (!data_->node_active[node]) continue;
				const auto& neighbors = data_->node_to_nodes[node];
				if (neighbors.empty()) continue;

				double sum_w = 0, sum_wd = 0;
				for (auto ai : neighbors) {
					if (!data_->node_active[ai]) continue;
					double d = vec3_length(vec3_sub(pos[ai], pos[node]));
					double w = (d > kEPS) ? 1.0 / d : 1.0;
					sum_w += w;
					sum_wd += w * dist_copy[ai];
				}
				if (sum_w > kEPS) {
					final_distances_[node] = (1.0 - omega_dis) * dist_copy[node]
						+ omega_dis * (sum_wd / sum_w);
				}
			}
		}

		if (params_.smooth_normals && current_layer_ >= params_.orthogonal_layers) {

			std::vector<double> node_omega(nn, 0.5);
			for (std::uint32_t node = 0; node < nn; ++node) {
				if (!data_->node_active[node]) continue;
				double ma = data_->node_manifold_angle[node];
				// Nodes at sharp features get lower omega (less smoothing)
				// Flat nodes (ma≈π) get higher omega (more smoothing)
				double ratio = std::abs(ma - kPI) / kPI; // 0=flat, ~1=sharp
				node_omega[node] = std::clamp(0.6 - 0.4 * ratio, 0.1, 0.6);
			}

			// 4b. Compute inverse-distance weights for each node's neighbors
			std::vector<double> link_over_dist_sum(nn, 0.0);
			std::vector<std::vector<double>> link_over_dist(nn);
			for (std::uint32_t node = 0; node < nn; ++node) {
				if (!data_->node_active[node]) continue;
				const auto& neighbors = data_->node_to_nodes[node];
				link_over_dist[node].resize(neighbors.size(), 0.0);
				for (std::size_t i = 0; i < neighbors.size(); ++i) {
					auto ai = neighbors[i];
					if (!data_->node_active[ai]) continue;
					double d = vec3_length(vec3_sub(pos[ai], pos[node]));
					if (d > kEPS) {
						link_over_dist[node][i] = 1.0 / d;
						link_over_dist_sum[node] += 1.0 / d;
					}
				}
			}

			// 4c. Iterative normal optimization loop (up to 500 iterations)
			//     Optimization: skip converged nodes, inline vector math
			constexpr int kMaxNormalIter = 500;
			constexpr double kConvergenceThreshold = 1e-7;
			const double cos_max_dev = std::cos(params_.max_normal_deviation_deg * kPI / 180.0);

			for (int iter = 1; iter <= kMaxNormalIter; ++iter) {
				double entity_group_max_change = 0.0;

				for (std::uint32_t node = 0; node < nn; ++node) {
					if (!data_->node_active[node]) continue;
					const double inv_w = link_over_dist_sum[node];
					if (inv_w <= 1e-10) continue;

					const auto& neighbors = data_->node_to_nodes[node];
					const auto& weights = link_over_dist[node];

					// Inline weighted average of neighbor normals
					double sx = 0, sy = 0, sz = 0;
					for (std::size_t i = 0; i < neighbors.size(); ++i) {
						const double w = weights[i];
						if (w <= 0) continue;
						const auto& fn = final_normals_[neighbors[i]];
						sx += fn[0] * w;
						sy += fn[1] * w;
						sz += fn[2] * w;
					}
					const double inv = 1.0 / inv_w;
					sx *= inv; sy *= inv; sz *= inv;
					double len = std::sqrt(sx * sx + sy * sy + sz * sz);
					if (len < kEPS) continue;
					sx /= len; sy /= len; sz /= len;

					// Blend: omega * averaged + (1-omega) * current
					const double omega = node_omega[node];
					const double om1 = 1.0 - omega;
					const auto& cur = final_normals_[node];
					double nx = omega * sx + om1 * cur[0];
					double ny = omega * sy + om1 * cur[1];
					double nz = omega * sz + om1 * cur[2];
					len = std::sqrt(nx * nx + ny * ny + nz * nz);
					if (len < kEPS) continue;
					nx /= len; ny /= len; nz /= len;

					// Check angle against all adjacent face normals
					bool rejected = false;
					for (auto fi : data_->node_to_faces[node]) {
						if (!data_->face_active[fi]) continue;
						const auto& fn = data_->face_normal[fi];
						if (nx * fn[0] + ny * fn[1] + nz * fn[2] < 0.0) { rejected = true; break; }
					}
					if (rejected) continue; // Skip this iter, retry next time

					// Check deviation from initial normal
					const auto& ini = initial_normals_[node];
					if (nx * ini[0] + ny * ini[1] + nz * ini[2] < cos_max_dev) continue;

					// Accept: blend new with current
					double ax = nx + cur[0], ay = ny + cur[1], az = nz + cur[2];
					len = std::sqrt(ax * ax + ay * ay + az * az);
					if (len < kEPS) continue;
					ax /= len; ay /= len; az /= len;

					double dx = ax - cur[0], dy = ay - cur[1], dz = az - cur[2];
					double change = std::sqrt(dx * dx + dy * dy + dz * dz);
					entity_group_max_change = std::max(entity_group_max_change, change);
					final_normals_[node] = { ax, ay, az };
				}

				if (entity_group_max_change < kConvergenceThreshold) {
					bl_log_info("BL normal smoothing converged at iter %d", iter);
					break;
				}
			}
		}


		// 5a. Compute node_smooth_height: project march vector onto face normals
		std::vector<double> node_smooth_height(nn, 0.0);
		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) continue;
			Vec3 march_vec = vec3_scale(final_normals_[node], final_distances_[node]);
			const auto& faces = data_->node_to_faces[node];
			double proj_sum = 0.0;
			int nf = 0;
			for (auto fi : faces) {
				if (!data_->face_active[fi]) continue;
				double proj = std::abs(vec3_dot(march_vec, data_->face_normal[fi]));
				proj_sum += proj;
				nf++;
			}
			node_smooth_height[node] = (nf > 0) ? proj_sum / nf : final_distances_[node];
		}

		// 5b. Iterative smoothing loop
		std::vector<double> smooth_offset(nn, 0.0);
		for (std::uint32_t node = 0; node < nn; ++node) {
			smooth_offset[node] = node_smooth_height[node];
		}

		constexpr int kMaxHeightIter = 100;
		for (int iter = 0; iter < kMaxHeightIter; ++iter) {
			bool break_smooth = true;
			for (std::uint32_t node = 0; node < nn; ++node) {
				if (!data_->node_active[node]) continue;
				const auto& neighbors = data_->node_to_nodes[node];

				double sum_w = 0, sum_wh = 0;
				for (auto ai : neighbors) {
					if (!data_->node_active[ai]) continue;
					double d = vec3_length(vec3_sub(pos[ai], pos[node]));
					double w = (d > kEPS) ? 1.0 / d : 1.0;
					sum_w += w;
					sum_wh += w * smooth_offset[ai];
				}
				if (sum_w <= kEPS) continue;
				double smooth_avg = sum_wh / sum_w;

				// Omega based on manifold angle (sharper features -> less smoothing)
				double ma = data_->node_manifold_angle[node];
				double ratio = std::abs(ma - kPI) / kPI;
				double omega = std::clamp(0.125 * std::exp(-ratio) * 0.636, 0.01, 0.5);

				double height = node_smooth_height[node] * omega
					+ smooth_avg * (1.0 - omega);

				if (std::abs(smooth_offset[node] - height) >= 0.01 * std::abs(smooth_offset[node])) {
					break_smooth = false;
				}
				smooth_offset[node] = height;
			}
			if (break_smooth) break;
		}

		// 5c. Apply smoothed heights back to final_distances
		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) continue;
			if (node_smooth_height[node] > 0.05) {
				final_distances_[node] = smooth_offset[node] / node_smooth_height[node]
					* final_distances_[node];
			}
		}


		// ==================================================================
		// 6. Squeeze by shortest distance (collision detection)
		// ==================================================================
		squeeze_thickness_by_shortest_distance();

		// ==================================================================
		// 7. Build final march vector
		// ==================================================================
		data_->march_direction.resize(nn);
		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) {
				data_->march_direction[node] = { 0, 0, 0 };
				continue;
			}
			data_->march_direction[node] = vec3_scale(
				final_normals_[node], final_distances_[node]);
		}
	}

	void BoundaryLayerSolver::smooth_normals_laplacian(std::vector<Vec3>& normals)
	{
		const auto nn = data_->surface_nodes.size();
		constexpr double omega = 0.5;

		for (int iter = 0; iter < params_.smooth_iterations; ++iter) {
			std::vector<Vec3> smoothed(nn);

			for (std::uint32_t node = 0; node < nn; ++node) {
				if (!data_->node_active[node] || data_->is_boundary_node[node]) {
					smoothed[node] = normals[node];
					continue;
				}

				// Weighted average of neighbors
				Vec3 sum = { 0, 0, 0 };
				double weight_sum = 0;

				for (auto neigh : data_->node_to_nodes[node]) {
					if (!data_->node_active[neigh]) continue;
					double w = 1.0;
					// Weight by node type: reduce smoothing at features
					if (data_->node_type[neigh] == NodeTypeBL::CONCAVE ||
						data_->node_type[neigh] == NodeTypeBL::CONVEX) {
						w = 0.3;
					}
					sum = vec3_add(sum, vec3_scale(normals[neigh], w));
					weight_sum += w;
				}

				if (weight_sum < kEPS) {
					smoothed[node] = normals[node];
					continue;
				}

				Vec3 avg = vec3_scale(sum, 1.0 / weight_sum);
				avg = vec3_normalized(avg);

				// Blend: omega * smoothed + (1-omega) * original
				Vec3 blended = vec3_add(
					vec3_scale(avg, omega),
					vec3_scale(normals[node], 1.0 - omega));
				blended = vec3_normalized(blended);

				// Check max deviation from initial normal
				double angle_dev = get_vec3d_angle(blended, initial_normals_[node]);
				double max_dev = params_.max_normal_deviation_deg * kPI / 180.0;

				if (angle_dev > max_dev) {
					// Clamp: interpolate back toward initial
					double t = max_dev / angle_dev;
					blended = vec3_add(
						vec3_scale(initial_normals_[node], 1.0 - t),
						vec3_scale(blended, t));
					blended = vec3_normalized(blended);
				}

				smoothed[node] = blended;
			}

			normals = smoothed;
		}
	}

	// ===========================================================================
	// Thickness computation
	// ===========================================================================

	void BoundaryLayerSolver::initialize_base_first_heights()
	{
		const auto nn = data_->surface_nodes.size();
		base_first_heights_.assign(nn, params_.first_height);

		if (!params_.use_aspect_first_height()) {
			return;
		}

		const auto& surface_nodes = data_->surface_nodes;
		double min_height = std::numeric_limits<double>::max();
		double max_height = 0.0;
		double sum_height = 0.0;
		std::size_t num_active = 0;

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) {
				base_first_heights_[node] = 0.0;
				continue;
			}

			const auto& neighbors = data_->node_to_nodes[node];
			double base_height = params_.first_height;
			double sum_length = 0.0;
			int valid_neighbor_count = 0;

			for (const auto neighbor : neighbors) {
				const double edge_length = vec3_length(
					vec3_sub(surface_nodes[neighbor], surface_nodes[node]));
				if (edge_length > kEPS && std::isfinite(edge_length)) {
					sum_length += edge_length;
					++valid_neighbor_count;
				}
			}

			if (valid_neighbor_count > 0) {
				const double average_length =
					sum_length / static_cast<double>(valid_neighbor_count);
				base_height = params_.first_height_aspect * average_length;
			}

			if (!std::isfinite(base_height) || base_height <= 0.0) {
				base_height = params_.first_height;
			}
			if (!params_.use_aspect_first_height() && params_.min_first_height > 0.0) {
				base_height = std::max(base_height, params_.min_first_height);
			}

			base_first_heights_[node] = base_height;
			min_height = std::min(min_height, base_height);
			max_height = std::max(max_height, base_height);
			sum_height += base_height;
			++num_active;
		}

		if (num_active > 0) {
			bl_log_info(
				"BL aspect first height: aspect=%.6g, min=%.6g, max=%.6g, avg=%.6g",
				params_.first_height_aspect,
				min_height,
				max_height,
				sum_height / static_cast<double>(num_active));
		}
	}

	double BoundaryLayerSolver::layer_height_for_node(
		std::uint32_t node,
		int layer_index
	) const noexcept
	{
		if (node >= base_first_heights_.size()) {
			return params_.layer_height(layer_index);
		}
		return params_.layer_height_from_base(base_first_heights_[node], layer_index);
	}

	void BoundaryLayerSolver::initialize_thickness_of_nodes()
	{
		const auto nn = data_->surface_nodes.size();
		initialize_base_first_heights();
		initial_distances_.resize(nn);
		final_distances_.resize(nn);

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) {
				initial_distances_[node] = 0.0;
				continue;
			}
			initial_distances_[node] = layer_height_for_node(node, 0);
		}

		final_distances_ = initial_distances_;
	}

	void BoundaryLayerSolver::smooth_thickness()
	{
		const auto nn = data_->surface_nodes.size();
		constexpr double omega = 0.5;
		constexpr int smooth_iters = 3;

		std::vector<double> smoothed(nn);

		for (int iter = 0; iter < smooth_iters; ++iter) {
			for (std::uint32_t node = 0; node < nn; ++node) {
				if (!data_->node_active[node]) {
					smoothed[node] = final_distances_[node];
					continue;
				}

				double sum = 0;
				int count = 0;
				for (auto neigh : data_->node_to_nodes[node]) {
					if (!data_->node_active[neigh]) continue;
					sum += final_distances_[neigh];
					++count;
				}

				if (count == 0) {
					smoothed[node] = final_distances_[node];
					continue;
				}

				double avg = sum / count;
				smoothed[node] = (1.0 - omega) * final_distances_[node] + omega * avg;
			}

			final_distances_ = smoothed;
		}
	}

	void BoundaryLayerSolver::squeeze_thickness_by_shortest_distance()
	{
		if (!params_.allow_squeeze) return;

		const auto nn = data_->surface_nodes.size();
		const auto nf = data_->surface_faces.size();
		const auto& current_pos = data_->layer_positions.back();
		FaceBVH front_bvh;
		std::vector<std::array<std::uint32_t, 4>> front_face_verts;
		std::vector<int> front_face_nv;
		std::vector<std::uint32_t> front_local_to_global;
		std::vector<int> front_global_to_local;
		build_active_face_bvh(
			current_pos,
			front_bvh,
			front_face_verts,
			front_face_nv,
			front_local_to_global,
			front_global_to_local
		);

		std::vector<bool> skip_front(front_face_verts.size(), false);
		std::size_t wall_limited = 0;
		std::size_t front_limited = 0;

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) continue;

			const auto& dir = final_normals_[node];
			double dist = final_distances_[node];
			if (dist < kEPS) continue;

			double best_allowed = dist;
			bool limited_by_wall = false;
			bool limited_by_front = false;

			// Check collision with non-BL faces using BVH
			if (!data_->non_bl_face_verts.empty()) {
				std::vector<bool> skip_face(data_->non_bl_face_verts.size(), false);
				double hit = data_->non_bl_bvh.ray_nearest_intersection(
					current_pos[node], dir, dist * 3.0, skip_face);

				if (hit > 0) {
					double max_allowed = hit * params_.proximity_factor;
					if (max_allowed < best_allowed) {
						best_allowed = max_allowed;
						limited_by_wall = true;
					}
				}
			}

			// Check collision with other BL front faces as well. This catches
			// opposite walls and disconnected components that the non-BL BVH
			// cannot see.
			if (!front_face_verts.empty()) {
				std::fill(skip_front.begin(), skip_front.end(), false);
				for (auto fi : data_->node_to_faces[node]) {
					if (fi >= static_cast<std::uint32_t>(nf)) {
						continue;
					}

					const int local_fi = front_global_to_local[fi];
					if (local_fi >= 0) {
						skip_front[static_cast<std::size_t>(local_fi)] = true;
					}

					for (auto nfi : data_->face_to_faces[fi]) {
						const int local_nfi = front_global_to_local[nfi];
						if (local_nfi >= 0) {
							skip_front[static_cast<std::size_t>(local_nfi)] = true;
						}
					}
				}

				const double origin_bias = std::max(dist * 1.0e-6, 1.0e-9);
				const Vec3 query_origin = vec3_add(current_pos[node], vec3_scale(dir, origin_bias));
				double hit = front_bvh.ray_nearest_intersection(
					query_origin, dir, dist * 3.0 + origin_bias, skip_front);

				if (hit > 0) {
					hit = std::max(hit - origin_bias, 0.0);
					double max_allowed = hit * params_.proximity_factor;
					if (max_allowed < best_allowed) {
						best_allowed = max_allowed;
						limited_by_front = true;
					}
				}
			}

			if (best_allowed < dist) {
				const double min_allowed = params_.use_aspect_first_height()
					? 0.0
					: params_.min_first_height;
				final_distances_[node] = std::min(dist, std::max(best_allowed, min_allowed));
				if (limited_by_wall) {
					++wall_limited;
				}
				if (limited_by_front) {
					++front_limited;
				}
			}
		}

	}

	// ===========================================================================
	// Front propagation
	// ===========================================================================

	bool BoundaryLayerSolver::propagate_front_one_layer(int layer)
	{
		const auto nn = data_->surface_nodes.size();
		const auto current_pos = data_->layer_positions.back();
		clear_problem_marks();

		deactivate_reached_front_entities(layer);

		std::unordered_set<std::uint32_t> modify_nodes;
		std::unordered_set<std::uint32_t> modify_faces;

		std::vector<Vec3> next_positions(nn, current_pos.empty() ? Vec3{ 0.0, 0.0, 0.0 } : current_pos.front());
		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) {
				next_positions[node] = current_pos[node];
				continue;
			}
			next_positions[node] = vec3_add(current_pos[node], data_->march_direction[node]);
		}

		bool valid = check_node_corridor_intersections(
			current_pos,
			next_positions,
			modify_nodes,
			modify_faces
		);
		if (!valid && params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
			return false;
		}

		data_->layer_positions.push_back(next_positions);
		auto& top_positions = data_->layer_positions.back();

		valid = check_top_face_quality(
			current_pos,
			top_positions,
			modify_nodes,
			modify_faces
		) && valid;
		if (!valid && params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
			return false;
		}

		valid = check_top_face_intersections(
			top_positions,
			modify_nodes,
			modify_faces
		) && valid;
		valid = check_side_face_intersections(
			current_pos,
			top_positions,
			modify_nodes,
			modify_faces
		) && valid;
		if (!valid && params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
			return false;
		}

		if (!modify_nodes.empty() || !modify_faces.empty()) {
			apply_problematic_area_treatment(layer, modify_nodes, modify_faces);
			enforce_max_consecutive_collapsed_layers(modify_faces);
		}

		recompute_node_thickness();
		clear_problem_marks();
		return true;
	}

	void BoundaryLayerSolver::deactivate_reached_front_entities(int layer)
	{
		const auto nf = data_->surface_faces.size();
		const auto nn = data_->surface_nodes.size();

		for (std::uint32_t fi = 0; fi < nf; ++fi) {
			if (data_->max_face_level[fi] != layer) {
				continue;
			}
			data_->face_active[fi] = false;
		}

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (data_->max_node_level[node] != layer) {
				continue;
			}
			data_->node_active[node] = false;
		}
	}

	bool BoundaryLayerSolver::check_node_corridor_intersections(
		const std::vector<Vec3>& current_positions,
		const std::vector<Vec3>& next_positions,
		std::unordered_set<std::uint32_t>& modify_nodes,
		std::unordered_set<std::uint32_t>& modify_faces
	) const
	{
		const auto nn = data_->surface_nodes.size();
		const auto nf = data_->surface_faces.size();

		FaceBVH front_bvh;
		std::vector<std::array<std::uint32_t, 4>> face_verts;
		std::vector<int> face_nv;
		std::vector<std::uint32_t> local_to_global;
		std::vector<int> global_to_local;
		build_active_face_bvh(
			current_positions,
			front_bvh,
			face_verts,
			face_nv,
			local_to_global,
			global_to_local
		);

		std::vector<std::uint32_t> candidates;
		bool all_ok = true;

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node] || data_->is_boundary_node[node]) {
				continue;
			}

			const Vec3 march = vec3_sub(next_positions[node], current_positions[node]);
			const double node_height = vec3_length(march);
			if (node_height < kEPS) {
				continue;
			}

			Vec3 node_normal = vec3_scale(march, 1.0 / node_height);
			double max_neighbor_dist = 0.0;
			for (auto neigh : data_->node_to_nodes[node]) {
				max_neighbor_dist = std::max(
					max_neighbor_dist,
					vec3_length(vec3_sub(current_positions[neigh], current_positions[node]))
				);
			}

			double dist = 2.0 * node_height + 0.5 * max_neighbor_dist;
			double range = 0.5 * dist;
			const double visibility_angle = std::min(compute_node_visibility_angle(node), kPI * 0.25);
			const double cos_vis = std::max(std::cos(visibility_angle), 1.0e-8);
			range /= cos_vis;
			const Vec3 center = vec3_add(current_positions[node], vec3_scale(node_normal, range));
			const double radius = range * std::sin(visibility_angle);

			AABB query_box;
			query_box.expand(current_positions[node]);
			query_box.expand(next_positions[node]);
			query_box.expand(center);
			query_box.pad(std::max(radius, kIntersectionEps));

			if (!face_verts.empty()) {
				front_bvh.query_overlap(query_box, candidates);
				for (auto local_face : candidates) {
					const auto face_id = local_to_global[local_face];
					if (face_id >= nf || !data_->face_active[face_id]) {
						continue;
					}

					bool linked = false;
					for (auto fi : data_->node_to_faces[node]) {
						if (fi == face_id) {
							linked = true;
							break;
						}
					}
					if (linked) {
						continue;
					}

					const auto& face = data_->surface_faces[face_id];
					bool hits = false;
					for (int tri_idx = 0; tri_idx < face_triangle_count(face); ++tri_idx) {
						const auto tri = face_triangle(current_positions, face, tri_idx);
						if (segment_triangle_intersect_strict(current_positions[node], next_positions[node], tri)) {
							hits = true;
							break;
						}
					}
					if (!hits) {
						continue;
					}

					std::ostringstream details;
					details << std::fixed << std::setprecision(10);
					details << "### Trigger\n\n";
					details << "- check: `check_node_corridor_intersections()`\n";
					details << "- predicate: `segment_triangle_intersect_strict(current_positions[node], next_positions[node], tri)`\n";
					details << "- node_id: `" << node << "`\n";
					details << "- hit_face_id: `" << face_id << "`\n";
					details << "- node_height: `" << node_height << "`\n";
					details << "- visibility_angle_rad: `" << visibility_angle << "`\n";
					details << "- visibility_angle_deg: `" << (visibility_angle * 180.0 / kPI) << "`\n";
					details << "- query_distance: `" << dist << "`\n";
					details << "\n### Segment\n\n";
					details << "- current_position: `" << format_vec3(current_positions[node]) << "`\n";
					details << "- next_position: `" << format_vec3(next_positions[node]) << "`\n";
					details << "- march_direction: `" << format_vec3(node_normal) << "`\n";
					export_first_problem_debug(
						"node_corridor_intersection",
						static_cast<std::int32_t>(face_id),
						{ copy_face_points(current_positions, face) },
						{ face.num_vertices },
						details.str(),
						{ current_positions[node], next_positions[node] }
					);

					all_ok = false;
					modify_nodes.insert(node);
					modify_faces.insert(face_id);
					if (params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
						if (!params_.retain_valid_layers) {
							throw std::runtime_error("Boundary layer front intersects existing faces.");
						}
						return false;
					}
				}
			}

			if (!data_->non_bl_face_verts.empty()) {
				data_->non_bl_bvh.query_overlap(query_box, candidates);
				for (auto wall_local : candidates) {
					SurfaceFace wall_face;
					wall_face.vertices = data_->non_bl_face_verts[wall_local];
					wall_face.num_vertices = data_->non_bl_face_nv[wall_local];

					bool hits = false;
					for (int tri_idx = 0; tri_idx < face_triangle_count(wall_face); ++tri_idx) {
						const auto tri = face_triangle(data_->non_bl_nodes, wall_face, tri_idx);
						if (segment_triangle_intersect_strict(current_positions[node], next_positions[node], tri)) {
							hits = true;
							break;
						}
					}
					if (!hits) {
						continue;
					}

					std::ostringstream details;
					details << std::fixed << std::setprecision(10);
					details << "### Trigger\n\n";
					details << "- check: `check_node_corridor_intersections()`\n";
					details << "- predicate: `segment_triangle_intersect_strict(current_positions[node], next_positions[node], wall_face_triangle)`\n";
					details << "- node_id: `" << node << "`\n";
					details << "- wall_face_local_id: `" << wall_local << "`\n";
					details << "- node_height: `" << node_height << "`\n";
					details << "- visibility_angle_rad: `" << visibility_angle << "`\n";
					details << "- visibility_angle_deg: `" << (visibility_angle * 180.0 / kPI) << "`\n";
					details << "- query_distance: `" << dist << "`\n";
					details << "\n### Segment\n\n";
					details << "- current_position: `" << format_vec3(current_positions[node]) << "`\n";
					details << "- next_position: `" << format_vec3(next_positions[node]) << "`\n";
					details << "- march_direction: `" << format_vec3(node_normal) << "`\n";
					export_first_problem_debug(
						"node_corridor_wall_intersection",
						-1,
						{ copy_face_points(data_->non_bl_nodes, wall_face) },
						{ wall_face.num_vertices },
						details.str(),
						{ current_positions[node], next_positions[node] }
					);

					all_ok = false;
					modify_nodes.insert(node);
					if (params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
						if (!params_.retain_valid_layers) {
							throw std::runtime_error("Boundary layer front intersects wall faces.");
						}
						return false;
					}
				}
			}
		}

		return all_ok;
	}

	bool BoundaryLayerSolver::check_top_face_intersections(
		const std::vector<Vec3>& top_positions,
		std::unordered_set<std::uint32_t>& modify_nodes,
		std::unordered_set<std::uint32_t>& modify_faces
	) const
	{
		const auto nf = data_->surface_faces.size();
		bool all_ok = true;

		FaceBVH top_bvh;
		std::vector<std::array<std::uint32_t, 4>> top_face_verts;
		std::vector<int> top_face_nv;
		std::vector<std::uint32_t> top_local_to_global;
		std::vector<int> top_global_to_local;
		build_active_face_bvh(
			top_positions,
			top_bvh,
			top_face_verts,
			top_face_nv,
			top_local_to_global,
			top_global_to_local
		);

		std::vector<std::uint32_t> candidates;
		std::vector<std::uint32_t> wall_candidates;

		for (std::size_t local_fi = 0; local_fi < top_local_to_global.size(); ++local_fi) {
			const auto face_id = top_local_to_global[local_fi];
			if (face_id >= nf || !data_->face_active[face_id]) {
				continue;
			}

			const auto& face = data_->surface_faces[face_id];
			const AABB box = make_face_aabb(top_positions, face);
			top_bvh.query_overlap(box, candidates);
			for (auto local_other : candidates) {
				if (local_other <= local_fi) {
					continue;
				}

				const auto other_id = top_local_to_global[local_other];
				if (other_id >= nf || !data_->face_active[other_id]) {
					continue;
				}

				const auto& other = data_->surface_faces[other_id];
				if (faces_share_vertex(face, other)) {
					continue;
				}

				if (!faces_intersect_strict(top_positions, face, top_positions, other)) {
					continue;
				}

				std::ostringstream details;
				details << "### Trigger\n\n";
				details << "- check: `check_top_face_intersections()`\n";
				details << "- predicate: `faces_intersect_strict(top_positions, face, top_positions, other)`\n";
				details << "- face_id_a: `" << face_id << "`\n";
				details << "- face_id_b: `" << other_id << "`\n";
				export_first_problem_debug(
					"top_face_intersection",
					static_cast<std::int32_t>(face_id),
					{
					  copy_face_points(top_positions, face),
					  copy_face_points(top_positions, other)
					},
					{ face.num_vertices, other.num_vertices },
					details.str(),
					{}
				);

				all_ok = false;
				modify_faces.insert(face_id);
				modify_faces.insert(other_id);
				for (int i = 0; i < face.num_vertices; ++i) {
					modify_nodes.insert(face.vertices[i]);
				}
				for (int i = 0; i < other.num_vertices; ++i) {
					modify_nodes.insert(other.vertices[i]);
				}

				if (params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
					if (!params_.retain_valid_layers) {
						throw std::runtime_error("Boundary layer top faces intersect.");
					}
					return false;
				}
			}

			if (!data_->non_bl_face_verts.empty()) {
				data_->non_bl_bvh.query_overlap(box, wall_candidates);
				for (auto wall_local : wall_candidates) {
					SurfaceFace wall_face;
					wall_face.vertices = data_->non_bl_face_verts[wall_local];
					wall_face.num_vertices = data_->non_bl_face_nv[wall_local];
					if (!faces_intersect_strict(top_positions, face, data_->non_bl_nodes, wall_face)) {
						continue;
					}

					std::ostringstream details;
					details << "### Trigger\n\n";
					details << "- check: `check_top_face_intersections()`\n";
					details << "- predicate: `faces_intersect_strict(top_positions, face, data_->non_bl_nodes, wall_face)`\n";
					details << "- face_id: `" << face_id << "`\n";
					details << "- wall_face_local_id: `" << wall_local << "`\n";
					export_first_problem_debug(
						"top_face_wall_intersection",
						static_cast<std::int32_t>(face_id),
						{
						  copy_face_points(top_positions, face),
						  copy_face_points(data_->non_bl_nodes, wall_face)
						},
						{ face.num_vertices, wall_face.num_vertices },
						details.str(),
						{}
					);

					all_ok = false;
					modify_faces.insert(face_id);
					for (int i = 0; i < face.num_vertices; ++i) {
						modify_nodes.insert(face.vertices[i]);
					}
					if (params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
						if (!params_.retain_valid_layers) {
							throw std::runtime_error("Boundary layer top faces intersect wall boundary.");
						}
						return false;
					}
					break;
				}
			}
		}

		return all_ok;
	}

	bool BoundaryLayerSolver::check_side_face_intersections(
		const std::vector<Vec3>& bottom_positions,
		const std::vector<Vec3>& top_positions,
		std::unordered_set<std::uint32_t>& modify_nodes,
		std::unordered_set<std::uint32_t>& modify_faces
	) const
	{
		const auto nf = data_->surface_faces.size();
		bool all_ok = true;

		FaceBVH bottom_bvh;
		std::vector<std::array<std::uint32_t, 4>> bottom_face_verts;
		std::vector<int> bottom_face_nv;
		std::vector<std::uint32_t> bottom_local_to_global;
		std::vector<int> bottom_global_to_local;
		build_active_face_bvh(
			bottom_positions,
			bottom_bvh,
			bottom_face_verts,
			bottom_face_nv,
			bottom_local_to_global,
			bottom_global_to_local
		);

		FaceBVH top_bvh;
		std::vector<std::array<std::uint32_t, 4>> top_face_verts;
		std::vector<int> top_face_nv;
		std::vector<std::uint32_t> top_local_to_global;
		std::vector<int> top_global_to_local;
		build_active_face_bvh(
			top_positions,
			top_bvh,
			top_face_verts,
			top_face_nv,
			top_local_to_global,
			top_global_to_local
		);

		struct SweptSideFace {
			ExplicitFace face;
			AABB box;
			std::array<std::uint32_t, 2> edge_nodes{ 0, 0 };
			std::vector<std::uint32_t> incident_faces;
		};

		std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::uint32_t>> edge_to_faces;
		for (std::uint32_t fi = 0; fi < nf; ++fi) {
			if (!data_->face_active[fi]) {
				continue;
			}
			const auto& face = data_->surface_faces[fi];
			for (int ei = 0; ei < face.num_vertices; ++ei) {
				const auto v0 = face.vertices[ei];
				const auto v1 = face.vertices[(ei + 1) % face.num_vertices];
				const auto key = std::minmax(v0, v1);
				auto& uses = edge_to_faces[{key.first, key.second}];
				if (std::find(uses.begin(), uses.end(), fi) == uses.end()) {
					uses.push_back(fi);
				}
			}
		}

		std::vector<SweptSideFace> side_faces;
		side_faces.reserve(edge_to_faces.size());
		std::vector<Vec3> side_positions;
		side_positions.reserve(edge_to_faces.size() * 4U);
		std::vector<std::array<std::uint32_t, 4>> side_face_verts;
		side_face_verts.reserve(edge_to_faces.size());
		std::vector<int> side_face_nv;
		side_face_nv.reserve(edge_to_faces.size());

		for (const auto& entry : edge_to_faces) {
			const auto v0 = entry.first.first;
			const auto v1 = entry.first.second;

			ExplicitFace side;
			side.num_vertices = 4;
			side.points[0] = bottom_positions[v0];
			side.points[1] = bottom_positions[v1];
			side.points[2] = top_positions[v1];
			side.points[3] = top_positions[v0];

			const double lift0 = vec3_length(vec3_sub(side.points[3], side.points[0]));
			const double lift1 = vec3_length(vec3_sub(side.points[2], side.points[1]));
			if (lift0 <= kIntersectionEps && lift1 <= kIntersectionEps) {
				continue;
			}
			if (explicit_face_area(side) <= 1e-14) {
				continue;
			}

			SweptSideFace side_face;
			side_face.face = side;
			side_face.box = make_explicit_face_aabb(side);
			side_face.edge_nodes = { v0, v1 };
			side_face.incident_faces = entry.second;
			side_faces.push_back(side_face);

			const auto base = static_cast<std::uint32_t>(side_positions.size());
			side_positions.push_back(side.points[0]);
			side_positions.push_back(side.points[1]);
			side_positions.push_back(side.points[2]);
			side_positions.push_back(side.points[3]);
			side_face_verts.push_back({ base, base + 1U, base + 2U, base + 3U });
			side_face_nv.push_back(4);
		}

		if (side_faces.empty()) {
			return all_ok;
		}

		const auto mark_face_and_nodes = [&](std::uint32_t face_id) {
			if (face_id >= nf) {
				return;
			}
			modify_faces.insert(face_id);
			const auto& face = data_->surface_faces[face_id];
			for (int vi = 0; vi < face.num_vertices; ++vi) {
				modify_nodes.insert(face.vertices[vi]);
			}
			};

		const auto mark_side = [&](const SweptSideFace& side) {
			modify_nodes.insert(side.edge_nodes[0]);
			modify_nodes.insert(side.edge_nodes[1]);
			for (const auto face_id : side.incident_faces) {
				mark_face_and_nodes(face_id);
			}
			};

		const auto is_incident_face =
			[](const SweptSideFace& side, std::uint32_t face_id) {
			return std::find(side.incident_faces.begin(), side.incident_faces.end(), face_id) !=
				side.incident_faces.end();
			};

		const auto shares_side_endpoint =
			[](const SweptSideFace& side, const SurfaceFace& face) {
			return face_contains_node(face, side.edge_nodes[0]) ||
				face_contains_node(face, side.edge_nodes[1]);
			};

		const auto handle_collision = [&](const char* message) -> bool {
			if (params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
				if (!params_.retain_valid_layers) {
					throw std::runtime_error(message);
				}
				return true;
			}
			return false;
			};

		std::vector<std::uint32_t> candidates;
		std::vector<std::uint32_t> wall_candidates;

		for (std::size_t side_index = 0; side_index < side_faces.size(); ++side_index) {
			const auto& side = side_faces[side_index];
			bool hit = false;

			if (!bottom_face_verts.empty()) {
				bottom_bvh.query_overlap(side.box, candidates);
				for (const auto local_face : candidates) {
					const auto face_id = bottom_local_to_global[local_face];
					if (face_id >= nf || !data_->face_active[face_id] || is_incident_face(side, face_id)) {
						continue;
					}
					const auto& face = data_->surface_faces[face_id];
					if (shares_side_endpoint(side, face)) {
						continue;
					}
					if (!explicit_face_intersects_strict(side.face, bottom_positions, face)) {
						continue;
					}

					std::ostringstream details;
					details << "### Trigger\n\n";
					details << "- check: `check_side_face_intersections()`\n";
					details << "- predicate: `explicit_face_intersects_strict(side.face, bottom_positions, face)`\n";
					details << "- hit_face_id: `" << face_id << "`\n";
					details << "- swept_edge_nodes: `" << side.edge_nodes[0] << ", " << side.edge_nodes[1] << "`\n";
					export_first_problem_debug(
						"side_face_current_front_intersection",
						static_cast<std::int32_t>(face_id),
						{
						  copy_face_points(side.face),
						  copy_face_points(bottom_positions, face)
						},
						{ side.face.num_vertices, face.num_vertices },
						details.str(),
						{}
					);

					all_ok = false;
					mark_side(side);
					mark_face_and_nodes(face_id);
					if (handle_collision("Boundary layer side wall intersects current front faces.")) {
						return false;
					}
					hit = true;
					break;
				}
			}
			if (hit) {
				continue;
			}

			if (!top_face_verts.empty()) {
				top_bvh.query_overlap(side.box, candidates);
				for (const auto local_face : candidates) {
					const auto face_id = top_local_to_global[local_face];
					if (face_id >= nf || !data_->face_active[face_id] || is_incident_face(side, face_id)) {
						continue;
					}
					const auto& face = data_->surface_faces[face_id];
					if (shares_side_endpoint(side, face)) {
						continue;
					}
					if (!explicit_face_intersects_strict(side.face, top_positions, face)) {
						continue;
					}

					std::ostringstream details;
					details << "### Trigger\n\n";
					details << "- check: `check_side_face_intersections()`\n";
					details << "- predicate: `explicit_face_intersects_strict(side.face, top_positions, face)`\n";
					details << "- hit_face_id: `" << face_id << "`\n";
					details << "- swept_edge_nodes: `" << side.edge_nodes[0] << ", " << side.edge_nodes[1] << "`\n";
					export_first_problem_debug(
						"side_face_top_front_intersection",
						static_cast<std::int32_t>(face_id),
						{
						  copy_face_points(side.face),
						  copy_face_points(top_positions, face)
						},
						{ side.face.num_vertices, face.num_vertices },
						details.str(),
						{}
					);

					all_ok = false;
					mark_side(side);
					mark_face_and_nodes(face_id);
					if (handle_collision("Boundary layer side wall intersects pushed top faces.")) {
						return false;
					}
					hit = true;
					break;
				}
			}
			if (hit) {
				continue;
			}

			if (!data_->non_bl_face_verts.empty()) {
				data_->non_bl_bvh.query_overlap(side.box, wall_candidates);
				for (const auto wall_local : wall_candidates) {
					SurfaceFace wall_face;
					wall_face.vertices = data_->non_bl_face_verts[wall_local];
					wall_face.num_vertices = data_->non_bl_face_nv[wall_local];
					if (!explicit_face_intersects_strict(side.face, data_->non_bl_nodes, wall_face)) {
						continue;
					}

					std::ostringstream details;
					details << "### Trigger\n\n";
					details << "- check: `check_side_face_intersections()`\n";
					details << "- predicate: `explicit_face_intersects_strict(side.face, data_->non_bl_nodes, wall_face)`\n";
					details << "- wall_face_local_id: `" << wall_local << "`\n";
					details << "- swept_edge_nodes: `" << side.edge_nodes[0] << ", " << side.edge_nodes[1] << "`\n";
					export_first_problem_debug(
						"side_face_wall_intersection",
						-1,
						{
						  copy_face_points(side.face),
						  copy_face_points(data_->non_bl_nodes, wall_face)
						},
						{ side.face.num_vertices, wall_face.num_vertices },
						details.str(),
						{}
					);

					all_ok = false;
					mark_side(side);
					if (handle_collision("Boundary layer side wall intersects non-BL wall faces.")) {
						return false;
					}
					hit = true;
					break;
				}
			}
		}

		FaceBVH side_bvh;
		side_bvh.build(side_positions, side_face_verts, side_face_nv);
		std::vector<std::uint32_t> side_candidates;
		for (std::size_t side_index = 0; side_index < side_faces.size(); ++side_index) {
			SurfaceFace side_face;
			side_face.vertices = side_face_verts[side_index];
			side_face.num_vertices = 4;

			side_bvh.query_overlap(side_faces[side_index].box, side_candidates);
			for (const auto local_other : side_candidates) {
				if (local_other <= side_index) {
					continue;
				}
				if (side_faces[side_index].edge_nodes[0] == side_faces[local_other].edge_nodes[0] ||
					side_faces[side_index].edge_nodes[0] == side_faces[local_other].edge_nodes[1] ||
					side_faces[side_index].edge_nodes[1] == side_faces[local_other].edge_nodes[0] ||
					side_faces[side_index].edge_nodes[1] == side_faces[local_other].edge_nodes[1]) {
					continue;
				}
				if (!explicit_faces_intersect_strict(side_faces[side_index].face, side_faces[local_other].face)) {
					continue;
				}

				std::ostringstream details;
				details << "### Trigger\n\n";
				details << "- check: `check_side_face_intersections()`\n";
				details << "- predicate: `explicit_faces_intersect_strict(side_a, side_b)`\n";
				details << "- side_index_a: `" << side_index << "`\n";
				details << "- side_index_b: `" << local_other << "`\n";
				details << "- swept_edge_nodes_a: `" << side_faces[side_index].edge_nodes[0] << ", "
					<< side_faces[side_index].edge_nodes[1] << "`\n";
				details << "- swept_edge_nodes_b: `" << side_faces[local_other].edge_nodes[0] << ", "
					<< side_faces[local_other].edge_nodes[1] << "`\n";
				export_first_problem_debug(
					"side_face_self_intersection",
					-1,
					{
					  copy_face_points(side_faces[side_index].face),
					  copy_face_points(side_faces[local_other].face)
					},
		{
		  side_faces[side_index].face.num_vertices,
		  side_faces[local_other].face.num_vertices
		},
					details.str(),
					{}
				);

				all_ok = false;
				mark_side(side_faces[side_index]);
				mark_side(side_faces[local_other]);
				if (handle_collision("Boundary layer side walls intersect each other.")) {
					return false;
				}
			}
		}

		return all_ok;
	}

	bool BoundaryLayerSolver::check_top_face_quality(
		const std::vector<Vec3>& bottom_positions,
		const std::vector<Vec3>& top_positions,
		std::unordered_set<std::uint32_t>& modify_nodes,
		std::unordered_set<std::uint32_t>& modify_faces
	) const
	{
		const auto nf = data_->surface_faces.size();
		bool all_ok = true;

		for (std::uint32_t fi = 0; fi < nf; ++fi) {
			if (!data_->face_active[fi]) {
				continue;
			}

			const auto& face = data_->surface_faces[fi];
			bool bad = false;

			if (face.num_vertices == 3) {
				const auto& b0 = bottom_positions[face.vertices[0]];
				const auto& b1 = bottom_positions[face.vertices[1]];
				const auto& b2 = bottom_positions[face.vertices[2]];
				const auto& t0 = top_positions[face.vertices[0]];
				const auto& t1 = top_positions[face.vertices[1]];
				const auto& t2 = top_positions[face.vertices[2]];

				const double base_skew = compute_skewness_tri(b0, b1, b2);
				const double top_skew = compute_skewness_tri(t0, t1, t2);
				const bool negative_jac = prism_has_negative_jacobian(b0, b1, b2, t0, t1, t2);
				const double volume_skew = negative_jac
					? 1.0
					: compute_prism_skewness(b0, b1, b2, t0, t1, t2);
				const bool top_skew_fail = top_skew > params_.max_skewness;
				const bool negative_jac_fail = negative_jac;
				const bool volume_skew_fail = !negative_jac && volume_skew > params_.max_skewness;

				bad = top_skew_fail || negative_jac_fail || volume_skew_fail;
				if (bad) {
					const Vec3 base_normal = calc_face_normal_from_three(b0, b1, b2);
					const Vec3 top_normal = calc_face_normal_from_three(t0, t1, t2);
					const double top_normal_alignment = std::clamp(
						std::abs(vec3_dot(base_normal, top_normal)),
						0.0,
						1.0
					);
					const double top_normal_misalignment = 1.0 - top_normal_alignment;

					const std::array<Vec3, 3> base = { b0, b1, b2 };
					const std::array<Vec3, 3> top = { t0, t1, t2 };
					std::array<double, 3> side_lengths{};
					std::array<double, 3> side_alignments{};
					std::array<double, 3> side_tilts{};
					std::array<double, 3> projected_heights{};
					double max_side_tilt = 0.0;
					for (int i = 0; i < 3; ++i) {
						const Vec3 side = vec3_sub(top[static_cast<std::size_t>(i)], base[static_cast<std::size_t>(i)]);
						side_lengths[static_cast<std::size_t>(i)] = vec3_length(side);
						if (side_lengths[static_cast<std::size_t>(i)] > kEPS) {
							side_alignments[static_cast<std::size_t>(i)] = std::clamp(
								std::abs(vec3_dot(side, base_normal)) / side_lengths[static_cast<std::size_t>(i)],
								0.0,
								1.0
							);
						}
						side_tilts[static_cast<std::size_t>(i)] = std::sqrt(std::max(
							0.0,
							1.0 - side_alignments[static_cast<std::size_t>(i)] * side_alignments[static_cast<std::size_t>(i)]
						));
						projected_heights[static_cast<std::size_t>(i)] = vec3_dot(side, base_normal);
						max_side_tilt = std::max(max_side_tilt, side_tilts[static_cast<std::size_t>(i)]);
					}

					const char* primary_reason = negative_jac_fail
						? "top_face_quality_negative_jacobian"
						: (volume_skew_fail
							? "top_face_quality_volume_skew"
							: "top_face_quality_top_skewness");

					std::ostringstream details;
					details << std::fixed << std::setprecision(10);
					details << "### Trigger\n\n";
					details << "- check: `check_top_face_quality()`\n";
					details << "- face_id: `" << fi << "`\n";
					details << "- cell_type: `prism`\n";
					details << "- primary_reason: `" << primary_reason << "`\n";
					details << "- external_bl_max_skewness: `" << params_.max_skewness << "`\n";
					details << "- top_skew_fail: `" << (top_skew_fail ? "true" : "false") << "`\n";
					details << "- negative_jacobian_fail: `" << (negative_jac_fail ? "true" : "false") << "`\n";
					details << "- volume_skew_fail: `" << (volume_skew_fail ? "true" : "false") << "`\n";
					details << "\n### Scalars\n\n";
					details << "- base_skew: `" << base_skew << "`\n";
					details << "- top_skew: `" << top_skew << "`\n";
					details << "- top_normal_alignment: `" << top_normal_alignment << "`\n";
					details << "- top_normal_misalignment: `" << top_normal_misalignment << "`\n";
					details << "- max_side_tilt: `" << max_side_tilt << "`\n";
					details << "- volume_skew: `" << volume_skew << "`\n";
					details << "\n### Formula\n\n";
					details << "- `top_skew = compute_skewness_tri(t0, t1, t2)`\n";
					details << "- `negative_jacobian = prism_has_negative_jacobian(b0, b1, b2, t0, t1, t2)`\n";
					details << "- `tri_skew = max(compute_skewness_tri(b0, b1, b2), compute_skewness_tri(t0, t1, t2))`\n";
					details << "- `side_i = t_i - b_i`\n";
					details << "- `alignment_i = abs(dot(side_i, base_normal)) / |side_i|`\n";
					details << "- `tilt_i = sqrt(1 - alignment_i^2)`\n";
					details << "- `top_normal_misalignment = 1 - abs(dot(base_normal, top_normal))`\n";
					details << "- `volume_skew = max(tri_skew, max_side_tilt, top_normal_misalignment)`\n";
					details << "- `volume_skew` gate uses external `bl_max_skewness`, i.e. `params_.max_skewness`\n";
					details << "\n### Side Metrics\n\n";
					for (int i = 0; i < 3; ++i) {
						details << "- side " << i
							<< ": len=`" << side_lengths[static_cast<std::size_t>(i)]
							<< "`, alignment=`" << side_alignments[static_cast<std::size_t>(i)]
							<< "`, tilt=`" << side_tilts[static_cast<std::size_t>(i)]
							<< "`, projected_height=`" << projected_heights[static_cast<std::size_t>(i)]
							<< "`\n";
					}

					export_first_problem_debug(
						primary_reason,
						static_cast<std::int32_t>(fi),
						{
						  make_debug_face(b0, b1, b2),
						  make_debug_face(t0, t1, t2),
						  make_debug_face(b0, b1, t1, t0),
						  make_debug_face(b1, b2, t2, t1),
						  make_debug_face(b2, b0, t0, t2)
						},
						{ 3, 3, 4, 4, 4 },
						details.str(),
						{}
					);
				}
			}
			else if (face.num_vertices == 4) {
				const auto& b0 = bottom_positions[face.vertices[0]];
				const auto& b1 = bottom_positions[face.vertices[1]];
				const auto& b2 = bottom_positions[face.vertices[2]];
				const auto& b3 = bottom_positions[face.vertices[3]];
				const auto& t0 = top_positions[face.vertices[0]];
				const auto& t1 = top_positions[face.vertices[1]];
				const auto& t2 = top_positions[face.vertices[2]];
				const auto& t3 = top_positions[face.vertices[3]];

				const double base_skew = compute_skewness_quad(b0, b1, b2, b3);
				const double top_skew = compute_skewness_quad(t0, t1, t2, t3);
				const double warping = compute_warping_quad(t0, t1, t2, t3);
				const bool negative_jac = hexa_has_negative_jacobian(b0, b1, b2, b3, t0, t1, t2, t3);
				const double volume_skew = negative_jac
					? 1.0
					: compute_hexa_skewness(b0, b1, b2, b3, t0, t1, t2, t3);
				const bool top_skew_fail = top_skew > params_.max_skewness;
				const bool warping_fail = warping > params_.max_warping_deg;
				const bool negative_jac_fail = negative_jac;
				const bool volume_skew_fail = !negative_jac && volume_skew > params_.max_skewness;

				bad = top_skew_fail || warping_fail || negative_jac_fail || volume_skew_fail;
				if (bad) {
					const Vec3 base_normal = vec3_normalized(vec3_add(
						vec3_cross(vec3_sub(b1, b0), vec3_sub(b2, b0)),
						vec3_cross(vec3_sub(b2, b0), vec3_sub(b3, b0))
					));
					const Vec3 top_normal = vec3_normalized(vec3_add(
						vec3_cross(vec3_sub(t1, t0), vec3_sub(t2, t0)),
						vec3_cross(vec3_sub(t2, t0), vec3_sub(t3, t0))
					));
					const double top_normal_alignment = std::clamp(
						std::abs(vec3_dot(base_normal, top_normal)),
						0.0,
						1.0
					);
					const double top_normal_misalignment = 1.0 - top_normal_alignment;

					const std::array<Vec3, 4> base = { b0, b1, b2, b3 };
					const std::array<Vec3, 4> top = { t0, t1, t2, t3 };
					std::array<double, 4> side_lengths{};
					std::array<double, 4> side_alignments{};
					std::array<double, 4> side_tilts{};
					std::array<double, 4> projected_heights{};
					double max_side_tilt = 0.0;
					for (int i = 0; i < 4; ++i) {
						const Vec3 side = vec3_sub(top[static_cast<std::size_t>(i)], base[static_cast<std::size_t>(i)]);
						side_lengths[static_cast<std::size_t>(i)] = vec3_length(side);
						if (side_lengths[static_cast<std::size_t>(i)] > kEPS) {
							side_alignments[static_cast<std::size_t>(i)] = std::clamp(
								std::abs(vec3_dot(side, base_normal)) / side_lengths[static_cast<std::size_t>(i)],
								0.0,
								1.0
							);
						}
						side_tilts[static_cast<std::size_t>(i)] = std::sqrt(std::max(
							0.0,
							1.0 - side_alignments[static_cast<std::size_t>(i)] * side_alignments[static_cast<std::size_t>(i)]
						));
						projected_heights[static_cast<std::size_t>(i)] = vec3_dot(side, base_normal);
						max_side_tilt = std::max(max_side_tilt, side_tilts[static_cast<std::size_t>(i)]);
					}

					const char* primary_reason = negative_jac_fail
						? "top_face_quality_negative_jacobian"
						: (volume_skew_fail
							? "top_face_quality_volume_skew"
							: (warping_fail
								? "top_face_quality_warping"
								: "top_face_quality_top_skewness"));

					std::ostringstream details;
					details << std::fixed << std::setprecision(10);
					details << "### Trigger\n\n";
					details << "- check: `check_top_face_quality()`\n";
					details << "- face_id: `" << fi << "`\n";
					details << "- cell_type: `hexa`\n";
					details << "- primary_reason: `" << primary_reason << "`\n";
					details << "- external_bl_max_skewness: `" << params_.max_skewness << "`\n";
					details << "- top_skew_fail: `" << (top_skew_fail ? "true" : "false") << "`\n";
					details << "- warping_fail: `" << (warping_fail ? "true" : "false") << "`\n";
					details << "- negative_jacobian_fail: `" << (negative_jac_fail ? "true" : "false") << "`\n";
					details << "- volume_skew_fail: `" << (volume_skew_fail ? "true" : "false") << "`\n";
					details << "\n### Scalars\n\n";
					details << "- base_skew: `" << base_skew << "`\n";
					details << "- top_skew: `" << top_skew << "`\n";
					details << "- warping_deg: `" << warping << "`\n";
					details << "- top_normal_alignment: `" << top_normal_alignment << "`\n";
					details << "- top_normal_misalignment: `" << top_normal_misalignment << "`\n";
					details << "- max_side_tilt: `" << max_side_tilt << "`\n";
					details << "- volume_skew: `" << volume_skew << "`\n";
					details << "\n### Formula\n\n";
					details << "- `top_skew = compute_skewness_quad(t0, t1, t2, t3)`\n";
					details << "- `warping = compute_warping_quad(t0, t1, t2, t3)`\n";
					details << "- `negative_jacobian = hexa_has_negative_jacobian(b0, b1, b2, b3, t0, t1, t2, t3)`\n";
					details << "- `quad_skew = max(compute_skewness_quad(b0, b1, b2, b3), compute_skewness_quad(t0, t1, t2, t3))`\n";
					details << "- `side_i = t_i - b_i`\n";
					details << "- `alignment_i = abs(dot(side_i, base_normal)) / |side_i|`\n";
					details << "- `tilt_i = sqrt(1 - alignment_i^2)`\n";
					details << "- `top_normal_misalignment = 1 - abs(dot(base_normal, top_normal))`\n";
					details << "- `volume_skew = max(quad_skew, max_side_tilt, top_normal_misalignment)`\n";
					details << "- `volume_skew` gate uses external `bl_max_skewness`, i.e. `params_.max_skewness`\n";
					details << "\n### Side Metrics\n\n";
					for (int i = 0; i < 4; ++i) {
						details << "- side " << i
							<< ": len=`" << side_lengths[static_cast<std::size_t>(i)]
							<< "`, alignment=`" << side_alignments[static_cast<std::size_t>(i)]
							<< "`, tilt=`" << side_tilts[static_cast<std::size_t>(i)]
							<< "`, projected_height=`" << projected_heights[static_cast<std::size_t>(i)]
							<< "`\n";
					}

					export_first_problem_debug(
						primary_reason,
						static_cast<std::int32_t>(fi),
						{
						  make_debug_face(b0, b1, b2, b3),
						  make_debug_face(t0, t1, t2, t3),
						  make_debug_face(b0, b1, t1, t0),
						  make_debug_face(b1, b2, t2, t1),
						  make_debug_face(b2, b3, t3, t2),
						  make_debug_face(b3, b0, t0, t3)
						},
						{ 4, 4, 4, 4, 4, 4 },
						details.str(),
						{}
					);
				}
			}

			if (!bad) {
				continue;
			}

			all_ok = false;
			modify_faces.insert(fi);
			for (int i = 0; i < face.num_vertices; ++i) {
				modify_nodes.insert(face.vertices[i]);
			}

			if (params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Stop) {
				if (!params_.retain_valid_layers) {
					throw std::runtime_error("Boundary layer generated invalid prism/hexa quality.");
				}
				return false;
			}
		}

		return all_ok;
	}

	void BoundaryLayerSolver::apply_problematic_area_treatment(
		int layer,
		std::unordered_set<std::uint32_t>& modify_nodes,
		std::unordered_set<std::uint32_t>& modify_faces
	)
	{
		auto& top = data_->layer_positions.back();
		const auto& prev = data_->layer_positions[data_->layer_positions.size() - 2];
		const auto& base = data_->layer_positions.front();
		const auto nn = data_->surface_nodes.size();

		if (params_.treatment == BoundaryLayerParams::ProblematicAreaTreatment::Exclude) {
			for (auto node : modify_nodes) {
				if (node >= nn) {
					continue;
				}
				data_->node_active[node] = false;
				data_->max_node_level[node] = 0;
				for (std::size_t ilayer = 1; ilayer < data_->layer_positions.size(); ++ilayer) {
					data_->layer_positions[ilayer][node] = base[node];
				}
				for (auto face : data_->node_to_faces[node]) {
					if (face >= data_->surface_faces.size()) {
						continue;
					}
					if (data_->face_active[face]) {
						data_->face_active[face] = false;
						data_->max_face_level[face] = 0;
						modify_faces.insert(face);
					}
				}
			}

			for (std::uint32_t node = 0; node < nn; ++node) {
				if (!data_->node_active[node]) {
					continue;
				}
				bool one_active = false;
				for (auto face : data_->node_to_faces[node]) {
					if (face < data_->face_active.size() && data_->face_active[face]) {
						one_active = true;
						break;
					}
				}
				if (one_active) {
					continue;
				}
				data_->node_active[node] = false;
				data_->max_node_level[node] = 0;
				for (std::size_t ilayer = 1; ilayer < data_->layer_positions.size(); ++ilayer) {
					data_->layer_positions[ilayer][node] = base[node];
				}
			}
			return;
		}

		for (auto node : modify_nodes) {
			if (node >= nn) {
				continue;
			}
			data_->node_active[node] = false;
			data_->max_node_level[node] = std::min(data_->max_node_level[node], layer);
			top[node] = prev[node];
			for (auto face : data_->node_to_faces[node]) {
				if (face >= data_->surface_faces.size()) {
					continue;
				}
				if (data_->face_active[face]) {
					data_->face_active[face] = false;
					data_->max_face_level[face] = std::min(data_->max_face_level[face], layer);
					modify_faces.insert(face);
				}
			}
		}

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) {
				continue;
			}
			bool one_active = false;
			for (auto face : data_->node_to_faces[node]) {
				if (face < data_->face_active.size() && data_->face_active[face]) {
					one_active = true;
					break;
				}
			}
			if (one_active) {
				continue;
			}
			data_->node_active[node] = false;
			data_->max_node_level[node] = std::min(data_->max_node_level[node], layer);
			top[node] = prev[node];
		}
	}

	void BoundaryLayerSolver::enforce_max_consecutive_collapsed_layers(
		const std::unordered_set<std::uint32_t>& seed_faces
	) {
		if (params_.max_consecutive_collapsed_layers <= 0 || seed_faces.empty()) {
			return;
		}

		std::vector<std::vector<int>> bucket(
			static_cast<std::size_t>(params_.max_consecutive_collapsed_layers) + 1U
		);
		for (auto face : seed_faces) {
			if (face < data_->surface_faces.size()) {
				bucket[0].push_back(static_cast<int>(face));
			}
		}

		int empty_count = 0;
		for (int b = 0; empty_count <= params_.max_consecutive_collapsed_layers;
			b = (b + 1) % (params_.max_consecutive_collapsed_layers + 1)) {
			auto& queue = bucket[static_cast<std::size_t>(b)];
			if (queue.empty()) {
				++empty_count;
				continue;
			}

			empty_count = 0;
			for (int face : queue) {
				if (face < 0 || face >= static_cast<int>(data_->surface_faces.size())) {
					continue;
				}
				const int level_f = data_->max_face_level[static_cast<std::size_t>(face)];

				std::unordered_set<std::uint32_t> related_faces;
				for (auto node : data_->surface_faces[static_cast<std::size_t>(face)].vertices) {
					if (node >= data_->node_to_faces.size()) {
						continue;
					}
					for (auto rface : data_->node_to_faces[node]) {
						related_faces.insert(rface);
					}
				}

				for (auto g : related_faces) {
					const int want = level_f + params_.max_consecutive_collapsed_layers;
					if (data_->max_face_level[g] <= want) {
						continue;
					}
					data_->max_face_level[g] = want;
					const int delay = want - level_f;
					const int nb = (b + delay) % (params_.max_consecutive_collapsed_layers + 1);
					bucket[static_cast<std::size_t>(nb)].push_back(static_cast<int>(g));
				}
			}

			queue.clear();
		}
	}

	void BoundaryLayerSolver::recompute_node_thickness()
	{
		const auto nn = data_->surface_nodes.size();
		data_->node_thickness.assign(nn, 0.0);

		if (data_->layer_positions.size() < 2) {
			return;
		}

		for (std::size_t layer = 1; layer < data_->layer_positions.size(); ++layer) {
			const auto& prev = data_->layer_positions[layer - 1];
			const auto& cur = data_->layer_positions[layer];
			for (std::uint32_t node = 0; node < nn; ++node) {
				data_->node_thickness[node] += vec3_length(vec3_sub(cur[node], prev[node]));
			}
		}
	}

	double BoundaryLayerSolver::compute_node_visibility_angle(std::uint32_t node) const
	{
		if (node >= data_->node_manifold_angle.size()) {
			return 0.0;
		}

		const double manifold = data_->node_manifold_angle[node];
		if (!std::isfinite(manifold) || manifold <= 0.0) {
			return 0.0;
		}

		return std::clamp(std::abs(manifold - kPI) * 0.5, 0.0, kHALF_PI - 1.0e-6);
	}

	// ===========================================================================
	// Quality metrics
	// ===========================================================================

	double BoundaryLayerSolver::compute_skewness_tri(
		const Vec3& a, const Vec3& b, const Vec3& c)
	{
		auto ab = vec3_sub(b, a);
		auto ac = vec3_sub(c, a);
		auto bc = vec3_sub(c, b);

		double lab = vec3_length(ab);
		double lac = vec3_length(ac);
		double lbc = vec3_length(bc);

		if (lab < kEPS || lac < kEPS || lbc < kEPS) return 1.0;

		// Equilateral skewness
		double area = 0.5 * vec3_length(vec3_cross(ab, ac));
		double s = (lab + lac + lbc) / 2.0;
		double r_in = area / s;
		double r_circ = (lab * lac * lbc) / (4.0 * area);

		if (r_circ < kEPS) return 1.0;

		// Ideal: r_in/r_circ = 0.5 for equilateral
		double quality = (r_in / r_circ) / 0.5;
		return 1.0 - std::clamp(quality, 0.0, 1.0);
	}

	double BoundaryLayerSolver::compute_skewness_quad(
		const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d)
	{
		// Average of two triangle skewnesses
		double s1 = compute_skewness_tri(a, b, c);
		double s2 = compute_skewness_tri(a, c, d);
		return std::max(s1, s2);
	}

	double BoundaryLayerSolver::compute_warping_quad(
		const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d)
	{
		// Warping: deviation of quad from planar
		auto n1 = calc_face_normal_from_three(a, b, c);
		auto n2 = calc_face_normal_from_three(a, c, d);
		double dot = vec3_dot(n1, n2);
		dot = std::clamp(dot, -1.0, 1.0);
		return std::acos(dot) * 180.0 / kPI;
	}

	double BoundaryLayerSolver::signed_corner_jacobian(
		const Vec3& e0, const Vec3& e1, const Vec3& e2)
	{
		const double denom = vec3_length(e0) * vec3_length(e1) * vec3_length(e2);
		if (denom < kEPS) {
			return 0.0;
		}
		return vec3_dot(e0, vec3_cross(e1, e2)) / denom;
	}

	bool BoundaryLayerSolver::prism_has_negative_jacobian(
		const Vec3& b0, const Vec3& b1, const Vec3& b2,
		const Vec3& t0, const Vec3& t1, const Vec3& t2)
	{
		constexpr double kProjectedHeightEps = 1.0e-8;
		const Vec3 base_normal = calc_face_normal_from_three(b0, b1, b2);
		if (vec3_length(base_normal) < kEPS) {
			return true;
		}

		const Vec3 top_normal = calc_face_normal_from_three(t0, t1, t2);
		if (vec3_length(top_normal) < kEPS || vec3_dot(base_normal, top_normal) <= kProjectedHeightEps) {
			return true;
		}

		const std::array<Vec3, 3> base = { b0, b1, b2 };
		const std::array<Vec3, 3> top = { t0, t1, t2 };
		double reference_sign = 0.0;

		for (int i = 0; i < 3; ++i) {
			const double projected_height = vec3_dot(vec3_sub(top[i], base[i]), base_normal);
			if (std::abs(projected_height) <= kProjectedHeightEps) {
				return true;
			}
			if (reference_sign == 0.0) {
				reference_sign = (projected_height > 0.0) ? 1.0 : -1.0;
				continue;
			}
			if (projected_height * reference_sign <= kProjectedHeightEps) {
				return true;
			}
		}

		return false;
	}

	bool BoundaryLayerSolver::hexa_has_negative_jacobian(
		const Vec3& b0, const Vec3& b1, const Vec3& b2, const Vec3& b3,
		const Vec3& t0, const Vec3& t1, const Vec3& t2, const Vec3& t3)
	{
		constexpr double kProjectedHeightEps = 1.0e-8;
		const Vec3 base_normal = vec3_normalized(vec3_add(
			vec3_cross(vec3_sub(b1, b0), vec3_sub(b2, b0)),
			vec3_cross(vec3_sub(b2, b0), vec3_sub(b3, b0))
		));
		if (vec3_length(base_normal) < kEPS) {
			return true;
		}

		const Vec3 top_normal = vec3_normalized(vec3_add(
			vec3_cross(vec3_sub(t1, t0), vec3_sub(t2, t0)),
			vec3_cross(vec3_sub(t2, t0), vec3_sub(t3, t0))
		));
		if (vec3_length(top_normal) < kEPS || vec3_dot(base_normal, top_normal) <= kProjectedHeightEps) {
			return true;
		}

		const std::array<Vec3, 4> base = { b0, b1, b2, b3 };
		const std::array<Vec3, 4> top = { t0, t1, t2, t3 };
		double reference_sign = 0.0;

		for (int i = 0; i < 4; ++i) {
			const double projected_height = vec3_dot(vec3_sub(top[i], base[i]), base_normal);
			if (std::abs(projected_height) <= kProjectedHeightEps) {
				return true;
			}
			if (reference_sign == 0.0) {
				reference_sign = (projected_height > 0.0) ? 1.0 : -1.0;
				continue;
			}
			if (projected_height * reference_sign <= kProjectedHeightEps) {
				return true;
			}
		}

		return false;
	}

	double BoundaryLayerSolver::compute_prism_skewness(
		const Vec3& b0, const Vec3& b1, const Vec3& b2,
		const Vec3& t0, const Vec3& t1, const Vec3& t2)
	{
		const Vec3 base_normal = calc_face_normal_from_three(b0, b1, b2);
		const Vec3 top_normal = calc_face_normal_from_three(t0, t1, t2);
		if (vec3_length(base_normal) < kEPS || vec3_length(top_normal) < kEPS) {
			return 1.0;
		}

		const double tri_skew = std::max(
			compute_skewness_tri(b0, b1, b2),
			compute_skewness_tri(t0, t1, t2)
		);

		double max_side_tilt = 0.0;
		const std::array<Vec3, 3> base = { b0, b1, b2 };
		const std::array<Vec3, 3> top = { t0, t1, t2 };
		for (int i = 0; i < 3; ++i) {
			const auto side = vec3_sub(top[i], base[i]);
			const double len_side = vec3_length(side);
			if (len_side < kEPS) {
				return 1.0;
			}
			const double normal_alignment = std::clamp(
				std::abs(vec3_dot(side, base_normal)) / len_side,
				0.0,
				1.0
			);
			const double side_tilt = std::sqrt(std::max(0.0, 1.0 - normal_alignment * normal_alignment));
			max_side_tilt = std::max(max_side_tilt, side_tilt);
		}

		const double top_normal_misalignment = 1.0 - std::clamp(std::abs(vec3_dot(base_normal, top_normal)), 0.0, 1.0);
		return std::max({ tri_skew, max_side_tilt, top_normal_misalignment });
	}

	double BoundaryLayerSolver::compute_hexa_skewness(
		const Vec3& b0, const Vec3& b1, const Vec3& b2, const Vec3& b3,
		const Vec3& t0, const Vec3& t1, const Vec3& t2, const Vec3& t3)
	{
		const Vec3 base_normal = vec3_normalized(vec3_add(
			vec3_cross(vec3_sub(b1, b0), vec3_sub(b2, b0)),
			vec3_cross(vec3_sub(b2, b0), vec3_sub(b3, b0))
		));
		const Vec3 top_normal = vec3_normalized(vec3_add(
			vec3_cross(vec3_sub(t1, t0), vec3_sub(t2, t0)),
			vec3_cross(vec3_sub(t2, t0), vec3_sub(t3, t0))
		));
		if (vec3_length(base_normal) < kEPS || vec3_length(top_normal) < kEPS) {
			return 1.0;
		}

		const double quad_skew = std::max(
			compute_skewness_quad(b0, b1, b2, b3),
			compute_skewness_quad(t0, t1, t2, t3)
		);

		const std::array<Vec3, 4> base = { b0, b1, b2, b3 };
		const std::array<Vec3, 4> top = { t0, t1, t2, t3 };
		double max_side_tilt = 0.0;
		for (int i = 0; i < 4; ++i) {
			const auto side = vec3_sub(top[i], base[i]);
			const double len_side = vec3_length(side);
			if (len_side < kEPS) {
				return 1.0;
			}
			const double normal_alignment = std::clamp(
				std::abs(vec3_dot(side, base_normal)) / len_side,
				0.0,
				1.0
			);
			const double side_tilt = std::sqrt(std::max(0.0, 1.0 - normal_alignment * normal_alignment));
			max_side_tilt = std::max(max_side_tilt, side_tilt);
		}

		const double top_normal_misalignment = 1.0 - std::clamp(std::abs(vec3_dot(base_normal, top_normal)), 0.0, 1.0);
		return std::max({ quad_skew, max_side_tilt, top_normal_misalignment });
	}

	std::vector<double> BoundaryLayerSolver::compute_distance_coefficient() const
	{
		const auto nn = data_->surface_nodes.size();
		std::vector<double> alphas(nn, 1.0);

		for (std::uint32_t node = 0; node < nn; ++node) {
			if (!data_->node_active[node]) continue;

			double betai = data_->node_manifold_angle[node];
			double grow = 1.0;

			if (betai > kPI) {
				// Convex node: angle > 180°, layers spread out
				grow = 1.0 / std::sin(kPI - betai / 2.0);
			}
			else if (betai > 0) {
				// Concave node: angle < 180°, layers compress -> need to grow MORE
				grow = 1.0 / std::sin(betai / 2.0);
			}

			grow = std::min(grow, 1.5); // Cap growth factor

			// Additional smoothing factor based on angle
			double alp = 1.0 + 0.5 * std::cos(betai / 2.0);
			alphas[node] = alp * grow;
		}

		return alphas;
	}

} // namespace sqmesh::mesh
