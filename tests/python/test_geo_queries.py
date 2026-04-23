# Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
# Licensed under the GNU Affero General Public License v3.0 or later.
# See LICENSE file in the project root for details.
#
import tempfile
import textwrap
import unittest
from pathlib import Path

import sqmesh.base as base
import sqmesh.geo as geo


def write_patch_stl(path: Path) -> None:
    path.write_text(
        textwrap.dedent(
            """\
            solid square_patch
              facet normal 0 0 1
                outer loop
                  vertex 0 0 0
                  vertex 1 0 0
                  vertex 1 1 0
                endloop
              endfacet
              facet normal 0 0 1
                outer loop
                  vertex 0 0 0
                  vertex 1 1 0
                  vertex 0 1 0
                endloop
              endfacet
            endsolid square_patch
            """
        ),
        encoding="utf-8",
    )


class GeoQueryTest(unittest.TestCase):
    def tearDown(self) -> None:
        try:
            base.shutdown_all()
        except base.SQMeshError:
            pass

    def test_stl_topology_and_view_queries(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir, base.Context() as context:
            stl_path = Path(temp_dir) / "patch.stl"
            write_patch_stl(stl_path)

            model = geo.import_stl(str(stl_path), context=context)
            summary = model.summary()
            self.assertEqual(summary.face_count, 1)
            self.assertEqual(summary.edge_count, 5)
            self.assertEqual(summary.vertex_count, 4)

            snapshot = model.topology_snapshot()
            self.assertEqual(snapshot.entity_count(geo.TopologyDimension.face), 1)
            self.assertEqual(snapshot.entity_count(geo.TopologyDimension.edge), 5)
            self.assertEqual(snapshot.entity_count(geo.TopologyDimension.vertex), 4)

            view = model.view()
            self.assertEqual(len(view.faces), 1)
            face = view.faces[0]
            boundary = face.boundary_loops()
            self.assertEqual(boundary.unknown_loop_count, 1)
            self.assertTrue(boundary.loops[0].closed)
            self.assertEqual(len(boundary.loops[0].vertex_ids), 4)

            feature_report = geo.feature_edges(model, context=context)
            self.assertEqual(len(feature_report.edges), 4)

            with self.assertRaises(base.UnsupportedError):
                face.uv_bounds()

            with self.assertRaises(base.UnsupportedError):
                face.sample(0.0, 0.0)

    def test_missing_file_maps_to_io_error(self) -> None:
        with base.Context() as context:
            with self.assertRaises(base.IoError):
                geo.import_stl("missing.stl", context=context)


if __name__ == "__main__":
    unittest.main()
