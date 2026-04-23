# Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
# Licensed under the GNU Affero General Public License v3.0 or later.
# See LICENSE file in the project root for details.
#
import tempfile
import unittest
from pathlib import Path

import sqmesh.base as base
import sqmesh.geo as geo
import sqmesh.mesh as mesh


class MeshWorkflowTest(unittest.TestCase):
    def tearDown(self) -> None:
        try:
            base.shutdown_all()
        except base.SQMeshError:
            pass

    def test_dummy_volume_mesher_and_msh_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir, base.Context() as context:
            model = geo.create_placeholder_model(context=context)
            output_mesh = mesh.create_volume_mesh(
                model,
                "dummy_mesher",
                context=context,
                parameters={
                    "minimum_length": 0.25,
                    "growth_rate": 1.2,
                    "element_type": "tetra",
                },
            )

            summary = output_mesh.summary()
            self.assertEqual(summary.node_count, 4)
            self.assertEqual(summary.face_count, 4)
            self.assertEqual(summary.cell_count, 1)
            self.assertEqual(output_mesh.nodes_count(), 4)
            self.assertEqual(output_mesh.cells_count(), 1)

            domain = output_mesh.domain_snapshot()
            self.assertEqual(domain.entity_group_count(), 3)
            self.assertEqual(domain.summary().cell_count, 1)

            report = output_mesh.quality_report()
            self.assertGreaterEqual(report.supported_element_count, 1)

            msh_path = Path(temp_dir) / "dummy.msh"
            output_mesh.export_msh(str(msh_path))
            imported_mesh = mesh.import_msh(str(msh_path), context=context)
            imported_summary = imported_mesh.summary()
            self.assertEqual(imported_summary.node_count, 4)
            self.assertEqual(imported_summary.face_count, 4)
            self.assertEqual(imported_summary.cell_count, 1)

    def test_python_dict_parameters_reject_unsupported_value_types(self) -> None:
        with base.Context() as context:
            model = geo.create_placeholder_model(context=context)
            with self.assertRaises(TypeError):
                mesh.create_volume_mesh(
                    model,
                    "dummy_mesher",
                    context=context,
                    parameters={"bad": object()},
                )

    def test_boolean_parameter_values_preserve_boolean_semantics(self) -> None:
        direct_value = mesh.ParameterValue(True)
        self.assertEqual(direct_value.type(), mesh.ParameterType.boolean)
        self.assertTrue(direct_value.is_boolean())
        self.assertFalse(direct_value.is_integer())
        self.assertEqual(direct_value.to_python(), True)

        converted = mesh._coerce_parameters({"enabled": True, "disabled": False})
        self.assertEqual(converted.try_get_boolean("enabled"), True)
        self.assertEqual(converted.try_get_boolean("disabled"), False)
        self.assertIsNone(converted.try_get_integer("enabled"))
        self.assertEqual(converted.get("enabled"), True)


if __name__ == "__main__":
    unittest.main()
