# Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
# Licensed under the GNU Affero General Public License v3.0 or later.
# See LICENSE file in the project root for details.
#
import unittest

import sqmesh.base as base
import sqmesh.geo as geo
import sqmesh.mesh as mesh


class BaseRuntimeTest(unittest.TestCase):
    def tearDown(self) -> None:
        try:
            base.shutdown_all()
        except base.SQMeshError:
            pass

    def test_context_lifecycle_and_current_bindings(self) -> None:
        with base.Context() as context:
            self.assertTrue(base.is_initialized())
            current = base.current_context()
            self.assertIsNotNone(current)
            self.assertEqual(current.handle, context.handle)
            self.assertNotEqual(base.current_session(context), base.INVALID_HANDLE)

        self.assertIsNone(base.current_context())
        self.assertFalse(base.is_initialized())

    def test_not_initialized_error_mapping(self) -> None:
        with self.assertRaises(base.NotInitializedError):
            geo.create_placeholder_model()

    def test_stale_wrappers_do_not_fall_back_to_new_current_context(self) -> None:
        first_context = base.Context()
        model = geo.create_placeholder_model(context=first_context)
        output_mesh = mesh.create_volume_mesh(
            model,
            "dummy_mesher",
            context=first_context,
            parameters={
                "minimum_length": 0.25,
                "growth_rate": 1.2,
                "element_type": "tetra",
            },
        )
        first_context.close()

        with base.Context() as second_context:
            self.assertEqual(base.current_context().handle, second_context.handle)

            with self.assertRaises(base.InvalidHandleError):
                base.current_session(first_context)

            with self.assertRaises(base.InvalidHandleError):
                model.summary()

            with self.assertRaises(base.InvalidHandleError):
                output_mesh.summary()

            with self.assertRaises(base.InvalidHandleError):
                mesh.create_volume_mesh(
                    model,
                    "dummy_mesher",
                    context=second_context,
                    parameters={
                        "minimum_length": 0.25,
                        "growth_rate": 1.2,
                        "element_type": "tetra",
                    },
                )

            with self.assertRaises(base.InvalidHandleError):
                mesh.create_surface_mesh(
                    model,
                    "auto_cfd_surface_mesher",
                    context=second_context,
                    parameters={"minimum_length": 1.0, "element_type": "tri"},
                )


if __name__ == "__main__":
    unittest.main()
