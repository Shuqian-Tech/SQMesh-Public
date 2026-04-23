// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/base/api.hpp"

#include <cstdlib>

int main()
{
  sqmesh::base::ContextHandle context = sqmesh::invalid_handle;
  return context == sqmesh::invalid_handle ? EXIT_SUCCESS : EXIT_FAILURE;
}
