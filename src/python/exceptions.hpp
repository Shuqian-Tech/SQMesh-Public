// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"

#include <stdexcept>
#include <string>

class SQMeshError : public std::runtime_error
{
public:
  using std::runtime_error::runtime_error;
};

class InvalidArgumentError : public SQMeshError
{
public:
  using SQMeshError::SQMeshError;
};

class NotInitializedError : public SQMeshError
{
public:
  using SQMeshError::SQMeshError;
};

class InvalidHandleError : public SQMeshError
{
public:
  using SQMeshError::SQMeshError;
};

class OwnerMismatchError : public SQMeshError
{
public:
  using SQMeshError::SQMeshError;
};

class IoError : public SQMeshError
{
public:
  using SQMeshError::SQMeshError;
};

class UnsupportedError : public SQMeshError
{
public:
  using SQMeshError::SQMeshError;
};

class InternalError : public SQMeshError
{
public:
  using SQMeshError::SQMeshError;
};

[[noreturn]] void throw_sqmesh_error(sqmesh::base::StatusCode status);
void throw_on_status(sqmesh::base::StatusCode status);
