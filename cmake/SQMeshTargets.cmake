include_guard(GLOBAL)

function(sqmesh_configure_target target_name)
  target_compile_features(${target_name} PUBLIC cxx_std_17)
  set_target_properties(${target_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_include_directories(
    ${target_name}
    PUBLIC
      $<BUILD_INTERFACE:${SQMESH_GENERATED_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:include>
  )

  if(MSVC)
    target_compile_options(${target_name} PRIVATE /W4 /permissive- /utf-8)
  else()
    target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()
