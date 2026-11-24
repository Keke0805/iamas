list(APPEND _iamas_operators_src
    iamas_geometry_utils.F
    iamas_matrix_operations.F
    iamas_rbf_interpolation.F
    iamas_spline_interpolation.F
    iamas_tensor_operations.F
    iamas_tracer_advection_helpers.F
    iamas_tracer_advection_mono.F
    iamas_tracer_advection_std.F
    iamas_vector_operations.F
    iamas_vector_reconstruction.F)

add_library(operators ${_iamas_operators_src})

iamas_fortran_target(operators)

add_library(${PROJECT_NAME}::operators ALIAS operators)

set_target_properties(operators PROPERTIES OUTPUT_NAME iamas_operators)
target_link_libraries(operators PUBLIC ${PROJECT_NAME}::framework)

install(TARGETS operators EXPORT ${PROJECT_NAME}Exports
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
