
set(iAMAS_FRAMEWORK_SOURCES
    iamas_abort.F
    iamas_attlist.F
    iamas_block_creator.F
    iamas_block_decomp.F
    iamas_bootstrapping.F
    iamas_c_interfacing.F
    iamas_constants.F
    iamas_decomp.F
    iamas_derived_types.F
    iamas_dmpar.F
    iamas_domain_routines.F
    iamas_field_accessor.F
    iamas_field_routines.F
    iamas_forcing.F
    iamas_framework.F
    iamas_hash.F
    iamas_io.F
    iamas_io_streams.F
    iamas_io_units.F
    iamas_kind_types.F
    iamas_log.F
    iamas_pool_routines.F
    iamas_sort.F
    iamas_stream_list.F
    iamas_stream_manager.F
    iamas_threading.F
    iamas_timekeeping.F
    iamas_timer.F
    pool_hash.c
    random_id.c
    regex_matching.c
    xml_stream_parser.c
    c_interface_file.c
)

add_library(framework ${iAMAS_FRAMEWORK_SOURCES})
set_iAMAS_DEBUG_flag(framework)
set(FRAMEWORK_COMPILE_DEFINITIONS
        USE_PIO2
        iAMAS_PIO_SUPPORT
        iamas=1
        iAMAS_NATIVE_TIMERS)
target_compile_definitions(framework PRIVATE ${FRAMEWORK_COMPILE_DEFINITIONS})

iamas_fortran_target(framework)
add_library(${PROJECT_NAME}::framework ALIAS framework)

set_target_properties(framework PROPERTIES OUTPUT_NAME iamas_framework)

set(FRAMEWORK_LINK_LIBRARIES
        ${PROJECT_NAME}::external::esmf
        ${PROJECT_NAME}::external::ezxml
        PIO::PIO_Fortran
        PIO::PIO_C
        PnetCDF::PnetCDF_Fortran
        NetCDF::NetCDF_Fortran
        NetCDF::NetCDF_C
        MPI::MPI_Fortran)

if (iAMAS_PROFILE)
    list(APPEND FRAMEWORK_LINK_LIBRARIES GPTL::GPTL)
endif ()
target_link_libraries(framework PUBLIC ${FRAMEWORK_LINK_LIBRARIES})

install(TARGETS framework EXPORT ${PROJECT_NAME}Exports
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
