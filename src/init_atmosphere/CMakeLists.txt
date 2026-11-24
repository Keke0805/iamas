# iAMAS/src/init_atmosphere
#
# Targets
#   iAMAS::core::init_atmosphere

## Generated includes
set(init_atm_core_inc
        block_dimension_routines.inc
        core_variables.inc
        define_packages.inc
        domain_variables.inc
        namelist_call.inc
        namelist_defines.inc
        setup_immutable_streams.inc
        structs_and_variables.inc)

## core_init_atosphere
set(init_atm_core_srcs
        iamas_atm_advection.F
        iamas_atmphys_constants.F
        iamas_atmphys_date_time.F
        iamas_atmphys_functions.F
        iamas_atmphys_initialize_real.F
        iamas_atmphys_utilities.F
        iamas_init_atm_bitarray.F
        iamas_init_atm_cases.F
        iamas_init_atm_core.F
        iamas_init_atm_core_interface.F
        iamas_init_atm_gwd.F
        iamas_init_atm_hinterp.F
        iamas_init_atm_llxy.F
        iamas_init_atm_queue.F
        iamas_init_atm_read_met.F
        iamas_init_atm_static.F
        iamas_init_atm_surface.F
        iamas_init_atm_fdda.F
        iamas_init_atm_emission.F #zhxia add this for emission_update. 0414/2025
        iamas_init_atm_chem.F #zhxia add this for reading chemical field from wrfchem. 06/2025
        iamas_init_atm_vinterp.F
        iamas_geotile_manager.F #qkyang add to speed init_static
        iamas_kd_tree.F #qkyang add to speed init_static
        iamas_parse_geoindex.F #qkyang add to speed init_static
        iamas_stack.F #qkyang add to speed init_static
        AddUpStreamBuildParam.F #qkyang add to Calculate upstream building parameters
        AddSlopeShadowParam.F #qkyang add to precompute slope and shadow static fields
        read_geogrid.c)

add_library(core_init_atmosphere ${init_atm_core_srcs})
if (${DO_PHYSICS})
    target_compile_definitions(core_init_atmosphere PRIVATE DO_PHYSICS)
endif ()
if (iAMAS_DOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
else ()
    target_compile_definitions(core_init_atmosphere PRIVATE SINGLE_PRECISION)
endif ()
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    target_compile_definitions(core_init_atmosphere PRIVATE iAMAS_DEBUG)
endif ()
if (${PIO_FOUND})
    if (${USE_PIO2})
        target_compile_definitions(core_init_atmosphere PRIVATE USE_PIO2)
    else ()
        target_compile_definitions(core_init_atmosphere PRIVATE USE_PIO1)
    endif ()
    target_compile_definitions(core_init_atmosphere PRIVATE iAMAS_PIO_SUPPORT)
endif ()
target_compile_definitions(core_init_atmosphere PRIVATE iamas=1)
target_compile_definitions(framework PRIVATE iAMAS_NATIVE_TIMERS)
iamas_core_target(CORE init_atmosphere TARGET core_init_atmosphere INCLUDES ${init_atm_core_inc})
