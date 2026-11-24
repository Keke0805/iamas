# Source files physics/
set(ATMOSPHERE_CORE_PHYSICS_SOURCES
    iamas_atmphys_camrad_init.F
    iamas_atmphys_constants.F
    iamas_atmphys_control.F
    iamas_atmphys_date_time.F
    iamas_atmphys_driver.F
    iamas_atmphys_driver_cloudiness.F
    iamas_atmphys_driver_convection.F
    iamas_atmphys_driver_gwdo.F
    iamas_atmphys_driver_lsm.F
    iamas_atmphys_driver_microphysics.F
    iamas_atmphys_driver_oml.F
    iamas_atmphys_driver_pbl.F
    iamas_atmphys_driver_radiation_lw.F
    iamas_atmphys_driver_radiation_sw.F
    iamas_atmphys_driver_sfclayer.F
    iamas_atmphys_finalize.F
    iamas_atmphys_functions.F
    iamas_atmphys_init.F
    iamas_atmphys_init_microphysics.F
    # iamas_atmphys_initialize_real.F
    iamas_atmphys_interface.F
    iamas_atmphys_landuse.F
    iamas_atmphys_lsm_noahinit.F
    iamas_atmphys_manager.F
    iamas_atmphys_o3climatology.F
    iamas_atmphys_packages.F
    iamas_atmphys_rrtmg_lwinit.F
    iamas_atmphys_rrtmg_swinit.F
    iamas_atmphys_todynamics.F
    iamas_atmphys_update.F
    iamas_atmphys_update_surface.F
    iamas_atmphys_vars.F
    iamas_atmphys_utilities.F
    iamas_atmphys_fdda.F
    iamas_init_topo_wind.F
)
list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_SOURCES PREPEND physics/)

# physics/schemes/
set(ATMOSPHERE_CORE_PHYSICS_SCHEMES_SOURCES
    libmassv.F
    module_bl_gwdo.F
    module_bl_mynn.F
    module_bl_ysu.F
    module_bl_shinhong.F
    module_cam_error_function.F
    module_cam_shr_kind_mod.F
    module_cam_support.F
    module_cu_gf.iamas.F
    module_cu_kfeta.F
    module_cu_mskf.F
    module_cu_ntiedtke.F
    module_cu_tiedtke.F
    module_mp_kessler.F
    module_mp_morr_two_moment.F
    module_mp_p3.F
    module_mp_radar.F
    module_mp_thompson.F
    module_mp_thompson_cldfra3.F
    module_mp_wsm6.F
    module_tompkins_cloud.F
    module_ra_cam.F
    module_ra_cam_support.F
    module_ra_rrtmg_lw.F
    module_ra_rrtmg_sw.F
    module_ra_rrtmg_vinterp.F
    module_sf_bem.F
    module_sf_bep.F
    module_sf_bep_bem.F
    module_sf_mynn.F
    module_sf_noah_seaice.F
    module_sf_noah_seaice_drv.F
    module_sf_noahdrv.F
    module_sf_noahlsm.F
    module_sf_noahlsm_glacial_only.F
    module_sf_oml.F
    module_sf_sfcdiags.F
    module_sf_sfclay.F
    module_sf_urban.F)

list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_SCHEMES_SOURCES PREPEND physics/schemes/)

set(ATMOSPHERE_CORE_PHYSICS_MMM_SOURCES
    bl_gwdo.F90
    bl_ysu.F90
    module_bl_shinhong.F90
    cu_ntiedtke.F90
    module_libmassv.F90
    mp_wsm6.F90
    mp_wsm6_effectRad.F90
    bl_mynn.F90
    bl_mynn_subroutines.F90
    mp_radar.F90
    mynn_shared.F90
    sf_mynn.F90
    sf_sfclayrev.F90)

list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_MMM_SOURCES PREPEND physics/physics_mmm/)

set(ATMOSPHERE_CORE_PHYSICS_NOAMP_UTILITY_SOURCES
    CheckNanMod.F90 ErrorHandleMod.F90 Machine.F90)
list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_NOAMP_UTILITY_SOURCES
     PREPEND physics/physics_noahmp/utility/)

set(ATMOSPHERE_CORE_PHYSICS_NOAMP_iAMAS_DRIVER_SOURCES
    BiochemVarInTransferMod.F90
    ConfigVarOutTransferMod.F90
    ForcingVarInTransferMod.F90
    NoahmpDriverMainMod.F90
    NoahmpIOVarFinalizeMod.F90
    NoahmpReadNamelistMod.F90
    PedoTransferSR2006Mod.F90
    BiochemVarOutTransferMod.F90
    EnergyVarInTransferMod.F90
    ForcingVarOutTransferMod.F90
    NoahmpIOVarInitMod.F90
    NoahmpReadTableMod.F90
    WaterVarInTransferMod.F90
    ConfigVarInTransferMod.F90
    EnergyVarOutTransferMod.F90
    NoahmpInitMainMod.F90
    NoahmpIOVarType.F90
    NoahmpSnowInitMod.F90
    WaterVarOutTransferMod.F90)
list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_NOAMP_iAMAS_DRIVER_SOURCES
     PREPEND physics/physics_noahmp/drivers/iamas/)

set(ATMOSPHERE_CORE_PHYSICS_NOAMP_SRC_SOURCES
    AtmosForcingMod.F90
    BalanceErrorCheckGlacierMod.F90
    BalanceErrorCheckMod.F90
    BiochemCropMainMod.F90
    BiochemNatureVegMainMod.F90
    BiochemVarInitMod.F90
    BiochemVarType.F90
    CanopyHydrologyMod.F90
    CanopyRadiationTwoStreamMod.F90
    CanopyWaterInterceptMod.F90
    CarbonFluxCropMod.F90
    CarbonFluxNatureVegMod.F90
    ConfigVarInitMod.F90
    ConfigVarType.F90
    ConstantDefineMod.F90
    CropGrowDegreeDayMod.F90
    CropPhotosynthesisMod.F90
    EnergyMainGlacierMod.F90
    EnergyMainMod.F90
    EnergyVarInitMod.F90
    EnergyVarType.F90
    ForcingVarInitMod.F90
    ForcingVarType.F90
    GeneralInitGlacierMod.F90
    GeneralInitMod.F90
    GlacierIceThermalPropertyMod.F90
    GlacierPhaseChangeMod.F90
    GlacierTemperatureMainMod.F90
    GlacierTemperatureSolverMod.F90
    GlacierThermalDiffusionMod.F90
    GroundAlbedoGlacierMod.F90
    GroundAlbedoMod.F90
    GroundRoughnessPropertyGlacierMod.F90
    GroundRoughnessPropertyMod.F90
    GroundThermalPropertyGlacierMod.F90
    GroundThermalPropertyMod.F90
    GroundWaterMmfMod.F90
    GroundWaterTopModelMod.F90
    HumiditySaturationMod.F90
    IrrigationFloodMod.F90
    IrrigationInfilPhilipMod.F90
    IrrigationMicroMod.F90
    IrrigationPrepareMod.F90
    IrrigationSprinklerMod.F90
    IrrigationTriggerMod.F90
    Makefile
    MatrixSolverTriDiagonalMod.F90
    NoahmpMainGlacierMod.F90
    NoahmpMainMod.F90
    NoahmpVarType.F90
    PhenologyMainMod.F90
    PrecipitationHeatAdvectGlacierMod.F90
    PrecipitationHeatAdvectMod.F90
    PsychrometricVariableGlacierMod.F90
    PsychrometricVariableMod.F90
    ResistanceAboveCanopyChen97Mod.F90
    ResistanceAboveCanopyMostMod.F90
    ResistanceBareGroundChen97Mod.F90
    ResistanceBareGroundMostMod.F90
    ResistanceCanopyStomataBallBerryMod.F90
    ResistanceCanopyStomataJarvisMod.F90
    ResistanceGroundEvaporationGlacierMod.F90
    ResistanceGroundEvaporationMod.F90
    ResistanceLeafToGroundMod.F90
    RunoffSubSurfaceDrainageMod.F90
    RunoffSubSurfaceEquiWaterTableMod.F90
    RunoffSubSurfaceGroundWaterMod.F90
    RunoffSubSurfaceShallowMmfMod.F90
    RunoffSurfaceBatsMod.F90
    RunoffSurfaceDynamicVicMod.F90
    RunoffSurfaceExcessDynamicVicMod.F90
    RunoffSurfaceFreeDrainMod.F90
    RunoffSurfaceTopModelEquiMod.F90
    RunoffSurfaceTopModelGrdMod.F90
    RunoffSurfaceTopModelMmfMod.F90
    RunoffSurfaceVicMod.F90
    RunoffSurfaceXinAnJiangMod.F90
    ShallowWaterTableMmfMod.F90
    SnowAgingBatsMod.F90
    SnowAlbedoBatsMod.F90
    SnowAlbedoClassMod.F90
    SnowCoverGlacierMod.F90
    SnowCoverGroundNiu07Mod.F90
    SnowfallBelowCanopyMod.F90
    SnowLayerCombineMod.F90
    SnowLayerDivideMod.F90
    SnowLayerWaterComboMod.F90
    SnowpackCompactionMod.F90
    SnowpackHydrologyGlacierMod.F90
    SnowpackHydrologyMod.F90
    SnowThermalPropertyMod.F90
    SnowWaterMainGlacierMod.F90
    SnowWaterMainMod.F90
    SoilHydraulicPropertyMod.F90
    SoilMoistureSolverMod.F90
    SoilSnowTemperatureMainMod.F90
    SoilSnowTemperatureSolverMod.F90
    SoilSnowThermalDiffusionMod.F90
    SoilSnowWaterPhaseChangeMod.F90
    SoilThermalPropertyMod.F90
    SoilWaterDiffusionRichardsMod.F90
    SoilWaterInfilGreenAmptMod.F90
    SoilWaterInfilPhilipMod.F90
    SoilWaterInfilSmithParlangeMod.F90
    SoilWaterMainMod.F90
    SoilWaterSupercoolKoren99Mod.F90
    SoilWaterSupercoolNiu06Mod.F90
    SoilWaterTranspirationMod.F90
    SurfaceAlbedoGlacierMod.F90
    SurfaceAlbedoMod.F90
    SurfaceEmissivityGlacierMod.F90
    SurfaceEmissivityMod.F90
    SurfaceEnergyFluxBareGroundMod.F90
    SurfaceEnergyFluxGlacierMod.F90
    SurfaceEnergyFluxVegetatedMod.F90
    SurfaceRadiationGlacierMod.F90
    SurfaceRadiationMod.F90
    TileDrainageEquiDepthMod.F90
    TileDrainageHooghoudtMod.F90
    TileDrainageSimpleMod.F90
    VaporPressureSaturationMod.F90
    WaterMainGlacierMod.F90
    WaterMainMod.F90
    WaterTableDepthSearchMod.F90
    WaterTableEquilibriumMod.F90
    WaterVarInitMod.F90
    WaterVarType.F90)
list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_NOAMP_SRC_SOURCES
     PREPEND physics/physics_noahmp/src/)

# diagnostics/
set(ATMOSPHERE_CORE_DIAGNOSTIC_SOURCES
    iamas_convective_diagnostics.F
    iamas_isobaric_diagnostics.F
    iamas_atm_diagnostic_template.F
    iamas_atm_diagnostics_manager.F
    iamas_atm_diagnostics_utils.F
    iamas_pv_diagnostics.F
    iamas_soundings.F
    iamas_aerosols_diagnostics.F
    iamas_column_diagnostics.F
    iamas_isohypsic_diagnostics.F
    iamas_cldfracmax_diagnostics.F)

list(TRANSFORM ATMOSPHERE_CORE_DIAGNOSTIC_SOURCES PREPEND diagnostics/)

# dynamics/
set(ATMOSPHERE_CORE_DYNAMICS_SOURCES iamas_atm_boundaries.F iamas_atm_iau.F
                                     iamas_les_init.F
                                     iamas_les_sgs.F
                                     iamas_atm_time_integration.F)
list(TRANSFORM ATMOSPHERE_CORE_DYNAMICS_SOURCES PREPEND dynamics/)

# utils/
set(ATMOSPHERE_CORE_UTILS_SOURCES atmphys_build_tables_thompson.F
                                  build_tables.F)
list(TRANSFORM ATMOSPHERE_CORE_UTILS_SOURCES PREPEND utils/)

# core_atosphere
set(ATMOSPHERE_CORE_SOURCES
    iamas_atm_dimensions.F iamas_atm_threading.F iamas_atm_core.F
    iamas_atm_core_interface.F
    # iamas_atm_halos.F
)

# Generated includes
set(ATMOSPHERE_CORE_INCLUDES
    block_dimension_routines.inc
    core_variables.inc
    define_packages.inc
    domain_variables.inc
    namelist_call.inc
    namelist_defines.inc
    setup_immutable_streams.inc
    structs_and_variables.inc)

# Chemistry
set(ATMOSPHERE_CORE_CHEMISTRY_SOURCES
    iamas_atmchem_binaer_init.F
    iamas_atmchem_gas_init.F
    iamas_atmchem_configs.F
    iamas_atmchem_constants.F
    iamas_atmchem_control.F
    iamas_atmchem_driver.F
    iamas_atmchem_driver_convection.F
    iamas_atmchem_driver_drydep.F
    iamas_atmchem_driver_optical.F
    iamas_atmchem_driver_wetscav.F
    iamas_atmchem_gocart_dust_init.F
    iamas_atmchem_interface.F
    iamas_atmchem_packages.F
    iamas_atmchem_todynamics.F
    iamas_atmchem_vars.F
    iamas_atmchem_utilities.F
    iamas_atmchem_driver_emission.F
    iamas_atmchem_driver_photolysis.F
    iamas_atmchem_driver_kinetics_thermodynamics.F)
list(TRANSFORM ATMOSPHERE_CORE_CHEMISTRY_SOURCES PREPEND chemistry/)

set(ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_SOURCES
    module_aer_drydep.F
    module_cal_cldfra.F
    module_data_cmu_bulkaqchem.F
    module_cmu_dvode_solver.F
    module_cmu_bulkaqchem.F
    module_ctrans_aqchem.F
    module_ctrans_general.F
    module_ctrans_kfeta.F
    module_cu_gd.F
    module_data_binaer_asect.F
    module_data_binaer_other.F
    module_data_gocart_dust.F
    module_data_mosaic_therm.F
    module_data_radm2.F
    module_data_rrtmgaeropt.F
    module_data_sorgam.F
    module_dep_simple.F
    module_gocart_drydep.F
    module_mixactivate.F
    module_mixactivate_wrappers.F
    module_mosaic_addemiss.F
    # module_mosaic_drydep.F
    module_mosaic_settling.F
    module_mosaic_wetscav.F
    module_optical_averaging.F
    module_phot_fastj.F
    module_fastj_mie.F
    module_fastj_data.F
    module_data_cbmz.F
    module_vertmx_wrf.F
    module_linoz.F
	module_sia_bins_simple.F
	module_soa_bins_simple.F
	module_mosaic_driver.F
	module_cmu_cloudchem.F
	module_mosaic_therm.F
	module_mosaic_csuesat.F
	module_mosaic_newnuc.F
	module_mosaic_coag.F
	module_mosaic_movesect.F)
if(USE_DNN)
  list(APPEND ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_SOURCES module_AIMACI_mosaic.F
       module_AIPC_saprc99.F)
endif()
if(USE_KPP)
  list(APPEND ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_SOURCES module_KPP_saprc99.F)
endif()

list(TRANSFORM ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_SOURCES
     PREPEND chemistry/schemes/)
# ============================================================================
# KPP Compilation and Code Generation Logic with Pure Configure-Time Execution
# Author: Zihan Xia, xiazh@mail.ustc.edu.cn
# Date: 2025/10/29
# ============================================================================
if(USE_KPP)
  message(STATUS "*** Enable KPP - symbolic chemistry Kinetics PreProcessor, Version 3.3.0 ***")

  # --- Part 1: Compile the KPP executable ---
  set(KPP_MAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/chemistry/schemes/KPP/kpp-3.3.0/src)
  set(KPP_EXECUTABLE ${CMAKE_BINARY_DIR}/bin/kpp)
  
  # define compile command
  set(KPP_MAKE_COMMAND
      ${CMAKE_MAKE_PROGRAM}
      "CC=${CMAKE_C_COMPILER}"
      "CC_FLAGS=${CMAKE_C_FLAGS}"
      "FLEX=${FLEX_EXECUTABLE}"
      "BISON=${BISON_EXECUTABLE}"
  )
  
  message(STATUS "Building KPP tool...")
  execute_process(
    COMMAND ${KPP_MAKE_COMMAND}
    WORKING_DIRECTORY ${KPP_MAKE_DIR}
    RESULT_VARIABLE KPP_BUILD_RESULT
  )
  if(NOT KPP_BUILD_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to build KPP tool.")
  endif()
  
  # copy kpp exe to specific dir
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy
            ${KPP_MAKE_DIR}/../bin/kpp
            ${KPP_EXECUTABLE}
  )

  # --- Part 2: Run KPP to generate mechanisms codes ---
  message(STATUS "Running KPP to generate mechanism files...")
  set(KPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/chemistry/schemes/KPP/kpp-3.3.0)
  set(KPP_MECHANISMS_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/chemistry/schemes/KPP/mechanisms)
  file(GLOB MECHANISM_DIRS LIST_DIRECTORIES true ${KPP_MECHANISMS_PARENT_DIR}/*)
  message(STATUS "Found KPP mechanisms: ${MECHANISM_DIRS}")

  set(ALL_KPP_GENERATED_SOURCES "")

  foreach(_mechanism_dir IN LISTS MECHANISM_DIRS)
    get_filename_component(_mechanism_name ${_mechanism_dir} NAME)
    message(STATUS "Configuring KPP for mechanism: ${_mechanism_name}")

    execute_process(
      COMMAND ${CMAKE_COMMAND} -E env "KPP_HOME=${KPP_DIR}"
              ${KPP_EXECUTABLE} ${_mechanism_name}.kpp
      WORKING_DIRECTORY ${_mechanism_dir}
      RESULT_VARIABLE KPP_RUN_RESULT
    )
    if(NOT KPP_RUN_RESULT EQUAL 0)
      message(FATAL_ERROR "KPP failed for mechanism '${_mechanism_name}'.")
    endif()

    file(GLOB MECHANISM_SOURCES
        "${_mechanism_dir}/${_mechanism_name}_*.f90"
    )
    list(FILTER MECHANISM_SOURCES EXCLUDE REGEX ".*${_mechanism_name}_(Initialize|Main)\\.f90$")
    message(STATUS "Mechanism '${_mechanism_name}' generated sources: ${MECHANISM_SOURCES}")

    list(APPEND ALL_KPP_GENERATED_SOURCES ${MECHANISM_SOURCES})
  endforeach()

  list(APPEND ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_SOURCES ${ALL_KPP_GENERATED_SOURCES})
  message(STATUS "Successfully added all KPP-generated sources to the build.")
endif()
# =================================================================

# ============================================================================
# Kokkos Compilation
# Author: Xin Qu, xinqu@mail.ustc.edu.cn
# Date: 2026/01/26
# ============================================================================
enable_language(CXX)
set(ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_KOKKOS_SOURCES
    module_optical_averaging_flcl.cpp module_mosaic_settling_flcl.cpp module_ctrans_general_flcl.cpp module_mixactivate_flcl.cpp)
list(TRANSFORM ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_KOKKOS_SOURCES
     PREPEND chemistry/schemes/)
add_library(flcl_chemistry STATIC
            ${ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_KOKKOS_SOURCES})
target_link_libraries(flcl_chemistry PUBLIC flcl::flcl)

set_target_properties(flcl_chemistry PROPERTIES POSITION_INDEPENDENT_CODE ON
                                                LINK_LANGUAGE CXX )
if(PLATFORM_FLAG MATCHES "sunway@swustc" OR PLATFORM_FLAG MATCHES "sunway@swustciamas" OR PLATFORM_FLAG MATCHES "sunway@swiamas")
    set_target_properties(flcl_chemistry PROPERTIES COMPILE_OPTIONS "-mhybrid-coding;-fpermissive")
endif()    
set(ATMOSPHERE_CORE_PHYSICS_SCHEMES_KOKKOS_SOURCES
    module_ra_rrtmg_lw_flcl.cpp 
    module_ra_rrtmg_sw_flcl.cpp)
list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_SCHEMES_KOKKOS_SOURCES
     PREPEND physics/schemes/)
if(PLATFORM_FLAG MATCHES "sunway@swustc" OR PLATFORM_FLAG MATCHES "sunway@swustciamas" OR PLATFORM_FLAG MATCHES "sunway@swiamas")
	set_source_files_properties(${ATMOSPHERE_CORE_PHYSICS_SCHEMES_KOKKOS_SOURCES} PROPERTIES COMPILE_OPTIONS "-mhybrid-coding;-fpermissive")
endif()
add_library(flcl_physics STATIC
    ${ATMOSPHERE_CORE_PHYSICS_SCHEMES_KOKKOS_SOURCES})
target_link_libraries(flcl_physics PUBLIC Kokkos::kokkos)
target_link_libraries(flcl_physics PUBLIC flcl::flcl)
set_target_properties(flcl_physics PROPERTIES POSITION_INDEPENDENT_CODE ON
                                                LINK_LANGUAGE CXX )

set(ATMOSPHERE_CORE_PHYSICS_KOKKOS_SOURCES
    iamas_atmphys_todynamics.cpp)
list(TRANSFORM ATMOSPHERE_CORE_PHYSICS_KOKKOS_SOURCES
     PREPEND physics/)
if(PLATFORM_FLAG MATCHES "sunway@swustc" OR PLATFORM_FLAG MATCHES "sunway@swustciamas" OR PLATFORM_FLAG MATCHES "sunway@swiamas")
	set_source_files_properties(${ATMOSPHERE_CORE_PHYSICS_KOKKOS_SOURCES} PROPERTIES COMPILE_OPTIONS "-mhybrid-coding;-fpermissive")
endif()
add_library(flcl_physics_tend STATIC
    ${ATMOSPHERE_CORE_PHYSICS_KOKKOS_SOURCES})
target_link_libraries(flcl_physics_tend PUBLIC Kokkos::kokkos)
target_link_libraries(flcl_physics_tend PUBLIC flcl::flcl)
set_target_properties(flcl_physics_tend PROPERTIES POSITION_INDEPENDENT_CODE ON
                                                LINK_LANGUAGE CXX )

set(ATMOSPHERE_CORE_DYNAMICS_KOKKOS_SOURCES
    iamas_atm_time_integration_scalars.cpp
    iamas_atm_time_integration_scalars_mono_precomm.cpp
    iamas_atm_time_integration_scalars_mono_afcomm.cpp
    iamas_atm_time_integration_large_step_variables.cpp
    iamas_atm_time_integration_solve_diagnostics.cpp
    iamas_atm_time_integration_dyn_tend.cpp
    iamas_atm_time_integration_acoustic_step.cpp)
list(TRANSFORM ATMOSPHERE_CORE_DYNAMICS_KOKKOS_SOURCES PREPEND dynamics/)
add_library(flcl_dynamics STATIC ${ATMOSPHERE_CORE_DYNAMICS_KOKKOS_SOURCES})
target_link_libraries(flcl_dynamics PUBLIC flcl::flcl)

set_target_properties(flcl_dynamics PROPERTIES POSITION_INDEPENDENT_CODE ON
                                                LINK_LANGUAGE CXX )
if(PLATFORM_FLAG MATCHES "sunway@swustc" OR PLATFORM_FLAG MATCHES "sunway@swustciamas" OR PLATFORM_FLAG MATCHES "sunway@swiamas")
    set_target_properties(flcl_dynamics PROPERTIES COMPILE_OPTIONS "-mhybrid-coding;-fpermissive")
endif()


if(flcl_FOUND)				       
	install(TARGETS flcl_chemistry flcl_dynamics flcl_physics flcl_physics_tend
        	EXPORT ${PROJECT_NAME}ExportsCore)
else()
	install(TARGETS flcl_chemistry flcl_dynamics flcl_physics flcl_physics_tend flcl
        	EXPORT ${PROJECT_NAME}ExportsCore)
endif()

# ============================================================================
# AI chemical schemes compile.
# Author: Zihan Xia, xiazh@mail.ustc.edu.cn
# ============================================================================
if(USE_DNN)
  set(AICHEM_SOURCES AIMACImosaic.cpp AIPCsaprc99.cpp)
  list(TRANSFORM AICHEM_SOURCES PREPEND chemistry/schemes/)
  add_library(aichem STATIC ${AICHEM_SOURCES})
  target_link_libraries(aichem PUBLIC torch)
  set_target_properties(aichem PROPERTIES POSITION_INDEPENDENT_CODE ON)
  install(TARGETS aichem EXPORT ${PROJECT_NAME}ExportsCore)
endif()
# ============================================================================

add_library(
  core_atmosphere
  ${ATMOSPHERE_CORE_SOURCES}
  # ${ATMOSPHERE_CORE_PHYSICS_NOAMP_UTILITY_SOURCES}
  # ${ATMOSPHERE_CORE_PHYSICS_NOAMP_iAMAS_DRIVER_SOURCES}
  # ${ATMOSPHERE_CORE_PHYSICS_NOAMP_SRC_SOURCES}
  ${ATMOSPHERE_CORE_PHYSICS_SOURCES}
  # ${ATMOSPHERE_CORE_PHYSICS_MMM_SOURCES}
  ${ATMOSPHERE_CORE_PHYSICS_SCHEMES_SOURCES}
  ${ATMOSPHERE_CORE_DIAGNOSTIC_SOURCES}
  ${ATMOSPHERE_CORE_DYNAMICS_SOURCES}
  ${ATMOSPHERE_CORE_CHEMISTRY_SOURCES}
  ${ATMOSPHERE_CORE_CHEMISTRY_SCHEMES_SOURCES})

set(CORE_ATMOSPHERE_COMPILE_DEFINITIONS iamas=1 iAMAS_NATIVE_TIMERS)
if(${DO_PHYSICS})
  list(APPEND CORE_ATMOSPHERE_COMPILE_DEFINITIONS DO_PHYSICS)
  list(APPEND ATMOSPHERE_COMPILE_DEFINITIONS DO_PHYSICS)
endif()
if(${DO_CHEMISTRY})
  list(APPEND CORE_ATMOSPHERE_COMPILE_DEFINITIONS DO_CHEMISTRY)
  list(APPEND ATMOSPHERE_COMPILE_DEFINITIONS DO_CHEMISTRY)
endif()
if(USE_DEBUG)
  list(APPEND CORE_ATMOSPHERE_COMPILE_DEFINITIONS USE_DEBUG)
endif()
if(USE_DNN)
  list(APPEND CORE_ATMOSPHERE_COMPILE_DEFINITIONS USE_DNN)
endif()
if(USE_KPP)
  list(APPEND CORE_ATMOSPHERE_COMPILE_DEFINITIONS USE_KPP)
endif()
if(${USE_SW})
    list(APPEND ATMOSPHERE_COMPILE_DEFINITIONS USE_SW)
endif()
if(${USE_SERIALIZATION})
    list(APPEND ATMOSPHERE_COMPILE_DEFINITIONS USE_SERIALIZATION)
endif()
target_compile_definitions(core_atmosphere
                           PRIVATE ${CORE_ATMOSPHERE_COMPILE_DEFINITIONS})
target_compile_definitions(flcl_physics
                            PUBLIC ${ATMOSPHERE_COMPILE_DEFINITIONS})
target_compile_definitions(flcl_chemistry
                            PUBLIC ${ATMOSPHERE_COMPILE_DEFINITIONS})
target_link_libraries(core_atmosphere PUBLIC flcl_chemistry flcl_dynamics flcl_physics flcl_physics_tend
                                             $<$<BOOL:${USE_DNN}>:aichem>)
set_iamas_debug_flag(core_atmosphere)
iamas_core_target(CORE atmosphere TARGET core_atmosphere INCLUDES
                 ${ATMOSPHERE_CORE_INCLUDES})
set_target_properties(atmosphere PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(atmosphere_build_tables ${ATMOSPHERE_CORE_UTILS_SOURCES})
set_target_properties(atmosphere_build_tables PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(atmosphere_build_tables PUBLIC core_atmosphere)
iamas_fortran_target(atmosphere_build_tables)
install(TARGETS atmosphere_build_tables EXPORT ${PROJECT_NAME}ExportsCore)
