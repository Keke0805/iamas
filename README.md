# iAMAS_v2.6.3

- Added support for source-area building turbulence parameterization
- Improvement of the Chengyu model
- Add support for the new aerosol dry deposition scheme (zhang et al.2001). We can set config_aer_drydep_opt=301 (or 302, or 311, or 312) to use this scheme.
- Fix a bug when we use shao scheme to compute dust emission.
- Add the TVD-FFSL hybrid scheme for chemical tracers transport
- Fix the bug when chem_opt is set to a scheme that does not contain "_aq". For non-_aq schemes, we use "vertmx" for aerosols and gases; for _aq schemes, we use "vertmx" for gases and "mixactivate" for aerosols.
- Replace seas with explicit Na/Cl and update related processes for Na, Cl, Ca, and CO3
- Add selectable dust source-function support for EROD and SSM
- Fix the bug that fail to generate precomupted coefficients for TVD-FFSL transport scheme
- ctrans_general mixactivate mosaic_settling OPT
- Fixed the issue of calculating Shortwave radiation heating rate in the shortwave radiation scheme.

# iAMAS_v2.6.2

- Integrate simple_SIA scheme for inorganic aerosol formation.
- Integrate simple_SOA scheme for organic aerosol formation.
- Add Kokkos code for the acoustics_step & physics_tend.
- Complete integration of mechanism-based MOSAIC scheme.
- Add new Kokkos code in Morrison Scheme.
- Improving Boundary Handling for Stable Long-term Regional Chemistry Simulations
- Fix LBC humidity vertical interpolation bug for top-to-bottom levels(following MPAS v8.3.1)
- Modify vertical velocity LBC from zero-gradient to constant w=0

# iAMAS_v2.6.1

- Add P3 from WRFv4.6.
- Add FDDA nudging support for theta, qv, u, and v.
- Fixed a critical bug that prevented simulation when the height of the layer exceeded 70
- Add Kokkos code for the dyn_tend,large_step_variables,solve_diagnostics subroutines

# iAMAS_v2.5.1

- Integrate the entire KPP (Kinetic PreProcessor) framework within iAMAS, supporting both general platforms and the Sunway platform (with contributions from Xinqu).
- Adapting the dynamic core for deep-atmosphere.
- Fix potential errors in filename parsing when using multiple files for IO.
- Use C for reading graph.info.part.* file, significantly improving read efficiency for high-resolution cases.
- Add support for the surface drag parameterization (config_topo_wind).
- Add rrtmg_sw/rrtmg_lw code in Kokkos.
- Optimize the computation speed of generating emission_update.nc using fine2coarse method.
- Add the Linoz module in atmosphere_model for the simulation of stratospheric ozone chemistry.

# iAMAS_v2.0.3

- Introduce the Fragment-PIO feature. Besides, add new function which could output region-of-interest data rather than global data.
- The pressure-level height diagnostic is now interpolated using ln(p) for increased accuracy. This significantly reduces the diagnostics height bias in the previous version.
- Reset the default for the lower air-temperature extrapolaton (config_extrap_airtemp) from "linear" to "lapse-rate" in the namelist.init_atmosphere. Besides, no extrapolation of the horizontal velocity in initilization.
- Add implementation of Local Climate Zone (LCZ) landuse classes using WUDAPT (51-61)
- Optimize array on driver_convcetion and driver_drydep
- Refactor emission module in iAMAS. From a technical implementation perspective, emissions are currently categorized into two types: (1) Initial emissions (e.g., EROD, DMS, AH) that are read only at the initial timestep, with subsequent calculations performed online. These are stored in init_final.nc. (2) Time-updated emissions (e.g., anthropogenic and biomass burning emissions), which are updated at specified intervals and stored in emission_update.nc.
- Add urban model: Multi - layer urban canopy (BEP) and indoor–outdoor exchange (BEM) models. config_urban_scheme = 0, 1, 2, 3 is now supported to turn on Non - urban, Urban, BEP, and BEP_BEM.
- Add support for parallel calculation to generate init_static.nc using graph.info.part.xxx
- Add support for the YSU boundary layer scheme in the BEP and BEP_BEM urban models.
- Add the ability to read simulation results (only gas species now) from wrfchem as initial values in iAMAS.
- Add a new diagnostic framework for aerosols. Dust/seas column density is calculated.
- Add support for compilation on the H800 computing platform by setting the -P parameter of the compilation command to "h800".
- The code structure of the long-wave radiation module is optimized so that the compiler can automatically vectorize part loops, and the module is accelerated by about 30% on the Intel platform.
- Add support for the Shin-Hong boundary layer scheme.
- Replace YAKL support with Kokkos support. (require: glibc >= 2.25)
- Update the optical, dynamics, and mosaic modules to Kokkos-based parallel code.

# iAMAS_v2.0.2

- Limit the maximum of dt using Morrison scheme to 90s, same as Thompson.
- Now use config_aer_cloud_droplet_on and config_aer_cloud_ice_on to turn on/off aerosol-cloud feedbacks. Both droplet and ice feedbacks are available for Thompson. Only droplet feedback is available for Morrison for now. I'm working on the cloud ice feedbacks.
- When config_aer_cloud_droplet_on=false, qndrop is set to cloud droplet number (similar to nc) to diagnose nc calculated from predefined Nt_c
- Fix a bug that causing qndrop calculated from mixactivate is not going back into memory in drydep driver. Causing qndrop values are 0. while sending to microphysics drivers.
- Fix issue #21 that causing runtime error in Debug Mode. Checking multiple gas indexes while accessing gas related
- config_urban_scheme=1 is now supported to turn on the urban model (only values 0 and 1 are currently allowed). Notes for Adding urban model: (1) Three parameters—coszr, declin, and hrang—are obtained from driver_radiation_sw. (2) Array slicing acceleration has been completed in the Land Surface Model. (3) The consistency of simulation results before and after array slicing has been verified, and debug logs, including time detection and parameter outputs, have been removed.
- Fix a bug that result in wrong mminlu in module_dep_simple.F.

# iAMAS_v2.0.1

- Fix bug where rrtmg_sw/rrtmg_lw was unnecessarily called multiple times to diagnose aerosol radiation fluxes when config_aer_ra_feedback=false.
- Fix a bug when we compile iAMAS in sunway@iamas node.
- Close the log information of running yakl.
- Fix the path error of the build_cmake script.
- Fix bug in run_atmosphere.sh to avoid errors while generating init_final.nc.
- Fix compilation failure when -B option was specified
- Modify file struct and delete core_landice, core_ocean, core_seaice, core_sw, core_test
- Change the name of the executable file to "atmosphere" and "init_atmosphere"
- YAKL optimization options upgraded to O3
- Adjust the log information when we choose mismatched schemes.
- Move the YAKL compiled files to build/build_yakl
- Add "config_input_co2_file" and "config_anthro_emis_file" in namelist configs to determine the file paths of input co2 file and anthro emis file.
- Add information in LICIENSE and version information in README.md

# iAMAS_v2.0.0

- Add stable support for CMake build system
- Integrate YAKL (Yet Another Kernel Library) framework for iAMAS, for swHPC, GPU application

# iAMAS_v1.7.0

- Introduce the framework for simulating aerosol chemistry with AI scheme

# iAMAS_v1.6.0

- Introduce the framework for simulating gas chemistry with AI scheme

# iAMAS_v1.4.0

- Introduce the framework for simulating biomass and anthropogenic species.

# iAMAS_v1.0.0

- Introduce the framework for simulating aerosol-cloud interaction

# iAMAS_v0.7.1

- Introduce the framework for simulating aerosol-radiation interaction

# iAMAS_v0.6.0

- Introduce the framework for simulating full-size range of dust from 0-40um

# iAMAS_v0.1.0

- Introduce the chemistry framework within the core_atmosphere module.
- Added foundational features to support chemistry simulations, including the coupling of chemical tendencies with dynamics.
- Establish a baseline for future development, enhancements, and extensions.

# iAMAS_v0.0.0

Migrate from MPAS-Atmosphere v7.0



The integrated Atmospheric Model Across Scales (iAMAS) is a model developed for simulating all the components related to atmospheric processes across scales, including troposphere, stratosphere, mesosphere, ionosphere, for research in climate, weather, and air quality. The similar framework is also extended to other planets than the Earth, such as Mars.



For information about building and running the model, please refer to the wiki page: https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/wikis/home
