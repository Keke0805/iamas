# 完整层级计时

说明：max(s) 为各 MPI rank 累计时间的最大值；total% = 该项 max / 总运行 max。包含子项的父项与子项不能相加；不同项的最大值可能来自不同 rank，不能相减得到精确独占时间。“待迁移”表示值得检查的数据并行计算，并非所有循环都能直接并行：垂直递推须保留层间顺序，优先并行网格柱/物种。MPI通信和主机管理单独列出，已有Kokkos的活化、沉降、辐射也作标记。


| section | 调用次数（全部rank） | max(s) | total% |
|---|---:|---:|---:|
| `total time` | 360 | 10565.3035 | 100.000% |
| &emsp;↳ `initialize` | 360 | 80.3186 | 0.760% |
| &emsp;&emsp;↳ `Stream input in init` | 360 | 54.8181 | 0.519% |
| &emsp;&emsp;&emsp;↳ `iAMAS_readStream` | 360 | 51.4709 | 0.487% |
| &emsp;&emsp;&emsp;↳ `Exch all halos in read stream` | 360 | 0.1377 | 0.001% |
| &emsp;&emsp;&emsp;↳ `Postread in read stream` | 360 | 0.0220 | 0.000% |
| &emsp;&emsp;↳ `Blocks init` | 360 | 20.7561 | 0.196% |
| &emsp;↳ `iAMAS_writeStream` | 720 | 146.7723 | 1.389% |
| &emsp;↳ `iAMAS_readStream` | 720 | 0.1964 | 0.002% |
| &emsp;↳ `Exch all halos in read stream` | 720 | 0.0015 | 0.000% |
| &emsp;↳ `Postread in read stream` | 720 | 0.0001 | 0.000% |
| &emsp;↳ `time integration` | 172800 | 8311.3049 | 78.666% |
| &emsp;&emsp;↳ `physics driver` | 172800 | 687.3390 | 6.506% |
| &emsp;&emsp;&emsp;↳ `phys_allocate` | 172800 | 0.0440 | 0.000% |
| &emsp;&emsp;&emsp;↳ `phys_state_prepare` | 172800 | 11.4475 | 0.108% |
| &emsp;&emsp;&emsp;↳ `cloudiness_driver` | 8640 | 84.1148 | 0.796% |
| &emsp;&emsp;&emsp;&emsp;↳ `calc_cldfraction` | 8640 | 83.9242 | 0.794% |
| &emsp;&emsp;&emsp;↳ `radiation_sw_driver` | 8640 | 253.6458 | 2.401% |
| &emsp;&emsp;&emsp;&emsp;↳ `RRTMG_sw` | 8640 | 250.1166 | 2.367% |
| &emsp;&emsp;&emsp;↳ `radiation_lw_driver` | 8640 | 124.8752 | 1.182% |
| &emsp;&emsp;&emsp;&emsp;↳ `RRTMG_lw` | 8640 | 123.2206 | 1.166% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `array_t` | 8640 | 0.0526 | 0.000% |
| &emsp;&emsp;&emsp;↳ `surface_layer_driver` | 172800 | 7.3764 | 0.070% |
| &emsp;&emsp;&emsp;&emsp;↳ `MYNN_sfclay` | 172800 | 5.0923 | 0.048% |
| &emsp;&emsp;&emsp;↳ `noah_land_surface` | 172800 | 5.6324 | 0.053% |
| &emsp;&emsp;&emsp;↳ `pbl_driver` | 172800 | 88.3572 | 0.836% |
| &emsp;&emsp;&emsp;&emsp;↳ `MYNN_pbl` | 172800 | 61.4733 | 0.582% |
| &emsp;&emsp;&emsp;↳ `convection_prep` | 172800 | 0.0035 | 0.000% |
| &emsp;&emsp;&emsp;↳ `convection_driver` | 172800 | 143.4797 | 1.358% |
| &emsp;&emsp;&emsp;&emsp;↳ `Grell-Freitas` | 172800 | 135.5745 | 1.283% |
| &emsp;&emsp;&emsp;↳ `convection_writeback` | 172800 | 0.0620 | 0.001% |
| &emsp;&emsp;&emsp;↳ `phys_release` | 172800 | 0.0347 | 0.000% |
| &emsp;&emsp;↳ `chemistry driver` | 172800 | 1843.4727 | 17.448% |
| &emsp;&emsp;&emsp;↳ `chemistry_from_iAMAS` | 172800 | 9.0740 | 0.086% |
| &emsp;&emsp;&emsp;↳ `optical driver` | 9000 | 794.3710 | 7.519% |
| &emsp;&emsp;&emsp;&emsp;↳ `optical_prep_sectional` | 9000 | 58.1430 | 0.550% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `opt_prep_bulk` | 9000 | 42.5609 | 0.403% |
| &emsp;&emsp;&emsp;&emsp;↳ `mieaer` | 9000 | 716.6056 | 6.783% |
| &emsp;&emsp;&emsp;↳ `photolysis driver` | 9000 | 0.0186 | 0.000% |
| &emsp;&emsp;&emsp;↳ `emission driver` | 172800 | 28.4805 | 0.270% |
| &emsp;&emsp;&emsp;&emsp;↳ `emis_natural` | 172800 | 13.5046 | 0.128% |
| &emsp;&emsp;&emsp;↳ `dry deposition` | 172800 | 774.2851 | 7.329% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_get_fields` | 172800 | 0.0390 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_init_arrays` | 172800 | 11.9335 | 0.113% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_gas_velocity` | 172800 | 1.4329 | 0.014% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_aer_resistance` | 172800 | 0.0296 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_aer_velocity` | 172800 | 9.5803 | 0.091% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_velocity_prepare` | 172800 | 0.6079 | 0.006% |
| &emsp;&emsp;&emsp;&emsp;↳ `vertmx` | 172800 | 22.2632 | 0.211% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `vertmx_workspace` | 172800 | 0.1847 | 0.002% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `vertmx_geometry` | 172800 | 1.7952 | 0.017% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `vertmx_solve_writeback` | 172800 | 20.3522 | 0.193% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_density_prepare` | 172800 | 0.4207 | 0.004% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_mixactivate` | 172800 | 513.1550 | 4.857% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_positive_clip` | 172800 | 13.6926 | 0.130% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_settling` | 172800 | 207.7632 | 1.966% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_budget_update` | 172800 | 0.0009 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;↳ `dd_scalar_writeback` | 172800 | 1.1950 | 0.011% |
| &emsp;&emsp;&emsp;↳ `chem convection` | 172800 | 267.5422 | 2.532% |
| &emsp;&emsp;&emsp;&emsp;↳ `assignment` | 172800 | 23.3751 | 0.221% |
| &emsp;&emsp;&emsp;&emsp;↳ `calculate1` | 172800 | 27.7295 | 0.262% |
| &emsp;&emsp;&emsp;&emsp;↳ `cup_ct` | 172800 | 100.9936 | 0.956% |
| &emsp;&emsp;&emsp;&emsp;↳ `calculate2` | 172800 | 11.6698 | 0.110% |
| &emsp;&emsp;&emsp;&emsp;↳ `calculate3` | 172800 | 63.1647 | 0.598% |
| &emsp;&emsp;&emsp;↳ `Kinetics thermodynamics driver` | 4680 | 0.0019 | 0.000% |
| &emsp;&emsp;&emsp;↳ `wet scavenging` | 172800 | 83.6126 | 0.791% |
| &emsp;&emsp;&emsp;&emsp;↳ `wet_state_prepare` | 172800 | 0.4276 | 0.004% |
| &emsp;&emsp;&emsp;&emsp;↳ `wet_scheme` | 172800 | 82.7595 | 0.783% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `wet_incloud` | 172800 | 0.0034 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `wet_precip_profiles` | 172800 | 2.0700 | 0.020% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `wet_aerosol_impaction` | 172800 | 18.8928 | 0.179% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `wet_gas_rain` | 172800 | 2.6359 | 0.025% |
| &emsp;&emsp;&emsp;&emsp;&emsp;↳ `wet_species_update` | 172800 | 43.4344 | 0.411% |
| &emsp;&emsp;&emsp;&emsp;↳ `wet_budget_update` | 172800 | 0.0007 | 0.000% |
| &emsp;&emsp;↳ `atm_rk_integration_setup` | 172800 | 51.9800 | 0.492% |
| &emsp;&emsp;↳ `atm_compute_moist_coefficients` | 172800 | 10.7698 | 0.102% |
| &emsp;&emsp;↳ `physics_get_tend` | 172800 | 107.4112 | 1.017% |
| &emsp;&emsp;↳ `chemistry_addtend` | 172800 | 143.4799 | 1.358% |
| &emsp;&emsp;&emsp;↳ `chem_src_get_fields` | 172800 | 0.0067 | 0.000% |
| &emsp;&emsp;&emsp;↳ `chem_src_clear` | 172800 | 15.2637 | 0.144% |
| &emsp;&emsp;&emsp;↳ `chem_src_get_emission` | 172800 | 0.0038 | 0.000% |
| &emsp;&emsp;&emsp;↳ `chem_src_update` | 172800 | 128.7701 | 1.219% |
| &emsp;&emsp;↳ `atm_compute_vert_imp_coefs` | 518400 | 12.1349 | 0.115% |
| &emsp;&emsp;↳ `atm_compute_dyn_tend` | 1555200 | 515.5767 | 4.880% |
| &emsp;&emsp;↳ `small_step_prep` | 1555200 | 32.7263 | 0.310% |
| &emsp;&emsp;↳ `atm_advance_acoustic_step` | 2073600 | 160.4536 | 1.519% |
| &emsp;&emsp;↳ `atm_divergence_damping_3d` | 2073600 | 34.3297 | 0.325% |
| &emsp;&emsp;↳ `atm_recover_large_step_variables` | 1555200 | 190.9633 | 1.807% |
| &emsp;&emsp;↳ `atm_compute_solve_diagnostics` | 1555200 | 293.6054 | 2.779% |
| &emsp;&emsp;↳ `atm_rk_dynamics_substep_finish` | 518400 | 24.5624 | 0.232% |
| &emsp;&emsp;↳ `atm_advance_for_scalars` | 345600 | 353.1606 | 3.343% |
| &emsp;&emsp;&emsp;↳ `advance_scalars_yakl` | 345600 | 353.0798 | 3.342% |
| &emsp;&emsp;↳ `atm_advance_for_chem` | 345600 | 1627.9864 | 15.409% |
| &emsp;&emsp;&emsp;↳ `advance_scalars_yakl` | 345600 | 1627.9022 | 15.408% |
| &emsp;&emsp;↳ `atm_advance_mono_for_scalars` | 172800 | 467.8636 | 4.428% |
| &emsp;&emsp;&emsp;↳ `mono_work_prep` | 172800 | 62.4081 | 0.591% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_get_field` | 172800 | 0.0023 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_source_update` | 172800 | 5.4925 | 0.052% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_source_halo` | 172800 | 55.9487 | 0.530% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_density_update` | 172800 | 2.0351 | 0.019% |
| &emsp;&emsp;&emsp;↳ `mono_precomm` | 172800 | 269.9993 | 2.556% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_init` | 172800 | 20.2530 | 0.192% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_ghost` | 172800 | 0.0142 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_vert_flux` | 172800 | 17.5295 | 0.166% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_horiz_flux` | 172800 | 76.7722 | 0.727% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_vert_upwind` | 172800 | 37.7083 | 0.357% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_horiz_lbc` | 172800 | 0.0780 | 0.001% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_div_limiter` | 172800 | 125.6270 | 1.189% |
| &emsp;&emsp;&emsp;↳ `mono_halo` | 172800 | 97.3863 | 0.922% |
| &emsp;&emsp;&emsp;↳ `mono_afcomm` | 172800 | 185.5178 | 1.756% |
| &emsp;&emsp;↳ `atm_advance_mono_for_chem` | 172800 | 1541.5949 | 14.591% |
| &emsp;&emsp;&emsp;↳ `mono_work_prep` | 172800 | 156.4489 | 1.481% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_get_field` | 172800 | 0.0023 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_source_update` | 172800 | 23.5940 | 0.223% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_source_halo` | 172800 | 132.0835 | 1.250% |
| &emsp;&emsp;&emsp;&emsp;↳ `prep_density_update` | 172800 | 2.0367 | 0.019% |
| &emsp;&emsp;&emsp;↳ `mono_precomm` | 172800 | 883.1815 | 8.359% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_init` | 172800 | 102.8972 | 0.974% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_ghost` | 172800 | 0.0528 | 0.000% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_vert_flux` | 172800 | 38.1988 | 0.362% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_horiz_flux` | 172800 | 228.2531 | 2.160% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_vert_upwind` | 172800 | 126.8981 | 1.201% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_horiz_lbc` | 172800 | 0.3414 | 0.003% |
| &emsp;&emsp;&emsp;&emsp;↳ `pc_div_limiter` | 172800 | 396.1890 | 3.750% |
| &emsp;&emsp;&emsp;↳ `mono_halo` | 172800 | 136.2345 | 1.289% |
| &emsp;&emsp;&emsp;↳ `mono_afcomm` | 172800 | 478.5480 | 4.529% |
| &emsp;&emsp;↳ `microphysics` | 172800 | 224.5350 | 2.125% |
| &emsp;&emsp;&emsp;↳ `mp_allocate` | 172800 | 2.5823 | 0.024% |
| &emsp;&emsp;&emsp;↳ `mp_state_prepare` | 172800 | 18.3648 | 0.174% |
| &emsp;&emsp;&emsp;↳ `Thompson` | 172800 | 149.7345 | 1.417% |
| &emsp;&emsp;&emsp;↳ `mp_precip_writeback` | 172800 | 6.3876 | 0.060% |
| &emsp;&emsp;&emsp;↳ `mp_state_writeback` | 172800 | 27.7301 | 0.262% |
| &emsp;&emsp;&emsp;↳ `mp_release` | 172800 | 0.0640 | 0.001% |
| &emsp;&emsp;&emsp;↳ `mp_radar_diagnostics` | 4320 | 3.4745 | 0.033% |
| &emsp;&emsp;&emsp;↳ `mp_relative_humidity` | 4320 | 0.1499 | 0.001% |
| &emsp;↳ `write output files` | 172800 | 1890.1668 | 17.890% |
| &emsp;&emsp;↳ `iAMAS_writeStream` | 8640 | 1860.5984 | 17.610% |
