# 12小时Kokkos迁移候选与时间占比

来源：job **9737442**；360 MPE；480步×90秒；模型检查通过，Error/Critical=0/0。

总运行 max **10565.3035秒**；时间积分 max **8311.3049秒**。下表使用本次长测数据，不引用6分钟短测比例。

当前构建仍为Kokkos Serial/MPE。“可迁移”表示源码中存在值得封装的数据并行计算，不表示只加parallel_for就能直接使用CPE；迁移后还需验证后端、数值结果和实际收益。

## 优先检查的计算阶段（按本次占比排序）

表中不把已有子计时的候选父项重复列入；“优先”仅按测得耗时排序，实施次序还要考虑工作量与正确性验证。占比不等于可节省比例，不给出相加后的预计加速。

| 过程 / 细分阶段 | max(s) | 占总运行 | 占时间积分 | 大概改法 |
|---|---:|---:|---:|---|
| `microphysics:Thompson`<br>Thompson | 149.7345 | 1.417% | 1.802% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `physics driver:convection_driver:Grell-Freitas`<br>Grell-Freitas | 135.5745 | 1.283% | 1.631% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `chemistry_addtend:chem_src_update`<br>排放源项累加 | 128.7701 | 1.219% | 1.549% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry driver:chem convection:cup_ct`<br>cup_ct | 100.9936 | 0.956% | 1.215% | 按柱隔离工作区；处理Fortran回调，保留柱内递推 |
| `physics driver:cloudiness_driver:calc_cldfraction`<br>calc_cldfraction | 83.9242 | 0.794% | 1.010% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:chem convection:calculate3`<br>calculate3 | 63.1647 | 0.598% | 0.760% | 先核对写入范围和阶段依赖，再并行独立单元/层 |
| `physics driver:pbl_driver:MYNN_pbl`<br>MYNN_pbl | 61.4733 | 0.582% | 0.740% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `atm_rk_integration_setup`<br>atm_rk_integration_setup | 51.9800 | 0.492% | 0.625% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:wet scavenging:wet_scheme:wet_species_update`<br>各物种浓度更新 | 43.4344 | 0.411% | 0.523% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry driver:optical driver:optical_prep_sectional:opt_prep_bulk`<br>光学预处理层/单元/粒径段/波段主循环 | 42.5609 | 0.403% | 0.512% | 并行网格柱/独立层；粒径和波段求和保留原顺序 |
| `atm_divergence_damping_3d`<br>atm_divergence_damping_3d | 34.3297 | 0.325% | 0.413% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `small_step_prep`<br>small_step_prep | 32.7263 | 0.310% | 0.394% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `microphysics:mp_state_writeback`<br>微物理状态和倾向回写 | 27.7301 | 0.262% | 0.334% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry driver:chem convection:calculate1`<br>calculate1 | 27.7295 | 0.262% | 0.334% | 先核对写入范围和阶段依赖，再并行独立单元/层 |
| `atm_rk_dynamics_substep_finish`<br>atm_rk_dynamics_substep_finish | 24.5624 | 0.232% | 0.296% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `atm_advance_mono_for_chem:mono_work_prep:prep_source_update`<br>加入源项并清零倾向 | 23.5940 | 0.223% | 0.284% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry driver:chem convection:assignment`<br>assignment | 23.3751 | 0.221% | 0.281% | 先核对写入范围和阶段依赖，再并行独立单元/层 |
| `chemistry driver:wet scavenging:wet_scheme:wet_aerosol_impaction`<br>雨下气溶胶碰并清除 | 18.8928 | 0.179% | 0.227% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `microphysics:mp_state_prepare`<br>微物理输入准备 | 18.3648 | 0.174% | 0.221% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry_addtend:chem_src_clear`<br>化学倾向清零 | 15.2637 | 0.144% | 0.184% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry driver:dry deposition:dd_positive_clip`<br>浓度非负限幅 | 13.6926 | 0.130% | 0.165% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry driver:emission driver:emis_natural`<br>自然源排放 | 13.5046 | 0.128% | 0.162% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `atm_compute_vert_imp_coefs`<br>atm_compute_vert_imp_coefs | 12.1349 | 0.115% | 0.146% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:dry deposition:dd_init_arrays`<br>数组初始化、预算快照 | 11.9335 | 0.113% | 0.144% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:chem convection:calculate2`<br>calculate2 | 11.6698 | 0.110% | 0.140% | 先核对写入范围和阶段依赖，再并行独立单元/层 |
| `physics driver:phys_state_prepare`<br>物理状态准备和转换 | 11.4475 | 0.108% | 0.138% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `atm_compute_moist_coefficients`<br>atm_compute_moist_coefficients | 10.7698 | 0.102% | 0.130% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:dry deposition:dd_aer_velocity`<br>气溶胶沉降速度 | 9.5803 | 0.091% | 0.115% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `chemistry driver:chemistry_from_iAMAS`<br>chemistry_from_iAMAS | 9.0740 | 0.086% | 0.109% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `microphysics:mp_precip_writeback`<br>降水及化学清除诊断回写 | 6.3876 | 0.060% | 0.077% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `physics driver:noah_land_surface`<br>Noah陆面计算接口 | 5.6324 | 0.053% | 0.068% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `atm_advance_mono_for_scalars:mono_work_prep:prep_source_update`<br>加入源项并清零倾向 | 5.4925 | 0.052% | 0.066% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `physics driver:surface_layer_driver:MYNN_sfclay`<br>MYNN_sfclay | 5.0923 | 0.048% | 0.061% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `microphysics:mp_radar_diagnostics`<br>可选雷达反射率 | 3.4745 | 0.033% | 0.042% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:wet scavenging:wet_scheme:wet_gas_rain`<br>气体雨洗 | 2.6359 | 0.025% | 0.032% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:wet scavenging:wet_scheme:wet_precip_profiles`<br>降水和云量廓线、垂直累积 | 2.0700 | 0.020% | 0.025% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `chemistry driver:dry deposition:dd_gas_velocity`<br>气体沉降速度（Wesely） | 1.4329 | 0.014% | 0.017% | 按网格柱并行；保留柱内求解顺序，工作数组私有化 |
| `chemistry driver:dry deposition:dd_scalar_writeback`<br>结果写回标量数组 | 1.1950 | 0.011% | 0.014% | 按位置、垂直层或物种拆分独立更新；保留物种映射和边界 |
| `chemistry driver:dry deposition:dd_velocity_prepare`<br>沉降速度整理与限幅 | 0.6079 | 0.006% | 0.007% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:dry deposition:dd_density_prepare`<br>密度及单位转换 | 0.4207 | 0.004% | 0.005% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `microphysics:mp_relative_humidity`<br>相对湿度诊断 | 0.1499 | 0.001% | 0.002% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `physics driver:convection_writeback`<br>物理对流后结果回写 | 0.0620 | 0.001% | 0.001% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `physics driver:radiation_lw_driver:RRTMG_lw:array_t`<br>array_t | 0.0526 | 0.000% | 0.001% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `physics driver:convection_prep`<br>物理对流前数据准备 | 0.0035 | 0.000% | 0.000% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:wet scavenging:wet_scheme:wet_incloud`<br>云内清除分支 | 0.0034 | 0.000% | 0.000% | 审计独立循环、临时数组和主机调用；保留存在依赖的顺序 |
| `chemistry driver:dry deposition:dd_budget_update`<br>沉降收支更新 | 0.0009 | 0.000% | 0.000% | 先确认分支启用且时间足够，再迁移数组差分/累加 |

## 待迁移模块的父级总计

这些阶段包含上表子项，不能与上表相加。父项未覆盖的部分可能仍有准备和未拆分计算，不能将父max减子max推算精确独占时间。

| 阶段 | max(s) | 占总运行 | 判断 |
|---|---:|---:|---|
| `chemistry driver:chem convection` | 267.5422 | 2.532% | 待迁移计算；包含已列细项，不能重复计入 |
| `chemistry_addtend` | 143.4799 | 1.358% | 待迁移计算；包含已列细项，不能重复计入 |
| `chemistry driver:wet scavenging` | 83.6126 | 0.791% | 待迁移计算；包含已列细项，不能重复计入 |
| `chemistry driver:wet scavenging:wet_scheme` | 82.7595 | 0.783% | 待迁移计算；包含已列细项，不能重复计入 |
| `chemistry driver:optical driver:optical_prep_sectional` | 58.1430 | 0.550% | 待迁移计算；包含已列细项，不能重复计入 |
| `chemistry driver:emission driver` | 28.4805 | 0.270% | 待迁移计算；包含已列细项，不能重复计入 |

## 混合阶段，需要继续区分计算和管理

含配置查询、数组准备、原有同步或日志的阶段，不将整个计时都当成可迁移数值计算。

| 阶段 | max(s) | 占总运行 | 判断 |
|---|---:|---:|---|
| `microphysics` | 224.5350 | 2.125% | 需继续拆分确认；包含已列细项，不能重复计入 |
| `atm_advance_mono_for_chem:mono_work_prep` | 156.4489 | 1.481% | 需继续拆分确认；包含已列细项，不能重复计入 |
| `physics driver:convection_driver` | 143.4797 | 1.358% | 需继续拆分确认；包含已列细项，不能重复计入 |
| `physics driver:pbl_driver` | 88.3572 | 0.836% | 需继续拆分确认；包含已列细项，不能重复计入 |
| `physics driver:cloudiness_driver` | 84.1148 | 0.796% | 需继续拆分确认；包含已列细项，不能重复计入 |
| `atm_advance_mono_for_scalars:mono_work_prep` | 62.4081 | 0.591% | 需继续拆分确认；包含已列细项，不能重复计入 |
| `physics driver:surface_layer_driver` | 7.3764 | 0.070% | 需继续拆分确认；包含已列细项，不能重复计入 |
| `atm_advance_mono_for_chem:mono_work_prep:prep_density_update` | 2.0367 | 0.019% | 需继续拆分确认；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_work_prep:prep_density_update` | 2.0351 | 0.019% | 需继续拆分确认；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:wet scavenging:wet_state_prepare` | 0.4276 | 0.004% | 需继续拆分确认；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:dry deposition:dd_aer_resistance` | 0.0296 | 0.000% | 需继续拆分确认；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:wet scavenging:wet_budget_update` | 0.0007 | 0.000% | 需继续拆分确认；该计时阶段（仍可能含未拆出的准备） |

## 已有Kokkos的计算和接口

这些部分可以继续优化，但不属于本次“从未使用Kokkos改成Kokkos”的新增迁移候选。接口中可能仍有Fortran准备。

| 阶段 | max(s) | 占总运行 | 判断 |
|---|---:|---:|---|
| `atm_advance_for_chem:advance_scalars_yakl` | 1627.9022 | 15.408% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:optical driver:mieaer` | 716.6056 | 6.783% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_compute_dyn_tend` | 515.5767 | 4.880% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:dry deposition:dd_mixactivate` | 513.1550 | 4.857% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_afcomm` | 478.5480 | 4.529% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_precomm:pc_div_limiter` | 396.1890 | 3.750% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_for_scalars:advance_scalars_yakl` | 353.0798 | 3.342% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_compute_solve_diagnostics` | 293.6054 | 2.779% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `physics driver:radiation_sw_driver:RRTMG_sw` | 250.1166 | 2.367% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_precomm:pc_horiz_flux` | 228.2531 | 2.160% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:dry deposition:dd_settling` | 207.7632 | 1.966% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_recover_large_step_variables` | 190.9633 | 1.807% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_afcomm` | 185.5178 | 1.756% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_acoustic_step` | 160.4536 | 1.519% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_precomm:pc_vert_upwind` | 126.8981 | 1.201% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_precomm:pc_div_limiter` | 125.6270 | 1.189% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `physics_get_tend` | 107.4112 | 1.017% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_precomm:pc_init` | 102.8972 | 0.974% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_precomm:pc_horiz_flux` | 76.7722 | 0.727% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_precomm:pc_vert_flux` | 38.1988 | 0.362% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_precomm:pc_vert_upwind` | 37.7083 | 0.357% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:dry deposition:vertmx:vertmx_solve_writeback` | 20.3522 | 0.193% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_precomm:pc_init` | 20.2530 | 0.192% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_precomm:pc_vert_flux` | 17.5295 | 0.166% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:dry deposition:vertmx:vertmx_geometry` | 1.7952 | 0.017% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_precomm:pc_horiz_lbc` | 0.3414 | 0.003% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:dry deposition:vertmx:vertmx_workspace` | 0.1847 | 0.002% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_precomm:pc_horiz_lbc` | 0.0780 | 0.001% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_precomm:pc_ghost` | 0.0528 | 0.000% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_precomm:pc_ghost` | 0.0142 | 0.000% | 已有Kokkos及接口；该计时阶段（仍可能含未拆出的准备） |

## 通信、文件I/O和主机管理

MPI交换、PIO写入、字段查询和内存管理不直接当成数值kernel迁移。该表也包含父子计时，不能简单求和。

| 阶段 | max(s) | 占总运行 | 判断 |
|---|---:|---:|---|
| `total time:write output files` | 1890.1668 | 17.890% | 文件I/O；包含已列细项，不能重复计入 |
| `total time:write output files:iAMAS_writeStream` | 1860.5984 | 17.610% | 文件I/O；该计时阶段（仍可能含未拆出的准备） |
| `total time:iAMAS_writeStream` | 146.7723 | 1.389% | 文件I/O；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_halo` | 136.2345 | 1.289% | 通信；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_work_prep:prep_source_halo` | 132.0835 | 1.250% | 通信；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_halo` | 97.3863 | 0.922% | 通信；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_work_prep:prep_source_halo` | 55.9487 | 0.530% | 通信；该计时阶段（仍可能含未拆出的准备） |
| `microphysics:mp_allocate` | 2.5823 | 0.024% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `total time:iAMAS_readStream` | 0.1964 | 0.002% | 文件I/O；该计时阶段（仍可能含未拆出的准备） |
| `microphysics:mp_release` | 0.0640 | 0.001% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `physics driver:phys_allocate` | 0.0440 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `chemistry driver:dry deposition:dd_get_fields` | 0.0390 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `physics driver:phys_release` | 0.0347 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `chemistry_addtend:chem_src_get_fields` | 0.0067 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `chemistry_addtend:chem_src_get_emission` | 0.0038 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_scalars:mono_work_prep:prep_get_field` | 0.0023 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `atm_advance_mono_for_chem:mono_work_prep:prep_get_field` | 0.0023 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |
| `total time:Exch all halos in read stream` | 0.0015 | 0.000% | 通信；该计时阶段（仍可能含未拆出的准备） |
| `total time:Postread in read stream` | 0.0001 | 0.000% | 主机管理；该计时阶段（仍可能含未拆出的准备） |

## 解释范围

- max是各rank累计耗时的最大值，不同阶段最大值可能来自不同rank；占总运行=max/根计时max。
- 初始化、周期性辐射、正式输出均按12小时配置执行；条件关闭的分支没有可用于评估其激活后成本的时间。
- 父子计时包含关系不可相加；表中的叶阶段只表示没有更细已记录子项，不保证内部没有接口准备。
- 垂直递推保持层间顺序，优先并行不同网格柱或独立物种；不能把所有垂直层直接同时求解。
- 迁移候选依据前次源码审计和本次新增计时标签；计时实验本身未修改数值算法。

未执行标签：emis_anthropogenic, emis_biomass。

[完整层级计时](hotspot-table.md) · [新增阶段耗时](新增阶段耗时.md) · [候选分类数据](kokkos-migration-summary.json) · [原始日志](log.atmosphere.0000.out)
