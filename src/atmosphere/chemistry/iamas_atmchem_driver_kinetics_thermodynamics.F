!************************************************************************
! This computer software was prepared by University of Science and Technology of China (USTC).
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute gases and aerosols concentrations for chem(:,:,:,:).
! * Author: Zihan Xia from AEMOL USTC, xiazh@mail.ustc.edu.cn
! Last update: zhxia 09/2025
!
! Contact:
! Chun Zhao, PhD
! Professor
! University of Science and Technology of China
! Hefei, Anhui, 230026
! Email: ChunZhao@ustc.edu.cn
!
! References:
! [1] Xia Z, Zhao C, Du Q, et al. Advancing sophisticated photochemistry simulation in atmospheric numerical models with Artificial Intelligence PhotoChemistry (AIPC) scheme using the feature‐mapping subspace self‐attention algorithm[J]. Journal of Advances in Modeling Earth Systems, 2025, 17(8): e2024MS004646.
! [2] Xia Z, Zhao C, Yang Z, et al. Toward a learnable Artificial Intelligence Model for Aerosol Chemistry and Interactions (AIMACI) based on the Multi-Head Self-Attention algorithm[J]. Atmospheric Chemistry and Physics, 2025, 25(12): 6197-6218.
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_driver_kinetics_thermodynamics
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_atmchem_vars
 use iamas_atmchem_configs
 use iamas_atmchem_utilities
 use iamas_atmphys_vars
 use module_linoz
 use module_sia_bins_simple
 use module_soa_bins_simple 
 use module_mosaic_driver
 use module_cmu_cloudchem
#ifdef USE_KPP
 use module_KPP_saprc99
#endif

#ifdef USE_DNN
 use module_AIPC_saprc99
 use module_AIMACI_mosaic
#endif

 implicit none
 private
 public:: driver_kinetics_thermodynamics

CONTAINS

!=================================================================================================================
 subroutine driver_kinetics_thermodynamics(curr_secs, itimestep, itimestepc, dtstepc, state, time_lev, diag, diag_physics, tend_chemistry, mesh, configs, its, ite)
!=================================================================================================================

 IMPLICIT NONE
 real(kind=RKIND),intent(in):: curr_secs
 integer, intent(in) :: itimestep, itimestepc
 real(kind=RKIND), intent(in) :: dtstepc
 type(iamas_pool_type),intent(inout):: state
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in   ):: configs
 type(iamas_pool_type),intent(in   ):: diag, diag_physics
 type(iamas_pool_type),intent(inout):: tend_chemistry
 integer, intent(in) :: its, ite, time_lev
!
! local variables
!
 real,pointer :: config_chemdt
 character(len=StrKIND),pointer :: chem_opt, phot_opt, gaschem_opt, aerchem_sia_opt, &
                                   aerchem_soa_opt, stra_gaschem_opt, cloudchem_opt
 real(kind=RKIND),dimension(:,:,:),pointer :: ph
 real(kind=RKIND),dimension(:,:,:),pointer :: tend_totchem
 real(kind=RKIND),dimension(:),pointer :: coszr
 integer, dimension(:), pointer :: iLev_DT
 integer, dimension(:), pointer :: is_trop_layer ! save the tropopause layer index. @zhxia
 real(kind=RKIND), dimension(:,:), pointer :: th_phy
 logical, dimension(:,:),allocatable :: is_trop ! bool array for mask. @zhxia 
 integer, pointer :: nCells, nVertLevels
 integer :: iCell, k, k_trop
 real(kind=RKIND),dimension(:,:)  ,pointer ::  relhum
 logical :: do_volatilize ! for simple SIA scheme
 logical,pointer :: config_do_chem_budget
 integer :: vbs_nbin(1) 

 ! linoz scheme vars
 character(len=StrKIND), pointer :: config_linoz_data_file

 call iamas_pool_get_config(configs,'config_chemdt',config_chemdt)
 call iamas_pool_get_config(configs,'config_chem_opt',chem_opt)
 call iamas_pool_get_config(configs,'config_phot_opt',phot_opt)
 call iamas_pool_get_config(configs,'config_gaschem_opt',gaschem_opt)
 call iamas_pool_get_config(configs,'config_aerchem_sia_opt',aerchem_sia_opt)
 call iamas_pool_get_config(configs,'config_aerchem_soa_opt',aerchem_soa_opt)
 call iamas_pool_get_config(configs,'config_stra_gaschem_opt',stra_gaschem_opt) !stratosphere
 call iamas_pool_get_config(configs,'config_cloudchem_opt',cloudchem_opt)
 call iamas_pool_get_config(configs,'config_do_chem_budget',config_do_chem_budget)

 if (config_do_chem_budget) then
     chem_old(1:num_chem,kts:kte,its:ite) = chem(1:num_chem,kts:kte,its:ite)
 end if

!--------------------------------------------
!! Photochemistry
!--------------------------------------------
 ! troposphere
 if (gaschem_opt .ne. 'no_gaschem') then
     !use diag var "iLev_DT" to set the level of tropopause. iLev_DT + 1 = bottom of stratosphere. zhxia@202510
     call iamas_pool_get_array(diag, 'iLev_DT', iLev_DT)
     call iamas_pool_get_dimension(mesh, 'nCells', nCells)
     call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
     call iamas_pool_get_array(diag, 'theta', th_phy)
     call iamas_pool_get_array(diag, 'is_trop_layer', is_trop_layer)

     allocate(is_trop(nVertLevels, nCells))
     ! Use the Hybrid PV-Theta Method to determine tropopause. zhxia@202602
     do iCell = 1, nCells
        
        k_trop = nVertLevels + 1
        do k = 1, nVertLevels
           if (th_phy(k, iCell) >= 380.0) then
              k_trop = min(k,iLev_DT(iCell))
              exit
           endif
        enddo

        is_trop_layer(iCell) = k_trop       
 
        do k = 1, nVertLevels
           is_trop(k, iCell) = (k <= is_trop_layer(iCell))
        end do
     enddo

     call iamas_pool_get_array(tend_chemistry,'ph'  , ph )
     call iamas_pool_get_array(diag_physics,'coszr'     ,coszr     )

     call chemistry_debug_message('photochemistry calculation')

     call chemistry_debug_message("config_gaschem_opt = "//trim(gaschem_opt))

     gaschem_opt_select: SELECT CASE(trim(gaschem_opt))
     CASE ("AIPC_saprc99") !zhxia@202507
#ifdef USE_DNN
         if  ((abs(config_chemdt - 3600.0) < 1e-6) .AND. (trim(phot_opt) /= 'no_phot'))  then
             if (trim(chem_opt) == 'bin04aer_gas_aq' .OR. &
                 trim(chem_opt) == 'bin04aer_gas') then 

                 if (itimestepc>=2) then !zhxia: itimestepc=1---->curr_secs=dt!=chemdt. But AI model was trained on chemdt=1hour.
                    call chemistry_debug_message('photochemistry driver: call AIPC_saprc99_driver')
       
                    call AIPC_saprc99_driver(                                                 &
                            chem = chem(1:num_chem,kts:kte,its:ite),                             &
                            rho_phy = rho_phy(kts:kte,its:ite), t_phy = t_phy(kts:kte,its:ite),  &
                            p_phy = pres_hyd(kts:kte,its:ite), qv = qv(kts:kte,its:ite),         &
                            coszen = coszr(its:ite), ph = ph(1:num_phr,kts:kte,its:ite),         &
                            is_trop = is_trop(kts:kte,its:ite),                                  &
                            its=its, ite=ite, kts=kts, kte=kte                                   &
                            )
        
                    call chemistry_debug_message('AIPC_saprc99_driver success')
                 else
                    call chemistry_debug_message('itimestepc=1---->curr_secs=dt!=chemdt. But AIPC model was trained on chemdt=1hour.')
                 endif
             else 
                 call chemistry_error_fatal('Error: AIPC_saprc99_driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use photochemistry scheme (AIPC). The rest of chem_opt needs to set config_gaschem_opt = no_gaschem.')
             end if
         else 
             call chemistry_error_fatal('Error: AIPC_saprc99_driver fail. Only config_chemdt = 3600 (1 hour) can use photochemistry scheme (AIPC_saprc99).')
         end if
#else
         call chemistry_error_fatal('--- Error: config_gaschem_opt=AIPC_saprc99 but atmosphere_model was not compiled with DNN=true')
#endif

     CASE ("KPP_saprc99") !zhxia@202512
#ifdef USE_KPP
         if (trim(chem_opt) == 'bin04aer_gas_aq' .OR. &
             trim(chem_opt) == 'bin04aer_gas') then 

             call chemistry_debug_message('photochemistry driver: call KPP_saprc99_driver')
       
             call KPP_saprc99_driver( itimestepc, mesh, state,                                         &
                            chem = chem(1:num_chem,kts:kte,its:ite), dtstepc = dtstepc,          &
                            rho_phy = rho_phy(kts:kte,its:ite), t_phy = t_phy(kts:kte,its:ite),  &
                            p_phy = pres_hyd(kts:kte,its:ite), qv = qv(kts:kte,its:ite),         &
                            ph = ph(1:num_phr,kts:kte,its:ite),                                  &
                            is_trop = is_trop(kts:kte,its:ite),                                  &
                            its=its, ite=ite, kts=kts, kte=kte                                   &
                            )
        
             call chemistry_debug_message('KPP_saprc99_driver success')
         else
             call chemistry_error_fatal('Error: KPP_saprc99_driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use photochemistry scheme (KPP_saprc99). The rest of chem_opt needs to set config_gaschem_opt = no_gaschem.')
         endif 
#else
         call chemistry_error_fatal('--- Error: config_gaschem_opt=KPP_saprc99 but atmosphere_model was not compiled with KPP=true')
#endif

     END SELECT gaschem_opt_select

     deallocate(is_trop)
 
 end if
 ! stratosphere
 if (stra_gaschem_opt .ne. 'no_stra_gaschem') then
     call iamas_pool_get_config(configs , 'config_linoz_data_file', config_linoz_data_file)

     call chemistry_debug_message('stratospheric photochemistry calculation')

     call chemistry_debug_message("config_stra_gaschem_opt = "//trim(stra_gaschem_opt))

     stra_gaschem_opt_select: SELECT CASE(trim(stra_gaschem_opt))
     CASE ("linoz") ! only for ozone. zhxia@202512
         if (trim(chem_opt) == 'bin04aer_gas_aq' .OR. &
             trim(chem_opt) == 'bin04aer_gas') then 
             call chemistry_debug_message('stratospheric photochemistry driver: call linoz_driver')

             ! Read the input data file for the Linoz stratospheric ozone chemistry at first time (only perform once).
             if (itimestepc==1) then

                 call linoz_read( state = state, linoz_data_file = config_linoz_data_file )
             end if

             ! Compute linoz 
             call linoz_driver(dtstepc = dtstepc, diag =diag, state = state, mesh = mesh, configs = configs        , &
                               its=its, ite=ite, kts=kts, kte=kte                                                    &
                               )
        
             call chemistry_debug_message('linoz_driver success')
         else 
             call chemistry_error_fatal('Error: linoz_driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use stratospheric photochemistry scheme (linoz). The rest of chem_opt needs to set config_stra_gaschem_opt = no_stra_gaschem.')
         end if
          
     END SELECT stra_gaschem_opt_select
 end if

!------------------------------------------
! Cloud aqueous chemistry
!------------------------------------------
 gas_aqfrac(1:num_gas,kts:kte,its:ite) = 0.0_RKIND

 if (cloudchem_opt .ne. 'no_cloudchem') then

     call chemistry_debug_message('cloud aqueous chemistry calculation')

     call chemistry_debug_message("config_cloudchem_opt = "//trim(cloudchem_opt))

     cloudchem_opt_select: SELECT CASE(trim(cloudchem_opt))
     CASE ("CMU_cloudchem")
         if (trim(chem_opt) == 'bin04aer_gas_aq') then
             call iamas_pool_get_array(diag_physics, 'cldfrac', cldfrac)
             call iamas_pool_get_array(tend_chemistry, 'ph', ph)

             if ((p_ph_no2 < param_first_scalar) .or. (p_ph_no2 > num_phr)) then
                 call chemistry_error_fatal('Error: CMU_cloudchem cannot find valid ph_no2 photolysis index.')
             endif

             call chemistry_debug_message('cloudchem driver: call cmu_cloudchem_driver')

             call cmu_cloudchem_driver(                                               &
                  id = 1, ktau = itimestep, ktauc = itimestepc, dtstepc = dtstepc,       &
                  num_moist = num_moist,                                                 &
                  t_phy = t_phy(kts:kte,its:ite), rho_phy = rho_phy(kts:kte,its:ite),    &
                  p_phy = pres_hyd(kts:kte,its:ite),                                     &
                  moist = moist(1:num_moist,kts:kte,its:ite),                            &
                  chem = chem(1:num_chem,kts:kte,its:ite),                               &
                  gas_aqfrac = gas_aqfrac(1:num_gas,kts:kte,its:ite),                    &
                  numgas_aqfrac = num_gas,                                               &
                  cldfra = cldfrac(kts:kte,its:ite),                                     &
                  ph_no2 = ph(p_ph_no2,kts:kte,its:ite),                                 &
                  its = its, ite = ite, kts = kts, kte = kte                             &
                  )

             call chemistry_debug_message('cmu_cloudchem_driver success')
         else
             call chemistry_error_fatal('Error: CMU_cloudchem currently requires chem_opt = bin04aer_gas_aq.')
         endif
     END SELECT cloudchem_opt_select

 end if
!------------------------------------------
! Aerosol chemistry and interactions
!------------------------------------------
! *** SIA ***
 if (aerchem_sia_opt .ne. 'no_sia') then

     call chemistry_debug_message('secondary inorganic aerosol (SIA) chemistry and interactions calculation')

     call chemistry_debug_message("config_aerchem_sia_opt = "//trim(aerchem_sia_opt))

     aerchem_sia_opt_select: SELECT CASE(trim(aerchem_sia_opt))
     CASE ("AIMACI_mosaic") !zhxia@202508
#ifdef USE_DNN
         if  (abs(config_chemdt - 3600.0) < 1e-6) then
             if (trim(chem_opt) == 'bin04aer_gas_aq' .or. trim(chem_opt) == 'bin04aer_gas') then
                 if (itimestepc>=2) then 
                    call chemistry_debug_message('aerosol driver: call AIMACI_mosaic_driver')
        
                    call AIMACI_mosaic_driver(                                                 &
                         chem = chem(1:num_chem,kts:kte,its:ite),                              &
                         rho_phy = rho_phy(kts:kte,its:ite), t_phy = t_phy(kts:kte,its:ite), &
                         p_phy = pres_hyd(kts:kte,its:ite), qv = qv(kts:kte,its:ite),          &
                         its=its, ite=ite, kts=kts, kte=kte                                    &
                         )
        
                    call chemistry_debug_message('AIMACI_mosaic_driver success')
                 else
                    call chemistry_debug_message('itimestepc=1---->curr_secs=dt!=chemdt. But AIMAIC model was trained on chemdt=1hour.')
                 endif
             else 
                 call chemistry_error_fatal('Error: AIMACI_mosaic_driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use aerosol chemistry scheme (AIMACI). The rest of chem_opt needs to set config_aerchem_sia_opt = no_sia.')
             endif
         else 
             call chemistry_error_fatal('Error: AIMACI_mosaic_driver fail. Only config_chemdt = 3600 (1 hour) can use aerosol chemistry scheme (AIMACI).')
         end if
#else
         call chemistry_error_fatal('--- Error: config_aerchem_sia_opt=AIMACI_mosaic but atmosphere_model was not compiled with DNN=true')
#endif
    ! xyruan add simple SIA scheme
      CASE ("simple_SIA")
         if (trim(chem_opt) == 'bin04aer_gas_aq' .or. trim(chem_opt) == 'bin04aer_gas') then
             call iamas_pool_get_array(diag ,'relhum' , relhum)

             call chemistry_debug_message('aerosol driver: call simple_SIA_driver')

             call simple_SIA_driver( dtstepc = dtstepc, t_phy = t_phy(kts:kte,its:ite), &
                                     relhum = relhum(kts:kte,its:ite), &
                                     chem = chem(1:num_chem,kts:kte,its:ite), &
                                     do_volatilize = .true.,                  &
                                     its=its, ite=ite, kts=kts, kte=kte )

             call chemistry_debug_message('simple_SIA_driver success')
         else
             call chemistry_error_fatal('Error: simple_SIA_driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use aerosol chemistry scheme (simple_SIA). The rest of chem_opt needs to set config_aerchem_sia_opt = no_sia.')
         endif
    ! xyruan add mechanism-based MOSAIC scheme
      CASE ("MOSAIC_mech")
         if (trim(chem_opt) == 'bin04aer_gas_aq' .or. trim(chem_opt) == 'bin04aer_gas') then 
             vbs_nbin = -99  !initialize vbs_nbin for non-vbs cases
             call chemistry_debug_message('aerosol driver: call MOSAIC_mech_driver')

             call MOSAIC_mech_driver(                                                        &
                  id = 1, curr_secs = curr_secs, ktau = itimestep, dtstep = dt_dyn,             &
                  ktauc = itimestepc, dtstepc = dtstepc, configs = configs,                     &
                  t_phy = t_phy(kts:kte,its:ite),  rho_phy = rho_phy(kts:kte,its:ite),          &
                  p_phy = pres_hyd(kts:kte,its:ite),                                            &
                  moist = moist(1:num_moist,kts:kte,its:ite),                                   &
                  chem = chem(1:num_chem,kts:kte,its:ite),                                      &
                  vbs_nbin = vbs_nbin,                                                          &
                  its=its, ite=ite, kts=kts, kte=kte                                            &
                  )
             
             call chemistry_debug_message('MOSAIC_mech_driver success')
         else
             call chemistry_error_fatal('Error: MOSAIC_mech_driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use aerosol chemistry scheme (MOSAIC_mech). The rest of chem_opt needs to set config_aerchem__sia_opt = no_sia.')
         endif
     END SELECT aerchem_sia_opt_select
 
 end if

! *** SOA ***
 if (aerchem_soa_opt .ne. 'no_soa') then

     call chemistry_debug_message('secondary organic aerosol (SOA) chemistry and interactions calculation')

     call chemistry_debug_message("config_aerchem_soa_opt = "//trim(aerchem_soa_opt))

     aerchem_soa_opt_select: SELECT CASE(trim(aerchem_soa_opt))
    ! xyruan add simple SOA scheme
     CASE ("simple_SOA")
         if (trim(chem_opt) == 'bin04aer_gas_aq' .or. trim(chem_opt) == 'bin04aer_gas') then

             call chemistry_debug_message('aerosol driver: call simple_SOA_driver')

             call simple_SOA_driver( dtstepc = dtstepc, t_phy = t_phy(kts:kte,its:ite), &
                                     rho_phy = rho_phy(kts:kte,its:ite), &
                                     chem = chem(1:num_chem,kts:kte,its:ite), &
                                     its=its, ite=ite, kts=kts, kte=kte )

             call chemistry_debug_message('simple_SOA_driver success')
         else
             call chemistry_error_fatal('Error: simple_SOA_driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use aerosol chemistry scheme (simple_SOA). The rest of chem_opt needs to set config_aerchem_soa_opt = no_soa.')
         endif
     END SELECT aerchem_soa_opt_select

 end if

 if (config_do_chem_budget) then
     call iamas_pool_get_array(tend_chemistry,'tend_totchem',tend_totchem)
     tend_totchem(1:num_chem,kts:kte,its:ite) = &
          (chem(1:num_chem,kts:kte,its:ite) - chem_old(1:num_chem,kts:kte,its:ite)) / dtstepc
 end if

end subroutine driver_kinetics_thermodynamics

END MODULE iamas_atmchem_driver_kinetics_thermodynamics
