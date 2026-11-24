! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_driver_microphysics
 use iamas_kind_types
 use iamas_pool_routines
 use iamas_timer, only : iamas_timer_start, iamas_timer_stop
 
 use iamas_atmphys_constants
 use iamas_atmphys_init_microphysics
 use iamas_atmphys_interface
 use iamas_atmphys_vars

!wrf physics:
 use module_mp_kessler
 use module_mp_thompson
 use module_mp_wsm6
 use module_mp_morr_two_moment
 use module_mp_p3

 implicit none
 private
 public:: allocate_microphysics,   &
          deallocate_microphysics, &
          driver_microphysics,     &
          microphysics_init


!iAMAS driver for parameterization of cloud microphysics processes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in iamas_atmphys_driver_microphysics:
! ------------------------------------------------
! allocate_microphysics     : allocate local arrays for parameterization of cloud microphysics.
! deallocate_microphysics   : deallocate local arrays for parameterization of cloud microphysics.
! microphysics_init         : initialization of individual cloud microphysics schemes.
! driver_microphysics       : main driver (called from iamas_atm_time_integration).
! precip_from_iAMAS          : initialize timestep local arrays for precipitation.
! precip_to_iAMAS            : copy local arrays to iAMAS arrays.
! compute_radar_reflectivity: compute radar reflectivities.
! compute_relhum            : compute relative humidity.
!
! WRF physics called from microphysics_driver:
! --------------------------------------------
! * module_mp_kessler : Kessler cloud microphysics.
! * module_mp_thompson: Thompson cloud microphysics.   
! * module_mp_wsm6    : WSM6 cloud microphysics.
! * module_mp_morrison: Morrison cloud microphysics.
! * module_mp_p3      : P3 cloud microphysics.(Predicted Particle Properties)
!
! comments:
! unlike all the other physics parameterizations, parameterizations of cloud microphysics schemes
! are called at the bottom of subroutine atm_srk3 in module atm_time_integration.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed call to the Thompson cloud microphysics scheme until scheme is updated to that in WRF revision 3.5.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added subroutine compute_relhum to calculate the relative humidity using the functions rslf and rsif from 
!   the Thompson cloud microphysics scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-12. 
! * removed the argument tend from the call to microphysics_from_iAMAS (not needed).
!   Laura D. Fowler (laura@ucar.edu) / 2013-11-07.
! * in call to subroutine wsm6, replaced the variable g (that originally pointed to gravity) with gravity, 
!   for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * moved the variable relhum from the diag_physics to the diag pool. Changed the argument
!   list for the subroutine compute_relhum accordingly.
!   Laura D. Fowler (laura@ucar.edu) / 2015-04-22.
! * added parameterization of the Thompson cloud microphysics from WRF version 3.8.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-28.
! * in subroutine compute_relhum, multiply relhum by 100. so that it has the same unit as in the initial
!   conditions.
!   Laura D. Fowler (laura@ucar.edu) / 2016-06-20.
! * added parameterization of the WSM6 cloud microphysics from WRF version 3.8.1. To initialize WSM6 as in its
!   original version, set the hail_option to 0.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-19.
! * since we removed the local variable microp_scheme from iamas_atmphys_vars.F, now defines microp_scheme as a
!   pointer to config_microp_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.

!--- initialization option for WSM6 from WRF version 3.8.1. this option could also be set as a namelist parameter.
 integer,parameter:: hail_opt = 0


 contains


!=================================================================================================================
 subroutine allocate_microphysics(configs)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 character(len=StrKIND),pointer:: p3_scheme

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call iamas_pool_get_config(configs,'config_mp_p3_scheme', p3_scheme)

!sounding variables:
 if(.not.allocated(rho_p)    ) allocate(rho_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(th_p)     ) allocate(th_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pi_p)     ) allocate(pi_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres_p)   ) allocate(pres_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(z_p)      ) allocate(z_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(dz_p)     ) allocate(dz_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(w_p)      ) allocate(w_p(ims:ime,kms:kme,jms:jme)   )

!mass mixing ratios:
 if(.not.allocated(qv_p)     ) allocate(qv_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qc_p)     ) allocate(qc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qr_p)     ) allocate(qr_p(ims:ime,kms:kme,jms:jme))

 !surface precipitation:
 if(.not.allocated(rainnc_p) ) allocate(rainnc_p(ims:ime,jms:jme) )
 if(.not.allocated(rainncv_p)) allocate(rainncv_p(ims:ime,jms:jme))

!jwfeng: 03/22/2020 add for wrfchem
#ifdef DO_CHEMISTRY
!jwfeng: 03/21/2021 add for Ice_nucleation feedback
 if(.not.allocated(dust_num_p) ) allocate(dust_num_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(sol_num_p)  ) allocate(sol_num_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(qndrop_p) ) allocate(qndrop_p(ims:ime,kms:kme,jms:jme) )
#endif

 microp_select: select case(microp_scheme)

    case ("mp_thompson","mp_wsm6","mp_morrison")
       !mass mixing ratios:
       if(.not.allocated(qi_p)        ) allocate(qi_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qs_p)        ) allocate(qs_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qg_p)        ) allocate(qg_p(ims:ime,kms:kme,jms:jme))

       !surface precipitation:
       if(.not.allocated(sr_p)        ) allocate(sr_p(ims:ime,jms:jme)        )
       if(.not.allocated(snownc_p)    ) allocate(snownc_p(ims:ime,jms:jme)    )
       if(.not.allocated(snowncv_p)   ) allocate(snowncv_p(ims:ime,jms:jme)   )
       if(.not.allocated(graupelnc_p) ) allocate(graupelnc_p(ims:ime,jms:jme) )
       if(.not.allocated(graupelncv_p)) allocate(graupelncv_p(ims:ime,jms:jme))

       !cloud water,cloud ice,and snow effective radii:
       if(.not.allocated(recloud_p) ) allocate(recloud_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(reice_p)   ) allocate(reice_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(resnow_p)  ) allocate(resnow_p(ims:ime,kms:kme,jms:jme) )

    microp2_select: select case(microp_scheme)

       case("mp_thompson")
          !number concentrations:
          if(.not.allocated(ntc_p)) allocate(ntc_p(ims:ime,jms:jme))
          if(.not.allocated(muc_p)) allocate(muc_p(ims:ime,jms:jme))
          if(.not.allocated(ni_p) ) allocate(ni_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(nr_p) ) allocate(nr_p(ims:ime,kms:kme,jms:jme))

          if(.not.allocated(rainprod_p)) allocate(rainprod_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(evapprod_p)) allocate(evapprod_p(ims:ime,kms:kme,jms:jme))

          !cloud water/ice sink:
          if(.not.allocated(qlsink_p) ) allocate(qlsink_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(ilsink_p) ) allocate(ilsink_p(ims:ime,kms:kme,jms:jme) )
          !preciptation rate of rain/ice/snow/graupel [kg/m2/s]
          if(.not.allocated( precr_p) ) allocate( precr_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( preci_p) ) allocate( preci_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( precs_p) ) allocate( precs_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( precg_p) ) allocate( precg_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( pnc_wcd_3d_p) ) allocate( pnc_wcd_3d_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( pnc_wau_3d_p) ) allocate( pnc_wau_3d_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( pnc_rcw_3d_p) ) allocate( pnc_rcw_3d_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( pni_wfz_3d_p) ) allocate( pni_wfz_3d_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( pnc_scw_3d_p) ) allocate( pnc_scw_3d_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( pnc_gcw_3d_p) ) allocate( pnc_gcw_3d_p(ims:ime,kms:kme,jms:jme) )

          qlsink_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          ilsink_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           precr_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           preci_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           precs_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           precg_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          pnc_wcd_3d_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          pnc_wau_3d_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          pnc_rcw_3d_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          pni_wfz_3d_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          pnc_scw_3d_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          pnc_gcw_3d_p(ims:ime,kms:kme,jms:jme) = 0._RKIND

       case("mp_morrison")
          !number concentrations:
          if(.not.allocated(ni_p) ) allocate(ni_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(nr_p) ) allocate(nr_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(ng_p) ) allocate(ng_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(ns_p) ) allocate(ns_p(ims:ime,kms:kme,jms:jme))

          if(.not.allocated(rqrcuten_p) ) allocate(rqrcuten_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(rqscuten_p) ) allocate(rqscuten_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(rqicuten_p) ) allocate(rqicuten_p(ims:ime,kms:kme,jms:jme))
          !cloud water/ice sink:
          if(.not.allocated(qlsink_p) ) allocate(qlsink_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(ilsink_p) ) allocate(ilsink_p(ims:ime,kms:kme,jms:jme) )
          !preciptation rate of rain/ice/snow/graupel [kg/m2/s]
          if(.not.allocated( precr_p) ) allocate( precr_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( preci_p) ) allocate( preci_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( precs_p) ) allocate( precs_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated( precg_p) ) allocate( precg_p(ims:ime,kms:kme,jms:jme) )

          qlsink_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
          ilsink_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           precr_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           preci_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           precs_p(ims:ime,kms:kme,jms:jme) = 0._RKIND
           precg_p(ims:ime,kms:kme,jms:jme) = 0._RKIND

          if(.not.allocated(rainprod_p)) allocate(rainprod_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(evapprod_p)) allocate(evapprod_p(ims:ime,kms:kme,jms:jme))

       case default

    end select microp2_select

    case ("mp_p3")
       !mass mixing ratios:
       if(.not.allocated(qi_p)        ) allocate(qi_p(ims:ime,kms:kme,jms:jme))

       !rime ice mass mixing ratio
       if(.not.allocated(qir_p) ) allocate(qir_p(ims:ime,kms:kme,jms:jme))
       !ice rime volume mixing ratio
       if(.not.allocated(qib_p) ) allocate(qib_p(ims:ime,kms:kme,jms:jme))

       !number concentrations:
       if(.not.allocated(ni_p) ) allocate(ni_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(nr_p) ) allocate(nr_p(ims:ime,kms:kme,jms:jme))

       !surface precipitation:
       if(.not.allocated(sr_p)        ) allocate(sr_p(ims:ime,jms:jme)        )
       if(.not.allocated(snownc_p)    ) allocate(snownc_p(ims:ime,jms:jme)    )
       if(.not.allocated(snowncv_p)   ) allocate(snowncv_p(ims:ime,jms:jme)   )

       !cloud water,cloud ice effective radii:
       if(.not.allocated(recloud_p) ) allocate(recloud_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(reice_p)   ) allocate(reice_p(ims:ime,kms:kme,jms:jme)  )

       !preivous timestep
       if(.not.allocated(qv_old_p) ) allocate(qv_old_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(th_old_p) ) allocate(th_old_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(refl10cm_p) ) allocate(refl10cm_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(vmi3d_p) ) allocate(vmi3d_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(di3d_p) ) allocate(di3d_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rhopo3d_p) ) allocate(rhopo3d_p(ims:ime,kms:kme,jms:jme))

       p3_select: select case(p3_scheme)
          case("P3_1CATEGORY_NC")
            if(.not.allocated(nc_p) ) allocate(nc_p(ims:ime,kms:kme,jms:jme))
          case("P3_2CATEGORY")
            if(.not.allocated(nc_p) ) allocate(nc_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(qi2_p) )allocate(qi2_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(ni2_p) )allocate(ni2_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(qir2_p) )allocate(qir2_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(qib2_p) )allocate(qib2_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(vmi3d_2_p) ) allocate(vmi3d_2_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(di3d_2_p) ) allocate(di3d_2_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(rhopo3d_2_p) ) allocate(rhopo3d_2_p(ims:ime,kms:kme,jms:jme))
          case("P3_1CAT_3MOM")
            if(.not.allocated(nc_p) ) allocate(nc_p(ims:ime,kms:kme,jms:jme))
            if(.not.allocated(zi_p) ) allocate(zi_p(ims:ime,kms:kme,jms:jme))
          case default
       end select p3_select

    case default

 end select microp_select

 end subroutine allocate_microphysics

!=================================================================================================================
 subroutine deallocate_microphysics(configs)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 character(len=StrKIND),pointer:: p3_scheme

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call iamas_pool_get_config(configs,'config_mp_p3_scheme', p3_scheme)

!sounding variables:
 if(allocated(rho_p)     ) deallocate(rho_p     )
 if(allocated(th_p)      ) deallocate(th_p      )
 if(allocated(pi_p)      ) deallocate(pi_p      )
 if(allocated(pres_p)    ) deallocate(pres_p    )
 if(allocated(z_p)       ) deallocate(z_p       )
 if(allocated(dz_p)      ) deallocate(dz_p      )
 if(allocated(w_p)       ) deallocate(w_p       )

!mass mixing ratios:
 if(allocated(qv_p)      ) deallocate(qv_p      )
 if(allocated(qc_p)      ) deallocate(qc_p      )
 if(allocated(qr_p)      ) deallocate(qr_p      )

 !surface precipitation:
 if(allocated(rainnc_p)  ) deallocate(rainnc_p  )
 if(allocated(rainncv_p) ) deallocate(rainncv_p )

!jwfeng: 03/22/2020 add for wrfchem
#ifdef DO_CHEMISTRY
 if(allocated(qndrop_p) ) deallocate(qndrop_p  )
 if(allocated(dust_num_p) ) deallocate(dust_num_p  )
 if(allocated(sol_num_p)  ) deallocate(sol_num_p  )
#endif

 microp_select: select case(microp_scheme)

    case ("mp_thompson","mp_wsm6","mp_morrison")
       !mass mixing ratios:
       if(allocated(qi_p)         ) deallocate(qi_p         )
       if(allocated(qs_p)         ) deallocate(qs_p         )
       if(allocated(qg_p)         ) deallocate(qg_p         )

       !surface precipitation:
       if(allocated(sr_p)         ) deallocate(sr_p         )
       if(allocated(snownc_p)     ) deallocate(snownc_p     )
       if(allocated(snowncv_p)    ) deallocate(snowncv_p    )
       if(allocated(graupelnc_p)  ) deallocate(graupelnc_p  )
       if(allocated(graupelncv_p) ) deallocate(graupelncv_p )

       !cloud water,cloud ice,and snow effective radii:
       if(allocated(recloud_p) ) deallocate(recloud_p )
       if(allocated(reice_p)   ) deallocate(reice_p   )
       if(allocated(resnow_p)  ) deallocate(resnow_p  )

    microp2_select: select case(microp_scheme)

       case("mp_thompson")
          !number concentrations:
          if(allocated(ntc_p)) deallocate(ntc_p)
          if(allocated(muc_p)) deallocate(muc_p)
          if(allocated(ni_p) ) deallocate(ni_p )
          if(allocated(nr_p) ) deallocate(nr_p )

          if(allocated(rainprod_p)) deallocate(rainprod_p)
          if(allocated(evapprod_p)) deallocate(evapprod_p)

          !cloud water/ice sink:
          if(allocated(qlsink_p) ) deallocate(qlsink_p  )
          if(allocated(ilsink_p) ) deallocate(ilsink_p  )
          !preciptation rate of rain/ice/snow/graupel [kg/m2/s]
          if(allocated( precr_p) ) deallocate( precr_p  )
          if(allocated( preci_p) ) deallocate( preci_p  )
          if(allocated( precs_p) ) deallocate( precs_p  )
          if(allocated( precg_p) ) deallocate( precg_p  )
          if(allocated(pnc_wcd_3d_p)) deallocate(pnc_wcd_3d_p)
          if(allocated(pnc_wau_3d_p)) deallocate(pnc_wau_3d_p)
          if(allocated(pnc_rcw_3d_p)) deallocate(pnc_rcw_3d_p)
          if(allocated(pni_wfz_3d_p)) deallocate(pni_wfz_3d_p)
          if(allocated(pnc_scw_3d_p)) deallocate(pnc_scw_3d_p)
          if(allocated(pnc_gcw_3d_p)) deallocate(pnc_gcw_3d_p)

       case("mp_morrison")
          !number concentrations:
          if(allocated(ni_p) ) deallocate(ni_p )
          if(allocated(nr_p) ) deallocate(nr_p )
          if(allocated(ng_p) ) deallocate(ng_p )
          if(allocated(ns_p) ) deallocate(ns_p )

          !cloud water/ice sink:
          if(allocated(qlsink_p) ) deallocate(qlsink_p  )
          if(allocated(ilsink_p) ) deallocate(ilsink_p  )
          !preciptation rate of rain/ice/snow/graupel [kg/m2/s]
          if(allocated( precr_p) ) deallocate( precr_p  )
          if(allocated( preci_p) ) deallocate( preci_p  )
          if(allocated( precs_p) ) deallocate( precs_p  )
          if(allocated( precg_p) ) deallocate( precg_p  )

          if(allocated(rqrcuten_p) ) deallocate(rqrcuten_p)
          if(allocated(rqscuten_p) ) deallocate(rqscuten_p)
          if(allocated(rqicuten_p) ) deallocate(rqicuten_p)

          if(allocated(rainprod_p)) deallocate(rainprod_p)
          if(allocated(evapprod_p)) deallocate(evapprod_p)

       case default

    end select microp2_select

    case ("mp_p3")
       !mass mixing ratios:
       if(allocated(qi_p)        ) deallocate(qi_p)

       !rime ice mass mixing ratio
       if(allocated(qir_p) ) deallocate(qir_p)
       !ice rime volume mixing ratio
       if(allocated(qib_p) ) deallocate(qib_p)

       !number concentrations:
       if(allocated(ni_p) ) deallocate(ni_p)
       if(allocated(nr_p) ) deallocate(nr_p)

       !surface precipitation:
       if(allocated(sr_p)        ) deallocate(sr_p        )
       if(allocated(snownc_p)    ) deallocate(snownc_p    )
       if(allocated(snowncv_p)   ) deallocate(snowncv_p   )

       !cloud water,cloud ice effective radii:
       if(allocated(recloud_p) ) deallocate(recloud_p)
       if(allocated(reice_p)   ) deallocate(reice_p  )

       !preivous timestep
       if(allocated(qv_old_p) ) deallocate(qv_old_p)
       if(allocated(th_old_p) ) deallocate(th_old_p)
       if(allocated(refl10cm_p) ) deallocate(refl10cm_p)
       if(allocated(vmi3d_p) ) deallocate(vmi3d_p)
       if(allocated(di3d_p) ) deallocate(di3d_p)
       if(allocated(rhopo3d_p) ) deallocate(rhopo3d_p)
       p3_select: select case(p3_scheme)
          case("P3_1CATEGORY_NC")
            !number concentrations:
            if(allocated(nc_p) ) deallocate(nc_p)
          case("P3_2CATEGORY")
            if(allocated(nc_p) ) deallocate(nc_p)
            if(allocated(qi2_p) )deallocate(qi2_p)
            if(allocated(ni2_p) )deallocate(ni2_p)
            if(allocated(qir2_p) )deallocate(qir2_p)
            if(allocated(qib2_p) )deallocate(qib2_p)
            if(allocated(vmi3d_2_p) ) deallocate(vmi3d_2_p)
            if(allocated(di3d_2_p) ) deallocate(di3d_2_p)
            if(allocated(rhopo3d_2_p) ) deallocate(rhopo3d_2_p)
          case("P3_1CAT_3MOM")
            if(allocated(nc_p) ) deallocate(nc_p)
            if(allocated(zi_p) ) deallocate(zi_p)
          case default
       end select p3_select


    case default

 end select microp_select

 end subroutine deallocate_microphysics

!=================================================================================================================
 subroutine microphysics_init(dminfo,configs,mesh,sfc_input,diag_physics)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: sfc_input

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics

!local pointer:
 character(len=StrKIND),pointer:: microp_scheme
 character(len=StrKIND),pointer:: p3_scheme
 integer :: i_err

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

  microp_select: select case(microp_scheme)

     case("mp_thompson")
        call thompson_init(l_mp_tables)
        call init_thompson_clouddroplets_foriAMAS(mesh,configs,sfc_input,diag_physics)

     case("mp_wsm6")
        call wsm6init(rho_a,rho_r,rho_s,cliq,cpv,hail_opt,.false.)

     case("mp_morrison")
        !morr_rimed_ice = 1 added temporarily by J.Gu
         call morr_two_moment_init(morr_rimed_ice = 1)

     case("mp_p3")
        call iamas_pool_get_config(configs,'config_mp_p3_scheme', p3_scheme)
        p3_select: select case(p3_scheme)
            case("P3_1CATEGORY")
                call p3_init('.',1,.false.,'WRF',i_err,.false.)
            case("P3_1CATEGORY_NC")
                call p3_init('.',1,.false.,'WRF',i_err,.false.)
            case("P3_2CATEGORY")
                call p3_init('.',2,.false.,'WRF',i_err,.false.)
            case("P3_1CAT_3MOM")
                call p3_init('.',1,.true.,'WRF',i_err,.false.)
           case default
        end select p3_select


     case default
    
  end select microp_select

 end subroutine microphysics_init

!=================================================================================================================
 subroutine driver_microphysics(clock,configs,mesh,state,time_lev,diag,diag_physics,tend,tend_physics,itimestep,its,ite)
!=================================================================================================================

!input arguments:
 type(iAMAS_Clock_type),intent(in):: clock
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh

 integer,intent(in):: time_lev
 integer,intent(in):: itimestep
 integer,intent(in):: its,ite
    
!inout arguments:
 type(iamas_pool_type),intent(inout):: state
 type(iamas_pool_type),intent(inout):: diag
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: tend
 type(iamas_pool_type),intent(inout):: tend_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 character(len=StrKIND),pointer:: p3_scheme
 real(kind=RKIND),dimension(:,:),pointer  :: delta_rc
#ifdef DO_CHEMISTRY
 real(kind=RKIND),dimension(:,:),pointer  :: tend_qndrop_mp
 integer,pointer :: index_qndrop
 character(len=StrKIND),pointer :: chem_opt
 logical :: update_qndrop
#endif

!local variables and arrays:
 logical:: log_microphysics
 integer:: i,icell,icount,istep,j,k,kk
 logical,pointer:: debug_on, config_aer_cloud_droplet_on, config_aer_cloud_ice_on, config_wetscav_on
!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write('')
!call iamas_log_write('---enter subroutine driver_microphysics:')

 call iamas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call iamas_pool_get_config(configs,'config_print_debug_info' ,debug_on )
 call iamas_pool_get_config(configs,'config_aer_cloud_droplet_on' ,config_aer_cloud_droplet_on )
 call iamas_pool_get_config(configs,'config_aer_cloud_ice_on' ,config_aer_cloud_ice_on )
 call iamas_pool_get_config(configs,'config_wetscav_on' ,config_wetscav_on )

 if (debug_on) call iamas_log_write("--- Doing allocate_microphysics")
!... allocation of microphysics arrays:
!$OMP MASTER
 call allocate_microphysics(configs)
 if (debug_on) call iamas_log_write("--- End allocate_microphysics")
!$OMP END MASTER
!$OMP BARRIER
!... initialization of precipitation related arrays:
 if (debug_on) call iamas_log_write("--- Doing precip_from_iAMAS")
 call precip_from_iAMAS(configs,diag_physics,its,ite)
 if (debug_on) call iamas_log_write("--- End precip_from_iAMAS")

!... initialization of soundings for non-hydrostatic dynamical cores.
 if (debug_on) call iamas_log_write("--- Doing microphysics_from_iAMAS")
 call microphysics_from_iAMAS(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,its,ite)
 if (debug_on) call iamas_log_write("--- End microphysics_from_iAMAS")

#ifdef DO_CHEMISTRY
 nullify(chem_opt)
 call iamas_pool_get_config(configs, 'config_chem_opt', chem_opt)
 update_qndrop = .true.
 if (associated(chem_opt)) update_qndrop = trim(chem_opt) /= 'no_chemistry'
 if (update_qndrop) then
 call iamas_pool_get_dimension(state,'index_qndrop'  ,index_qndrop)
 if (index_qndrop .gt. 0) then
    call iamas_pool_get_array(tend, 'tend_qndrop_mp', tend_qndrop_mp)
    do j = jts, jte
    do k = kts, kte
    do i = its, ite
       tend_qndrop_mp(k, i) = qndrop_p(i, k, j)
    end do
    end do
    end do
 end if
 end if
#endif
call iamas_pool_get_array(diag, 'delta_rc', delta_rc)
do j = jts, jte
do k = kts, kte
do i = its, ite
   if (trim(microp_scheme) == 'mp_kessler') then
      delta_rc(k, i) = qc_p(i, k, j)
   else
      delta_rc(k, i) = qc_p(i, k, j) + qi_p(i, k, j)
   endif
end do
end do
end do

!... call to different cloud microphysics schemes:
 microp_select: select case(microp_scheme)
    
    case ("mp_kessler")
       if (debug_on) call iamas_log_write("--- Doing kessler")
       call iamas_timer_start('Kessler')
       call kessler( &
                  t        = th_p      , qv    = qv_p  , qc     = qc_p     ,                &
                  qr       = qr_p      , rho   = rho_p , pii    = pi_p     ,                &
                  dt_in    = dt_microp , z     = z_p   , xlv    = xlv      ,                &
                  cp       = cp        , ep2   = ep_2  , svp1   = svp1     ,                &
                  svp2     = svp2      , svp3  = svp3  , svpt0  = svpt0    ,                &
                  rhowater = rho_w     , dz8w  = dz_p  , rainnc = rainnc_p ,                &
                  rainncv  = rainncv_p ,                                                    &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde   , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme   , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                   )
       call iamas_timer_stop('Kessler')
       if (debug_on) call iamas_log_write("--- End kessler")

    case ("mp_thompson")
       if (debug_on) call iamas_log_write("--- Doing mp_thompson")
       istep = 1
       call iamas_timer_start('Thompson')
       do while (istep .le. n_microp)
       if (debug_on) call iamas_log_write("------- istep = $i",intArgs=(/istep/))
       call mp_gt_driver( &
                  th        = th_p        , qv         = qv_p         , qc         = qc_p         , &
                  qr        = qr_p        , qi         = qi_p         , qs         = qs_p         , &
                  qg        = qg_p        , ni         = ni_p         , nr         = nr_p         , &
#ifdef DO_CHEMISTRY
                  nc        = qndrop_p    , nifa       = dust_num_p   , nwfa       = sol_num_p    , &
#endif
                  pii       = pi_p        , p          = pres_p       , dz         = dz_p         , & 
                  w         = w_p         , dt_in      = dt_microp    , itimestep  = itimestep    , &
                  rainnc    = rainnc_p    , rainncv    = rainncv_p    , snownc     = snownc_p     , &
                  snowncv   = snowncv_p   , graupelnc  = graupelnc_p  , graupelncv = graupelncv_p , &
                  sr        = sr_p        , rainprod   = rainprod_p   , evapprod   = evapprod_p   , &
                  re_cloud  = recloud_p   , re_ice     = reice_p      , re_snow    = resnow_p     , &
                  has_reqc  = has_reqc    , has_reqi   = has_reqi     , has_reqs   = has_reqs     , &
                  ntc       = ntc_p       , muc        = muc_p        ,                             &
!jwfeng: 03/22/2020 add for iAMAS-Chem
                  qlsink    = qlsink_p    , ilsink     = ilsink_p     ,                             &
                  precr     = precr_p     , preci      = preci_p      ,                             &
                  precs     = precs_p     , precg      = precg_p      ,                             &
!jwfeng: 06/26/2023
                  pnc_wcd_3d= pnc_wcd_3d_p, pnc_wau_3d = pnc_wau_3d_p , pnc_rcw_3d = pnc_rcw_3d_p , &
                  pni_wfz_3d= pni_wfz_3d_p, pnc_scw_3d = pnc_scw_3d_p , pnc_gcw_3d = pnc_gcw_3d_p , &
#ifdef DO_CHEMISTRY
                  use_aer_IN= config_aer_cloud_ice_on                                             , &
#endif
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde           , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme           , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                           )
       istep = istep + 1
       enddo
       call iamas_timer_stop('Thompson')
       if (debug_on) call iamas_log_write("--- End mp_thompson")

    case ("mp_wsm6")
       if (debug_on) call iamas_log_write("--- Doing wsm6")
       call iamas_timer_start('WSM6')
       call wsm6( &
                  th       = th_p        , q          = qv_p         , qc       = qc_p      , &
                  qr       = qr_p        , qi         = qi_p         , qs       = qs_p      , &
                  qg       = qg_p        , den        = rho_p        , pii      = pi_p      , &
                  p        = pres_p      , delz       = dz_p         , delt     = dt_microp , &
                  g        = gravity     , cpd        = cp           , cpv      = cpv       , &
                  rd       = R_d         , rv         = R_v          , t0c      = svpt0     , &
                  ep1      = ep_1        , ep2        = ep_2         , qmin     = epsilon   , &
                  xls      = xls         , xlv0       = xlv          , xlf0     = xlf       , &
                  den0     = rho_a       , denr       = rho_w        , cliq     = cliq      , &
                  cice     = cice        , psat       = psat         , rain     = rainnc_p  , &
                  rainncv  = rainncv_p   , snow       = snownc_p     , snowncv  = snowncv_p , &
                  graupel  = graupelnc_p , graupelncv = graupelncv_p , sr       = sr_p      , &
                  re_cloud = recloud_p   , re_ice     = reice_p      , re_snow  = resnow_p  , &
                  has_reqc = has_reqc    , has_reqi   = has_reqi     , has_reqs = has_reqs  , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde     , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme     , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte       &
               )
       call iamas_timer_stop('WSM6')
       if (debug_on) call iamas_log_write("--- End wsm6")

    case ("mp_morrison")
       if (debug_on) call iamas_log_write("--- Doing mp_morrison")
       istep = 1
       call iamas_timer_start('Morrison')
       do while(istep .le. n_microp)
          if (debug_on) call iamas_log_write("------- istep = $i",intArgs=(/istep/))
          call mp_morr_two_moment(itimestep = itimestep, th = th_p, qv = qv_p, qc = qc_p, qr = qr_p, qi = qi_p, qs = qs_p, qg = qg_p, &
                                ni = ni_p, ns = ns_p, nr = nr_p, ng = ng_p, &
                                pii = pi_p, p = pres_p, dt_in = dt_microp, dz = dz_p, w = w_p, & !rho and ht are not used in this scheme / J.Gu !Important
                                rainnc = rainnc_p, rainncv = rainncv_p, sr = sr_p, &
                                snownc = snownc_p, snowncv = snowncv_p, graupelnc = graupelnc_p, graupelncv = graupelncv_p, &
                                qrcuten = rqrcuten_p, qscuten = rqscuten_p, qicuten = rqicuten_p, &
#ifdef DO_CHEMISTRY
                                f_qndrop = config_aer_cloud_droplet_on, qndrop = qndrop_p, &
                                wetscav_on = config_wetscav_on, rainprod = rainprod_p, evapprod = evapprod_p, &
                                qlsink = qlsink_p, precr = precr_p, preci = preci_p, precs = precs_p, precg = precg_p, &
#endif
                                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde     , &
                                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme     , &
                                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte       &
                               )
          istep = istep + 1
       end do
       call iamas_timer_stop('Morrison')
       if (debug_on) call iamas_log_write("--- End mp_morrison")

    case ("mp_p3")
       if (debug_on) call iamas_log_write("--- Doing P3")
       call iamas_pool_get_config(configs,'config_mp_p3_scheme', p3_scheme)
       call iamas_timer_start('P3')
       p3_select: select case(p3_scheme)
          case("P3_1CATEGORY")
              if (debug_on) call iamas_log_write("------ Calling P3 one category")
              call mp_p3_wrapper_wrf(itimestep = itimestep, TH_3d = th_p, Qv_3d = qv_p, Qc_3d = qc_p, Qr_3d = qr_p, &
                                     Qnr_3d = nr_p, Qi1_3d = qi_p, Qir1_3d = qir_p, Qni1_3d = ni_p, Qib1_3d = qib_p, &
                                     th_old_3d = th_old_p, qv_old_3d = qv_old_p, pii = pi_p, p = pres_p, dt = dt_microp, dz = dz_p, &
                                     w = w_p, rainnc = rainnc_p, rainncv = rainncv_p, sr = sr_p, &
                                     snownc = snownc_p, snowncv = snowncv_p, n_icecat = 1, &
                                     ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde     , &
                                     ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme     , &
                                     its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     , &
                                     diag_zdbz_3d = refl10cm_p, diag_effc_3d = recloud_p, diag_effi_3d = reice_p, &
                                     diag_vmi_3d = vmi3d_p, diag_di_3d = di3d_p, diag_rhopo_3d = rhopo3d_p) ! for P3

          case("P3_1CATEGORY_NC")
              if (debug_on) call iamas_log_write("------ Calling P3 one category nc")
              call mp_p3_wrapper_wrf(itimestep = itimestep, TH_3d = th_p, Qv_3d = qv_p, Qc_3d = qc_p, Qr_3d = qr_p, &
                                     Qnr_3d = nr_p, Qi1_3d = qi_p, Qir1_3d = qir_p, Qni1_3d = ni_p, Qib1_3d = qib_p, &
                                     th_old_3d = th_old_p, qv_old_3d = qv_old_p, nc_3d = nc_p, pii = pi_p, p = pres_p, dt = dt_microp, dz = dz_p, &
                                     w = w_p, rainnc = rainnc_p, rainncv = rainncv_p, sr = sr_p, &
                                     snownc = snownc_p, snowncv = snowncv_p, n_icecat = 1, &
                                     ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde     , &
                                     ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme     , &
                                     its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     , &
                                     diag_zdbz_3d = refl10cm_p, diag_effc_3d = recloud_p, diag_effi_3d = reice_p, &
                                     diag_vmi_3d = vmi3d_p, diag_di_3d = di3d_p, diag_rhopo_3d = rhopo3d_p) ! for P3

          case("P3_2CATEGORY")
              if (debug_on) call iamas_log_write("------ Calling P3 two category")
              call mp_p3_wrapper_wrf_2cat(itimestep = itimestep, TH_3d = th_p, Qv_3d = qv_p, Qc_3d = qc_p, Qr_3d = qr_p, Qnr_3d = nr_p, &
                                     Qi1_3d = qi_p, Qir1_3d = qir_p, Qni1_3d = ni_p, Qib1_3d = qib_p, &
                                     Qi2_3d = qi2_p, Qir2_3d = qir2_p, Qni2_3d = ni2_p, Qib2_3d = qib2_p, &
                                     th_old_3d = th_old_p, qv_old_3d = qv_old_p, nc_3d = nc_p, pii = pi_p, p = pres_p, dt = dt_microp, dz = dz_p, &
                                     w = w_p, rainnc = rainnc_p, rainncv = rainncv_p, sr = sr_p, &
                                     snownc = snownc_p, snowncv = snowncv_p, n_icecat = 2, &
                                     ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde     , &
                                     ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme     , &
                                     its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     , &
                                     diag_zdbz_3d = refl10cm_p, diag_effc_3d = recloud_p, diag_effi_3d = reice_p, &
                                     diag_vmi_3d = vmi3d_p, diag_di_3d = di3d_p, diag_rhopo_3d = rhopo3d_p, &
                                     diag_vmi2_3d = vmi3d_2_p, diag_di2_3d = di3d_2_p, diag_rhopo2_3d = rhopo3d_2_p) ! for P3
          case("P3_1CAT_3MOM")
              if (debug_on) call iamas_log_write("------ Calling P3 one category 3 moment")
              call mp_p3_wrapper_wrf(itimestep = itimestep, TH_3d = th_p, Qv_3d = qv_p, Qc_3d = qc_p, Qr_3d = qr_p, &
                                     Qnr_3d = nr_p, Qi1_3d = qi_p, Qir1_3d = qir_p, Qni1_3d = ni_p, Qib1_3d = qib_p, &
                                     th_old_3d = th_old_p, qv_old_3d = qv_old_p, nc_3d = nc_p, pii = pi_p, p = pres_p, dt = dt_microp, dz = dz_p, &
                                     w = w_p, rainnc = rainnc_p, rainncv = rainncv_p, sr = sr_p, &
                                     snownc = snownc_p, snowncv = snowncv_p, n_icecat = 1, &
                                     ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde     , &
                                     ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme     , &
                                     its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     , &
                                     diag_zdbz_3d = refl10cm_p, diag_effc_3d = recloud_p, diag_effi_3d = reice_p, &
                                     diag_vmi_3d = vmi3d_p, diag_di_3d = di3d_p, diag_rhopo_3d = rhopo3d_p, &
                                     Qzi1_3d = zi_p &
                               ! not currently output
                               !     ,diag_dhmax_3d = dhmax3d
                               !     ,diag_lami_3d  = lami3d
                               !     ,diag_mui_3d = mui3d
                                     ) ! for P3
         case default
       end select p3_select

       call iamas_timer_stop('P3')
       if (debug_on) call iamas_log_write("--- End P3")


    case default

 end select microp_select
#ifdef DO_CHEMISTRY
       if (update_qndrop) then
       if (index_qndrop .gt. 0) then
          do j = jts, jte
          do k = kts, kte
          do i = its, ite
             tend_qndrop_mp(k, i) = (qndrop_p(i, k, j) - tend_qndrop_mp(k, i))/dt_dyn
          end do
          end do
          end do
       end if
       end if
#endif
 do j = jts, jte
 do k = kts, kte
 do i = its, ite
    ! Change in cloud condensate due to microphysics: after minus before.
    if (trim(microp_scheme) == 'mp_kessler') then
       delta_rc(k, i) = qc_p(i, k, j) - delta_rc(k,i)
    else
       delta_rc(k, i) = (qc_p(i, k, j) + qi_p(i, k, j)) - delta_rc(k,i)
    endif
 end do
 end do
 end do

!... calculate the 10cm radar reflectivity and relative humidity, if needed:
 if (l_diags) then
 
    !ensure that we only call compute_radar_reflectivity() if we are using an MPS that supports
    !the computation of simulated radar reflectivity:
    if(trim(microp_scheme) == "mp_wsm6"     .or. &
       trim(microp_scheme) == "mp_thompson" .or. &
       trim(microp_scheme) == "mp_morrison") then
       call compute_radar_reflectivity(configs,diag_physics,its,ite)
    else
       call iamas_log_write('*** NOTICE: NOT computing simulated radar reflectivity')
       call iamas_log_write('            since WSM6 or Thompson or Morrison microphysics scheme was not selected')
    end if

    !calculate the relative humidity over water if the temperature is strictly greater than 0.C,
    !over ice otherwise.
    call compute_relhum(diag,its,ite)

 end if

!... copy updated precipitation from the wrf-physics grid back to the geodesic-dynamics grid:
 if (debug_on) call iamas_log_write("--- Doing precip_to_iAMAS")
 call precip_to_iAMAS(clock,configs,diag_physics,its,ite)
 if (debug_on) call iamas_log_write("--- End precip_to_iAMAS")

!... copy updated cloud microphysics variables from the wrf-physics grid back to the geodesic-
!    dynamics grid:
 if (debug_on) call iamas_log_write("--- Doing microphysics_to_iAMAS")
 call microphysics_to_iAMAS(configs,mesh,state,time_lev,diag,diag_physics,tend,tend_physics,itimestep,its,ite)
 if (debug_on) call iamas_log_write("--- End microphysics_to_iAMAS")

!... deallocation of all microphysics arrays:
!$OMP BARRIER
!$OMP MASTER
 if (debug_on) call iamas_log_write("--- Doing deallocate_microphysics")
 call deallocate_microphysics(configs)
 if (debug_on) call iamas_log_write("--- End deallocate_microphysics")
!$OMP END MASTER

!call iamas_log_write('---enter subroutine driver_microphysics:')
!call iamas_log_write('')
 end subroutine driver_microphysics

!=================================================================================================================
 subroutine precip_from_iAMAS(configs,diag_physics,its,ite)
!=================================================================================================================

!input variables:
 type(iamas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!output variables:
 type(iamas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,pointer:: nCellsSolve
 real,dimension(:),pointer:: graupelncv,rainncv,snowncv,sr 

!local variables and arrays:
 integer:: i,j,k

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call iamas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call iamas_pool_get_array(diag_physics,'sr'        ,sr        )

!variables common to all cloud microphysics schemes:
 do j = jts, jte
 do i = its, ite
    rainncv_p(i,j) = 0._RKIND
    rainnc_p(i,j)  = 0._RKIND
 enddo
 enddo

 do i = its,ite
    rainncv(i) = 0._RKIND
 enddo

!variables specific to different cloud microphysics schemes:
 microp_select: select case(microp_scheme)

    case ("mp_thompson","mp_wsm6","mp_morrison")
       call iamas_pool_get_array(diag_physics,'graupelncv',graupelncv)
       call iamas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
       do j = jts, jte
       do i = its, ite
          snowncv_p(i,j)    = 0._RKIND
          graupelncv_p(i,j) = 0._RKIND
          snownc_p(i,j)     = 0._RKIND
          graupelnc_p(i,j)  = 0._RKIND
          sr_p(i,j)         = 0._RKIND
       enddo
       enddo

       do i = its,ite
          snowncv(i)    = 0._RKIND
          graupelncv(i) = 0._RKIND
          sr(i)         = 0._RKIND
       enddo 

    microp2_select: select case(microp_scheme)
       case("mp_thompson","mp_morrison")
          !jwfeng: 03/22/2020 add for wrfchem
          do j = jts, jte
          do k = kts, kte
          do i = its, ite
             qlsink_p(i,k,j) = 0._RKIND
             ilsink_p(i,k,j) = 0._RKIND
              precr_p(i,k,j) = 0._RKIND
              preci_p(i,k,j) = 0._RKIND
              precs_p(i,k,j) = 0._RKIND
              precg_p(i,k,j) = 0._RKIND
              if(microp_scheme == "mp_thompson") then
         pnc_wcd_3d_p(i,k,j) = 0._RKIND
         pnc_wau_3d_p(i,k,j) = 0._RKIND
         pnc_rcw_3d_p(i,k,j) = 0._RKIND
         pni_wfz_3d_p(i,k,j) = 0._RKIND
         pnc_scw_3d_p(i,k,j) = 0._RKIND
         pnc_gcw_3d_p(i,k,j) = 0._RKIND
              endif

          enddo
          enddo
          enddo
       case default
    end select microp2_select

    case ("mp_p3")
       call iamas_pool_get_array(diag_physics,'snowncv',snowncv)
       do j = jts, jte
       do i = its, ite
          snowncv_p(i,j)    = 0._RKIND
          snownc_p(i,j)     = 0._RKIND
          sr_p(i,j)         = 0._RKIND
       enddo
       enddo

       do i = its,ite
          snowncv(i)    = 0._RKIND
          sr(i)         = 0._RKIND
       enddo 



    case default

 end select microp_select

 end subroutine precip_from_iAMAS

!=================================================================================================================
 subroutine precip_to_iAMAS(clock,configs,diag_physics,its,ite)
!=================================================================================================================

 use iamas_timekeeping

!input arguments:
 type(iAMAS_Clock_type), intent(in):: clock
 type(iamas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,dimension(:),pointer:: i_rainnc

 real(kind=RKIND),pointer:: config_bucket_rainnc
 real(kind=RKIND),dimension(:),pointer:: precipw
 real(kind=RKIND),dimension(:),pointer:: graupelnc,rainnc,rainnc_hr,snownc
 real(kind=RKIND),dimension(:),pointer:: graupelncv,rainncv,snowncv,sr

!jwfeng: 03/22/2020 add for iAMAS-Chem
 real(kind=RKIND),dimension(:,:),pointer:: qlsink,ilsink
 real(kind=RKIND),dimension(:,:),pointer:: precr,preci,precs,precg

!jwfeng: 06/26/2023 add for iAMAS-Chem
 real(kind=RKIND),dimension(:,:),pointer:: pnc_wcd_3d, pnc_wau_3d , pnc_rcw_3d
 real(kind=RKIND),dimension(:,:),pointer:: pni_wfz_3d, pnc_scw_3d , pnc_gcw_3d

 logical,pointer:: debug_on

!local variables and arrays:
 integer:: i,j,k
 real(kind=RKIND):: rho_a

 !czhao obtain time information for hourly integration
 type(iAMAS_Time_Type):: currTime
 character(len=StrKIND)   :: timeStamp
 integer:: ierr
 real(kind=RKIND):: utc_h
 integer:: yr,h,m,s,s_n,s_d,DoY 

!----------------------------------------------------------------------------------------------------
  !czhao obtain time for hourly integration 
!  call iamas_log_write('czhao check time clock at the begin')

   currTime = iamas_get_clock_time(clock,iAMAS_NOW,ierr)
   call iamas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,H=h,M=m, &
                    S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
   utc_h = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0

! call iamas_log_write('czhao check h  = $i', intArgs=(/h/))
! call iamas_log_write('czhao check m  = $i', intArgs=(/m/))
! call iamas_log_write('czhao check utc_h=$r', realArgs=(/utc_h/))
  
 call iamas_pool_get_config(configs,'config_microp_scheme',microp_scheme       )
 call iamas_pool_get_config(configs,'config_bucket_rainnc',config_bucket_rainnc)
 call iamas_pool_get_config(configs,'config_print_debug_info' ,debug_on )
 call iamas_pool_get_array(diag_physics,'i_rainnc'  ,i_rainnc  )
 call iamas_pool_get_array(diag_physics,'precipw'   ,precipw   )
 call iamas_pool_get_array(diag_physics,'rainnc'    ,rainnc    )
 call iamas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call iamas_pool_get_array(diag_physics,'rainnc_hr' ,rainnc_hr )
 call iamas_pool_get_array(diag_physics,'sr'        ,sr        )

 do i = its,ite
    precipw(i) = 0._RKIND    
 enddo

!variables common to all cloud microphysics schemes:
 do j = jts,jte
 do i = its,ite

    !precipitable water:
    do k = kts,kte
       rho_a = rho_p(i,k,j) / (1._RKIND + qv_p(i,k,j))
       precipw(i) = precipw(i) + qv_p(i,k,j) * rho_a * dz_p(i,k,j)
    enddo


    !time-step precipitation:
    rainncv(i) = rainnc_p(i,j)
    
    !accumulated precipitation:
    rainnc(i) = rainnc(i) + rainncv(i)

    !czhao obtain hourly rainnc
    if (utc_h .eq. h) then 
    rainnc_hr(i) = 0.0
    endif
    rainnc_hr(i) = rainnc_hr(i) + rainncv(i)


   if(l_acrain .and. config_bucket_rainnc.gt.0._RKIND .and. &
      rainnc(i).gt.config_bucket_rainnc) then
      i_rainnc(i) = i_rainnc(i) + 1
      rainnc(i)   = rainnc(i) - config_bucket_rainnc
   endif
 
 enddo
 enddo

!variables specific to different cloud microphysics schemes:
 microp_select_init: select case(microp_scheme)

    case ("mp_thompson","mp_wsm6","mp_morrison")
       call iamas_pool_get_array(diag_physics,'graupelnc' ,graupelnc )
       call iamas_pool_get_array(diag_physics,'graupelncv',graupelncv)
       call iamas_pool_get_array(diag_physics,'snownc'    ,snownc    )
       call iamas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
       do j = jts,jte
       do i = its,ite
          !time-step precipitation:
          snowncv(i)    = snownc_p(i,j)
          graupelncv(i) = graupelnc_p(i,j)
          sr(i) = (snownc_p(i,j) + graupelnc_p(i,j)) / (rainnc_p(i,j)+1.e-12)

          !accumulated precipitation:
          snownc(i)    = snownc(i) + snowncv(i)
          graupelnc(i) = graupelnc(i) + graupelncv(i)
       enddo
       enddo


    microp2_select_init: select case(microp_scheme)

       case("mp_thompson","mp_morrison")
          !jwfeng: 03/22/2020
          !cloud water/ice sink:
          call iamas_pool_get_array(diag_physics,'qlsink'    ,qlsink    )
          call iamas_pool_get_array(diag_physics,'ilsink'    ,ilsink    )
          call iamas_pool_get_array(diag_physics,'precr'     ,precr     )
          call iamas_pool_get_array(diag_physics,'preci'     ,preci     )
          call iamas_pool_get_array(diag_physics,'precs'     ,precs     )
          call iamas_pool_get_array(diag_physics,'precg'     ,precg     )
          do j = jts,jte
          do i = its,ite
          do k = kts,kte
             qlsink(k,i) = qlsink_p(i,k,j)
             ilsink(k,i) = ilsink_p(i,k,j)
              precr(k,i) =  precr_p(i,k,j)
              preci(k,i) =  preci_p(i,k,j)
              precs(k,i) =  precs_p(i,k,j)
              precg(k,i) =  precg_p(i,k,j)
          enddo
          enddo
          enddo
          if(debug_on) call iamas_log_write( 'max/min of qlsink = $r/$r',realArgs=(/maxval(qlsink(kts:kte,its:ite)),minval(qlsink(kts:kte,its:ite))/))
          if(debug_on) call iamas_log_write( 'max/min of ilsink = $r/$r',realArgs=(/maxval(ilsink(kts:kte,its:ite)),minval(ilsink(kts:kte,its:ite))/))
          if(debug_on) call iamas_log_write( 'max/min of precr  = $r/$r',realArgs=(/maxval(precr(kts:kte,its:ite)),minval(precr(kts:kte,its:ite))  /))
          if(debug_on) call iamas_log_write( 'max/min of preci  = $r/$r',realArgs=(/maxval(preci(kts:kte,its:ite)),minval(preci(kts:kte,its:ite))  /))
          if(debug_on) call iamas_log_write( 'max/min of precs  = $r/$r',realArgs=(/maxval(precs(kts:kte,its:ite)),minval(precs(kts:kte,its:ite))  /))
          if(debug_on) call iamas_log_write( 'max/min of precg  = $r/$r',realArgs=(/maxval(precg(kts:kte,its:ite)),minval(precg(kts:kte,its:ite))  /))
          if(debug_on) call iamas_log_write( 'its/ite, jts/jte, kts/kte = $i/$i, $i/$i, $i/$i', intArgs=(/its,ite,jts,jte,kts,kte/))
          if(microp_scheme == "mp_thompson") then
             call iamas_pool_get_array(diag_physics,'pnc_wcd_3d',pnc_wcd_3d)
             call iamas_pool_get_array(diag_physics,'pnc_wau_3d',pnc_wau_3d)
             call iamas_pool_get_array(diag_physics,'pnc_rcw_3d',pnc_rcw_3d)
             call iamas_pool_get_array(diag_physics,'pni_wfz_3d',pni_wfz_3d)
             call iamas_pool_get_array(diag_physics,'pnc_scw_3d',pnc_scw_3d)
             call iamas_pool_get_array(diag_physics,'pnc_gcw_3d',pnc_gcw_3d)
             do j = jts,jte
             do i = its,ite
             do k = kts,kte
                 pnc_wcd_3d(k,i) = pnc_wcd_3d_p(i,k,j)
                 pnc_wau_3d(k,i) = pnc_wau_3d_p(i,k,j)
                 pnc_rcw_3d(k,i) = pnc_rcw_3d_p(i,k,j)
                 pni_wfz_3d(k,i) = pni_wfz_3d_p(i,k,j)
                 pnc_scw_3d(k,i) = pnc_scw_3d_p(i,k,j)
                 pnc_gcw_3d(k,i) = pnc_gcw_3d_p(i,k,j)
             enddo
             enddo
             enddo
          endif

       case default

    end select microp2_select_init

    case ("mp_p3")
       call iamas_pool_get_array(diag_physics,'snownc' ,snownc )
       call iamas_pool_get_array(diag_physics,'snowncv',snowncv)
       do j = jts,jte
       do i = its,ite
          !time-step precipitation:
          snowncv(i)    = snownc_p(i,j)
          sr(i) = snownc_p(i,j) / (rainnc_p(i,j)+1.e-12)

          !accumulated precipitation:
          snownc(i)    = snownc(i) + snowncv(i)
       enddo
       enddo

    case default

 end select microp_select_init

 end subroutine precip_to_iAMAS

!=================================================================================================================
 subroutine compute_radar_reflectivity(configs,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 real(kind=RKIND),dimension(:),pointer:: refl10cm_max,refl10cm_1km,refl10cm_1km_max

!local variables and arrays:
 integer:: i,j,k,kp
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qc1d,qr1d,qs1d,qg1d,t1d,p1d,nr1d,dBZ1d,zp
 real(kind=RKIND),dimension(:),allocatable:: ns1d,ng1d !added for Morrison
 real(kind=RKIND):: w1,w2

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call iamas_pool_get_array(diag_physics,'refl10cm_max',refl10cm_max)
 call iamas_pool_get_array(diag_physics,'refl10cm_1km',refl10cm_1km)
 call iamas_pool_get_array(diag_physics,'refl10cm_1km_max',refl10cm_1km_max)

 microp_select: select case(microp_scheme)

    case ("mp_kessler")
       call physics_error_fatal('--- calculation of radar reflectivity is not available' // &
                                 'with kessler cloud microphysics')

    case ("mp_wsm6")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j)+0.5*dz_p(i,1,j) ! height AGL
          enddo

          call refl10cm_wsm6(qv1d,qr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp)   ) deallocate(zp   )

    case ("mp_thompson")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qc1d) ) allocate(qc1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(nr1d) ) allocate(nr1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qc1d(k)  = qc_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             nr1d(k)  = nr_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j)+0.5*dz_p(i,1,j) ! height AGL
          enddo

          call calc_refl10cm(qv1d,qc1d,qr1d,nr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qc1d) ) deallocate(qc1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(nr1d) ) deallocate(nr1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp)   ) deallocate(zp   )

    case ("mp_morrison")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       !if(.not.allocated(qc1d) ) allocate(qc1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(nr1d) ) allocate(nr1d(kts:kte) )
       if(.not.allocated(ns1d) ) allocate(ns1d(kts:kte) )
       if(.not.allocated(ng1d) ) allocate(ng1d(kts:kte) )
       if(.not.allocated(dBZ1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             !qc1d(k)  = qc_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             nr1d(k)  = nr_p(i,k,j)
             ns1d(k)  = ns_p(i,k,j)
             ng1d(k)  = ng_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j)+0.5*dz_p(i,1,j) ! height AGL
          enddo

          !call refl10cm_hm(qv1d,qc1d,qr1d,nr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)
          call refl10cm_hm(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       !if(allocated(qc1d) ) deallocate(qc1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(nr1d) ) deallocate(nr1d )
       if(allocated(ns1d) ) deallocate(ns1d )
       if(allocated(ng1d) ) deallocate(ng1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp)   ) deallocate(zp   )

    case default

 end select microp_select

 end subroutine compute_radar_reflectivity

!=================================================================================================================
 subroutine compute_relhum(diag,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag

!local pointers:
 real(kind=RKIND),dimension(:,:),pointer:: relhum

!local variables and arrays:
 integer:: i,j,k

 real(kind=RKIND):: tempc, rh
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qvs1d,t1d,p1d

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_array(diag,'relhum',relhum)
 do k = kts,kte
    do i = its,ite
       relhum(k,i) = 0._RKIND
    enddo
 enddo

 if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
 if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
 if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
 if(.not.allocated(qvs1d)) allocate(qvs1d(kts:kte))

 do j = jts,jte
 do i = its,ite
    do k = kts,kte
       p1d(k)  = pres_p(i,k,j)
       t1d(k)  = th_p(i,k,j) * pi_p(i,k,j)
       tempc   = t1d(k) - 273.16_RKIND
       qvs1d(k) = rslf(p1d(k),t1d(k))
       if(tempc .le. 0._RKIND) qvs1d(k) = rsif(p1d(k),t1d(k))
       qv1d(k) = qv_p(i,k,j)
       relhum(k,i) = qv1d(k) / qvs1d(k)
       relhum(k,i) = relhum(k,i) * 100._RKIND
    enddo
 enddo
 enddo

 if(allocated(p1d)  ) deallocate(p1d  )
 if(allocated(t1d)  ) deallocate(t1d  )
 if(allocated(qv1d) ) deallocate(qv1d )
 if(allocated(qvs1d)) deallocate(qvs1d)

 end subroutine compute_relhum

!=================================================================================================================
 end module iamas_atmphys_driver_microphysics
!=================================================================================================================
