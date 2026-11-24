! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_driver_convection
 use iamas_kind_types
 use iamas_pool_routines
 use iamas_timer, only : iamas_timer_start, iamas_timer_stop

 use iamas_atmphys_constants
 use iamas_atmphys_utilities
 use iamas_atmphys_vars

!wrf physics:
 use module_cu_gf
 use module_cu_kfeta
 use module_cu_mskf
 use module_cu_tiedtke
 use module_cu_ntiedtke

 implicit none
 private
 public:: allocate_convection,     &
          deallocate_convection,   &
          init_convection,         &
          driver_convection,       &
          update_convection_step1, &
          update_convection_step2


!iAMAS driver for parameterization of convection.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in iamas_atmphys_driver_convection:
! ----------------------------------------------
! allocate_convection     : allocate local arrays for parameterization of convection.
! deallocate_convection   : deallocate local arrays for parameterization of convection.
! init_convection         : initialization of individual convection scheme.
! driver_convection       : main driver (called from subroutine physics_driver).
! convection_from_iAMAS    : initialize local arrays.
! convection_to_iAMAS      : copy local arrays to iAMAS arrays.
! update_convection_step1 : updates lifetime of deep convective clouds in Kain-Fritsch scheme.
! update_convection_step2 : updates accumulated precipitation output from convection schemes.
!
! WRF physics called from driver_convection:
! ------------------------------------------
! * module_cu_kfeta   : Kain-Fritsch convection scheme.
! * module_cu_tiedtke : Tiedtke convection scheme.

! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to the subroutines kf_eta_cps
!   and tiedtke. 
! * removed call to the updated Kain-Fritsch convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added the mean distance between cell centers in the call to the Kain-Fritsch parameterization of convection.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-22. 
! * in call to subroutine kf_eta_cps, replaced the variable g (that originally pointed to gravity) with gravity,
!   for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * Modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * in ./physics_wrf, updated the Tiedtke convection scheme to that in WRF version 3.6.1. The call to subroutine 
!   cu_tiedtke has been updated accordingly to include the sensible heat flux.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-11.
! * renamed "kain_fritsch" with "cu_kain_fritsch" and "tiedtke" with "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * added the implementation of the Grell-Freitas convection scheme (option cu_grell_freitas).
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
!   code cleaner, now that we are using physics packages.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-01.
! * in the call to subroutine cu_grell_freitas, added the array pratec for the calculation of rainc.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-14.
! * in the call to subroutine cu_grell_freitas, removed the argument ktop_deep.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-20.
! * updated the call to the Tiedtke parameterization of convection to WRF version 3.8.1. This option is
!   triggered using the option cu_tiedtke.
!   Laura D. Fowler (laura@ucar.edu) / 2016-08-18.
! * added the call to the "new" Tiedtke parameterization of convection from WRF version 3.8.1. This option is
!   triggered using the option cu_ntiedtke.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-20.
! * for the kain_fritsch parameterization of convection, change the definition of dx_p to match that used in the
!   Grell-Freitas and "new Tiedtke" parameterization.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.
! * since we removed the local variable convection_scheme from iamas_atmphys_vars.F, now defines convection_scheme
!   as a pointer to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2917-02-16.


 contains


!=================================================================================================================
 subroutine allocate_convection(configs)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs

!local variables and pointers:
 character(len=StrKIND),pointer:: convection_scheme

 integer:: i,k,j

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(.not.allocated(cu_act_flag)      ) allocate(cu_act_flag(ims:ime,jms:jme)             )
 if(.not.allocated(rthcuten_p)       ) allocate(rthcuten_p(ims:ime,kms:kme,jms:jme)      )
 if(.not.allocated(rqvcuten_p)       ) allocate(rqvcuten_p(ims:ime,kms:kme,jms:jme)      )
 if(.not.allocated(rqccuten_p)       ) allocate(rqccuten_p(ims:ime,kms:kme,jms:jme)      )
 if(.not.allocated(rqicuten_p)       ) allocate(rqicuten_p(ims:ime,kms:kme,jms:jme)      )
 if(.not.allocated(rthcuten_shalc_p) ) allocate(rthcuten_shalc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvcuten_shalc_p) ) allocate(rqvcuten_shalc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqccuten_shalc_p) ) allocate(rqccuten_shalc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqicuten_shalc_p) ) allocate(rqicuten_shalc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(pratec_p)         ) allocate(pratec_p(ims:ime,jms:jme)                )
 if(.not.allocated(raincv_p)         ) allocate(raincv_p(ims:ime,jms:jme)                )

 do i = its,ite
 do j = jts,jte
    pratec_p(i,j) = 0._RKIND
    raincv_p(i,j) = 0._RKIND
 enddo
 enddo

 do i = its,ite
 do k = kts,kte
 do j = jts,jte
    rthcuten_p(i,k,j) = 0._RKIND
    rqvcuten_p(i,k,j) = 0._RKIND
    rqccuten_p(i,k,j) = 0._RKIND
    rqicuten_p(i,k,j) = 0._RKIND
    rthcuten_shalc_p(i,k,j) = 0._RKIND
    rqvcuten_shalc_p(i,k,j) = 0._RKIND
    rqccuten_shalc_p(i,k,j) = 0._RKIND
    rqicuten_shalc_p(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas", "cu_kain_fritsch", "cu_kain_fritsch_sca")
       if(.not.allocated(cubot_p)        ) allocate(cubot_p(ims:ime,jms:jme)            )
       if(.not.allocated(cutop_p)        ) allocate(cutop_p(ims:ime,jms:jme)            )

       if(.not.allocated(qccu_p)         ) allocate(qccu_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(qicu_p)         ) allocate(qicu_p(ims:ime,kms:kme,jms:jme)      )
       if(.not.allocated(qccu_shalc_p)   ) allocate(qccu_shalc_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qicu_shalc_p)   ) allocate(qicu_shalc_p(ims:ime,kms:kme,jms:jme))

       if(.not.allocated(rinc_deepc_p)   ) allocate(rinc_deepc_p(ims:ime,        jms:jme))
       if(.not.allocated(rinc_shalc_p)   ) allocate(rinc_shalc_p(ims:ime,        jms:jme))

       if(.not.allocated(kbcon_deepc_p)  ) allocate(kbcon_deepc_p(ims:ime,jms:jme)      )
       if(.not.allocated(kbcon_shalc_p)  ) allocate(kbcon_shalc_p(ims:ime,jms:jme)      )
       if(.not.allocated(ktop_deepc_p)   ) allocate(ktop_deepc_p(ims:ime,jms:jme)       )
       if(.not.allocated(ktop_shalc_p)   ) allocate(ktop_shalc_p(ims:ime,jms:jme)       )
       if(.not.allocated(xmb_deepc_p)    ) allocate(xmb_deepc_p(ims:ime,jms:jme)        )
       if(.not.allocated(xmb_shalc_p)    ) allocate(xmb_shalc_p(ims:ime,jms:jme)        )

       cu_kfgf_select: select case(convection_scheme)
          case ("cu_grell_freitas")
             if(.not.allocated(k22_deepc_p)    ) allocate(k22_deepc_p(ims:ime,jms:jme)        )
             if(.not.allocated(k22_shalc_p)    ) allocate(k22_shalc_p(ims:ime,jms:jme)        )
             if(.not.allocated(ht_p)           ) allocate(ht_p(ims:ime,jms:jme)               )
             if(.not.allocated(kpbl_p)         ) allocate(kpbl_p(ims:ime,jms:jme)             )

             if(.not.allocated(dx_p)           ) allocate(dx_p(ims:ime,jms:jme)               )
             if(.not.allocated(area_p)         ) allocate(area_p(ims:ime,jms:jme)             )
             if(.not.allocated(gsw_p)          ) allocate(gsw_p(ims:ime,jms:jme)              )
             if(.not.allocated(hfx_p)          ) allocate(hfx_p(ims:ime,jms:jme)              )
             if(.not.allocated(qfx_p)          ) allocate(qfx_p(ims:ime,jms:jme)              )
             if(.not.allocated(xland_p)        ) allocate(xland_p(ims:ime,jms:jme)            )

             if(.not.allocated(rthblten_p)     ) allocate(rthblten_p(ims:ime,kms:kme,jms:jme) )
             if(.not.allocated(rthdynten_p)    ) allocate(rthdynten_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(rthraten_p)     ) allocate(rthraten_p(ims:ime,kms:kme,jms:jme) )
             if(.not.allocated(rqvblten_p)     ) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme) )
             if(.not.allocated(rqvdynten_p)    ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(rucuten_p)      ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)  )
             if(.not.allocated(rvcuten_p)      ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)  )

          case ("cu_kain_fritsch", "cu_kain_fritsch_sca")
             if(.not.allocated(dx_p)             ) allocate(dx_p(ims:ime,jms:jme)                    )
             if(.not.allocated(area_p)           ) allocate(area_p(ims:ime,jms:jme)                  )
             if(.not.allocated(nca_p)            ) allocate(nca_p(ims:ime,jms:jme)                   )
             if(.not.allocated(timec_p)          ) allocate(timec_p(ims:ime,jms:jme)                 )
             if(.not.allocated(w0avg_p)          ) allocate(w0avg_p(ims:ime,kms:kme,jms:jme)         )
             if(.not.allocated(rqrcuten_p)       ) allocate(rqrcuten_p(ims:ime,kms:kme,jms:jme)      )
             if(.not.allocated(rqscuten_p)       ) allocate(rqscuten_p(ims:ime,kms:kme,jms:jme)      )
             if(.not.allocated(rqrcuten_shalc_p) ) allocate(rqrcuten_shalc_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(rqscuten_shalc_p) ) allocate(rqscuten_shalc_p(ims:ime,kms:kme,jms:jme))

             if(.not.allocated(cldfrac_deepc_p)  ) allocate(cldfrac_deepc_p(ims:ime,kms:kme,jms:jme) )
             if(.not.allocated(cldfrac_shalc_p)  ) allocate(cldfrac_shalc_p(ims:ime,kms:kme,jms:jme) )
             !jwfeng:2021-06-03
             if(.not.allocated(umf_p)            ) allocate(umf_p(ims:ime,kms:kme,jms:jme)           )
             if(.not.allocated(uer_p)            ) allocate(uer_p(ims:ime,kms:kme,jms:jme)           )
             if(.not.allocated(udr_p)            ) allocate(udr_p(ims:ime,kms:kme,jms:jme)           )
             if(.not.allocated(dmf_p)            ) allocate(dmf_p(ims:ime,kms:kme,jms:jme)           )
             if(.not.allocated(der_p)            ) allocate(der_p(ims:ime,kms:kme,jms:jme)           )
             if(.not.allocated(ddr_p)            ) allocate(ddr_p(ims:ime,kms:kme,jms:jme)           )
             if(.not.allocated( wu_p)            ) allocate( wu_p(ims:ime,kms:kme,jms:jme)           )
             if(.not.allocated(ldb_p)            ) allocate(ldb_p(ims:ime,        jms:jme)           )
             if(.not.allocated(let_p)            ) allocate(let_p(ims:ime,        jms:jme)           )
             if(.not.allocated(ksl_p)            ) allocate(ksl_p(ims:ime,        jms:jme)           )
             if(.not.allocated(lfs_p)            ) allocate(lfs_p(ims:ime,        jms:jme)           )
             if(.not.allocated(klc_p)            ) allocate(klc_p(ims:ime,        jms:jme)           )
             if(.not.allocated(ishall_p)         ) allocate(ishall_p(ims:ime,     jms:jme)           )
             if(.not.allocated(check_kfeta_p)    ) allocate(check_kfeta_p(ims:ime,jms:jme)           )

             do i = its,ite
             do j = jts,jte
                cubot_p(i,j) = DBLE(kte+1)
                cutop_p(i,j) = DBLE(kts)
             enddo
             enddo

             do i = its,ite
             do k = kts,kte
             do j = jts,jte
                rqrcuten_p(i,k,j) = 0._RKIND
                rqscuten_p(i,k,j) = 0._RKIND
                rqrcuten_shalc_p(i,k,j) = 0._RKIND
                rqscuten_shalc_p(i,k,j) = 0._RKIND
             enddo
             enddo
             enddo

             cu_kf_select: select case(convection_scheme)
                case ("cu_kain_fritsch_sca")
                   if(.not.allocated(hpbl_p)   ) allocate(hpbl_p(ims:ime,jms:jme)           )
                   if(.not.allocated(ust_p)    ) allocate(ust_p(ims:ime,jms:jme)            )
                   if(.not.allocated(wstar_p)  ) allocate(wstar_p(ims:ime,jms:jme)          )
                   if(.not.allocated(zol_p)    ) allocate(zol_p(ims:ime,jms:jme)            )
                   if(.not.allocated(w_deepc_p)) allocate(w_deepc_p(ims:ime,kms:kme,jms:jme))
                case default
             end select cu_kf_select

          case default

       end select cu_kfgf_select

    case ("cu_tiedtke","cu_ntiedtke")
       if(.not.allocated(hfx_p)        ) allocate(hfx_p(ims:ime,jms:jme)                )
       if(.not.allocated(qfx_p)        ) allocate(qfx_p(ims:ime,jms:jme)                )
       if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)              )
       if(.not.allocated(rucuten_p)    ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(rvcuten_p)    ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)    )

       do i = its,ite
       do j = jts,jte
          hfx_p(i,j)   = 0._RKIND
          qfx_p(i,j)   = 0._RKIND
          xland_p(i,j) = 0._RKIND
       enddo
       enddo

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rucuten_p(i,k,j)     = 0._RKIND
          rvcuten_p(i,k,j)     = 0._RKIND
       enddo
       enddo
       enddo

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             if(.not.allocated(rqvdynten_p)  ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme)  )
             if(.not.allocated(rqvdynblten_p)) allocate(rqvdynblten_p(ims:ime,kms:kme,jms:jme))

             do i = its,ite
             do k = kts,kte
             do j = jts,jte
                rqvdynten_p(i,k,j)   = 0._RKIND
                rqvdynblten_p(i,k,j) = 0._RKIND
             enddo
             enddo
             enddo

          case("cu_ntiedtke")
             if(.not.allocated(dx_p)     ) allocate(dx_p(ims:ime,jms:jme)             )
             if(.not.allocated(rqvften_p)) allocate(rqvften_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(rthften_p)) allocate(rthften_p(ims:ime,kms:kme,jms:jme))

             do i = its,ite
             do k = kts,kte
             do j = jts,jte
                rqvften_p(i,k,j) = 0._RKIND
                rthften_p(i,k,j) = 0._RKIND
             enddo
             enddo
             enddo

          case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine allocate_convection

!=================================================================================================================
 subroutine deallocate_convection(configs)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(allocated(cu_act_flag)     ) deallocate(cu_act_flag     )
 if(allocated(rthcuten_p)      ) deallocate(rthcuten_p      )
 if(allocated(rqvcuten_p)      ) deallocate(rqvcuten_p      )
 if(allocated(rqccuten_p)      ) deallocate(rqccuten_p      )
 if(allocated(rqicuten_p)      ) deallocate(rqicuten_p      )
 if(allocated(rthcuten_shalc_p)) deallocate(rthcuten_shalc_p)
 if(allocated(rqvcuten_shalc_p)) deallocate(rqvcuten_shalc_p)
 if(allocated(rqccuten_shalc_p)) deallocate(rqccuten_shalc_p)
 if(allocated(rqicuten_shalc_p)) deallocate(rqicuten_shalc_p)
 if(allocated(pratec_p)        ) deallocate(pratec_p        )
 if(allocated(raincv_p)        ) deallocate(raincv_p        )

 convection_select: select case(convection_scheme)
    case ("cu_grell_freitas", "cu_kain_fritsch", "cu_kain_fritsch_sca")
       if(allocated(cubot_p)        ) deallocate(cubot_p        )
       if(allocated(cutop_p)        ) deallocate(cutop_p        )

       if(allocated(qccu_p)         ) deallocate(qccu_p         )
       if(allocated(qicu_p)         ) deallocate(qicu_p         )
       if(allocated(qccu_shalc_p)   ) deallocate(qccu_shalc_p   )
       if(allocated(qicu_shalc_p)   ) deallocate(qicu_shalc_p   )

       if(allocated(rinc_deepc_p)   ) deallocate(rinc_deepc_p   )
       if(allocated(rinc_shalc_p)   ) deallocate(rinc_shalc_p   )

       if(allocated(kbcon_deepc_p)  ) deallocate(kbcon_deepc_p  )
       if(allocated(kbcon_shalc_p)  ) deallocate(kbcon_shalc_p  )
       if(allocated(ktop_deepc_p)   ) deallocate(ktop_deepc_p   )
       if(allocated(ktop_shalc_p)   ) deallocate(ktop_shalc_p   )
       if(allocated(xmb_deepc_p)    ) deallocate(xmb_deepc_p    )
       if(allocated(xmb_shalc_p)    ) deallocate(xmb_shalc_p    )

       cu_kfgf_select: select case(convection_scheme)
          case ("cu_grell_freitas")
             if(allocated(ht_p)           ) deallocate(ht_p           )
             if(allocated(kpbl_p)         ) deallocate(kpbl_p         )
             if(allocated(k22_deepc_p)    ) deallocate(k22_deepc_p    )
             if(allocated(k22_shalc_p)    ) deallocate(k22_shalc_p    )

             if(allocated(dx_p)           ) deallocate(dx_p           )
             if(allocated(area_p)         ) deallocate(area_p         )
             if(allocated(gsw_p)          ) deallocate(gsw_p          )
             if(allocated(hfx_p)          ) deallocate(hfx_p          )
             if(allocated(qfx_p)          ) deallocate(qfx_p          )
             if(allocated(xland_p)        ) deallocate(xland_p        )
             if(allocated(rthblten_p)     ) deallocate(rthblten_p     )
             if(allocated(rthdynten_p)    ) deallocate(rthdynten_p    )
             if(allocated(rthraten_p)     ) deallocate(rthraten_p     )
             if(allocated(rqvblten_p)     ) deallocate(rqvblten_p     )
             if(allocated(rqvdynten_p)    ) deallocate(rqvdynten_p    )
             if(allocated(rucuten_p)      ) deallocate(rucuten_p      )
             if(allocated(rvcuten_p)      ) deallocate(rvcuten_p      )

          case ("cu_kain_fritsch", "cu_kain_fritsch_sca")
             if(allocated(dx_p)               ) deallocate(dx_p               )
             if(allocated(area_p)             ) deallocate(area_p             )
             if(allocated(nca_p)              ) deallocate(nca_p              )
             if(allocated(timec_p)            ) deallocate(timec_p            )
             if(allocated(w0avg_p)            ) deallocate(w0avg_p            )
             if(allocated(rqrcuten_p)         ) deallocate(rqrcuten_p         )
             if(allocated(rqscuten_p)         ) deallocate(rqscuten_p         )
             if(allocated(rqrcuten_shalc_p)   ) deallocate(rqrcuten_shalc_p   )
             if(allocated(rqscuten_shalc_p)   ) deallocate(rqscuten_shalc_p   )

             if(allocated(cldfrac_deepc_p)    ) deallocate(cldfrac_deepc_p    )
             if(allocated(cldfrac_shalc_p)    ) deallocate(cldfrac_shalc_p    )
             !jwfeng:2021-06-03
             if(allocated(umf_p)        ) deallocate(umf_p        )
             if(allocated(uer_p)        ) deallocate(uer_p        )
             if(allocated(udr_p)        ) deallocate(udr_p        )
             if(allocated(dmf_p)        ) deallocate(dmf_p        )
             if(allocated(der_p)        ) deallocate(der_p        )
             if(allocated(ddr_p)        ) deallocate(ddr_p        )
             if(allocated( wu_p)        ) deallocate( wu_p        )
             if(allocated(ldb_p)        ) deallocate(ldb_p        )
             if(allocated(let_p)        ) deallocate(let_p        )
             if(allocated(ksl_p)        ) deallocate(ksl_p        )
             if(allocated(lfs_p)        ) deallocate(lfs_p        )
             if(allocated(klc_p)        ) deallocate(klc_p        )
             if(allocated(ishall_p)     ) deallocate(ishall_p     )
             if(allocated(check_kfeta_p)) deallocate(check_kfeta_p)

             cu_kf_select: select case(convection_scheme)
                case ("cu_kain_fritsch_sca")
                   if(allocated(hpbl_p)         ) deallocate(hpbl_p         )
                   if(allocated(ust_p)          ) deallocate(ust_p          )
                   if(allocated(wstar_p)        ) deallocate(wstar_p        )
                   if(allocated(zol_p)          ) deallocate(zol_p          )
                   if(allocated(w_deepc_p)      ) deallocate(w_deepc_p      )
                case default
             end select cu_kf_select

       end select cu_kfgf_select

    case ("cu_tiedtke","cu_ntiedtke")
       if(allocated(hfx_p)        ) deallocate(hfx_p        )
       if(allocated(qfx_p)        ) deallocate(qfx_p        )
       if(allocated(xland_p)      ) deallocate(xland_p      )
       if(allocated(rucuten_p)    ) deallocate(rucuten_p    )
       if(allocated(rvcuten_p)    ) deallocate(rvcuten_p    )

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             if(allocated(rqvdynten_p)  ) deallocate(rqvdynten_p  )
             if(allocated(rqvdynblten_p)) deallocate(rqvdynblten_p)

          case ("cu_ntiedtke")
             if(allocated(dx_p)     ) deallocate(dx_p)
             if(allocated(rqvften_p)) deallocate(rqvften_p)
             if(allocated(rthften_p)) deallocate(rthften_p)

          case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine deallocate_convection

!=================================================================================================================
 subroutine init_convection(mesh,configs,diag_physics)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics

!local pointers:
 logical,pointer:: config_do_restart
 character(len=StrKIND),pointer:: convection_scheme

 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: nca

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_dimension(mesh,'nCells',nCells)
 call iamas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
 call iamas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 convection_select: select case(convection_scheme)

    case ("cu_kain_fritsch")
       call iamas_pool_get_array(diag_physics,'nca',nca)
       if(.not. config_do_restart) then
          do iCell = 1, nCells
             nca(iCell) = -100._RKIND
          enddo
       endif
       call kf_lutab(svp1,svp2,svp3,svpt0)

    case ("cu_kain_fritsch_sca")
       call iamas_pool_get_array(diag_physics,'nca',nca)
       if(.not. config_do_restart) then
          do iCell = 1, nCells
             nca(iCell) = -100._RKIND
          enddo
       endif
       call mskf_lutab(svp1,svp2,svp3,svpt0)

    case default

 end select convection_select

 end subroutine init_convection

!=================================================================================================================
 subroutine driver_convection(itimestep,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: sfc_input

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: tend_physics

!local variables and arrays:
 logical:: log_convection

 integer:: i,j
 integer:: icount,initflag

 real(kind=RKIND):: dx

!local pointers:
 logical,pointer:: config_do_restart
 integer,pointer:: gfconv_closure_deep,gfconv_closure_shallow
 integer,pointer:: nClosures
 character(len=StrKIND),pointer:: convection_scheme
 real(kind=RKIND),pointer:: len_disp

!variables specific to Kain_Fritsch parameterization:
 logical:: warm_rain,adapt_step_flag
 integer:: ktau
 real(kind=RKIND):: curr_secs
 real(kind=RKIND):: cudt
 real(kind=RKIND):: cudtacttime

!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write('')
!call iamas_log_write('--- enter subroutine driver_convection:')

 call iamas_pool_get_config(configs,'config_gfconv_closure_deep',gfconv_closure_deep)
 call iamas_pool_get_config(configs,'config_gfconv_closure_shallow',gfconv_closure_shallow)
 call iamas_pool_get_config(configs,'config_len_disp'         ,len_disp         )
 call iamas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
 call iamas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

!initialize instantaneous precipitation, and copy convective tendencies from the dynamics to
!the physics grid:
 call convection_from_iAMAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)

!... convert the convection time-step to minutes:
 cudt = dt_cu/60.

!... call to convection schemes:
 curr_secs   = -1
 cudtacttime = -1
 adapt_step_flag = .false.
 do j = jts, jte
 do i = its, ite
    cu_act_flag(i,j) = .false.
 enddo
 enddo

!... initialization of initflag needed in the calls cu_tiedtke and cu_ntiedtke:
 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas")
       call iamas_timer_start('Grell-Freitas')
       call iamas_pool_get_dimension(mesh,'nClosures',nClosures)

       call cu_grell_freitas( &
             itimestep       = itimestep              , dt            = dt_dyn          , &
             dxCell          = dx_p                   , areaCell      = area_p          , &
             u               = u_p                    , v             = v_p             , &
             w               = w_p                    , t             = t_p             , &
             q               = qv_p                   , rho           = rho_p           , &
             p               = pres_hyd_p             , pi            = pi_p            , &
             p8w             = pres2_hyd_p            , dz8w          = dz_p            , &
             ht              = ht_p                   , xland         = xland_p         , &
             gsw             = gsw_p                  , xlv           = xlv             , &
             cp              = cp                     , g             = gravity         , &
             r_v             = R_v                    , hfx           = hfx_p           , &
             qfx             = qfx_p                  , rthblten      = rthblten_p      , &
             rqvblten        = rqvblten_p             , rthften       = rthdynten_p     , &
             rqvften         = rqvdynten_p            , rthraten      = rthraten_p      , &
             kpbl            = kpbl_p                 , raincv        = raincv_p        , &
             pratec          = pratec_p               , htop          = cutop_p         , &
             hbot            = cubot_p                , k22_deepc     = k22_deepc_p     , &
             kbcon_deepc     = kbcon_deepc_p          , ktop_deepc    = ktop_deepc_p    , &
             k22_shalc       = k22_shalc_p            , kbcon_shalc   = kbcon_shalc_p   , &
             ktop_shalc      = ktop_shalc_p           , xmb_total     = xmb_deepc_p     , &
             xmb_shallow     = xmb_shalc_p            , gdc           = qccu_p          , &
             gdc2            = qicu_p                 , rthcuten      = rthcuten_p      , &
             rqvcuten        = rqvcuten_p             , rqccuten      = rqccuten_p      , &
             rqicuten        = rqicuten_p             , gdcs          = qccu_shalc_p    , &
             gdc2s           = qicu_shalc_p           , rthcutens     = rthcuten_shalc_p, &
             rqvcutens       = rqvcuten_shalc_p       , rqccutens     = rqccuten_shalc_p, &
             rqicutens       = rqicuten_shalc_p       , rinc          = rinc_deepc_p    , &
             rincs           = rinc_shalc_p           , nclos         = nClosures       , &
             ichoice_deep    = gfconv_closure_deep                                      , &
             ichoice_shallow = gfconv_closure_shallow                                   , &
             ishallow_g3     = ishallow                                                 , &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde      , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme      , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte        &
                            )
       call iamas_timer_stop('Grell-Freitas')

    case ("cu_kain_fritsch")
       if(itimestep == 1) then
          ktau = itimestep
       else
          ktau = itimestep + 1
       endif

       call iamas_timer_start('Kain-Fritsch')
       call  kf_eta_cps ( &
             pcps            = pres_hyd_p      , t         = t_p          , &
             dt              = dt_dyn          , ktau      = ktau         , &
             dxCell          = dx_p            , areaCell  = area_p       , &
             cudt            = cudt            , curr_secs = curr_secs    , & 
             adapt_step_flag = adapt_step_flag , rho       = rho_p        , &
             raincv          = raincv_p        , pratec    = pratec_p     , &
             nca             = nca_p           , u         = u_p          , &
             v               = v_p             , th        = th_p         , &
             pi              = pi_p            , w         = w_p          , &
             dz8w            = dz_p            , w0avg     = w0avg_p      , &
             xlv0            = xlv0            , xlv1      = xlv1         , &
             xls0            = xls0            , xls1      = xls1         , &
             cp              = cp              , r         = r_d          , &
             g               = gravity         , ep1       = ep_1         , &
             ep2             = ep_2            , svp1      = svp1         , &
             svp2            = svp2            , svp3      = svp3         , &
             svpt0           = svpt0           , stepcu    = n_cu         , &
             cu_act_flag     = cu_act_flag     , warm_rain = warm_rain    , &
             cutop           = cutop_p         , cubot     = cubot_p      , &
             qv              = qv_p            , f_qv      = f_qv         , &
             f_qc            = f_qc            , f_qr      = f_qr         , &
             f_qi            = f_qi            , f_qs      = f_qs         , &
             rinc            = rinc_deepc_p    , rincs     = rinc_shalc_p , &
             rthcuten        = rthcuten_p      , rqvcuten  = rqvcuten_p   , &
             rqccuten        = rqccuten_p      , rqrcuten  = rqrcuten_p   , &
             rqicuten        = rqicuten_p      , rqscuten  = rqscuten_p   , &
             !jwfeng: 2021-06-03
             !Passing vars to iAMAS pool for iAMAS-chem
             umf3d           = umf_p           , uer3d     = uer_p        , &
             udr3d           = udr_p           , dmf3d     = dmf_p        , &
             der3d           = der_p           , ddr3d     = ddr_p        , &
              wu3d           =  wu_p           , ldb2d     = ldb_p        , &
             let2d           = let_p           , ksl2d     = ksl_p        , &
             lfs2d           = lfs_p           , klc2d     = klc_p        , &
             ishall2d        = ishall_p        , check2d   = check_kfeta_p, &
             ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde         , &
             ims=ims, ime=ime, jms=jms, jme=jme, kms=kds, kme=kme         , &
             its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte           &
                        )
       call iamas_timer_stop('Kain-Fritsch')

    case ("cu_kain_fritsch_sca")
       if(itimestep == 1) then
          ktau = itimestep
       else
          ktau = itimestep + 1
       endif

       call mskf_cps( &
             pcps            = pres_hyd_p      , t               = t_p             , &
             dt              = dt_dyn          , ktau            = ktau            , &
             dxCell          = dx_p            , areaCell        = area_p          , &
             cudt            = cudt            , adapt_step_flag = adapt_step_flag , &
             rho             = rho_p           , raincv          = raincv_p        , &
             pratec          = pratec_p        , nca             = nca_p           , &
             u               = u_p             , v               = v_p             , &
             th              = th_p            , pi              = pi_p            , &
             w               = w_p             , dz8w            = dz_p            , &
             w0avg           = w0avg_p         , xlv0            = xlv0            , &
             xlv1            = xlv1            , xls0            = xls0            , &
             xls1            = xls1            , cp              = cp              , &
             r               = r_d             , g               = gravity         , &
             ep1             = ep_1            , ep2             = ep_2            , &
             svp1            = svp1            , svp2            = svp2            , &
             svp3            = svp3            , svpt0           = svpt0           , &
             stepcu          = n_cu            , cu_act_flag     = cu_act_flag     , &
             warm_rain       = warm_rain       , cutop           = cutop_p         , &
             cubot           = cubot_p         , qv              = qv_p            , &
             f_qv            = f_qv            , f_qc            = f_qc            , &
             f_qr            = f_qr            , f_qi            = f_qi            , &
             f_qs            = f_qs            , kf_edrates      = kf_edrates      , &
             w_up            = w_deepc_p       , timec_kf        = timec_p         , &
             zol             = zol_p           , wstar           = wstar_p         , &
             ust             = ust_p           , pblh            = hpbl_p          , &
             cldfra_dp_kf    = cldfrac_deepc_p , cldfra_sh_kf    = cldfrac_shalc_p , &
             qc_kf           = qccu_p          , qi_kf           = qicu_p          , &
             qccu_shalc      = qccu_shalc_p    , qicu_shalc      = qicu_shalc_p    , &
             rinc            = rinc_deepc_p    , rincs           = rinc_shalc_p    , &
             xmb_deepc       = xmb_deepc_p     , xmb_shalc       = xmb_shalc_p     , &
             udr_kf          = udr_p           , ddr_kf          = ddr_p           , &
             uer_kf          = uer_p           , der_kf          = der_p           , &
             rthcuten        = rthcuten_p      , rqvcuten        = rqvcuten_p      , &
             rqccuten        = rqccuten_p      , rqrcuten        = rqrcuten_p      , &
             rqicuten        = rqicuten_p      , rqscuten        = rqscuten_p      , &
             kbcon_deepc     = kbcon_deepc_p   , ktop_deepc      = ktop_deepc_p    , &
             kbcon_shalc     = kbcon_shalc_p   , ktop_shalc      = ktop_shalc_p    , &
             rthcutens       = rthcuten_shalc_p, rqvcutens       = rqvcuten_shalc_p, &
             rqccutens       = rqccuten_shalc_p, rqicutens       = rqicuten_shalc_p, &
             rqrcutens       = rqrcuten_shalc_p, rqscutens       = rqscuten_shalc_p, &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                    )

    case("cu_tiedtke")
       call iamas_timer_start('Tiedtke')
       call cu_tiedtke( &
             pcps        = pres_hyd_p    , p8w       = pres2_hyd_p   ,               &
             znu         = znu_hyd_p     , t3d       = t_p           ,               &
             dt          = dt_dyn        , itimestep = initflag      ,               &
             stepcu      = n_cu          , raincv    = raincv_p      ,               &
             pratec      = pratec_p      , qfx       = qfx_p         ,               &
             u3d         = u_p           , v3d       = v_p           ,               &
             w           = w_p           , qv3d      = qv_p          ,               &
             qc3d        = qc_p          , qi3d      = qi_p          ,               &
             pi3d        = pi_p          , rho3d     = rho_p         ,               &
             qvften      = rqvdynten_p   , qvpblten  = rqvdynblten_p ,               &
             dz8w        = dz_p          , xland     = xland_p       ,               &
             cu_act_flag = cu_act_flag   , f_qv      = f_qv          ,               &
             f_qc        = f_qc          , f_qr      = f_qr          ,               &
             f_qi        = f_qi          , f_qs      = f_qs          ,               &
             rthcuten    = rthcuten_p    , rqvcuten  = rqvcuten_p    ,               &
             rqccuten    = rqccuten_p    , rqicuten  = rqicuten_p    ,               &
             rucuten     = rucuten_p     , rvcuten   = rvcuten_p     ,               &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                      )
       call iamas_timer_stop('Tiedtke')

    case("cu_ntiedtke")
       call iamas_timer_start('New_Tiedtke')
       call cu_ntiedtke( &
             pcps      = pres_hyd_p  , p8w         = pres2_hyd_p ,                   &
             t3d       = t_p         , dz8w        = dz_p        ,                   &
             dt        = dt_dyn      , itimestep   = initflag    ,                   &
             stepcu    = n_cu        , raincv      = raincv_p    ,                   & 
             pratec    = pratec_p    , qfx         = qfx_p       ,                   &
             hfx       = hfx_p       , xland       = xland_p     ,                   &
             dx        = dx_p        , u3d         = u_p         ,                   &
             v3d       = v_p         , w           = w_p         ,                   &
             qv3d      = qv_p        , qc3d        = qc_p        ,                   &
             qi3d      = qi_p        , pi3d        = pi_p        ,                   &
             rho3d     = rho_p       , qvften      = rqvften_p   ,                   &
             thften    = rthften_p   , cu_act_flag = cu_act_flag ,                   &
             f_qv      = f_qv        , f_qc        = f_qc        ,                   &
             f_qr      = f_qr        , f_qi        = f_qi        ,                   &
             f_qs      = f_qs        , rthcuten    = rthcuten_p  ,                   &
             rqvcuten  = rqvcuten_p  , rqccuten    = rqccuten_p  ,                   &
             rqicuten  = rqicuten_p  , rucuten     = rucuten_p   ,                   &
             rvcuten   = rvcuten_p   ,                                               &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                              )
       call iamas_timer_stop('New_Tiedtke')

    case default

 end select convection_select

!copy instantaneous and accumulated precipitation, convective tendencies, and "other" arrays
!specific to convection parameterization back to the dynamics grid:
 call convection_to_iAMAS(configs,diag_physics,tend_physics,its,ite)

!call iamas_log_write('--- end subroutine driver_convection.')
!call iamas_log_write('')

 end subroutine driver_convection

!=================================================================================================================
 subroutine convection_from_iAMAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: sfc_input
 type(iamas_pool_type),intent(in):: diag_physics
 type(iamas_pool_type),intent(in):: tend_physics

 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in):: dt_dyn

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 integer,dimension(:),pointer:: kpbl,k22_shalc,k22_deepc
 integer,dimension(:),pointer:: kbcon_deepc,kbcon_shalc,ktop_deepc,ktop_shalc
 real(kind=RKIND),dimension(:),pointer  :: areaCell,meshDensity
 real(kind=RKIND),dimension(:),pointer  :: nca,timec,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:),pointer  :: gsw,hfx,qfx,xland
 real(kind=RKIND),dimension(:),pointer  :: hpbl,ust,wstar,zol
 real(kind=RKIND),dimension(:),pointer  :: xmb_deepc,xmb_shalc
 real(kind=RKIND),dimension(:),pointer  :: rinc_deepc,rinc_shalc
 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 real(kind=RKIND),dimension(:,:),pointer:: qc_cu,qi_cu,qc_cu_shalc,qi_cu_shalc
 real(kind=RKIND),dimension(:,:),pointer:: w0avg,w_deepc
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_deepc,cldfrac_shalc
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten 
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten_shalc,rqvcuten_shalc,rqccuten_shalc,rqicuten_shalc
 real(kind=RKIND),dimension(:,:),pointer:: rqrcuten_shalc,rqscuten_shalc
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rthdynten,rthratenlw,rthratensw
 real(kind=RKIND),dimension(:,:),pointer:: rqvblten,rqvdynten,rucuten,rvcuten

!local variables:
 integer:: i,j,k
 integer:: iEdge
 real(kind=RKIND),pointer:: len_disp

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call iamas_pool_get_array(diag_physics,'cuprec',cuprec)
 call iamas_pool_get_array(diag_physics,'raincv',raincv)

 call iamas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call iamas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call iamas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call iamas_pool_get_array(tend_physics,'rqicuten',rqicuten)
 call iamas_pool_get_array(tend_physics,'rthcuten_shalc',rthcuten_shalc)
 call iamas_pool_get_array(tend_physics,'rqvcuten_shalc',rqvcuten_shalc)
 call iamas_pool_get_array(tend_physics,'rqccuten_shalc',rqccuten_shalc)
 call iamas_pool_get_array(tend_physics,'rqicuten_shalc',rqicuten_shalc)

 do j = jts,jte
 do i = its,ite
    raincv_p(i,j)     = raincv(i)
    pratec_p(i,j)     = cuprec(i)
    do k = kts,kte
       rthcuten_p(i,k,j) = rthcuten(k,i)
       rqvcuten_p(i,k,j) = rqvcuten(k,i)
       rqccuten_p(i,k,j) = rqccuten(k,i)
       rqicuten_p(i,k,j) = rqicuten(k,i)
       rthcuten_shalc_p(i,k,j) = rthcuten_shalc(k,i)
       rqvcuten_shalc_p(i,k,j) = rqvcuten_shalc(k,i)
       rqccuten_shalc_p(i,k,j) = rqccuten_shalc(k,i)
       rqicuten_shalc_p(i,k,j) = rqicuten_shalc(k,i)
    enddo
 enddo
 enddo
 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas", "cu_kain_fritsch", "cu_kain_fritsch_sca")
       call iamas_pool_get_array(diag_physics,'cubot'        ,cubot        )
       call iamas_pool_get_array(diag_physics,'cutop'        ,cutop        )

       call iamas_pool_get_array(diag_physics,'qc_cu'        ,qc_cu        )
       call iamas_pool_get_array(diag_physics,'qi_cu'        ,qi_cu        )
       call iamas_pool_get_array(diag_physics,'qc_cu_shalc'  ,qc_cu_shalc  )
       call iamas_pool_get_array(diag_physics,'qi_cu_shalc'  ,qi_cu_shalc  )
       call iamas_pool_get_array(diag_physics,'rinc_deepc'   ,rinc_deepc   )
       call iamas_pool_get_array(diag_physics,'rinc_shalc'   ,rinc_shalc   )

       call iamas_pool_get_array(diag_physics,'kbcon_deepc'  ,kbcon_deepc  )
       call iamas_pool_get_array(diag_physics,'kbcon_shalc'  ,kbcon_shalc  )
       call iamas_pool_get_array(diag_physics,'ktop_deepc'   ,ktop_deepc   )
       call iamas_pool_get_array(diag_physics,'ktop_shalc'   ,ktop_shalc   )
       call iamas_pool_get_array(diag_physics,'xmb_deepc'    ,xmb_deepc    )
       call iamas_pool_get_array(diag_physics,'xmb_shalc'    ,xmb_shalc    )

       do j = jts,jte
       do i = its,ite
          cubot_p(i,j)       = cubot(i)
          cutop_p(i,j)       = cutop(i)
          rinc_deepc_p(i,j)  = rinc_deepc(i)
          rinc_shalc_p(i,j)  = rinc_shalc(i)
          kbcon_deepc_p(i,j) = kbcon_deepc(i)
          kbcon_shalc_p(i,j) = kbcon_shalc(i)
          ktop_deepc_p(i,j)  = ktop_deepc(i)
          ktop_shalc_p(i,j)  = ktop_shalc(i)
          xmb_deepc_p(i,j)   = xmb_deepc(i)
          xmb_shalc_p(i,j)   = xmb_shalc(i)
          do k = kts,kte
             qccu_p(i,k,j)       = qc_cu(k,i)
             qicu_p(i,k,j)       = qi_cu(k,i)
             qccu_shalc_p(i,k,j) = qc_cu_shalc(k,i)
             qicu_shalc_p(i,k,j) = qi_cu_shalc(k,i)
          enddo
       enddo
       enddo

       cu_kfgf_select: select case(convection_scheme)

          case ("cu_grell_freitas")
             call iamas_pool_get_config(configs,'config_len_disp',len_disp)

             call iamas_pool_get_array(mesh,'areaCell'   ,areaCell   )
             call iamas_pool_get_array(mesh,'meshDensity',meshDensity)
             call iamas_pool_get_array(mesh,'zgrid'      ,zgrid      )

             call iamas_pool_get_array(sfc_input,'xland',xland)

             call iamas_pool_get_array(diag_physics,'gsw'          ,gsw          )
             call iamas_pool_get_array(diag_physics,'hfx'          ,hfx          )
             call iamas_pool_get_array(diag_physics,'qfx'          ,qfx          )
             call iamas_pool_get_array(diag_physics,'kpbl'         ,kpbl         )
             call iamas_pool_get_array(diag_physics,'k22_deepc'    ,k22_deepc    )
             call iamas_pool_get_array(diag_physics,'k22_shalc'    ,k22_shalc    )

             call iamas_pool_get_array(tend_physics,'rqvblten'  ,rqvblten  )
             call iamas_pool_get_array(tend_physics,'rqvdynten' ,rqvdynten )
             call iamas_pool_get_array(tend_physics,'rthblten'  ,rthblten  )
             call iamas_pool_get_array(tend_physics,'rthdynten' ,rthdynten )
             call iamas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
             call iamas_pool_get_array(tend_physics,'rthratensw',rthratensw)
             call iamas_pool_get_array(tend_physics,'rucuten'   ,rucuten   )
             call iamas_pool_get_array(tend_physics,'rvcuten'   ,rvcuten   )

             do j = jts,jte
             do i = its,ite
                dx_p(i,j) = len_disp / meshDensity(i)**0.25
                area_p(i,j)          = areaCell(i)
                ht_p(i,j)            = zgrid(1,i)
                xland_p(i,j)         = xland(i)
                gsw_p(i,j)           = gsw(i)
                hfx_p(i,j)           = hfx(i)
                qfx_p(i,j)           = qfx(i)
                kpbl_p(i,j)          = kpbl(i)

                k22_deepc_p(i,j)     = k22_deepc(i)
                k22_shalc_p(i,j)     = k22_shalc(i)

                do k = kts,kte
                   rthblten_p(i,k,j)  = rthblten(k,i)
                   rthdynten_p(i,k,j) = rthdynten(k,i)
                   rthraten_p(i,k,j)  = rthratenlw(k,i) + rthratensw(k,i)
                   rqvblten_p(i,k,j)  = rqvblten(k,i)
                   rqvdynten_p(i,k,j) = rqvdynten(k,i)
                   rucuten_p(i,k,j)   = rucuten(k,i)
                   rvcuten_p(i,k,j)   = rvcuten(k,i)
                enddo
             enddo
             enddo

          case ("cu_kain_fritsch", "cu_kain_fritsch_sca")
             call iamas_pool_get_config(configs,'config_len_disp',len_disp)

             call iamas_pool_get_array(mesh,'areaCell',areaCell)
             call iamas_pool_get_array(mesh,'meshDensity',meshDensity)

             call iamas_pool_get_array(diag_physics,'nca'     ,nca     )
             call iamas_pool_get_array(diag_physics,'w0avg'   ,w0avg   )
             call iamas_pool_get_array(diag_physics,'timec'   ,timec   )
             call iamas_pool_get_array(diag_physics,'cldfrac_deepc',cldfrac_deepc)
             call iamas_pool_get_array(diag_physics,'cldfrac_shalc',cldfrac_shalc)

             call iamas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
             call iamas_pool_get_array(tend_physics,'rqscuten',rqscuten)
             call iamas_pool_get_array(tend_physics,'rqrcuten_shalc',rqrcuten_shalc)
             call iamas_pool_get_array(tend_physics,'rqscuten_shalc',rqscuten_shalc)

             do j = jts,jte
             do i = its,ite
                !area of grid-cell:       
                dx_p(i,j)    = len_disp / meshDensity(i)**0.25
                area_p(i,j)  = areaCell(i)
                timec_p(i,j) = timec(i)

                do k = kts,kte
                   rqrcuten_p(i,k,j) = rqrcuten(k,i)
                   rqscuten_p(i,k,j) = rqscuten(k,i)
                   rqrcuten_shalc_p(i,k,j) = rqrcuten_shalc(k,i)
                   rqscuten_shalc_p(i,k,j) = rqscuten_shalc(k,i)
                enddo

                !decreases the characteristic time period that convection remains active. When nca_p
                !becomes less than the convective timestep, convective tendencies and precipitation
                !are reset to zero (note that this is also done in subroutine kf_eta_cps).
                nca_p(i,j) = nca(i)
                !jwfeng:2021-06-03
                ldb_p(i,j) = 0
                let_p(i,j) = 0
                ksl_p(i,j) = 0
                lfs_p(i,j) = 0
                klc_p(i,j) = 0
                ishall_p(i,j) = 0
                check_kfeta_p(i,j) = 0

                do k = kts,kte
                   w0avg_p(i,k,j) = w0avg(k,i)
                   cldfrac_deepc_p(i,k,j) = cldfrac_deepc(k,i)
                   cldfrac_shalc_p(i,k,j) = cldfrac_shalc(k,i)
                   !jwfeng:2021-06-03
                   umf_p(i,k,j) = 0.0
                   uer_p(i,k,j) = 0.0
                   udr_p(i,k,j) = 0.0
                   dmf_p(i,k,j) = 0.0
                   der_p(i,k,j) = 0.0
                   ddr_p(i,k,j) = 0.0
                    wu_p(i,k,j) = 0.0
                enddo
             enddo
             enddo
             cu_kf_select: select case(convection_scheme)
                case ("cu_kain_fritsch_sca")
                   call iamas_pool_get_array(diag_physics,'hpbl'   , hpbl   )
                   call iamas_pool_get_array(diag_physics,'ust'    , ust    )
                   call iamas_pool_get_array(diag_physics,'wstar'  , wstar  )
                   call iamas_pool_get_array(diag_physics,'zol'    , zol    )
                   call iamas_pool_get_array(diag_physics,'w_deepc', w_deepc)
                   do j = jts,jte
                   do i = its,ite
                      hpbl_p(i,j)  = hpbl(i)
                      ust_p(i,j)   = ust(i)
                      wstar_p(i,j) = wstar(i)
                      zol_p(i,j)   = zol(i)
                      do k = kts,kte
                         w_deepc_p(i,k,j)   = w_deepc(k,i)
                      enddo
                   enddo
                   enddo

                case default

             end select cu_kf_select

          case default

       end select cu_kfgf_select

    case ("cu_tiedtke","cu_ntiedtke")
       call iamas_pool_get_array(sfc_input,'xland',xland)
       call iamas_pool_get_array(diag_physics,'hfx',hfx)
       call iamas_pool_get_array(diag_physics,'qfx',qfx)

       call iamas_pool_get_array(tend_physics,'rqvblten' ,rqvblten )
       call iamas_pool_get_array(tend_physics,'rqvdynten',rqvdynten)
       call iamas_pool_get_array(tend_physics,'rucuten'  ,rucuten  )
       call iamas_pool_get_array(tend_physics,'rvcuten'  ,rvcuten  )

       do j = jts,jte
       do i = its,ite
          xland_p(i,j) = xland(i)
          hfx_p(i,j)   = hfx(i)
          qfx_p(i,j)   = qfx(i)
       enddo

       do k = kts,kte
       do i = its,ite
          rucuten_p(i,k,j) = rucuten(k,i)
          rvcuten_p(i,k,j) = rvcuten(k,i)
       enddo
       enddo
       enddo

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             do j = jts,jte
             do k = kts,kte
             do i = its,ite
                rqvdynblten_p(i,k,j) = rqvblten(k,i) 
                rqvdynten_p(i,k,j)   = rqvdynten(k,i)
             enddo
             enddo
             enddo
 
          case("cu_ntiedtke")
             call iamas_pool_get_config(configs,'config_len_disp',len_disp)
             call iamas_pool_get_array(mesh,'meshDensity',meshDensity)

             call iamas_pool_get_array(tend_physics,'rqvdynten' ,rqvdynten )
             call iamas_pool_get_array(tend_physics,'rqvblten'  ,rqvblten  )
             call iamas_pool_get_array(tend_physics,'rthdynten' ,rthdynten )
             call iamas_pool_get_array(tend_physics,'rthblten'  ,rthblten  )
             call iamas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
             call iamas_pool_get_array(tend_physics,'rthratensw',rthratensw)
             
             do j = jts,jte
             do i = its,ite
                dx_p(i,j) = len_disp / meshDensity(i)**0.25
             enddo

             do k = kts,kte
             do i = its,ite
                rqvften_p(i,k,j) = rqvdynten(k,i) + rqvblten(k,i)
                rthften_p(i,k,j) = (rthdynten(k,i) + rthblten(k,i) &
                                 + rthratenlw(k,i) + rthratensw(k,i)) * pi_p(i,k,j) 
             enddo
             enddo
             enddo

          case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine convection_from_iAMAS

!=================================================================================================================
 subroutine convection_to_iAMAS(configs,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: tend_physics

!local variables:
 integer:: i,k,j

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 integer,dimension(:),pointer:: k22_shalc,k22_deepc
 integer,dimension(:),pointer:: kbcon_deepc,kbcon_shalc,ktop_deepc,ktop_shalc
 real(kind=RKIND),dimension(:),pointer  :: nca,timec,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:),pointer  :: xmb_deepc,xmb_shalc
 real(kind=RKIND),dimension(:),pointer  :: rinc_deepc,rinc_shalc,rinc_overlap
 real(kind=RKIND),dimension(:,:),pointer:: qc_cu,qi_cu,qi_cu_shalc,qc_cu_shalc
 real(kind=RKIND),dimension(:,:),pointer:: w0avg,w_deepc
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac_deepc,cldfrac_shalc
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten_shalc,rqvcuten_shalc,rqccuten_shalc,rqicuten_shalc
 real(kind=RKIND),dimension(:,:),pointer:: rqrcuten_shalc,rqscuten_shalc
 real(kind=RKIND),dimension(:,:),pointer:: rucuten,rvcuten
 !jwfeng: 2021-06-03
 real(kind=RKIND),dimension(:,:),pointer:: umf,uer,udr,dmf,der,ddr,wu
 integer,dimension(:),pointer:: ldb,let,ksl,lfs,klc,ishall,check_kfeta

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call iamas_pool_get_array(diag_physics,'cuprec',cuprec)
 call iamas_pool_get_array(diag_physics,'raincv',raincv)

 call iamas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call iamas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call iamas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call iamas_pool_get_array(tend_physics,'rqicuten',rqicuten)
 call iamas_pool_get_array(tend_physics,'rthcuten_shalc',rthcuten_shalc)
 call iamas_pool_get_array(tend_physics,'rqvcuten_shalc',rqvcuten_shalc)
 call iamas_pool_get_array(tend_physics,'rqccuten_shalc',rqccuten_shalc)
 call iamas_pool_get_array(tend_physics,'rqicuten_shalc',rqicuten_shalc)

 do j = jts,jte
 do i = its,ite
    raincv(i)     = raincv_p(i,j)
    cuprec(i)     = pratec_p(i,j)
    do k = kts, kte
       rthcuten(k,i)    = rthcuten_p(i,k,j)
       rqvcuten(k,i)    = rqvcuten_p(i,k,j)
       rqccuten(k,i)    = rqccuten_p(i,k,j)
       rqicuten(k,i)    = rqicuten_p(i,k,j)
       rthcuten_shalc(k,i)    = rthcuten_shalc_p(i,k,j)
       rqvcuten_shalc(k,i)    = rqvcuten_shalc_p(i,k,j)
       rqccuten_shalc(k,i)    = rqccuten_shalc_p(i,k,j)
       rqicuten_shalc(k,i)    = rqicuten_shalc_p(i,k,j)
    enddo
 enddo
 enddo

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas", "cu_kain_fritsch", "cu_kain_fritsch_sca")
       call iamas_pool_get_array(diag_physics,'cubot'        ,cubot        )
       call iamas_pool_get_array(diag_physics,'cutop'        ,cutop        )
       call iamas_pool_get_array(diag_physics,'qc_cu'        ,qc_cu        )
       call iamas_pool_get_array(diag_physics,'qi_cu'        ,qi_cu        )
       call iamas_pool_get_array(diag_physics,'qc_cu_shalc'  ,qc_cu_shalc  )
       call iamas_pool_get_array(diag_physics,'qi_cu_shalc'  ,qi_cu_shalc  )
       call iamas_pool_get_array(diag_physics,'rinc_deepc'   ,rinc_deepc   )
       call iamas_pool_get_array(diag_physics,'rinc_shalc'   ,rinc_shalc   )
       call iamas_pool_get_array(diag_physics,'rinc_overlap' ,rinc_overlap )
       call iamas_pool_get_array(diag_physics,'kbcon_deepc'  ,kbcon_deepc  )
       call iamas_pool_get_array(diag_physics,'kbcon_shalc'  ,kbcon_shalc  )
       call iamas_pool_get_array(diag_physics,'ktop_deepc'   ,ktop_deepc   )
       call iamas_pool_get_array(diag_physics,'ktop_shalc'   ,ktop_shalc   )
       call iamas_pool_get_array(diag_physics,'xmb_deepc'    ,xmb_deepc    )
       call iamas_pool_get_array(diag_physics,'xmb_shalc'    ,xmb_shalc    )

       do j = jts,jte
       do i = its,ite
          cubot(i)       = cubot_p(i,j)
          cutop(i)       = cutop_p(i,j)
          kbcon_deepc(i) = kbcon_deepc_p(i,j)
          kbcon_shalc(i) = kbcon_shalc_p(i,j)
          ktop_deepc(i)  = ktop_deepc_p(i,j)
          ktop_shalc(i)  = ktop_shalc_p(i,j)
          xmb_deepc(i)   = xmb_deepc_p(i,j)
          xmb_shalc(i)   = xmb_shalc_p(i,j)
          rinc_deepc(i)  = rinc_deepc_p(i,j)
          rinc_shalc(i)  = rinc_shalc_p(i,j)
          if (rinc_deepc(i) .gt. 0._RKIND .and. rinc_shalc(i) .gt. 0._RKIND) then
             rinc_overlap(i) = 1._RKIND
          else
             rinc_overlap(i) = 0._RKIND
          end if

          do k = kts, kte
             qc_cu(k,i)       = qccu_p(i,k,j)
             qi_cu(k,i)       = qicu_p(i,k,j)
             qc_cu_shalc(k,i) = qccu_shalc_p(i,k,j)
             qi_cu_shalc(k,i) = qicu_shalc_p(i,k,j)
          enddo
       enddo
       enddo

       cu_kfgf_select: select case(convection_scheme)

          case ("cu_grell_freitas")
             call iamas_pool_get_array(diag_physics,'k22_deepc'    ,k22_deepc    )
             call iamas_pool_get_array(diag_physics,'k22_shalc'    ,k22_shalc    )

             call iamas_pool_get_array(tend_physics,'rucuten',rucuten)
             call iamas_pool_get_array(tend_physics,'rvcuten',rvcuten)

             do j = jts,jte
             do i = its,ite
                k22_deepc(i)     = k22_deepc_p(i,j)
                k22_shalc(i)     = k22_shalc_p(i,j)

                do k = kts,kte
                   rucuten(k,i)  = rucuten_p(i,k,j)
                   rvcuten(k,i)  = rvcuten_p(i,k,j)
                enddo
             enddo
             enddo

          case ("cu_kain_fritsch", "cu_kain_fritsch_sca")
             call iamas_pool_get_array(diag_physics,'nca'  ,nca  )
             call iamas_pool_get_array(diag_physics,'w0avg',w0avg)
             call iamas_pool_get_array(diag_physics,'timec',timec)
             call iamas_pool_get_array(diag_physics,'cldfrac_deepc',cldfrac_deepc)
             call iamas_pool_get_array(diag_physics,'cldfrac_shalc',cldfrac_shalc)
             !jwfeng: 2021-06-03
             call iamas_pool_get_array(diag_physics,'umf'  ,umf  )
             call iamas_pool_get_array(diag_physics,'uer'  ,uer  )
             call iamas_pool_get_array(diag_physics,'udr'  ,udr  )
             call iamas_pool_get_array(diag_physics,'dmf'  ,dmf  )
             call iamas_pool_get_array(diag_physics,'der'  ,der  )
             call iamas_pool_get_array(diag_physics,'ddr'  ,ddr  )
             call iamas_pool_get_array(diag_physics, 'wu'  , wu  )
             call iamas_pool_get_array(diag_physics,'ldb'  ,ldb  )
             call iamas_pool_get_array(diag_physics,'let'  ,let  )
             call iamas_pool_get_array(diag_physics,'ksl'  ,ksl  )
             call iamas_pool_get_array(diag_physics,'lfs'  ,lfs  )
             call iamas_pool_get_array(diag_physics,'klc'  ,klc  )
             call iamas_pool_get_array(diag_physics,'ishall', ishall )
             call iamas_pool_get_array(diag_physics,'check_kfeta', check_kfeta)

             call iamas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
             call iamas_pool_get_array(tend_physics,'rqscuten',rqscuten)
             call iamas_pool_get_array(tend_physics,'rqrcuten_shalc',rqrcuten_shalc)
             call iamas_pool_get_array(tend_physics,'rqscuten_shalc',rqscuten_shalc)

             do j = jts,jte
             do i = its,ite
                nca(i)   = nca_p(i,j)
                timec(i) = timec_p(i,j)
                !jwfeng: 2021-06-03
                ldb(i)   = ldb_p(i,j)
                let(i)   = let_p(i,j)
                ksl(i)   = ksl_p(i,j)
                lfs(i)   = lfs_p(i,j)
                klc(i)   = klc_p(i,j)
                ishall(i) = ishall_p(i,j)
                check_kfeta(i) = check_kfeta_p(i,j)
                do k = kts, kte
                   w0avg(k,i)          = w0avg_p(i,k,j)
                   cldfrac_deepc(k,i)  = cldfrac_deepc_p(i,k,j)
                   cldfrac_shalc(k,i)  = cldfrac_shalc_p(i,k,j)
                   rqrcuten(k,i)       = rqrcuten_p(i,k,j)
                   rqscuten(k,i)       = rqscuten_p(i,k,j)
                   rqrcuten_shalc(k,i) = rqrcuten_shalc_p(i,k,j)
                   rqscuten_shalc(k,i) = rqscuten_shalc_p(i,k,j)
                   !jwfeng: 2021-06-03
                   umf(k,i)      = umf_p(i,k,j)
                   uer(k,i)      = uer_p(i,k,j)
                   udr(k,i)      = udr_p(i,k,j)
                   dmf(k,i)      = dmf_p(i,k,j)
                   der(k,i)      = der_p(i,k,j)
                   ddr(k,i)      = ddr_p(i,k,j)
                    wu(k,i)      =  wu_p(i,k,j)
                enddo
             enddo
             enddo

             cu_kf_select: select case(convection_scheme)
                case ("cu_kain_fritsch_sca")

                   call iamas_pool_get_array(diag_physics,'w_deepc', w_deepc)

                   do j = jts,jte
                   do i = its,ite
                      do k = kts, kte
                         w_deepc(k,i) = w_deepc_p(i,k,j)
                      enddo
                   enddo
                   enddo

                case default
             end select cu_kf_select

          case default
       end select cu_kfgf_select

    case ("cu_tiedtke","cu_ntiedtke")
       call iamas_pool_get_array(tend_physics,'rucuten',rucuten)
       call iamas_pool_get_array(tend_physics,'rvcuten',rvcuten)

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          rucuten(k,i) = rucuten_p(i,k,j)
          rvcuten(k,i) = rvcuten_p(i,k,j)
       enddo
       enddo
       enddo

    case default

 end select convection_select

 end subroutine convection_to_iAMAS

!=================================================================================================================
 subroutine update_convection_step1(configs,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: tend_physics

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten_shalc,rqvcuten_shalc,rqccuten_shalc,rqicuten_shalc
 real(kind=RKIND),dimension(:,:),pointer:: rqrcuten_shalc,rqscuten_shalc

!local variables and arrays:
 integer:: i,k

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 convection_select: select case(convection_scheme)

    case ("cu_kain_fritsch", "cu_kain_fritsch_sca")
      call iamas_pool_get_array(diag_physics,'nca'   ,nca   )
      call iamas_pool_get_array(diag_physics,'cubot' ,cubot )
      call iamas_pool_get_array(diag_physics,'cutop' ,cutop )
      call iamas_pool_get_array(diag_physics,'cuprec',cuprec)
      call iamas_pool_get_array(diag_physics,'raincv',raincv)

      call iamas_pool_get_array(tend_physics,'rthcuten',rthcuten)
      call iamas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
      call iamas_pool_get_array(tend_physics,'rqccuten',rqccuten)
      call iamas_pool_get_array(tend_physics,'rqicuten',rqicuten)
      call iamas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
      call iamas_pool_get_array(tend_physics,'rqscuten',rqscuten)
      call iamas_pool_get_array(tend_physics,'rthcuten_shalc',rthcuten_shalc)
      call iamas_pool_get_array(tend_physics,'rqvcuten_shalc',rqvcuten_shalc)
      call iamas_pool_get_array(tend_physics,'rqccuten_shalc',rqccuten_shalc)
      call iamas_pool_get_array(tend_physics,'rqicuten_shalc',rqicuten_shalc)
      call iamas_pool_get_array(tend_physics,'rqrcuten_shalc',rqrcuten_shalc)
      call iamas_pool_get_array(tend_physics,'rqscuten_shalc',rqscuten_shalc)

       do i = its, ite
          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          if(nca(i) .gt. 0.) then
             nca(i) = nca(i) - dt_dyn
             
             if(nca(i) .lt. 0.5*dt_dyn) then
                do k = kts,kte
                   rthcuten(k,i) = 0._RKIND
                   rqvcuten(k,i) = 0._RKIND
                   rqccuten(k,i) = 0._RKIND
                   rqrcuten(k,i) = 0._RKIND
                   rqicuten(k,i) = 0._RKIND
                   rqscuten(k,i) = 0._RKIND
                   rthcuten_shalc(k,i) = 0._RKIND
                   rqvcuten_shalc(k,i) = 0._RKIND
                   rqccuten_shalc(k,i) = 0._RKIND
                   rqrcuten_shalc(k,i) = 0._RKIND
                   rqicuten_shalc(k,i) = 0._RKIND
                   rqscuten_shalc(k,i) = 0._RKIND
                enddo
                raincv(i) = 0._RKIND
                cuprec(i) = 0._RKIND
                cubot(i)  = kte+1
                cutop(i)  = kts
             endif
          endif
       enddo

    case default

 end select convection_select

 end subroutine update_convection_step1

!=================================================================================================================
 subroutine update_convection_step2(clock,configs,diag_physics,its,ite)
!=================================================================================================================

 use iamas_timekeeping

!input arguments:
 type(iAMAS_Clock_type), intent(in):: clock
 type(iamas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(iamas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer,dimension(:),pointer:: i_rainc

 real(kind=RKIND),pointer:: bucket_rainc
 real(kind=RKIND),dimension(:),pointer:: cuprec,rainc,rainc_hr

!local variables and arrays:
 integer:: i

 !czhao obtain time information for hourly integration
 type(iAMAS_Time_Type):: currTime
 character(len=StrKIND)   :: timeStamp
 integer:: ierr
 real(kind=RKIND):: utc_h
 integer:: yr,h,m,s,s_n,s_d,DoY

!----------------------------------------------------------------------------
  !czhao obtain time for hourly integration 
!  call iamas_log_write('czhao check time clock at the begin')

   currTime = iamas_get_clock_time(clock,iAMAS_NOW,ierr)
   call iamas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,H=h,M=m, &
                    S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
   utc_h = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0

! call iamas_log_write('czhao check h  = $i', intArgs=(/h/))
! call iamas_log_write('czhao check m  = $i', intArgs=(/m/))
! call iamas_log_write('czhao check utc_h=$r', realArgs=(/utc_h/))

 call iamas_pool_get_config(configs,'config_bucket_rainc',bucket_rainc)

 call iamas_pool_get_array(diag_physics,'i_rainc',i_rainc)
 call iamas_pool_get_array(diag_physics,'cuprec' ,cuprec )
 call iamas_pool_get_array(diag_physics,'rainc'  ,rainc  )
 call iamas_pool_get_array(diag_physics,'rainc_hr' ,rainc_hr )


!update the accumulated precipitation at the end of each dynamic time-step:
 do i = its, ite
    rainc(i) = rainc(i) + cuprec(i) * dt_dyn

    !czhao obtain hourly rainnc
    if (utc_h .eq. h) then
    rainc_hr(i) = 0.0
    endif
    rainc_hr(i) = rainc_hr(i) + cuprec(i) * dt_dyn 

    if(l_acrain .and. bucket_rainc.gt.0._RKIND .and. &
       rainc(i).gt.bucket_rainc) then
       i_rainc(i) = i_rainc(i) + 1
       rainc(i)   = rainc(i) - bucket_rainc
    endif
 enddo

 end subroutine update_convection_step2

!=================================================================================================================
 end module iamas_atmphys_driver_convection
!=================================================================================================================
