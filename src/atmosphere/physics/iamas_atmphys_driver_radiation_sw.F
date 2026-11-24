! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_driver_radiation_sw
 use iamas_kind_types
 use iamas_pool_routines
 use iamas_timer, only : iamas_timer_start, iamas_timer_stop

 use iamas_atmphys_constants
 use iamas_atmphys_manager, only: gmt,curr_julday,julday,year
 use iamas_atmphys_camrad_init
 use iamas_atmphys_rrtmg_swinit
 use iamas_atmphys_vars
 use iamas_atmphys_utilities, only:physics_message
 
!wrf physics:
 use module_ra_cam
 use module_ra_rrtmg_sw

 use module_data_binaer_asect, only:ncomp_aer, lptr_dust_aer, lptr_cl_aer, lptr_na_aer

 implicit none
 private
 public:: allocate_radiation_sw,   &
          deallocate_radiation_sw, &
          driver_radiation_sw,     &
          init_radiation_sw,       &
          radconst

!iAMAS driver for parameterization of shortwave radiation codes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in iamas_atmphys_driver_radiation_sw:
! ------------------------------------------------
! allocate_radiation_sw  : allocate local arrays for parameterization of sw radiation codes.
! deallocate_radiation_sw: deallocate local arrays for parameterization of sw radiation codes.
! init_radiation_sw      : initialization of individual sw radiation codes.
! driver_radiation_sw    : main driver (called from subroutine physics_driver).
! radiation_sw_from_iAMAS : initialize local arrays.
! radiation_sw_to_iAMAS   : copy local arrays to iAMAS arrays.
! radconst               : calculate solar declination,...
!
! WRF physics called from driver_radiation_sw:
! --------------------------------------------
! * module_ra_cam        : CAM short wave radiation code.
! * module_ra_rrtmg_sw   : RRTMG short wave radiation code.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutines
!   rrtmg_sw and camrad.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added structure diag in the call to subroutine init_radiation_sw and call to subroutine
!   camradinit for initialization of variable mxaerl.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-01.
! * modified the call to subroutine rrtmg_swrad to include the option of using the same ozone
!   climatology as the one used in the CAM radiation codes.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-17.
! * in call to subroutine rrtmg_swrad, replaced the variable g (that originally pointed to
!   gravity) with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * cleaned-up the call to rrtmg_swrad after cleaning up subroutine rrtmg_swrad in module_ra_rrtmg_sw.F.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-05.
! * added the cloud radii for cloud water, cloud ice, and snow calculated in the Thompson cloud microphysics
!   scheme in the call to subroutine rrtmg_swrad.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-07.
! * removed qr_p, and qg_p in the call to rrtmg_swrad since not used in the calculation of the cloud optical
!   properties.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-08.
! * in the call to rrtmg_swrad, substituted the variables qv_p, qc_p, qi_p, and qs_p with qvrad_p, qcrad_p,
!   qirad_p, and qsrad_p initialized in subroutine cloudiness_from_iAMAS.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-09.
! * in subroutines radiation_sw_from_iAMAS, revised the initialization of re_cloud, re_ice, re_snow, to
!   handle the case when the cloud microphysics parameterization is turned off, i.e. config_microp_scheme='off'.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-10.
! * since we removed the local variable radt_sw_scheme from iamas_atmphys_vars.F, now defines radt_sw_scheme
!   as a pointer to config_radt_sw_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2917-02-16.


 contains


!=================================================================================================================
 subroutine allocate_radiation_sw(configs,xtime_s)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs
 real(kind=RKIND),intent(in):: xtime_s

!local pointers:
 character(len=StrKIND),pointer:: radt_sw_scheme

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme)
!added by qkyang to support urban model+++++++++++++++++++++
!  if(.not.allocated(hrang_p)      ) allocate(hrang_p(ims:ime,jms:jme)              )
!added by qkyang to support urban model---------------------
 if(.not.allocated(f_ice)        ) allocate(f_ice(ims:ime,kms:kme,jms:jme)        )
 if(.not.allocated(f_rain)       ) allocate(f_rain(ims:ime,kms:kme,jms:jme)       )

 if(.not.allocated(xlat_p)       ) allocate(xlat_p(ims:ime,jms:jme)               )
 if(.not.allocated(xlon_p)       ) allocate(xlon_p(ims:ime,jms:jme)               )

 if(.not.allocated(sfc_albedo_p) ) allocate(sfc_albedo_p(ims:ime,jms:jme)         )
 if(.not.allocated(snow_p)       ) allocate(snow_p(ims:ime,jms:jme)               )
 if(.not.allocated(tsk_p)        ) allocate(tsk_p(ims:ime,jms:jme)                )
 if(.not.allocated(xice_p)       ) allocate(xice_p(ims:ime,jms:jme)               )
 if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)              )

 if(.not.allocated(coszr_p)      ) allocate(coszr_p(ims:ime,jms:jme)              )
 if(.not.allocated(gsw_p)        ) allocate(gsw_p(ims:ime,jms:jme)                )
 if(.not.allocated(swcf_p)       ) allocate(swcf_p(ims:ime,jms:jme)               )
 if(.not.allocated(swdnb_p)      ) allocate(swdnb_p(ims:ime,jms:jme)              )
 if(.not.allocated(swdnbc_p)     ) allocate(swdnbc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swdnt_p)      ) allocate(swdnt_p(ims:ime,jms:jme)              )
 if(.not.allocated(swdntc_p)     ) allocate(swdntc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swupb_p)      ) allocate(swupb_p(ims:ime,jms:jme)              )
 if(.not.allocated(swupbc_p)     ) allocate(swupbc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swupt_p)      ) allocate(swupt_p(ims:ime,jms:jme)              )
 if(.not.allocated(swuptc_p)     ) allocate(swuptc_p(ims:ime,jms:jme)             )

 if(.not.allocated(rthratensw_p) ) allocate(rthratensw_p(ims:ime,kms:kme,jms:jme) )

!jwfeng from USTC AEMOL: 2020/10/09
#ifdef DO_CHEMISTRY
 if(.not.allocated(tauaersw_p)      ) allocate(tauaersw_p(ims:ime,kms:kme,jms:jme,1:14)      )
 if(.not.allocated(gaersw_p)        ) allocate(gaersw_p(ims:ime,kms:kme,jms:jme,1:14)        )
 if(.not.allocated(waersw_p)        ) allocate(waersw_p(ims:ime,kms:kme,jms:jme,1:14)        )
 if(.not.allocated(tauaersw_dust_p) ) allocate(tauaersw_dust_p(ims:ime,kms:kme,jms:jme,1:14) )
 if(.not.allocated(gaersw_dust_p)   ) allocate(gaersw_dust_p(ims:ime,kms:kme,jms:jme,1:14)   )
 if(.not.allocated(waersw_dust_p)   ) allocate(waersw_dust_p(ims:ime,kms:kme,jms:jme,1:14)   )
 if(.not.allocated(tauaersw_seas_p) ) allocate(tauaersw_seas_p(ims:ime,kms:kme,jms:jme,1:14) )
 if(.not.allocated(gaersw_seas_p)   ) allocate(gaersw_seas_p(ims:ime,kms:kme,jms:jme,1:14)   )
 if(.not.allocated(waersw_seas_p)   ) allocate(waersw_seas_p(ims:ime,kms:kme,jms:jme,1:14)   )
 if(.not.allocated(tauaersw_other_p)) allocate(tauaersw_other_p(ims:ime,kms:kme,jms:jme,1:14))
 if(.not.allocated(gaersw_other_p)  ) allocate(gaersw_other_p(ims:ime,kms:kme,jms:jme,1:14)  )
 if(.not.allocated(waersw_other_p)  ) allocate(waersw_other_p(ims:ime,kms:kme,jms:jme,1:14)  )
 if(.not.allocated(swdnb_dust_p)    ) allocate(swdnb_dust_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swdnbc_dust_p)   ) allocate(swdnbc_dust_p(ims:ime,jms:jme)                )
 if(.not.allocated(swdnt_dust_p)    ) allocate(swdnt_dust_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swdntc_dust_p)   ) allocate(swdntc_dust_p(ims:ime,jms:jme)                )
 if(.not.allocated(swupb_dust_p)    ) allocate(swupb_dust_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swupbc_dust_p)   ) allocate(swupbc_dust_p(ims:ime,jms:jme)                )
 if(.not.allocated(swupt_dust_p)    ) allocate(swupt_dust_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swuptc_dust_p)   ) allocate(swuptc_dust_p(ims:ime,jms:jme)                )
 if(.not.allocated(swdnb_seas_p)    ) allocate(swdnb_seas_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swdnbc_seas_p)   ) allocate(swdnbc_seas_p(ims:ime,jms:jme)                )
 if(.not.allocated(swdnt_seas_p)    ) allocate(swdnt_seas_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swdntc_seas_p)   ) allocate(swdntc_seas_p(ims:ime,jms:jme)                )
 if(.not.allocated(swupb_seas_p)    ) allocate(swupb_seas_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swupbc_seas_p)   ) allocate(swupbc_seas_p(ims:ime,jms:jme)                )
 if(.not.allocated(swupt_seas_p)    ) allocate(swupt_seas_p(ims:ime,jms:jme)                 )
 if(.not.allocated(swuptc_seas_p)   ) allocate(swuptc_seas_p(ims:ime,jms:jme)                )
 if(.not.allocated(swdnb_other_p)   ) allocate(swdnb_other_p(ims:ime,jms:jme)                )
 if(.not.allocated(swdnbc_other_p)  ) allocate(swdnbc_other_p(ims:ime,jms:jme)               )
 if(.not.allocated(swdnt_other_p)   ) allocate(swdnt_other_p(ims:ime,jms:jme)                )
 if(.not.allocated(swdntc_other_p)  ) allocate(swdntc_other_p(ims:ime,jms:jme)               )
 if(.not.allocated(swupb_other_p)   ) allocate(swupb_other_p(ims:ime,jms:jme)                )
 if(.not.allocated(swupbc_other_p)  ) allocate(swupbc_other_p(ims:ime,jms:jme)               )
 if(.not.allocated(swupt_other_p)   ) allocate(swupt_other_p(ims:ime,jms:jme)                )
 if(.not.allocated(swuptc_other_p)  ) allocate(swuptc_other_p(ims:ime,jms:jme)               )

 if(.not.allocated(rthratensw_dust_p) ) allocate(rthratensw_dust_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rthratensw_seas_p) ) allocate(rthratensw_seas_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(rthratensw_other_p)) allocate(rthratensw_other_p(ims:ime,kms:kme,jms:jme))
#endif

 radiation_sw_select: select case (trim(radt_sw_scheme))

    case("rrtmg_sw")
       if(.not.allocated(recloud_p)    ) allocate(recloud_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(reice_p)      ) allocate(reice_p(ims:ime,kms:kme,jms:jme)         )
       if(.not.allocated(resnow_p)     ) allocate(resnow_p(ims:ime,kms:kme,jms:jme)        )

       if(.not.allocated(alswvisdir_p) ) allocate(alswvisdir_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswvisdif_p) ) allocate(alswvisdif_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswnirdir_p) ) allocate(alswnirdir_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswnirdif_p) ) allocate(alswnirdif_p(ims:ime,jms:jme)            )
       if(.not.allocated(swvisdir_p)   ) allocate(swvisdir_p(ims:ime,jms:jme)              )
       if(.not.allocated(swvisdif_p)   ) allocate(swvisdif_p(ims:ime,jms:jme)              )
       if(.not.allocated(swnirdir_p)   ) allocate(swnirdir_p(ims:ime,jms:jme)              )
       if(.not.allocated(swnirdif_p)   ) allocate(swnirdif_p(ims:ime,jms:jme)              )

       if(.not.allocated(swdnflx_p)    ) allocate(swdnflx_p(ims:ime,kms:kme+1,jms:jme)     )
       if(.not.allocated(swdnflxc_p)   ) allocate(swdnflxc_p(ims:ime,kms:kme+1,jms:jme)    )
       if(.not.allocated(swupflx_p)    ) allocate(swupflx_p(ims:ime,kms:kme+1,jms:jme)     )
       if(.not.allocated(swupflxc_p)   ) allocate(swupflxc_p(ims:ime,kms:kme+1,jms:jme)    )

       if(.not.allocated(pin_p)        ) allocate(pin_p(num_oznlevels)                     )
       if(.not.allocated(o3clim_p)     ) allocate(o3clim_p(ims:ime,1:num_oznlevels,jms:jme))

    case("cam_sw")
       if(.not.allocated(glw_p)        ) allocate(glw_p(ims:ime,jms:jme)                )
       if(.not.allocated(lwcf_p)       ) allocate(lwcf_p(ims:ime,jms:jme)               )
       if(.not.allocated(lwdnb_p)      ) allocate(lwdnb_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwdnbc_p)     ) allocate(lwdnbc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwdnt_p)      ) allocate(lwdnt_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwdntc_p)     ) allocate(lwdntc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwupb_p)      ) allocate(lwupb_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwupbc_p)     ) allocate(lwupbc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwupt_p)      ) allocate(lwupt_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwuptc_p)     ) allocate(lwuptc_p(ims:ime,jms:jme)             )
       if(.not.allocated(olrtoa_p)     ) allocate(olrtoa_p(ims:ime,jms:jme)             )
       if(.not.allocated(sfc_emiss_p)  ) allocate(sfc_emiss_p(ims:ime,jms:jme)          )
       if(.not.allocated(rthratenlw_p) ) allocate(rthratenlw_p(ims:ime,kms:kme,jms:jme) )

       if(.not.allocated(cemiss_p)     ) allocate(cemiss_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(taucldc_p)    ) allocate(taucldc_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(taucldi_p)    ) allocate(taucldi_p(ims:ime,kms:kme,jms:jme)    )

       if(.not.allocated(pin_p)        ) allocate(pin_p(num_oznlevels)                  )
       if(.not.allocated(ozmixm_p) ) &
          allocate(ozmixm_p(ims:ime,1:num_oznlevels,jms:jme,num_months) )
       
       if(.not.allocated(m_hybi_p)     ) allocate(m_hybi_p(num_aerlevels)               )
       if(.not.allocated(m_psn_p)      ) allocate(m_psn_p(ims:ime,jms:jme)              )
       if(.not.allocated(m_psp_p)      ) allocate(m_psp_p(ims:ime,jms:jme)              )
       if(.not.allocated(aerosolcn_p)  ) &
          allocate(aerosolcn_p(ims:ime,1:num_aerlevels,jms:jme,num_aerosols) )
       if(.not.allocated(aerosolcp_p)  ) &
          allocate(aerosolcp_p(ims:ime,1:num_aerlevels,jms:jme,num_aerosols) )

       !allocate these arrays on the first time step, only:
       if(xtime_s .lt. 1.e-12) then

          if(.not.allocated(emstot_p) ) allocate(emstot_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(abstot_p) ) &
             allocate(abstot_p(ims:ime,kms:kme,cam_abs_dim2,jms:jme) )
          if(.not.allocated(absnxt_p) ) &
             allocate(absnxt_p(ims:ime,kms:kme,cam_abs_dim1,jms:jme) )

       endif

    case default

 end select radiation_sw_select

 end subroutine allocate_radiation_sw

!=================================================================================================================
 subroutine deallocate_radiation_sw(configs)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: radt_sw_scheme

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme)
!added by qkyang to support urban model+++++++++++++++++++++
!  if(allocated(hrang_p)      ) deallocate(hrang_p      )
!added by qkyang to support urban model---------------------
 if(allocated(f_ice)        ) deallocate(f_ice        )
 if(allocated(f_rain)       ) deallocate(f_rain       )
 if(allocated(xlat_p)       ) deallocate(xlat_p       )
 if(allocated(xlon_p)       ) deallocate(xlon_p       )
 if(allocated(sfc_albedo_p) ) deallocate(sfc_albedo_p )
 if(allocated(snow_p)       ) deallocate(snow_p       )
 if(allocated(tsk_p)        ) deallocate(tsk_p        )
 if(allocated(xice_p)       ) deallocate(xice_p       )
 if(allocated(xland_p)      ) deallocate(xland_p      )
 if(allocated(coszr_p)      ) deallocate(coszr_p      )
 if(allocated(gsw_p)        ) deallocate(gsw_p        )
 if(allocated(swcf_p)       ) deallocate(swcf_p       )
 if(allocated(swdnb_p)      ) deallocate(swdnb_p      )
 if(allocated(swdnbc_p)     ) deallocate(swdnbc_p     )
 if(allocated(swdnt_p)      ) deallocate(swdnt_p      )
 if(allocated(swdntc_p)     ) deallocate(swdntc_p     )
 if(allocated(swupb_p)      ) deallocate(swupb_p      )
 if(allocated(swupbc_p)     ) deallocate(swupbc_p     )
 if(allocated(swupt_p)      ) deallocate(swupt_p      )
 if(allocated(swuptc_p)     ) deallocate(swuptc_p     )
 
 if(allocated(rthratensw_p) ) deallocate(rthratensw_p )
!jwfeng from USTC AEMOL: 2020/10/09
#ifdef DO_CHEMISTRY
 if(allocated(tauaersw_p)      ) deallocate(tauaersw_p      )
 if(allocated(gaersw_p)        ) deallocate(gaersw_p        )
 if(allocated(waersw_p)        ) deallocate(waersw_p        )
 if(allocated(tauaersw_dust_p) ) deallocate(tauaersw_dust_p )
 if(allocated(gaersw_dust_p)   ) deallocate(gaersw_dust_p   )
 if(allocated(waersw_dust_p)   ) deallocate(waersw_dust_p   )
 if(allocated(tauaersw_seas_p) ) deallocate(tauaersw_seas_p )
 if(allocated(gaersw_seas_p)   ) deallocate(gaersw_seas_p   )
 if(allocated(waersw_seas_p)   ) deallocate(waersw_seas_p   )
 if(allocated(tauaersw_other_p)) deallocate(tauaersw_other_p)
 if(allocated(gaersw_other_p)  ) deallocate(gaersw_other_p  )
 if(allocated(waersw_other_p)  ) deallocate(waersw_other_p  )
 if(allocated(swdnb_dust_p)    ) deallocate(swdnb_dust_p    )
 if(allocated(swdnbc_dust_p)   ) deallocate(swdnbc_dust_p   )
 if(allocated(swdnt_dust_p)    ) deallocate(swdnt_dust_p    )
 if(allocated(swdntc_dust_p)   ) deallocate(swdntc_dust_p   )
 if(allocated(swupb_dust_p)    ) deallocate(swupb_dust_p    )
 if(allocated(swupbc_dust_p)   ) deallocate(swupbc_dust_p   )
 if(allocated(swupt_dust_p)    ) deallocate(swupt_dust_p    )
 if(allocated(swuptc_dust_p)   ) deallocate(swuptc_dust_p   )
 if(allocated(swdnb_seas_p)    ) deallocate(swdnb_seas_p    )
 if(allocated(swdnbc_seas_p)   ) deallocate(swdnbc_seas_p   )
 if(allocated(swdnt_seas_p)    ) deallocate(swdnt_seas_p    )
 if(allocated(swdntc_seas_p)   ) deallocate(swdntc_seas_p   )
 if(allocated(swupb_seas_p)    ) deallocate(swupb_seas_p    )
 if(allocated(swupbc_seas_p)   ) deallocate(swupbc_seas_p   )
 if(allocated(swupt_seas_p)    ) deallocate(swupt_seas_p    )
 if(allocated(swuptc_seas_p)   ) deallocate(swuptc_seas_p   )
 if(allocated(swdnb_other_p)   ) deallocate(swdnb_other_p   )
 if(allocated(swdnbc_other_p)  ) deallocate(swdnbc_other_p  )
 if(allocated(swdnt_other_p)   ) deallocate(swdnt_other_p   )
 if(allocated(swdntc_other_p)  ) deallocate(swdntc_other_p  )
 if(allocated(swupb_other_p)   ) deallocate(swupb_other_p   )
 if(allocated(swupbc_other_p)  ) deallocate(swupbc_other_p  )
 if(allocated(swupt_other_p)   ) deallocate(swupt_other_p   )
 if(allocated(swuptc_other_p)  ) deallocate(swuptc_other_p  )

 if(allocated(rthratensw_dust_p) ) deallocate(rthratensw_dust_p )
 if(allocated(rthratensw_seas_p) ) deallocate(rthratensw_seas_p )
 if(allocated(rthratensw_other_p)) deallocate(rthratensw_other_p)
#endif

 radiation_sw_select: select case (trim(radt_sw_scheme))

    case("rrtmg_sw")
       if(allocated(recloud_p)    ) deallocate(recloud_p    )
       if(allocated(reice_p)      ) deallocate(reice_p      )
       if(allocated(resnow_p)     ) deallocate(resnow_p     )

       if(allocated(alswvisdir_p) ) deallocate(alswvisdir_p )
       if(allocated(alswvisdif_p) ) deallocate(alswvisdif_p )
       if(allocated(alswnirdir_p) ) deallocate(alswnirdir_p )
       if(allocated(alswnirdif_p) ) deallocate(alswnirdif_p )

       if(allocated(swdnflx_p)    ) deallocate(swdnflx_p    )
       if(allocated(swdnflxc_p)   ) deallocate(swdnflxc_p   )
       if(allocated(swupflx_p)    ) deallocate(swupflx_p    )
       if(allocated(swupflxc_p)   ) deallocate(swupflxc_p   )

       if(allocated(pin_p)        ) deallocate(pin_p        )
       if(allocated(o3clim_p)     ) deallocate(o3clim_p     )

    case("cam_sw")
       if(allocated(pin_p)        ) deallocate(pin_p        )
       if(allocated(m_hybi_p)     ) deallocate(m_hybi_p     )

       if(allocated(xlat_p)       ) deallocate(xlat_p       )
       if(allocated(xlon_p)       ) deallocate(xlon_p       )

       if(allocated(glw_p)        ) deallocate(glw_p        )
       if(allocated(lwcf_p)       ) deallocate(lwcf_p       )
       if(allocated(lwdnb_p)      ) deallocate(lwdnb_p      )
       if(allocated(lwdnbc_p)     ) deallocate(lwdnbc_p     )
       if(allocated(lwdnt_p)      ) deallocate(lwdnt_p      )
       if(allocated(lwdntc_p)     ) deallocate(lwdntc_p     )
       if(allocated(lwupb_p)      ) deallocate(lwupb_p      )
       if(allocated(lwupbc_p)     ) deallocate(lwupbc_p     )
       if(allocated(lwupt_p)      ) deallocate(lwupt_p      )
       if(allocated(lwuptc_p)     ) deallocate(lwuptc_p     )
       if(allocated(olrtoa_p)     ) deallocate(olrtoa_p     )
       if(allocated(sfc_emiss_p)  ) deallocate(sfc_emiss_p  )
       if(allocated(rthratenlw_p) ) deallocate(rthratenlw_p )

       if(allocated(cemiss_p)     ) deallocate(cemiss_p     )
       if(allocated(ozmixm_p)     ) deallocate(ozmixm_p     )
       if(allocated(taucldc_p)    ) deallocate(taucldc_p    )
       if(allocated(taucldi_p)    ) deallocate(taucldi_p    )

       if(allocated(m_psn_p)      ) deallocate(m_psn_p      )
       if(allocated(m_psp_p)      ) deallocate(m_psp_p      )
       if(allocated(aerosolcn_p)  ) deallocate(aerosolcn_p  )
       if(allocated(aerosolcp_p)  ) deallocate(aerosolcp_p  )

    case default

 end select radiation_sw_select

 end subroutine deallocate_radiation_sw

!=================================================================================================================
 subroutine radiation_sw_from_iAMAS(configs,mesh,state,time_lev,diag_physics,diag_chemistry,atm_input, &
                                   sfc_input,xtime_s,its,ite)
!=================================================================================================================

!input arguments:
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: state
 type(iamas_pool_type),intent(in):: diag_physics
!jwfeng from USTC AEMOL: 2020/10/09
 type(iamas_pool_type),intent(in):: diag_chemistry
 type(iamas_pool_type),intent(in):: atm_input
 type(iamas_pool_type),intent(in):: sfc_input

 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

 real(kind=RKIND),intent(in):: xtime_s

!local variables:
 integer:: i,j,k,h,n
 character(len=StrKIND):: msg

!local pointers:
 logical,pointer:: config_o3climatology
 character(len=StrKIND),pointer:: radt_sw_scheme
 character(len=StrKIND),pointer:: microp_scheme
 character(len=StrKIND),pointer:: config_chem_opt
 logical,pointer:: config_microp_re

 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: skintemp,snow,xice,xland
 real(kind=RKIND),dimension(:),pointer    :: m_ps,pin
 real(kind=RKIND),dimension(:),pointer    :: sfc_albedo,sfc_emiss
 real(kind=RKIND),dimension(:,:),pointer  :: cldfrac,m_hybi,o3clim
 real(kind=RKIND),dimension(:,:),pointer  :: re_cloud,re_ice,re_snow
 real(kind=RKIND),dimension(:,:,:),pointer:: aerosols,ozmixm
!jwfeng from USTC AEMOL: 2020/10/09
#ifdef DO_CHEMISTRY
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaersw,gaersw,waersw
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaersw_dust,gaersw_dust,waersw_dust
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaersw_seas,gaersw_seas,waersw_seas
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaersw_other,gaersw_other,waersw_other
#endif

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_o3climatology' ,config_o3climatology)
 call iamas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme      )
 call iamas_pool_get_config(configs,'config_microp_scheme' ,microp_scheme       )
 call iamas_pool_get_config(configs,'config_microp_re'     ,config_microp_re    )
 call iamas_pool_get_config(configs,'config_chem_opt'      ,config_chem_opt     )

 call iamas_pool_get_array(mesh,'latCell',latCell)
 call iamas_pool_get_array(mesh,'lonCell',lonCell)

 call iamas_pool_get_array(state,'aerosols',aerosols,time_lev)

 call iamas_pool_get_array(sfc_input,'skintemp',skintemp)
 call iamas_pool_get_array(sfc_input,'snow'    ,snow    )
 call iamas_pool_get_array(sfc_input,'xice'    ,xice    )
 call iamas_pool_get_array(sfc_input,'xland'   ,xland   )

 call iamas_pool_get_array(atm_input,'pin'     ,pin     )
 call iamas_pool_get_array(atm_input,'ozmixm'  ,ozmixm  )

 call iamas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call iamas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call iamas_pool_get_array(diag_physics,'cldfrac'   ,cldfrac   )
 call iamas_pool_get_array(diag_physics,'o3clim'    ,o3clim    )
 call iamas_pool_get_array(diag_physics,'m_hybi'    ,m_hybi    )
 call iamas_pool_get_array(diag_physics,'m_ps'      ,m_ps      )

!jwfeng from USTC AEMOL: 2020/10/09
#ifdef DO_CHEMISTRY
 if(trim(config_chem_opt) .ne. "no_chemistry") then
    call iamas_pool_get_array(diag_chemistry,'tauaersw',tauaersw)
    call iamas_pool_get_array(diag_chemistry,'gaersw'  ,gaersw  )
    call iamas_pool_get_array(diag_chemistry,'waersw'  ,waersw  )
    call iamas_pool_get_array(diag_chemistry,'tauaersw_dust',tauaersw_dust)
    call iamas_pool_get_array(diag_chemistry,'gaersw_dust'  ,gaersw_dust  )
    call iamas_pool_get_array(diag_chemistry,'waersw_dust'  ,waersw_dust  )
    call iamas_pool_get_array(diag_chemistry,'tauaersw_seas',tauaersw_seas)
    call iamas_pool_get_array(diag_chemistry,'gaersw_seas'  ,gaersw_seas  )
    call iamas_pool_get_array(diag_chemistry,'waersw_seas'  ,waersw_seas  )
 endif
#endif

 do j = jts,jte
 do i = its,ite
    xlat_p(i,j)       = latCell(i) / degrad
    xlon_p(i,j)       = lonCell(i) / degrad

    sfc_albedo_p(i,j) = sfc_albedo(i)
    snow_p(i,j)       = snow(i)
    tsk_p(i,j)        = skintemp(i)
    xice_p(i,j)       = xice(i)
    xland_p(i,j)      = xland(i)
 enddo
 enddo
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac_p(i,k,j) = cldfrac(k,i)
 enddo
 enddo
 enddo

!jwfeng from USTC AEMOL: 2020/10/09
#ifdef DO_CHEMISTRY
 if(trim(config_chem_opt) .ne. "no_chemistry") then
 do h = 1,14
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
     tauaersw_p(i,k,j,h)       = tauaersw(h,k,i)
     gaersw_p(i,k,j,h)         = gaersw(h,k,i)
     waersw_p(i,k,j,h)         = waersw(h,k,i)
     tauaersw_dust_p(i,k,j,h)  = tauaersw_dust(h,k,i)
     gaersw_dust_p(i,k,j,h)    = gaersw_dust(h,k,i)
     waersw_dust_p(i,k,j,h)    = waersw_dust(h,k,i)
     tauaersw_seas_p(i,k,j,h)  = tauaersw_seas(h,k,i)
     gaersw_seas_p(i,k,j,h)    = gaersw_seas(h,k,i)
     waersw_seas_p(i,k,j,h)    = waersw_seas(h,k,i)
     tauaersw_other_p(i,k,j,h) = 0.0_RKIND
     gaersw_other_p(i,k,j,h)   = 0.0_RKIND
     waersw_other_p(i,k,j,h)   = 1.0_RKIND
 enddo
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    swdnb_other_p(i,j)    = 0.0_RKIND
    swdnbc_other_p(i,j)   = 0.0_RKIND
    swdnt_other_p(i,j)    = 0.0_RKIND
    swdntc_other_p(i,j)   = 0.0_RKIND
    swupb_other_p(i,j)    = 0.0_RKIND
    swupbc_other_p(i,j)   = 0.0_RKIND
    swupt_other_p(i,j)    = 0.0_RKIND
    swuptc_other_p(i,j)   = 0.0_RKIND
    swdnb_dust_p(i,j)     = 0.0_RKIND
    swdnbc_dust_p(i,j)    = 0.0_RKIND
    swdnt_dust_p(i,j)     = 0.0_RKIND
    swdntc_dust_p(i,j)    = 0.0_RKIND
    swupb_dust_p(i,j)     = 0.0_RKIND
    swupbc_dust_p(i,j)    = 0.0_RKIND
    swupt_dust_p(i,j)     = 0.0_RKIND
    swuptc_dust_p(i,j)    = 0.0_RKIND
    swdnb_seas_p(i,j)     = 0.0_RKIND
    swdnbc_seas_p(i,j)    = 0.0_RKIND
    swdnt_seas_p(i,j)     = 0.0_RKIND
    swdntc_seas_p(i,j)    = 0.0_RKIND
    swupb_seas_p(i,j)     = 0.0_RKIND
    swupbc_seas_p(i,j)    = 0.0_RKIND
    swupt_seas_p(i,j)     = 0.0_RKIND
    swuptc_seas_p(i,j)    = 0.0_RKIND
 enddo

 do k = kts,kte
 do i = its,ite
    rthratensw_dust_p(i,k,j)  = 0.0_RKIND
    rthratensw_seas_p(i,k,j)  = 0.0_RKIND
    rthratensw_other_p(i,k,j) = 0.0_RKIND
 enddo
 enddo
 enddo
 endif
#endif

!initialization:
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    f_ice(i,k,j)  = 0.0_RKIND
    f_rain(i,k,j) = 0.0_RKIND
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    coszr_p(i,j)    = 0.0_RKIND
    !added by qkyang to support urban model+++++++++++++++++++++
!     hrang_p(i,j)    = 0.0_RKIND
    !added by qkyang to support urban model---------------------
    gsw_p(i,j)      = 0.0_RKIND
    swcf_p(i,j)     = 0.0_RKIND
    swdnb_p(i,j)    = 0.0_RKIND
    swdnbc_p(i,j)   = 0.0_RKIND
    swdnt_p(i,j)    = 0.0_RKIND
    swdntc_p(i,j)   = 0.0_RKIND
    swupb_p(i,j)    = 0.0_RKIND
    swupbc_p(i,j)   = 0.0_RKIND
    swupt_p(i,j)    = 0.0_RKIND
    swuptc_p(i,j)   = 0.0_RKIND
 enddo

 do k = kts,kte
 do i = its,ite
    rthratensw_p(i,k,j) = 0.0_RKIND
 enddo
 enddo
 enddo

 radiation_sw_select: select case (trim(radt_sw_scheme))

    case("rrtmg_sw")

       microp_select: select case(microp_scheme)
          case("mp_thompson","mp_wsm6")
             if(config_microp_re) then
                call iamas_pool_get_array(diag_physics,'re_cloud',re_cloud)
                call iamas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
                call iamas_pool_get_array(diag_physics,'re_snow' ,re_snow )

                do j = jts,jte
                do k = kts,kte
                do i = its,ite
                   recloud_p(i,k,j) = re_cloud(k,i)
                   reice_p(i,k,j)   = re_ice(k,i)
                   resnow_p(i,k,j)  = re_snow(k,i)
                enddo
                enddo
                enddo
             else
                ! These are set in module iamas_atmphys_manager and should not be set again
                !has_reqc = 0
                !has_reqi = 0
                !has_reqs = 0
                do j = jts,jte
                do k = kts,kte
                do i = its,ite
                   recloud_p(i,k,j) = 0._RKIND
                   reice_p(i,k,j)   = 0._RKIND
                   resnow_p(i,k,j)  = 0._RKIND
                enddo
                enddo
                enddo
             endif

          case default          
       end select microp_select

       do j = jts,jte
       do k = kts,kte+2
       do i = its,ite
          swdnflx_p(i,k,j)  = 0.0_RKIND
          swdnflxc_p(i,k,j) = 0.0_RKIND
          swupflx_p(i,k,j)  = 0.0_RKIND
          swupflxc_p(i,k,j) = 0.0_RKIND
       enddo
       enddo
       enddo

       !ozone volume mixing ratio:
       if(config_o3climatology) then
          do k = 1, num_oznLevels
             pin_p(k) = pin(k)
          enddo
          do j = jts,jte
          do k = 1, num_oznLevels
             do i = its,ite
                o3clim_p(i,k,j) = o3clim(k,i)
             enddo
          enddo
          enddo
       else
          do k = 1, num_oznLevels
             pin_p(k) = 0.0_RKIND
          enddo
          do j = jts,jte
          do k = 1, num_oznLevels
             do i = its,ite
                o3clim_p(i,k,j) = 0.0_RKIND
             enddo
          enddo
          enddo
       endif

    case("cam_sw")
       do j = jts,jte
       do i = its,ite
          sfc_emiss_p(i,j)  = sfc_emiss(i)

          olrtoa_p(i,j)     = 0.0_RKIND
          glw_p(i,j)        = 0.0_RKIND
          lwcf_p(i,j)       = 0.0_RKIND
          lwdnb_p(i,j)      = 0.0_RKIND
          lwdnbc_p(i,j)     = 0.0_RKIND
          lwdnt_p(i,j)      = 0.0_RKIND
          lwdntc_p(i,j)     = 0.0_RKIND
          lwupb_p(i,j)      = 0.0_RKIND
          lwupbc_p(i,j)     = 0.0_RKIND
          lwupt_p(i,j)      = 0.0_RKIND
          lwuptc_p(i,j)     = 0.0_RKIND
       enddo
       do k = kts,kte
       do i = its,ite
          rthratenlw_p(i,k,j) = 0.0_RKIND
          cemiss_p(i,k,j)     = 0.0_RKIND
          taucldc_p(i,k,j)    = 0.0_RKIND
          taucldi_p(i,k,j)    = 0.0_RKIND
       enddo
       enddo
       enddo
       !infrared absorption:
       if(xtime_s .lt. 1.e-12) then
          do j = jts,jte
          do n = 1,cam_abs_dim1
          do k = kts,kte
          do i = its,ite
             absnxt_p(i,k,n,j) = 0.0_RKIND
          enddo
          enddo
          enddo
          do n = 1,cam_abs_dim2
          do k = kts,kte+1
          do i = its,ite
             abstot_p(i,k,n,j) = 0.0_RKIND
          enddo
          enddo
          enddo
          do k = kts,kte+1
          do i = its,ite
             emstot_p(i,k,j) = 0.0_RKIND
          enddo
          enddo
          enddo
       endif
       !ozone mixing ratio:
       do k = 1, num_oznlevels
          pin_p(k) = pin(k)
       enddo
       do n = 1, num_months
          do j = jts,jte
          do k = 1, num_oznlevels
          do i = its,ite
             ozmixm_p(i,k,j,n) = ozmixm(n,k,i)
          enddo
          enddo
          enddo
       enddo
       !aerosol mixing ratio:
       do k = 1, num_aerlevels
          m_hybi_p(k) = m_hybi(k,1)
       enddo
       do i = its,ite
       do j = jts,jte
          m_psp_p(i,j) = m_ps(i)
          m_psn_p(i,j) = m_ps(i)
       enddo
       enddo
       do n = 1,num_aerosols
       do j = jts,jte
       do k = 1, num_aerlevels
       do i = its,ite
          aerosolcp_p(i,k,j,n) = aerosols(n,k,i)
          aerosolcn_p(i,k,j,n) = aerosols(n,k,i)
       enddo
       enddo
       enddo
       enddo

    case default

 end select radiation_sw_select

 end subroutine radiation_sw_from_iAMAS

!=================================================================================================================
 subroutine radiation_sw_to_iAMAS(diag_physics,tend_physics,its,ite)
!=================================================================================================================

!Added by qkyang to support urban model+++++++++++++++++++++
 real(kind=RKIND),pointer  :: declin
!  real(kind=RKIND),dimension(:),pointer  :: hrang
!Added by qkyang to support urban model--------------------- 

!input arguments:
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: tend_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,j,k,n

!local pointers:
 real(kind=RKIND),dimension(:),pointer  :: coszr,gsw,swcf,swdnb,swdnbc,swdnt,swdntc, &
                                           swupb,swupbc,swupt,swuptc
 real(kind=RKIND),dimension(:,:),pointer:: rthratensw
#ifdef DO_CHEMISTRY
 real(kind=RKIND),dimension(:),pointer  :: swdnb_dust,swdnbc_dust,swdnt_dust,swdntc_dust, &
                                           swupb_dust,swupbc_dust,swupt_dust,swuptc_dust
 real(kind=RKIND),dimension(:),pointer  :: swdnb_seas,swdnbc_seas,swdnt_seas,swdntc_seas, &
                                           swupb_seas,swupbc_seas,swupt_seas,swuptc_seas
 real(kind=RKIND),dimension(:),pointer  :: swdnb_other,swdnbc_other,swdnt_other,swdntc_other, &
                                           swupb_other,swupbc_other,swupt_other,swuptc_other
 real(kind=RKIND),dimension(:,:),pointer:: rthratensw_dust,rthratensw_seas,rthratensw_other
#endif

!-----------------------------------------------------------------------------------------------------------------

!Added by qkyang to support urban model+++++++++++++++++++++
 call iamas_pool_get_array(diag_physics,'declin'    ,declin    )
!  call iamas_pool_get_array(diag_physics,'hrang'     ,hrang     )
!Added by qkyang to support urban model--------------------- 
 call iamas_pool_get_array(diag_physics,'coszr'     ,coszr     )
 call iamas_pool_get_array(diag_physics,'gsw'       ,gsw       )
 call iamas_pool_get_array(diag_physics,'swcf'      ,swcf      )
 call iamas_pool_get_array(diag_physics,'swdnb'     ,swdnb     )
 call iamas_pool_get_array(diag_physics,'swdnbc'    ,swdnbc    )
 call iamas_pool_get_array(diag_physics,'swdnt'     ,swdnt     )
 call iamas_pool_get_array(diag_physics,'swdntc'    ,swdntc    )
 call iamas_pool_get_array(diag_physics,'swupb'     ,swupb     )
 call iamas_pool_get_array(diag_physics,'swupbc'    ,swupbc    )
 call iamas_pool_get_array(diag_physics,'swupt'     ,swupt     )
 call iamas_pool_get_array(diag_physics,'swuptc'    ,swuptc    )
 call iamas_pool_get_array(tend_physics,'rthratensw',rthratensw)

#ifdef DO_CHEMISTRY
 call iamas_pool_get_array(diag_physics,'swdnb_seas' ,swdnb_seas )
 call iamas_pool_get_array(diag_physics,'swdnbc_seas',swdnbc_seas)
 call iamas_pool_get_array(diag_physics,'swdnt_seas' ,swdnt_seas )
 call iamas_pool_get_array(diag_physics,'swdntc_seas',swdntc_seas)
 call iamas_pool_get_array(diag_physics,'swupb_seas' ,swupb_seas )
 call iamas_pool_get_array(diag_physics,'swupbc_seas',swupbc_seas)
 call iamas_pool_get_array(diag_physics,'swupt_seas' ,swupt_seas )
 call iamas_pool_get_array(diag_physics,'swuptc_seas',swuptc_seas)
 call iamas_pool_get_array(tend_physics,'rthratensw_seas',rthratensw_seas)
 call iamas_pool_get_array(diag_physics,'swdnb_dust' ,swdnb_dust )
 call iamas_pool_get_array(diag_physics,'swdnbc_dust',swdnbc_dust)
 call iamas_pool_get_array(diag_physics,'swdnt_dust' ,swdnt_dust )
 call iamas_pool_get_array(diag_physics,'swdntc_dust',swdntc_dust)
 call iamas_pool_get_array(diag_physics,'swupb_dust' ,swupb_dust )
 call iamas_pool_get_array(diag_physics,'swupbc_dust',swupbc_dust)
 call iamas_pool_get_array(diag_physics,'swupt_dust' ,swupt_dust )
 call iamas_pool_get_array(diag_physics,'swuptc_dust',swuptc_dust)
 call iamas_pool_get_array(tend_physics,'rthratensw_dust',rthratensw_dust)
 call iamas_pool_get_array(diag_physics,'swdnb_other' ,swdnb_other )
 call iamas_pool_get_array(diag_physics,'swdnbc_other',swdnbc_other)
 call iamas_pool_get_array(diag_physics,'swdnt_other' ,swdnt_other )
 call iamas_pool_get_array(diag_physics,'swdntc_other',swdntc_other)
 call iamas_pool_get_array(diag_physics,'swupb_other' ,swupb_other )
 call iamas_pool_get_array(diag_physics,'swupbc_other',swupbc_other)
 call iamas_pool_get_array(diag_physics,'swupt_other' ,swupt_other )
 call iamas_pool_get_array(diag_physics,'swuptc_other',swuptc_other)
 call iamas_pool_get_array(tend_physics,'rthratensw_other',rthratensw_other)
#endif

 do j = jts,jte

 do i = its,ite
    coszr(i) = coszr_p(i,j)
!     hrang(i) = hrang_p(i,j)
    gsw(i)    = gsw_p(i,j)
    swcf(i)   = swcf_p(i,j)
    swdnb(i)  = swdnb_p(i,j)
    swdnbc(i) = swdnbc_p(i,j)
    swdnt(i)  = swdnt_p(i,j)
    swdntc(i) = swdntc_p(i,j)
    swupb(i)  = swupb_p(i,j)
    swupbc(i) = swupbc_p(i,j)
    swupt(i)  = swupt_p(i,j)
    swuptc(i) = swuptc_p(i,j)
 enddo

 do k = kts,kte
 do i = its,ite
    rthratensw(k,i) = rthratensw_p(i,k,j)
 enddo
 enddo

 enddo

#ifdef DO_CHEMISTRY
 do j = jts,jte
 do i = its,ite
    swdnb_dust(i)  = swdnb_dust_p(i,j)
    swdnbc_dust(i) = swdnbc_dust_p(i,j)
    swdnt_dust(i)  = swdnt_dust_p(i,j)
    swdntc_dust(i) = swdntc_dust_p(i,j)
    swupb_dust(i)  = swupb_dust_p(i,j)
    swupbc_dust(i) = swupbc_dust_p(i,j)
    swupt_dust(i)  = swupt_dust_p(i,j)
    swuptc_dust(i) = swuptc_dust_p(i,j)
    swdnb_seas(i)  = swdnb_seas_p(i,j)
    swdnbc_seas(i) = swdnbc_seas_p(i,j)
    swdnt_seas(i)  = swdnt_seas_p(i,j)
    swdntc_seas(i) = swdntc_seas_p(i,j)
    swupb_seas(i)  = swupb_seas_p(i,j)
    swupbc_seas(i) = swupbc_seas_p(i,j)
    swupt_seas(i)  = swupt_seas_p(i,j)
    swuptc_seas(i) = swuptc_seas_p(i,j)
    swdnb_other(i)  = swdnb_other_p(i,j)
    swdnbc_other(i) = swdnbc_other_p(i,j)
    swdnt_other(i)  = swdnt_other_p(i,j)
    swdntc_other(i) = swdntc_other_p(i,j)
    swupb_other(i)  = swupb_other_p(i,j)
    swupbc_other(i) = swupbc_other_p(i,j)
    swupt_other(i)  = swupt_other_p(i,j)
    swuptc_other(i) = swuptc_other_p(i,j)
 enddo
 do k = kts,kte
 do i = its,ite
    rthratensw_dust(k,i)  = rthratensw_dust_p(i,k,j)
    rthratensw_seas(k,i)  = rthratensw_seas_p(i,k,j)
    rthratensw_other(k,i) = rthratensw_other_p(i,k,j)
 enddo
 enddo
 enddo
#endif
 end subroutine radiation_sw_to_iAMAS

!=================================================================================================================
 subroutine init_radiation_sw(dminfo,configs,mesh,atm_input,diag,diag_physics,state,time_lev)
!=================================================================================================================

!input arguments:
 type(dm_info), intent(in):: dminfo
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in),optional:: mesh
 type(iamas_pool_type),intent(in),optional:: diag
 type(iamas_pool_type),intent(in),optional:: diag_physics

 integer,intent(in),optional:: time_lev

!inout arguments:
 type(iamas_pool_type),intent(inout),optional:: atm_input
 type(iamas_pool_type),intent(inout),optional:: state

!local pointers:
 character(len=StrKIND),pointer:: radt_sw_scheme

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme)

!call to shortwave radiation scheme:
 radiation_sw_select: select case (trim(radt_sw_scheme))

    case ("rrtmg_sw")
       call rrtmg_initsw_foriAMAS(dminfo)

    case("cam_sw")
       call camradinit(dminfo,mesh,atm_input,diag,diag_physics,state,time_lev)

    case default

 end select radiation_sw_select

 end subroutine init_radiation_sw

!=================================================================================================================
 subroutine driver_radiation_sw(itimestep,configs,mesh,state,time_lev,diag_physics,diag_chemistry,atm_input, &
                                sfc_input,tend_physics,xtime_s,its,ite)
!=================================================================================================================

!Added by qkyang to support urban model+++++++++++++++++++++
!  real(kind=RKIND),pointer  :: declin
 real(kind=RKIND),dimension(:),pointer  :: hrang
!Added by qkyang to support urban model--------------------- 
!Added by qkyang to support slope and shading effects>>>>>>>>>>>>>>>>>
integer,pointer:: config_slope_rad
integer,pointer:: config_topo_shading
real(kind=RKIND),pointer:: config_shadlen
!Added by qkyang to support slope and shading effects<<<<<<<<<<<<<<<<
!input arguments:
 integer,intent(in):: itimestep
 type(iamas_pool_type),intent(in):: configs
 type(iamas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(iamas_pool_type),intent(inout):: state
 type(iamas_pool_type),intent(inout):: diag_physics
!jwfeng from USTC AEMOL: 2020/10/09
 type(iamas_pool_type),intent(inout):: diag_chemistry
 type(iamas_pool_type),intent(inout):: atm_input
 type(iamas_pool_type),intent(inout):: sfc_input
 type(iamas_pool_type),intent(inout):: tend_physics


!local pointers:
 logical,pointer:: config_o3climatology
 character(len=StrKIND),pointer:: radt_sw_scheme
 character(len=StrKIND),pointer:: config_chem_opt
 integer,pointer:: config_aer_opt 
 logical,pointer:: config_aer_ra_feedback 

!local variables:
 integer:: o3input
 real(kind=RKIND):: radt,xtime_m

!Added by qkyang to support urban model+++++++++++++++++++++
!  call iamas_pool_get_array(diag_physics,'declin'    ,declin    )
 call iamas_pool_get_array(diag_physics,'hrang'     ,hrang     )
!Added by qkyang to support urban model--------------------- 

!Added by qkyang to support slope and shading effects>>>>>>>>>>>>>>>>>
 call iamas_pool_get_config(configs,'config_slope_rad',config_slope_rad)
 call iamas_pool_get_config(configs,'config_topo_shading',config_topo_shading)
 call iamas_pool_get_config(configs,'config_shadlen',config_shadlen)
!Added by qkyang to support slope and shading effects<<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------------------------------------------------
!call iamas_log_write(' --- enter subroutine driver_radiation_sw: $i',intArgs=(/itimestep/))

 call iamas_pool_get_config(configs,'config_o3climatology' ,config_o3climatology)
 call iamas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme      )
 call iamas_pool_get_config(configs,'config_aer_opt',config_aer_opt    )
 call iamas_pool_get_config(configs,'config_aer_ra_feedback',config_aer_ra_feedback    )
 call iamas_pool_get_config(configs,'config_chem_opt',config_chem_opt    )

 xtime_m = xtime_s/60.

!copy iAMAS arrays to local arrays:
 call radiation_sw_from_iAMAS(configs,mesh,state,time_lev,diag_physics,diag_chemistry,atm_input,sfc_input,xtime_s,its,ite)

! This should be OMP MASTER with barrier afterwards or OMP SINGLE, since declin and solcon
! are global variables in iamas_atmphys_vars.F and race conditions may occur otherwise!
!$OMP SINGLE
!... calculates solar declination:
!call radconst(declin,solcon,julday,degrad,dpd)
 call radconst(declin,solcon,curr_julday,degrad,dpd)
!call iamas_log_write('     ITIMESTEP   = $i', intArgs=(/itimestep/))
!call iamas_log_write('     YEAR        = $i', intArgs=(/year/))
!call iamas_log_write('     JULDAY      = $i', intArgs=(/julday/))
!call iamas_log_write('     GMT         = $r', realArgs=(/gmt/))
!call iamas_log_write('     XTIME_M     = $r', realArgs=(/xtime_m/))
!call iamas_log_write('     CURR_JULDAY = $r', realArgs=(/curr_julday/))
!call iamas_log_write('     SOLCON      = $r', realArgs=(/solcon/))
!call iamas_log_write('     DECLIN      = $r', realArgs=(/declin/))
!$OMP END SINGLE

!... convert the radiation time_step to minutes:
 radt = dt_radtsw/60.
    
!call to shortwave radiation scheme:
 radiation_sw_select: select case (trim(radt_sw_scheme))

    case ("rrtmg_sw")
       o3input = 0
       if(config_o3climatology) o3input = 2

       call iamas_timer_start('RRTMG_sw')
#ifdef DO_CHEMISTRY
    if(trim(config_chem_opt) .ne. "no_chemistry" .and. config_aer_ra_feedback) then
       if (ncomp_aer(1) .ge. 2) then
       call rrtmg_swrad( &
              p3d        = pres_hyd_p   , p8w        = pres2_hyd_p   , pi3d     = pi_p     , &
              t3d        = t_p          , t8w        = t2_p          , dz8w     = dz_p     , &
!             qv3d       = qv_p         , qc3d       = qc_p          , qi3d     = qi_p     , &
!             qs3d       = qs_p         , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              qv3d       = qvrad_p      , qc3d       = qcrad_p       , qi3d     = qirad_p  , &
              qs3d       = qsrad_p      , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              albedo     = sfc_albedo_p , xland      = xland_p       , xice     = xice_p   , &
              snow       = snow_p       , coszr      = coszr_p       , xtime    = xtime_m  , &
              hrang      = hrang(its:ite)      , &!qkyang added hrang   =hrang(its:ite)     ,
              gmt        = gmt          , julday     = julday        , radt     = radt     , &
              degrad     = degrad       , declin     = declin        , solcon   = solcon   , &
              xlat       = xlat_p       , xlong      = xlon_p        , icloud   = icloud   , &
              o3input    = o3input      , noznlevels = num_oznlevels , pin      = pin_p    , &
              o3clim     = o3clim_p     , gsw        = gsw_p         , swcf     = swcf_p   , &
              config_aer_opt = config_aer_opt, config_aer_ra_feedback = config_aer_ra_feedback,&
              config_chem_opt = config_chem_opt                                            , &
              rthratensw = rthratensw_other_p,has_reqc   = has_reqc  , has_reqi = has_reqi , &
              has_reqs   = has_reqs     , re_cloud   = recloud_p     , re_ice   = reice_p  , &
              re_snow    = resnow_p     , swupt      = swupt_other_p , swuptc   = swuptc_other_p , &
              swdnt      = swdnt_other_p , swdntc    = swdntc_other_p, swupb    = swupb_other_p  , &
              swupbc     = swupbc_other_p, swdnb     = swdnb_other_p , swdnbc   = swdnbc_other_p , &
              tauaersw   = tauaersw_other_p , gaersw = gaersw_other_p, waersw   = waersw_other_p , & !jwfeng
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,        &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,        &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                       )
       end if
       call rrtmg_swrad( &
              p3d        = pres_hyd_p   , p8w        = pres2_hyd_p   , pi3d     = pi_p     , &
              t3d        = t_p          , t8w        = t2_p          , dz8w     = dz_p     , &
!             qv3d       = qv_p         , qc3d       = qc_p          , qi3d     = qi_p     , &
!             qs3d       = qs_p         , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              qv3d       = qvrad_p      , qc3d       = qcrad_p       , qi3d     = qirad_p  , &
              qs3d       = qsrad_p      , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              albedo     = sfc_albedo_p , xland      = xland_p       , xice     = xice_p   , &
              snow       = snow_p       , coszr      = coszr_p       , xtime    = xtime_m  , &
              hrang      = hrang(its:ite)      , &!qkyang added hrang   =hrang(its:ite)     ,
              gmt        = gmt          , julday     = julday        , radt     = radt     , &
              degrad     = degrad       , declin     = declin        , solcon   = solcon   , &
              xlat       = xlat_p       , xlong      = xlon_p        , icloud   = icloud   , &
              o3input    = o3input      , noznlevels = num_oznlevels , pin      = pin_p    , &
              o3clim     = o3clim_p     , gsw        = gsw_p         , swcf     = swcf_p   , &
              config_aer_opt = config_aer_opt, config_aer_ra_feedback = config_aer_ra_feedback,&
              config_chem_opt = config_chem_opt                                            , &
              rthratensw = rthratensw_p , has_reqc   = has_reqc      , has_reqi = has_reqi , &
              has_reqs   = has_reqs     , re_cloud   = recloud_p     , re_ice   = reice_p  , &
              re_snow    = resnow_p     , swupt      = swupt_p       , swuptc   = swuptc_p , &
              swdnt      = swdnt_p      , swdntc     = swdntc_p      , swupb    = swupb_p  , &
              swupbc     = swupbc_p     , swdnb      = swdnb_p       , swdnbc   = swdnbc_p , &
              tauaersw   = tauaersw_p   , gaersw     = gaersw_p      , waersw   = waersw_p , & !jwfeng
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,        &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,        &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                       )
       ! jwfeng: 2022/07/25 Calculate optical params of dust aerosols.
       if (lptr_dust_aer(1,1,1) .ge. 1) then
       call rrtmg_swrad( &
              p3d        = pres_hyd_p   , p8w        = pres2_hyd_p   , pi3d     = pi_p     , &
              t3d        = t_p          , t8w        = t2_p          , dz8w     = dz_p     , &
!             qv3d       = qv_p         , qc3d       = qc_p          , qi3d     = qi_p     , &
!             qs3d       = qs_p         , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              qv3d       = qvrad_p      , qc3d       = qcrad_p       , qi3d     = qirad_p  , &
              qs3d       = qsrad_p      , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              albedo     = sfc_albedo_p , xland      = xland_p       , xice     = xice_p   , &
              snow       = snow_p       , coszr      = coszr_p       , xtime    = xtime_m  , &
              hrang      = hrang(its:ite)      , &!qkyang added hrang   =hrang(its:ite)     ,
              gmt        = gmt          , julday     = julday        , radt     = radt     , &
              degrad     = degrad       , declin     = declin        , solcon   = solcon   , &
              xlat       = xlat_p       , xlong      = xlon_p        , icloud   = icloud   , &
              o3input    = o3input      , noznlevels = num_oznlevels , pin      = pin_p    , &
              o3clim     = o3clim_p     , gsw        = gsw_p         , swcf     = swcf_p   , &
              config_aer_opt = config_aer_opt, config_aer_ra_feedback = config_aer_ra_feedback,&
              config_chem_opt = config_chem_opt                                            , &
              rthratensw = rthratensw_dust_p , has_reqc = has_reqc   , has_reqi = has_reqi , &
              has_reqs = has_reqs       ,re_cloud = recloud_p     , re_ice = reice_p       , &
              re_snow  = resnow_p       ,swupt    = swupt_dust_p  , swuptc = swuptc_dust_p , &
              swdnt    = swdnt_dust_p   ,swdntc   = swdntc_dust_p , swupb  = swupb_dust_p  , &
              swupbc   = swupbc_dust_p  ,swdnb    = swdnb_dust_p  , swdnbc = swdnbc_dust_p , &
              tauaersw = tauaersw_dust_p,gaersw   = gaersw_dust_p , waersw = waersw_dust_p , & !jwfeng
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,        &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,        &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                       )
       end if
       ! jwfeng: 2022/07/25 Calculate optical params of seas aerosols.
       if (lptr_cl_aer(1,1,1) .ge. 1 .and. lptr_na_aer(1,1,1) .ge. 1) then
       call rrtmg_swrad( &
              p3d        = pres_hyd_p   , p8w        = pres2_hyd_p   , pi3d     = pi_p     , &
              t3d        = t_p          , t8w        = t2_p          , dz8w     = dz_p     , &
!             qv3d       = qv_p         , qc3d       = qc_p          , qi3d     = qi_p     , &
!             qs3d       = qs_p         , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              qv3d       = qvrad_p      , qc3d       = qcrad_p       , qi3d     = qirad_p  , &
              qs3d       = qsrad_p      , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              albedo     = sfc_albedo_p , xland      = xland_p       , xice     = xice_p   , &
              snow       = snow_p       , coszr      = coszr_p       , xtime    = xtime_m  , &
              hrang      = hrang(its:ite)      , &!qkyang added hrang   =hrang(its:ite)     ,
              gmt        = gmt          , julday     = julday        , radt     = radt     , &
              degrad     = degrad       , declin     = declin        , solcon   = solcon   , &
              xlat       = xlat_p       , xlong      = xlon_p        , icloud   = icloud   , &
              o3input    = o3input      , noznlevels = num_oznlevels , pin      = pin_p    , &
              o3clim     = o3clim_p     , gsw        = gsw_p         , swcf     = swcf_p   , &
              config_aer_opt = config_aer_opt, config_aer_ra_feedback = config_aer_ra_feedback,&
              config_chem_opt = config_chem_opt                                            , &
              rthratensw = rthratensw_seas_p , has_reqc = has_reqc   , has_reqi = has_reqi , &
              has_reqs = has_reqs       ,re_cloud = recloud_p     , re_ice = reice_p       , &
              re_snow  = resnow_p       ,swupt    = swupt_seas_p  , swuptc = swuptc_seas_p , &
              swdnt    = swdnt_seas_p   ,swdntc   = swdntc_seas_p , swupb  = swupb_seas_p  , &
              swupbc   = swupbc_seas_p  ,swdnb    = swdnb_seas_p  , swdnbc = swdnbc_seas_p , &
              tauaersw = tauaersw_seas_p,gaersw   = gaersw_seas_p , waersw = waersw_seas_p , & !jwfeng
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,        &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,        &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                       )
       end if
    else
       call rrtmg_swrad( &
              p3d        = pres_hyd_p   , p8w        = pres2_hyd_p   , pi3d     = pi_p     , &
              t3d        = t_p          , t8w        = t2_p          , dz8w     = dz_p     , &
!             qv3d       = qv_p         , qc3d       = qc_p          , qi3d     = qi_p     , &
!             qs3d       = qs_p         , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              qv3d       = qvrad_p      , qc3d       = qcrad_p       , qi3d     = qirad_p  , &
              qs3d       = qsrad_p      , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              albedo     = sfc_albedo_p , xland      = xland_p       , xice     = xice_p   , &
              snow       = snow_p       , coszr      = coszr_p       , xtime    = xtime_m  , &
              hrang      = hrang(its:ite)      , &!qkyang added hrang   =hrang(its:ite)     ,
              gmt        = gmt          , julday     = julday        , radt     = radt     , &
              degrad     = degrad       , declin     = declin        , solcon   = solcon   , &
              xlat       = xlat_p       , xlong      = xlon_p        , icloud   = icloud   , &
              o3input    = o3input      , noznlevels = num_oznlevels , pin      = pin_p    , &
              o3clim     = o3clim_p     , gsw        = gsw_p         , swcf     = swcf_p   , &
              config_aer_opt = config_aer_opt, config_aer_ra_feedback = config_aer_ra_feedback,&
              config_chem_opt = config_chem_opt                                            , &
              rthratensw = rthratensw_p , has_reqc   = has_reqc      , has_reqi = has_reqi , &
              has_reqs   = has_reqs     , re_cloud   = recloud_p     , re_ice   = reice_p  , &
              re_snow    = resnow_p     , swupt      = swupt_p       , swuptc   = swuptc_p , &
              swdnt      = swdnt_p      , swdntc     = swdntc_p      , swupb    = swupb_p  , &
              swupbc     = swupbc_p     , swdnb      = swdnb_p       , swdnbc   = swdnbc_p , &
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,        &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,        &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                       )
    endif
#else
       call rrtmg_swrad( &
              p3d        = pres_hyd_p   , p8w        = pres2_hyd_p   , pi3d     = pi_p     , &
              t3d        = t_p          , t8w        = t2_p          , dz8w     = dz_p     , &
!             qv3d       = qv_p         , qc3d       = qc_p          , qi3d     = qi_p     , &
!             qs3d       = qs_p         , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              qv3d       = qvrad_p      , qc3d       = qcrad_p       , qi3d     = qirad_p  , &
              qs3d       = qsrad_p      , cldfra3d   = cldfrac_p     , tsk      = tsk_p    , &
              albedo     = sfc_albedo_p , xland      = xland_p       , xice     = xice_p   , &
              snow       = snow_p       , coszr      = coszr_p       , xtime    = xtime_m  , &
              hrang      = hrang(its:ite)      , &!qkyang added hrang   =hrang(its:ite)     ,
              gmt        = gmt          , julday     = julday        , radt     = radt     , &
              degrad     = degrad       , declin     = declin        , solcon   = solcon   , &
              xlat       = xlat_p       , xlong      = xlon_p        , icloud   = icloud   , &
              o3input    = o3input      , noznlevels = num_oznlevels , pin      = pin_p    , &
              o3clim     = o3clim_p     , gsw        = gsw_p         , swcf     = swcf_p   , &
              config_aer_opt = config_aer_opt, config_aer_ra_feedback = config_aer_ra_feedback,&
              config_chem_opt = config_chem_opt                                            , &
              rthratensw = rthratensw_p , has_reqc   = has_reqc      , has_reqi = has_reqi , &
              has_reqs   = has_reqs     , re_cloud   = recloud_p     , re_ice   = reice_p  , &
              re_snow    = resnow_p     , swupt      = swupt_p       , swuptc   = swuptc_p , &
              swdnt      = swdnt_p      , swdntc     = swdntc_p      , swupb    = swupb_p  , &
              swupbc     = swupbc_p     , swdnb      = swdnb_p       , swdnbc   = swdnbc_p , &
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,        &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,        &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                       )
#endif
       call iamas_timer_stop('RRTMG_sw')

    case ("cam_sw")
       call iamas_timer_start('CAMRAD_sw')
       call camrad( dolw = .false. , dosw = .true. ,                                         &
                p_phy         = pres_hyd_p    , p8w           = pres2_hyd_p   ,              &
                pi_phy        = pi_p          , t_phy         = t_p           ,              &
                z             = zmid_p        , dz8w          = dz_p          ,              &            
                rthratenlw    = rthratenlw_p  , rthratensw    = rthratensw_p  ,              &
                swupt         = swupt_p       , swuptc        = swuptc_p      ,              &
                swdnt         = swdnt_p       , swdntc        = swdntc_p      ,              &
                lwupt         = lwupt_p       , lwuptc        = lwuptc_p      ,              &
                lwdnt         = lwdnt_p       , lwdntc        = lwdntc_p      ,              &
                swupb         = swupb_p       , swupbc        = swupbc_p      ,              &
                swdnb         = swdnb_p       , swdnbc        = swdnbc_p      ,              &
                lwupb         = lwupb_p       , lwupbc        = lwupbc_p      ,              &
                lwdnb         = lwdnb_p       , lwdnbc        = lwdnbc_p      ,              &
                swcf          = swcf_p        , lwcf          = lwcf_p        ,              &
                gsw           = gsw_p         , glw           = glw_p         ,              &
                olr           = olrtoa_p      , cemiss        = cemiss_p      ,              &
                taucldc       = taucldc_p     , taucldi       = taucldi_p     ,              & 
                coszr         = coszr_p       , hrang         = hrang(its:ite)       ,              &!qkyang added hrang   =hrang(its:ite)     ,
                albedo        = sfc_albedo_p  ,              & 
                emiss         = sfc_emiss_p   , tsk           = tsk_p         ,              & 
                xlat          = xlat_p        , xlong         = xlon_p        ,              &
                rho_phy       = rho_p         , qv3d          = qv_p          ,              & 
                qc3d          = qc_p          , qr3d          = qr_p          ,              &
                qi3d          = qi_p          , qs3d          = qs_p          ,              &
                qg3d          = qg_p          , f_qv          = f_qv          ,              &
                f_qc          = f_qc          , f_qr          = f_qr          ,              &
                f_qi          = f_qi          , f_qs          = f_qs          ,              &
                f_qg          = f_qg          , f_ice_phy     = f_ice         ,              &
                f_rain_phy    = f_rain        , cldfra        = cldfrac_p     ,              &
                xland         = xland_p       , xice          = xice_p        ,              &
                num_months    = num_months    , levsiz        = num_oznlevels ,              & 
                pin0          = pin_p         , ozmixm        = ozmixm_p      ,              &
                paerlev       = num_aerlevels , naer_c        = num_aerosols  ,              &
                m_psp         = m_psp_p       , m_psn         = m_psn_p       ,              &
                aerosolcp     = aerosolcp_p   , aerosolcn     = aerosolcn_p   ,              &
                m_hybi0       = m_hybi_p      , snow          = snow_p        ,              &
                cam_abs_dim1  = cam_abs_dim1  , cam_abs_dim2  = cam_abs_dim2  ,              &
                gmt           = gmt           , yr            = year          ,              &
                julday        = julday        , julian        = curr_julday   ,              &
                dt            = dt_dyn        , xtime         = xtime_m       ,              &
                declin        = declin        , solcon        = solcon        ,              &
                radt          = radt          , degrad        = degrad        ,              &
                n_cldadv      = 3             , abstot_3d     = abstot_p      ,              &
                absnxt_3d     = absnxt_p      , emstot_3d     = emstot_p      ,              &
                doabsems      = doabsems      ,                                              &
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,      &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,      &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte        &
                  )
       call iamas_timer_stop('CAMRAD_sw')

    case default

 end select radiation_sw_select

#ifdef DO_CHEMISTRY
! RRTMG leaves aerosol-component diagnostics undefined in a pure-physics run.
! Keep the chemistry-only "other" SW diagnostics finite and identically zero.
 if (trim(config_chem_opt) == "no_chemistry") then
    swdnb_other_p(its:ite,jts:jte) = 0.0_RKIND
    swdnt_other_p(its:ite,jts:jte) = 0.0_RKIND
    swupb_other_p(its:ite,jts:jte) = 0.0_RKIND
    swupt_other_p(its:ite,jts:jte) = 0.0_RKIND
 end if
#endif

!copy local arrays to iAMAS grid:
 call radiation_sw_to_iAMAS(diag_physics,tend_physics,its,ite)

!Added by qkyang to support slope and shading effects>>>>>>>>>>>>>>>>>
if (config_slope_rad == 1 .or. config_topo_shading == 1) then
   call apply_slope_rad_correction(mesh,sfc_input,diag_physics,its,ite, &
                                   config_slope_rad,config_topo_shading,config_shadlen)
endif
!Added by qkyang to support slope and shading effects<<<<<<<<<<<<<<<<

!call iamas_log_write('--- end subroutine driver_radiation_sw.')

 end subroutine driver_radiation_sw

!Added by qkyang to support slope and shading effects>>>>>>>>>>>>>>>>>
!=================================================================================================================
subroutine apply_slope_rad_correction(mesh,sfc_input,diag_physics,its,ite, &
                                      config_slope_rad,config_topo_shading,config_shadlen)
!=================================================================================================================
 ! This routine applies a terrain slope correction to surface shortwave flux.
 ! Design choice:
 ! - Keep the implementation local to the SW driver (minimal module spread).
 ! - Reuse already available diagnostic radiation quantities.
 !
 ! Physical idea:
 ! 1) Read precomputed static terrain geometry (slope/hx/hy and horizon angle).
 ! 2) Convert gradient to slope azimuth when slope correction is enabled.
 ! 3) Recompute effective cosine of solar zenith over the sloped surface.
 ! 4) Scale direct-beam part of shortwave radiation, keep diffuse part less sensitive.
 !
 ! Important scope:
 ! - This is slope correction only (no terrain shadow mask propagation).

 integer,intent(in):: its,ite
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: sfc_input
 type(iamas_pool_type),intent(inout):: diag_physics
integer,intent(in):: config_slope_rad,config_topo_shading
real(kind=RKIND),intent(in):: config_shadlen

 integer:: i,iBin,nShadowBins
 real(kind=RKIND):: pi
 real(kind=RKIND):: hx,hy,slope,slp_azi
 real(kind=RKIND):: csza,sinz,denom,argu,sol_azi,csza_slp
 real(kind=RKIND):: solar_elev,azi_norm
 real(kind=RKIND):: ratio,log_ratio,diffuse_frac,corr_fac

 integer,dimension(:),pointer:: shadowmask
 real(kind=RKIND),dimension(:),pointer:: csza_slp_out
 real(kind=RKIND),dimension(:),pointer:: latCell
 real(kind=RKIND),dimension(:),pointer:: ter_slope,ter_hx,ter_hy
 real(kind=RKIND),dimension(:,:),pointer:: ter_horizon_angle
 real(kind=RKIND),dimension(:),pointer:: coszr,hrang,gsw,swdnb,swdnt
 real(kind=RKIND),pointer:: declin_diag



 ! Mesh latitude is needed by solar-azimuth geometry.
 call iamas_pool_get_array(mesh,'latCell',latCell)
 call iamas_pool_get_array(sfc_input,'ter_slope',ter_slope)
 call iamas_pool_get_array(sfc_input,'ter_hx',ter_hx)
 call iamas_pool_get_array(sfc_input,'ter_hy',ter_hy)
 call iamas_pool_get_array(sfc_input,'ter_horizon_angle',ter_horizon_angle)

 ! Radiation diagnostics already computed by the active SW scheme.
 ! - declin: solar declination
 ! - coszr : cosine solar zenith angle over horizontal surface
 ! - hrang : solar hour angle
 ! - gsw   : net shortwave flux at surface (target to be corrected)
 ! - swdnb/swdnt: used to diagnose diffuse fraction proxy
 ! - csza_slp/csza_slp_out: cosine solar zenith angle over sloped surface
 ! - When the slope faces the sun, csza_slp > 0.
 ! - The more perpendicular sunlight is to the slope, the larger the value,
 ! - with a theoretical maximum of 1 (direct sunlight onto the slope).
 ! - shadowmask: A shadow mask typically refers to the shady side of a mountain.
 call iamas_pool_get_array(diag_physics,'declin',declin_diag)
 call iamas_pool_get_array(diag_physics,'coszr',coszr)
 call iamas_pool_get_array(diag_physics,'hrang',hrang)
 call iamas_pool_get_array(diag_physics,'gsw',gsw)
 call iamas_pool_get_array(diag_physics,'swdnb',swdnb)
 call iamas_pool_get_array(diag_physics,'swdnt',swdnt)
 call iamas_pool_get_array(diag_physics,'shadowmask',shadowmask)
 call iamas_pool_get_array(diag_physics,'csza_slp',csza_slp_out)

 pi = acos(-1._RKIND)

 ! Loop over local solve range for this thread/block slice.
 do i = its,ite
   shadowmask(i) = 0
   csza_slp_out(i) = 0._RKIND

    ! Skip night and numerically tiny flux cases.
    if(coszr(i) <= 1.e-6_RKIND .or. gsw(i) <= 0._RKIND) cycle

    ! Clamp baseline cosine zenith into physical range.
    csza = max(0._RKIND,min(1._RKIND,coszr(i))) !coszr = sin(xxlat) * sin(declin) + cos(xxlat) * cos(declin) * cos(hrang(i))
    ! Very low sun angle is numerically sensitive and physically uncertain.
    if(csza < 1.e-2_RKIND) cycle

    ! Compute solar azimuth from spherical geometry:
    !   csza = cos(zenith)
    ! and use hrang sign to place azimuth in the correct half plane.
    sinz = sqrt(max(0._RKIND,1._RKIND-csza*csza))
    denom = sinz*cos(latCell(i))
    if(abs(denom) > 1.e-12_RKIND) then
       argu = (csza*sin(latCell(i)) - sin(declin_diag))/denom
       ! Numerical safety for acos domain.
       argu = max(-1._RKIND,min(1._RKIND,argu))
       sol_azi = sign(acos(argu),sin(hrang(i))) + pi
    else
       ! Polar and near-singular fallback.
       sol_azi = pi
    endif

   ! Default geometry for "no slope" path.
   slope = 0._RKIND
   slp_azi = 0._RKIND

   if(config_slope_rad == 1) then
      ! Read static terrain geometry precomputed in init_atm_static.
      slope = max(0._RKIND, ter_slope(i))
      hx = ter_hx(i)
      hy = ter_hy(i)
      if(slope > 1.e-4_RKIND) then
         ! Slope azimuth needs to be in the same coordinate system as sol_azi.
         ! sol_azi is defined as: 0=South, pi/2=West, pi=North, -pi/2=East.
         ! atan2(hx,hy) gives the "uphill" direction with 0=North, pi/2=East.
         ! The "downhill" (slope face) direction with 0=South is numerically identical to the "uphill" direction with 0=North.
         ! Therefore, NO +pi is needed here.
         slp_azi = atan2(hx,hy)
      else
         slope = 0._RKIND
      endif
   endif

    ! Effective cosine zenith over sloped terrain.
    ! Equivalent to projecting sun direction onto slope normal.
    csza_slp = csza*cos(slope) + sinz*sin(slope)*cos(sol_azi-slp_azi)
    ! Negative projection means self-shaded slope for direct beam.
    if(csza_slp <= 1.e-4_RKIND) csza_slp = 0._RKIND

   if(config_topo_shading == 1 .and. config_shadlen > 0._RKIND) then
      ! Use static directional horizon angle:
      ! if horizon(dir) >= solar_elevation then direct beam is blocked.
      nShadowBins = size(ter_horizon_angle,1)
      solar_elev = asin(csza)
      
      ! Match solar azimuth to the bin definition in AddSlopeShadowParam.
      ! In AddSlopeShadowParam: 0=North, pi/2=East, pi=South, 3pi/2=West
      ! sol_azi here is defined as: 0=South, pi/2=West, pi=North, 3pi/2=East
      ! Therefore, matched_azi = sol_azi + pi
      azi_norm = modulo(sol_azi + pi, 2._RKIND*pi)
      
      iBin = 1 + int(azi_norm / (2._RKIND*pi) * real(nShadowBins, RKIND))
      iBin = max(1, min(nShadowBins, iBin))
      if(ter_horizon_angle(iBin,i) >= solar_elev) then
         shadowmask(i) = 1
         csza_slp = 0._RKIND
      endif
   endif

    ! Diffuse fraction proxy follows the same WRF slope_rad style:
    ! TOA vs surface downwelling ratio is used as a cloudiness/transmittance indicator.
    ! More diffuse conditions => smaller terrain orientation impact.
    ratio = max(1.0001_RKIND,swdnt(i)/max(swdnb(i),1.e-3_RKIND))
    log_ratio = log(ratio)
    if(log_ratio > 1.e-6_RKIND) then
       diffuse_frac = min(1._RKIND,1._RKIND/max(0.1_RKIND,2.1_RKIND-2.8_RKIND*log(log_ratio)))
    else
       diffuse_frac = 1._RKIND
    endif

    ! Correction factor:
    ! - diffuse component unchanged
   ! - direct component scaled by csza_slp / csza
   ! This keeps behavior consistent with the WRF slope_rad formulation, and
   ! naturally supports topo_shading (csza_slp=0 when blocked).
    corr_fac = 1._RKIND
    if(diffuse_frac < 1._RKIND) then
       corr_fac = diffuse_frac + (1._RKIND-diffuse_frac)*csza_slp/max(csza,1.e-6_RKIND)
       ! Prevent negative surface flux from geometry correction.
       corr_fac = max(0._RKIND,corr_fac)
    endif

    ! Apply correction to net surface shortwave flux.
    gsw(i) = gsw(i)*corr_fac
    csza_slp_out(i) = csza_slp
 end do
end subroutine apply_slope_rad_correction
!Added by qkyang to support slope and shading effects<<<<<<<<<<<<<<<<

!=================================================================================================================
 subroutine radconst(declin,solcon,julian,degrad,dpd)
!=================================================================================================================

!input arguments:
!integer,intent(in):: julian
 real(kind=RKIND),intent(in):: julian
 real(kind=RKIND),intent(in):: degrad,dpd

!output arguments:
 real(kind=RKIND),intent(out):: declin,solcon

!local variables:
 real(kind=RKIND):: obecl,sinob,sxlong,arg,decdeg,djul,rjul,eccfac

!-----------------------------------------------------------------------------------------------------------------

   declin=0.
   solcon=0.

!obecl : obliquity = 23.5 degree.
        
   obecl=23.5*degrad
   sinob=sin(obecl)
        
!calculate longitude of the sun from vernal equinox:        

   if(julian.ge.80.)sxlong=dpd*(julian-80.)
   if(julian.lt.80.)sxlong=dpd*(julian+285.)
   sxlong=sxlong*degrad
   arg=sinob*sin(sxlong)
   declin=asin(arg)
   decdeg=declin/degrad

!solar constant eccentricity factor (paltridge and platt 1976)

   djul=julian*360./365.
   rjul=djul*degrad
   eccfac=1.000110+0.034221*cos(rjul)+0.001280*sin(rjul)+0.000719*  &
          cos(2*rjul)+0.000077*sin(2*rjul)
   solcon=solcon_0*eccfac

 end subroutine radconst

!=================================================================================================================
 end module iamas_atmphys_driver_radiation_sw
!=================================================================================================================
