!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute photolysis rates for chem(:,:,:,:).
! * Author: Zihan Xia from AEMOL USTC, xiazh@mail.ustc.edu.cn
! Last update: zhxia 05/2024
!
! Contact:
! Chun Zhao, PhD
! Professor
! University of Science and Technology of China
! Hefei, Anhui, 230026
! Email: ChunZhao@ustc.edu.cn
!
! References:
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_driver_photolysis
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_atmchem_vars
 use iamas_atmchem_configs
 use iamas_atmchem_utilities
 use iamas_atmphys_vars

 use module_phot_fastj

 implicit none
 private
 public:: driver_photolysis

 !integer,pointer :: nswbands.
 !define in module_data_rrtmgaeropt.F 
 INTEGER nswbands,nlwbands
 PARAMETER (nswbands = 14,nlwbands=16)

CONTAINS

!=================================================================================================================
 subroutine driver_photolysis(curr_secs, diag_chemistry, tend_chemistry, configs, its, ite)
!=================================================================================================================

 IMPLICIT NONE
 real(kind=RKIND),intent(in   ):: curr_secs
 type(iamas_pool_type),intent(in   ):: configs
 type(iamas_pool_type),intent(in   ):: diag_chemistry
 type(iamas_pool_type),intent(inout):: tend_chemistry
 integer,intent(in):: its,ite

! local variables
 character(len=StrKIND),pointer :: phot_opt, chem_opt
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw, gaersw, waersw, bscoefsw
 real(kind=RKIND),dimension(:,:,:),pointer :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer
 real(kind=RKIND),dimension(:,:,:),pointer :: ph

 call iamas_pool_get_config(configs,'config_phot_opt',phot_opt)
 call iamas_pool_get_config(configs,'config_chem_opt',chem_opt)

 if (phot_opt .ne. 'no_phot') then
     !call iamas_pool_get_array(diag_chemistry,'tauaersw'  , tauaersw )
     !call iamas_pool_get_array(diag_chemistry,'gaersw'    , gaersw   )
     !call iamas_pool_get_array(diag_chemistry,'waersw'    , waersw   )
     !call iamas_pool_get_array(diag_chemistry,'bscoefsw'  , bscoefsw )
     call iamas_pool_get_array(diag_chemistry,'l2aer'     , l2aer    )
     call iamas_pool_get_array(diag_chemistry,'l3aer'     , l3aer    )
     call iamas_pool_get_array(diag_chemistry,'l4aer'     , l4aer    )
     call iamas_pool_get_array(diag_chemistry,'l5aer'     , l5aer    )
     call iamas_pool_get_array(diag_chemistry,'l6aer'     , l6aer    )
     call iamas_pool_get_array(diag_chemistry,'l7aer'     , l7aer    )
     call iamas_pool_get_array(tend_chemistry,'ph'        , ph       )

     call chemistry_debug_message('photolysis rate calculation')

     call chemistry_debug_message("config_phot_opt = "//trim(phot_opt))

     phot_opt_select: SELECT CASE(trim(phot_opt))
     CASE ("photfastj")
         if (trim(chem_opt) == 'bin04aer_gas_aq' .OR. &
            trim(chem_opt) == 'bin04aer_gas') then
            call chemistry_debug_message('photolysis driver: call fastj_driver')

            call fastj_driver(                                                                             &
                 id = 1, curr_secs = curr_secs, dtstep = dt_dyn,                                             &
                 configs = configs, gmt = gmt, julday = INT(curr_julday),                                  &
                 t_phy = t_phy(kts:kte,its:ite), moist = moist(1:num_moist,kts:kte,its:ite),               &
                 p8w = pres2_hyd(kts:kte+1,its:ite), p_phy = pres_hyd(kts:kte,its:ite),                    &
                 chem = chem(1:num_chem,kts:kte,its:ite), rho_phy = rho_phy(kts:kte,its:ite),              &
                 dz8w = dz(kts:kte,its:ite), xlat = xlat(its:ite), xlong = xlong(its:ite),                 &
                 z_at_w = zgrid(kts:kte+1,its:ite),                                                        &
                 ph = ph(1:num_phr,kts:kte,its:ite),                                                       &
                 ivgtyp = ivgtyp(its:ite),                                                                 &
                 !tauaersw = tauaersw(1:nswbands,kts:kte,its:ite), gaersw = gaersw(1:nswbands,kts:kte,its:ite),           &
                 !waersw = waersw(1:nswbands,kts:kte,its:ite), bscoefsw = bscoefsw(1:nswbands,kts:kte,its:ite),           &              
                 !don't use in this code.   zhxia,10/2024
                 l2aer = l2aer(1:nswbands,kts:kte,its:ite), l3aer = l3aer(1:nswbands,kts:kte,its:ite),     &
                 l4aer = l4aer(1:nswbands,kts:kte,its:ite), l5aer = l5aer(1:nswbands,kts:kte,its:ite),     &
                 l6aer = l6aer(1:nswbands,kts:kte,its:ite), l7aer = l7aer(1:nswbands,kts:kte,its:ite),     &
                 its=its, ite=ite, kts=kts, kte=kte                                                        &
                 )

            call chemistry_debug_message('fastj success')
         else
            call chemistry_error_fatal('Error: photolysis driver fail. Only chem_opt = bin04aer_gas_aq or bin04aer_gas can use photolysis scheme (fastj). The rest of chem_opt needs to set config_phot_opt = no_phot.')
         endif

     END SELECT phot_opt_select

 else
     call iamas_log_write("Don't do photolysis rate computation")
 end if

end subroutine driver_photolysis

END MODULE iamas_atmchem_driver_photolysis

