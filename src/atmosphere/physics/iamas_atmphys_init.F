! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
!=================================================================================================================
 module iamas_atmphys_init
!qkyang added to support Shin-Hong PBL>>>>>>>
 use iamas_atmphys_vars, only:its,ite,kts,kte
!qkyang added to support Shin-Hong PBL<<<<<<
 use iamas_kind_types
 use iamas_pool_routines
 use iamas_timekeeping

 use iamas_atmphys_driver_convection, only: init_convection
 use iamas_atmphys_driver_lsm,only: init_lsm
 use iamas_atmphys_driver_microphysics
 use iamas_atmphys_driver_radiation_lw, only: init_radiation_lw
 use iamas_atmphys_driver_radiation_sw, only: init_radiation_sw
 use iamas_atmphys_driver_sfclayer

 use iamas_atmphys_landuse
 use iamas_atmphys_o3climatology
 use iamas_init_topo_wind, only: init_topo_wind_fields !qkyang added to support topo_wind=1,2
 use module_tompkins_cloud, only: tompkins_cloud_init

 implicit none
 private
 public:: physics_init

!iAMAS main initialization subroutine for all physics parameterizations.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in iamas_atmphys_init:
! ---------------------------------
! physics_init    : call initialization of individual physics parameterizations.
! init_dir_forphys: needed for initialization of "reconstruct" subroutines.
! r3_normalize    : needed for initialization of "reconstruct" subroutines.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added structure diag in calls to subroutine init_radiation_lw and init_radiation_sw.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-01.
! * added call to subroutine init_o3climatology. reads monthly-mean climatological ozone data and interpolates
!   ozone data to the iAMAS grid.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-03.
! * added the calculation of the mean distance between cell centers.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-22.
! * added initialization of variable xicem.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-24.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * added initialization of the accumulated surface pressure. Added initialization of the tendency and the
!   accumulated tendency of the surface pressure.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * changed the argument list in the call to subroutine microphysics_init, needed to include the Thompson
!   parameterization of cloud microphysics.
!   Laura D. Fowler (laura@ucar.edu) / 2015-03-28.
! * modified the initialization of i_rainc and i_rainnc, now that the convection and cloud microphysics
!   parameterizations are in "packages."
!   Laura D. Fowler (laura@ucar.edu) / 2106-04-13.
! * removed the calculation of the variable dcEdge_m which is no longer needed in the different physics
!   parameterizations.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.


 contains


!=================================================================================================================
 subroutine physics_init(dminfo,clock,configs,mesh,diag,tend,state,time_lev,diag_physics, &
                         atm_input,sfc_input)
!=================================================================================================================
!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>
 character(len=StrKIND),pointer:: pbl_scheme
 real(kind=RKIND),dimension(:,:),pointer::tke_pbl
 integer:: k
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<<
!input arguments:
 type(dm_info),intent(in):: dminfo
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: configs
 type(iAMAS_Clock_type),intent(in):: clock

 integer,intent(in):: time_lev

!inout arguments:
 type(iamas_pool_type),intent(inout):: state
 type(iamas_pool_type),intent(inout):: diag
 type(iamas_pool_type),intent(inout):: tend
 type(iamas_pool_type),intent(inout):: diag_physics
 type(iamas_pool_type),intent(inout):: atm_input
 type(iamas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: config_do_restart,       &
                   config_o3climatology

 character(len=StrKIND),pointer::            &
                   config_radt_cld_scheme,   &
                   config_lsm_scheme,        &
                   config_microp_scheme,     &
                   config_convection_scheme, &
                   config_sfclayer_scheme,   &
                   config_radt_lw_scheme,    &
                   config_radt_sw_scheme

 integer,pointer:: nCellsSolve,nLags
 integer,dimension(:),pointer  :: i_rainc,i_rainnc
 integer,dimension(:),pointer  :: i_acswdnb,i_acswdnbc,i_acswdnt,i_acswdntc, &
                                  i_acswupb,i_acswupbc,i_acswupt,i_acswuptc, &
                                  i_aclwdnb,i_aclwdnbc,i_aclwdnt,i_aclwdntc, &
                                  i_aclwupb,i_aclwupbc,i_aclwupt,i_aclwuptc

 real(kind=RKIND),dimension(:),pointer  :: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                           acswupb,acswupbc,acswupt,acswuptc,  &
                                           aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                           aclwupb,aclwupbc,aclwupt,aclwuptc
 real(kind=RKIND),dimension(:),pointer  :: nsteps_accum,ndays_accum,tday_accum, &
                                           tyear_accum,tyear_mean
 real(kind=RKIND),dimension(:),pointer  :: sst,sstsk,tmn,xice,xicem
 real(kind=RKIND),dimension(:,:),pointer:: tlag
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer:: q_tom, bdw_tom

 real(kind=RKIND),dimension(:),pointer  :: t_oml, t_oml_initial, t_oml_200m_initial
 real(kind=RKIND),dimension(:),pointer  :: h_oml, h_oml_initial, hu_oml, hv_oml
 real(kind=RKIND), pointer :: config_oml_hml0
 integer,pointer:: nCells,index_qv
 logical,pointer:: config_oml1d
 
 

!local variables and arrays:
 type(iAMAS_Time_Type):: currTime

 logical:: init_done
 integer:: ierr,julday 
 integer:: iCell,iLag,iEdge,nEdges_m

!-----------------------------------------------------------------------------------------------------------------
! call iamas_log_write('')
! call iamas_log_write('--- enter subroutine physics_init:')

!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>
 call iamas_pool_get_config(configs,'config_pbl_scheme',pbl_scheme       )
 call iamas_pool_get_array(diag_physics,'tke_pbl',tke_pbl  )
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<<

 call iamas_pool_get_config(configs,'config_do_restart'       ,config_do_restart       )
 call iamas_pool_get_config(configs,'config_o3climatology'    ,config_o3climatology    )
 call iamas_pool_get_config(configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call iamas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )
 call iamas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
 call iamas_pool_get_config(configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )
 call iamas_pool_get_config(configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call iamas_pool_get_config(configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )

 call iamas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call iamas_pool_get_dimension(mesh,'nLags'      ,nLags      )

 call iamas_pool_get_array(diag_physics,'i_acswdnb'   ,i_acswdnb   )
 call iamas_pool_get_array(diag_physics,'i_acswdnbc'  ,i_acswdnbc  )
 call iamas_pool_get_array(diag_physics,'i_acswdnt'   ,i_acswdnt   )
 call iamas_pool_get_array(diag_physics,'i_acswdntc'  ,i_acswdntc  )
 call iamas_pool_get_array(diag_physics,'i_acswupb'   ,i_acswupb   )
 call iamas_pool_get_array(diag_physics,'i_acswupbc'  ,i_acswupbc  )
 call iamas_pool_get_array(diag_physics,'i_acswupt'   ,i_acswupt   )
 call iamas_pool_get_array(diag_physics,'i_acswuptc'  ,i_acswuptc  )
 call iamas_pool_get_array(diag_physics,'i_aclwdnb'   ,i_aclwdnb   )
 call iamas_pool_get_array(diag_physics,'i_aclwdnbc'  ,i_aclwdnbc  )
 call iamas_pool_get_array(diag_physics,'i_aclwdnt'   ,i_aclwdnt   )
 call iamas_pool_get_array(diag_physics,'i_aclwdntc'  ,i_aclwdntc  )
 call iamas_pool_get_array(diag_physics,'i_aclwupb'   ,i_aclwupb   )
 call iamas_pool_get_array(diag_physics,'i_aclwupbc'  ,i_aclwupbc  )
 call iamas_pool_get_array(diag_physics,'i_aclwupt'   ,i_aclwupt   )
 call iamas_pool_get_array(diag_physics,'i_aclwuptc'  ,i_aclwuptc  )

 call iamas_pool_get_array(diag_physics,'acswdnb'     ,acswdnb     )
 call iamas_pool_get_array(diag_physics,'acswdnbc'    ,acswdnbc    )
 call iamas_pool_get_array(diag_physics,'acswdnt'     ,acswdnt     )
 call iamas_pool_get_array(diag_physics,'acswdntc'    ,acswdntc    )
 call iamas_pool_get_array(diag_physics,'acswupb'     ,acswupb     )
 call iamas_pool_get_array(diag_physics,'acswupbc'    ,acswupbc    )
 call iamas_pool_get_array(diag_physics,'acswupt'     ,acswupt     )
 call iamas_pool_get_array(diag_physics,'acswuptc'    ,acswuptc    )
 call iamas_pool_get_array(diag_physics,'aclwdnb'     ,aclwdnb     )
 call iamas_pool_get_array(diag_physics,'aclwdnbc'    ,aclwdnbc    )
 call iamas_pool_get_array(diag_physics,'aclwdnt'     ,aclwdnt     )
 call iamas_pool_get_array(diag_physics,'aclwdntc'    ,aclwdntc    )
 call iamas_pool_get_array(diag_physics,'aclwupb'     ,aclwupb     )
 call iamas_pool_get_array(diag_physics,'aclwupbc'    ,aclwupbc    )
 call iamas_pool_get_array(diag_physics,'aclwupt'     ,aclwupt     )
 call iamas_pool_get_array(diag_physics,'aclwuptc'    ,aclwuptc    )

 call iamas_pool_get_array(diag_physics,'nsteps_accum',nsteps_accum)
 call iamas_pool_get_array(diag_physics,'ndays_accum' ,ndays_accum )
 call iamas_pool_get_array(diag_physics,'tday_accum'  ,tday_accum  )
 call iamas_pool_get_array(diag_physics,'tyear_accum' ,tyear_accum )
 call iamas_pool_get_array(diag_physics,'tyear_mean'  ,tyear_mean  )
 call iamas_pool_get_array(diag_physics,'tlag'        ,tlag        )
 call iamas_pool_get_array(diag_physics,'sstsk'       ,sstsk       )
 call iamas_pool_get_array(diag_physics,'xicem'       ,xicem       )

 call iamas_pool_get_array(sfc_input,'sst' ,sst )
 call iamas_pool_get_array(sfc_input,'tmn' ,tmn )
 call iamas_pool_get_array(sfc_input,'xice',xice)

 call iamas_pool_get_array(diag_physics,'t_oml'             ,t_oml)
 call iamas_pool_get_array(diag_physics,'t_oml_initial'     ,t_oml_initial)
 call iamas_pool_get_array(diag_physics,'t_oml_200m_initial',t_oml_200m_initial)
 call iamas_pool_get_array(diag_physics,'h_oml'             ,h_oml)
 call iamas_pool_get_array(diag_physics,'h_oml_initial'     ,h_oml_initial)
 call iamas_pool_get_array(diag_physics,'hu_oml'            ,hu_oml)
 call iamas_pool_get_array(diag_physics,'hv_oml'            ,hv_oml)
 call iamas_pool_get_config(configs,'config_oml1d'          ,config_oml1d  )
 call iamas_pool_get_config(configs,'config_oml_hml0'       ,config_oml_hml0  )
 call iamas_pool_get_dimension(mesh,'nCells',nCells)

 call iamas_pool_get_config(configs,'config_radt_cld_scheme'   ,config_radt_cld_scheme   )
 currTime = iamas_get_clock_time(clock,iAMAS_NOW,ierr)
 call iamas_get_time(curr_time=currTime,DoY=julday,ierr=ierr)

!initialization of east-north directions to convert u-tendencies from cell centers to cell
!edges:
 call init_dirs_forphys(mesh)
 if(.not. config_do_restart .and. trim(config_radt_cld_scheme) .eq. 'cld_fraction_tompkins') then
    call iamas_pool_get_array(state,'scalars',scalars,time_lev)
    call iamas_pool_get_array(diag,'q_tom', q_tom)
    call iamas_pool_get_array(diag,'bdw_tom', bdw_tom)
    call iamas_pool_get_dimension(state,'index_qv',index_qv)
    call tompkins_cloud_init(ite, kte, scalars(index_qv,kts:kte,its:ite), q_tom(kts:kte,its:ite), bdw_tom(kts:kte,its:ite))
 endif
!qkyang added to support Shin-Hong PBL>>>>>>>>>>>>
!initialization of shortwave radiation processes:
if(.not. config_do_restart) then
    if(trim(pbl_scheme) .eq. 'bl_shinhong') then
    do k = kts,kte
    do iCell = 1, nCellsSolve
       tke_pbl(k,iCell)  = 0.005_RKIND!Based on tke_pbl(i,k,j) = epsq2l/2.,epsq2l = 0.01 in module_bl_shinhong.F
    enddo
    enddo
    !test to output tke_pbl(1,1),tke_pbl(2,1),tke_pbl(2,2),tke_pbl(2,3)
   !  write(*,*) 'tke_pbl(1,1),tke_pbl(2,1),tke_pbl(2,2),tke_pbl(2,3) = ',are all 0.005
    endif
 endif
!qkyang added to support Shin-Hong PBL<<<<<<<<<<<<<

!initialization of counters i_rainc and i_rainnc. i_rainc and i_rainnc track the number of
!times the accumulated convective (rainc) and grid-scale (rainnc) rain exceed the prescribed
!threshold value:
 if(.not. config_do_restart .and. config_convection_scheme.ne.'off') then
    call iamas_pool_get_array(diag_physics,'i_rainc',i_rainc)
    do iCell = 1, nCellsSolve
       i_rainc(iCell)  = 0
    enddo
 endif
 if(.not. config_do_restart .and. config_microp_scheme.ne.'off') then
    call iamas_pool_get_array(diag_physics,'i_rainnc',i_rainnc)
    do iCell = 1, nCellsSolve
       i_rainnc(iCell) = 0
    enddo
 endif

!initialization of counters i_acsw* and i_aclw*. i_acsw* and i_aclw* track the number of times
!the accumulated long and short-wave radiation fluxes exceed their prescribed theshold values.
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       i_acswdnb(iCell)  = 0
       i_acswdnbc(iCell) = 0
       i_acswdnt(iCell)  = 0
       i_acswdntc(iCell) = 0
       i_acswupb(iCell)  = 0
       i_acswupbc(iCell) = 0
       i_acswupt(iCell)  = 0
       i_acswuptc(iCell) = 0

       i_aclwdnb(iCell)  = 0
       i_aclwdnbc(iCell) = 0
       i_aclwdnt(iCell)  = 0
       i_aclwdntc(iCell) = 0
       i_aclwupb(iCell)  = 0
       i_aclwupbc(iCell) = 0
       i_aclwupt(iCell)  = 0
       i_aclwuptc(iCell) = 0

       acswdnb(iCell)  = 0._RKIND
       acswdnbc(iCell) = 0._RKIND
       acswdnt(iCell)  = 0._RKIND
       acswdntc(iCell) = 0._RKIND
       acswupb(iCell)  = 0._RKIND
       acswupbc(iCell) = 0._RKIND
       acswupt(iCell)  = 0._RKIND
       acswuptc(iCell) = 0._RKIND

       aclwdnb(iCell)  = 0._RKIND
       aclwdnbc(iCell) = 0._RKIND
       aclwdnt(iCell)  = 0._RKIND
       aclwdntc(iCell) = 0._RKIND
       aclwupb(iCell)  = 0._RKIND
       aclwupbc(iCell) = 0._RKIND
       aclwupt(iCell)  = 0._RKIND
       aclwuptc(iCell) = 0._RKIND
    enddo
 endif

!initialization of xicem:
 if(.not.config_do_restart) then
!    call iamas_log_write('--- initialization of xicem:')
    do iCell = 1, nCellsSolve
       xicem(iCell) = xice(iCell)
    enddo
 endif

!initialization of the local sea-surface temperature when a diurnal cycle of the
!sea-surface temperature is applied. This avoids having the array sstsk equal to
!zero over land:
 if(.not. config_do_restart) then
!    call iamas_log_write('--- initialization of sstsk:')
    do iCell = 1, nCellsSolve
       sstsk(iCell) = sst(iCell)
    enddo
 endif

! initialized the 1D ocean mixed-layer model  (code from wrf module_sf_oml)
 if (config_oml1d) then
   if (.not. config_do_restart) then
      call iamas_log_write('--- initialization of 1D ocean mixed layer model ')
      do iCell = 1, nCellsSolve
        t_oml(iCell) = sst(iCell)
        t_oml_initial(iCell) = sst(iCell)
      end do
      if (config_oml_hml0 .gt. 0) then
        do iCell = 1, nCellsSolve
          h_oml(iCell) = config_oml_hml0
          h_oml_initial(iCell) = config_oml_hml0
          hu_oml(iCell) = 0.
          hv_oml(iCell) = 0.
          t_oml_200m_initial(iCell) = sst(iCell) - 5.
        end do
      else if (config_oml_hml0 .eq. 0) then
! initializing with climatological mixed layer depth only
        do iCell = 1, nCellsSolve
          h_oml(iCell) = h_oml_initial(iCell)
          hu_oml(iCell) = 0.
          hv_oml(iCell) = 0.
          t_oml_200m_initial(iCell) = sst(iCell) - 5.
        end do
      else
        do iCell = 1, nCellsSolve
          h_oml(iCell) = h_oml_initial(iCell)
          ! WRF COMMENT:
          ! fill in near coast area with SST: 200 K was set as missing value in ocean pre-processing code
          if( (t_oml_200m_initial(iCell) > 200.) .and. (t_oml_200m_initial(iCell) <= 200.) )  &
               t_oml_200m_initial(iCell) = sst(iCell)
        end do
      end if
   end if
 end if

!initialization of temperatures needed for updating the deep soil temperature:
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       nsteps_accum(iCell) = 0._RKIND
       ndays_accum(iCell)  = 0._RKIND
       tday_accum(iCell)   = 0._RKIND
       tyear_accum(iCell)  = 0._RKIND
       tyear_mean(iCell)   = tmn(iCell)
       do iLag = 1, nLags
          tlag(iLag,iCell) = tmn(iCell)
       enddo
    enddo
 endif

!read the input files that contain the monthly-mean ozone climatology on fixed pressure levels:
 if(config_o3climatology .and. (.not. config_do_restart)) &
    call init_o3climatology(mesh,atm_input)

!initialization of global surface properties. set here for now, but may be moved when time
!manager is implemented:
 call landuse_init_foriAMAS(dminfo,julday,mesh,configs,diag_physics,sfc_input)

!initialization of parameterized convective processes:
 if(config_convection_scheme .ne. 'off') &
    call init_convection(mesh,configs,diag_physics)

!initialization of cloud microphysics processes:
 if(config_microp_scheme .ne. 'off') &
    call microphysics_init(dminfo,configs,mesh,sfc_input,diag_physics)

!initialization of surface layer processes:
 if(config_sfclayer_scheme .ne. 'off') call init_sfclayer(configs)

!initialization of land-surface model:
!if(.not. config_do_restart) then
!   if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,configs,diag_physics,sfc_input)
!endif
 if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,configs,diag_physics,sfc_input)


!qkyang added to support topo_wind=1,2>>>>>>>>>>>>>>
!initialization of topographic wind parameters: 
   call init_topo_wind_fields(configs, mesh, sfc_input, diag_physics)
!qkyang added end to support topo_wind=1,2<<<<<<<<<<<

!initialization of shortwave radiation processes:
 init_done = .false.
 if(config_radt_sw_scheme.ne.'off') then
    if(trim(config_radt_sw_scheme) .eq. 'cam_sw') then
       call init_radiation_sw(dminfo,configs,mesh,atm_input,diag,diag_physics,state,time_lev)
       init_done = .true.
    else
       call init_radiation_sw(dminfo,configs)
    endif
 endif

!initialization of longwave radiation processes: if we run the CAM radiation codes, the initia
!lization of the longwave and shortwave parameterizations is the same, and needs to be called
!only once:
 if(config_radt_lw_scheme.ne.'off') then    
    if(trim(config_radt_lw_scheme) .eq. 'cam_lw') then
       if(.not. init_done) then
          call init_radiation_lw(dminfo,configs,mesh,atm_input,diag,diag_physics,state,time_lev)
       else
!          call iamas_log_write('')
!          call iamas_log_write('--- camrad lw initialization done above')
       endif
    else
       call init_radiation_lw(dminfo,configs)
    endif
 endif
 
! call iamas_log_write('')
! call iamas_log_write('--- end subroutine physics_init')
! call iamas_log_write('')

 end subroutine physics_init

!=================================================================================================================
 subroutine init_dirs_forphys(mesh)
!=================================================================================================================

!inout arguments:
!----------------
 type(iamas_pool_type),intent(in):: mesh

!local pointers:
 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------

 call iamas_pool_get_dimension(mesh,'nCells',nCells)

 call iamas_pool_get_array(mesh,'latCell',latCell)
 call iamas_pool_get_array(mesh,'lonCell',lonCell)
 call iamas_pool_get_array(mesh,'east'   ,east   )
 call iamas_pool_get_array(mesh,'north'  ,north  )

!Compute unit vectors in east and north directions for each cell:
 do iCell = 1, nCells

    east(1,iCell) = -sin(lonCell(iCell))
    east(2,iCell) =  cos(lonCell(iCell))
    east(3,iCell) =  0.0
    call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

    north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
    north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
    north(3,iCell) =  cos(latCell(iCell))
    call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

 end do

 end subroutine init_dirs_forphys

!=================================================================================================================
 subroutine r3_normalize(ax, ay, az)
!=================================================================================================================
!normalizes the vector (ax, ay, az)

 real (kind=RKIND), intent(inout) :: ax, ay, az
 real (kind=RKIND) :: mi

!-----------------------------------------------------------------------------------------------------------------

 mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

 end subroutine r3_normalize

!=================================================================================================================
 end module iamas_atmphys_init
!=================================================================================================================
