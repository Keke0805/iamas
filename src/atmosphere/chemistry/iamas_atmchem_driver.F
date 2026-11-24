!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! General driver for all chemistry procecces.
! * Author: Jiawang Feng and Zihan Xia from AEMOL USTC
! Last update: zhxia 9/2024
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! References: 
!     The code structure is similar to the physics drivers
!     in iAMAS-v5.2:src/core_atmosphere/physics/iamas_atmphys_driver.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_driver
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines

 use iamas_atmchem_driver_optical
 use iamas_atmchem_driver_photolysis 
 use iamas_atmchem_driver_emission  
 use iamas_atmchem_driver_drydep
 use iamas_atmchem_driver_convection
 use iamas_atmchem_driver_wetscav
 use iamas_atmchem_driver_kinetics_thermodynamics

 use iamas_atmchem_interface
 use iamas_atmchem_vars
 use iamas_atmchem_binaer_init
 use iamas_atmchem_utilities
 use iamas_timer

 implicit none
 private
 public:: chemistry_driver,driver_chem_diag


!iAMAS top chemistry driver.

! subroutine chemistry_driver is the top chemistry driver from which separate drivers.

! <Description> zhxia 09/2024
! The main sections of chem_driver are:
!
! (1) Initialization of variables as needed for chemistry
!
! (2) Calls to drivers for the various chemistry tasks:
!        optical_driver
!        photolysis_driver
!        emissions_driver
!        dry_dep_driver
!        convection_driver (convective tracer transport)
!        kinetics_thermodynamics_driver (gas & aerosol)
!        wetscav_driver
 
 contains

!=================================================================================================================
 subroutine driver_chem_diag(domain)
!=================================================================================================================
 use module_data_binaer_asect

!inout arguments:
 type(domain_type),intent(inout):: domain

!local arguments:
 type(iamas_pool_type),pointer :: mesh
 type(iamas_pool_type),pointer :: state
 type(iamas_pool_type),pointer :: diag_chemistry
 type(block_type),pointer:: block

 integer:: thread
 integer,pointer:: nThreads, nVertLevels
 integer,dimension(:),pointer:: cellSolveThreadStart, cellSolveThreadEnd
 integer :: its,ite

 integer:: i,k,n,time_lev
 integer:: l_dust,l_cl,l_na
 integer:: iphase,itype
 real(kind=RKIND),dimension(:,:),pointer :: pm2_5_dust, pm2_5_seas
 real(kind=RKIND),dimension(:,:),pointer :: pm10_dust, pm10_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: chem

 itype = 1
 iphase = ai_phase

block => domain % blocklist
do while(associated(block))
    call iamas_pool_get_subpool(block%structs,'mesh'          ,mesh          )
    call iamas_pool_get_subpool(block%structs,'state'         ,state         )
    call iamas_pool_get_subpool(block%structs,'diag_chemistry',diag_chemistry)

    call iamas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)
    call iamas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
    call iamas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

    call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
    call iamas_pool_get_array(diag_chemistry,'pm2_5_dust',pm2_5_dust)
    call iamas_pool_get_array(diag_chemistry,'pm2_5_seas',pm2_5_seas)
    call iamas_pool_get_array(diag_chemistry,'pm10_dust',pm10_dust)
    call iamas_pool_get_array(diag_chemistry,'pm10_seas',pm10_seas)
    time_lev = 2
    call iamas_pool_get_array(state,'chem',chem,time_lev)
    
!jwfeng: Please CHECK if I have used OMP right
!$OMP PARALLEL DO
    do thread=1,nThreads
        its = cellSolveThreadStart(thread)
        ite = cellSolveThreadEnd(thread)
        do k = 1, nVertLevels
        do i = its,ite
           pm2_5_dust(k,i) = 0.0_RKIND
           pm2_5_seas(k,i) = 0.0_RKIND
           pm10_dust(k,i) = 0.0_RKIND
           pm10_seas(k,i) = 0.0_RKIND
        enddo
        enddo

        do k = 1, nVertLevels
        do i = its,ite
        do n = 1, nsize_aer(itype)
           l_dust = lptr_dust_aer(n,itype,iphase)
           if (l_dust > param_first_scalar .and. dhi_sect(n,itype) < 10.001e-4 ) then
               pm10_dust(k,i) = chem(l_dust,k,i) + pm10_dust(k,i)
           end if
           if (l_dust > param_first_scalar .and. dhi_sect(n,itype) < 2.5001e-4 ) then
               pm2_5_dust(k,i) = chem(l_dust,k,i) + pm2_5_dust(k,i)
           end if

           l_cl = lptr_cl_aer(n,itype,iphase)
           l_na = lptr_na_aer(n,itype,iphase)
           if (l_cl > param_first_scalar .and. dhi_sect(n,itype) < 10.001e-4 ) then
               pm10_seas(k,i) = chem(l_cl,k,i) + pm10_seas(k,i)
           end if
           if (l_na > param_first_scalar .and. dhi_sect(n,itype) < 10.001e-4 ) then
               pm10_seas(k,i) = chem(l_na,k,i) + pm10_seas(k,i)
           end if
           if (l_cl > param_first_scalar .and. dhi_sect(n,itype) < 2.5001e-4 ) then
               pm2_5_seas(k,i) = chem(l_cl,k,i) + pm2_5_seas(k,i)
           end if
           if (l_na > param_first_scalar .and. dhi_sect(n,itype) < 2.5001e-4 ) then
               pm2_5_seas(k,i) = chem(l_na,k,i) + pm2_5_seas(k,i)
           end if
        enddo
        enddo
        enddo
    end do
!$OMP END PARALLEL DO

    block => block % next
end do

 end subroutine driver_chem_diag

!=================================================================================================================
 subroutine chemistry_driver(domain,itimestep,xtime_d,xtime_s,curr_secs)
!=================================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s,xtime_d,curr_secs !zhxia add curr_secs for total simulation secons from start time.

!inout arguments:
 type(domain_type),intent(inout):: domain

!local pointers:
 type(iamas_pool_type),pointer::  configs,        &
                                 mesh,           &
                                 state,          &
                                 emis_state,     &
                                 sfc_input,      &
                                 diag,           &
                                 diag_physics,   &
                                 diag_chemistry, &
                                 tend,           &
                                 tend_chemistry

!local variables:
 type(block_type),pointer:: block

 integer:: time_lev
 integer:: thread

 integer,pointer:: nThreads
 integer,dimension(:),pointer:: cellSolveThreadStart, cellSolveThreadEnd
 logical,pointer:: config_drydep_on, config_wetscav_on
 real,pointer:: config_photdt, config_dt, config_chemdt !zhxia add chemdt for chemistry.
 logical :: do_photstep, do_chemstep !zhxia add do_chemstep for chemistry.
 logical,pointer:: config_emis_update !zhxia add config_emis_update for chemistry
 character(len=StrKIND),pointer:: config_dust_opt, config_aer_op_opt, config_ctrans_opt, config_phot_opt, config_gaschem_opt, config_stra_gaschem_opt !zhxia add phot_opt, gaschem_opt, and config_stra_gaschem_opt for chemistry.
 character(len=StrKIND),pointer:: config_aerchem_sia_opt, config_aerchem_soa_opt !xyruan add for aerosol chemistry
 character(len=StrKIND):: debugString

! zhxia add these for updating the current seconds, current julian day and current year
 type (iAMAS_Clock_type), pointer :: clock
 type(iAMAS_Time_Type) :: currTime
 character(len=StrKIND) :: timeStamp
 integer :: ierr
 real(kind=RKIND) :: utc_s,utc_h
 integer :: h,m,s,s_n,s_d,DoY,yr,month

! zhxia add these for chemical computation
 real(kind=RKIND) :: dtstepc
 real(kind=RKIND), save :: last_chem_time
 integer, save :: itimestepc =0
!=================================================================================================================

call iamas_timer_start('chemistry driver')
call chemistry_message('----------------------------------------------------------------------')
call chemistry_message('----------------------Starting chemistry driver-----------------------')
call chemistry_message('----------------------------------------------------------------------')

! zhxia add: Update the current seconds, current julian day and current year
  clock => domain % clock
  currTime    = iamas_get_clock_time(clock,iAMAS_NOW,ierr)
  call iamas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,MM=month,H=h,M=m,&
                     S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
  utc_h       = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
  utc_s       = real(s + s_n/s_d)
  curr_year   = yr
  curr_month  = month
  curr_julday = real(DoY-1) + utc_h / 24.0 !start from 0
  LeapYear    = isLeapYear(yr)


block => domain % blocklist
do while(associated(block))

    call iamas_pool_get_subpool(block%structs,'mesh'          ,mesh          )
    call iamas_pool_get_subpool(block%structs,'state'         ,state         )
    call iamas_pool_get_subpool(block%structs,'emis_state'    ,emis_state    )
    call iamas_pool_get_subpool(block%structs,'diag'          ,diag          )
    call iamas_pool_get_subpool(block%structs,'diag_physics'  ,diag_physics  )
    call iamas_pool_get_subpool(block%structs,'diag_chemistry',diag_chemistry)
    call iamas_pool_get_subpool(block%structs,'sfc_input'     ,sfc_input     )
    call iamas_pool_get_subpool(block%structs,'tend'          ,tend          )
    call iamas_pool_get_subpool(block%structs,'tend_chemistry',tend_chemistry)

    call iamas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

    call iamas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
    call iamas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

!	jwfeng from USTC AEMOL:  02/2020
!	Preparation for chemistry_driver

!	jwfeng:  03/2020
    call iamas_pool_get_config(block%configs,'config_aer_op_opt' ,config_aer_op_opt )
    call iamas_pool_get_config(block%configs,'config_phot_opt' ,config_phot_opt) !zhxia add this for photolysis
    call iamas_pool_get_config(block%configs,'config_gaschem_opt' ,config_gaschem_opt) !zhxia add this for photochemistry
    call iamas_pool_get_config(block%configs,'config_aerchem_sia_opt' ,config_aerchem_sia_opt) !xyruan add this for inorganic aerosol chemistry and interactions
    call iamas_pool_get_config(block%configs,'config_aerchem_soa_opt' ,config_aerchem_soa_opt) !xyruan add this for organic aerosol chemistry and interactions
    call iamas_pool_get_config(block%configs,'config_stra_gaschem_opt' ,config_stra_gaschem_opt) !zhxia add this for stratospheric photochemistry
    call iamas_pool_get_config(block%configs,'config_drydep_on'  ,config_drydep_on  )
    call iamas_pool_get_config(block%configs,'config_ctrans_opt' ,config_ctrans_opt )
    call iamas_pool_get_config(block%configs,'config_wetscav_on' ,config_wetscav_on )
    call iamas_pool_get_config(block%configs,'config_dust_opt'   ,config_dust_opt   )
    call iamas_pool_get_config(block%configs,'config_emis_update',config_emis_update) !zhxia add this for emission_update.
    call iamas_pool_get_config(block%configs,'config_photdt'     ,config_photdt     )
    call iamas_pool_get_config(block%configs,'config_dt'         ,config_dt         )
    call iamas_pool_get_config(block%configs,'config_chemdt'     ,config_chemdt     ) !zhxia add this for chemistry integration
!-----------------------------------------------------------------------------------------!
! Allocate arrays shared by all chemistry parameterizations:
    call allocate_forall_chemistry
    call chemistry_debug_message('forall chemistry allocate success')
    time_lev = 1
    call iamas_timer_start("chemistry_from_iAMAS")
!$OMP PARALLEL DO
    do thread=1,nThreads
       call chemistry_from_iAMAS(block%configs,mesh,state,time_lev,sfc_input,diag,diag_physics, tend_chemistry, &
                                cellSolveThreadStart(thread), cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO
    call iamas_timer_stop("chemistry_from_iAMAS")
!-----------------------------------------------------------------------------------------!
! Optical start 
!-----------------------------------------------------------------------------------------!
    do_photstep = .false.
    if(itimestep==1) then  !First step always do photstep
       do_photstep = .true.
    else
       if ( (config_photdt<=0) .or. &
            ( &
              ( curr_secs+real(config_dt,8)+0.01 ) >= &
              ( int( curr_secs/real(config_photdt,8)+1,8 )*real(config_photdt,8) ) &
            ) &
          ) then
          do_photstep = .true.
       endif
    endif
    if(do_photstep) then
        call chemistry_message('--- time to run the dust optical driver config_aer_op_opt   = '//trim(config_aer_op_opt))
        call iamas_timer_start('optical driver')

!$OMP PARALLEL DO
        do thread=1,nThreads
           call driver_optical(curr_secs, diag, diag_chemistry, block%configs, &
                               cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
        end do
!$OMP END PARALLEL DO
        call chemistry_debug_message('optical driver success')
! Optical end 
        call iamas_timer_stop('optical driver')
    else
        call chemistry_message('--- time to run the dust optical driver config_aer_op_opt   = Dont do photstep at this timestep')
    endif
!-----------------------------------------------------------------------------------------!
! Photolysis start. (zhxia 09/2024)
!-----------------------------------------------------------------------------------------!
    if(do_photstep) then
        call chemistry_message('--- time to run the photolysis driver config_phot_opt = '//trim(config_phot_opt))
        call iamas_timer_start('photolysis driver')
        
!$OMP PARALLEL DO
        do thread=1,nThreads
           call driver_photolysis(curr_secs, diag_chemistry, tend_chemistry, block%configs, &
                                  cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
        end do
!$OMP END PARALLEL DO
        call chemistry_debug_message('photolysis driver success')
! Photolysis end
        call iamas_timer_stop('photolysis driver')
    else
        call chemistry_message('--- time to run the photolysis driver config_phot_opt = Dont do photstep at this timestep')
    endif
!-----------------------------------------------------------------------------------------!
! Emisssion start (zhxia 04/2025)
!-----------------------------------------------------------------------------------------!
    call chemistry_message('--- time to run the emission driver')
    call chemistry_message('config_dust_opt = '//trim(config_dust_opt)//'  config_emis_update = '//merge("true ", "false", config_emis_update))
    call iamas_timer_start('emission driver')

!$OMP PARALLEL DO
    do thread=1,nThreads
       call driver_emission(itimestep,block%configs,mesh,state,emis_state,sfc_input , &
                            diag,diag_physics,diag_chemistry,tend_chemistry         , &
                            cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO

    call chemistry_debug_message('emission driver success')

! Emisssion end 
    call iamas_timer_stop('emission driver')
!-----------------------------------------------------------------------------------------!
! Dry deposition start
!-----------------------------------------------------------------------------------------!
    call iamas_timer_start('dry deposition')
    write(debugString,'(a,l)') '--- time to run the dry deposition driver config_drydep_on  = ', config_drydep_on
    call chemistry_message(trim(debugString))
if (config_drydep_on) then

!$OMP PARALLEL DO
    do thread=1,nThreads
       call driver_drydep(itimestep,block%configs,mesh,state,                    &
                          sfc_input,diag,diag_physics,tend,tend_chemistry,       &
                          cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO
    call chemistry_debug_message('drydep driver success')
! Dry deposition end 
end if
    call iamas_timer_stop('dry deposition')
!-----------------------------------------------------------------------------------------!
! chem convection start
!-----------------------------------------------------------------------------------------!
    call iamas_timer_start('chem convection')
if (trim(config_ctrans_opt) .ne. 'no_ctrans') then
    call chemistry_message('--- time to run the chem ctrans driver config_ctrans_opt    = '//trim(config_ctrans_opt))
    
!$OMP PARALLEL DO
    do thread=1,nThreads
       call driver_convection( &
            state, diag_physics, time_lev, &
            block%configs, diag_chemistry, tend_chemistry, &
            cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO
    call chemistry_debug_message('ctrans driver success')
! Chem convection end
end if
    call iamas_timer_stop('chem convection')
!-----------------------------------------------------------------------------------------!
! Kinetics thermodynamics start (zhxia 09/2024) 
! History:
! (1)zhxia 09/2024: Add photochemistry driver with AIPC scheme.
! (2)zhxia 10/2024: Add aerosol chemistry driver with AIMACI scheme.
! (3)zhxia 09/2025: Combine gas and aerosol into one driver called Kinetics thermodynamics driver.
!-----------------------------------------------------------------------------------------!
    do_chemstep = .false.
    if (itimestep==1) then
       do_chemstep = .true.
       dtstepc = config_dt
       last_chem_time = curr_secs
       itimestepc = itimestepc + 1
    else
       if ( (config_chemdt<=0) .or. &
            ( &
              ( curr_secs+real(config_dt,8)+0.01 ) >= &
              ( int( curr_secs/real(config_chemdt,8)+1,8 )*real(config_chemdt,8) ) &
            ) &
          ) then  
          do_chemstep = .true.
          dtstepc = curr_secs - last_chem_time
          last_chem_time = curr_secs
          itimestepc = itimestepc + 1
       endif
    endif

    if (do_chemstep) then
        call chemistry_message('--- time to run the kinetics thermodynamics driver: config_gaschem_opt = '//trim(config_gaschem_opt))
        call chemistry_message('--- time to run the kinetics thermodynamics driver: config_aerchem_sia_opt =' //trim(config_aerchem_sia_opt))
        call chemistry_message('--- time to run the kinetics thermodynamics driver: config_aerchem_soa_opt =' //trim(config_aerchem_soa_opt))
        call chemistry_message('--- time to run the kinetics thermodynamics driver: config_stra_gaschem_opt = '//trim(config_stra_gaschem_opt))

        call iamas_timer_start('Kinetics thermodynamics driver')

!$OMP PARALLEL DO
        do thread = 1, nThreads
            call driver_kinetics_thermodynamics(curr_secs, itimestep, itimestepc, dtstepc, state, time_lev,  &
                           diag, diag_physics, tend_chemistry, mesh, block%configs, &
                           cellSolveThreadStart(thread), cellSolveThreadEnd(thread))
        end do
!$OMP END PARALLEL DO
        call chemistry_debug_message('Kinetics thermodynamics driver success')
! Kinetics thermodynamics end
        call iamas_timer_stop('Kinetics thermodynamics driver')

    else
        call chemistry_message('--- time to run the kinetics thermodynamics driver: Dont do chemdt at this timestep')

    end if
!-----------------------------------------------------------------------------------------!
! wet scavenging start
!-----------------------------------------------------------------------------------------!
    call iamas_timer_start('wet scavenging')
if (config_wetscav_on) then
    write(debugString, '(a,l)') '--- time to run the wet scavenging driver config_wetscav_on = ', config_wetscav_on
    call chemistry_message(trim(debugString))
! wet scavenging start 
!$OMP PARALLEL DO
    do thread=1,nThreads
       call driver_wetscav(state, diag_physics, time_lev, block%configs, tend_chemistry, &
                          cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    call chemistry_debug_message('wetscav driver success')
    end do
!$OMP END PARALLEL DO
! wet scavenging end 
else
    write(debugString, '(a,l)') '--- time to run the wet scavenging driver config_wetscav_on = ', config_wetscav_on
    call chemistry_message(trim(debugString))
end if
    call iamas_timer_stop('wet scavenging')
!-----------------------------------------------------------------------------------------!
! Deallocate arrays shared by all chemistry parameterizations:
!-----------------------------------------------------------------------------------------!
    call deallocate_forall_chemistry
    call chemistry_debug_message('forall chemistry deallocate success')
    block => block % next
end do 

call iamas_timer_stop('chemistry driver')

 end subroutine chemistry_driver

!=================================================================================================================
 end module iamas_atmchem_driver
!=================================================================================================================
