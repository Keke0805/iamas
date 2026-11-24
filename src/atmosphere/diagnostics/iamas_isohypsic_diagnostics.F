! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module isohypsic_diagnostics

    use iamas_dmpar
    use iamas_kind_types
    use iamas_derived_types
    use iamas_pool_routines
    use iamas_constants
    use iamas_log, only : iamas_log_write

    type (iAMAS_pool_type), pointer :: mesh
    type (iAMAS_pool_type), pointer :: state
    type (iAMAS_pool_type), pointer :: diag

    type (iAMAS_clock_type), pointer :: clock

    public :: isohypsic_diagnostics_setup, &
              isohypsic_diagnostics_compute

    private

    logical :: need_t, need_q, need_u, need_v, need_p, need_w, need_qc


    contains


    !-----------------------------------------------------------------------
    !  routine isohypsic_diagnostics_setup
    !
    !> \brief Set up the isohypsic diagnostics module
    !> \author zyzhang
    !> \date   20 June 2025
    !> \details
    !>  This routine sets up the isohypsic diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine isohypsic_diagnostics_setup(all_pools, simulation_clock)

        use iamas_derived_types, only : iAMAS_pool_type, iAMAS_clock_type
        use iamas_pool_routines, only : iamas_pool_get_subpool

        implicit none

        type (iAMAS_pool_type), pointer :: all_pools
        type (iAMAS_clock_type), pointer :: simulation_clock

        clock => simulation_clock

        call iamas_pool_get_subpool(all_pools, 'mesh', mesh)
        call iamas_pool_get_subpool(all_pools, 'state', state)
        call iamas_pool_get_subpool(all_pools, 'diag', diag)
   
    end subroutine isohypsic_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine isohypsic_diagnostics_compute
    !
    !> \brief Compute isohypsic diagnostic before model output is written
    !> \author zyzhang
    !> \date   20 June 2025
    !> \details
    !>  Compute isohypsic diagnostic before model output is written. Code called
    !>  from here was previously in iamas_atm_interp_diagnostics.F.
    !
    !-----------------------------------------------------------------------
    subroutine isohypsic_diagnostics_compute()

        use iamas_atm_diagnostics_utils, only : iAMAS_field_will_be_written

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.

        need_t = iAMAS_field_will_be_written('temperature_isohypsic')
        need_any_diags = need_any_diags .or. need_t
        need_u = iAMAS_field_will_be_written('u_isohypsic')
        need_any_diags = need_any_diags .or. need_u
        need_v = iAMAS_field_will_be_written('v_isohypsic')
        need_any_diags = need_any_diags .or. need_v
        need_q = iAMAS_field_will_be_written('relhum_isohypsic')
        need_any_diags = need_any_diags .or. need_q
        need_p = iAMAS_field_will_be_written('pressure_isohypsic')
        need_any_diags = need_any_diags .or. need_p
        need_w = iAMAS_field_will_be_written('w_isohypsic')
        need_any_diags = need_any_diags .or. need_w
        need_qc = iAMAS_field_will_be_written('qc_isohypsic')
        need_any_diags = need_any_diags .or. need_qc

        if (need_any_diags) then
            call hinterp_diagnostics(mesh, state, 1, diag)
        end if
   
    end subroutine isohypsic_diagnostics_compute


   !==================================================================================================
    subroutine hinterp_diagnostics(mesh, state, time_lev, diag)
   !==================================================================================================
   
     !input arguments:
      type (iamas_pool_type), intent(in)  :: mesh
      type (iamas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev              ! which time level to use from state
     
     !inout arguments:
      type (iamas_pool_type), intent(inout) :: diag

     !local variables:
      integer :: iCell,iVert,iVertD,k,kk
      integer, pointer :: nCells, nCellsSolve, nVertLevels, nVertices,vertexDegree, nIsohypsicLevels
      integer :: nVertLevelsP1
      integer, pointer :: index_qv, num_scalars, index_qc
      integer, dimension(:,:), pointer :: cellsOnVertex
     
      type (field2DReal), pointer:: height_field
     
      real (kind=RKIND), dimension(:), pointer :: areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      
      real (kind=RKIND), dimension(:,:), pointer :: exner, zgrid
      real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p
      real (kind=RKIND), dimension(:,:), pointer :: relhum, theta_m, vorticity
      real (kind=RKIND), dimension(:,:), pointer :: umeridional, uzonal, vvel
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: qc

      real (kind=RKIND), dimension(:), pointer :: isohypsic_levels
      real (kind=RKIND), dimension(:,:), pointer :: temperature_isohypsic
      real (kind=RKIND), dimension(:,:), pointer :: u_isohypsic
      real (kind=RKIND), dimension(:,:), pointer :: v_isohypsic
      real (kind=RKIND), dimension(:,:), pointer :: relhum_isohypsic
      real (kind=RKIND), dimension(:,:), pointer :: pressure_isohypsic
      real (kind=RKIND), dimension(:,:), pointer :: w_isohypsic
      real (kind=RKIND), dimension(:,:), pointer :: qc_isohypsic

     !---------------------------------------------------------

      real (kind=RKIND), dimension(:,:), allocatable :: height, heightCp1, height2, height_v
      real (kind=RKIND), dimension(:,:), allocatable :: temperature, pressure

     !------------- local interpolated fields ------------------

      integer :: nIntH
      real (kind=RKIND), dimension(:,:), allocatable :: field_in,height_in
      real (kind=RKIND), dimension(:,:), allocatable :: field_interp,height_interp

     !----------------------------------------------------------

      call iamas_pool_get_dimension(mesh, 'nCells', nCells)
      call iamas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call iamas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call iamas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call iamas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
      call iamas_pool_get_dimension(mesh, 'nIsohypsicLevels', nIsohypsicLevels)
      call iamas_pool_get_dimension(state, 'index_qv', index_qv)
      call iamas_pool_get_dimension(state, 'index_qc', index_qc)
      call iamas_pool_get_dimension(state, 'num_scalars', num_scalars)

      nVertLevelsP1 = nVertLevels + 1

      call iamas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call iamas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call iamas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

      call iamas_pool_get_array(mesh, 'zgrid', zgrid)
      call iamas_pool_get_array(state, 'w', vvel, time_lev)
      call iamas_pool_get_array(state, 'theta_m', theta_m, time_lev)
      call iamas_pool_get_array(state, 'scalars', scalars, time_lev)
      qc => scalars(index_qc,:,:)

      call iamas_pool_get_field(mesh, 'zgrid', height_field)
      call iamas_dmpar_exch_halo_field(height_field)
       
      call iamas_pool_get_array(diag, 'exner', exner)
      call iamas_pool_get_array(diag, 'pressure_base', pressure_b)
      call iamas_pool_get_array(diag, 'pressure_p', pressure_p)
      call iamas_pool_get_array(diag, 'vorticity', vorticity)
      call iamas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
      call iamas_pool_get_array(diag, 'uReconstructZonal', uzonal)
      call iamas_pool_get_array(diag, 'relhum', relhum)

      call iamas_pool_get_array(diag, 'isohypsic_levels', isohypsic_levels)
      call iamas_pool_get_array(diag, 'temperature_isohypsic', temperature_isohypsic)
      call iamas_pool_get_array(diag, 'u_isohypsic', u_isohypsic)
      call iamas_pool_get_array(diag, 'v_isohypsic', v_isohypsic)
      call iamas_pool_get_array(diag, 'relhum_isohypsic', relhum_isohypsic)
      call iamas_pool_get_array(diag, 'pressure_isohypsic', pressure_isohypsic)
      call iamas_pool_get_array(diag, 'w_isohypsic', w_isohypsic)
      call iamas_pool_get_array(diag, 'qc_isohypsic', qc_isohypsic)

      if(.not.allocated(height) ) allocate(height(nVertLevels,nCells) )
      if(.not.allocated(heightCp1) ) allocate(heightCp1(nVertLevels,nCells+1) )
      if(.not.allocated(height2)   ) allocate(height2(nVertLevelsP1,nCells)   )
      if(.not.allocated(height_v)  ) allocate(height_v(nVertLevels,nVertices) )
      if(.not.allocated(temperature) ) allocate(temperature(nVertLevels,nCells)   )
      if(.not.allocated(pressure) ) allocate(pressure(nVertLevels,nCells)   )

      isohypsic_levels(1) = 100.0
      isohypsic_levels(2) = 130.0
      isohypsic_levels(3) = 150.0
      isohypsic_levels(4) = 170.0
      isohypsic_levels(5) = 190.0
      isohypsic_levels(6) = 210.0
      isohypsic_levels(7) = 230.0

     !calculation of height at cell centers (at mass points):
      do iCell = 1, nCells
      do k = 1, nVertLevels
         height(k,iCell) = 0.5*(zgrid(k,iCell) - zgrid(1,iCell) + zgrid(k+1,iCell) - zgrid(1,iCell))
         heightCp1(k,iCell) = height(k,iCell)
      enddo
      enddo
      do iCell = nCells+1, nCells+1
      do k = 1, nVertLevels
         heightCp1(k,iCell)   = 0.5*(zgrid(k,iCell) - zgrid(1,iCell) + zgrid(k+1,iCell) - zgrid(1,iCell))
      enddo
      enddo
     !calculation of height at cell centers (at vertical velocity points)
      do iCell = 1, nCells
      do k = 1, nVertLevelsP1
         height2(k,iCell) = zgrid(k,iCell) - zgrid(1,iCell)
      enddo
      enddo
     !calculation of height at cell vertices (at mass points):
      do iVert = 1, nVertices
        height_v(:,iVert) = 0._RKIND
    
        do k = 1, nVertLevels
        do iVertD = 1, vertexDegree
           height_v(k,iVert) = height_v(k,iVert) &
                   + kiteAreasOnVertex(iVertD,iVert)*heightCp1(k,cellsOnVertex(iVertD,iVert))
        enddo
        height_v(k,iVert) = height_v(k,iVert) / areaTriangle(iVert)
        enddo
      enddo

      if (need_t) then

       !calculation of temperature at cell centers:
        do iCell = 1, nCells
        do k = 1, nVertLevels
           temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell) 
        enddo
        enddo

        allocate(field_in(nCells, nVertLevels))
        allocate(height_in(nCells, nVertLevels))
        allocate(field_interp(nCells, nIsohypsicLevels))
        allocate(height_interp(nCells, nIsohypsicLevels))

        do iCell = 1, nCells
        do k = 1, nIsohypsicLevels
           height_interp(iCell,k) = isohypsic_levels(k)
        end do
        end do

        do iCell = 1, nCells
        do k = 1, nVertLevels
           field_in(iCell,k) = temperature(k,iCell)
        end do
        end do
        do iCell = 1, nCells
        do k = 1, nVertLevels
           height_in(iCell,k) = height(k,iCell)
        end do
        end do

        call interp_tofixed_height(nCells, nVertLevels, nIsohypsicLevels, height_in, field_in, height_interp, field_interp)

        do iCell = 1, nCells
        do k = 1, nIsohypsicLevels
           temperature_isohypsic(k,iCell) = field_interp(iCell,k)
        end do
        end do

        deallocate(field_in)
        deallocate(field_interp)
        deallocate(height_in)
        deallocate(height_interp)

      end if

      if (need_p) then

       !calculation of pressure at cell centers:
        do iCell = 1,nCells
        do k = 1,nVertLevels
           pressure(k,iCell)    = pressure_p(k,iCell) + pressure_b(k,iCell)
        enddo
        enddo

        allocate(field_in(nCells, nVertLevels))
        allocate(height_in(nCells, nVertLevels))
        allocate(field_interp(nCells, nIsohypsicLevels))
        allocate(height_interp(nCells, nIsohypsicLevels))

        do k=1,nIsohypsicLevels
            height_interp(1:nCells,k) = isohypsic_levels(k)
         end do
 
         do iCell=1,nCells
         do k=1,nVertLevels
            field_in(iCell,k) = pressure(k,iCell)
         end do
         end do
         do iCell=1,nCells
         do k=1,nVertLevels
            height_in(iCell,k) = height(k,iCell)
         end do
         end do
 
         call interp_tofixed_height(nCells, nVertLevels, nIsohypsicLevels, height_in, field_in, height_interp, field_interp)
 
         do k=1,nIsohypsicLevels
            pressure_isohypsic(k,1:nCells) = field_interp(1:nCells,k)
         end do
 
         deallocate(field_in)
         deallocate(field_interp)
         deallocate(height_in)
         deallocate(height_interp)

      end if

      if (need_u) then
 
         allocate(field_in(nCells, nVertLevels))
         allocate(height_in(nCells, nVertLevels))
         allocate(field_interp(nCells, nIsohypsicLevels))
         allocate(height_interp(nCells, nIsohypsicLevels))
 
         do k=1,nIsohypsicLevels
             height_interp(:,k) = isohypsic_levels(k)
          end do
  
          do iCell=1,nCells
          do k=1,nVertLevels
             field_in(iCell,k) = uzonal(k,iCell)
          end do
          end do
          do iCell=1,nCells
          do k=1,nVertLevels
             height_in(iCell,k) = height(k,iCell)
          end do
          end do
  
          call interp_tofixed_height(nCells, nVertLevels, nIsohypsicLevels, height_in, field_in, height_interp, field_interp)
  
          do k=1,nIsohypsicLevels
             u_isohypsic(k,1:nCells) = field_interp(1:nCells,k)
          end do
  
          deallocate(field_in)
          deallocate(field_interp)
          deallocate(height_in)
          deallocate(height_interp)
 
      end if

      if (need_v) then
 
        allocate(field_in(nCells, nVertLevels))
        allocate(height_in(nCells, nVertLevels))
        allocate(field_interp(nCells, nIsohypsicLevels))
        allocate(height_interp(nCells, nIsohypsicLevels))

        do k=1,nIsohypsicLevels
            height_interp(:,k) = isohypsic_levels(k)
         end do
 
         do iCell=1,nCells
         do k=1,nVertLevels
            field_in(iCell,k) = umeridional(k,iCell)
         end do
         end do
         do iCell=1,nCells
         do k=1,nVertLevels
            height_in(iCell,k) = height(k,iCell)
         end do
         end do
 
         call interp_tofixed_height(nCells, nVertLevels, nIsohypsicLevels, height_in, field_in, height_interp, field_interp)
 
         do k=1,nIsohypsicLevels
            v_isohypsic(k,1:nCells) = field_interp(1:nCells,k)
         end do
 
         deallocate(field_in)
         deallocate(field_interp)
         deallocate(height_in)
         deallocate(height_interp)

      end if

      if (need_q) then
 
        allocate(field_in(nCells, nVertLevels))
        allocate(height_in(nCells, nVertLevels))
        allocate(field_interp(nCells, nIsohypsicLevels))
        allocate(height_interp(nCells, nIsohypsicLevels))

        do k=1,nIsohypsicLevels
            height_interp(:,k) = isohypsic_levels(k)
         end do
 
         do iCell=1,nCells
         do k=1,nVertLevels
            field_in(iCell,k) = relhum(k,iCell)
         end do
         end do
         do iCell=1,nCells
         do k=1,nVertLevels
            height_in(iCell,k) = height(k,iCell)
         end do
         end do
 
         call interp_tofixed_height(nCells, nVertLevels, nIsohypsicLevels, height_in, field_in, height_interp, field_interp)
 
         do k=1,nIsohypsicLevels
            relhum_isohypsic(k,1:nCells) = field_interp(1:nCells,k)
         end do
 
         deallocate(field_in)
         deallocate(field_interp)
         deallocate(height_in)
         deallocate(height_interp)

      end if

      if (need_w) then
 
         allocate(field_in(nCells, nVertLevelsP1))
         allocate(height_in(nCells, nVertLevelsP1))
         allocate(field_interp(nCells, nIsohypsicLevels))
         allocate(height_interp(nCells, nIsohypsicLevels))
 
         do k=1,nIsohypsicLevels
             height_interp(:,k) = isohypsic_levels(k)
          end do
  
          do iCell=1,nCells
          do k=1,nVertLevelsP1
             field_in(iCell,k) = vvel(k,iCell)
          end do
          end do
          do iCell=1,nCells
          do k=1,nVertLevelsP1
             height_in(iCell,k) = height2(k,iCell)
          end do
          end do
  
          call interp_tofixed_height(nCells, nVertLevelsP1, nIsohypsicLevels, height_in, field_in, height_interp, field_interp)
  
          do k=1,nIsohypsicLevels
             w_isohypsic(k,1:nCells) = field_interp(1:nCells,k)
          end do
  
          deallocate(field_in)
          deallocate(field_interp)
          deallocate(height_in)
          deallocate(height_interp)
 
       end if

       if (need_qc) then
 
         allocate(field_in(nCells, nVertLevels))
         allocate(height_in(nCells, nVertLevels))
         allocate(field_interp(nCells, nIsohypsicLevels))
         allocate(height_interp(nCells, nIsohypsicLevels))
 
         do k=1,nIsohypsicLevels
             height_interp(:,k) = isohypsic_levels(k)
          end do
  
          do iCell=1,nCells
          do k=1,nVertLevels
             field_in(iCell,k) = qc(k,iCell)
          end do
          end do
          do iCell=1,nCells
          do k=1,nVertLevels
             height_in(iCell,k) = height(k,iCell)
          end do
          end do
  
          call interp_tofixed_height(nCells, nVertLevels, nIsohypsicLevels, height_in, field_in, height_interp, field_interp)
  
          do k=1,nIsohypsicLevels
             qc_isohypsic(k,1:nCells) = field_interp(1:nCells,k)
          end do
  
          deallocate(field_in)
          deallocate(field_interp)
          deallocate(height_in)
          deallocate(height_interp)
 
       end if

      if(allocated(temperature) ) deallocate(temperature )
      if(allocated(height2)     ) deallocate(height2     )
      if(allocated(height)      ) deallocate(height      )
      if(allocated(heightCp1)   ) deallocate(heightCp1   )
      if(allocated(height_v)    ) deallocate(height_v    )
      if(allocated(pressure)    ) deallocate(pressure    )     
  
     end subroutine hinterp_diagnostics

   !==================================================================================================
    subroutine interp_tofixed_height(ncol,nlev_in,nlev_out,heig_in,field_in,heig_out,field_out)
   !==================================================================================================

   !input arguments:
    integer,intent(in):: ncol,nlev_in,nlev_out
   
    real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: heig_in,field_in
    real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: heig_out
     
   !output arguments:
    real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out
     
   !local variables:
   ! integer:: i1,i2,icol,k,kk
    integer:: icol,k,kk
    integer:: kkstart,kount
    integer,dimension(ncol):: kupper
     
    real(kind=RKIND):: dpl,dpu

   !--------------------------------------------------------------------------------------------------

    do icol=1,ncol
       kupper(icol) = 1
    enddo

    do k=1,nlev_out

      kkstart = nlev_in
      do icol=1,ncol
         kkstart = min0(kkstart,kupper(icol))
      enddo
      kount = 0

      do kk=kkstart,nlev_in-1
         do icol=1,ncol
            if(heig_out(icol,k).gt.heig_in(icol,kk).and.heig_out(icol,k).le.heig_in(icol,kk+1))then
               kupper(icol) = kk
               kount = kount + 1
            endif
         enddo

         if (kount.eq.ncol) then
            do icol=1,ncol
               dpu = heig_out(icol,k) - heig_in(icol,kupper(icol))
               dpl = heig_in(icol,kupper(icol)+1) - heig_out(icol,k)
               field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                                 + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
            enddo
            goto 35
         end if
      enddo

      do icol = 1, ncol
         if(heig_out(icol,k) .lt. heig_in(icol,1)) then
            field_out(icol,k) = field_in(icol,1)
         elseif(heig_out(icol,k) .gt. heig_in(icol,nlev_in)) then
            !field_out(icol,k) = field_in(icol,1)*heig_out(icol,k)/heig_in(icol,1)
            field_out(icol,k) = field_in(icol,nlev_in)
         else
            dpu = heig_out(icol,k) - heig_in(icol,kupper(icol))
            dpl = heig_in(icol,kupper(icol)+1) - heig_out(icol,k)
            field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                              + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
         endif
      enddo

    35 continue

    enddo

    end subroutine interp_tofixed_height

end module isohypsic_diagnostics
