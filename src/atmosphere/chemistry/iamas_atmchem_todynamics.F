!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to pass all chemistry tendencies from chemistry to dynamic modules.
! * Author: Jiawang Feng from AEMOL USTC
! Last update: jwfeng 11/2019
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! References: 
!     iAMAS-v5.2:src/core_atmosphere/physics/iamas_atmphys_todynamics.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module iamas_atmchem_todynamics
 use iamas_timer, only: iamas_timer_start, iamas_timer_stop
 use iamas_kind_types
 use iamas_derived_types
 use iamas_pool_routines
 use iamas_dmpar
 use iamas_atmchem_vars
 use iamas_atmphys_vars

 implicit none
 private
 public:: chemistry_get_tend


!Interface between the chemistry schemes and the non-hydrostatic dynamical core.
!
! subroutines in iamas_atmchem_todynamics:
! ---------------------------------------
! chemistry_get_tend: add and mass-weigh tendencies before being added to dynamics tendencies.
! tend_toEdges      : interpolate wind-tendencies from centers to edges of grid-cells.
! ----------------------------------------

 contains
 

!=================================================================================================================
 subroutine chemistry_get_tend( mesh, state, tend, tend_chemistry )
!=================================================================================================================
 use iamas_atm_dimensions
!input variables:
 type(iamas_pool_type),intent(in):: mesh
 type(iamas_pool_type),intent(in):: state
 type(iamas_pool_type),intent(in):: tend_chemistry
!inout variables:
 type(iamas_pool_type),intent(inout):: tend

!local variables:
 integer        :: i,iCell,k,n
 integer,pointer:: nCells,nCellsSolve

 real(kind=RKIND),dimension(:,:),  pointer :: mass          ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:,:),pointer :: tend_chem

!=================================================================================================================
 call iamas_timer_start('chem_src_get_fields')
 call iamas_pool_get_dimension(mesh, 'nCells',         nCells        )
 call iamas_pool_get_dimension(mesh, 'nCellsSolve',    nCellsSolve   )
 call iamas_pool_get_array(state,'rho_zz',   mass, 2)
 call iamas_pool_get_array(tend ,'chem_tend',tend_chem)
 call iamas_timer_stop('chem_src_get_fields')

 call iamas_timer_start('chem_src_clear')
 tend_chem(:,:,:) = 0._RKIND
 call iamas_timer_stop('chem_src_clear')
 call chemistry_get_tend_work(state,nCells,nCellsSolve, &
                           tend_chemistry, &
                           mass,tend_chem  &
                           )

 end subroutine chemistry_get_tend

 !==================================================================================================
 subroutine chemistry_get_tend_work(state,nCells,nCellsSolve, &
                           tend_chemistry, &
                           mass,tend_chem  &
                           )

!==================================================================================================

    use iamas_atm_dimensions
 
    implicit none
    type(iamas_pool_type),intent(in) :: state
    type(iamas_pool_type),intent(in) :: tend_chemistry
    integer, intent(in) :: nCells, nCellsSolve
    real(kind=RKIND),dimension(nVertLevels,nCells+1), intent(in) :: mass
    real(kind=RKIND),dimension(num_chem,nVertLevels,nCells+1), intent(inout) :: tend_chem
! local variables
    real(kind=RKIND),dimension(:,:,:),pointer :: tend_emis
    integer :: i, k, h
 
    ! add this for update chem with emission. @zhxia 0418/2025
 call iamas_timer_start('chem_src_get_emission')
    call iamas_pool_get_array(tend_chemistry,'tend_emis'     ,tend_emis     )
 call iamas_timer_stop('chem_src_get_emission')
 call iamas_timer_start('chem_src_update')
    do h = 1, num_chem
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_chem(h,k,i) = tend_chem(h,k,i) + tend_emis(h,k,i)*mass(k,i)
    enddo
    enddo
    enddo
 call iamas_timer_stop('chem_src_update')
 end subroutine chemistry_get_tend_work
 

!=================================================================================================================
 end module iamas_atmchem_todynamics
!=================================================================================================================
