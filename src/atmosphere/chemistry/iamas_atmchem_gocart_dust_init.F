!************************************************************************
!
! Module to initialize varaiables for gocart dust scheme.
! * Author: Jiawang Feng and Chun Zhao / AEMOL USTC
! Last update: Dec 2019
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! Please report any bugs or problems to Chun Zhao
!
! References: 
!   WRF-Chem-v3.5:chem/module_gocart_dust_init.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

Module iamas_atmchem_gocart_dust_init

    use module_data_gocart_dust
    use iamas_kind_types
    use iamas_pool_routines


    implicit none
    private
    public:: gocart_dust_init_foriAMAS

    logical,public,save :: has_done_gocart_init = .false.         !Greenwich mean time hour of model start (hr)

    contains

    subroutine gocart_dust_init_foriAMAS(configs)
        integer :: istat,lc,slcats,iindex
        character(len=StrKIND),pointer:: mminsl
        character(4) :: sltype
        type(iamas_pool_type),intent(in):: configs

        call iamas_pool_get_config(configs,'input_soil_data',mminsl)

    !read in soil properties from soilparm.tbl:
        open(16,file='SOILPARM.TBL',form='FORMATTED',status='OLD',iostat=istat)

        read(16,*)
        read(16,2000) sltype
        2000 format(a4)
        read(16,*)slcats,iindex

    !prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008:
        if(size(bb    ) < slcats .or. &
           size(drysmc) < slcats .or. &
           size(f11   ) < slcats .or. &
           size(maxsmc) < slcats .or. &
           size(refsmc) < slcats .or. &
           size(satpsi) < slcats .or. &
           size(satdk ) < slcats .or. &
           size(satdw ) < slcats .or. &
           size(wltsmc) < slcats .or. &
           size(qtz   ) < slcats) then
           write(0,*) 'table sizes too small for value of slcats in gocart'
        endif

        if(sltype .eq. mminsl) then
           do lc = 1, slcats
              read(16,*) iindex,bb(lc),drysmc(lc),f11(lc),maxsmc(lc),refsmc(lc),satpsi(lc), &
                         satdk(lc),satdw(lc),wltsmc(lc),qtz(lc)
           enddo
        endif

        close(16)
    end subroutine

END Module iamas_atmchem_gocart_dust_init
