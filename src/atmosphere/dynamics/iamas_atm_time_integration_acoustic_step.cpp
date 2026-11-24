/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Xinqu develop acoustic_step subroutine in kokkos

    * Author: Xin Qu
    Last update: Nov 2025

    Contact:
    Chun Zhao, PhD
    Professor
    University of Science and Technology of China
    Hefei, Anhui, 230026
    Email: ChunZhao@ustc.edu.cn

    Please report any bugs or problems to Chun Zhao

    Additional information:
    aemol.ustc.edu.cn: resources
**/
#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>
#include <flcl-cxx.hpp>
#include <cstring>
#include <iostream>
typedef double real;
using real1d = Kokkos::Experimental::OffsetView<real*, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real2d = Kokkos::Experimental::OffsetView<real**, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real3d = Kokkos::Experimental::OffsetView<real***, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real4d = Kokkos::Experimental::OffsetView<real****, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real5d = Kokkos::Experimental::OffsetView<real*****, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int2d = Kokkos::Experimental::OffsetView<int**, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int5d = Kokkos::Experimental::OffsetView<int*****, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int1d = Kokkos::Experimental::OffsetView<int*, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using int3d = Kokkos::Experimental::OffsetView<int***, Kokkos::LayoutLeft,flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;

extern "C" void atm_advance_acoustic_step_kokkos(int &maxEdges,int &nVertLevels,int &nCells, int &nEdges, int &nCellsSolve, int &cellStart, int &cellEnd, int &vertexStart, 
    int &vertexEnd, int &edgeStart, int &edgeEnd, int &cellSolveStart, int &cellSolveEnd, int &vertexSolveStart, int &vertexSolveEnd, int &edgeSolveStart, int &edgeSolveEnd, 
    flcl_ndarray_t* rho_zz_p, flcl_ndarray_t* theta_m_p, flcl_ndarray_t* ru_p_f, flcl_ndarray_t* rw_p_f, flcl_ndarray_t* rtheta_pp_p, flcl_ndarray_t* rtheta_pp_old_p, 
    flcl_ndarray_t* zz_p, flcl_ndarray_t* exner_p, flcl_ndarray_t* cqu_p, flcl_ndarray_t* ruAvg_p, flcl_ndarray_t* wwAvg_p, 
    flcl_ndarray_t* rho_pp_p, flcl_ndarray_t* cofwt_p, flcl_ndarray_t* coftz_p, flcl_ndarray_t* zxu_p, flcl_ndarray_t* a_tri_p, flcl_ndarray_t* alpha_tri_p, flcl_ndarray_t* gamma_tri_p, 
    flcl_ndarray_t* dss_p, flcl_ndarray_t* tend_ru_p, flcl_ndarray_t* tend_rho_p, flcl_ndarray_t* tend_rt_p, 
    flcl_ndarray_t* tend_rw_p, flcl_ndarray_t* zgrid_p, flcl_ndarray_t* cofwr_p, flcl_ndarray_t* cofwz_p, flcl_ndarray_t* w_p, flcl_ndarray_t* ru_f, flcl_ndarray_t* ru_save_p, 
    flcl_ndarray_t* rw_f, flcl_ndarray_t* rw_save_p, flcl_ndarray_t* fzm_p, flcl_ndarray_t* fzp_p, flcl_ndarray_t* rdzw_p, flcl_ndarray_t* dcEdge_p, flcl_ndarray_t* invDcEdge_p, 
    flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* cofrz_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* cellsOnEdge_p, 
    flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnCell_sign_p, 
    real &dts, int &small_step, flcl_ndarray_t* epssm_p, real &cf1, real &cf2, real &cf3, 
    flcl_ndarray_t* specZoneMaskEdge_p, flcl_ndarray_t* specZoneMaskCell_p){
    auto rho_zz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_zz_p); auto rho_zz = real2d(rho_zz_raw,{1,1});
    auto theta_m_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*theta_m_p); auto theta_m = real2d(theta_m_raw,{1,1});
    auto ru_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_p_f); auto ru_p = real2d(ru_p_raw,{1,1});
    auto rw_p_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_p_f); auto rw_p = real2d(rw_p_raw,{1,1});
    auto rtheta_pp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rtheta_pp_p); auto rtheta_pp = real2d(rtheta_pp_raw,{1,1});
    auto rtheta_pp_old_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rtheta_pp_old_p); auto rtheta_pp_old = real2d(rtheta_pp_old_raw,{1,1});
    auto zz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zz_p); auto zz = real2d(zz_raw,{1,1});
    auto exner_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*exner_p); auto exner = real2d(exner_raw,{1,1});
    auto cqu_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*cqu_p); auto cqu = real2d(cqu_raw,{1,1});
    auto ruAvg_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ruAvg_p); auto ruAvg = real2d(ruAvg_raw,{1,1});
    auto wwAvg_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*wwAvg_p); auto wwAvg = real2d(wwAvg_raw,{1,1});
    auto rho_pp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rho_pp_p); auto rho_pp = real2d(rho_pp_raw,{1,1});
    auto cofwt_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*cofwt_p); auto cofwt = real2d(cofwt_raw,{1,1});
    auto coftz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*coftz_p); auto coftz = real2d(coftz_raw,{1,1});
    auto zxu_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zxu_p); auto zxu = real2d(zxu_raw,{1,1});
    auto a_tri_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*a_tri_p); auto a_tri = real2d(a_tri_raw,{1,1});
    auto alpha_tri_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*alpha_tri_p); auto alpha_tri = real2d(alpha_tri_raw,{1,1});
    auto gamma_tri_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*gamma_tri_p); auto gamma_tri = real2d(gamma_tri_raw,{1,1});
    auto dss_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*dss_p); auto dss = real2d(dss_raw,{1,1});
    auto tend_ru_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_ru_p); auto tend_ru = real2d(tend_ru_raw,{1,1});
    auto tend_rho_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_rho_p); auto tend_rho = real2d(tend_rho_raw,{1,1});
    auto tend_rt_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_rt_p); auto tend_rt = real2d(tend_rt_raw,{1,1});
    auto tend_rw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_rw_p); auto tend_rw = real2d(tend_rw_raw,{1,1});
    auto zgrid_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*zgrid_p); auto zgrid = real2d(zgrid_raw,{1,1});
    auto cofwr_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*cofwr_p); auto cofwr = real2d(cofwr_raw,{1,1});
    auto cofwz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*cofwz_p); auto cofwz = real2d(cofwz_raw,{1,1});
    auto w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*w_p); auto w = real2d(w_raw,{1,1});
    auto ru_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_f); auto ru = real2d(ru_raw,{1,1});
    auto ru_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ru_save_p); auto ru_save = real2d(ru_save_raw,{1,1});
    auto rw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_f); auto rw = real2d(rw_raw,{1,1});
    auto rw_save_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rw_save_p); auto rw_save = real2d(rw_save_raw,{1,1});
    auto fzm_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fzm_p); auto fzm = real1d(fzm_raw,{1});
    auto fzp_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fzp_p); auto fzp = real1d(fzp_raw,{1});
    auto rdzw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rdzw_p); auto rdzw = real1d(rdzw_raw,{1});
    auto dcEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dcEdge_p); auto dcEdge = real1d(dcEdge_raw,{1});
    auto invDcEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invDcEdge_p); auto invDcEdge = real1d(invDcEdge_raw,{1});
    auto invAreaCell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p); auto invAreaCell = real1d(invAreaCell_raw,{1});
    auto cofrz_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*cofrz_p); auto cofrz = real1d(cofrz_raw,{1});
    auto dvEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p); auto dvEdge = real1d(dvEdge_raw,{1});
    auto nEdgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p); auto nEdgesOnCell = int1d(nEdgesOnCell_raw,{1});
    auto cellsOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p); auto cellsOnEdge = int2d(cellsOnEdge_raw,{1,1});
    auto edgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p); auto edgesOnCell = int2d(edgesOnCell_raw,{1,1});
    auto edgesOnCell_sign_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p); auto edgesOnCell_sign = real2d(edgesOnCell_sign_raw,{1,1});
    auto epssm_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*epssm_p); auto epssm = real2d(epssm_raw,{1,1});
    auto specZoneMaskEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*specZoneMaskEdge_p); auto specZoneMaskEdge = real1d(specZoneMaskEdge_raw,{1});
    auto specZoneMaskCell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*specZoneMaskCell_p); auto specZoneMaskCell = real1d(specZoneMaskCell_raw,{1});
    real resm_s[nVertLevels * (nCells + 1)];
    auto resm = real2d{(real*)resm_s,{1,1},{nVertLevels + 1, nCells + 2}};
    
    const real rgas = 287.0;
    const real cp = 7. * rgas / 2.0;
    const real cv = cp - rgas;
    const real gravity = 9.80616;
    real rcv = rgas / (cp - rgas);
    real c2 = cp * rcv;
    // for(int iCell = 1; iCell <= nCells + 1; iCell++){
    Kokkos::parallel_for("Init_Resm", Kokkos::RangePolicy<>(1, nCells + 2), 
        KOKKOS_LAMBDA(const int iCell) {
        for(int k = 1; k <= nVertLevels; k++){
            resm(k, iCell) = (1.0 - epssm(k, iCell)) / (1.0 + epssm(k, iCell));
        }
    // }
    });
    Kokkos::fence();
    real rdts = 1./dts;
    if(small_step != 1){ //not needed on first small step
        /**
         * forward-backward acoustic step integration.
         * begin by updating the horizontal velocity u, 
         * and accumulating the contribution from the updated u to the other tendencies.
         * we are looping over all edges, but only computing on edges of owned cells. This will include updates of
         * all owned edges plus some edges that are owned by other blocks.  We perform these redundant computations
         * so that we do not have to communicate updates of u to update the cell variables (rho, w, and theta). 
         * MGD this loop will not be very load balanced with if-test below
         */
        // for(int iEdge = edgeStart; iEdge <= edgeEnd; iEdge++){ // MGD do we really just need edges touching owned cells? 
        Kokkos::parallel_for("Update_Ru_Forward", Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), 
        KOKKOS_LAMBDA(const int iEdge) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            // update edges for block-owned cells
            if(cell1 <= nCellsSolve || cell2 <= nCellsSolve){
                for(int k = 1; k <= nVertLevels; k++){
                    real pgrad = ((rtheta_pp(k,cell2)-rtheta_pp(k,cell1))*invDcEdge(iEdge) )/(.5*(zz(k,cell2)+zz(k,cell1)));
                    pgrad = cqu(k,iEdge)*0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad;
                    pgrad = pgrad + 0.5*zxu(k,iEdge)*gravity*(rho_pp(k,cell1)+rho_pp(k,cell2));
                    ru_p(k,iEdge) = ru_p(k,iEdge) + dts*(tend_ru(k,iEdge) - (1.0 - specZoneMaskEdge(iEdge))*pgrad);
                }

                // accumulate ru_p for use later in scalar transport
                for(int k = 1; k <= nVertLevels; k++){
                    ruAvg(k,iEdge) = ruAvg(k,iEdge) + ru_p(k,iEdge);
                }
            } // end test for block-owned cells
        // } // end loop over edges
        });
        Kokkos::fence();
    }else{ // this is all that us needed for ru_p update for first acoustic step in RK substep
        // for(int iEdge = edgeStart; iEdge <= edgeEnd; iEdge++){ // MGD do we really just need edges touching owned cells?
        Kokkos::parallel_for("Update_Ru_First", Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), 
        KOKKOS_LAMBDA(const int iEdge) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            // update edges for block-owned cells
            if(cell1 <= nCellsSolve || cell2 <= nCellsSolve){
                for(int k = 1; k <= nVertLevels; k++){
                    ru_p(k,iEdge) = dts*tend_ru(k,iEdge);
                }

                for(int k = 1; k <= nVertLevels; k++){
                    ruAvg(k,iEdge) = ru_p(k,iEdge);
                }
            } //end test for block-owned cells
        // } //end loop over edges
        });
        Kokkos::fence();
    } // test for first acoustic step
    Kokkos::parallel_for("Init_Theta_Old", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), 
    KOKKOS_LAMBDA(const int iCell) {
        if(small_step == 1){ // initialize here on first small timestep.
            for(int k = 1; k <= nVertLevels; k++){
                rtheta_pp_old(k,iCell) = 0.0;
            }
        }else{
            for(int k = 1; k <= nVertLevels; k++){
                rtheta_pp_old(k,iCell) = rtheta_pp(k,iCell);
            }
        }
    });
    Kokkos::fence();

    // for(int iCell = cellSolveStart; iCell <= cellSolveEnd; iCell++){ // loop over all owned cells to solve
    Kokkos::parallel_for("Vertical_Implicit_Solve", Kokkos::RangePolicy<>(cellSolveStart, cellSolveEnd + 1), 
        KOKKOS_LAMBDA(const int iCell) {
        real rs_s[nVertLevels];
        real ts_s[nVertLevels];
        auto rs = real1d((real*)rs_s, {1},{nVertLevels + 1});
        auto ts = real1d((real*)ts_s, {1},{nVertLevels + 1});
        if(small_step == 1){ // initialize here on first small timestep.
            for(int k = 1; k <= nVertLevels; k++){
                wwAvg(k,iCell) = 0.0;         
                rho_pp(k,iCell) = 0.0;            
                rtheta_pp(k,iCell) = 0.0;            
                rw_p(k,iCell) = 0.0;
            }
            wwAvg(nVertLevels + 1,iCell) = 0.0;
            rw_p(nVertLevels + 1,iCell) = 0.0;
        }
        if(specZoneMaskCell(iCell) == 0.0){ // not specified zone, compute...
            for(int k = 1; k <= nVertLevels; k++){
                ts(k) = 0.0;
                rs(k) = 0.0;
            }
            for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                int iEdge = edgesOnCell(i,iCell);
                int cell1 = cellsOnEdge(1,iEdge);
                int cell2 = cellsOnEdge(2,iEdge);
                for(int k = 1; k <= nVertLevels; k++){
                    real flux = edgesOnCell_sign(i,iCell)*dts*dvEdge(iEdge)*ru_p(k,iEdge) * invAreaCell(iCell);
                    rs(k) = rs(k)-flux;
                    ts(k) = ts(k)-flux*0.5*(theta_m(k,cell2)+theta_m(k,cell1));
                }
            }
            /**
             * vertically implicit acoustic and gravity wave integration.
             * this follows Klemp et al MWR 2007, with the addition of an implicit Rayleigh damping of w
             * serves as a gravity-wave absorbing layer, from Klemp et al 2008.
             */
            for(int k = 1; k <= nVertLevels; k++){
                rs(k) = rho_pp(k,iCell) + dts*tend_rho(k,iCell) + rs(k)                  
                            - cofrz(k)*resm(k,iCell)*(rw_p(k+1,iCell)-rw_p(k,iCell)); 
                ts(k) = rtheta_pp(k,iCell) + dts*tend_rt(k,iCell) + ts(k)                
                                - resm(k,iCell)*rdzw(k)*( coftz(k+1,iCell)*rw_p(k+1,iCell)     
                                                -coftz(k,iCell)*rw_p(k,iCell));
            }
            for(int k = 2; k <= nVertLevels; k++){
                wwAvg(k,iCell) = wwAvg(k,iCell) + 0.5*(1.0-epssm(k,iCell))*rw_p(k,iCell);
            }
            for(int k = 2; k <= nVertLevels; k++){
                rw_p(k,iCell) = rw_p(k,iCell) +  dts*tend_rw(k,iCell)                       
                       - cofwz(k,iCell)*((zz(k  ,iCell)*ts(k)                           
                                     -zz(k-1,iCell)*ts(k-1))                            
                               +resm(k,iCell)*(zz(k  ,iCell)*rtheta_pp(k  ,iCell)                
                                     -zz(k-1,iCell)*rtheta_pp(k-1,iCell)))              
                       - cofwr(k,iCell)*((rs(k)+rs(k-1))                                
                               +resm(k,iCell)*(rho_pp(k,iCell)+rho_pp(k-1,iCell)))               
                       + cofwt(k  ,iCell)*(ts(k  )+resm(k,iCell)*rtheta_pp(k  ,iCell))           
                       + cofwt(k-1,iCell)*(ts(k-1)+resm(k,iCell)*rtheta_pp(k-1,iCell));
            }
            // tridiagonal solve sweeping up and then down the column
            for(int k = 2; k <= nVertLevels; k++){
                rw_p(k,iCell) = (rw_p(k,iCell)-a_tri(k,iCell)*rw_p(k-1,iCell))*alpha_tri(k,iCell);
            }
            for(int k = nVertLevels; k >= 1; k--){
                rw_p(k,iCell) = rw_p(k,iCell) - gamma_tri(k,iCell)*rw_p(k+1,iCell);
            }
            //  the implicit Rayleigh damping on w (gravity-wave absorbing) 
            for(int k = 2; k <= nVertLevels; k++){
                rw_p(k,iCell) = (rw_p(k,iCell) + (rw_save(k  ,iCell) - rw(k  ,iCell)) -dts*dss(k,iCell)*               
                        (fzm(k)*zz (k,iCell)+fzp(k)*zz (k-1,iCell))        
                        *(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell))       
                                 *w(k,iCell)    )/(1.0+dts*dss(k,iCell)) 
                         - (rw_save(k  ,iCell) - rw(k  ,iCell));
            }
            // accumulate (rho*omega)' for use later in scalar transport
            for(int k = 2; k <= nVertLevels; k++){
                wwAvg(k,iCell) = wwAvg(k,iCell) + 0.5*(1.0+epssm(k,iCell))*rw_p(k,iCell);
            }
            // update rho_pp and theta_pp given updated rw_p
            for(int k = 1; k <= nVertLevels; k++){
                rho_pp(k,iCell) = rs(k) - cofrz(k) *(rw_p(k+1,iCell)-rw_p(k  ,iCell));
                rtheta_pp(k,iCell) = ts(k) - rdzw(k)*(coftz(k+1,iCell)*rw_p(k+1,iCell)  
                                -coftz(k  ,iCell)*rw_p(k  ,iCell));
            }
        }else{ // specifed zone in regional_iAMAS
            for(int k = 1; k <= nVertLevels; k++){
                rho_pp(k,iCell) = rho_pp(k,iCell) + dts*tend_rho(k,iCell);
                rtheta_pp(k,iCell) = rtheta_pp(k,iCell) + dts*tend_rt(k,iCell);
                rw_p(k,iCell) = rw_p(k,iCell) + dts*tend_rw(k,iCell);
                wwAvg(k,iCell) = wwAvg(k,iCell) + 0.5*(1.0+epssm(k,iCell))*rw_p(k,iCell);
            }
        }
    // } // end of loop over cells
    });
    Kokkos::fence();
}
