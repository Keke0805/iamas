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

extern "C" void physics_get_tend_phase1_kokkos(int &nCells, int &nVertLevels, int &nEdges, int &num_scalars, int &nCellsSolve, int &nEdgesSolve,
    flcl_ndarray_t* tend_u_p,flcl_ndarray_t* rublten_Edge_p,flcl_ndarray_t* mass_edge_p,flcl_ndarray_t* tend_th_p,flcl_ndarray_t* tend_scalars_p,flcl_ndarray_t* rqvblten_p,flcl_ndarray_t* rqcblten_p,flcl_ndarray_t* rqiblten_p,
    int &index_qv,int &index_qc,int &index_qi,flcl_ndarray_t* mass_p,flcl_ndarray_t* rthblten_p){
    auto tend_u_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_u_p); auto tend_u = real2d(tend_u_raw,{1,1});
    auto rublten_Edge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rublten_Edge_p); auto rublten_Edge = real2d(rublten_Edge_raw,{1,1});
    auto mass_edge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*mass_edge_p); auto mass_edge = real2d(mass_edge_raw,{1,1});
    auto tend_th_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_th_p); auto tend_th = real2d(tend_th_raw,{1,1});
    auto tend_scalars_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*tend_scalars_p); auto tend_scalars = real3d(tend_scalars_raw,{1,1,1});
    auto rqvblten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rqvblten_p); auto rqvblten = real2d(rqvblten_raw,{1,1});
    auto rqcblten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rqcblten_p); auto rqcblten = real2d(rqcblten_raw,{1,1});
    auto rqiblten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rqiblten_p); auto rqiblten = real2d(rqiblten_raw,{1,1});
    auto mass_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*mass_p); auto mass = real2d(mass_raw,{1,1});
    auto rthblten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rthblten_p); auto rthblten = real2d(rthblten_raw,{1,1});
    Kokkos::parallel_for("physics_get_tend_phase1_kokkos", Kokkos::RangePolicy<>(1, nEdgesSolve + 1), 
        KOKKOS_LAMBDA(const int i) {
        for(int k = 1; k <= nVertLevels; k++){
            tend_u(k,i)=tend_u(k,i)+rublten_Edge(k,i)*mass_edge(k,i);
        }
    // }
    });
    Kokkos::fence();
    // for(int i = 1; i <= nCellsSolve; i++){
    Kokkos::parallel_for("physics_get_tend_phase1_kokkos", Kokkos::RangePolicy<>(1, nCellsSolve + 1), 
        KOKKOS_LAMBDA(const int i) {
        for(int k = 1; k <= nVertLevels; k++){
            tend_th(k,i) = tend_th(k,i) + rthblten(k,i)*mass(k,i);
            tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvblten(k,i)*mass(k,i);
            tend_scalars(index_qc,k,i) = tend_scalars(index_qc,k,i) + rqcblten(k,i)*mass(k,i);
            tend_scalars(index_qi,k,i) = tend_scalars(index_qi,k,i) + rqiblten(k,i)*mass(k,i);
        }
    // }
    });
    Kokkos::fence();
}

extern "C" void physics_get_tend_phase2_kokkos(int &nCells, int &nVertLevels, int &nEdges, int &num_scalars, int &nCellsSolve,
    flcl_ndarray_t* tend_th_p,flcl_ndarray_t* tend_scalars_p,
    int &index_qv,int &index_qc,int &index_qi,flcl_ndarray_t* mass_p,flcl_ndarray_t* rthcuten_p,flcl_ndarray_t* rqvcuten_p,flcl_ndarray_t* rqccuten_p,flcl_ndarray_t* rqicuten_p){
    auto tend_th_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_th_p); auto tend_th = real2d(tend_th_raw,{1,1});
    auto tend_scalars_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*tend_scalars_p); auto tend_scalars = real3d(tend_scalars_raw,{1,1,1});
    auto mass_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*mass_p); auto mass = real2d(mass_raw,{1,1});
    auto rthcuten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rthcuten_p); auto rthcuten = real2d(rthcuten_raw,{1,1});
    auto rqvcuten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rqvcuten_p); auto rqvcuten = real2d(rqvcuten_raw,{1,1});
    auto rqccuten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rqccuten_p); auto rqccuten = real2d(rqccuten_raw,{1,1});
    auto rqicuten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rqicuten_p); auto rqicuten = real2d(rqicuten_raw,{1,1});
    // for(int i = 1; i <= nCellsSolve; i++){
    Kokkos::parallel_for("physics_get_tend_phase2_kokkos", Kokkos::RangePolicy<>(1, nCellsSolve + 1), 
        KOKKOS_LAMBDA(const int i) {
        for(int k = 1; k <= nVertLevels; k++){
            tend_th(k,i) = tend_th(k,i) + rthcuten(k,i)*mass(k,i);
            tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvcuten(k,i)*mass(k,i);
            tend_scalars(index_qc,k,i) = tend_scalars(index_qc,k,i) + rqccuten(k,i)*mass(k,i);
            tend_scalars(index_qi,k,i) = tend_scalars(index_qi,k,i) + rqicuten(k,i)*mass(k,i);
        }
    // }
    });
    Kokkos::fence();
}

extern "C" void physics_get_tend_phase3_kokkos(int &nVertLevels, int &nCells, int &nEdges, int &num_scalars, int &nCellsSolve,
    flcl_ndarray_t* tend_scalars_p,flcl_ndarray_t* tend_th_p,
    int &index_qv,flcl_ndarray_t* mass_p,flcl_ndarray_t* rthfddaten_p,flcl_ndarray_t* rqvfddaten_p){
    auto tend_th_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_th_p); auto tend_th = real2d(tend_th_raw,{1,1});
    auto tend_scalars_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*tend_scalars_p); auto tend_scalars = real3d(tend_scalars_raw,{1,1,1});
    auto mass_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*mass_p); auto mass = real2d(mass_raw,{1,1});
    auto rthfddaten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rthfddaten_p); auto rthfddaten = real2d(rthfddaten_raw,{1,1});
    auto rqvfddaten_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rqvfddaten_p); auto rqvfddaten = real2d(rqvfddaten_raw,{1,1});
    // for(int i = 1; i <= nCellsSolve; i++){
    Kokkos::parallel_for("physics_get_tend_phase3_kokkos", Kokkos::RangePolicy<>(1, nCellsSolve + 1), 
        KOKKOS_LAMBDA(const int i) {
        for(int k = 1; k <= nVertLevels; k++){
            tend_th(k,i) = tend_th(k,i) + rthfddaten(k,i)*mass(k,i);
            //tend_qv(k,i) = tend_qv(k,i) + rqvfddaten(k,i)*mass(k,i)
            tend_scalars(index_qv,k,i) = tend_scalars(index_qv,k,i) + rqvfddaten(k,i)*mass(k,i);
        }
    // }
    });
    Kokkos::fence();
}

extern "C" void physics_get_tend_phase4_kokkos(int &nCells, int &nVertLevels, int &nEdges, int &num_scalars, int &nCellsSolve,
    flcl_ndarray_t* tend_th_p,flcl_ndarray_t* tend_scalars_p,
    int &index_qv,flcl_ndarray_t* theta_m_p,flcl_ndarray_t* tend_theta_p,flcl_ndarray_t* scalars_p){
    auto tend_th_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_th_p); auto tend_th = real2d(tend_th_raw,{1,1});
    auto tend_scalars_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*tend_scalars_p); auto tend_scalars = real3d(tend_scalars_raw,{1,1,1});
    auto theta_m_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*theta_m_p); auto theta_m = real2d(theta_m_raw,{1,1});
    auto tend_theta_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tend_theta_p); auto tend_theta = real2d(tend_theta_raw,{1,1});
    auto scalars_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*scalars_p); auto scalars = real3d(scalars_raw,{1,1,1});
    const real R_v = 461.6;
    const real R_d = 287.0;
    // for(int i = 1; i <= nCellsSolve; i++){
    Kokkos::parallel_for("physics_get_tend_phase4_kokkos", Kokkos::RangePolicy<>(1, nCellsSolve + 1), 
        KOKKOS_LAMBDA(const int i) {
        for(int k = 1; k <= nVertLevels; k++){
            real coeff = (1. + R_v/R_d * scalars(index_qv,k,i));
            tend_th(k,i) = coeff * tend_th(k,i) + R_v/R_d * theta_m(k,i) * tend_scalars(index_qv,k,i) / coeff;
            tend_theta(k,i) = tend_theta(k,i) + tend_th(k,i);
        }
    // }
    });
    Kokkos::fence();
}

extern "C" void tend_toEdges_kokkos(int &nCells, int &nVertLevels,int &nCellsSolve, int &nEdges, flcl_ndarray_t* Ux_tend_p, flcl_ndarray_t* Uy_tend_p, flcl_ndarray_t* U_tend_p, flcl_ndarray_t* edgeNormalVectors_p, flcl_ndarray_t* east_p, flcl_ndarray_t* north_p, flcl_ndarray_t* cellsOnEdge_p){
    auto Ux_tend_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*Ux_tend_p); auto Ux_tend = real2d(Ux_tend_raw,{1,1});
    auto Uy_tend_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*Uy_tend_p); auto Uy_tend = real2d(Uy_tend_raw,{1,1});
    auto U_tend_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*U_tend_p); auto U_tend = real2d(U_tend_raw,{1,1});
    auto edgeNormalVectors_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgeNormalVectors_p); auto edgeNormalVectors = real2d(edgeNormalVectors_raw,{1,1});
    auto east_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*east_p); auto east = real2d(east_raw,{1,1});
    auto north_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*north_p); auto north = real2d(north_raw,{1,1});
    auto cellsOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p); auto cellsOnEdge = int2d(cellsOnEdge_raw,{1,1});
    // for(int iEdge = 1; iEdge <= nEdges; iEdge++){
    Kokkos::parallel_for("tend_toEdges_kokkos", Kokkos::RangePolicy<>(1, nEdges + 1), 
        KOKKOS_LAMBDA(const int iEdge) {
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        for(int k = 1; k <= nVertLevels; k++){
            U_tend(k,iEdge) = Ux_tend(k,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   
                                +  edgeNormalVectors(2,iEdge) * east(2,cell1)   
                                +  edgeNormalVectors(3,iEdge) * east(3,cell1))  
                            + Uy_tend(k,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)   
                                +  edgeNormalVectors(2,iEdge) * north(2,cell1)   
                                +  edgeNormalVectors(3,iEdge) * north(3,cell1))  
                            + Ux_tend(k,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   
                                +  edgeNormalVectors(2,iEdge) * east(2,cell2)   
                                +  edgeNormalVectors(3,iEdge) * east(3,cell2))  
                            + Uy_tend(k,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)   
                                +  edgeNormalVectors(2,iEdge) * north(2,cell2)   
                                +  edgeNormalVectors(3,iEdge) * north(3,cell2));
        }
    // }
    });
    Kokkos::fence();
}