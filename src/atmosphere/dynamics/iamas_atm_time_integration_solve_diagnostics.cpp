/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Xinqu develop solve_diagnostics subroutine in kokkos

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

extern "C" void atm_compute_solve_diagnostics_kokkos(int &nVertLevels, int &maxEdges2, int &maxEdges, int &nCells, int &nEdges, int &nVertices, 
        int &vertexDegree, real &dt, real &config_apvm_upwinding, 
        flcl_ndarray_t* fVertex_p, flcl_ndarray_t* fEdge_p, flcl_ndarray_t* invAreaTriangle_p, flcl_ndarray_t* invAreaCell_p, flcl_ndarray_t* dvEdge_p, flcl_ndarray_t* dcEdge_p, flcl_ndarray_t* invDvEdge_p, flcl_ndarray_t* invDcEdge_p, 
        flcl_ndarray_t* weightsOnEdge_p, flcl_ndarray_t* kiteAreasOnVertex_p, flcl_ndarray_t* h_edge_p, flcl_ndarray_t* h_p, flcl_ndarray_t* u_p, flcl_ndarray_t* v_p, flcl_ndarray_t* vorticity_p, flcl_ndarray_t* ke_p, flcl_ndarray_t* pv_edge_p, flcl_ndarray_t* pv_vertex_p, flcl_ndarray_t* pv_cell_p, 
        flcl_ndarray_t* gradPVn_p, flcl_ndarray_t* gradPVt_p, flcl_ndarray_t* divergence_p, flcl_ndarray_t* cellsOnEdge_p, flcl_ndarray_t* cellsOnVertex_p, flcl_ndarray_t* verticesOnEdge_p, flcl_ndarray_t* edgesOnCell_p, flcl_ndarray_t* edgesOnEdge_p, 
        flcl_ndarray_t* edgesOnVertex_p, flcl_ndarray_t* kiteForCell_p, flcl_ndarray_t* verticesOnCell_p, flcl_ndarray_t* edgesOnVertex_sign_p, flcl_ndarray_t* edgesOnCell_sign_p, flcl_ndarray_t* nEdgesOnCell_p, flcl_ndarray_t* nEdgesOnEdge_p, 
        int &cellStart, int &cellEnd, int &vertexStart, int &vertexEnd, int &edgeStart, int &edgeEnd, flcl_ndarray_t* ke_edge_p, flcl_ndarray_t* ke_vertex_p,
        int *rk_step){
    // auto ke_edge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ke_edge_p); auto ke_edge = real2d(ke_edge_raw,{1,1});
    // auto ke_vertex_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ke_vertex_p); auto ke_vertex = real2d(ke_vertex_raw,{1,1});
    auto fVertex_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fVertex_p); auto fVertex = real1d(fVertex_raw,{1});
    auto fEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fEdge_p); auto fEdge = real1d(fEdge_raw,{1});
    auto invAreaTriangle_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaTriangle_p); auto invAreaTriangle = real1d(invAreaTriangle_raw,{1});
    auto invAreaCell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invAreaCell_p); auto invAreaCell = real1d(invAreaCell_raw,{1});
    auto dvEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dvEdge_p); auto dvEdge = real1d(dvEdge_raw,{1});
    auto dcEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*dcEdge_p); auto dcEdge = real1d(dcEdge_raw,{1});
    auto invDvEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invDvEdge_p); auto invDvEdge = real1d(invDvEdge_raw,{1});
    auto invDcEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*invDcEdge_p); auto invDcEdge = real1d(invDcEdge_raw,{1});
    auto weightsOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*weightsOnEdge_p); auto weightsOnEdge = real2d(weightsOnEdge_raw,{1,1});
    auto kiteAreasOnVertex_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kiteAreasOnVertex_p); auto kiteAreasOnVertex = real2d(kiteAreasOnVertex_raw,{1,1});
    auto h_edge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*h_edge_p); auto h_edge = real2d(h_edge_raw,{1,1});
    auto h_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*h_p); auto h = real2d(h_raw,{1,1});
    auto u_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*u_p); auto u = real2d(u_raw,{1,1});
    auto v_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*v_p); auto v = real2d(v_raw,{1,1});
    auto vorticity_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*vorticity_p); auto vorticity = real2d(vorticity_raw,{1,1});
    auto ke_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ke_p); auto ke = real2d(ke_raw,{1,1});
    auto pv_edge_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pv_edge_p); auto pv_edge = real2d(pv_edge_raw,{1,1});
    auto pv_vertex_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pv_vertex_p); auto pv_vertex = real2d(pv_vertex_raw,{1,1});
    auto pv_cell_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*pv_cell_p); auto pv_cell = real2d(pv_cell_raw,{1,1});
    auto gradPVn_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*gradPVn_p); auto gradPVn = real2d(gradPVn_raw,{1,1});
    auto gradPVt_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*gradPVt_p); auto gradPVt = real2d(gradPVt_raw,{1,1});
    auto divergence_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*divergence_p); auto divergence = real2d(divergence_raw,{1,1});
    auto cellsOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnEdge_p); auto cellsOnEdge = int2d(cellsOnEdge_raw,{1,1});
    auto cellsOnVertex_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*cellsOnVertex_p); auto cellsOnVertex = int2d(cellsOnVertex_raw,{1,1});
    auto verticesOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*verticesOnEdge_p); auto verticesOnEdge = int2d(verticesOnEdge_raw,{1,1});
    auto edgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnCell_p); auto edgesOnCell = int2d(edgesOnCell_raw,{1,1});
    auto edgesOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnEdge_p); auto edgesOnEdge = int2d(edgesOnEdge_raw,{1,1});
    auto edgesOnVertex_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*edgesOnVertex_p); auto edgesOnVertex = int2d(edgesOnVertex_raw,{1,1});
    auto kiteForCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*kiteForCell_p); auto kiteForCell = int2d(kiteForCell_raw,{1,1});
    auto verticesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t **>(*verticesOnCell_p); auto verticesOnCell = int2d(verticesOnCell_raw,{1,1});
    auto edgesOnVertex_sign_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnVertex_sign_p); auto edgesOnVertex_sign = real2d(edgesOnVertex_sign_raw,{1,1});
    auto edgesOnCell_sign_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*edgesOnCell_sign_p); auto edgesOnCell_sign = real2d(edgesOnCell_sign_raw,{1,1});
    auto nEdgesOnCell_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnCell_p); auto nEdgesOnCell = int1d(nEdgesOnCell_raw,{1});
    auto nEdgesOnEdge_raw = flcl::view_from_ndarray<flcl::flcl_view_i32_c_t *>(*nEdgesOnEdge_p); auto nEdgesOnEdge = int1d(nEdgesOnEdge_raw,{1});
    const bool hollingsworth = true;
    real ke_edge_s[nVertLevels * (nEdges + 1)];
    real ke_vertex_s[nVertLevels * (nVertices + 1)];
    auto ke_edge = real2d((real*)ke_edge_s,{1,1},{nVertLevels+1,nEdges+2});
    auto ke_vertex = real2d((real*)ke_vertex_s,{1,1},{nVertLevels+1,nVertices+2});
    for(int k = 1; k <= nVertLevels; ++k){
        ke_edge(k,nEdges+1) = 0.0;
        ke_vertex(k,nVertices+1) = 0.0;
    }
    Kokkos::fence();
    // Compute height on cell edges at velocity locations
    Kokkos::parallel_for(
        "ComputeKE_Initial", Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), KOKKOS_LAMBDA(int iEdge) {
    // for(int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge){
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        for(int k = 1; k <= nVertLevels; ++k){
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2));
        }
        // the first openmp barrier below is set so that ke_edge is computed
        // it would be good to move this somewhere else?
        real efac = dcEdge(iEdge)*dvEdge(iEdge);

        for(int k = 1; k <= nVertLevels; ++k){
            ke_edge(k,iEdge) = efac*u(k,iEdge)*u(k,iEdge);
        }
    // }
    });
    Kokkos::fence();
    // Compute circulation and relative vorticity at each vertex
    Kokkos::parallel_for(
        "Compute_vorticity", Kokkos::RangePolicy<>(vertexStart, vertexEnd + 1), KOKKOS_LAMBDA(int iVertex) {
    // for(int iVertex = vertexStart; iVertex <= vertexEnd; ++iVertex){
        for(int k = 1; k <= nVertLevels; ++k){
            vorticity(k,iVertex) = 0.0;
        }
        for(int i = 1; i <= vertexDegree; i++){
            int iEdge = edgesOnVertex(i,iVertex);
            real s = edgesOnVertex_sign(i,iVertex) * dcEdge(iEdge);
            for(int k = 1; k <= nVertLevels; ++k){
                vorticity(k,iVertex) = vorticity(k,iVertex) + s * u(k,iEdge);
            }
        }
        for(int k = 1; k <= nVertLevels; ++k){
            vorticity(k,iVertex) = vorticity(k,iVertex) * invAreaTriangle(iVertex);
        }
    // }
    });
    Kokkos::fence();
    // Compute the divergence at each cell center
    // Kokkos::parallel_for(
    //     "Compute_divergence", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
    Kokkos::parallel_for(
        "Compute_divergence", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellStart; iCell <= cellEnd; ++iCell){
        for(int k = 1; k <= nVertLevels; ++k){
            divergence(k,iCell) = 0.0;
        }
        for(int i = 1; i <= nEdgesOnCell(iCell); i++){
            int iEdge = edgesOnCell(i,iCell);
            real s = edgesOnCell_sign(i,iCell) * dvEdge(iEdge);
            for(int k = 1; k <= nVertLevels; ++k){
                divergence(k,iCell) = divergence(k,iCell) + s * u(k,iEdge);
            }
        }
        real r = invAreaCell(iCell);
        for(int k = 1; k <= nVertLevels; ++k){
            divergence(k,iCell) = divergence(k,iCell) * r;
        }
    // }
    });
    Kokkos::fence();
    // Compute kinetic energy in each cell (Ringler et al JCP 2009)
    // Replace 2.0 with 2 in exponentiation to avoid outside chance that
    // compiler will actually allow "float raised to float" operation

    Kokkos::parallel_for(
        "Compute_ke", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
    // for(int iCell = cellStart; iCell <= cellEnd; ++iCell){
        for(int k = 1; k <= nVertLevels; ++k){
            ke(k,iCell) = 0.0;
        }
        for(int i = 1; i <= nEdgesOnCell(iCell); i++){
            int iEdge = edgesOnCell(i,iCell);
            for(int k = 1; k <= nVertLevels; ++k){
                ke(k,iCell) = ke(k,iCell) + 0.25 * ke_edge(k,iEdge);
            }
        }
        for(int k = 1; k <= nVertLevels; ++k){
            ke(k,iCell) = ke(k,iCell) * invAreaCell(iCell);
        }
    // }
    });
    Kokkos::fence();
    if (hollingsworth){
        // Compute ke at cell vertices - AG's new KE construction, part 1
        // *** approximation here because we don't have inner triangle areas

        // Replace 2.0 with 2 in exponentiation to avoid outside chance that
        // compiler will actually allow "float raised to float" operation

        Kokkos::parallel_for(
            "ComputeKE_Standard", Kokkos::RangePolicy<>(vertexStart, vertexEnd + 1), KOKKOS_LAMBDA(int iVertex) {
        // for(int iVertex = vertexStart; iVertex <= vertexEnd; ++iVertex){
            real r = 0.25 * invAreaTriangle(iVertex);
            for(int k = 1; k <= nVertLevels; ++k){
                ke_vertex(k,iVertex) = (  ke_edge(k,edgesOnVertex(1,iVertex))+ke_edge(k,edgesOnVertex(2,iVertex))+ke_edge(k,edgesOnVertex(3,iVertex)) )*r;
            }
        // }
        });
        Kokkos::fence();
        // adjust ke at cell vertices - AG's new KE construction, part 2
        real ke_fact = 1.0 - .375;
        Kokkos::parallel_for(
            "ComputeKE_Hollingsworth", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellStart; iCell <= cellEnd; ++iCell){
            for(int k = 1; k <= nVertLevels; ++k){
                ke(k,iCell) = ke_fact * ke(k,iCell);
            }
        // }
        });
        Kokkos::fence();
        Kokkos::parallel_for(
            "ComputeOnCells", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellStart; iCell <= cellEnd; ++iCell){
            real r = invAreaCell(iCell);
            for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                int iVertex = verticesOnCell(i,iCell);
                int j = kiteForCell(i,iCell);
                for(int k = 1; k <= nVertLevels; ++k){
                    ke(k,iCell) = ke(k,iCell) + (1.-ke_fact)*kiteAreasOnVertex(j,iVertex) * ke_vertex(k,iVertex) * r;
                }
            }
        // }
        });
    }
    Kokkos::fence();
    // Compute v (tangential) velocities following Thuburn et al JCP 2009
    // The tangential velocity is only used to compute the Smagorinsky coefficient
    bool reconstruct_v = true;
    if(rk_step != nullptr) {  
        if(*rk_step != 3) { 
            reconstruct_v = false;
        }
    }
    if(reconstruct_v){
        Kokkos::parallel_for(
            "Compute_v", Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), KOKKOS_LAMBDA(int iEdge) {
        // for(int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge){
            for(int k = 1; k <= nVertLevels; ++k){
                v(k,iEdge) = 0.0;
            }
            for(int i = 1; i <= nEdgesOnEdge(iEdge); i++){
                int eoe = edgesOnEdge(i,iEdge);
                for(int k = 1; k <= nVertLevels; ++k){
                    v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe);
                }
            }
        // }
        });
    }
    Kokkos::fence();
    /**
     * Compute height at vertices, pv at vertices, and average pv to edge locations
     *  ( this computes pv_vertex at all vertices bounding real cells )
     * Avoid dividing h_vertex by areaTriangle and move areaTriangle into
     * numerator for the pv_vertex calculation
     */
    Kokkos::parallel_for(
        "Compute_pv_vertex", Kokkos::RangePolicy<>(vertexStart, vertexEnd + 1), KOKKOS_LAMBDA(int iVertex) {
    // for(int iVertex = vertexStart; iVertex <= vertexEnd; ++iVertex){
        for(int k = 1; k <= nVertLevels; ++k){
            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex));
        }
    // }
    });
    Kokkos::fence();
    // Compute pv at the edges
    //   ( this computes pv_edge at all edges bounding real cells )
    Kokkos::parallel_for(
        "Compute_pv_edge", Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), KOKKOS_LAMBDA(int iEdge) {
    // for(int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge){
        for(int k = 1; k <= nVertLevels; ++k){
            pv_edge(k,iEdge) =  0.5 * (pv_vertex(k,verticesOnEdge(1,iEdge)) + pv_vertex(k,verticesOnEdge(2,iEdge)));
        }
    // }
    });
    Kokkos::fence();
    if (config_apvm_upwinding > 0.0){
        /**
         * Compute pv at cell centers
         *  ( this computes pv_cell for all real cells )
         * only needed for APVM upwinding
         */
        Kokkos::parallel_for(
            "Compute_pv_cell", Kokkos::RangePolicy<>(cellStart, cellEnd + 1), KOKKOS_LAMBDA(int iCell) {
        // for(int iCell = cellStart; iCell <= cellEnd; ++iCell){
            for(int k = 1; k <= nVertLevels; ++k){
                pv_cell(k,iCell) = 0.0;
            }
            real r = invAreaCell(iCell);
            for(int i = 1; i <= nEdgesOnCell(iCell); i++){
                int iVertex = verticesOnCell(i,iCell);
                int j = kiteForCell(i,iCell);
                for(int k = 1; k <= nVertLevels; ++k){
                    pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(j,iVertex) * pv_vertex(k,iVertex) * r;
                }
            }
        // }
        });
        Kokkos::fence();
        /**
         * Modify PV edge with upstream bias. 
         * Compute gradient of PV in the tangent direction
         *  ( this computes gradPVt at all edges bounding real cells )
         * Compute gradient of PV in normal direction
         * (tdr: 2009-10-02: this is not correct because the pv_cell in the halo is not correct)
         * Modify PV edge with upstream bias.
         * Merged loops for calculating gradPVt, gradPVn and pv_edge
         * Also precomputed inverses of dvEdge and dcEdge to avoid repeated divisions
         */
        real r = config_apvm_upwinding * dt;
        Kokkos::parallel_for(
            "Compute_pv_upstream", Kokkos::RangePolicy<>(edgeStart, edgeEnd + 1), KOKKOS_LAMBDA(int iEdge) {
        // for(int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge){
            real r1 = 1.0 * invDvEdge(iEdge);
            real r2 = 1.0 * invDcEdge(iEdge);
            for(int k = 1; k <= nVertLevels; ++k){
                gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) * r1;
                gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) * r2;
                pv_edge(k,iEdge) = pv_edge(k,iEdge) - r * (v(k,iEdge) * gradPVt(k,iEdge) + u(k,iEdge) * gradPVn(k,iEdge));
            }
        // }
        });
    }
    Kokkos::fence();
}
