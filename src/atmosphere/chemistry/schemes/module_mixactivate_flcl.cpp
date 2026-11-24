// =============================================================================

#include <math.h>
#include <algorithm>
#include <Kokkos_Core.hpp>
#include <flcl-cxx.hpp>
#include <utility>

typedef double real;

typedef Kokkos::View<bool**, Kokkos::LayoutLeft> bool2d;

typedef Kokkos::View<bool***, Kokkos::LayoutLeft> bool3d;

typedef Kokkos::View<int*, Kokkos::LayoutLeft> int1d;

typedef Kokkos::View<int**, Kokkos::LayoutLeft> int2d;

typedef Kokkos::View<int***, Kokkos::LayoutLeft> int3d;

typedef Kokkos::View<int****, Kokkos::LayoutLeft> int4d;

typedef Kokkos::View<real*, Kokkos::LayoutLeft> real1d;

typedef Kokkos::View<real**, Kokkos::LayoutLeft> real2d;

typedef Kokkos::View<real***, Kokkos::LayoutLeft> real3d;

typedef Kokkos::View<real****, Kokkos::LayoutLeft> real4d;

typedef Kokkos::View<real*****, Kokkos::LayoutLeft> real5d;

constexpr int CHEM_WORKSPACE_CHUNK_SIZE = 64;

//module_atmchem_constants
constexpr real g = 9.80616;
constexpr real r_v = 461.6;
constexpr real rhowater = 1000.;
constexpr real xlv = 2.5e6;
constexpr real mwdry   = 28.966;
//mpas_constants
constexpr real rgas    = 287.0;
constexpr real rvord   = r_v/rgas;
constexpr real cp   = 7.*rgas/2.;
constexpr real r_d    = 287.0;
constexpr real ep_2    = r_d/r_v;

template <class ChemView, class NpvView, class MaerView, class MaerTotView, class Aer2View>
KOKKOS_INLINE_FUNCTION void loadaer(const ChemView &chem,int k,real cs,const NpvView &npv,
                          const real2d &dlo_sect,const real2d &dhi_sect,int maxd_acomp,const int1d &ncomp,
                          const int3d &numptr_aer, const real2d &dens_aer,
                          const int4d &massptr_aer,
                          const MaerView &maerosol, const MaerView &maerosolcw,
                          const MaerTotView &maerosol_tot, const MaerTotView &maerosol_totcw,
                          const Aer2View &naerosol,const Aer2View &naerosolcw,
                          const Aer2View &vaerosol,const Aer2View &vaerosolcw,
                          int m,int n,int ai_phase,int cw_phase,int i);

template <class Aer2View, class FluxFullView, class WorkReal2View,
          class WorkReal1View, class WorkBool2View, class WorkInt1View>
KOKKOS_INLINE_FUNCTION void activate(real wbar, real sigw, real wdiab, real wminf, real wmaxf, real tair, real rhoair,
                      int msectional, int maxd_atype, int ntype_aer, int maxd_asize, const int1d &nsize_aer,
                      const Aer2View &na, const Aer2View &volc,const real2d &dlo_sect,const real2d &dhi_sect,const real2d &sigman, const Aer2View &hygro,
                      const Aer2View &fn, const Aer2View &fs, const Aer2View &fm, const Aer2View &fluxn, const Aer2View &fluxs, const Aer2View &fluxm, const FluxFullView &flux_fullact,
                      const WorkReal2View &ylo, const WorkReal2View &yhi, const WorkReal2View &ymean,
                      const WorkReal2View &alnsign, const WorkReal2View &am, const WorkReal2View &lnhygro,
                      const WorkReal2View &f1, const WorkReal2View &sm, const WorkReal2View &lnsm,
                      const WorkReal2View &eta, const WorkReal2View &asub, const WorkReal2View &bsub,
                      const WorkReal2View &sumflxn, const WorkReal2View &sumflxs, const WorkReal2View &sumflxm,
                      const WorkReal2View &sumfn, const WorkReal2View &sumfs, const WorkReal2View &sumfm,
                      const WorkReal2View &fnold, const WorkReal2View &fsold, const WorkReal2View &fmold,
                      const WorkReal2View &rlo, const WorkReal2View &rhi, const WorkReal2View &rmean,
                      const WorkReal2View &etafactor2,
                      const WorkReal1View &totn, const WorkReal1View &gmrad,
                      const WorkReal1View &gmradsq, const WorkReal1View &gmlnsig,
                      const WorkReal1View &gmsm, const WorkReal1View &sumns,
                      const WorkBool2View &bin_is_empty, const WorkBool2View &bin_is_narrow,
                      const WorkInt1View &tmp_nsize_aer,
                      int grid_id, int ktau, int ii,  int kk,real smax_prescribed);

KOKKOS_INLINE_FUNCTION real erf_alt(real &x);
KOKKOS_INLINE_FUNCTION real erfc_num_recipes(real &x);

template <class Real2View>
KOKKOS_INLINE_FUNCTION void maxsat(real &zeta, const Real2View &eta,
                        int &maxd_atype,int &ntype_aer,int &maxd_asize,const int1d &nsize_aer,
                        const Real2View &sm, const Real2View &alnsign,
                        const Real2View &f1,real &smax,int i);

template <class Real2View, class Real1View, class Int1View>
KOKKOS_INLINE_FUNCTION void maxsat_1(real &zeta, const Real2View &eta,
                        int &maxd_atype,int &ntype_aer,int &maxd_asize,const Int1View &nsize_aer,
                        const Real1View &sm, const Real1View &gmlnsig,
                        const Real2View &f1,real &smax,int i);

template <class Real2View>
KOKKOS_INLINE_FUNCTION void maxsat_init(int &maxd_atype, int &ntype_aer,
           int &maxd_asize, const int1d &nsize_aer,
           const Real2View &alnsign, const Real2View &f1,int i);

template <class Real2View, class Real1View, class Int1View>
KOKKOS_INLINE_FUNCTION void maxsat_init_1(int &maxd_atype, int &ntype_aer,
           int &maxd_asize, const Int1View &nsize_aer,
           const Real1View &gmlnsig, const Real2View &f1,int i);

template <class QView, class SrcView, class EkkView, class OverlapView, class QactView>
KOKKOS_INLINE_FUNCTION void explmix_t(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,
                       const QactView &qactold,int lnum,int lnumcw,int nnew,int nsav,int i);

template <class QView, class SrcView, class EkkView, class OverlapView>
KOKKOS_INLINE_FUNCTION void explmix_t_3d(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,
                       const QView &qactold,int lnum,int lnumcw,int nnew,int nsav,int i);

template <class QView, class SrcView, class EkkView, class OverlapView>
KOKKOS_INLINE_FUNCTION void explmix_f(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,int i);

template <class QView, class SrcView, class EkkView, class OverlapView>
KOKKOS_INLINE_FUNCTION void explmix_f_3d(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,int lnumcw,int nnew,int nsav,int i);

struct mixactivate_workspace {
    real2d qndrop, lcldfra, lcldfra_old, wtke, zn, zs;
    real2d cs, csbot, csbot_cscen, dz, qndrop_new, ekd, ekk, srcn;
    real2d surfrate, overlapp, overlapm, ekkp, ekkm, source, flux_fullact;
    real3d naerosol, naerosolcw, maerosol_tot, maerosol_totcw, vaerosol, vaerosolcw;
    real4d maerosol, maerosolcw, raercol;
    real3d fn, fs, fm, fluxn, fluxs, fluxm;
    real3d ylo, yhi, ymean, alnsign, am, lnhygro, f1, sm_activate, lnsm;
    real3d eta, asub, bsub, sumflxn, sumflxs, sumflxm;
    real3d sumfn, sumfs, sumfm, fnold, fsold, fmold, rlo, rhi, rmean, etafactor2;
    real2d totn, gmrad, gmradsq, gmlnsig, gmsm, sumns;
    bool3d bin_is_empty, bin_is_narrow;
    int2d tmp_nsize_aer;
    real4d nact, mact, ccnfact;
    real3d npv, hygro_aer, ccn, amcube, argfactor;
    real2d rhodz;

    void init(int kte, int ncol, int num_chem,
              int maxd_acomp, int maxd_asize, int maxd_atype) {
        qndrop = real2d("mix_ws_qndrop", kte+1, ncol);
        lcldfra = real2d("mix_ws_lcldfra", kte, ncol);
        lcldfra_old = real2d("mix_ws_lcldfra_old", kte, ncol);
        wtke = real2d("mix_ws_wtke", kte, ncol);
        zn = real2d("mix_ws_zn", kte+1, ncol);
        zs = real2d("mix_ws_zs", kte+1, ncol);
        cs = real2d("mix_ws_cs", kte, ncol);
        csbot = real2d("mix_ws_csbot", kte+1, ncol);
        csbot_cscen = real2d("mix_ws_csbot_cscen", kte+1, ncol);
        dz = real2d("mix_ws_dz", kte, ncol);
        qndrop_new = real2d("mix_ws_qndrop_new", kte+1, ncol);
        ekd = real2d("mix_ws_ekd", kte, ncol);
        ekk = real2d("mix_ws_ekk", kte+1, ncol);
        srcn = real2d("mix_ws_srcn", kte, ncol);
        surfrate = real2d("mix_ws_surfrate", num_chem, ncol);
        overlapp = real2d("mix_ws_overlapp", kte, ncol);
        overlapm = real2d("mix_ws_overlapm", kte, ncol);
        ekkp = real2d("mix_ws_ekkp", kte+1, ncol);
        ekkm = real2d("mix_ws_ekkm", kte, ncol);
        source = real2d("mix_ws_source", kte+1, ncol);
        flux_fullact = real2d("mix_ws_flux_fullact", kte, ncol);
        naerosol = real3d("mix_ws_naerosol", maxd_asize, maxd_atype, ncol);
        naerosolcw = real3d("mix_ws_naerosolcw", maxd_asize, maxd_atype, ncol);
        maerosol_tot = real3d("mix_ws_maerosol_tot", maxd_asize, maxd_atype, ncol);
        maerosol_totcw = real3d("mix_ws_maerosol_totcw", maxd_asize, maxd_atype, ncol);
        vaerosol = real3d("mix_ws_vaerosol", maxd_asize, maxd_atype, ncol);
        vaerosolcw = real3d("mix_ws_vaerosolcw", maxd_asize, maxd_atype, ncol);
        maerosol = real4d("mix_ws_maerosol", maxd_acomp, maxd_asize, maxd_atype, ncol);
        maerosolcw = real4d("mix_ws_maerosolcw", maxd_acomp, maxd_asize, maxd_atype, ncol);
        raercol = real4d("mix_ws_raercol", kte+1, num_chem, 2, ncol);
        fn = real3d("mix_ws_fn", maxd_asize, maxd_atype, ncol);
        fs = real3d("mix_ws_fs", maxd_asize, maxd_atype, ncol);
        fm = real3d("mix_ws_fm", maxd_asize, maxd_atype, ncol);
        fluxn = real3d("mix_ws_fluxn", maxd_asize, maxd_atype, ncol);
        fluxs = real3d("mix_ws_fluxs", maxd_asize, maxd_atype, ncol);
        fluxm = real3d("mix_ws_fluxm", maxd_asize, maxd_atype, ncol);
        ylo = real3d("mix_ws_activate_ylo", maxd_asize, maxd_atype, ncol);
        yhi = real3d("mix_ws_activate_yhi", maxd_asize, maxd_atype, ncol);
        ymean = real3d("mix_ws_activate_ymean", maxd_asize, maxd_atype, ncol);
        alnsign = real3d("mix_ws_activate_alnsign", maxd_asize, maxd_atype, ncol);
        am = real3d("mix_ws_activate_am", maxd_asize, maxd_atype, ncol);
        lnhygro = real3d("mix_ws_activate_lnhygro", maxd_asize, maxd_atype, ncol);
        f1 = real3d("mix_ws_activate_f1", maxd_asize, maxd_atype, ncol);
        sm_activate = real3d("mix_ws_activate_sm", maxd_asize, maxd_atype, ncol);
        lnsm = real3d("mix_ws_activate_lnsm", maxd_asize, maxd_atype, ncol);
        eta = real3d("mix_ws_activate_eta", maxd_asize, maxd_atype, ncol);
        asub = real3d("mix_ws_activate_asub", maxd_asize, maxd_atype, ncol);
        bsub = real3d("mix_ws_activate_bsub", maxd_asize, maxd_atype, ncol);
        sumflxn = real3d("mix_ws_activate_sumflxn", maxd_asize, maxd_atype, ncol);
        sumflxs = real3d("mix_ws_activate_sumflxs", maxd_asize, maxd_atype, ncol);
        sumflxm = real3d("mix_ws_activate_sumflxm", maxd_asize, maxd_atype, ncol);
        sumfn = real3d("mix_ws_activate_sumfn", maxd_asize, maxd_atype, ncol);
        sumfs = real3d("mix_ws_activate_sumfs", maxd_asize, maxd_atype, ncol);
        sumfm = real3d("mix_ws_activate_sumfm", maxd_asize, maxd_atype, ncol);
        fnold = real3d("mix_ws_activate_fnold", maxd_asize, maxd_atype, ncol);
        fsold = real3d("mix_ws_activate_fsold", maxd_asize, maxd_atype, ncol);
        fmold = real3d("mix_ws_activate_fmold", maxd_asize, maxd_atype, ncol);
        rlo = real3d("mix_ws_activate_rlo", maxd_asize, maxd_atype, ncol);
        rhi = real3d("mix_ws_activate_rhi", maxd_asize, maxd_atype, ncol);
        rmean = real3d("mix_ws_activate_rmean", maxd_asize, maxd_atype, ncol);
        etafactor2 = real3d("mix_ws_activate_etafactor2", maxd_asize, maxd_atype, ncol);
        totn = real2d("mix_ws_activate_totn", maxd_atype, ncol);
        gmrad = real2d("mix_ws_activate_gmrad", maxd_atype, ncol);
        gmradsq = real2d("mix_ws_activate_gmradsq", maxd_atype, ncol);
        gmlnsig = real2d("mix_ws_activate_gmlnsig", maxd_atype, ncol);
        gmsm = real2d("mix_ws_activate_gmsm", maxd_atype, ncol);
        sumns = real2d("mix_ws_activate_sumns", maxd_atype, ncol);
        bin_is_empty = bool3d("mix_ws_activate_bin_is_empty", maxd_asize, maxd_atype, ncol);
        bin_is_narrow = bool3d("mix_ws_activate_bin_is_narrow", maxd_asize, maxd_atype, ncol);
        tmp_nsize_aer = int2d("mix_ws_activate_tmp_nsize_aer", maxd_atype, ncol);
        nact = real4d("mix_ws_nact", kte, maxd_asize, maxd_atype, ncol);
        mact = real4d("mix_ws_mact", kte, maxd_asize, maxd_atype, ncol);
        npv = real3d("mix_ws_npv", maxd_asize, maxd_atype, ncol);
        hygro_aer = real3d("mix_ws_hygro_aer", maxd_asize, maxd_atype, ncol);
        ccn = real3d("mix_ws_ccn", 6, kte, ncol);
        ccnfact = real4d("mix_ws_ccnfact", 6, maxd_asize, maxd_atype, ncol);
        amcube = real3d("mix_ws_amcube", maxd_asize, maxd_atype, ncol);
        argfactor = real3d("mix_ws_argfactor", maxd_asize, maxd_atype, ncol);
        rhodz = real2d("mix_ws_rhodz", kte, ncol);
    }
};

static mixactivate_workspace *mixactivate_ws = nullptr;
static bool is_mixactivate_ws_inited = false;
static int mixactivate_ws_kte = 0, mixactivate_ws_ncol = 0, mixactivate_ws_num_chem = 0;
static int mixactivate_ws_maxd_acomp = 0, mixactivate_ws_maxd_asize = 0;
static int mixactivate_ws_maxd_atype = 0;

static bool mixactivate_workspace_matches(int kte, int ncol, int num_chem,
                                          int maxd_acomp, int maxd_asize,
                                          int maxd_atype) {
    return is_mixactivate_ws_inited && mixactivate_ws != nullptr &&
           mixactivate_ws_kte == kte && mixactivate_ws_ncol == ncol &&
           mixactivate_ws_num_chem == num_chem &&
           mixactivate_ws_maxd_acomp == maxd_acomp &&
           mixactivate_ws_maxd_asize == maxd_asize &&
           mixactivate_ws_maxd_atype == maxd_atype;
}

static mixactivate_workspace &get_mixactivate_workspace(
        int kte, int ncol, int num_chem,
        int maxd_acomp, int maxd_asize, int maxd_atype) {
    if (!mixactivate_workspace_matches(kte, ncol, num_chem,
                                      maxd_acomp, maxd_asize, maxd_atype)) {
        if (mixactivate_ws != nullptr) {
            Kokkos::fence();
            delete mixactivate_ws;
            mixactivate_ws = nullptr;
        }
        mixactivate_ws = new mixactivate_workspace();
        mixactivate_ws->init(kte, ncol, num_chem,
                             maxd_acomp, maxd_asize, maxd_atype);
        mixactivate_ws_kte = kte;
        mixactivate_ws_ncol = ncol;
        mixactivate_ws_num_chem = num_chem;
        mixactivate_ws_maxd_acomp = maxd_acomp;
        mixactivate_ws_maxd_asize = maxd_asize;
        mixactivate_ws_maxd_atype = maxd_atype;
        is_mixactivate_ws_inited = true;
    }
    return *mixactivate_ws;
}

extern "C" void mixactivate_workspace_finalize() {
    if (mixactivate_ws != nullptr) {
        if (Kokkos::is_initialized()) {
            Kokkos::fence();
        }
        delete mixactivate_ws;
        mixactivate_ws = nullptr;
    }
    is_mixactivate_ws_inited = false;
    mixactivate_ws_kte = 0;
    mixactivate_ws_ncol = 0;
    mixactivate_ws_num_chem = 0;
    mixactivate_ws_maxd_acomp = 0;
    mixactivate_ws_maxd_asize = 0;
    mixactivate_ws_maxd_atype = 0;
}

extern "C" void c_mixactivate(int& msectional,
           flcl_ndarray_t* chem_p,int& num_chem,flcl_ndarray_t* qv_p,flcl_ndarray_t* qc_p,flcl_ndarray_t* qi_p,flcl_ndarray_t* qs_p,flcl_ndarray_t* qndrop3d_p,flcl_ndarray_t* tend_qndrop_activation_p,
           flcl_ndarray_t* temp_p,flcl_ndarray_t* p_phy_p,flcl_ndarray_t* w_p,flcl_ndarray_t* ddvel_p,int& idrydep_onoff,
           int& maxd_acomp,int& maxd_asize,int& maxd_atype,int& maxd_aphase,
           flcl_ndarray_t* ncomp_aer_p,flcl_ndarray_t* nsize_aer_p,int& ntype_aer,int& nphase_aer,
           flcl_ndarray_t* numptr_aer_p,flcl_ndarray_t* massptr_aer_p,flcl_ndarray_t* dlo_sect_p,flcl_ndarray_t* dhi_sect_p,flcl_ndarray_t* sigmag_aer_p,flcl_ndarray_t* dpvolmean_aer_p,
           flcl_ndarray_t* dens_aer_p,flcl_ndarray_t* mw_aer_p,
           flcl_ndarray_t* waterptr_aer_p,flcl_ndarray_t* hygro_p,int& ai_phase,int& cw_phase,
           int& its,int& ite,int& kts,int& kte,
           flcl_ndarray_t* rho_p,flcl_ndarray_t* zm_p,flcl_ndarray_t* dz8w_p,flcl_ndarray_t* p_at_w_p,flcl_ndarray_t* t_at_w_p,flcl_ndarray_t* kvh_p,
           flcl_ndarray_t* cldfra_p,flcl_ndarray_t* cldfra_old_p,flcl_ndarray_t* qsrflx_p,
           flcl_ndarray_t* ccn1_p,flcl_ndarray_t* ccn2_p,flcl_ndarray_t* ccn3_p,flcl_ndarray_t* ccn4_p,flcl_ndarray_t* ccn5_p,flcl_ndarray_t* ccn6_p,flcl_ndarray_t* nsource_p,
           int& grid_id,int& ktau,real& dtstep,
           bool& f_qv,bool& f_qc,bool& f_qi,bool& f_qs,flcl_ndarray_t* f_ice_phy_p,flcl_ndarray_t* f_rain_phy_p,int& sw_physics){
        //receive
        //maxd_asize = 10
        //maxd_atype = 1

        //input variables
        auto chem_all                   = flcl::view_from_ndarray<real***>(*chem_p);
        auto qv_all                     = flcl::view_from_ndarray<real**>(*qv_p);
        auto qc_all                     = flcl::view_from_ndarray<real**>(*qc_p);
        auto qi_all                     = flcl::view_from_ndarray<real**>(*qi_p);
        auto qs_all                     = flcl::view_from_ndarray<real**>(*qs_p);
        auto qndrop3d_all               = flcl::view_from_ndarray<real**>(*qndrop3d_p);
        auto tend_qndrop_activation_all = flcl::view_from_ndarray<real**>(*tend_qndrop_activation_p);

        auto temp_all                   = flcl::view_from_ndarray<real**>(*temp_p);
        auto p_phy_all                  = flcl::view_from_ndarray<real**>(*p_phy_p);
        auto w_all                      = flcl::view_from_ndarray<real**>(*w_p);
        auto ddvel_all                  = flcl::view_from_ndarray<real**>(*ddvel_p);

        auto ncomp_aer              = flcl::view_from_ndarray<int*>(*ncomp_aer_p);
        auto nsize_aer              = flcl::view_from_ndarray<int*>(*nsize_aer_p);

        auto numptr_aer             = flcl::view_from_ndarray<int***>(*numptr_aer_p);
        auto massptr_aer            = flcl::view_from_ndarray<int****>(*massptr_aer_p);
        auto dlo_sect               = flcl::view_from_ndarray<real**>(*dlo_sect_p);
        auto dhi_sect               = flcl::view_from_ndarray<real**>(*dhi_sect_p);
        auto sigmag_aer             = flcl::view_from_ndarray<real**>(*sigmag_aer_p);
        auto dpvolmean_aer          = flcl::view_from_ndarray<real**>(*dpvolmean_aer_p);

        auto dens_aer               = flcl::view_from_ndarray<real**>(*dens_aer_p);
        auto mw_aer                 = flcl::view_from_ndarray<real**>(*mw_aer_p);

        auto waterptr_aer           = flcl::view_from_ndarray<int**>(*waterptr_aer_p);
        auto hygro_all              = flcl::view_from_ndarray<real****>(*hygro_p);

        auto rho_all                    = flcl::view_from_ndarray<real**>(*rho_p);
        auto zm_all                     = flcl::view_from_ndarray<real**>(*zm_p);
        auto dz8w_all                   = flcl::view_from_ndarray<real**>(*dz8w_p);
        auto p_at_w_all                 = flcl::view_from_ndarray<real**>(*p_at_w_p);
        auto t_at_w_all                 = flcl::view_from_ndarray<real**>(*t_at_w_p);
        auto kvh_all                    = flcl::view_from_ndarray<real**>(*kvh_p);

        auto cldfra_all                 = flcl::view_from_ndarray<real**>(*cldfra_p);
        auto cldfra_old_all             = flcl::view_from_ndarray<real**>(*cldfra_old_p);
        auto qsrflx_all                 = flcl::view_from_ndarray<real**>(*qsrflx_p);

        auto ccn1_all                   = flcl::view_from_ndarray<real**>(*ccn1_p);
        auto ccn2_all                   = flcl::view_from_ndarray<real**>(*ccn2_p);
        auto ccn3_all                   = flcl::view_from_ndarray<real**>(*ccn3_p);
        auto ccn4_all                   = flcl::view_from_ndarray<real**>(*ccn4_p);
        auto ccn5_all                   = flcl::view_from_ndarray<real**>(*ccn5_p);
        auto ccn6_all                   = flcl::view_from_ndarray<real**>(*ccn6_p);
        auto nsource_all                = flcl::view_from_ndarray<real**>(*nsource_p);

        auto f_ice_phy_all              = flcl::view_from_ndarray<real**>(*f_ice_phy_p);
        auto f_rain_phy_all             = flcl::view_from_ndarray<real**>(*f_rain_phy_p);

        //local variables
        const int ncol = ite - its + 1;
#ifdef USE_SERIALIZATION
        const int workspace_ncol = 1;
#else
        const int workspace_ncol = CHEM_WORKSPACE_CHUNK_SIZE;
#endif
        auto &workspace = get_mixactivate_workspace(
            kte, workspace_ncol, num_chem,
            maxd_acomp, maxd_asize, maxd_atype);

        auto qndrop      = workspace.qndrop;
        auto lcldfra     = workspace.lcldfra;
        auto lcldfra_old = workspace.lcldfra_old;
        auto wtke        = workspace.wtke;
        auto zn          = workspace.zn;
        auto zs          = workspace.zs;

        constexpr real zkmin = 0.01;
        constexpr real zkmax = 100.0;

        auto cs             = workspace.cs;
        auto csbot          = workspace.csbot;
        auto csbot_cscen    = workspace.csbot_cscen;
        auto dz             = workspace.dz;

        real wdiab;
        constexpr real wmixmin = 0.2;

        auto qndrop_new     = workspace.qndrop_new;
        auto ekd            = workspace.ekd;
        auto ekk            = workspace.ekk;
        auto srcn           = workspace.srcn;

        constexpr real sq2pi = 2.5066282746;
        real dtinv;

        int km1,kp1;
        real wbar,wmix,wmin,wmax;
        real dum;
        real tmpa, tmpb, tmpc, tmpc1, tmpc2, tmpd, tmpe, tmpf;
        real tmpcourno;
        real dact;
        real fluxntot;
        real fac_srflx;
        real depvel_drop, depvel_tmp;
        constexpr real depvel_uplimit = 1.0;                      //! upper limit for dep vels (m/s)

        auto surfrate = workspace.surfrate;                     //! surface exchange rate (/s)
        real surfratemax;                                         //! max surfrate for all species treated here
        real surfrate_drop;                                       //! surfade exchange rate for droplelts
        real dtmin,tinv,dtt;
        int nsubmix,nsubmix_bnd;
        int i,k,m,n,nsub;
        real dtmix;
        real alogarg;
        real pi;
        int nnew,nsav,ntemp;
        auto overlapp        = workspace.overlapp;              //! cloud overlap
        auto overlapm        = workspace.overlapm;              //! cloud overlap
        auto ekkp            = workspace.ekkp;                  //! zn*zs*density*diffusivity
        auto ekkm            = workspace.ekkm;                  //! zn*zs*density*diffusivity

        int lnum,lnumcw,l,lmass,lmasscw,lsfc,lsfccw,ltype,lsig,lwater;
        auto naerosol        = workspace.naerosol;
        auto naerosolcw      = workspace.naerosolcw;
        auto maerosol        = workspace.maerosol;
        auto maerosolcw      = workspace.maerosolcw;
        auto maerosol_tot    = workspace.maerosol_tot;
        auto maerosol_totcw  = workspace.maerosol_totcw;
        auto vaerosol        = workspace.vaerosol;
        auto vaerosolcw      = workspace.vaerosolcw;
        auto raercol         = workspace.raercol;
        auto source          = workspace.source;

        auto fn              = workspace.fn;
        auto fs              = workspace.fs;
        auto fm              = workspace.fm;
        auto fluxn           = workspace.fluxn;
        auto fluxs           = workspace.fluxs;
        auto fluxm           = workspace.fluxm;
        auto flux_fullact    = workspace.flux_fullact;

        auto ylo             = workspace.ylo;
        auto yhi             = workspace.yhi;
        auto ymean           = workspace.ymean;
        auto alnsign         = workspace.alnsign;
        auto am              = workspace.am;
        auto lnhygro         = workspace.lnhygro;
        auto f1              = workspace.f1;
        auto sm_activate     = workspace.sm_activate;
        auto lnsm            = workspace.lnsm;
        auto eta             = workspace.eta;
        auto asub            = workspace.asub;
        auto bsub            = workspace.bsub;
        auto sumflxn         = workspace.sumflxn;
        auto sumflxs         = workspace.sumflxs;
        auto sumflxm         = workspace.sumflxm;
        auto sumfn           = workspace.sumfn;
        auto sumfs           = workspace.sumfs;
        auto sumfm           = workspace.sumfm;
        auto fnold           = workspace.fnold;
        auto fsold           = workspace.fsold;
        auto fmold           = workspace.fmold;
        auto rlo             = workspace.rlo;
        auto rhi             = workspace.rhi;
        auto rmean           = workspace.rmean;
        auto etafactor2      = workspace.etafactor2;
        auto totn            = workspace.totn;
        auto gmrad           = workspace.gmrad;
        auto gmradsq         = workspace.gmradsq;
        auto gmlnsig         = workspace.gmlnsig;
        auto gmsm            = workspace.gmsm;
        auto sumns           = workspace.sumns;
        auto bin_is_empty    = workspace.bin_is_empty;
        auto bin_is_narrow   = workspace.bin_is_narrow;
        auto tmp_nsize_aer   = workspace.tmp_nsize_aer;

        auto nact            = workspace.nact;
        auto mact            = workspace.mact;
        auto npv             = workspace.npv;
        real scale;

        auto hygro_aer       = workspace.hygro_aer;
        constexpr int psat=6;                                                //! number of supersaturations to calc ccn concentration
        auto ccn             = workspace.ccn;

        constexpr real surften = 0.076;                                      //! surface tension of water w/respect to air (N/m)
        auto ccnfact         = workspace.ccnfact;
        auto amcube          = workspace.amcube;
        auto argfactor       = workspace.argfactor;
        real aten;                                                           //! surface tension parameter
        real t0;                                                             //! reference temperature
        real sm;                                                             //! critical supersaturation
        real arg;

        constexpr int icheck_colmass = 0;                                   //! icheck_colmass > 0 turns on mass/number conservation checking
                                                                            //! values of 1, 10, 100 produce less to more diagnostics
        int4d colmass_worst_ij;
        int1d colmass_maxworst_i;
        real5d colmass_bgn;
        real5d colmass_end;
        real5d colmass_sfc;
        real3d colmass_worst;
        real colmass_maxworst_r;
        auto rhodz = workspace.rhodz;
        real rhodzsum;

        real smax_prescribed = 1.0; //activate's variable

        if (icheck_colmass > 0) {
            colmass_worst_ij  = int4d("colmass_worst_ij", 2, maxd_acomp+1, maxd_asize, maxd_atype);
            colmass_maxworst_i= int1d("colmass_maxworst_i", 3);
            colmass_bgn   = real5d("colmass_bgn", maxd_acomp+1, maxd_asize, maxd_atype,maxd_aphase, workspace_ncol);
            colmass_end   = real5d("colmass_end", maxd_acomp+1, maxd_asize, maxd_atype,maxd_aphase, workspace_ncol);
            colmass_sfc   = real5d("colmass_sfc", maxd_acomp+1, maxd_asize, maxd_atype,maxd_aphase, workspace_ncol);
            colmass_worst     = real3d("colmass_worst", maxd_acomp+1, maxd_asize, maxd_atype);

            for (int ac = 0; ac < maxd_acomp+1; ++ac) {
                for (int as = 0; as < maxd_asize; ++as) {
                    for (int at = 0; at < maxd_atype; ++at) {
                        colmass_worst_ij(0,ac,as,at) = -1;
                        colmass_worst_ij(1,ac,as,at) = -1;
                        colmass_worst(ac,as,at) = 0.0;
                    }
                }
            }
        }

        arg = 1.0;
        // if (abs(0.8427-ERF_ALT(arg))/0.8427>0.001) then
        //     call chemistry_error_fatal('dropmixnuc: Error function error')
        // endif
        arg = 0.0;
        // if (abs(ERF_ALT(arg)) > 1.0e-6) then  !Jiawang Feng : Precison issues may happen in some machines
        //     call chemistry_error_fatal('dropmixnuc: Error function error')
        // endif

        pi = 4.*Kokkos::atan(1.0);
        dtinv=1./dtstep;

        depvel_drop =  0.1; //! prescribed here rather than getting it from dry_dep_driver
        if (idrydep_onoff <= 0) depvel_drop =  0.0;
        depvel_drop =  Kokkos::min(depvel_drop,depvel_uplimit);
        const real zkmin_local = zkmin;
        const real zkmax_local = zkmax;
        const real wmixmin_local = wmixmin;

        t0=273.15;                   //!wig, 1-Mar-2009: Added .15
        aten=2.*surften/(r_v*t0*rhowater);

        // The serial backend reuses one workspace column. Parallel backends
        // reuse a fixed 64-column workspace, including for the tail chunk.
#ifdef USE_SERIALIZATION
        for (int column = 0; column < ncol; ++column) {
        const int current_chunk_len = 1;
        const auto column_range = std::make_pair(column, column + 1);
#else
        for (int chunk_start = 0; chunk_start < ncol;
             chunk_start += CHEM_WORKSPACE_CHUNK_SIZE) {
        const int current_chunk_len = std::min(
            CHEM_WORKSPACE_CHUNK_SIZE, ncol - chunk_start);
        const auto column_range = std::make_pair(
            chunk_start, chunk_start + current_chunk_len);
#endif

        auto chem = Kokkos::subview(
            chem_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
        auto qv = Kokkos::subview(qv_all, Kokkos::ALL(), column_range);
        auto qc = Kokkos::subview(qc_all, Kokkos::ALL(), column_range);
        auto qi = Kokkos::subview(qi_all, Kokkos::ALL(), column_range);
        auto qs = Kokkos::subview(qs_all, Kokkos::ALL(), column_range);
        auto qndrop3d = Kokkos::subview(
            qndrop3d_all, Kokkos::ALL(), column_range);
        auto tend_qndrop_activation = Kokkos::subview(
            tend_qndrop_activation_all, Kokkos::ALL(), column_range);
        auto temp = Kokkos::subview(temp_all, Kokkos::ALL(), column_range);
        auto p_phy = Kokkos::subview(p_phy_all, Kokkos::ALL(), column_range);
        auto w = Kokkos::subview(w_all, Kokkos::ALL(), column_range);
        auto ddvel = Kokkos::subview(ddvel_all, Kokkos::ALL(), column_range);
        auto hygro = Kokkos::subview(
            hygro_all, Kokkos::ALL(), column_range,
            Kokkos::ALL(), Kokkos::ALL());
        auto rho = Kokkos::subview(rho_all, Kokkos::ALL(), column_range);
        auto zm = Kokkos::subview(zm_all, Kokkos::ALL(), column_range);
        auto dz8w = Kokkos::subview(dz8w_all, Kokkos::ALL(), column_range);
        auto kvh = Kokkos::subview(kvh_all, Kokkos::ALL(), column_range);
        auto cldfra = Kokkos::subview(
            cldfra_all, Kokkos::ALL(), column_range);
        auto cldfra_old = Kokkos::subview(
            cldfra_old_all, Kokkos::ALL(), column_range);
        auto qsrflx = Kokkos::subview(
            qsrflx_all, Kokkos::ALL(), column_range);
        auto ccn1 = Kokkos::subview(ccn1_all, Kokkos::ALL(), column_range);
        auto ccn2 = Kokkos::subview(ccn2_all, Kokkos::ALL(), column_range);
        auto ccn3 = Kokkos::subview(ccn3_all, Kokkos::ALL(), column_range);
        auto ccn4 = Kokkos::subview(ccn4_all, Kokkos::ALL(), column_range);
        auto ccn5 = Kokkos::subview(ccn5_all, Kokkos::ALL(), column_range);
        auto ccn6 = Kokkos::subview(ccn6_all, Kokkos::ALL(), column_range);
        auto nsource = Kokkos::subview(
            nsource_all, Kokkos::ALL(), column_range);

        //OVERALL_MAIN_I_LOOP:
        Kokkos::parallel_for("mixactivate_main_i",
            Kokkos::RangePolicy<>(0, current_chunk_len),
            KOKKOS_LAMBDA(int i) {
            int km1,kp1;
            real wbar,wmix,wdiab,wmin,wmax;
            real tmpa, tmpb, tmpc, tmpc1, tmpc2, tmpd, tmpe, tmpf;
            real tmpcourno;
            real dact;
            real fluxntot;
            real fac_srflx;
            real depvel_tmp;
            real surfratemax;
            real surfrate_drop;
            real dtmin,tinv,dtt;
            int nsubmix,nsubmix_bnd;
            real dtmix;
            real alogarg;
            int nnew,nsav,ntemp;
            int lnum,lnumcw,lmass,lmasscw,lwater;
            real scale;
            real sm;
            real arg;
            real rhodzsum = 0.0;

            for (int l = 0; l < num_chem; ++l) {
                qsrflx(l,i) = 0.0;
            }
            for (int n = 0; n < ntype_aer; ++n) {
                for (int m = 0; m < nsize_aer(n); ++m) {
                    real alogsig_local = Kokkos::log(sigmag_aer(m,n));
                    real dgnum_local = dpvolmean_aer(m,n) *
                        Kokkos::exp(-1.5*alogsig_local*alogsig_local);
                    real exp45logsig_local =
                        Kokkos::exp(4.5*alogsig_local*alogsig_local);
                    npv(m,n,i) = 6.0/(pi*Kokkos::pow(0.01*dgnum_local,3)*
                        exp45logsig_local);
                    argfactor(m,n,i) = 2.0/(3.0*Kokkos::sqrt(2.0)*alogsig_local);
                    amcube(m,n,i) = 3.0/(4.0*pi*exp45logsig_local*npv(m,n,i));
                }
            }

            // !------------------------------------------------------------
            // !         load number nucleated into qndrop on cloud boundaries
            // !------------------------------------------------------------

            for(int k = 1; k<kte; k++){
                zs(k,i)=1./(zm(k,i)-zm(k-1,i));
            }
            zs(0,i) = zs(1,i);
            zs(kte,i) = 0;

            for (int k = 0; k<kte; k++){

                real qcld_local;
                if(f_qi) qcld_local=qc(k,i)+qi(k,i);
                else qcld_local=qc(k,i);

                if(qcld_local < -1.|| qcld_local > 1.){
                    Kokkos::abort("mixactivate: qcld out of range");
                }
                if(qcld_local > 1.e-20){
                    lcldfra(k,i)=cldfra(k,i)*qc(k,i)/qcld_local;
                    lcldfra_old(k,i)=cldfra_old(k,i)*qc(k,i)/qcld_local;
                }
                else{
                    lcldfra(k,i)=0.;
                    lcldfra_old(k,i)=0.;
                }
                qndrop(k,i)=qndrop3d(k,i);
                cs(k,i)=rho(k,i);              //! air density (kg/m3)
                dz(k,i)=dz8w(k,i);
                for(int n = 0; n < ntype_aer; n++){
                    for(int m = 0; m < nsize_aer(n); m++){
                        nact(k,m,n,i)=0.;
                        mact(k,m,n,i)=0.;
                    }
                }
                zn(k,i)=1./(cs(k,i)*dz(k,i));
                if(k>0){
                    ekd(k,i)=kvh(k,i);
                    ekd(k,i)=Kokkos::max(ekd(k,i),zkmin_local);
                    ekd(k,i)=Kokkos::min(ekd(k,i),zkmax_local);
                }
                else ekd(k,i)=0;
                if(k == 0){
                    wtke(k,i)=sq2pi*depvel_drop;
                }
                else wtke(k,i)=sq2pi*ekd(k,i)/dz(k,i);
                wtke(k,i)=Kokkos::max(wtke(k,i),wmixmin_local);
                nsource(k,i)=0.;
            }
            nsource(nsource.extent(0)-1,i) = 0.;
            qndrop(kte,i)      = 0.;
            zn(kte,i)          = 0.;

            for(int k=1; k<kte; k++){
                real tmpa_local = dz(k-1,i);
                real tmpb_local = dz(k,i);
                real tmpc_local = tmpa_local/(tmpa_local + tmpb_local);
                csbot(k,i) = cs(k-1,i)*(1.0-tmpc_local) + cs(k,i)*tmpc_local;
                csbot_cscen(k,i) = csbot(k,i)/cs(k,i);
            }
            csbot(kts - 1,i) = cs(kts - 1,i);
	        csbot_cscen(kts - 1,i) = 1.0;
	        csbot(kte,i) = cs(kte - 1,i);
            csbot_cscen(kte,i) = 1.0;

            //!  calculate surface rate and mass mixing ratio for aerosol
            surfratemax = 0.0;
            // nsav=1;
            // nnew=2;
            nsav=0;
            nnew=1;
            surfrate_drop=depvel_drop/dz(0,i);
            surfratemax = Kokkos::max( surfratemax, surfrate_drop );
            for(int n = 0; n<ntype_aer; n++){
                for(int m = 0; m<nsize_aer(n); m++){
                    lnum=numptr_aer(m,n,ai_phase - 1);
                    lnumcw=numptr_aer(m,n,cw_phase - 1);
                    if(lnum > 0){
                        depvel_tmp = Kokkos::max( 0.0, Kokkos::min( ddvel(lnum-1,i), depvel_uplimit ) );
                        surfrate(lnum-1,i)=depvel_tmp/dz(0,i);
                        surfrate(lnumcw-1,i)=surfrate_drop;
                        surfratemax = Kokkos::max( surfratemax, surfrate(lnum-1,i) );
                        scale = 1.;
                        for(int k = 0; k < kte; k++){
                            raercol(k,lnumcw-1,nsav,i)=chem(lnumcw-1,k,i)*scale;
                            raercol(k,lnum-1,nsav,i)=chem(lnum-1,k,i)*scale;
                        }
                    }
                    for(int l = 0; l < ncomp_aer(n); l++){
                        lmass=massptr_aer(l,m,n,ai_phase - 1);
                        lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                        scale = 1.e-9;
                        depvel_tmp = Kokkos::max( 0.0, Kokkos::min( ddvel(lmass-1,i), depvel_uplimit ) );
                        surfrate(lmass-1,i)=depvel_tmp/dz(0,i);
                        surfrate(lmasscw-1,i)=surfrate_drop;
                        surfratemax = Kokkos::max( surfratemax, surfrate(lmass-1,i) );
                        for(int k = 0; k < kte; k++){
                            raercol(k,lmasscw-1,nsav,i)=chem(lmasscw-1,k,i)*scale;
                            raercol(k,lmass-1,nsav,i)=chem(lmass-1,k,i)*scale;
                        }
                    }
                    lwater=waterptr_aer(m,n);
                    if(lwater > 0){
                        depvel_tmp = Kokkos::max( 0.0, Kokkos::min( ddvel(lwater-1,i), depvel_uplimit ) );
                        surfrate(lwater-1,i)=depvel_tmp/dz(0,i);
                        surfratemax = Kokkos::max( surfratemax, surfrate(lwater-1,i) );
                        for(int k = 0; k < kte; k++){
                            raercol(k,lwater-1,nsav,i)=chem(lwater-1,k,i);
                        }
                    }
                }// size
            }// type

            //! mass conservation checking
            if (icheck_colmass > 0){
                for(int ac = 0; ac<maxd_acomp+1;ac++){
                    for(int as = 0; as<maxd_asize; as++){
                        for(int at = 0; at<maxd_atype; at++){
                            for(int phase = 0; phase<maxd_aphase; phase++){
                                colmass_bgn(ac,as,at,phase,i) = 0.0;
                                colmass_end(ac,as,at,phase,i) = 0.0;
                                colmass_sfc(ac,as,at,phase,i) = 0.0;
                            }
                        }
                    }
                }
                for(int k = 0; k<kte; k++){
                    rhodz(k,i) = 1.0/zn(k,i);
                    rhodzsum = rhodzsum + rhodz(k,i);
                }
                for(int n = 0; n<ntype_aer; n++){
                    for(int m = 0; m<nsize_aer(n); m++){
                        lnum=numptr_aer(m,n,ai_phase - 1);
                        lnumcw=numptr_aer(m,n,cw_phase - 1);
                        if(lnum>0){
                            for(int k = 0; k<kte; k++){
                                colmass_bgn(0,m,n,0,i) = colmass_bgn(0,m,n,0,i)+chem(lnum-1,k,i  )*rhodz(k,i);
                                colmass_bgn(0,m,n,1,i) = colmass_bgn(0,m,n,1,i)+chem(lnumcw-1,k,i)*rhodz(k,i);
                            }
                        }
                        for(int l = 0; l < ncomp_aer(n); l++){
                            lmass=massptr_aer(l,m,n,ai_phase - 1);
                            lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                            for(int k = 0; k<kte; k++){
                                colmass_bgn(l+1,m,n,0,i) = colmass_bgn(l+1,m,n,0,i)+chem(lmass-1,k,i  )*rhodz(k,i);
                                colmass_bgn(l+1,m,n,1,i) = colmass_bgn(l+1,m,n,1,i)+chem(lmasscw-1,k,i)*rhodz(k,i);
                            }
                        }
                    }   //size
                }   // type
            }   //icheck_colmass > 0

            // !------------------------------------------------------------
            // !        droplet nucleation/aerosol activation
            // !------------------------------------------------------------
            //! k-loop for growing/shrinking cloud calcs .............................
            //GROW_SHRINK_MAIN_K_LOOP:
            for(int k = 0; k < kte; k++){
                km1=Kokkos::max(k-1,0);
                kp1=Kokkos::min(k+1,kte-1);
                tmpc1 = Kokkos::max( (lcldfra(k,i)-lcldfra_old(k,i)), 0.0 );

                if(k > 0){
                    tmpcourno = dtstep*Kokkos::max(w(k,i),0.0)/dz(k,i);
                    tmpc2 = Kokkos::max( (lcldfra(k,i)-lcldfra(km1,i)), 0.0 ) * tmpcourno;
                    tmpc2 = Kokkos::min( tmpc2, 1.0 );
                }
                else tmpc2 = 0.0;

                if ((tmpc1 > 0.001) || (tmpc2 > 0.001)){
                    wbar=w(k,i)+wtke(k,i);
                    wmix=0.;
                    wmin=0.;
                    wmax=50.;
                    wdiab=0;

                    for(int n = 0; n < ntype_aer; n++){
                        for(int m = 0; m < nsize_aer(n); m++){
                            loadaer(raercol,k,
                                cs(k,i), npv, dlo_sect,dhi_sect,
                                maxd_acomp, ncomp_aer,
                                numptr_aer,
                                dens_aer,
                                massptr_aer,
                                maerosol, maerosolcw,
                                maerosol_tot, maerosol_totcw,
                                naerosol, naerosolcw,
                                vaerosol, vaerosolcw,m,n,ai_phase,cw_phase,i);
                            hygro_aer(m,n,i)=hygro(k,i,m,n);
                        }
                    }

                    activate(wbar,wmix,wdiab,wmin,wmax,temp(k,i),cs(k,i),
                        msectional, maxd_atype, ntype_aer, maxd_asize, nsize_aer,
                        naerosol, vaerosol,
                        dlo_sect,dhi_sect,sigmag_aer,hygro_aer,
                        fn,fs,fm,fluxn,fluxs,fluxm,flux_fullact,
                        ylo,yhi,ymean,alnsign,am,lnhygro,f1,sm_activate,lnsm,
                        eta,asub,bsub,sumflxn,sumflxs,sumflxm,sumfn,sumfs,sumfm,
                        fnold,fsold,fmold,rlo,rhi,rmean,etafactor2,
                        totn,gmrad,gmradsq,gmlnsig,gmsm,sumns,
                        bin_is_empty,bin_is_narrow,tmp_nsize_aer,
                        grid_id,ktau,i,k,smax_prescribed);

                    for(int n = 0; n < ntype_aer; n++){
                        for(int m = 0; m < nsize_aer(n); m++){
                            lnum   = numptr_aer(m,n,ai_phase - 1);
                            lnumcw = numptr_aer(m,n,cw_phase - 1);
                            if (tmpc1 > 0.0) dact = tmpc1*fn(m,n,i)*raercol(k,lnum-1,nsav,i);
                            else dact = 0.0;
                            if (tmpc2 > 0.0) dact = dact + tmpc2*fn(m,n,i)*raercol(km1,lnum-1,nsav,i);
                            dact = Kokkos::min( dact, 0.99*raercol(k,lnum-1,nsav,i) );
                            raercol(k,lnumcw-1,nsav,i) = raercol(k,lnumcw-1,nsav,i)+dact;
                            raercol(k,lnum-1,nsav,i) = raercol(k,lnum-1,nsav,i)-dact;
                            tend_qndrop_activation(k,i) = dact*dtinv;
                            qndrop(k,i) = qndrop(k,i)+dact;
                            nsource(k,i) = nsource(k,i)+dact*dtinv;
                            for(int l = 0;l < ncomp_aer(n);l++){
                                lmass   = massptr_aer(l,m,n,ai_phase - 1);
                                lmasscw = massptr_aer(l,m,n,cw_phase - 1);
                                if (tmpc1 > 0.0) dact = tmpc1*fm(m,n,i)*raercol(k,lmass-1,nsav,i);          //! interstitial only
                                else dact = 0.0;
                                if (tmpc2 > 0.0) dact = dact + tmpc2*fm(m,n,i)*raercol(km1,lmass-1,nsav,i); //! interstitial only
                                dact = Kokkos::min( dact, 0.99*raercol(k,lmass-1,nsav,i) );
                                raercol(k,lmasscw-1,nsav,i)  =  raercol(k,lmasscw-1,nsav,i)+dact;
                                raercol(k,lmass-1,nsav,i)  =  raercol(k,lmass-1,nsav,i)-dact;
                            }
                        }
                    }
                }   //end ((tmpc1 > 0.001) || (tmpc2 > 0.001))

                if(lcldfra(k,i) < lcldfra_old(k,i) && lcldfra_old(k,i) > 1.e-20){
                    nsource(k,i)=nsource(k,i)+qndrop(k,i)*(lcldfra(k,i)-lcldfra_old(k,i))*dtinv;
                    qndrop(k,i)=qndrop(k,i)*(1.+lcldfra(k,i)-lcldfra_old(k,i));
                    tmpc = (lcldfra(k,i)-lcldfra_old(k,i))/lcldfra_old(k,i);
                    for(int n = 0; n < ntype_aer; n++){
                        for(int m = 0; m < nsize_aer(n); m++){
                            lnum=numptr_aer(m,n,ai_phase - 1);
                            lnumcw=numptr_aer(m,n,cw_phase - 1);
                            if(lnum > 0){
                                dact=raercol(k,lnumcw-1,nsav,i)*tmpc;
                                raercol(k,lnumcw-1,nsav,i)=raercol(k,lnumcw-1,nsav,i)+dact;
                                raercol(k,lnum-1,nsav,i)=raercol(k,lnum-1,nsav,i)-dact;
                            }
                            for(int l = 0 ; l < ncomp_aer(n);l++){
                                lmass=massptr_aer(l,m,n,ai_phase - 1);
                                lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                                dact=raercol(k,lmasscw-1,nsav,i)*tmpc;
                                raercol(k,lmasscw-1,nsav,i)=raercol(k,lmasscw-1,nsav,i)+dact;
                                raercol(k,lmass-1,nsav,i)=raercol(k,lmass-1,nsav,i)-dact;
                            }
                        }
                    }
                }
            }//end GROW_SHRINK_MAIN_K_LOOP:

        // ! end of k-loop for growing/shrinking cloud calcs ......................


        // ! ......................................................................
        // ! start of main k-loop for calc of old cloud activation tendencies ..........
        // ! this loop does "set up" for the nsubmix loop
        // !
        // ! rce-comment
        // !    changed this part of code to use current cloud fraction (lcldfra) exclusively
        //OLD_CLOUD_MAIN_K_LOOP:
        for(int k = 0; k < kte; k++){
            km1=Kokkos::max(k-1,kts-1);
            kp1=Kokkos::min(k+1,kte-1);
            flux_fullact(k,i) = 0.0;

            if(lcldfra(k,i) > 0.01){

                //old cloud
                if(lcldfra(k,i)-lcldfra(km1,i) > 0.01 || k == 0){
                    wdiab=0;
                    wmix=wtke(k,i);   //! spectrum of updrafts
                    wbar=w(k,i);    //! spectrum of updrafts
                    wmax=50.;
                    ekd(k,i)=wtke(k,i)*dz(k,i)/sq2pi;
                    alogarg=Kokkos::max(1.e-20,1/lcldfra(k,i)-1.);
                    wmin=wbar+wmix*0.25*sq2pi*Kokkos::log(alogarg);
                    for(int n = 0; n < ntype_aer; n++){
                        for(int m = 0; m < nsize_aer(n); m++){
                                loadaer(raercol,km1,
                                    cs(k,i), npv,dlo_sect,dhi_sect,
                                    maxd_acomp,ncomp_aer,
                                    numptr_aer,
                                    dens_aer,
                                    massptr_aer,
                                    maerosol, maerosolcw,
                                    maerosol_tot, maerosol_totcw,
                                    naerosol, naerosolcw,
                                    vaerosol, vaerosolcw,m,n,ai_phase,cw_phase,i);
                                hygro_aer(m,n,i)=hygro(k,i,m,n);
                        }
                    }

                    activate(wbar,wmix,wdiab,wmin,wmax,temp(k,i),cs(k,i),
                          msectional, maxd_atype, ntype_aer, maxd_asize, nsize_aer,
                          naerosol, vaerosol,
                          dlo_sect,dhi_sect, sigmag_aer,hygro_aer,
                          fn,fs,fm,fluxn,fluxs,fluxm,flux_fullact,
                          ylo,yhi,ymean,alnsign,am,lnhygro,f1,sm_activate,lnsm,
                          eta,asub,bsub,sumflxn,sumflxs,sumflxm,sumfn,sumfs,sumfm,
                          fnold,fsold,fmold,rlo,rhi,rmean,etafactor2,
                          totn,gmrad,gmradsq,gmlnsig,gmsm,sumns,
                          bin_is_empty,bin_is_narrow,tmp_nsize_aer,
                          grid_id,ktau,i,k,smax_prescribed);

                    // ! rce-comment
                    // !    the activation-fraction fluxes (fluxn, fluxm) from subr activate assume that
                    // !       wbar << wmix, which is valid for global-model scale but not mesoscale
                    // !    for wrf-chem application, divide these by flux_fullact to get a
                    // !       "flux-weighted-average" activation fraction, then multiply by (ekd(k,i)*zs(k,i))
                    // !       which is the local "turbulent vertical-mixing velocity"

                    if (k > 0){
                        if (flux_fullact(k,i) > 1.0e-20){
                        tmpa = ekd(k,i)*zs(k,i);
                        tmpf = flux_fullact(k,i);
                            for(int n = 0; n < ntype_aer; n++){
                                for(int m = 0; m < nsize_aer(n); m++){
                                    tmpb = Kokkos::max( fluxn(m,n,i), 0.0 ) / Kokkos::max( fluxn(m,n,i), tmpf );
                                    fluxn(m,n,i) = tmpa*tmpb;
                                    tmpb = Kokkos::max( fluxm(m,n,i), 0.0 ) / Kokkos::max( fluxm(m,n,i), tmpf );
                                    fluxm(m,n,i) = tmpa*tmpb;
                                }
                            }
                        }
                        else{
                            for(int n = 0; n < ntype_aer; n++){
                                for(int m = 0; m < nsize_aer(n); m++){
                                    fluxn(m,n,i) = 0.0;
                                    fluxm(m,n,i) = 0.0;
                                }
                            }
                        }
                    }

                    if(k > 0) tmpc = lcldfra(k,i)-lcldfra(km1,i);
                    else tmpc=lcldfra(k,i);

                    tmpe = csbot_cscen(k,i)/(dz(k,i));
                    fluxntot=0.;

                    for(int n = 0; n < ntype_aer; n++){
                        for(int m = 0; m < nsize_aer(n); m++){
                            fluxn(m,n,i)=fluxn(m,n,i)*tmpc;
                            fluxm(m,n,i)=fluxm(m,n,i)*tmpc;
                            lnum=numptr_aer(m,n,ai_phase - 1);
                            fluxntot=fluxntot+fluxn(m,n,i)*raercol(km1,lnum-1,nsav,i);
                            nact(k,m,n,i)=nact(k,m,n,i)+fluxn(m,n,i)*tmpe;
                            mact(k,m,n,i)=mact(k,m,n,i)+fluxm(m,n,i)*tmpe;
                        }
                    }

                    flux_fullact(k,i) = flux_fullact(k,i)*tmpe;
                    nsource(k,i)=nsource(k,i)+fluxntot*zs(k,i);
                    fluxntot=fluxntot*cs(k,i);
                }
            }

            else{               //lcldfra(k,i) <= 0.01
                if(qndrop(k,i) > 10000.e6){
                    //print *,'k,i,lcldfra,qndrop=',k,i,lcldfra(k,i),qndrop(k,i)
                    //print *,'cldfra,ql,qi',cldfra(k,i),qc(k,i),qi(k,i)
                }
                nsource(k,i)=nsource(k,i)-qndrop(k,i)*dtinv;
                qndrop(k,i)=0.;

                for(int n = 0; n < ntype_aer; n++){
                    for(int m = 0; m < nsize_aer(n); m++){
                        lnum=numptr_aer(m,n,ai_phase - 1);
                        lnumcw=numptr_aer(m,n,cw_phase - 1);
                        if(lnum > 0){
                            raercol(k,lnum-1,nsav,i)=raercol(k,lnum-1,nsav,i)+raercol(k,lnumcw-1,nsav,i);
                            raercol(k,lnumcw-1,nsav,i)=0.;
                        }
                        for(int l=0; l < ncomp_aer(n); l++){
                            lmass=massptr_aer(l,m,n,ai_phase - 1);
                            lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                            raercol(k,lmass-1,nsav,i)=raercol(k,lmass-1,nsav,i)+raercol(k,lmasscw-1,nsav,i);
                            raercol(k,lmasscw-1,nsav,i)=0.;
                        }
                    }
                }
            }//lcldfra(k,i)
        }//end OLD_CLOUD_MAIN_K_LOOP
    //////////////////////////////////////////////////////////////
        ntemp=nsav;
        nsav=nnew;
        nnew=ntemp;

        //!load new droplets in layers above, below cloud
        dtmin=dtstep;
        ekk(0,i)=0.0;

        // ! rce-comment -- ekd(k,i) is eddy-diffusivity at k/k-1 interface
        // !   want ekk(k,i) = ekd(k,i) * (density at k/k-1 interface)
        for(int k = 1; k < kte; k++) ekk(k,i)=ekd(k,i)*csbot(k,i);

        ekk(kte,i)=0.0;

        for(int k = 0; k < kte; k++){
            ekkp(k,i)=zn(k,i)*ekk(k+1,i)*zs(k+1,i);
            ekkm(k,i)=zn(k,i)*ekk(k,i)*zs(k,i);
            tinv=ekkp(k,i)+ekkm(k,i);
            if(k == 0) tinv=tinv+surfratemax;
            if(tinv > 1.e-6){
                dtt=1./tinv;
                dtmin=Kokkos::min(dtmin,dtt);
            }
        }
        dtmix=0.9*dtmin;
        nsubmix=dtstep/dtmix+1;
        if(nsubmix>100) nsubmix_bnd=100;
        else nsubmix_bnd=nsubmix;

        dtmix=dtstep/nsubmix;
        fac_srflx = -1.0/(zn(0,i)*nsubmix);

        for(int k = 0; k < kte; k++){
            kp1=Kokkos::min(k+1,kte-1);
            km1=Kokkos::max(k-1,0);
            if(lcldfra(kp1,i) > 0) overlapp(k,i)=Kokkos::min(lcldfra(k,i)/lcldfra(kp1,i),1.);
            else overlapp(k,i)=1.;
            if(lcldfra(km1,i) > 0) overlapm(k,i)=Kokkos::min(lcldfra(k,i)/lcldfra(km1,i),1.);
            else overlapm(k,i)=1.;
        }

        // ! ......................................................................
        // ! start of nsubmix-loop for calc of old cloud activation tendencies ....
        //OLD_CLOUD_NSUBMIX_LOOP:
        for(int nsub = 0; nsub < nsubmix; nsub++){
            for(int k = 0; k < kte; k++){
                qndrop_new(k,i)=qndrop(k,i);
                srcn(k,i)=0.0;}
            ntemp=nsav;
            nsav=nnew;
            nnew=ntemp;
            for(int n = 0; n < ntype_aer; n++){
                for(int m = 0; m < nsize_aer(n); m++){
                    lnum=numptr_aer(m,n,ai_phase - 1);
                    for(int k = 0; k < kte-1; k++){
                        srcn(k+1,i)=srcn(k+1,i)+nact(k+1,m,n,i)*(raercol(k,lnum-1,nsav,i));
                    }
                    srcn(0,i)=srcn(0,i)+nact(0,m,n,i)*(raercol(0,lnum-1,nsav,i));
                }
            }

            explmix_f(qndrop,srcn,ekkp,ekkm,overlapp,overlapm,
                 qndrop_new,surfrate_drop,kts,kte,dtmix,i);
            for(int n = 0; n < ntype_aer; n++){
                for(int m = 0; m < nsize_aer(n); m++){
                    lnum=numptr_aer(m,n,ai_phase - 1);
                    lnumcw=numptr_aer(m,n,cw_phase - 1);
                    if(lnum>0){
                        for(int k = 0; k < kte-1; k++){
                            source(k+1,i)=nact(k+1,m,n,i)*(raercol(k,lnum-1,nsav,i));
                        }
                        source(0,i)= nact(0,m,n,i)*(raercol(0,lnum-1,nsav,i));

                        explmix_f_3d(raercol,source,ekkp,ekkm,overlapp,overlapm,
                            raercol,surfrate(lnumcw-1,i),kts,kte,dtmix,lnumcw-1,nnew,nsav,i);

                        explmix_t_3d(raercol,source,ekkp,ekkm,overlapp,overlapm,
                            raercol,surfrate(lnum-1,i),kts,kte,dtmix,
                            raercol,lnum-1,lnumcw-1,nnew,nsav,i);

                        qsrflx(lnum-1,i) = qsrflx(lnum-1,i) + fac_srflx * raercol(kts-1,lnum-1,nsav,i)*surfrate(lnum-1,i);
                        qsrflx(lnumcw-1,i) = qsrflx(lnumcw-1,i) + fac_srflx * raercol(kts-1,lnumcw-1,nsav,i)*surfrate(lnumcw-1,i);
                        if (icheck_colmass > 0){
                            tmpf = dtmix*rhodz(kts-1,i);
                            colmass_sfc(0,m,n,0,i) = colmass_sfc(0,m,n,0,i) + raercol(kts-1,lnum-1,nsav,i)*surfrate(lnum-1,i)*tmpf;
                            colmass_sfc(0,m,n,1,i) = colmass_sfc(0,m,n,1,i) + raercol(kts-1,lnumcw-1,nsav,i)*surfrate(lnumcw-1,i)*tmpf;
                        }
                    }

                    for(int l = 0;l < ncomp_aer(n);l++){
                        lmass=massptr_aer(l,m,n,ai_phase - 1);
                        lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                        for(int k = 0; k < kte-1; k++) source(k+1,i)=mact(k+1,m,n,i)*(raercol(k,lmass-1,nsav,i));
                        source(kts-1,i)= mact(kts-1,m,n,i)*(raercol(kts-1,lmass-1,nsav,i));

                        explmix_f_3d(raercol,source,ekkp,ekkm,overlapp,overlapm,
                            raercol,surfrate(lmasscw-1,i),kts,kte,dtmix,lmasscw-1,nnew,nsav,i);

                        explmix_t_3d(raercol,source,ekkp,ekkm,overlapp,overlapm,
                            raercol,surfrate(lmass-1,i),kts,kte,dtmix,
                            raercol,lmass-1,lmasscw-1,nnew,nsav,i);

                        qsrflx(lmass-1,i) = qsrflx(lmass-1,i) + fac_srflx * raercol(kts-1,lmass-1,nsav,i)*surfrate(lmass-1,i);
                        qsrflx(lmasscw-1,i) = qsrflx(lmasscw-1,i) + fac_srflx * raercol(kts-1,lmasscw-1,nsav,i)*surfrate(lmasscw-1,i);
                        if (icheck_colmass > 0){
                            tmpf = dtmix*rhodz(kts-1,i)*1.0e9;
                            colmass_sfc(l+1,m,n,0,i) = colmass_sfc(l+1,m,n,0,i) + raercol(kts-1,lmass-1,nsav,i)*surfrate(lmass-1,i)*tmpf;
                            colmass_sfc(l+1,m,n,1,i) = colmass_sfc(l+1,m,n,1,i) + raercol(kts-1,lmasscw-1,nsav,i)*surfrate(lmasscw-1,i)*tmpf;
                        }
                    }
                    lwater=waterptr_aer(m,n);
                    if(lwater>0){
                        for(int k = 0; k < kte; k++) source(k,i)=0.;
                         explmix_t(raercol,source,ekkp,ekkm,overlapp,overlapm,
                            raercol,surfrate(lwater-1,i),kts,kte,dtmix,source,lwater-1,lwater-1,nnew,nsav,i);
                    }
                }//size
            }//type
        }//OLD_CLOUD_NSUBMIX_LOOP

        for(int k = 0; k < kte;k++){
            if(lcldfra(k,i) == 0.){
                qndrop(k,i)=0.;
                for(int n = 0; n < ntype_aer; n++){
                    for(int m = 0; m < nsize_aer(n); m++){
                        lnum=numptr_aer(m,n,ai_phase - 1);
                        lnumcw=numptr_aer(m,n,cw_phase - 1);
                        if(lnum > 0){
                            raercol(k,lnum-1,nnew,i)=raercol(k,lnum-1,nnew,i)+raercol(k,lnumcw-1,nnew,i);
                            raercol(k,lnumcw-1,nnew,i)=0.;
                        }
                        for(int l = 0; l < ncomp_aer(n); l++){
                            lmass=massptr_aer(l,m,n,ai_phase - 1);
                            lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                            raercol(k,lmass-1,nnew,i)=raercol(k,lmass-1,nnew,i)+raercol(k,lmasscw-1,nnew,i);
                            raercol(k,lmasscw-1,nnew,i)=0.;
                        }
                    }
                }
            }
        }

        //!        droplet number
        for(int k = 0; k < kte;k++){
            qndrop3d(k,i) = Kokkos::max(qndrop(k,i),1.e-6);

            if(qc(k,i) < -1. || qc(k,i) > 1.){
                Kokkos::abort("mixactivate: qc out of range");
            }
            if(qi(k,i) < -1. || qi(k,i) > 1.){
                Kokkos::abort("mixactivate: qi out of range");
            }
            if(qv(k,i) < -1. || qv(k,i) > 1.){
                Kokkos::abort("mixactivate: qv out of range");
            }
            cldfra_old(k,i) = cldfra(k,i);
        }

        for(int p = 0; p < psat; p++){
            for(int k = 0; k < kte; k++){
                ccn(p,k,i) = 0.;
            }
        }
        for(int n = 0; n < ntype_aer; n++){
            for(int m = 0; m < nsize_aer(n); m++){
                lnum=numptr_aer(m,n,ai_phase - 1);
                lnumcw=numptr_aer(m,n,cw_phase - 1);
                if(lnum > 0){
                    scale = 1.;
                    for(int k = 0; k < kte; k++){
                        chem(lnumcw-1,k,i)= raercol(k,lnumcw-1,nnew,i)*scale;
                        chem(lnum-1,k,i)= raercol(k,lnum-1,nnew,i)*scale;
                    }
                }
                for(int l = 0; l < ncomp_aer(n); l++){
                    lmass=massptr_aer(l,m,n,ai_phase - 1);
                    lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                    scale = 1.e9;
                    for(int k = 0; k < kte; k++){
                        chem(lmasscw-1,k,i)=raercol(k,lmasscw-1,nnew,i)*scale; //! ug/kg
                        chem(lmass-1,k,i)=raercol(k,lmass-1,nnew,i)*scale; //! ug/kg
                    }
                }
                lwater=waterptr_aer(m,n);

                for(int k = 0; k < kte; k++){
                    if(lwater>0) chem(lwater-1,k,i)=raercol(k,lwater-1,nnew,i);   //! don't convert units

                    sm=2.*aten*Kokkos::sqrt(aten/(27.*hygro(k,i,m,n)*amcube(m,n,i)));
                    for(int l = 0; l < psat; l++){
                        real supersat_local;
                        if (l == 0) supersat_local = 0.02;
                        else if (l == 1) supersat_local = 0.05;
                        else if (l == 2) supersat_local = 0.1;
                        else if (l == 3) supersat_local = 0.2;
                        else if (l == 4) supersat_local = 0.5;
                        else supersat_local = 1.0;
                        arg=argfactor(m,n,i)*Kokkos::log(sm/(0.01*supersat_local));
                        if(arg<2){
                            if(arg<-2) ccnfact(l,m,n,i)=1.e-6;    //! convert from #/m3 to #/cm3
                            else ccnfact(l,m,n,i)=1.e-6*0.5*erfc_num_recipes(arg);
                        }
                        else ccnfact(l,m,n,i) = 0.;
                        ccn(l,k,i)=ccn(l,k,i)+(raercol(k,lnum-1,nnew,i)+raercol(k,lnumcw-1,nnew,i))*cs(k,i)*ccnfact(l,m,n,i);
                    }
                }
            }
        }

        for(int l = 0; l < psat; l++){
            for(int k = 0; k < kte; k++){
                if(l==0)ccn1(k,i)=ccn(l,k,i);
                if(l==1)ccn2(k,i)=ccn(l,k,i);
                if(l==2)ccn3(k,i)=ccn(l,k,i);
                if(l==3)ccn4(k,i)=ccn(l,k,i);
                if(l==4)ccn5(k,i)=ccn(l,k,i);
                if(l==5)ccn6(k,i)=ccn(l,k,i);
            }
        }
        if (icheck_colmass > 0){
            for(int n = 0; n < ntype_aer; n++){
                for(int m = 0; m < nsize_aer(n); m++){
                    lnum=numptr_aer(m,n,ai_phase - 1);
                    lnumcw=numptr_aer(m,n,cw_phase - 1);
                    if(lnum>0){
                        for(int k = 0; k < kte; k++){
                            colmass_end(0,m,n,0,i) += chem(lnum-1,k,i  )*rhodz(k,i);
                            colmass_end(0,m,n,1,i) += chem(lnumcw-1,k,i)*rhodz(k,i);
                        }
                    }
                    for(int l = 0; l < ncomp_aer(n); l++){
                        lmass=massptr_aer(l,m,n,ai_phase - 1);
                        lmasscw=massptr_aer(l,m,n,cw_phase - 1);
                        for(int k = 0; k < kte; k++){
                            colmass_end(l+1,m,n,0,i) += chem(lmass-1,k,i  )*rhodz(k,i);
                            colmass_end(l+1,m,n,1,i) += chem(lmasscw-1,k,i)*rhodz(k,i);
                        }
                    }
                }
            }

            for(int n = 0; n < ntype_aer; n++){
                for(int m = 0; m < nsize_aer(n); m++){
                    for(int l = 0; l <= ncomp_aer(n); l++){
                        tmpa = ( colmass_bgn(l,m,n,0,i) + colmass_bgn(l,m,n,1,i) )/rhodzsum;
                        tmpb = ( colmass_end(l,m,n,0,i) + colmass_end(l,m,n,1,i) )/rhodzsum;
                        tmpc = ( colmass_sfc(l,m,n,0,i) + colmass_sfc(l,m,n,1,i) )/rhodzsum;
                        tmpd = (tmpb + tmpc) - tmpa;
                        tmpe = Kokkos::max( tmpa, 1.0e-20 );
                        if (Kokkos::abs(tmpd) < 1.0e5*tmpe) tmpf = tmpd/tmpe;
                        else if (tmpf < 0.0) tmpf = -1.0e5;
                        else tmpf = 1.0e5;
                        if (Kokkos::abs(tmpf) > Kokkos::abs(colmass_worst(l,m,n))){
                            colmass_worst(l,m,n) = tmpf;
                            colmass_worst_ij(0,l,m,n) = i+1;
                            colmass_worst_ij(1,l,m,n) = 1;
                        }
                    }
                }
            }
       }//! (icheck_colmass > 0)
        }); //OVERALL_MAIN_I_LOOP ! end of main loop over i
        }
        Kokkos::fence();

    if (icheck_colmass > 0){
        //if (icheck_colmass >= 100) call chemistry_message('mixactivate colmass worst errors bgn - type, size, comp, err, i, j')
        colmass_maxworst_r = 0.0;
        for(int i = 0; i < 3;i++) colmass_maxworst_i(i) = -1;

        for(int n = 0; n < ntype_aer; n++){
            for(int m = 0; m < nsize_aer(n); m++){
                for(int l=0; l <= ncomp_aer(n);l++){
                    // if (icheck_colmass >= 100) &
                    //     write(*,'(3i3,1p,e10.2,2i4)') n, m, l, &
                    //     colmass_worst(l,m,n), colmass_worst_ij(1:2,l,m,n)
                    if (Kokkos::abs(colmass_worst(l,m,n)) > Kokkos::abs(colmass_maxworst_r)){
                        colmass_maxworst_r = colmass_worst(l,m,n);
                        colmass_maxworst_i(0) = n+1;
                        colmass_maxworst_i(1) = m+1;
                        colmass_maxworst_i(2) = l;
                    }
                }
            }
        }
        //  if ((icheck_colmass >= 10) .or. (abs(colmass_maxworst_r) >= 1.0e-6)) &
        //      write(*,'(a,3i3,1p,e10.2)') 'mixactivate colmass maxworst', &
        //      colmass_maxworst_i(1:3), colmass_maxworst_r
    }//if (icheck_colmass > 0)
}

template <class ChemView, class NpvView, class MaerView, class MaerTotView, class Aer2View>
KOKKOS_INLINE_FUNCTION void loadaer(const ChemView &chem,int k,real cs,const NpvView &npv,
                          const real2d &dlo_sect,const real2d &dhi_sect,int maxd_acomp,const int1d &ncomp,
                          const int3d &numptr_aer, const real2d &dens_aer,
                          const int4d &massptr_aer,
                          const MaerView &maerosol, const MaerView &maerosolcw,
                          const MaerTotView &maerosol_tot, const MaerTotView &maerosol_totcw,
                          const Aer2View &naerosol,const Aer2View &naerosolcw,
                          const Aer2View &vaerosol,const Aer2View &vaerosolcw,
                          int m,int n,int ai_phase,int cw_phase,int i){

    int lnum,lnumcw,l,ltype,lmass,lmasscw,lsfc,lsfccw;
    real num_at_dhi, num_at_dlo;
    real npv_at_dhi, npv_at_dlo;
    real pi = 3.1415926526;
    real specvol; //! inverse aerosol material density (m3/kg)


    lnum=numptr_aer(m,n,ai_phase-1);
    lnumcw=numptr_aer(m,n,cw_phase-1);
    maerosol_tot(m,n,i)=0.;
    maerosol_totcw(m,n,i)=0.;
    vaerosol(m,n,i)=0.;
    vaerosolcw(m,n,i)=0.;
    for(int l = 0; l < ncomp(n); l++){
        lmass=massptr_aer(l,m,n,ai_phase-1);
        lmasscw=massptr_aer(l,m,n,cw_phase-1);
        maerosol(l,m,n,i)=chem(k,lmass-1,0,i)*cs;
        maerosol(l,m,n,i)=Kokkos::max(maerosol(l,m,n,i),0.);
        maerosolcw(l,m,n,i)=chem(k,lmasscw-1,0,i)*cs;

        maerosolcw(l,m,n,i)=Kokkos::max(maerosolcw(l,m,n,i),0.);
        maerosol_tot(m,n,i)=maerosol_tot(m,n,i)+maerosol(l,m,n,i);
        maerosol_totcw(m,n,i)=maerosol_totcw(m,n,i)+maerosolcw(l,m,n,i);
        // //! [ 1.e-3 factor because dens_aer is (g/cm3), specvol is (m3/kg) ]
        specvol=1.0e-3/dens_aer(l,n);
        ////////////////////
        vaerosol(m,n,i)=vaerosol(m,n,i)+maerosol(l,m,n,i)*specvol;
        vaerosolcw(m,n,i)=vaerosolcw(m,n,i)+maerosolcw(l,m,n,i)*specvol;
    }
    if(lnum > 0){
        // !            aerosol number predicted
        // ! [ 1.0e6 factor because because dhi_ & dlo_sect are (cm), vaerosol is (m3) ]

        npv_at_dhi = 6.0e6/(pi*dhi_sect(m,n)*dhi_sect(m,n)*dhi_sect(m,n));
        npv_at_dlo = 6.0e6/(pi*dlo_sect(m,n)*dlo_sect(m,n)*dlo_sect(m,n));

        naerosol(m,n,i)=chem(k,lnum-1,0,i)*cs;
        naerosolcw(m,n,i)=chem(k,lnumcw-1,0,i)*cs;
        num_at_dhi = vaerosol(m,n,i)*npv_at_dhi;
        num_at_dlo = vaerosol(m,n,i)*npv_at_dlo;

        naerosol(m,n,i) = Kokkos::max( num_at_dhi, Kokkos::min( num_at_dlo, naerosol(m,n,i) ) );
        num_at_dhi = vaerosolcw(m,n,i)*npv_at_dhi;
        num_at_dlo = vaerosolcw(m,n,i)*npv_at_dlo;
        naerosolcw(m,n,i) = Kokkos::max( num_at_dhi, Kokkos::min( num_at_dlo, naerosolcw(m,n,i) ) );
    }
    else{
        naerosol(m,n,i)=vaerosol(m,n,i)*npv(m,n,i);
        naerosol(m,n,i)=Kokkos::max(naerosol(m,n,i),0.);
        naerosolcw(m,n,i)=vaerosolcw(m,n,i)*npv(m,n,i);
        naerosolcw(m,n,i)=Kokkos::max(naerosolcw(m,n,i),0.);
    }
}

template <class Aer2View, class FluxFullView, class WorkReal2View,
          class WorkReal1View, class WorkBool2View, class WorkInt1View>
KOKKOS_INLINE_FUNCTION void activate(real wbar, real sigw, real wdiab, real wminf, real wmaxf, real tair, real rhoair,
                      int msectional, int maxd_atype, int ntype_aer, int maxd_asize, const int1d &nsize_aer,
                      const Aer2View &na, const Aer2View &volc,const real2d &dlo_sect,const real2d &dhi_sect,const real2d &sigman, const Aer2View &hygro,
                      const Aer2View &fn, const Aer2View &fs, const Aer2View &fm, const Aer2View &fluxn, const Aer2View &fluxs, const Aer2View &fluxm, const FluxFullView &flux_fullact,
                      const WorkReal2View &ylo, const WorkReal2View &yhi, const WorkReal2View &ymean,
                      const WorkReal2View &alnsign, const WorkReal2View &am, const WorkReal2View &lnhygro,
                      const WorkReal2View &f1, const WorkReal2View &sm, const WorkReal2View &lnsm,
                      const WorkReal2View &eta, const WorkReal2View &asub, const WorkReal2View &bsub,
                      const WorkReal2View &sumflxn, const WorkReal2View &sumflxs, const WorkReal2View &sumflxm,
                      const WorkReal2View &sumfn, const WorkReal2View &sumfs, const WorkReal2View &sumfm,
                      const WorkReal2View &fnold, const WorkReal2View &fsold, const WorkReal2View &fmold,
                      const WorkReal2View &rlo, const WorkReal2View &rhi, const WorkReal2View &rmean,
                      const WorkReal2View &etafactor2,
                      const WorkReal1View &totn, const WorkReal1View &gmrad,
                      const WorkReal1View &gmradsq, const WorkReal1View &gmlnsig,
                      const WorkReal1View &gmsm, const WorkReal1View &sumns,
                      const WorkBool2View &bin_is_empty, const WorkBool2View &bin_is_narrow,
                      const WorkInt1View &tmp_nsize_aer,
                      int grid_id, int ktau, int ii,  int kk,real smax_prescribed){

    //loacl

    int nx=200;
    int iquasisect_option, isectional;
    real integ,integf;
    real surften = 0.076;                                   //! surface tension of water w/respect to air (N/m)
    real p0 = 1013.25e2;                                    //! reference pressure (Pa)
    real t0 = 273.15;                                       //! reference temperature (K)
    real ycut, lnycut, betayy, betayy2, gammayy, phiyy;

    real pres;          //! pressure (Pa)
    real path;          //! mean free path (m)
    real diff;          //! diffusivity (m2/s)
    real conduct;       //! thermal conductivity (Joule/m/sec/deg)
    real diff0,conduct0;
    real es;            //! saturation vapor pressure
    real qs;            //! water vapor saturation mixing ratio
    real dqsdt;         //! change in qs with temperature
    real dqsdp;         //! change in qs with pressure
    real gg;            //! thermodynamic function (m2/s)
    real sqrtg;         //! Kokkos::sqrt(gg)

    real zeta;

    real lnsmax;        //! ln(smax)
    real alpha;
    real gamma;
    real beta;
    real gaus;
    bool top;           //! true if cloud top, false if cloud base or new cloud

    real aten;                                          //! surface tension parameter
    real sumflx_fullact;
    real wold,gold;
    real alogten,alog2,alog3,alogaten;
    real alogam;

    real ccc;
    real dumaa,dumbb;
    real wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb;
    real dfmin,dfmax,fnew,fold,fnmin,fnbar,fsbar,fmbar;
    real alw,sqrtalw;
    real smax;
    real x,arg;
    real xmincoeff,xcut;
    real z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf;
    real etafactor1,etafactor2max;
    int m,n,nw,nwmax;

    bool lable3000 = true;
    bool lable4000 = true;
    bool lable4100 = true;
    bool lable4800 = true;
    bool lable5000 = true;
    bool lable6000 = true;

    real ep_2 = 287.0/461.6;
    real svp1 = 0.6112;
    real svp2 = 17.67;
    real svp3 = 29.65;


    //!      numerical integration parameters
    real eps = 0.3;
    real fmax = 0.99;
    real sds = 3.;

    //!      mathematical constants
    real third, twothird, sixth, zero, one, two, three;
    real sq2  = 1.4142135624;
    real sqpi = 1.7724538509;
    real pi   = 3.1415926536;

    // ! rce 08-jul-2005
    // ! if either (na(n,m,ii) < nsmall) or (volc(n,m,ii) < vsmall)
    // ! then treat bin/mode (n,m) as being empty, and set its fn/fs/fm=0.0
    // !     (for single precision, gradual underflow starts around 1.0e-38,
    // !      and strange things can happen when in that region)
    real nsmall = 1.0e-20;                                          //! aer number conc in #/m3
    real vsmall = 1.0e-37;                                          //! aer volume conc in m3/m3
    bool all_bins_empty;
    int idiagaa, ipass_nwloop;
    int idiag_dndy_neg, idiag_fnsm_prob;


    // ! The flag for cloud top is no longer used so set it to false. This is an
    // ! antiquated feature related to radiation enhancing mass fluxes at cloud
    // ! top. It is currently, as of Feb. 2009, set to false in the CAM version
    // ! as well.
    top = false;

    // !.......................................................................
    // !
    // !   start calc. of modal or sectional activation properties (start of section 1)
    // !
    // !.......................................................................
    idiag_dndy_neg = 1;         //! set this to 0 to turn off   warnings about dn/dy < 0
    idiag_fnsm_prob = 1;        //! set this to 0 to turn off   warnings about fn/fs/fm misbehavior
    iquasisect_option = 2;

    if(msectional > 0) isectional = iquasisect_option;
    else isectional = 0;


    //!BSINGH - For WRFCuP
    // if ( present( smax_prescribed ) )  smax_prescribed = 1 is no_present
    if(smax_prescribed != 1.0){
        if(smax_prescribed <= 0.0){
            for(int n=0; n < ntype_aer; n++){
                for(int m=0; m < nsize_aer(n); m++){
                    fluxn(m,n,ii)=0.;
                    fn(m,n,ii)=0.;
                    fluxs(m,n,ii)=0.;
                    fs(m,n,ii)=0.;
                    fluxm(m,n,ii)=0.;
                    fm(m,n,ii)=0.;
                }
            }
            flux_fullact(kk,ii)=0.;
            return;
        }
    }

    // //!BSINGH - ENDS

    for(int n = 0; n < ntype_aer; n++){
       if(ntype_aer == 1 && nsize_aer(n) == 1 && na(0,0,ii) < 1.e-20){
            fn(0,0,ii)=0.;
            fs(0,0,ii)=0.;
            fm(0,0,ii)=0.;
            fluxn(0,0,ii)=0.;
            fluxs(0,0,ii)=0.;
            fluxm(0,0,ii)=0.;
            flux_fullact(kk,ii)=0.;
            return;
       }
    }

    zero = 0.0;
    one = 1.0;
    two = 2.0;
    three = 3.0;
    third = 1.0/3.0;
    twothird = 2.0/3.0; //!wig, 1-Mar-2009: Corrected value from 2/6
    sixth = 1.0/6.0;


    pres=r_d*rhoair*tair;
    diff0=0.211e-4*(p0/pres)*Kokkos::pow((tair/t0),1.94);
    conduct0=(5.69+0.017*(tair-t0))*4.186e2*1.e-5;  //! convert to J/m/s/deg
    es=1000.*svp1*Kokkos::exp( svp2*(tair-t0)/(tair-svp3) );
    qs=ep_2*es/(pres-es);
    dqsdt=xlv/(r_v*tair*tair)*qs;
    alpha=g*(xlv/(cp*r_v*tair*tair)-1./(r_d*tair));
    gamma=(1+xlv/cp*dqsdt)/(rhoair*qs);
    gg=1./(rhowater/(diff0*rhoair*qs)+xlv*rhowater/(conduct0*tair)*(xlv/(r_v*tair)-1.));
    sqrtg=Kokkos::sqrt(gg);
    beta=4.*pi*rhowater*gg*gamma;
    aten=2.*surften/(r_v*tair*rhowater);
    alogaten=Kokkos::log(aten);
    alog2=Kokkos::log(two);
    alog3=Kokkos::log(three);
    ccc=4.*pi*third;
    etafactor2max=1.e10/Kokkos::pow(alpha*wmaxf,1.5); //! this should make eta big if na is very small.


    all_bins_empty = true;
    for(int n = 0; n < ntype_aer; n++){
        totn(n,ii)=0.;
        gmrad(n,ii)=0.;
        gmradsq(n,ii)=0.;
        sumns(n,ii)=0.;
        for(int m = 0; m < nsize_aer(n); m++){
            alnsign(m,n,ii)=Kokkos::log(sigman(m,n));

            bin_is_empty(m,n,ii) = true;

            if (volc(m,n,ii) > vsmall && na(m,n,ii) > nsmall){
                bin_is_empty(m,n,ii) = false;
                all_bins_empty = false;
                lnhygro(m,n,ii)=Kokkos::log(hygro(m,n,ii));
                am(m,n,ii)=Kokkos::exp(-1.5*alnsign(m,n,ii)*alnsign(m,n,ii))*
                    Kokkos::pow(3.*volc(m,n,ii)/(4.*pi*na(m,n,ii)),third);

                if (isectional > 0){
                    totn(n,ii)=totn(n,ii)+na(m,n,ii);
                    alogam=Kokkos::log(am(m,n,ii));
                    gmrad(n,ii)=gmrad(n,ii)+na(m,n,ii)*alogam;
                    gmradsq(n,ii)=gmradsq(n,ii)+na(m,n,ii)*alogam*alogam;
                }
                etafactor2(m,n,ii)=1./(na(m,n,ii)*beta*sqrtg);

                if(hygro(m,n,ii) > 1.e-10) sm(m,n,ii)=2.*aten/(3.*am(m,n,ii))*Kokkos::sqrt(aten/(3.*hygro(m,n,ii)*am(m,n,ii)));
                else sm(m,n,ii)=100.;
            }
            else{
                sm(m,n,ii)=1.;
                etafactor2(m,n,ii)=etafactor2max;
            }

            lnsm(m,n,ii)=Kokkos::log(sm(m,n,ii));
            if ((isectional == 3) || (isectional == 4))
                sumns(n,ii)=sumns(n,ii)+na(m,n,ii)/Kokkos::pow(sm(m,n,ii),twothird);
        }   //size
    }       //type

    //!  if all bins are empty, set all activation fractions to zero and exit
    if ( all_bins_empty ){
        for(int n=0; n<ntype_aer; n++){
            for(int m=0; m<nsize_aer(n); m++){
                fluxn(m,n,ii)=0.;
                fn(m,n,ii)=0.;
                fluxs(m,n,ii)=0.;
                fs(m,n,ii)=0.;
                fluxm(m,n,ii)=0.;
                fm(m,n,ii)=0.;
            }
        }
        flux_fullact(kk,ii)=0.;
        return;
    }

    if (isectional <= 0){
        maxsat_init(maxd_atype, ntype_aer,
            maxd_asize, nsize_aer, alnsign, f1,ii);
        //goto 30000
        lable3000 = false;
    }
    if(lable3000 == true){
        for(int n = 0; n < ntype_aer; n++){
            gmrad(n,ii)=gmrad(n,ii)/Kokkos::max(totn(n,ii),1e-20);
            for(int nn = 0; nn < maxd_atype; nn++)
                gmlnsig(nn,ii)=gmradsq(n,ii)/totn(n,ii)-gmrad(n,ii)*gmrad(n,ii);
            gmlnsig(n,ii)=Kokkos::sqrt( Kokkos::max( 1.e-4, gmlnsig(n,ii) ) );
            gmrad(n,ii)=Kokkos::exp(gmrad(n,ii));
            if ((isectional == 3) || (isectional == 4)){
                gmsm(n,ii)=totn(n,ii)/sumns(n,ii);
                gmsm(n,ii)=gmsm(n,ii)*gmsm(n,ii)*gmsm(n,ii);
                gmsm(n,ii)=Kokkos::sqrt(gmsm(n,ii));
            }
            else gmsm(n,ii)=2.*aten/(3.*gmrad(n,ii))*Kokkos::sqrt(aten/(3.*hygro(nsize_aer(n)-1,n,ii)*gmrad(n,ii)));
        }
        for(int max = 0; max < maxd_atype; max++){
            tmp_nsize_aer(max,ii) = 1;
        }

        maxsat_init_1(maxd_atype, ntype_aer,
              maxd_asize, tmp_nsize_aer, gmlnsig, f1,ii);

        for(int n = 0; n < ntype_aer; n++){
            for(int m = 0; m < nsize_aer(n); m++){
                rlo(m,n,ii) = 0.5*0.01*dlo_sect(m,n);
                rhi(m,n,ii) = 0.5*0.01*dhi_sect(m,n);
                ylo(m,n,ii) = Kokkos::pow(rlo(m,n,ii)/rhi(m,n,ii),3);
                yhi(m,n,ii) = 1.0;
                bin_is_narrow(m,n,ii) = false;
                if ((rhi(m,n,ii)/rlo(m,n,ii)) <= 1.01) bin_is_narrow(m,n,ii) = true;
                if ( bin_is_empty(m,n,ii) ){
                    rmean(m,n,ii) = Kokkos::sqrt(rlo(m,n,ii)*rhi(m,n,ii));
                    ymean(m,n,ii) = Kokkos::pow(rmean(m,n,ii)/rhi(m,n,ii),3);
                    continue;
                }
                rmean(m,n,ii) = Kokkos::pow(volc(m,n,ii)/(ccc*na(m,n,ii)),third);
                rmean(m,n,ii) = Kokkos::max( rlo(m,n,ii), Kokkos::min( rhi(m,n,ii), rmean(m,n,ii) ) );
                ymean(m,n,ii) = Kokkos::pow(rmean(m,n,ii)/rhi(m,n,ii),3);
                if ( bin_is_narrow(m,n,ii) ) continue;
                if ((rhi(m,n,ii)/rmean(m,n,ii)) <= 1.01){
                    bin_is_narrow(m,n,ii) = true;
                    rlo(m,n,ii) = Kokkos::min( rmean(m,n,ii), (rhi(m,n,ii)/1.01) );
                    ylo(m,n,ii) = Kokkos::pow(rlo(m,n,ii)/rhi(m,n,ii),3);
                    continue;
                }
                else if ((rmean(m,n,ii)/rlo(m,n,ii)) <= 1.01){
                    bin_is_narrow(m,n,ii) = true;
                    rhi(m,n,ii) = Kokkos::max( rmean(m,n,ii), (rlo(m,n,ii)*1.01) );
                    ylo(m,n,ii) = Kokkos::pow(rlo(m,n,ii)/rhi(m,n,ii),3);
                    ymean(m,n,ii) = Kokkos::pow(rmean(m,n,ii)/rhi(m,n,ii),3);
                    continue;
                }
                // ! if rmean is somewhat close to either rlo or rhi, then dn/dy will be
                // !    negative near the upper or lower bin boundary
                // ! in these cases, assume that all the particles are in a subset of the full bin,
                // !    and adjust rlo or rhi so that rmean will be near the center of this subset
                // ! note that the bin is made narrower LOCALLY/TEMPORARILY,
                // !    just for the purposes of the activation calculation

                gammayy = (ymean(m,n,ii)-ylo(m,n,ii)) / (yhi(m,n,ii)-ylo(m,n,ii));
                if (gammayy < 0.34){
                   dumaa = ylo(m,n,ii) + (yhi(m,n,ii)-ylo(m,n,ii))*(gammayy/0.34);
                   rhi(m,n,ii) = rhi(m,n,ii)*Kokkos::pow(dumaa,third);
                   ylo(m,n,ii) = Kokkos::pow(rlo(m,n,ii)/rhi(m,n,ii),3);
                   ymean(m,n,ii) = Kokkos::pow(rmean(m,n,ii)/rhi(m,n,ii),3);
                }
                else if (gammayy >= 0.66){
                   dumaa = ylo(m,n,ii) + (yhi(m,n,ii)-ylo(m,n,ii))*((gammayy-0.66)/0.34);
                   ylo(m,n,ii) = dumaa;
                   rlo(m,n,ii) = rhi(m,n,ii)*Kokkos::pow(dumaa,third);
                }
                if ((rhi(m,n,ii)/rlo(m,n,ii)) <= 1.01){
                   bin_is_narrow(m,n,ii) = true;
                   continue;
                }

                betayy = ylo(m,n,ii)/yhi(m,n,ii);
                betayy2 = betayy*betayy;
                bsub(m,n,ii) = (12.0*ymean(m,n,ii) - 6.0*(1.0+betayy)) /(4.0*(1.0-betayy2*betayy) - 3.0*(1.0-betayy2)*(1.0+betayy));
                asub(m,n,ii) = (1.0 - bsub(m,n,ii)*(1.0-betayy2)*0.5) / (1.0-betayy);

                // if ( asub(m,n,ii)+bsub(m,n,ii)*ylo(m,n,ii) < 0. ){
                //     if (idiag_dndy_neg > 0){
                //         print *,'dndy<0 at lower boundary'
                //         print *,'n,m=',n,m
                //         print *,'na=',na(m,n,ii),' volc=',volc(m,n,ii)
                //         print *,'volc/(na*pi*4/3)=', (volc(m,n,ii)/(na(m,n,ii)*ccc))
                //         print *,'rlo(m,n,ii),rhi(m,n,ii)=',rlo(m,n,ii),rhi(m,n,ii)
                //         print *,'dlo_sect/2,dhi_sect/2=',   &
                //                 (0.005*dlo_sect(m,n)),(0.005*dhi_sect(m,n))
                //         print *,'asub,bsub,ylo,yhi=',asub(m,n,ii),bsub(m,n,ii),ylo(m,n,ii),yhi(m,n,ii)
                //         print *,'asub+bsub*ylo=',   &
                //                 (asub(m,n,ii)+bsub(m,n,ii)*ylo(m,n,ii))
                //         print *,'subr activate error 11 - i,k =', ii,  kk
                //     }
                // }

                // if ( asub(m,n,ii)+bsub(m,n,ii)*yhi(m,n,ii) .lt. 0. ) then
                //   if (idiag_dndy_neg .gt. 0) then
                //     print *,'dndy<0 at upper boundary'
                //     print *,'n,m=',n,m
                //     print *,'na=',na(m,n,ii),' volc=',volc(m,n,ii)
                //     print *,'volc/(na*pi*4/3)=', (volc(m,n,ii)/(na(m,n,ii)*ccc))
                //     print *,'rlo(m,n,ii),rhi(m,n,ii)=',rlo(m,n,ii),rhi(m,n,ii)
                //     print *,'dlo_sect/2,dhi_sect/2=',   &
                //              (0.005*dlo_sect(m,n)),(0.005*dhi_sect(m,n))
                //     print *,'asub,bsub,ylo,yhi=',asub(m,n,ii),bsub(m,n,ii),ylo(m,n,ii),yhi(m,n,ii)
                //     print *,'asub+bsub*yhi=',   &
                //              (asub(m,n,ii)+bsub(m,n,ii)*yhi(m,n,ii))
                //     print *,'subr activate error 12 - i,k =', ii, kk
                //   endif
                // endif
            }
        }
    }   //lable3000

    if(sigw <= 1.e-5) lable5000 = false;
    if(lable5000 == true){
        ipass_nwloop = 1;
        idiagaa = 0;

        while(lable4000){
            lable4000 = false;
            if(top){
                wmax=0.;
                wmin=Kokkos::min(zero,-wdiab);
            }
            else{
                wmax=Kokkos::min(wmaxf,wbar+sds*sigw);
                wmin=Kokkos::max(wminf,-wdiab);
            }
            wmin=Kokkos::max(wmin,wbar-sds*sigw);
            w=wmin;
            dwmax=eps*sigw;
            dw=dwmax;
            dfmax=0.2;
            dfmin=0.1;
            if(wmax<=w){
                for(int n=0; n<ntype_aer; n++){
                    for(int m=0; m<nsize_aer(n); m++){
                        fluxn(m,n,ii)=0.;
                        fn(m,n,ii)=0.;
                        fluxs(m,n,ii)=0.;
                        fs(m,n,ii)=0.;
                        fluxm(m,n,ii)=0.;
                        fm(m,n,ii)=0.;
                    }
                }
                flux_fullact(kk,ii)=0.;
                return;
            }
            for(int n=0; n<ntype_aer; n++){
                for(int m=0; m<nsize_aer(n); m++){
                    sumflxn(m,n,ii)=0.;
                    sumfn(m,n,ii)=0.;
                    fnold(m,n,ii)=0.;
                    sumflxs(m,n,ii)=0.;
                    sumfs(m,n,ii)=0.;
                    fsold(m,n,ii)=0.;
                    sumflxm(m,n,ii)=0.;
                    sumfm(m,n,ii)=0.;
                    fmold(m,n,ii)=0.;
                }
            }
            sumflx_fullact=0.;
            fold=0;
            gold=0;
            wold=w;

            nwmax = 200;
            dwmin = (wmax - wmin)/(nwmax-1);
            dwmin = Kokkos::min( dwmax, dwmin );
            dwmin = Kokkos::max( 0.01,  dwmin );

            // if (idiagaa.gt.0) then
            //  write(*,94700) ktau, grid_id, ii, kk, nwmax
            //  write(*,94710) 'wbar,sigw,wdiab=', wbar, sigw, wdiab
            //  write(*,94710) 'wmin,wmax,dwmin,dwmax=', wmin, wmax, dwmin, dwmax
            //  write(*,94720) -1, w, wold, dw
            // endif
            // 94700    format( / 'activate 47000 - ktau,id,ii,jj,kk,nwmax=', 6i5 )
            // 94710    format( 'activate 47000 - ', a, 6(1x,f11.5) )
            // 94720    format( 'activate 47000 - nw,w,wold,dw=', i5, 3(1x,f11.5) )

            for(int nw = 0; nw < nwmax; nw++){              //lable 47000
                lable4100 = true;
                while (lable4100){
                    lable4100 = false;

                    wnuc=w+wdiab;                               //lable 41000

                    // if (idiagaa.gt.0) write(*,94720) nw, w, wold, dw

                    alw=alpha*wnuc;
                    sqrtalw=Kokkos::sqrt(alw);
                    zeta=2.*sqrtalw*aten/(3.*sqrtg);
                    etafactor1=2.*alw*sqrtalw;

                    if (isectional > 0){
                        for(int n = 0;n < ntype_aer; n++){
                            if(totn(n,ii) > 1.e-10) eta(0,n,ii)=etafactor1/(totn(n,ii)*beta*sqrtg);
                            else eta(0,n,ii)=1.e10;
                        }

                        if (smax_prescribed != 1.0) smax = smax_prescribed;
                        else {
                            for(int max = 0; max < maxd_atype; max++){
                                tmp_nsize_aer(max,ii) = 1;
                            }
                            maxsat_1(zeta,eta,maxd_atype,ntype_aer,
                              maxd_asize,tmp_nsize_aer,gmsm,gmlnsig,f1,smax,ii);
                        }

                        lnsmax=Kokkos::log(smax);
                        x=2*(Kokkos::log(gmsm(0,ii))-lnsmax)/(3*sq2*gmlnsig(0,ii));
                        fnew=0.5*(1.-erf_alt(x));
                    }
                    else{
                        for(int n = 0; n < ntype_aer; n++){
                            for(int m = 0; m < nsize_aer(n); m++){
                                eta(m,n,ii)=etafactor1*etafactor2(m,n,ii);
                            }
                        }

                        if (smax_prescribed != 1.0) smax = smax_prescribed;
                        else{
                            maxsat(zeta,eta,maxd_atype,ntype_aer,
                                maxd_asize,nsize_aer,sm,alnsign,f1,smax,ii);
                        }

                        lnsmax=Kokkos::log(smax);
                        x=2*(lnsm(nsize_aer(0)-1,0,ii)-lnsmax)/(3*sq2*alnsign(nsize_aer(0)-1,0,ii));
                        fnew=0.5*(1.-erf_alt(x));
                    }

                    dwnew = dw;
                    if(fnew-fold > dfmax && nw > 0){
                        if (dw > 1.01*dwmin){
                            dw=0.7*dw;
                            dw=Kokkos::max(dw,dwmin);
                            w=wold+dw;
                            lable4100 = true;
                            continue;
                        }
                        else dwnew = dwmin;
                    }

                    if(fnew-fold < dfmin) dwnew=Kokkos::min(1.5*dw,dwmax);

                    fold=fnew;
                    z=(w-wbar)/(sigw*sq2);
                    gaus=Kokkos::exp(-z*z);
                    fnmin=1.;
                    xmincoeff=alogaten-2.*third*(lnsmax-alog2)-alog3;

                    for(int n = 0; n < ntype_aer; n++){
                        for(int m = 0; m < nsize_aer(n); m++){
                            if ( bin_is_empty(m,n,ii) ){
                                fn(m,n,ii)=0.;
                                fs(m,n,ii)=0.;
                                fm(m,n,ii)=0.;
                            }
                            else if ((isectional == 2) || (isectional == 4)){
                                xcut=xmincoeff-third*lnhygro(m,n,ii);
                                lnycut = 3.0 * ( xcut - Kokkos::log(rhi(m,n,ii)) );
                                lnycut = Kokkos::min( lnycut, Kokkos::log(yhi(m,n,ii)*1.0e5) );
                                ycut=Kokkos::exp(lnycut);

                                if(ycut > yhi(m,n,ii)){
                                    fn(m,n,ii)=0.;
                                    fs(m,n,ii)=0.;
                                    fm(m,n,ii)=0.;
                                }
                                else if(ycut < ylo(m,n,ii)){
                                    fn(m,n,ii)=1.;
                                    fs(m,n,ii)=1.;
                                    fm(m,n,ii)=1.;
                                }
                                else if ( bin_is_narrow(m,n,ii) ){
                                    if (ycut > ymean(m,n,ii)){
                                        fn(m,n,ii)=0.;
                                        fs(m,n,ii)=0.;
                                        fm(m,n,ii)=0.;
                                    }
                                    else{
                                        fn(m,n,ii)=1.;
                                        fs(m,n,ii)=1.;
                                        fm(m,n,ii)=1.;
                                    }
                                }
                                else{
                                    phiyy=ycut/yhi(m,n,ii);
                                    fn(m,n,ii) = asub(m,n,ii)*(1.0-phiyy) + 0.5*bsub(m,n,ii)*(1.0-phiyy*phiyy);
                                    // if (fn(m,n,ii).lt.zero .or. fn(m,n,ii).gt.one) then
                                    // if (idiag_fnsm_prob .gt. 0) then
                                    //     print *,'fn(',m,n,',ii)=',fn(m,n,ii),' outside 0,1 - activate err21'
                                    //     print *,'na,volc       =', na(m,n,ii), volc(m,n,ii)
                                    //     print *,'asub,bsub     =', asub(m,n,ii), bsub(m,n,ii)
                                    //     print *,'yhi,ycut      =', yhi(m,n,ii), ycut
                                    // endif
                                    // endif

                                    if (fn(m,n,ii) <= zero){
                                        fn(m,n,ii)=zero;
                                        fs(m,n,ii)=zero;
                                        fm(m,n,ii)=zero;
                                    }
                                    else if (fn(m,n,ii) >= one){
                                        fn(m,n,ii)=one;
                                        fs(m,n,ii)=one;
                                        fm(m,n,ii)=one;
                                    }
                                    else{
                                        fm(m,n,ii) = (yhi(m,n,ii)/ymean(m,n,ii)) * (0.5*asub(m,n,ii)*(1.0-phiyy*phiyy) + third*bsub(m,n,ii)*(1.0-phiyy*phiyy*phiyy));
                                        // if (fm(m,n,ii).lt.fn(m,n,ii) .or. fm(m,n,ii).gt.one) then
                                        //     if (idiag_fnsm_prob .gt. 0) then
                                        //         print *,'fm(',m,n,',ii)=',fm(m,n,ii),' outside fn,1 - activate err22'
                                        //         print *,'na,volc,fn    =', na(m,n,ii), volc(m,n,ii), fn(m,n,ii)
                                        //         print *,'asub,bsub     =', asub(m,n,ii), bsub(m,n,ii)
                                        //         print *,'yhi,ycut     =', yhi(m,n,ii), ycut
                                        //     endif
                                        // endif
                                        if (fm(m,n,ii) <= fn(m,n,ii)){
                                            fm(m,n,ii)=fn(m,n,ii);
                                            fs(m,n,ii)=fn(m,n,ii);
                                        }
                                        else if (fm(m,n,ii) >= one){
                                            fm(m,n,ii)=one;
                                            fs(m,n,ii)=one;
                                            fn(m,n,ii)=one;
                                        }
                                        else{
                                            dumaa = fn(m,n,ii)*(yhi(m,n,ii)/ymean(m,n,ii));
                                            fm(m,n,ii) = Kokkos::min( fm(m,n,ii), dumaa );
                                            dumaa = 1.0 + (fn(m,n,ii)-1.0)*(ylo(m,n,ii)/ymean(m,n,ii)) ;
                                            fm(m,n,ii) = Kokkos::min( fm(m,n,ii), dumaa );
                                            //! 10-nov-2005 rce - now calculate fs and bound it by fn, fm
                                            betayy = ylo(m,n,ii)/yhi(m,n,ii);
                                            dumaa = Kokkos::pow(phiyy,twothird);
                                            dumbb = Kokkos::pow(betayy,twothird);
                                            fs(m,n,ii) = (asub(m,n,ii)*(1.0-phiyy*dumaa) +
                                                    0.625*bsub(m,n,ii)*(1.0-phiyy*phiyy*dumaa)) /
                                                    (asub(m,n,ii)*(1.0-betayy*dumbb) +
                                                    0.625*bsub(m,n,ii)*(1.0-betayy*betayy*dumbb));
                                            fs(m,n,ii)=Kokkos::max(fs(m,n,ii),fn(m,n,ii));
                                            fs(m,n,ii)=Kokkos::min(fs(m,n,ii),fm(m,n,ii));
                                        }
                                    }
                                }
                            }
                            else{
                                x=2*(lnsm(m,n,ii)-lnsmax)/(3*sq2*alnsign(m,n,ii));
                                fn(m,n,ii)=0.5*(1.-erf_alt(x));
                                arg=x-sq2*alnsign(m,n,ii);
                                fs(m,n,ii)=0.5*(1.-erf_alt(arg));
                                arg=x-1.5*sq2*alnsign(m,n,ii);
                                fm(m,n,ii)=0.5*(1.-erf_alt(arg));
                            }

                            fnmin=Kokkos::min(fn(m,n,ii),fnmin);
                            // !               integration is second order accurate
                            // !               assumes linear variation of f*gaus with w
                            wb=(w+wold);
                            fnbar=(fn(m,n,ii)*gaus+fnold(m,n,ii)*gold);
                            fsbar=(fs(m,n,ii)*gaus+fsold(m,n,ii)*gold);
                            fmbar=(fm(m,n,ii)*gaus+fmold(m,n,ii)*gold);
                            if((top && w < 0.) || (! top && w > 0.)){
                                sumflxn(m,n,ii)=sumflxn(m,n,ii)+sixth*(wb*fnbar
                                +(fn(m,n,ii)*gaus*w+fnold(m,n,ii)*gold*wold))*dw;
                                sumflxs(m,n,ii)=sumflxs(m,n,ii)+sixth*(wb*fsbar
                                    +(fs(m,n,ii)*gaus*w+fsold(m,n,ii)*gold*wold))*dw;
                                sumflxm(m,n,ii)=sumflxm(m,n,ii)+sixth*(wb*fmbar
                                    +(fm(m,n,ii)*gaus*w+fmold(m,n,ii)*gold*wold))*dw;
                            }
                            sumfn(m,n,ii)=sumfn(m,n,ii)+0.5*fnbar*dw;
                            fnold(m,n,ii)=fn(m,n,ii);
                            sumfs(m,n,ii)=sumfs(m,n,ii)+0.5*fsbar*dw;
                            fsold(m,n,ii)=fs(m,n,ii);
                            sumfm(m,n,ii)=sumfm(m,n,ii)+0.5*fmbar*dw;
                            fmold(m,n,ii)=fm(m,n,ii);
                        }//end m
                    }//end n

                    sumflx_fullact = sumflx_fullact + sixth*(wb*(gaus+gold) + (gaus*w + gold*wold))*dw;
                    gold=gaus;
                    wold=w;
                    dw=dwnew;
                    if(nw > 0 && (w > wmax || fnmin > fmax)){
                        lable4800 = false;
                        nw = nwmax;
                        break;
                    }
                    w=w+dw;
                }
            }//end 47000
            if(lable4800){
                // write(0,*)'do loop is too short in activate'
                // print *,'do loop is too short in activate'
                // write(0,*)'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
                // print *,'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
                // write(0,*)'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
                // print *,'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
                // write(0,*)'wnuc=',wnuc
                // print *,'wnuc=',wnuc
                // do n=1,ntype_aer
                //     write(0,*)'ntype=',n
                //     print *,'ntype=',n
                //     write(0,*)'na=',(na(m,n,ii),m=1,nsize_aer(n))
                //     print *,'na=',(na(m,n,ii),m=1,nsize_aer(n))
                //     write(0,*)'fn=',(fn(m,n,ii),m=1,nsize_aer(n))
                //     print *,'fn=',(fn(m,n,ii),m=1,nsize_aer(n))
                // end do
                // !   dump all subr parameters to allow testing with standalone code
                // !   (build a driver that will read input and call activate)
                // write(0,*)'top,wbar,sigw,wdiab,tair,rhoair,ntype_aer='
                // print *,'top,wbar,sigw,wdiab,tair,rhoair,ntype_aer='
                // write(0,*) top,wbar,sigw,wdiab,tair,rhoair,ntype_aer
                // print *, top,wbar,sigw,wdiab,tair,rhoair,ntype_aer
                // write(0,*)'na='
                // print *,'na='
                // write(0,*) na
                // print *, na
                // write(0,*)'volc='
                // print *,'volc='
                // write(0,*) volc
                // print *, volc
                // write(0,*)'sigman='
                // print *,'sigman='
                // write(0,*) sigman
                // print *, sigman
                // write(0,*)'hygro='
                // print *,'hygro='
                // write(0,*) hygro
                // print *, hygro

                // write(0,*)'subr activate error 31 - i,k =', ii,  kk
                // print *,'subr activate error 31 - i,k =', ii,  kk
                if (ipass_nwloop == 1){
                    ipass_nwloop = 2;
                    idiagaa = 2;
                    lable4000 = true;
                    continue;
                }
                Kokkos::abort("STOP: activate before 48000");
            }
        }//end 4000

        if(!top && w < wmaxf){
            wnuc=w+wdiab;

            z1=(w-wbar)/(sigw*sq2);
            z2=(wmaxf-wbar)/(sigw*sq2);
            integ=sigw*0.5*sq2*sqpi*(erfc_num_recipes(z1)-erfc_num_recipes(z2));

            //!            consider only upward flow into cloud base when estimating flux
            wf1=Kokkos::max(w,zero);
            zf1=(wf1-wbar)/(sigw*sq2);
            gf1=Kokkos::exp(-zf1*zf1);
            wf2=Kokkos::max(wmaxf,zero);
            zf2=(wf2-wbar)/(sigw*sq2);
            gf2=Kokkos::exp(-zf2*zf2);
            gf=(gf1-gf2);
            integf=wbar*sigw*0.5*sq2*sqpi*(erfc_num_recipes(zf1)-erfc_num_recipes(zf2))+sigw*sigw*gf;

            for(int n = 0; n < ntype_aer; n++){
                for(int m = 0; m < nsize_aer(n); m++){
                    sumflxn(m,n,ii)=sumflxn(m,n,ii)+integf*fn(m,n,ii);
                    sumfn(m,n,ii)=sumfn(m,n,ii)+fn(m,n,ii)*integ;
                    sumflxs(m,n,ii)=sumflxs(m,n,ii)+integf*fs(m,n,ii);
                    sumfs(m,n,ii)=sumfs(m,n,ii)+fs(m,n,ii)*integ;
                    sumflxm(m,n,ii)=sumflxm(m,n,ii)+integf*fm(m,n,ii);
                    sumfm(m,n,ii)=sumfm(m,n,ii)+fm(m,n,ii)*integ;
                }
            }

            sumflx_fullact = sumflx_fullact + integf;
        }

        for(int n = 0; n < ntype_aer; n++){
            for(int m = 0; m < nsize_aer(n); m++){
                fn(m,n,ii)=sumfn(m,n,ii)/(sq2*sqpi*sigw);
                fluxn(m,n,ii)=sumflxn(m,n,ii)/(sq2*sqpi*sigw);
                if(fn(m,n,ii) > 1.01){
                    // if (idiag_fnsm_prob > 0){
                    //     print *,'fn=',fn(m,n,ii),' > 1 - activate err41'
                    //     print *,'w,m,n,na,am=',w,m,n,na(m,n,ii),am(m,n,ii)
                    //     print *,'integ,sumfn,sigw=',integ,sumfn(m,n,ii),sigw
                    //     print *,'subr activate error - i,k =', ii,  kk
                    // }
                    fluxn(m,n,ii) = fluxn(m,n,ii)/fn(m,n,ii);
                }

                fs(m,n,ii)=sumfs(m,n,ii)/(sq2*sqpi*sigw);
                fluxs(m,n,ii)=sumflxs(m,n,ii)/(sq2*sqpi*sigw);
                if(fs(m,n,ii) > 1.01){
                    // if (idiag_fnsm_prob .gt. 0) then
                    //     print *,'fs=',fs(m,n,ii),' > 1 - activate err42'
                    //     print *,'m,n,isectional=',m,n,isectional
                    //     print *,'alnsign(m,n,ii)=',alnsign(m,n,ii)
                    //     print *,'rcut,rlo(m,n,ii),rhi(m,n,ii)',Kokkos::exp(xcut),rlo(m,n,ii),rhi(m,n,ii)
                    //     print *,'w,m,na,am=',w,m,na(m,n,ii),am(m,n,ii)
                    //     print *,'integ,sumfs,sigw=',integ,sumfs(m,n,ii),sigw
                    fluxs(m,n,ii) = fluxs(m,n,ii)/fs(m,n,ii);
                }

                fm(m,n,ii)=sumfm(m,n,ii)/(sq2*sqpi*sigw);
                fluxm(m,n,ii)=sumflxm(m,n,ii)/(sq2*sqpi*sigw);

                if(fm(m,n,ii) > 1.01){
                    // if (idiag_fnsm_prob .gt. 0) then
                    //     print *,'fm(',m,n,',ii)=',fm(m,n,ii),' > 1 - activate err43'
                    // endif
                    fluxm(m,n,ii) = fluxm(m,n,ii)/fm(m,n,ii);
                }

            }
        }
        flux_fullact(kk,ii) = sumflx_fullact/(sq2*sqpi*sigw);
        lable6000 = false;
    }   //end lable5000

    if(lable6000){
        wnuc=wbar+wdiab;
        if(wnuc <= 0.){
            for(int n=0; n<ntype_aer; n++){
                for(int m=0; m<nsize_aer(n); m++){
                    fn(m,n,ii)=0;
                    fluxn(m,n,ii)=0;
                    fs(m,n,ii)=0;
                    fluxs(m,n,ii)=0;
                    fm(m,n,ii)=0;
                    fluxm(m,n,ii)=0;
                }
            }
            flux_fullact(kk,ii)=0.;
            return;
        }

        w=wbar;
        alw=alpha*wnuc;
        sqrtalw=Kokkos::sqrt(alw);
        zeta=2.*sqrtalw*aten/(3.*sqrtg);

        if (isectional > 0){
            for(int n = 0;n < ntype_aer; n++){
                if(totn(n,ii) > 1.e-10) eta(0,n,ii)=2*alw*sqrtalw/(totn(n,ii)*beta*sqrtg);
                else eta(0,n,ii)=1.e10;
            }

            if ( smax_prescribed != 1.0) smax = smax_prescribed;
            else {
                for(int max = 0; max < maxd_atype; max++){
                    tmp_nsize_aer(max,ii) = 1;
                }
                maxsat_1(zeta,eta,maxd_atype,ntype_aer,
                      maxd_asize,tmp_nsize_aer,gmsm,gmlnsig,f1,smax,ii);
            }
        }
        else{
            for(int n = 0; n < ntype_aer; n++){
                for(int m = 0; m < nsize_aer(n); m++){
                    if(na(m,n,ii) > 1.e-10) eta(m,n,ii)=2*alw*sqrtalw/(na(m,n,ii)*beta*sqrtg);
                    else eta(m,n,ii)=1.e10;
                }
            }

            if ( smax_prescribed != 1.0) smax = smax_prescribed;
            else{
                maxsat(zeta,eta,maxd_atype,ntype_aer,
                    maxd_asize,nsize_aer,sm,alnsign,f1,smax,ii);
            }
        }

        lnsmax=Kokkos::log(smax);
        xmincoeff=alogaten-2.*third*(lnsmax-alog2)-alog3;

        for(int n = 0; n < ntype_aer; n++){
            for(int m = 0; m < nsize_aer(n); m++){
                if ( bin_is_empty(m,n,ii) ){
                   fn(m,n,ii)=0.;
                   fs(m,n,ii)=0.;
                   fm(m,n,ii)=0.;
                }
               else if ((isectional == 2) || (isectional == 4)){
                    xcut=xmincoeff-third*lnhygro(m,n,ii);
                    lnycut = 3.0 * ( xcut - Kokkos::log(rhi(m,n,ii)) );
                    lnycut = Kokkos::min( lnycut, Kokkos::log(yhi(m,n,ii)*1.0e5) );
                    ycut=Kokkos::exp(lnycut);
                    if(ycut > yhi(m,n,ii)){
                        fn(m,n,ii)=0.;
                        fs(m,n,ii)=0.;
                        fm(m,n,ii)=0.;
                    }
                    else if(ycut < ylo(m,n,ii)){
                        fn(m,n,ii)=1.;
                        fs(m,n,ii)=1.;
                        fm(m,n,ii)=1.;
                    }
                    else if ( bin_is_narrow(m,n,ii) ){
                        if (ycut > ymean(m,n,ii)){
                            fn(m,n,ii)=0.;
                            fs(m,n,ii)=0.;
                            fm(m,n,ii)=0.;
                        }
                        else{
                            fn(m,n,ii)=1.;
                            fs(m,n,ii)=1.;
                            fm(m,n,ii)=1.;
                        }
                    }
                    else{
                        phiyy=ycut/yhi(m,n,ii);
                        fn(m,n,ii) = asub(m,n,ii)*(1.0-phiyy) + 0.5*bsub(m,n,ii)*(1.0-phiyy*phiyy);
                        // if (fn(m,n,ii).lt.zero .or. fn(m,n,ii).gt.one) then
                        //     if (idiag_fnsm_prob .gt. 0) then
                        //         print *,'fn(',m,n,',ii)=',fn(m,n,ii),' outside 0,1 - activate err21'
                        //         print *,'na,volc       =', na(m,n,ii), volc(m,n,ii)
                        //         print *,'asub,bsub     =', asub(m,n,ii), bsub(m,n,ii)
                        //         print *,'yhi,ycut      =', yhi(m,n,ii), ycut
                        //     endif
                        // endif
                        if (fn(m,n,ii) <= zero){
                            fn(m,n,ii)=zero;
                            fs(m,n,ii)=zero;
                            fm(m,n,ii)=zero;
                        }
                        else if (fn(m,n,ii) >= one){
                            fn(m,n,ii)=one;
                            fs(m,n,ii)=one;
                            fm(m,n,ii)=one;
                        }
                        else{
                            fm(m,n,ii) = (yhi(m,n,ii)/ymean(m,n,ii)) * (0.5*asub(m,n,ii)*(1.0-phiyy*phiyy) +
                                  third*bsub(m,n,ii)*(1.0-phiyy*phiyy*phiyy));
                            // if (fm(m,n,ii).lt.fn(m,n,ii) .or. fm(m,n,ii).gt.one) then
                            //     if (idiag_fnsm_prob .gt. 0) then
                            //         print *,'fm(',m,n,',ii)=',fm(m,n,ii),' outside fn,1 - activate err22'
                            //         print *,'na,volc,fn    =', na(m,n,ii), volc(m,n,ii), fn(m,n,ii)
                            //         print *,'asub,bsub     =', asub(m,n,ii), bsub(m,n,ii)
                            //         print *,'yhi,ycut      =', yhi(m,n,ii), ycut
                            //     endif
                            // endif

                            if (fm(m,n,ii) <= fn(m,n,ii)){
                                fm(m,n,ii)=fn(m,n,ii);
                                fs(m,n,ii)=fn(m,n,ii);
                            }
                            else if (fm(m,n,ii) >= one){
                                fm(m,n,ii)=one;
                                fs(m,n,ii)=one;
                                fn(m,n,ii)=one;
                            }
                            else{
                                dumaa = fn(m,n,ii)*(yhi(m,n,ii)/ymean(m,n,ii));
                                fm(m,n,ii) = Kokkos::min( fm(m,n,ii), dumaa );
                                dumaa = 1.0 + (fn(m,n,ii)-1.0)*(ylo(m,n,ii)/ymean(m,n,ii));
                                fm(m,n,ii) = Kokkos::min( fm(m,n,ii), dumaa );

                                betayy = ylo(m,n,ii)/yhi(m,n,ii);
                                dumaa = Kokkos::pow(phiyy,twothird);
                                dumbb = Kokkos::pow(betayy,twothird);
                                fs(m,n,ii) =
                                    (asub(m,n,ii)*(1.0-phiyy*dumaa) +
                                        0.625*bsub(m,n,ii)*(1.0-phiyy*phiyy*dumaa)) /
                                    (asub(m,n,ii)*(1.0-betayy*dumbb) +
                                        0.625*bsub(m,n,ii)*(1.0-betayy*betayy*dumbb));
                                fs(m,n,ii)=Kokkos::max(fs(m,n,ii),fn(m,n,ii));
                                fs(m,n,ii)=Kokkos::min(fs(m,n,ii),fm(m,n,ii));
                            }
                        }
                    }
               }
               else{
                    x=2*(lnsm(m,n,ii)-lnsmax)/(3*sq2*alnsign(m,n,ii));
                    fn(m,n,ii)=0.5*(1.-erf_alt(x));
                    arg=x-sq2*alnsign(m,n,ii);
                    fs(m,n,ii)=0.5*(1.-erf_alt(arg));
                    arg=x-1.5*sq2*alnsign(m,n,ii);
                    fm(m,n,ii)=0.5*(1.-erf_alt(arg));
               }
                if((top && wbar < 0.) || (!top && wbar > 0.)){
                    fluxn(m,n,ii)=fn(m,n,ii)*w;
                    fluxs(m,n,ii)=fs(m,n,ii)*w;
                    fluxm(m,n,ii)=fm(m,n,ii)*w;
                }
                else{
                   fluxn(m,n,ii)=0;
                   fluxs(m,n,ii)=0;
                   fluxm(m,n,ii)=0;
                }
            }//! m=0,nsize_aer(n)
        }//! n=0,ntype_aer

        if((top && wbar < 0.) || (!top && wbar > 0.)) flux_fullact(kk,ii) = w;
        else flux_fullact(kk,ii) = 0.0;
    }
}

KOKKOS_INLINE_FUNCTION real erf_alt(real &x){
    real erf_alt;
    erf_alt = 1. - erfc_num_recipes(x);
    return erf_alt;
}

KOKKOS_INLINE_FUNCTION real erfc_num_recipes(real &x) {
    real  erfc_num_recipes, erfc_dbl, dum, t, zz;

    zz = Kokkos::abs(x);
    t = 1.0/(1.0 + 0.5*zz);

    dum =  ( -zz*zz - 1.26551223 + t*(1.00002368 + t*(0.37409196 +
          t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +
                                           t*(-1.13520398 +
          t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))));

    erfc_dbl = t * Kokkos::exp(dum);

    if (x < 0.0) erfc_dbl = 2.0 - erfc_dbl;

    erfc_num_recipes = erfc_dbl;

    return erfc_num_recipes;
}

template <class Real2View>
KOKKOS_INLINE_FUNCTION void maxsat_init(int &maxd_atype, int &ntype_aer,
           int &maxd_asize, const int1d &nsize_aer,
           const Real2View &alnsign, const Real2View &f1,int i){
    int m,n;
    for(n = 0; n < ntype_aer; n++){
        for(m = 0; m < nsize_aer(n); m++){
            f1(m,n,i)=0.5*Kokkos::exp(2.5*alnsign(m,n,i)*alnsign(m,n,i));
        }
    }
}

template <class Real2View, class Real1View, class Int1View>
KOKKOS_INLINE_FUNCTION void maxsat_init_1(int &maxd_atype, int &ntype_aer,
           int &maxd_asize, const Int1View &nsize_aer,
           const Real1View &gmlnsig, const Real2View &f1,int i){
    int m,n;
    for(n = 0; n < ntype_aer; n++){
        for(m = 0; m < nsize_aer(n,i); m++){
            f1(m,n,i)=0.5*Kokkos::exp(2.5*gmlnsig(n,i)*gmlnsig(n,i));
        }
    }
}

template <class Real2View>
KOKKOS_INLINE_FUNCTION void maxsat(real &zeta, const Real2View &eta,
                        int &maxd_atype,int &ntype_aer,int &maxd_asize,const int1d &nsize_aer,
                        const Real2View &sm, const Real2View &alnsign,
                        const Real2View &f1,real &smax,int i){

    real g1, g2;
    real thesum;
    int m; //! size index
    int n; //! type index

    bool lable1 = true;

    for(int n = 0; n < ntype_aer; n++){
        for(int m = 0; m < nsize_aer(n); m++){
            if(zeta > 1.e5*eta(m,n,i) || sm(m,n,i)*sm(m,n,i) > 1.e5*eta(m,n,i)) smax=1.e-20;
            else {
                lable1 = false;
                break;
            }
            if( ! lable1) break;
        }
    }
    if( lable1 ) return;
    else {
        thesum=0;
        for(int n = 0; n < ntype_aer; n++){
            for(int m = 0; m < nsize_aer(n); m++){
                if(eta(m,n,i) > 1.e-20){
                    g1=Kokkos::sqrt(zeta/eta(m,n,i));
                    g1=g1*g1*g1;
                    g2=sm(m,n,i)/Kokkos::sqrt(eta(m,n,i)+3.*zeta);
                    g2=Kokkos::sqrt(g2);
                    g2=g2*g2*g2;
                    thesum=thesum +
                        (f1(m,n,i)*g1+(1.+0.25*alnsign(m,n,i))*g2)/(sm(m,n,i)*sm(m,n,i));
                }
                else thesum=1.e20;
            }
        }

        smax=1./Kokkos::sqrt(thesum);

        return;
    }
}

template <class Real2View, class Real1View, class Int1View>
KOKKOS_INLINE_FUNCTION void maxsat_1(real &zeta, const Real2View &eta,
                        int &maxd_atype,int &ntype_aer,int &maxd_asize,const Int1View &nsize_aer,
                        const Real1View &gm, const Real1View &gmlnsig,
                        const Real2View &f1,real &smax,int i){

    real g1, g2;
    real thesum;
    int m; //! size index
    int n; //! type index

    bool lable1 = true;

    for(int n = 0; n < ntype_aer; n++){
        for(int m = 0; m < nsize_aer(n,i); m++){
            if(zeta > 1.e5*eta(m,n,i) || gm(n,i)*gm(n,i) > 1.e5*eta(m,n,i)) smax=1.e-20;
            else {
                lable1 = false;
                break;
            }
            if( ! lable1) break;
        }
    }
    if( lable1 ) return;
    else {
        thesum=0;
        for(int n = 0; n < ntype_aer; n++){
            for(int m = 0; m < nsize_aer(n,i); m++){
                if(eta(m,n,i) > 1.e-20){
                    g1=Kokkos::sqrt(zeta/eta(m,n,i));
                    g1=g1*g1*g1;
                    g2=gm(n,i)/Kokkos::sqrt(eta(m,n,i)+3.*zeta);
                    g2=Kokkos::sqrt(g2);
                    g2=g2*g2*g2;
                    thesum=thesum +
                        (f1(m,n,i)*g1+(1.+0.25*gmlnsig(n,i))*g2)/(gm(n,i)*gm(n,i));
                }
                else thesum=1.e20;
            }
        }

        smax=1./Kokkos::sqrt(thesum);

        return;
    }
}

//is_unact = true

template <class QView, class SrcView, class EkkView, class OverlapView, class QactView>
KOKKOS_INLINE_FUNCTION void explmix_t(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,
                       const QactView &qactold,int lnum,int lnumcw,int nnew,int nsav,int i){
    q(0,lnum,nnew,i) = qold(0,lnum,nsav,i) + dt*(-src(0,i) + ekkp(0,i)*(qold(1,lnum,nsav,i)-qold(0,lnum,nsav,i)+qactold(1,i)*(1.0-overlapp(0,i)))
                     + ekkm(0,i)*(qold(0,lnum,nsav,i)-qold(0,lnum,nsav,i)+qactold(0,i)*(1.0-overlapm(0,i))));
    q(0,lnum,nnew,i) = Kokkos::max(q(0,lnum,nnew,i),0.);
    for(int k = 1; k < kte-1; k++){
        q(k,lnum,nnew,i) = qold(k,lnum,nsav,i) + dt*(-src(k,i) + ekkp(k,i)*(qold(k+1,lnum,nsav,i)-qold(k,lnum,nsav,i)+qactold(k+1,i)*(1.0-overlapp(k,i)))
                         + ekkm(k,i)*(qold(k-1,lnum,nsav,i)-qold(k,lnum,nsav,i)+qactold(k-1,i)*(1.0-overlapm(k,i))));
        q(k,lnum,nnew,i) = Kokkos::max(q(k,lnum,nnew,i),0.);
    }
    if(kte > 1){
        q(kte-1,lnum,nnew,i) = qold(kte-1,lnum,nsav,i) + dt*(-src(kte-1,i) + ekkp(kte-1,i)*(qold(kte-1,lnum,nsav,i)-qold(kte-1,lnum,nsav,i)+qactold(kte-1,i)*(1.0-overlapp(kte-1,i)))
                             + ekkm(kte-1,i)*(qold(kte-2,lnum,nsav,i)-qold(kte-1,lnum,nsav,i)+qactold(kte-2,i)*(1.0-overlapm(kte-1,i))));
        q(kte-1,lnum,nnew,i) = Kokkos::max(q(kte-1,lnum,nnew,i),0.);
    }
    q(0,lnum,nnew,i)=q(0,lnum,nnew,i)-surfrate*qold(0,lnum,nsav,i)*dt;
    q(0,lnum,nnew,i)=Kokkos::max(q(0,lnum,nnew,i),0.);
}

template <class QView, class SrcView, class EkkView, class OverlapView>
KOKKOS_INLINE_FUNCTION void explmix_t_3d(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,
                       const QView &qactold,int lnum,int lnumcw,int nnew,int nsav,int i){

    q(0,lnum,nnew,i) = qold(0,lnum,nsav,i) + dt*(-src(0,i) + ekkp(0,i)*(qold(1,lnum,nsav,i)-qold(0,lnum,nsav,i)+qactold(1,lnumcw,nsav,i)*(1.0-overlapp(0,i)))
                     + ekkm(0,i)*(qold(0,lnum,nsav,i)-qold(0,lnum,nsav,i)+qactold(0,lnumcw,nsav,i)*(1.0-overlapm(0,i))));
    q(0,lnum,nnew,i) = Kokkos::max(q(0,lnum,nnew,i),0.);
    for(int k = 1; k < kte-1; k++){
        q(k,lnum,nnew,i) = qold(k,lnum,nsav,i) + dt*(-src(k,i) + ekkp(k,i)*(qold(k+1,lnum,nsav,i)-qold(k,lnum,nsav,i)+qactold(k+1,lnumcw,nsav,i)*(1.0-overlapp(k,i)))
                         + ekkm(k,i)*(qold(k-1,lnum,nsav,i)-qold(k,lnum,nsav,i)+qactold(k-1,lnumcw,nsav,i)*(1.0-overlapm(k,i))));
        q(k,lnum,nnew,i) = Kokkos::max(q(k,lnum,nnew,i),0.);
    }
    if(kte > 1){
        q(kte-1,lnum,nnew,i) = qold(kte-1,lnum,nsav,i) + dt*(-src(kte-1,i) + ekkp(kte-1,i)*(qold(kte-1,lnum,nsav,i)-qold(kte-1,lnum,nsav,i)+qactold(kte-1,lnumcw,nsav,i)*(1.0-overlapp(kte-1,i)))
                             + ekkm(kte-1,i)*(qold(kte-2,lnum,nsav,i)-qold(kte-1,lnum,nsav,i)+qactold(kte-2,lnumcw,nsav,i)*(1.0-overlapm(kte-1,i))));
        q(kte-1,lnum,nnew,i) = Kokkos::max(q(kte-1,lnum,nnew,i),0.);
    }
    q(0,lnum,nnew,i)=q(0,lnum,nnew,i)-surfrate*qold(0,lnum,nsav,i)*dt;
    q(0,lnum,nnew,i)=Kokkos::max(q(0,lnum,nnew,i),0.);
}

//is_unact = flase
template <class QView, class SrcView, class EkkView, class OverlapView>
KOKKOS_INLINE_FUNCTION void explmix_f(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,int i){

    // k=0: km1=0
    q(0,i) = qold(0,i) + dt*(src(0,i) + ekkp(0,i)*(overlapp(0,i)*qold(1,i)-qold(0,i)) +
                         ekkm(0,i)*(overlapm(0,i)*qold(0,i)-qold(0,i)));
    q(0,i) = Kokkos::max(q(0,i),0.);

    for(int k = 1; k < kte-1; k++){
        q(k,i) = qold(k,i) + dt*(src(k,i) + ekkp(k,i)*(overlapp(k,i)*qold(k+1,i)-qold(k,i)) +
                              ekkm(k,i)*(overlapm(k,i)*qold(k-1,i)-qold(k,i)));
        q(k,i) = Kokkos::max(q(k,i),0.);
    }
    // k=kte-1: kp1=kte-1
    if(kte > 1){
        q(kte-1,i) = qold(kte-1,i) + dt*(src(kte-1,i) + ekkp(kte-1,i)*(overlapp(kte-1,i)*qold(kte-1,i)-qold(kte-1,i)) +
                                     ekkm(kte-1,i)*(overlapm(kte-1,i)*qold(kte-2,i)-qold(kte-1,i)));
        q(kte-1,i) = Kokkos::max(q(kte-1,i),0.);
    }
    q(0,i)=q(0,i)-surfrate*qold(0,i)*dt;
    q(0,i)=Kokkos::max(q(0,i),0.);
}

template <class QView, class SrcView, class EkkView, class OverlapView>
KOKKOS_INLINE_FUNCTION void explmix_f_3d(const QView &q, const SrcView &src, const EkkView &ekkp, const EkkView &ekkm, const OverlapView &overlapp, const OverlapView &overlapm,
                       const QView &qold, real surfrate,int kts,int kte,real dt,int lnumcw,int nnew,int nsav,int i){

    q(0,lnumcw,nnew,i) = qold(0,lnumcw,nsav,i) + dt*(src(0,i) + ekkp(0,i)*(overlapp(0,i)*qold(1,lnumcw,nsav,i)-qold(0,lnumcw,nsav,i)) +
                        ekkm(0,i)*(overlapm(0,i)*qold(0,lnumcw,nsav,i)-qold(0,lnumcw,nsav,i)));
    q(0,lnumcw,nnew,i) = Kokkos::max(q(0,lnumcw,nnew,i),0.);
    for(int k = 1; k < kte-1; k++){
        q(k,lnumcw,nnew,i) = qold(k,lnumcw,nsav,i) + dt*(src(k,i) + ekkp(k,i)*(overlapp(k,i)*qold(k+1,lnumcw,nsav,i)-qold(k,lnumcw,nsav,i)) +
                            ekkm(k,i)*(overlapm(k,i)*qold(k-1,lnumcw,nsav,i)-qold(k,lnumcw,nsav,i)));
        q(k,lnumcw,nnew,i) = Kokkos::max(q(k,lnumcw,nnew,i),0.);
    }
    if(kte > 1){
        q(kte-1,lnumcw,nnew,i) = qold(kte-1,lnumcw,nsav,i) + dt*(src(kte-1,i) + ekkp(kte-1,i)*(overlapp(kte-1,i)*qold(kte-1,lnumcw,nsav,i)-qold(kte-1,lnumcw,nsav,i)) +
                                ekkm(kte-1,i)*(overlapm(kte-1,i)*qold(kte-2,lnumcw,nsav,i)-qold(kte-1,lnumcw,nsav,i)));
        q(kte-1,lnumcw,nnew,i) = Kokkos::max(q(kte-1,lnumcw,nnew,i),0.);
    }
    q(0,lnumcw,nnew,i)=q(0,lnumcw,nnew,i)-surfrate*qold(0,lnumcw,nsav,i)*dt;
    q(0,lnumcw,nnew,i)=Kokkos::max(q(0,lnumcw,nnew,i),0.);
}
