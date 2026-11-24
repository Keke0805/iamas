/*
    This computer software was prepared by University of Science and Technology of
    China (USTC)
    NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
    EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

    Module to Compute Aerosol Optical Properties
    * Author: Shen Du and Xin Qu
    Last update: OCT 2024

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
// #include "YAKL.h"
#include <Kokkos_Core.hpp>
#include "flcl-cxx.hpp"
#include <algorithm>
#include <complex>
#include <iostream>
#include <array>
#include <cstring>
#include <stdbool.h>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <chrono>
#include <math.h>
#include "data_rrtmgaeropt.hpp"
#include <fstream>
#include <utility>

// using yakl::Array;
// using Kokkos::SArray;
// using yakl::memHost;
// using yakl::memDevice;
// using yakl::styleFortran;
using std::complex;
using std::string;
// using yakl::fortran::parallel_for;
// using yakl::fortran::Bounds;
// using yakl::fortran::LBnd;
// using yakl::COLON;

struct data_rrtmgaeropt data;

const int ncoef     = data.ncoef;
const int prefr     = data.prefr;
const int prefi     = data.prefi;
const int nswbands  = data.nswbands;
const int nlwbands  = data.nlwbands;
const real pi       = 3.14159265;
static bool ini_fit = true;

static int nrefr, nrefi;
static std::complex<double> crefw;
static real rmin, rmax;
static real xr;
static real pie, third;

#define POW3(x) ((x) * (x) * (x))
#define POW2(x) ((x) * (x))

#define thesum1(x) thesum1[x]
#define thesum2(x) thesum2[x]

typedef Kokkos::View<real**, Kokkos::LayoutLeft> real2d;

typedef Kokkos::View<real***, Kokkos::LayoutLeft> real3d;

typedef Kokkos::View<std::complex<double>***, Kokkos::LayoutLeft,flcl::HostMemorySpace> real_com3d;
typedef Kokkos::View<std::complex<double>****, Kokkos::LayoutLeft,flcl::HostMemorySpace> real_com4d;

constexpr int CHEM_WORKSPACE_CHUNK_SIZE = 64;

struct optical_averaging_workspace {
    real2d dz;
    real3d swsizeaer, swextaer, swwaer, swgaer, swtauaer;
    real3d swabsaer, swbscoef, lwextaer, lwtauaer;
    real3d l2, l3, l4, l5, l6, l7;
    real3d radius_wet_col, number_bin_col, radius_core_col;
    real_com3d swrefindx_col1, swrefindx_core_col1, swrefindx_shell_col1;
    real_com4d swrefindx_col, swrefindx_core_col, swrefindx_shell_col;
    real_com4d lwrefindx_col, lwrefindx_core_col, lwrefindx_shell_col;

    optical_averaging_workspace(int kte, int ncol, int nbin,
                                int nswbands, int nlwbands)
        : dz("optical_ws_dz", kte, ncol),
          swsizeaer("optical_ws_swsizeaer", nswbands, kte, ncol),
          swextaer("optical_ws_swextaer", nswbands, kte, ncol),
          swwaer("optical_ws_swwaer", nswbands, kte, ncol),
          swgaer("optical_ws_swgaer", nswbands, kte, ncol),
          swtauaer("optical_ws_swtauaer", nswbands, kte, ncol),
          swabsaer("optical_ws_swabsaer", nswbands, kte, ncol),
          swbscoef("optical_ws_swbscoef", nswbands, kte, ncol),
          lwextaer("optical_ws_lwextaer", nlwbands, kte, ncol),
          lwtauaer("optical_ws_lwtauaer", nlwbands, kte, ncol),
          l2("optical_ws_l2", nswbands, kte, ncol),
          l3("optical_ws_l3", nswbands, kte, ncol),
          l4("optical_ws_l4", nswbands, kte, ncol),
          l5("optical_ws_l5", nswbands, kte, ncol),
          l6("optical_ws_l6", nswbands, kte, ncol),
          l7("optical_ws_l7", nswbands, kte, ncol),
          radius_wet_col("optical_ws_radius_wet_col", nbin, kte, ncol),
          number_bin_col("optical_ws_number_bin_col", nbin, kte, ncol),
          radius_core_col("optical_ws_radius_core_col", nbin, kte, ncol),
          swrefindx_col1("optical_ws_swrefindx_col1", nbin, kte, ncol),
          swrefindx_core_col1("optical_ws_swrefindx_core_col1", nbin, kte, ncol),
          swrefindx_shell_col1("optical_ws_swrefindx_shell_col1", nbin, kte, ncol),
          swrefindx_col("optical_ws_swrefindx_col", nbin, kte, ncol, nswbands),
          swrefindx_core_col("optical_ws_swrefindx_core_col", nbin, kte, ncol, nswbands),
          swrefindx_shell_col("optical_ws_swrefindx_shell_col", nbin, kte, ncol, nswbands),
          lwrefindx_col("optical_ws_lwrefindx_col", nbin, kte, ncol, nlwbands),
          lwrefindx_core_col("optical_ws_lwrefindx_core_col", nbin, kte, ncol, nlwbands),
          lwrefindx_shell_col("optical_ws_lwrefindx_shell_col", nbin, kte, ncol, nlwbands) {}
};

static optical_averaging_workspace *optical_ws = nullptr;
static int optical_ws_kte = -1;
static int optical_ws_ncol = -1;
static int optical_ws_nbin = -1;
static int optical_ws_nswbands = -1;
static int optical_ws_nlwbands = -1;

static bool optical_workspace_matches(int kte, int ncol, int nbin,
                                      int nswbands, int nlwbands) {
    return optical_ws != nullptr &&
           optical_ws_kte == kte && optical_ws_ncol == ncol &&
           optical_ws_nbin == nbin &&
           optical_ws_nswbands == nswbands &&
           optical_ws_nlwbands == nlwbands;
}

static optical_averaging_workspace &get_optical_workspace(
        int kte, int ncol, int nbin, int nswbands, int nlwbands) {
    if (!optical_workspace_matches(kte, ncol, nbin, nswbands, nlwbands)) {
        if (optical_ws != nullptr) {
            Kokkos::fence();
            delete optical_ws;
            optical_ws = nullptr;
        }
        optical_ws = new optical_averaging_workspace(
            kte, ncol, nbin, nswbands, nlwbands);
        optical_ws_kte = kte;
        optical_ws_ncol = ncol;
        optical_ws_nbin = nbin;
        optical_ws_nswbands = nswbands;
        optical_ws_nlwbands = nlwbands;
    }
    return *optical_ws;
}

extern "C" void optical_averaging_workspace_finalize() {
    if (optical_ws != nullptr) {
        if (Kokkos::is_initialized()) {
            Kokkos::fence();
        }
        delete optical_ws;
        optical_ws = nullptr;
    }
    optical_ws_kte = -1;
    optical_ws_ncol = -1;
    optical_ws_nbin = -1;
    optical_ws_nswbands = -1;
    optical_ws_nlwbands = -1;
}

// clang-format off
KOKKOS_INLINE_FUNCTION void binterp_s(real (*table)[prefi][prefr][ncoef],int ns,int &im,int &jm,real &x,real &y,\
                int &ix, int &jy, real &t, real &u,real out[]);

KOKKOS_INLINE_FUNCTION void binterp_l(real (*table)[prefi][prefr][ncoef],int ns,int &im,int &jm,real &x,real &y,\
                int &ix, int &jy, real &t, real &u,real out[]);

KOKKOS_INLINE_FUNCTION void fitcurv(real *rs,real *yin,real (*coef)[prefi][prefr][ncoef],\
                int nr,int ni,int ns,const int maxm);

KOKKOS_INLINE_FUNCTION void fitcurv_nolog(real *rs,real *yin,real (*coef)[prefi][prefr][ncoef],\
                int nr,int ni,int ns,const int maxm);

KOKKOS_INLINE_FUNCTION void chebft(real (*c)[prefi][prefr][ncoef],\
                int nr,int ni,int ns,int n,real *f);

KOKKOS_INLINE_FUNCTION void miev0(real &xx,complex<double> &crefin,bool &perfct,real &mimcut,bool &anyang,int &numang,real *xmu,int &nmom,\
            int &ipolzn,int &momdim,bool *prnt,real *qext,real *qsca,real *gqsc,int n_,real (*pmom)[8],complex<double> &sforw,  \
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback);

KOKKOS_INLINE_FUNCTION void ckinmi(int &numang,const int &maxang, real &xx, bool &perfct, complex<double> &crefin,int &momdim, int &nmom, int &ipolzn,bool &anyang, \
            real *xmu,bool *calcmo,int &npquan);

KOKKOS_INLINE_FUNCTION void errmsg(string messag, bool fatal);

KOKKOS_INLINE_FUNCTION void wrtbad(string varnam,bool &erflag);

KOKKOS_INLINE_FUNCTION void small1(real &xx,int &numang,real *xmu,real *qext,real *qsca,real *gqsc,int n_,complex<double> &sforw,\
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback,\
            complex<double> *a,complex<double> *b);

KOKKOS_INLINE_FUNCTION void small2(real &xx,complex<double> &cior,bool &calcqe,int &numang,real *xmu,real *qext,real *qsca,real *gqsc,int n_,complex<double> &sforw,\
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback,\
            complex<double> *a,complex<double> *b);

KOKKOS_INLINE_FUNCTION void testmi(real *qext,real *qsca,real *gqsc,int n_,complex<double> &sforw,  \
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback, \
            real (*pmom)[8],int &momdim,bool &ok);

KOKKOS_INLINE_FUNCTION void tstbad(string varnam,real relerr,bool &ok);

KOKKOS_INLINE_FUNCTION void biga(complex<double> &cior,real &xx,int &ntrm,bool &noabs,bool &yesang,real *rbiga,complex<double> *cbiga);

KOKKOS_INLINE_FUNCTION real f1(real mre);

KOKKOS_INLINE_FUNCTION real f2(real mre);

KOKKOS_INLINE_FUNCTION complex<double> confra(int n,complex<double> &zinv,real &xx);

KOKKOS_INLINE_FUNCTION void lpcoef(int &ntrm,int &nmom,int &ipolzn,int &momdim,bool *calcmo,int &npquan,complex<double> *a,complex<double> *b, \
            real (*pmom)[8]);

KOKKOS_INLINE_FUNCTION void lpco1t(int &nmom,int &ipolzn,int &momdim,bool *calcmo,complex<double> *a,complex<double> *b,real (*pmom)[8]);

KOKKOS_INLINE_FUNCTION void lpco2t(int &nmom,int &ipolzn,int &momdim,bool *calcmo,complex<double> *a,complex<double> *b,real (*pmom)[8]);

KOKKOS_INLINE_FUNCTION void miprnt(bool *prnt,real &xx,bool &perfct,complex<double> &crefin,int &numang,real *xmu,real *qext,real *qsca,real *gqsc,int n_,\
            int &nmom,int &ipolzn,int &momdim,bool *calcmo,real (*pmom)[8],complex<double> &sforw,complex<double> &sback,\
            complex<double> *tforw,complex<double> *tback,complex<double> *s1,complex<double> *s2);


KOKKOS_INLINE_FUNCTION void mieaer_init(real &refrmin,real &refrmax,real &refimin,real &refimax,int naerosols,real &refr,real &refi,real &drefr,real &drefi,
		real &pext,real &pscat,real &pabs,real &pasm,real &ppmom2,real &ppmom3,real &ppmom4,real &ppmom5,real &ppmom6,real &ppmom7,
		real &sback2,real &bma,real &bpa,real &rmin,real &rmax,std::complex<double> &sforw,std::complex<double> &sback,
		std::complex<double> *tforw,std::complex<double> *tback,real (*pmom)[8],real *p2,real *p3,real *p4,real *p5,real *p6,real *p7,
		bool &perfct,bool &anyang,bool *prnt,real *xmu,int &numang,std::complex<double> *s1,std::complex<double> *s2,real &mimcut,
		real &thesize,real *qext,real *qsca,real *gqsc,real *qext4,real *qsca4,real *qabs4,real *asymm,real *sb2,std::complex<double> &crefd,
		std::complex<double> &crefw,const int nsiz,int &nmom,int &ipolzn,int &momdim,real *rs,int &nrefr,int &nrefi,real &pie,real &xr);
// KOKKOS_INLINE_FUNCTION void mieaer(const int& id,int iclm,int jclm,const int& nbin_a,
//     const real2d& number_bin_col,const real2d& radius_wet_col,
//     const real_com3d& swrefindx_col,
//     const real_com3d& lwrefindx_col,
//     const reald& dz,
//     const real& curr_secs,const int& kts,const int& kte,
//     const real2d& swsizeaer, const real2d& swextaer,
//     const real2d& swwaer, const real2d& swgaer,
//     const real2d& swtauaer, const real2d& swabsaer,
//     const real2d& lwextaer, const real2d& lwtauaer,
//     const real2d& l2, const real2d& l3,
//     const real2d& l4, const real2d& l5,
//     const real2d& l6, const real2d& l7,
//     const real2d& swbscoef);
KOKKOS_INLINE_FUNCTION void mieaer(const int& id,int iclm,const int& nbin_a,
    const real3d& number_bin_col,const real3d& radius_wet_col,
    const real_com4d& swrefindx_col,
    const real_com4d& lwrefindx_col,
    const real2d& dz,
    const real& curr_secs,const int& kts,const int& kte,
    const real3d& swsizeaer, const real3d& swextaer,
    const real3d& swwaer, const real3d& swgaer,
    const real3d& swtauaer, const real3d& swabsaer,
    const real3d& lwextaer, const real3d& lwtauaer,
    const real3d& l2, const real3d& l3,
    const real3d& l4, const real3d& l5,
    const real3d& l6, const real3d& l7,
    const real3d& swbscoef);
// clang-format on
void call_int() {
    const int ltype = 1;
    const int nsiz  = 200;
    const int nlog  = 30;
    int nmom        = 7;
    int ipolzn      = 0;
    int momdim      = 7;

    // static int nrefr, nrefi;
    int nr, ni;
    std::complex<double> sforw, sback, tforw[2], tback[2];
    real pmom[8][8];
    real p2[nsiz], p3[nsiz], p4[nsiz], p5[nsiz];
    real p6[nsiz], p7[nsiz];
    bool perfct = false, anyang = false, prnt[2];
    prnt[0]     = false;
    prnt[1]     = false;
    real xmu[1] = {1.0};
    int numang  = 0;
    std::complex<double> s1[1], s2[1];
    real mimcut = 0.0;
    real thesize;
    real qext[nsiz];
    real qsca[nsiz];
    real gqsc[nsiz];
    real qext4[nsiz], qsca4[nsiz], qabs4[nsiz];
    real asymm[nsiz];
    real sb2[nsiz];
    std::complex<double> crefin;
    std::complex<double> crefd;
    // static std::complex<double> crefw;
    // static real rmin, rmax;
    real bma, bpa;
    real refr;                             // real part of refractive index
    real refi;                             // imaginary part of refractive index
    real refrmin;                          // minimum of real part of refractive index
    real refrmax;                          // maximum of real part of refractive index
    real refimin;                          // minimum of imag part of refractive index
    real refimax;                          // maximum of imag part of refractive index
    real drefr;                            // increment in real part of refractive index
    real drefi;                            // increment in imag part of refractive index
    std::complex<real> specrefndx(ltype);  // refractive indices
    const int naerosols = 5;

    real weighte, weights, weighta;
    // real x;
    //  real thesum; // for normalizing things
    // real sizem; // size in microns
    int m, j, nc, klevel;
    real pext;    // parameterized specific extinction (cm2/g)
    real pscat;   // scattering cross section
    real pabs;    // parameterized specific extinction (cm2/g)
    real pasm;    // parameterized asymmetry factor
    real ppmom2;  // 2 Lengendre expansion coefficient (numbered 0,1,2,...)
    real ppmom3;  // 3     ...
    real ppmom4;  // 4     ...
    real ppmom5;  // 5     ...
    real ppmom6;  // 6     ...
    real ppmom7;  // 7     ...
    real sback2;  // JCB 2007/02/01 sback*conjg(sback)
    // static real xrmin, xrmax, xr;
    real rs[nsiz];   // surface mode radius (cm)
    real xrad;       // normalized aerosol radius
    real ch[ncoef];  // chebychev polynomial

    int i, k, l, ns;
    // static real pie, third;
    int ibin;
    std::string msg;
    int kcallmieaer = 0, kcallmieaer2 = 0;
    pie   = 4.0 * std::atan(1.0);
    third = 1.0 / 3.0;
    rmin  = data.rmmin;
    rmax  = data.rmmax;

    mieaer_init(refrmin, refrmax, refimin, refimax, naerosols, refr, refi, drefr, drefi, pext, pscat, pabs, pasm, ppmom2, ppmom3, ppmom4, ppmom5, ppmom6, ppmom7,
                sback2, bma, bpa, rmin, rmax, sforw, sback, tforw, tback, pmom, p2, p3, p4, p5, p6, p7, perfct, anyang, prnt, xmu, numang, s1, s2, mimcut, thesize,
                qext, qsca, gqsc, qext4, qsca4, qabs4, asymm, sb2, crefd, crefw, nsiz, nmom, ipolzn, momdim, rs, nrefr, nrefi, pie, xr);
}

extern "C" void mieaer_kokkos(int &kts, int &kte, int &nbin_o, int &nswbands, int &nlwbands,
                              int &id, flcl_ndarray_t *dz_p, real &curr_secs, flcl_ndarray_t *swsizeaer_p, flcl_ndarray_t *swextaer_p, flcl_ndarray_t *swwaer_p, flcl_ndarray_t *swgaer_p,
                              flcl_ndarray_t *tauaerlw_p, flcl_ndarray_t *lwtauaer_p,
                              flcl_ndarray_t *l2aer_p, flcl_ndarray_t *l2_p, flcl_ndarray_t *l3aer_p, flcl_ndarray_t *l3_p, flcl_ndarray_t *l4aer_p, flcl_ndarray_t *l4_p,
                              flcl_ndarray_t *l5aer_p, flcl_ndarray_t *l5_p, flcl_ndarray_t *l6aer_p, flcl_ndarray_t *l6_p, flcl_ndarray_t *l7aer_p, flcl_ndarray_t *l7_p,
                              flcl_ndarray_t *number_bin_col_p, flcl_ndarray_t *number_bin_p, flcl_ndarray_t *radius_wet_col_p, flcl_ndarray_t *radius_wet_p,
                              flcl_ndarray_t *swrefindx_col_p, flcl_ndarray_t *swrefindx_col1_p,
                              flcl_ndarray_t *swrefindx_p, flcl_ndarray_t *lwrefindx_col_p, flcl_ndarray_t *lwrefindx_p,
                              flcl_ndarray_t *radius_core_col_p, flcl_ndarray_t *radius_core_p,
                              flcl_ndarray_t *swrefindx_core_col_p, flcl_ndarray_t *swrefindx_core_p, flcl_ndarray_t *swrefindx_core_col1_p,
                              flcl_ndarray_t *swrefindx_shell_col_p, flcl_ndarray_t *swrefindx_shell_p,
                              flcl_ndarray_t *swrefindx_shell_col1_p, flcl_ndarray_t *lwrefindx_core_col_p, flcl_ndarray_t *lwrefindx_core_p,
                              flcl_ndarray_t *lwrefindx_shell_p,
                              int &option_method, int &option_mie, flcl_ndarray_t *tauaersw_p, flcl_ndarray_t *swtauaer_p,
                              flcl_ndarray_t *absaersw_p, flcl_ndarray_t *swabsaer_p, flcl_ndarray_t *extaersw_p, flcl_ndarray_t *gaersw_p,
                              flcl_ndarray_t *waersw_p, flcl_ndarray_t *bscoefsw_p, flcl_ndarray_t *swbscoef_p,
                              flcl_ndarray_t *extaerlw_p, flcl_ndarray_t *lwextaer_p, int &dust_diag,
                              flcl_ndarray_t *tauaersw_dust_p, flcl_ndarray_t *absaersw_dust_p, flcl_ndarray_t *gaersw_dust_p, flcl_ndarray_t *waersw_dust_p,
                              flcl_ndarray_t *tauaerlw_dust_p, int &seas_diag, flcl_ndarray_t *tauaersw_seas_p, flcl_ndarray_t *absaersw_seas_p,
                              flcl_ndarray_t *gaersw_seas_p, flcl_ndarray_t *waersw_seas_p, flcl_ndarray_t *tauaerlw_seas_p,
                              int &its, int &ite,
                              flcl_ndarray_t *dz8w_p, flcl_ndarray_t *lwrefindx_shell_col_p, int &iloop) {
    // (dz)\(iclm, ([^,()]+)\)
    // (swsizeaer|swextaer|swwaer|swgaer|swtauaer|swabsaer|swbscoef|lwextaer|lwtauaer|l2|l3|l4|l5|l6|l7|radius_wet_col|number_bin_col|radius_core_col|swrefindx_col1|swrefindx_core_col1|swrefindx_shell_col1)\(iclm, ([^,()]+), ([^(),]+)\)
    // (swrefindx_col|swrefindx_core_col|swrefindx_shell_col|lwrefindx_col|lwrefindx_core_col|lwrefindx_shell_col)\(iclm, ([^,()]+), ([^(),]+), ([^(),]+)\)

    const int total_ncol = ite - its + 1;
    if (total_ncol <= 0) return;
#ifdef USE_SERIALIZATION
    const int workspace_ncol = 1;
#else
    const int workspace_ncol = CHEM_WORKSPACE_CHUNK_SIZE;
#endif

    auto &workspace = get_optical_workspace(
        kte, workspace_ncol, nbin_o, nswbands, nlwbands);
    auto dz                      = workspace.dz;
    auto swsizeaer               = workspace.swsizeaer;
    auto swextaer                = workspace.swextaer;
    auto swwaer                  = workspace.swwaer;
    auto swgaer                  = workspace.swgaer;
    auto swtauaer                = workspace.swtauaer;
    auto swabsaer                = workspace.swabsaer;
    auto swbscoef                = workspace.swbscoef;
    auto lwextaer                = workspace.lwextaer;
    auto lwtauaer                = workspace.lwtauaer;
    auto l2                      = workspace.l2;
    auto l3                      = workspace.l3;
    auto l4                      = workspace.l4;
    auto l5                      = workspace.l5;
    auto l6                      = workspace.l6;
    auto l7                      = workspace.l7;
    auto radius_wet_col          = workspace.radius_wet_col;
    auto number_bin_col          = workspace.number_bin_col;
    auto radius_core_col         = workspace.radius_core_col;
    auto swrefindx_col1          = workspace.swrefindx_col1;
    auto swrefindx_core_col1     = workspace.swrefindx_core_col1;
    auto swrefindx_shell_col1    = workspace.swrefindx_shell_col1;
    auto swrefindx_col           = workspace.swrefindx_col;
    auto swrefindx_core_col      = workspace.swrefindx_core_col;
    auto swrefindx_shell_col     = workspace.swrefindx_shell_col;
    auto lwrefindx_col           = workspace.lwrefindx_col;
    auto lwrefindx_core_col      = workspace.lwrefindx_core_col;
    auto lwrefindx_shell_col     = workspace.lwrefindx_shell_col;


    auto dz8w_all = flcl::view_from_ndarray<real**>(*dz8w_p);

    auto tauaersw_all = flcl::view_from_ndarray<real***>(*tauaersw_p);
    auto absaersw_all = flcl::view_from_ndarray<real***>(*absaersw_p);
    auto extaersw_all = flcl::view_from_ndarray<real***>(*extaersw_p);
    auto gaersw_all = flcl::view_from_ndarray<real***>(*gaersw_p);
    auto waersw_all = flcl::view_from_ndarray<real***>(*waersw_p);
    auto bscoefsw_all = flcl::view_from_ndarray<real***>(*bscoefsw_p);
    auto tauaerlw_all = flcl::view_from_ndarray<real***>(*tauaerlw_p);
    auto extaerlw_all = flcl::view_from_ndarray<real***>(*extaerlw_p);
    auto l2aer_all = flcl::view_from_ndarray<real***>(*l2aer_p);
    auto l3aer_all = flcl::view_from_ndarray<real***>(*l3aer_p);
    auto l4aer_all = flcl::view_from_ndarray<real***>(*l4aer_p);
    auto l5aer_all = flcl::view_from_ndarray<real***>(*l5aer_p);
    auto l6aer_all = flcl::view_from_ndarray<real***>(*l6aer_p);
    auto l7aer_all = flcl::view_from_ndarray<real***>(*l7aer_p);
    auto radius_wet_all = flcl::view_from_ndarray<real***>(*radius_wet_p);
    auto number_bin_all = flcl::view_from_ndarray<real***>(*number_bin_p);
    auto radius_core_all = flcl::view_from_ndarray<real***>(*radius_core_p);
    auto tauaersw_seas_all = flcl::view_from_ndarray<real***>(*tauaersw_seas_p);
    auto absaersw_seas_all = flcl::view_from_ndarray<real***>(*absaersw_seas_p);
    auto gaersw_seas_all = flcl::view_from_ndarray<real***>(*gaersw_seas_p);
    auto waersw_seas_all = flcl::view_from_ndarray<real***>(*waersw_seas_p);
    auto tauaersw_dust_all = flcl::view_from_ndarray<real***>(*tauaersw_dust_p);
    auto absaersw_dust_all = flcl::view_from_ndarray<real***>(*absaersw_dust_p);
    auto gaersw_dust_all = flcl::view_from_ndarray<real***>(*gaersw_dust_p);
    auto waersw_dust_all = flcl::view_from_ndarray<real***>(*waersw_dust_p);
    auto tauaerlw_dust_all = flcl::view_from_ndarray<real***>(*tauaerlw_dust_p);
    auto tauaerlw_seas_all = flcl::view_from_ndarray<real***>(*tauaerlw_seas_p);

    auto swrefindx_all = flcl::view_from_ndarray<std::complex<double>****>(*swrefindx_p);
    auto swrefindx_core_all = flcl::view_from_ndarray<std::complex<double>****>(*swrefindx_core_p);
    auto swrefindx_shell_all = flcl::view_from_ndarray<std::complex<double>****>(*swrefindx_shell_p);
    auto lwrefindx_all = flcl::view_from_ndarray<std::complex<double>****>(*lwrefindx_p);
    auto lwrefindx_core_all = flcl::view_from_ndarray<std::complex<double>****>(*lwrefindx_core_p);
    auto lwrefindx_shell_all = flcl::view_from_ndarray<std::complex<double>****>(*lwrefindx_shell_p);


    if (ini_fit) {
        // printf("===========initial===========\n");
        call_int();
        ini_fit = false;
    }
    // Serial execution reuses one scratch column. Parallel execution reuses a
    // fixed 64-column workspace, including for the tail chunk.
#ifdef USE_SERIALIZATION
    for (int column = 0; column < total_ncol; ++column) {
    const int current_chunk_len = 1;
    const auto column_range = std::make_pair(column, column + 1);
#else
    for (int chunk_start = 0; chunk_start < total_ncol;
         chunk_start += CHEM_WORKSPACE_CHUNK_SIZE) {
    const int current_chunk_len = std::min(
        CHEM_WORKSPACE_CHUNK_SIZE, total_ncol - chunk_start);
    const auto column_range = std::make_pair(
        chunk_start, chunk_start + current_chunk_len);
#endif

    auto dz8w = Kokkos::subview(dz8w_all, Kokkos::ALL(), column_range);

    auto tauaersw = Kokkos::subview(
        tauaersw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto absaersw = Kokkos::subview(
        absaersw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto extaersw = Kokkos::subview(
        extaersw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto gaersw = Kokkos::subview(
        gaersw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto waersw = Kokkos::subview(
        waersw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto bscoefsw = Kokkos::subview(
        bscoefsw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto tauaerlw = Kokkos::subview(
        tauaerlw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto extaerlw = Kokkos::subview(
        extaerlw_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto l2aer = Kokkos::subview(
        l2aer_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto l3aer = Kokkos::subview(
        l3aer_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto l4aer = Kokkos::subview(
        l4aer_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto l5aer = Kokkos::subview(
        l5aer_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto l6aer = Kokkos::subview(
        l6aer_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto l7aer = Kokkos::subview(
        l7aer_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto radius_wet = Kokkos::subview(
        radius_wet_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto number_bin = Kokkos::subview(
        number_bin_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto radius_core = Kokkos::subview(
        radius_core_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto tauaersw_seas = Kokkos::subview(
        tauaersw_seas_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto absaersw_seas = Kokkos::subview(
        absaersw_seas_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto gaersw_seas = Kokkos::subview(
        gaersw_seas_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto waersw_seas = Kokkos::subview(
        waersw_seas_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto tauaersw_dust = Kokkos::subview(
        tauaersw_dust_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto absaersw_dust = Kokkos::subview(
        absaersw_dust_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto gaersw_dust = Kokkos::subview(
        gaersw_dust_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto waersw_dust = Kokkos::subview(
        waersw_dust_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto tauaerlw_dust = Kokkos::subview(
        tauaerlw_dust_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
    auto tauaerlw_seas = Kokkos::subview(
        tauaerlw_seas_all, Kokkos::ALL(), Kokkos::ALL(), column_range);

    auto swrefindx = Kokkos::subview(
        swrefindx_all, Kokkos::ALL(), Kokkos::ALL(),
        column_range, Kokkos::ALL());
    auto swrefindx_core = Kokkos::subview(
        swrefindx_core_all, Kokkos::ALL(), Kokkos::ALL(),
        column_range, Kokkos::ALL());
    auto swrefindx_shell = Kokkos::subview(
        swrefindx_shell_all, Kokkos::ALL(), Kokkos::ALL(),
        column_range, Kokkos::ALL());
    auto lwrefindx = Kokkos::subview(
        lwrefindx_all, Kokkos::ALL(), Kokkos::ALL(),
        column_range, Kokkos::ALL());
    auto lwrefindx_core = Kokkos::subview(
        lwrefindx_core_all, Kokkos::ALL(), Kokkos::ALL(),
        column_range, Kokkos::ALL());
    auto lwrefindx_shell = Kokkos::subview(
        lwrefindx_shell_all, Kokkos::ALL(), Kokkos::ALL(),
        column_range, Kokkos::ALL());

    Kokkos::parallel_for(
        "mieaer_kokkos", current_chunk_len, KOKKOS_LAMBDA(int idx) {
            int iclm = idx;
            real fv;
            std::complex<real> aa;
            std::complex<real> bb;
            for (int k = 0; k < kte; k++) {
                dz(k,iclm) = dz8w(k,iclm);
            }

            for (int k = 0; k < kte; k++) {
                for (int isize = 0; isize < nbin_o; isize++) {
                    number_bin_col(isize, k, iclm)       = number_bin(isize,k,iclm);
                    radius_wet_col(isize, k, iclm)       = radius_wet(isize,k,iclm);
                    swrefindx_col1(isize, k, iclm)       = swrefindx(isize,k,iclm, 2);
                    radius_core_col(isize, k, iclm)      = radius_core(isize,k,iclm);
                    swrefindx_core_col1(isize, k, iclm)  = swrefindx_core(isize,k,iclm, 2);
                    swrefindx_shell_col1(isize, k, iclm) = swrefindx_shell(isize,k,iclm, 2);

                    for (int i = 0; i < nswbands; i++) {
                        swrefindx_col(isize, k, iclm,i)       = swrefindx(isize,k,iclm, i);
                        swrefindx_core_col(isize, k, iclm,i)  = swrefindx_core(isize,k,iclm, i);
                        swrefindx_shell_col(isize, k, iclm,i) = swrefindx_shell(isize,k,iclm, i);
                    }

                    for (int i = 0; i < nlwbands; i++) {
                        lwrefindx_col(isize, k, iclm,i)       = lwrefindx(isize,k,iclm, i);
                        lwrefindx_core_col(isize, k, iclm,i)  = lwrefindx_core(isize,k,iclm, i);
                        lwrefindx_shell_col(isize, k, iclm,i) = lwrefindx_shell(isize,k,iclm, i);
                    }

                    if (option_method == 3 && option_mie == 2) {
                        swrefindx_col1(isize, k, iclm) = swrefindx_shell(isize,k,iclm, 2);
                        for (int i = 0; i < nswbands; i++) {
                            swrefindx_col(isize, k, iclm,i) = swrefindx_shell(isize,k,iclm, i);
                        }
                    }

                    if (radius_wet_col(isize, k, iclm) < 1e-20) {
                        radius_core_col(isize, k, iclm) = 0.0;
                    } else if (POW3(radius_core_col(isize, k, iclm)) != 0) {
                        if (radius_core_col(isize, k, iclm) / POW3(radius_wet_col(isize, k, iclm)) <= 0.0001) {
                            radius_core_col(isize, k, iclm) = 0.0;
                        }
                    }
                }
            }
            if (option_method == 2) {
                for (int k = 0; k < kte; k++) {
                    for (int isize = 0; isize < nbin_o; isize++) {
                        for (int ns = 0; ns < nswbands; ns++) {
                            fv                                      = POW3(radius_core_col(isize, k, iclm) / radius_wet_col(isize, k, iclm));
                            aa                                      = (POW2(swrefindx_core_col( isize, k, iclm,ns)) + real(2.0) * POW2(swrefindx_shell(isize,k,iclm, ns)));
                            bb                                      = fv * (POW2(swrefindx_core_col( isize, k, iclm,ns)) - POW2(swrefindx_shell(isize,k,iclm, ns)));
                            swrefindx_col(isize, k, iclm,ns) = swrefindx_shell(isize,k,iclm, ns) * sqrt((aa + real(2.0) * bb) / (aa - bb));

                            if (ns == 2) {
                                swrefindx_col1(isize, k, iclm) = swrefindx_shell(isize,k,iclm, ns) * sqrt((aa + real(2.0) * bb) / (aa - bb));
                            }
                        }
                    }
                }
            }

            if (option_method <= 2) {
                for (int k = 0; k < kte; k++) {
                    for (int isize = 0; isize < nbin_o; isize++) {
                        radius_core_col(isize, k, iclm) = 0.0;
                        swrefindx_core_col1(isize, k, iclm).real(0.0);
                        swrefindx_core_col1(isize, k, iclm).imag(0.0);
                        for (int i = 0; i < nswbands; i++) {
                            swrefindx_core_col(isize, k, iclm,i).real(0.0);
                            swrefindx_core_col(isize, k, iclm,i).imag(0.0);
                        }
                    }
                }
            }

            if (option_mie == 1) {
                mieaer(id, iclm, nbin_o,
                       number_bin_col, radius_wet_col, swrefindx_col, lwrefindx_col,
                       dz, curr_secs, kts, kte,
                       swsizeaer, swextaer, swwaer, swgaer, swtauaer, swabsaer, lwextaer, lwtauaer,
                       l2, l3, l4, l5, l6, l7, swbscoef);
            }

    //         // if(option_mie >= 2 && option_method <= 2){
    //         //  slave___module_optical_averaging_MOD_mieaer_sc(debug_on,id,iclm,jclm,nbin_o,
    //         //                                          number_bin_col_p,radius_wet_col_p,swrefindx_col1_p,radius_core_col_p,
    //         //                                          swrefindx_core_col1_p,dz_p, curr_secs,kte,swsizeaer_p,swextaer_p,
    //         //                                          swwaer_p,swgaer_p,swtauaer_p,l2_p,l3_p,l4_p,l5_p,l6_p,l7_p,swbscoef_p);
    //         // }

    //         // if(option_mie >= 2 && option_method == 3){
    //         //  slave___module_optical_averaging_MOD_mieaer_sc(debug_on,id,iclm,jclm,nbin_o,
    //         //                                          number_bin_col_p,radius_wet_col_p,swrefindx_shell_col1_p,radius_core_col_p,
    //         //                                          swrefindx_core_col1_p,dz_p, curr_secs,kte,swsizeaer_p,swextaer_p,
    //         //                                          swwaer_p,swgaer_p,swtauaer_p,l2_p,l3_p,l4_p,l5_p,l6_p,l7_p,swbscoef_p);
    //         // }

            if (iloop == 1) {
                for (int k = 0; k < kte; k++) {
                    for (int i = 0; i < nswbands; i++) {
                        tauaersw(i, k, iclm) = swtauaer(i, k, iclm);
                        absaersw(i, k, iclm) = swabsaer(i, k, iclm);
                        extaersw(i, k, iclm) = swextaer(i, k, iclm);
                        gaersw(i, k, iclm)   = swgaer(i, k, iclm);
                        waersw(i, k, iclm)   = swwaer(i, k, iclm);
                        bscoefsw(i, k, iclm) = swbscoef(i, k, iclm);
                        l2aer(i, k, iclm)    = l2(i, k, iclm);
                        l3aer(i, k, iclm)    = l3(i, k, iclm);
                        l4aer(i, k, iclm)    = l4(i, k, iclm);
                        l5aer(i, k, iclm)    = l5(i, k, iclm);
                        l6aer(i, k, iclm)    = l6(i, k, iclm);
                        l7aer(i, k, iclm)    = l7(i, k, iclm);
                    }

                    for (int i = 0; i < nlwbands; i++) {
                        tauaerlw(i, k, iclm) = lwtauaer(i, k, iclm);
                        extaerlw(i, k, iclm) = lwextaer(i, k, iclm);
                    }
                }
            }

            if (dust_diag == 1) {
                for (int k = 0; k < kte; k++) {
                    for (int i = 0; i < nswbands; i++) {
                        tauaersw_dust(i, k, iclm) = swtauaer(i, k, iclm);
                        absaersw_dust(i, k, iclm) = swabsaer(i, k, iclm);
                        gaersw_dust(i, k, iclm)   = swgaer(i, k, iclm);
                        waersw_dust(i, k, iclm)   = swwaer(i, k, iclm);
                    }

                    for (int i = 0; i < nlwbands; i++) {
                        tauaerlw_dust(i, k, iclm) = lwtauaer(i, k, iclm);
                    }
                }
            }

            if (seas_diag == 1) {
                for (int k = 0; k < kte; k++) {
                    for (int i = 0; i < nswbands; i++) {
                        tauaersw_seas(i, k, iclm) = swtauaer(i, k, iclm);
                        absaersw_seas(i, k, iclm) = swabsaer(i, k, iclm);
                        gaersw_seas(i, k, iclm)   = swgaer(i, k, iclm);
                        waersw_seas(i, k, iclm)   = swwaer(i, k, iclm);
                    }

                    for (int i = 0; i < nlwbands; i++) {
                        tauaerlw_seas(i, k, iclm) = lwtauaer(i, k, iclm);
                    }
                }
            }
        });
    }
    Kokkos::fence();
}

// clang-format off
KOKKOS_INLINE_FUNCTION void mieaer_init(real &refrmin,real &refrmax,real &refimin,real &refimax,const int naerosols,real &refr,real &refi,real &drefr,real &drefi,
		real &pext,real &pscat,real &pabs,real &pasm,real &ppmom2,real &ppmom3,real &ppmom4,real &ppmom5,real &ppmom6,real &ppmom7,
		real &sback2,real &bma,real &bpa,real &rmin,real &rmax,std::complex<double> &sforw,std::complex<double> &sback,
		std::complex<double> *tforw,std::complex<double> *tback,real (*pmom)[8],real *p2,real *p3,real *p4,real *p5,real *p6,real *p7,
		bool &perfct,bool &anyang,bool *prnt,real *xmu,int &numang,std::complex<double> *s1,std::complex<double> *s2,real &mimcut,
		real &thesize,real *qext,real *qsca,real *gqsc,real *qext4,real *qsca4,real *qabs4,real *asymm,real *sb2,std::complex<double> &crefd,
		std::complex<double> &crefw,const int nsiz,int &nmom,int &ipolzn,int &momdim,real *rs,int &nrefr,int &nrefi,real &pie,real &xr){
	//----------------------------------------------------------------------
        //shortwave
        //---------------------------------------------------------------------
        //wavelength loop
	int ns,l,k,i;
        for(ns=0;ns<data.nswbands;ns++){
            data.wavmidsw[ns] = 0.5*(data.wavmin[ns] + data.wavmax[ns])*1.0e-4;

            //parameterize aerosol radiative properties in terms of
            //relative humidity, surface mode wet radius, aerosol species,
            //and wavelength
            //first find min,max of real and imaginary parts of refractive index
            //real and imaginary parts of water refractive index

            data.crefwsw[ns]=std::complex<real>(data.refrwsw[ns],data.refiwsw[ns]);
            refrmin=std::real(data.crefwsw[ns]);
            refrmax=std::real(data.crefwsw[ns]);

            //change imaginary part of the refractive index from positive to negative
            refimin=-std::imag(data.crefwsw[ns]);
            refimax=-std::imag(data.crefwsw[ns]);

            for (int l = 0; l < naerosols; l++) {
            switch (l) {
                case 0:
                     refr = data.refrsw_dust[ns];
                     refi = -data.refisw_dust[ns];
                    break;
                case 1:
                     refr = data.refrsw_bc[ns];
                     refi = -data.refisw_bc[ns];
                    break;
                case 2:
                    refr = data.refrsw_oc[ns];
                    refi = -data.refisw_oc[ns];
                    break;
                case 3:
                    refr = data.refrsw_seas[ns];
                    refi = -data.refisw_seas[ns];
                    break;
                case 4:
                    refr = data.refrsw_sulf[ns];
                    refi = -data.refisw_sulf[ns];
                    break;
                }

                refrmin = std::min(refrmin, refr);
                refrmax = std::max(refrmax, refr);
                refimin = std::min(refimin, refi);
                refimax = std::max(refimax, refi);
            }

            drefr = refrmax - refrmin;
            if (drefr > 1e-4) {
                nrefr = data.prefr;
                drefr = drefr / (nrefr - 1);
            } else {
                nrefr = 1;
            }

            drefi = refimax - refimin;
            if (drefi > 1e-4) {
                nrefi = data.prefi;
                drefi = drefi / (nrefi - 1);
            } else {
                nrefi = 1;
            }

            bma = 0.5 * std::log(rmax / rmin); // JCB
            bpa = 0.5 * std::log(rmax * rmin); // JCB

            for (int nr = 0; nr < nrefr; nr++) {
                for(int ni = 0;ni < nrefi;ni++){
                    data.refrtabsw[ns][nr] = refrmin + nr  * drefr;
                    data.refitabsw[ns][ni] = refimin / 0.2 * pow(0.2, static_cast<real>(ni+1)); // slightly different from Ghan and Zaveri
                    if (ni == nrefi-1) {
                        data.refitabsw[ns][ni] = -1.0e-20; // JCB change
                    }
                    crefd=std::complex<real>(data.refrtabsw[ns][nr], data.refitabsw[ns][ni]);

                    // Mie calculations of optical efficiencies
                    for (int n = 0; n < nsiz; n++) {
                        xr = cos(pie * (static_cast<real>(n) + 0.5) / static_cast<real>(nsiz));
                        rs[n] = exp(xr * bma + bpa);

                        // Size parameter and weighted refractive index
                        thesize = 2. * pie * rs[n] / data.wavmidsw[ns];
                        thesize = std::min(static_cast<double>(thesize), 10000.0);

                        miev0(thesize, crefd, perfct, mimcut, anyang, numang, xmu, nmom, ipolzn, momdim, prnt, qext, qsca, gqsc,n, pmom, sforw, sback, s1, s2, tforw, tback);
                        qext4[n] = std::max(static_cast<double>(qext[n]), 1.e-20);
                        qsca4[n] = std::min(static_cast<double>(qsca[n]), static_cast<double>(qext[n]));
                        qabs4[n] = qext4[n] - qsca4[n];
                        qabs4[n] = std::max(static_cast<double>(qabs4[n]), 1.e-20); // avoid 0
                        asymm[n] = gqsc[n] / qsca4[n]; // assume always greater than zero
                        asymm[n] = std::max(static_cast<double>(asymm[n]), 1.e-20); // avoid negative

                        // Coefficients of phase function expansion; note modification by JCB of miev0 coefficients
                        p2[n] = pmom[0][2] / pmom[0][0] * 5.0;
                        p3[n] = pmom[0][3] / pmom[0][0] * 7.0;
                        p4[n] = pmom[0][4] / pmom[0][0] * 9.0;
                        p5[n] = pmom[0][5] / pmom[0][0] * 11.0;
                        p6[n] = pmom[0][6] / pmom[0][0] * 13.0;
                        p7[n] = pmom[0][7] / pmom[0][0] * 15.0;

                        // Backscattering efficiency, Bohren and Huffman, page 122
                        // As stated by Bohren and Huffman, this is 4*pie times what is should be
                        // May need to be smoothed - a very rough function - for the time being we won't apply smoothing
                        // and let the integration over the size distribution be the smoothing
                        sb2[n] = 4.0 * std::real(sback * std::conj(sback)) / (thesize * thesize); // JCB 2007/02/01

                        // PMA makes sure it is greater than zero
                        qext4[n] = std::max(static_cast<double>(qext4[n]), 1.e-20);
                        qabs4[n] = std::max(static_cast<double>(qabs4[n]), 1.e-20);
                        qsca4[n] = std::max(static_cast<double>(qsca4[n]), 1.e-20);
                        asymm[n] = std::max(static_cast<double>(asymm[n]), 1.e-20);
                        sb2[n] = std::max(static_cast<double>(sb2[n]), 1.e-20);
                    }
                    fitcurv(rs,qext4,data.extpsw,nr,ni,ns,nsiz);
                    fitcurv(rs,qabs4,data.abspsw,nr,ni,ns,nsiz);
                    fitcurv(rs,qsca4,data.ascatpsw,nr,ni,ns,nsiz); // scattering efficiency
                    fitcurv(rs,asymm,data.asmpsw,nr,ni,ns,nsiz);
                    fitcurv(rs,sb2,data.sbackpsw,nr,ni,ns,nsiz); // backscattering efficiency
                    fitcurv_nolog(rs,p2,data.pmom2psw,nr,ni,ns,nsiz);
                    fitcurv_nolog(rs,p3,data.pmom3psw,nr,ni,ns,nsiz);
                    fitcurv_nolog(rs,p4,data.pmom4psw,nr,ni,ns,nsiz);
                    fitcurv_nolog(rs,p5,data.pmom5psw,nr,ni,ns,nsiz);
                    fitcurv_nolog(rs,p6,data.pmom6psw,nr,ni,ns,nsiz);
                    fitcurv_nolog(rs,p7,data.pmom7psw,nr,ni,ns,nsiz);
                }
            }
        }


        //----------------------------------------------------------------------
        //longwave
        //---------------------------------------------------------------------
        //wavelength loop
        for (int ns = 0; ns < data.nlwbands; ns++) {
            // wavelength for longwave 1/cm --> cm
            data.wavmidlw[ns] = 0.5 * (1.0 / data.wavenumber1_longwave[ns] + 1.0 / data.wavenumber2_longwave[ns]);

            data.crefwlw[ns] = std::complex<real>(data.refrwlw[ns], data.refiwlw[ns]);
            refrmin = std::real(data.crefwlw[ns]);
            refrmax = std::real(data.crefwlw[ns]);
            refimin = -std::imag(data.crefwlw[ns]);
            refimax = -std::imag(data.crefwlw[ns]);

            // aerosol species loop (dust, BC, OC, Sea Salt, and Sulfate)
            for (int l = 0; l < naerosols; l++) {
                switch (l) {
                    case 0:
                        refr = data.refrlw_dust[ns];
                        refi = -data.refilw_dust[ns];
                        break;
                    case 1:
                        refr = data.refrlw_bc[ns];
                        refi = -data.refilw_bc[ns];
                        break;
                    case 2:
                        refr = data.refrlw_oc[ns];
                        refi = -data.refilw_oc[ns];
                        break;
                    case 3:
                        refr = data.refrlw_seas[ns];
                        refi = -data.refilw_seas[ns];
                        break;
                    case 4:
                        refr = data.refrlw_sulf[ns];
                        refi = -data.refilw_sulf[ns];
                        break;
                }
                refrmin = std::min(refrmin, refr);
                refrmax = std::max(refrmax, refr);
                refimin = std::min(refimin, refi);
                refimax = std::max(refimax, refi);
            }


            drefr = refrmax - refrmin;
            if (drefr > 1e-4) {
                nrefr = data.prefr;
                drefr = drefr / (nrefr - 1);
            }else{
                nrefr = 1;
            }

            drefi = refimax - refimin;
            if (drefi > 1e-4) {
                nrefi = data.prefi;
                drefi = drefi / (nrefi - 1);
            }else{
                nrefi = 1;
            }

            bma = 0.5 * std::log(rmax / rmin);
            bpa = 0.5 * std::log(rmax * rmin);

            for (int nr = 0; nr < nrefr; nr++){
                for(int ni = 0;ni < nrefi;ni++){
                    data.refrtablw[ns][nr] = refrmin + nr * drefr;
                    data.refitablw[ns][ni] = refimin / 0.2 *std::pow(0.2, ni + 1);
                    if (ni == nrefi - 1) {
                        data.refitablw[ns][ni] = -1e-21;
                    }
                    crefd= std::complex<real>(data.refrtablw[ns][nr], data.refitablw[ns][ni]);

                    for (int n = 0; n < nsiz; n++) {
                        xr = cos(pie * (static_cast<real>(n) + 0.5) / static_cast<real>(nsiz));
                        rs[n] = std::exp(xr * bma + bpa);

                        thesize = 2.0 * pie * rs[n] / data.wavmidlw[ns];
                        thesize = std::min(static_cast<double>(thesize), 10000.0);

                       miev0(thesize,crefd,perfct,mimcut,anyang,   \
                                numang,xmu,nmom,ipolzn,momdim,prnt,   \
                                qext,qsca,gqsc,n,pmom,sforw,sback,s1,  \
                                s2,tforw,tback );

                        qext4[n] = qext[n];
                        qext4[n] = std::max(static_cast<double>(qext4[n]), 1e-20);
                        qsca4[n] = std::min(static_cast<double>(qsca[n]), static_cast<double>(qext[n]));
                        qsca4[n] = std::max(static_cast<double>(qsca4[n]), 1e-20);
                        qabs4[n] = qext4[n] - qsca4[n];
                        qabs4[n] = std::max(static_cast<double>(qabs4[n]), 1e-20);
                        asymm[n] = gqsc[n]/qsca4[n];
                        asymm[n] = std::max(static_cast<double>(asymm[n]), 1e-20);
                    }
                    fitcurv(rs, qext4, data.extplw, nr,ni,ns, nsiz);
                    fitcurv(rs, qabs4, data.absplw, nr,ni,ns, nsiz);
                    fitcurv(rs, qsca4, data.ascatplw, nr,ni,ns, nsiz);
                    fitcurv(rs, asymm, data.asmplw, nr,ni,ns, nsiz);
                }
            }
        }

}
// clang-format on
KOKKOS_INLINE_FUNCTION void mieaer(const int &id, int iclm, const int &nbin_a,
                                   const real3d &number_bin_col, const real3d &radius_wet_col,
                                   const real_com4d &swrefindx_col,
                                   const real_com4d &lwrefindx_col,
                                   const real2d &dz,
                                   const real &curr_secs, const int &kts, const int &kte,
                                   const real3d &swsizeaer, const real3d &swextaer,
                                   const real3d &swwaer, const real3d &swgaer,
                                   const real3d &swtauaer, const real3d &swabsaer,
                                   const real3d &lwextaer, const real3d &lwtauaer,
                                   const real3d &l2, const real3d &l3,
                                   const real3d &l4, const real3d &l5,
                                   const real3d &l6, const real3d &l7,
                                   const real3d &swbscoef) {
    real xrmin = std::log(rmin);
    real xrmax = std::log(rmax);
    // ######################################################################
    // parameterization of mie calculation for shortwave
    // #####################################################################

    // reald thesum1("thesum1",yakl::Bnd {0,kte-1});
    real thesum1(kte);
    real thesum2(kte);

    for (int klevel = 0; klevel < kte; klevel++) {
        thesum1(klevel) = 0.0;
        for (int m = 0; m < nbin_a; m++) {
            thesum1(klevel) += number_bin_col(m, klevel, iclm);
        }
    }

    for (int klevel = 0; klevel < kte; klevel++) {
        for (int ns = 0; ns < nswbands; ns++) {
            real cext[ncoef];
            real casm[ncoef];
            real cscat[ncoef];
            real cpmom2[ncoef];
            real cpmom3[ncoef];
            real cpmom4[ncoef];
            real cpmom5[ncoef];
            real cpmom6[ncoef];
            real cpmom7[ncoef];
            real cpsback2p[ncoef];
            int itab, jtab;
            real ttab, utab;
            swtauaer(ns, klevel, iclm)  = 0.0;
            swabsaer(ns, klevel, iclm)  = 0.0;
            swwaer(ns, klevel, iclm)    = 0.0;
            swgaer(ns, klevel, iclm)    = 0.0;
            swsizeaer(ns, klevel, iclm) = 0.0;
            swextaer(ns, klevel, iclm)  = 0.0;
            l2(ns, klevel, iclm)        = 0.0;
            l3(ns, klevel, iclm)        = 0.0;
            l4(ns, klevel, iclm)        = 0.0;
            l5(ns, klevel, iclm)        = 0.0;
            l6(ns, klevel, iclm)        = 0.0;
            l7(ns, klevel, iclm)        = 0.0;
            swbscoef(ns, klevel, iclm)  = 0.0;

            // if (thesum <= 1e-21) continue;
            if (thesum1(klevel) > 1e-21) {
                for (int m = 0; m < nbin_a; m++) {
                    real sizem = radius_wet_col(m, klevel, iclm);

                    if (radius_wet_col(m, klevel, iclm) <= rmin) {
                        radius_wet_col(m, klevel, iclm) = rmin;
                        // write( msg, '(a, 5i4,1x, e11.4)' )	&
                        //'mieaer: radius_wet set to rmin,'  //	&
                        //'id,i,j,k,m,rm(m,k)', id, iclm, klevel, m, radius_wet_col(m+1,klevel+1);
                        // call chemistry_message( msg )
                    }
                    if (radius_wet_col(m, klevel, iclm) > rmax) {
                        radius_wet_col(m, klevel, iclm) = rmax;
                        if (number_bin_col(m, klevel, iclm) >= 1.e-10) {
                            // write( msg, '(a, 5i4,1x, 2e11.4)' )	&
                            //'mieaer: radius_wet set to rmax,'  //	&
                            //'id,i,j,k,m,rm(m,k),number', &
                            // id, iclm, klevel, m, radius_wet_col(m+1,klevel+1);,number_bin_col(m+1,klevel+1);
                            // call chemistry_message( msg )
                        }
                    }

                    real x                      = std::log(radius_wet_col(m, klevel, iclm));
                    std::complex<double> crefin = swrefindx_col(m, klevel, iclm,ns);
                    real refr                   = std::real(crefin);
                    real refi                   = -std::imag(crefin);
                    real xrad                   = x;
                    real thesize                = 2.0 * pie * std::exp(x) / data.wavmidsw[ns];

                    xrad = (2 * xrad - xrmax - xrmin) / (xrmax - xrmin);

                    if (std::abs(refr) > 10.0 || std::abs(refr) <= 0.001) {
                        // write ( msg, '(a,1x, e14.5)' )  &
                        //'mieaer /refr/ outside sw range 1e-3 - 10 ' //  &
                        //'refr= ', refr
                        // call chemistry_error_fatal( msg )
                    }
                    if (std::abs(refi) > 10.0) {
                        // write ( msg, '(a,1x, e14.5)' )  &
                        //'mieaer /refi/ sw >10 '  //  &
                        //'refi', refi
                        // call chemistry_error_fatal( msg )
                    }

                    itab = 0;

                    binterp_s(data.extpsw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cext);
                    binterp_s(data.ascatpsw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cscat);
                    binterp_s(data.asmpsw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, casm);
                    binterp_s(data.pmom2psw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cpmom2);
                    binterp_s(data.pmom3psw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cpmom3);
                    binterp_s(data.pmom4psw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cpmom4);
                    binterp_s(data.pmom5psw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cpmom5);
                    binterp_s(data.pmom6psw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cpmom6);
                    binterp_s(data.pmom7psw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cpmom7);
                    binterp_s(data.sbackpsw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cpsback2p);

                    real ch[ncoef];
                    ch[0]       = 1.0;
                    ch[1]       = xrad;
                    real pext   = 0.5 * cext[0];
                    real pscat  = 0.5 * cscat[0];
                    real pasm   = 0.5 * casm[0];
                    real ppmom2 = 0.5 * cpmom2[0];
                    real ppmom3 = 0.5 * cpmom3[0];
                    real ppmom4 = 0.5 * cpmom4[0];
                    real ppmom5 = 0.5 * cpmom5[0];
                    real ppmom6 = 0.5 * cpmom6[0];
                    real ppmom7 = 0.5 * cpmom7[0];
                    real sback2 = 0.5 * cpsback2p[0];

                    for (int nc = 1; nc < ncoef; nc++) {
                        if (nc > 1) {
                            ch[nc] = 2.0 * xrad * ch[nc - 1] - ch[nc - 2];
                        }

                        pext   = pext + ch[nc] * cext[nc];
                        pscat  = pscat + ch[nc] * cscat[nc];
                        pasm   = pasm + ch[nc] * casm[nc];
                        ppmom2 = ppmom2 + ch[nc] * cpmom2[nc];
                        ppmom3 = ppmom3 + ch[nc] * cpmom3[nc];
                        ppmom4 = ppmom4 + ch[nc] * cpmom4[nc];
                        ppmom5 = ppmom5 + ch[nc] * cpmom5[nc];
                        ppmom6 = ppmom6 + ch[nc] * cpmom6[nc];
                        ppmom7 = ppmom7 + ch[nc] * cpmom7[nc];
                        sback2 = sback2 + ch[nc] * cpsback2p[nc];
                    }

                    pext  = std::exp(pext);
                    pscat = std::exp(pscat);
                    pasm  = std::exp(pasm);
                    if (ppmom2 <= 0.0) ppmom2 = 0.0;
                    if (ppmom3 <= 0.0) ppmom3 = 0.0;
                    if (ppmom4 <= 0.0 || sizem <= 0.03e-04) ppmom4 = 0.0;
                    if (ppmom5 <= 0.0 || sizem <= 0.03e-04) ppmom5 = 0.0;
                    if (ppmom6 <= 0.0 || sizem <= 0.03e-04) ppmom6 = 0.0;
                    if (ppmom7 <= 0.0 || sizem <= 0.03e-04) ppmom7 = 0.0;
                    sback2 = std::exp(sback2);
                    if (sback2 <= 0.0) sback2 = 0.0;

                    // weights:s
                    pscat                                    = std::min(static_cast<double>(pscat), static_cast<double>(pext));  // czhao
                    real weighte                             = pext * pie * std::pow(std::exp(x), 2);                            // JCB, extinction cross section
                    real weights                             = pscat * pie * std::pow(std::exp(x), 2);                           // JCB, scattering cross section
                    swtauaer(ns, klevel, iclm) = swtauaer(ns, klevel, iclm) + weighte * number_bin_col(m, klevel, iclm);
                    // must be multiplied by deltaZ
                    // czhao add absorption optical depth
                    swabsaer(ns, klevel, iclm)  = swabsaer(ns, klevel, iclm) + (weighte - weights) * number_bin_col(m, klevel, iclm);
                    swsizeaer(ns, klevel, iclm) = swsizeaer(ns, klevel, iclm) + exp(x) * 10000.0 * number_bin_col(m, klevel, iclm);
                    swwaer(ns, klevel, iclm)    = swwaer(ns, klevel, iclm) + weights * number_bin_col(m, klevel, iclm);
                    swgaer(ns, klevel, iclm)    = swgaer(ns, klevel, iclm) + pasm * weights * number_bin_col(m, klevel, iclm);
                    // need weighting by scattering cross section ?  JCB 2004/02/09
                    l2(ns, klevel, iclm)       = l2(ns, klevel, iclm) + weights * ppmom2 * number_bin_col(m, klevel, iclm);
                    l3(ns, klevel, iclm)       = l3(ns, klevel, iclm) + weights * ppmom3 * number_bin_col(m, klevel, iclm);
                    l4(ns, klevel, iclm)       = l4(ns, klevel, iclm) + weights * ppmom4 * number_bin_col(m, klevel, iclm);
                    l5(ns, klevel, iclm)       = l5(ns, klevel, iclm) + weights * ppmom5 * number_bin_col(m, klevel, iclm);
                    l6(ns, klevel, iclm)       = l6(ns, klevel, iclm) + weights * ppmom6 * number_bin_col(m, klevel, iclm);
                    l7(ns, klevel, iclm)       = l7(ns, klevel, iclm) + weights * ppmom7 * number_bin_col(m, klevel, iclm);
                    swbscoef(ns, klevel, iclm) = swbscoef(ns, klevel, iclm) + pie * exp(x) * exp(x) * sback2 * number_bin_col(m, klevel, iclm);
                }
                // take averages - weighted by cross section - new code JCB 2004/02/09
                swsizeaer(ns, klevel, iclm) = swsizeaer(ns, klevel, iclm) / thesum1(klevel);
                swgaer(ns, klevel, iclm)    = swgaer(ns, klevel, iclm) / swwaer(ns, klevel, iclm);  // JCB removed *3 factor 2/9/2004
                // because factor is applied in subroutine opmie, file zz01fastj_mod.f
                l2(ns, klevel, iclm) = l2(ns, klevel, iclm) / swwaer(ns, klevel, iclm);
                l3(ns, klevel, iclm) = l3(ns, klevel, iclm) / swwaer(ns, klevel, iclm);
                l4(ns, klevel, iclm) = l4(ns, klevel, iclm) / swwaer(ns, klevel, iclm);
                l5(ns, klevel, iclm) = l5(ns, klevel, iclm) / swwaer(ns, klevel, iclm);
                l6(ns, klevel, iclm) = l6(ns, klevel, iclm) / swwaer(ns, klevel, iclm);
                l7(ns, klevel, iclm) = l7(ns, klevel, iclm) / swwaer(ns, klevel, iclm);
                // backscatter coef, divide by 4*Pie to get units of (km*ster)^-1 JCB 2007/02/01
                swbscoef(ns, klevel, iclm) = swbscoef(ns, klevel, iclm) * 1.0e5;  // units are now (km)^-1
                swextaer(ns, klevel, iclm) = swtauaer(ns, klevel, iclm) * 1.0e5;  // now true extincion, units (km)^-1
                // this must be last
                swwaer(ns, klevel, iclm) = swwaer(ns, klevel, iclm) / swtauaer(ns, klevel, iclm);
            }
            //}
            //  });
        }
    }

    // before returning, multiply tauaer by depth of individual cells.
    // tauaer is in cm-1, dz in m; multiply dz by 100 to convert from m to cm.
    for (int ns = 0; ns < data.nswbands; ns++) {
        for (int klevel = 0; klevel < kte; klevel++) {
            swtauaer(ns, klevel, iclm) = swtauaer(ns, klevel, iclm) * dz( klevel,iclm) * 100;
            swabsaer(ns, klevel, iclm) = swabsaer(ns, klevel, iclm) * dz( klevel,iclm) * 100;
        }
    }

    // ######################################################################
    // parameterization of mie calculation for longwave
    // #####################################################################
    // reald thesum2("thesum2",yakl::Bnd {0,kte-1});
    for (int klevel = 0; klevel < kte; klevel++) {
        thesum2(klevel) = 0.0;
        for (int m = 0; m < nbin_a; m++) {
            thesum2(klevel) += number_bin_col(m, klevel, iclm);
        }
    }

    for (int klevel = 0; klevel < kte; klevel++) {
        for (int ns = 0; ns < nlwbands; ns++) {
            // parallel_for("mieaer2",Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {kte, nlwbands}),
            //     KOKKOS_LAMBDA(int klevel, int ns){
            real cabs[ncoef];
            int itab, jtab;
            real ttab, utab;
            // Aerosol optical properties
            lwtauaer(ns, klevel, iclm) = 0.0;
            lwextaer(ns, klevel, iclm) = 0.0;

            if (thesum2(klevel) > 1e-21) {
                // Loop over the bins
                for (int m = 0; m < nbin_a; m++) {
                    // Size
                    real sizem                  = radius_wet_col(m, klevel, iclm);  // Radius in cm
                    real x                      = std::log(static_cast<double>(sizem));           // Radius in cm
                    std::complex<double> crefin = lwrefindx_col(m, klevel, iclm,ns);
                    real refr                   = crefin.real();
                    real refi                   = -crefin.imag();
                    real xrad                   = x;
                    real thesize                = 2.0 * pie * std::exp(x) / data.wavmidlw[ns];

                    // Normalize size parameter
                    xrad = (2 * xrad - xrmax - xrmin) / (xrmax - xrmin);
                    // cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    // ec  diagnostics
                    // ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                    //  retain this diagnostic code
                    if (std::abs(refr) > 10.0 || std::abs(refr) <= 0.001) {
                        // write ( msg, '(a,1x, e14.5)' )  &
                        //'mieaer /refr/ outside lw range 1e-3 - 10 ' //  &
                        //'refr= ', refr
                        // call chemistry_error_fatal( msg )
                    }
                    if (std::abs(refi) > 10.0) {
                        // write ( msg, '(a,1x, e14.5)' )  &
                        //'mieaer /refi/ lw >10 '  //  &
                        //'refi', refi
                        // call chemistry_error_fatal( msg )
                    }

                    // Interpolate coefficients linear in refractive index

                    itab = 0;

                    binterp_l(data.absplw, ns, nrefr, nrefi, refr, refi, itab, jtab, ttab, utab, cabs);

                    // Chebyshev polynomials
                    real ch[ncoef];
                    ch[0] = 1.0;
                    ch[1] = xrad;
                    for (int nc = 2; nc < ncoef; nc++) {
                        ch[nc] = 2.0 * xrad * ch[nc - 1] - ch[nc - 2];
                    }

                    // Parameterized optical properties
                    real pabs = 0.5 * cabs[0];
                    for (int nc = 1; nc < ncoef; nc++) {
                        pabs += ch[nc] * cabs[nc];
                    }
                    pabs = std::exp(pabs);

                    // Weights
                    real weighta = pabs * pie * std::pow(std::exp(x), 2);  // Extinction cross section
                    lwtauaer(ns, klevel, iclm) += weighta * number_bin_col(m, klevel, iclm);

                }  // End of nbin_a loop

                // Take averages - weighted by cross section
                lwextaer(ns, klevel, iclm) = lwtauaer(ns, klevel, iclm) * 1.0e5;  // True extinction, units (km)^-1
            }
            //  }); // End of wavelength loop
        }
    }

    for (int ns = 0; ns < data.nlwbands; ns++) {
        for (int klevel = 0; klevel < kte; klevel++) {
            lwtauaer(ns, klevel, iclm) = lwtauaer(ns, klevel, iclm) * dz( klevel,iclm) * 100;
        }
    }
}
// clang-format off
KOKKOS_INLINE_FUNCTION void binterp_s(real (*table)[prefi][prefr][ncoef],int ns,int &im,int &jm,real &x,real &y,\
                int &ix, int &jy, real &t, real &u,real out[]) {
    //bilinear interpolation of table

    int ip1, jp1;
    int i, j, k;
    real dx, dy, tu, tuc, tcu, tcuc;


  if (ix <= 0) {
    if (im > 1) {
      for (i = 0; i < im; ++i) {
        if (x < data.refrtabsw[ns][i])
          break;
      }
      ix = std::max(i, 1);
      ip1 = std::min(ix + 1, im);
      dx = data.refrtabsw[ns][ip1 - 1] - data.refrtabsw[ns][ix - 1];
      if (std::abs(dx) > 1.e-20) {
        t = (x - data.refrtabsw[ns][ix - 1]) /
            (data.refrtabsw[ns][ix] - data.refrtabsw[ns][ix - 1]);
      } else {
        t = 0;
      }
    } else {
      ix = 1;
      ip1 = 1;
      t = 0;
    }
    if (jm > 1) {
      for (j = 0; j < jm; ++j) {
        if (y < data.refitabsw[ns][j])
          break;
      }
      jy = std::max(j, 1);
      jp1 = std::min(jy + 1, jm);
      dy = data.refitabsw[ns][jp1 - 1] - data.refitabsw[ns][jy - 1];
      if (std::abs(dy) > 1.e-20) {
        u = (y - data.refitabsw[ns][jy - 1]) / dy;
      } else {
        u = 0;
      }
    } else {
      jy = 1;
      jp1 = 1;
      u = 0;
    }
  }

  jp1 = std::min(jy + 1, jm);
  ip1 = std::min(ix + 1, im);
  tu = t * u;
  tuc = t - tu;
  tcuc = 1 - tuc - u;
  tcu = u - tu;
  for (k = 0; k < ncoef; ++k) {
    out[k] = tcuc * table[ns][jy - 1][ix - 1][k] +
             tuc * table[ns][jy - 1][ip1 - 1][k] +
             tu * table[ns][jp1 - 1][ip1 - 1][k] +
             tcu * table[ns][jp1 - 1][ix - 1][k];
  }
  return;
}

KOKKOS_INLINE_FUNCTION void binterp_l(real (*table)[prefi][prefr][ncoef],int ns,int &im,int &jm,real &x,real &y,\
                int &ix, int &jy, real &t, real &u,real out[]) {
    //bilinear interpolation of table

    int ip1, jp1;
    int i, j, k;
    real dx, dy, tu, tuc, tcu, tcuc;


  if (ix <= 0) {
    if (im > 1) {
      for (i = 0; i < im; ++i) {
        if (x < data.refrtablw[ns][i])
          break;
      }
      ix = std::max(i, 1);
      ip1 = std::min(ix + 1, im);
      dx = data.refrtablw[ns][ip1 - 1] - data.refrtablw[ns][ix - 1];
      if (std::abs(dx) > 1.e-20) {
        t = (x - data.refrtablw[ns][ix - 1]) /
            (data.refrtablw[ns][ix] - data.refrtablw[ns][ix - 1]);
      } else {
        t = 0;
      }
    } else {
      ix = 1;
      ip1 = 1;
      t = 0;
    }
    if (jm > 1) {
      for (j = 0; j < jm; ++j) {
        if (y < data.refitablw[ns][j])
          break;
      }
      jy = std::max(j, 1);
      jp1 = std::min(jy + 1, jm);
      dy = data.refitablw[ns][jp1 - 1] - data.refitablw[ns][jy - 1];
      if (std::abs(dy) > 1.e-20) {
        u = (y - data.refitablw[ns][jy - 1]) / dy;
      } else {
        u = 0;
      }
    } else {
      jy = 1;
      jp1 = 1;
      u = 0;
    }
  }

  jp1 = std::min(jy + 1, jm);
  ip1 = std::min(ix + 1, im);
  tu = t * u;
  tuc = t - tu;
  tcuc = 1 - tuc - u;
  tcu = u - tu;
  for (k = 0; k < ncoef; ++k) {
    out[k] = tcuc * table[ns][jy - 1][ix - 1][k] +
             tuc * table[ns][jy - 1][ip1 - 1][k] +
             tu * table[ns][jp1 - 1][ip1 - 1][k] +
             tcu * table[ns][jp1 - 1][ix - 1][k];
  }
  return;
}



KOKKOS_INLINE_FUNCTION  void fitcurv(real *rs,real *yin,real (*coef)[prefi][prefr][ncoef],\
                int nr,int ni,int ns,const int maxm){
 //  fit y(x) using Chebychev polynomials
     real x[maxm],y[maxm];
     int m = 0;
     real xmin,xmax;
 //    char msg[300];
     for (m = 0; m < maxm; ++m) {
         x[m] = std::log(rs[m]);
         y[m] = std::log(yin[m]);
     }

     xmin = x[0];
     xmax = x[maxm - 1];
     for (m = 0; m < maxm; ++m) {
         x[m] = (2 * x[m] - xmax - xmin) / (xmax - xmin);
     }

     chebft(coef, nr, ni, ns, maxm, y);
     return;
}

KOKKOS_INLINE_FUNCTION void fitcurv_nolog(real *rs,real *yin,real (*coef)[prefi][prefr][ncoef],\
                int nr,int ni,int ns,const int maxm){
 //  fit y(x) using Chebychev polynomials
     real x[maxm],y[maxm];
     int m = 0;
     real xmin,xmax;
 //    char msg[300];
     for (m = 0; m < maxm; ++m) {
         x[m] = std::log(rs[m]);
         y[m] = yin[m]; // note, no "log" here
     }

     xmin = x[0];
     xmax = x[maxm - 1];
     for (m = 0; m < maxm; ++m) {
         x[m] = (2 * x[m] - xmax - xmin) / (xmax - xmin);
     }

     chebft(coef,nr, ni, ns, maxm, y);
     return;
}

KOKKOS_INLINE_FUNCTION  void chebft(real (*c)[prefi][prefr][ncoef],\
                int nr,int ni,int ns,int n,real *f){
 /**
     given a function f with values at zeroes x_k of Chebychef polynomial
     T_n(x), calculate coefficients c_j such that
     f(x)=sum(k=1,n) c_k t_(k-1)(y) - 0.5*c_1
     where y=(x-0.5*(xmax+xmin))/(0.5*(xmax-xmin))
     See Numerical Recipes, pp. 148-150.
 **/

 // local variables
     real fac,thesum;
     int j,k;
     fac = 2.0 / n;
     for (j = 1; j <= ncoef; ++j) {
         thesum = 0;
         for (k = 1; k <= n; ++k) {
             thesum += f[k-1] * std::cos((pi * (j - 1)) * ((k - 0.5) / n));
         }
         c[ns][ni][nr][j - 1] = fac * thesum;
     }
     return;
}



KOKKOS_INLINE_FUNCTION void miev0(real &xx,complex<double> &crefin,bool &perfct,real &mimcut,bool &anyang,int &numang,real *xmu,int &nmom,\
            int &ipolzn,int &momdim,bool *prnt,real *qext,real *qsca,real *gqsc,int n_,real pmom[8][8],complex<double> &sforw,  \
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback){
/**
    computes mie scattering and extinction efficiencies; asymmetry
    factor;  forward- and backscatter amplitude;  scattering
    amplitudes for incident polarization parallel and perpendicular
    to the plane of scattering, as functions of scattering angle;
    coefficients in the legendre polynomial expansions of either the
    unpolarized phase function or the polarized phase matrix;
    and some quantities needed in polarized radiative transfer.

    calls :  biga, ckinmi, small1, small2, testmi, miprnt,
               lpcoef, errmsg

    i n t e r n a l   v a r i a b l e s
    -----------------------------------

  an,bn           mie coefficients  little-a-sub-n, little-b-sub-n
                     ( ref. 1, eq. 16 )
  anm1,bnm1       mie coefficients  little-a-sub-(n-1),
                     little-b-sub-(n-1);  used in -gqsc- sum
  anp             coeffs. in s+ expansion ( ref. 2, p. 1507 )
  bnp             coeffs. in s- expansion ( ref. 2, p. 1507 )
  anpm            coeffs. in s+ expansion ( ref. 2, p. 1507 )
                     when  mu  is replaced by  - mu
  bnpm            coeffs. in s- expansion ( ref. 2, p. 1507 )
                     when  mu  is replaced by  - mu
  calcmo(k)       true, calculate moments for k-th phase quantity
                     (derived from -ipolzn-; used only in 'lpcoef')
  cbiga(n)        bessel function ratio capital-a-sub-n (ref. 2, eq. 2)
                     ( complex version )
  cior            complex index of refraction with negative
                     imaginary part (van de hulst convention)
  cioriv          1 / cior
  coeff           ( 2n + 1 ) / ( n ( n + 1 ) )
  fn              realing point version of index in loop performing
                     mie series summation
  lita,litb(n)    mie coefficients -an-, -bn-, saved in arrays for
                     use in calculating legendre moments *pmom*
  maxtrm          max. possible no. of terms in mie series
  mm              + 1 and  - 1,  alternately.
  mim             magnitude of imaginary refractive index
  mre             real part of refractive index
  maxang          max. possible value of input variable -numang-
  nangd2          (numang+1)/2 ( no. of angles in 0-90 deg; anyang=f )
  noabs           true, sphere non-absorbing (determined by -mimcut-)
  np1dn           ( n + 1 ) / n
  npquan          highest-numbered phase quantity for which moments are
                     to be calculated (the largest digit in -ipolzn-
                     if  ipolzn .ne. 0)
  ntrm            no. of terms in mie series
  pass1           true on first entry, false thereafter; for self-test
  pin(j)          angular function little-pi-sub-n ( ref. 2, eq. 3 )
                     at j-th angle
  pinm1(j)        little-pi-sub-(n-1) ( see -pin- ) at j-th angle
  psinm1          ricatti-bessel function psi-sub-(n-1), argument -xx-
  psin            ricatti-bessel function psi-sub-n of argument -xx-
                     ( ref. 1, p. 11 ff. )
  rbiga(n)        bessel function ratio capital-a-sub-n (ref. 2, eq. 2)
                     ( real version, for when imag refrac index = 0 )
  rioriv          1 / mre
  rn              1 / n
  rtmp            (real) temporary variable
  sp(j)           s+  for j-th angle  ( ref. 2, p. 1507 )
  sm(j)           s-  for j-th angle  ( ref. 2, p. 1507 )
  sps(j)          s+  for (numang+1-j)-th angle ( anyang=false )
  sms(j)          s-  for (numang+1-j)-th angle ( anyang=false )
  taun            angular function little-tau-sub-n ( ref. 2, eq. 4 )
                     at j-th angle
  tcoef           n ( n+1 ) ( 2n+1 ) (for summing tforw,tback series)
  twonp1          2n + 1
  yesang          true if scattering amplitudes are to be calculated
  zetnm1          ricatti-bessel function  zeta-sub-(n-1) of argument
                     -xx-  ( ref. 2, eq. 17 )
  zetn            ricatti-bessel function  zeta-sub-n of argument -xx-
 */
 /*********************************************************************
 ***********  i / o specifications for subroutine miev0  **************
 *********************************************************************/
    const int maxang = 501;
    const int mxang2 = maxang/2 + 1;
/**
   note --  maxtrm = 10100  is neces-
   sary to do some of the test probs,
   but 1100 is sufficient for most
   conceivable applications
   jwfeng: 2022/03/22, changed to 1500
   because 14 bands need more terms.
*/
 constexpr int maxtrm = 1500;
  constexpr real onethr = 1. / 3.;
  string msg = "";
  static bool anysav, calcmo[4], noabs, ok, persav, yesang;
  int npquan;
  int i, j, n, nmosav, iposav, numsav, ntrm, nangd2;
  real mim, mimsav, mre, mm, np1dn;
  bool noabs_;
  real rioriv, xmusav, xxsav, sq, fn, rn, twonp1, tcoef, coeff;
  real xinv, psinm1, chinm1, psin, chin, rtmp, taun;
  real pin[maxang], pinm1[maxang], rbiga[maxtrm];
  complex<double> an, bn, anm1, bnm1, anp, bnp, anpm, bnpm, cresav, cior,
      cioriv, ctmp, zet, zetnm1, zetn;
  complex<double> lita[maxtrm], litb[maxtrm], cbiga[maxtrm], sp[maxang],
      sm[maxang], sps[mxang2], sms[mxang2];
  static bool pass1 = true;

  static int flag = 0;
  union biga_c_r {
    real rbiga[maxtrm];
    complex<double> cbiga[maxtrm];
    biga_c_r(){};
  } c_r;

  if (pass1) {
    // save certain user input values
    xxsav = xx;
    cresav = crefin;
    mimsav = mimcut;
    persav = perfct;
    anysav = anyang;
    nmosav = nmom;
    iposav = ipolzn;
    numsav = numang;
    xmusav = xmu[0];
    // reset input values for test case
    xx = 10.0;
    crefin = complex<double>(1.5, -0.1);
    perfct = false;
    mimcut = 0.0;
    anyang = true;
    numang = 1;
    xmu[0] = -0.7660444;
    nmom = 1;
    ipolzn = -1;
  }
  /** check input and calculate
      certain variables from input **/
  // label10:
  // bool bkpt = true;
  // while (bkpt) {}
  do {
    flag = flag + 1;
    ckinmi(numang, maxang, xx, perfct, crefin, momdim, nmom, ipolzn, anyang,
           xmu, calcmo, npquan);

    if (!perfct) {
      cior = crefin;
      if (std::imag(cior) > 0.0)
      cior = std::conj(cior);
      mre = std::real(cior);
    //   std::cout<<"cior"<<cior<<std::endl;
    //   std::cout<<"mre2"<<mre<<std::endl;
      mim = -std::imag(cior);
      if (mim <= mimcut) {
        noabs = true;
      } else {
        noabs = false;
      }
      cioriv = 1.0 / cior;
    //   std::cout<<"mre2"<<mre<<std::endl;
      rioriv = 1.0 / mre;
    }
    if (perfct && xx <= 0.1) {
      /** use totally-reflecting
        small-particle limit   **/
      small1(xx, numang, xmu, qext, qsca, gqsc, n_, sforw, sback, s1, s2, tforw,
             tback, lita, litb);
      ntrm = 2;
      // goto label200;
    } else if (!perfct && xx * std::max(1.0, std::abs(cior)) <= 0.0) {
      /**
       * use general-refractive-index
          small-particle limit
          ( ref. 2, p. 1508 )
      **/
      noabs_ = !noabs;
      small2(xx, cior, noabs_, numang, xmu, qext, qsca, gqsc, n_, sforw, sback,
             s1, s2, tforw, tback, lita, litb);
      ntrm = 2;
      //   goto label200;
    } else {
      nangd2 = (numang + 1) / 2;
      yesang = numang > 0;
      /**  estimate number of terms in mie series
       ( ref. 2, p. 1508 ) **/
      if (xx <= 8.0) {
        ntrm = xx + 4.0 * std::pow(xx, onethr) + 1.0;
      } else if (xx < 4200.0) {
        ntrm = xx + 4.05 * std::pow(xx, onethr) + 2.0;
      } else {
        ntrm = xx + 4.0 * std::pow(xx, onethr) + 2.0;
      }
      if (ntrm + 1 > maxtrm) {
        // std::ostringstream oss;
        // oss << "ntrm = " << ntrm << " xx = " << xx;
        // msg = oss.str();
        // chemistry_error_message(trim(msg))
        errmsg("miev0--parameter maxtrm too small", true);
      }
      if (!perfct)
        biga(cior, xx, ntrm, noabs, yesang, c_r.rbiga, c_r.cbiga);
      xinv = 1.0 / xx;
      psinm1 = std::sin(xx);
      chinm1 = std::cos(xx);
      psin = psinm1 * xinv - chinm1;
      chin = chinm1 * xinv + psinm1;
      zetnm1 = complex<double>(psinm1, chinm1);
      zetn = complex<double>(psin, chin);
      // initialize previous coeffi-
      // cients for -gqsc- series
      anm1 = complex<double>(0.0, 0.0);
      bnm1 = complex<double>(0.0, 0.0);
      // initialize angular function little-pi
      // and sums for s+, s- ( ref. 2, p. 1507 )
      if (anyang) {
        for (j = 0; j < numang; ++j) {
          pinm1[j] = 0.0;
          pin[j] = 1.0;
          sp[j] = complex<double>(0.0, 0.0);
          sm[j] = complex<double>(0.0, 0.0);
        }
      } else {
        for (j = 0; j < nangd2; ++j) {
          pinm1[j] = 0.0;
          pin[j] = 1.0;
          sp[j] = complex<double>(0.0, 0.0);
          sm[j] = complex<double>(0.0, 0.0);
          sps[j] = complex<double>(0.0, 0.0);
          sms[j] = complex<double>(0.0, 0.0);
        }
      }
      // initialize mie sums for efficiencies, etc.
      qsca[n_] = 0.0;
      gqsc[n_] = 0.0;
      sforw = complex<double>(0., 0.);
      sback = complex<double>(0., 0.);
      tforw[0] = complex<double>(0., 0.);
      tback[0] = complex<double>(0., 0.);
      /**********  loop to sum mie series  ***************/
      mm = 1.0;
    for(int n = 0;n<ntrm;n++){
        /** compute various numerical coefficients **/
        fn = n + 1;
        rn = 1.0 / fn;
        np1dn = 1.0 + rn;
        twonp1 = 2 * (n + 1) + 1;
        coeff = twonp1 / (fn * (n + 2));
        tcoef = twonp1 * (fn * (n + 2));
        /** calculate mie series coefficients  **/
        if (perfct) {
          /** totally-reflecting case  **/
          an = complex<double>((fn * xinv) * psin - psinm1, 0.0) /
               (complex<double>(fn * xinv, 0.0) * zetn - zetnm1);
          bn = complex<double>(psin, 0.0) / zetn;
        } else if (noabs) {
          /** no-absorption case  **/
          an = (complex<double>((rioriv * c_r.rbiga[n] + (fn * xinv)) * psin -
                                psinm1,
                0.0)) /
               (complex<double>(rioriv * c_r.rbiga[n] + (fn * xinv), 0.0) *
                    zetn -
                zetnm1);
          bn = (complex<double>((mre * c_r.rbiga[n] + (fn * xinv)) * psin -
                                psinm1,
                0.0)) /
               (complex<double>(mre * c_r.rbiga[n] + (fn * xinv), 0.0) * zetn -
                zetnm1);
        } else {
          /** absorptive case  **/
          an = ((cioriv * c_r.cbiga[n] + complex<double>(fn * xinv, 0.0)) *
                    complex<double>(psin, 0.0) -
                complex<double>(psinm1, 0.0)) /
               ((cioriv * c_r.cbiga[n] + complex<double>(fn * xinv, 0.0)) *
                    zetn -
                zetnm1);
          bn = ((cior * c_r.cbiga[n] + complex<double>(fn * xinv, 0.0)) *
                    complex<double>(psin, 0.0) -
                complex<double>(psinm1, 0.0)) /
               ((cior * c_r.cbiga[n] + complex<double>(fn * xinv, 0.0)) * zetn -
                zetnm1);
          qsca[n_] += twonp1 * (std::norm(an) + std::norm(bn));
        }
        /** save mie coefficients for *pmom* calculation **/
        lita[n] = an;
        litb[n] = bn;
        /** increment mie sums for non-angle-
                 dependent quantities **/
        sforw = sforw + complex<double>(twonp1, 0.0) * (an + bn);
        tforw[0] += complex<double>(tcoef, 0.0) * (an - bn);
        sback = sback + complex<double>(mm * twonp1, 0.0) * (an - bn);
        tback[0] += complex<double>(mm * tcoef, 0.0) * (an + bn);
        gqsc[n_] +=
            (fn - rn) * std::real(anm1 * std::conj(an) + bnm1 * std::conj(bn)) +
            coeff * std::real(an * std::conj(bn));
        if (yesang) {
          /** put mie coefficients in form
          needed for computing s+, s-
          ( ref. 2, p. 1507 )      **/
          anp = complex<double>(coeff, 0.0) * (an + bn);
          bnp = complex<double>(coeff, 0.0) * (an - bn);
          /** increment mie sums for s+, s-
          while upward recursing
          angular functions little pi
          and little tau       **/
          if (anyang) {
            /** arbitrary angles
            vectorizable loop **/
            for (j = 0; j < numang; ++j) {
              rtmp = (xmu[j] * pin[j]) - pinm1[j];
              taun = fn * rtmp - pinm1[j];
              sp[j] += anp * complex<double>(pin[j] + taun, 0.);
              sm[j] += bnp * complex<double>(pin[j] - taun, 0.);
              pinm1[j] = pin[j];
              pin[j] = (xmu[j] * pin[j]) + np1dn * rtmp;
            }

          } else {
            /** angles symmetric about 90 degrees **/
            anpm = complex<double>(mm, 0.0) * anp;
            bnpm = complex<double>(mm, 0.0) * bnp;
            /** vectorizable loop **/
            for (j = 0; j < nangd2; ++j) {
              rtmp = (xmu[j] * pin[j]) - pinm1[j];
              taun = fn * rtmp - pinm1[j];
              sp[j] += anp * complex<double>(pin[j] + taun, 0.0);
              sms[j] += bnpm * complex<double>(pin[j] + taun, 0.0);
              sm[j] += bnp * complex<double>(pin[j] - taun, 0.0);
              sps[j] += anpm * complex<double>(pin[j] - taun, 0.0);
              pinm1[j] = pin[j];
              pin[j] = (xmu[j] * pin[j]) + np1dn * rtmp;
            }
          }
        }
        /** update relevant quantities for next
        pass through loop **/
        mm = -mm;
        anm1 = an;
        bnm1 = bn;
        /** upward recurrence for ricatti-bessel
        functions ( ref. 1, eq. 17 ) **/
        zet = complex<double>(twonp1 * xinv, 0.0) * zetn - zetnm1;
        zetnm1 = zetn;
        zetn = zet;
        psinm1 = psin;
        psin = std::real(zetn);
      }

      /**********  end loop to sum mie series  ***************/
      qext[n_] = 2. / std::pow(xx, 2) * std::real(sforw);
      if (perfct || noabs) {
        qsca[n_] = qext[n_];
      } else {
        qsca[n_] = 2. / std::pow(xx, 2) * qsca[n_];
      }
      gqsc[n_] = 4. / std::pow(xx, 2) * gqsc[n_];
      sforw = 0.5 * sforw;
      sback = 0.5 * sback;
      tforw[1] = 0.5 * (sforw + 0.25 * tforw[0]);
      tforw[0] = 0.5 * (sforw - 0.25 * tforw[0]);
      tback[1] = 0.5 * (sback + 0.25 * tback[0]);
      tback[0] = 0.5 * (-sback + 0.25 * tback[0]);
      if (yesang) {
        /** recover scattering amplitudes
        from s+, s- ( ref. 1, eq. 11 ) **/
        if (anyang) {
          /**vectorizable loop **/
          for (j = 0; j < numang; ++j) {
            s1[j] = 0.5 * (sp[j] + sm[j]);
            s2[j] = 0.5 * (sp[j] - sm[j]);
          }
        } else {
          /**vectorizable loop **/
          for (j = 0; j < nangd2; ++j) {
            s1[j] = 0.5 * (sp[j] + sm[j]);
            s2[j] = 0.5 * (sp[j] - sm[j]);
          }
          /**vectorizable loop **/
          for (j = 0; j < nangd2; ++j) {
            s1[numang - 1 - j] = 0.5 * (sps[j] + sms[j]);
            s2[numang - 1 - j] = 0.5 * (sps[j] - sms[j]);
          }
        }
      }
    }
    /**calculate legendre moments **/
    if (nmom > 0)
      lpcoef(ntrm, nmom, ipolzn, momdim, calcmo, npquan, lita, litb, pmom);
    if (std::imag(crefin) > 0.0) {
      /** take complex conjugates
        of scattering amplitudes **/
      sforw = std::conj(sforw);
      sback = std::conj(sback);
      for (i = 0; i < 2; ++i) {
        tforw[i] = std::conj(tforw[i]);
        tback[i] = std::conj(tback[i]);
      }
      for (j = 0; j < numang; ++j) {
        s1[j] = std::conj(s1[j]);
        s2[j] = std::conj(s2[j]);
      }
    }
    if (pass1) {
      /** compare test case results with
        correct answers and abort if bad **/
      testmi(qext, qsca, gqsc, n_, sforw, sback, s1, s2, tforw, tback, pmom,
             momdim, ok);
      if (!ok) {
        prnt[0] = false;
        prnt[1] = false;
        miprnt(prnt, xx, perfct, crefin, numang, xmu, qext, qsca, gqsc, n_,
               nmom, ipolzn, momdim, calcmo, pmom, sforw, sback, tforw, tback,
               s1, s2);
        errmsg("miev0 -- self-test failed", true);
      }
      /** restore user input values **/
      xx = xxsav;
      crefin = cresav;
      mimcut = mimsav;
      perfct = persav;
      anyang = anysav;
      nmom = nmosav;
      ipolzn = iposav;
      numang = numsav;
      xmu[0] = xmusav;
      pass1 = false;
      // goto label10;
    }
  } while (!pass1 && flag<2);

  if (prnt[0] || prnt[1]) {
    miprnt(prnt, xx, perfct, crefin, numang, xmu, qext, qsca, gqsc, n_, nmom,
           ipolzn, momdim, calcmo, pmom, sforw, sback, tforw, tback, s1, s2);
  }
  return;
}

KOKKOS_INLINE_FUNCTION void ckinmi(int &numang,const int &maxang, real &xx, bool &perfct, complex<double> &crefin,int &momdim, int &nmom, int &ipolzn,bool &anyang, \
            real *xmu,bool *calcmo,int &npquan){
    /*******check for bad input to 'miev0' and calculate -calcmo,npquan-******/
    int i,l,j,ip;
    string str;
    bool inperr = false;
    std::ostringstream oss;
    if(numang > maxang){
        errmsg( "miev0--parameter maxang too small", true );
        inperr = true;
    }
    if(numang < 0) wrtbad( "numang", inperr );
    if(xx < 0) wrtbad( "xx", inperr );
    if(!perfct && std::real(crefin) <= 0.) wrtbad( "crefin", inperr );
    if(momdim < 1.) wrtbad( "momdim", inperr );
    if(nmom != 0){
        if ( nmom < 0 || nmom > momdim ) wrtbad("nmom",inperr);
        if ( std::abs(ipolzn) > 4444 )   wrtbad( "ipolzn", inperr );
        npquan = 0;
        for(l = 0; l < 4; ++l){
            calcmo[l] = false;
        }
        if(ipolzn != 0){
                            /**
                             parse out -ipolzn- into its digits
                             to find which phase quantities are
                             to have their moments calculated
                             */

            oss << std::setw(4) << std::abs(ipolzn);
            str= oss.str();
            for(j = 0; j < 4; ++j){
                ip = str[j] - '0';
                if(ip >= 1 && ip <= 4) calcmo[ip - 1] = true;
                if(ip == 0 || (ip >= 5 && ip <= 9)) wrtbad( "ipolzn", inperr );
                npquan = std::max( npquan, ip );
            }
        }
    }
    if(anyang){
                            /**
                             allow for slight imperfections in
                             computation of cosine
                             */
        for(i = 0; i< numang; ++i){
            if ( xmu[i] < -1.00001 || xmu[i] > 1.00001 ) wrtbad( "xmu", inperr );
        }
        for(i = 0; i< ( numang + 1 ) / 2; ++i){
            if ( xmu[i] < -0.00001 || xmu[i] > 1.00001 ) wrtbad( "xmu", inperr );
        }
    }
    if(inperr) errmsg( "miev0--input error(s).  aborting...", true );
    if(xx > 20000.0 || std::real(crefin) > 10.0 || std::abs(std::imag(crefin)) > 10.0) errmsg(   \
           "miev0--xx or crefin outside tested range", false );
    return;
}
KOKKOS_INLINE_FUNCTION void errmsg(string messag, bool fatal){
    /** print out a warning or error message;  abort if error **/
    static bool once = false;
    static int maxmsg = 100,nummsg = 0;
    if(fatal) //chemistry_error_fatal( 'optical averaging mie fatal error '//messag )
    nummsg += 1;
    if(nummsg > maxmsg){
        //if ( .not.once )  write ( *,99 )
        if(!once) //chemistry_message( 'optical averaging mie: too many warning messages -- no longer printing ' )
        once = true;
    }else{
        //chemistry_message( 'optical averaging mie warning '  // messag ) ;
    }
    return;
}

KOKKOS_INLINE_FUNCTION void wrtbad(string varnam,bool &erflag){
    /**
        write names of erroneous variables
    input :   varnam = name of erroneous variable to be written
                ( character, any length )
    output :  erflag = logical flag, set true by this routine
    **/
    static int maxmsg = 50,nummsg = 0;
    nummsg = nummsg + 1;
//      write ( *, '(3a)' )  ' ****  input variable  ', varnam,   &
//                           '  in error  ****'
    //chemistry_message( 'optical averaging mie input variable in error ' // varnam );
    erflag = true;
    if(nummsg == maxmsg) errmsg ( "too many input variable errors.  aborting...$", true );
    return;
}

KOKKOS_INLINE_FUNCTION void small1(real &xx,int &numang,real *xmu,real *qext,real *qsca,real *gqsc,int n_,complex<double> &sforw,\
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback,\
            complex<double> *a,complex<double> *b){
    /**
    small-particle limit of mie quantities in totally reflecting
    limit ( mie series truncated after 2 terms )
    a,b       first two mie coefficients, with numerator and
              denominator expanded in powers of -xx- ( a factor
              of xx**3 is missing but is restored before return
              to calling program )  ( ref. 2, p. 1508 )
    **/
    int j;
    const real twothr = 2./3.;
    const real fivthr = 5./3.;
    const real fivnin = 5./9.;
    real rtmp;
    complex<double> ctmp;
    a[0] = complex<double> ( 0., twothr * ( 1. - 0.2 * std::pow(xx,2) ) )   \
             / complex<double> ( 1. - 0.5 * std::pow(xx,2), twothr * std::pow(xx,3) );
    b[0] = complex<double> ( 0., - ( 1. - 0.1 * std::pow(xx,2)  ) / 3. )   \
             / complex<double> ( 1.0 + 0.5 * std::pow(xx,2), - std::pow(xx,3) / 3. );
    a[1] = complex<double> ( 0.,   std::pow(xx,2) / 30. );
    b[1] = complex<double> ( 0., - std::pow(xx,2) / 45. );
    qsca[n_] = 6. * std::pow(xx,4) * ( std::norm(a[0]) + std::norm(b[0])    \
                            + fivthr * ( std::norm(a[1]) + std::norm(b[1]) )) ;
    qext[n_] = qsca[n_];
    gqsc[n_] = 6. * std::pow(xx,4) * std::real( a[0] * std::conj( a[1] + b[0] )   \
                          + ( b[0] + complex<double>(fivnin, 0.) * a[1] ) * std::conj( b[1] ) );
    rtmp = 1.5 * std::pow(xx,3);
    sforw      = complex<double>(rtmp, 0.) * ( a[0] + b[0] + complex<double>(fivthr, 0.) * ( a[1] + b[1] ) );
    sback      = complex<double>(rtmp, 0.) * ( a[0] - b[0] - complex<double>(fivthr, 0.) * ( a[1] - b[1] ) );
    tforw[0]   = complex<double>(rtmp, 0.) * ( b[0] + complex<double>(fivthr, 0.) * ( 2.*b[1] - a[1] ) );
    tforw[1]   = complex<double>(rtmp, 0.) * ( a[0] + complex<double>(fivthr, 0.) * ( 2.*a[1] - b[1] ) );
    tback[0]   = complex<double>(rtmp, 0.) * ( b[0] - complex<double>(fivthr, 0.) * ( 2.*b[1] + a[1] ) );
    tback[1]   = complex<double>(rtmp, 0.) * ( a[0] - complex<double>(fivthr, 0.) * ( 2.*a[1] + b[1] ) );
    for(j = 0;j < numang; ++j){
        s1[j] = complex<double>(rtmp, 0.) * ( a[0] + b[0] * complex<double>(xmu[j], 0.) + complex<double>(fivthr, 0.) *   \
                    ( a[1] * complex<double>(xmu[j], 0.) + b[1] * complex<double>( 2.*std::pow(xmu[j],2) - 1. ,0.)) );
        s2[j] = complex<double>(rtmp, 0.) * ( b[0] + a[0] * complex<double>(xmu[j], 0.) + complex<double>(fivthr, 0.) *   \
                    ( b[1] * complex<double>(xmu[j], 0.) + a[1] * complex<double>( 2.*std::pow(xmu[j],2) - 1. ,0.)) );
    }
                                                            /** recover actual mie coefficients **/
    a[0] = complex<double>(std::pow(xx,3),0.) * a[0];
    a[1] = complex<double>(std::pow(xx,3),0.) * a[1];
    b[0] = complex<double>(std::pow(xx,3),0.) * b[0];
    b[1] = complex<double>(std::pow(xx,3),0.) * b[1];
    return;
}

KOKKOS_INLINE_FUNCTION void small2(real &xx,complex<double> &cior,bool &calcqe,int &numang,real *xmu,real *qext,real *qsca,real *gqsc,int n_,complex<double> &sforw,\
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback,\
            complex<double> *a,complex<double> *b){
    /**
    small-particle limit of mie quantities in totally reflecting
    limit ( mie series truncated after 2 terms )
    a,b       first two mie coefficients, with numerator and
              denominator expanded in powers of -xx- ( a factor
              of xx**3 is missing but is restored before return
              to calling program )  ( ref. 2, p. 1508 )
    **/
    int j;
    const real twothr = 2./3.;
    const real fivthr = 5./3.;
    real rtmp;
    complex<double>  ctmp, ciorsq;
    ciorsq = std::pow(cior,2);
    ctmp = complex<double>( 0.0, twothr ) * ( ciorsq - 1.0 );
    complex<double> tmp1 = ( std::pow(ciorsq,2)/175. - 0.275 * ciorsq + 0.25 ) * std::pow(xx,4);
    complex<double> tmp2 = std::pow(xx,3) * ctmp * ( 1.0 - 0.1 * std::pow(xx,2) );
    a[0] = ctmp * ( 1.0 - 0.1 * std::pow(xx,2) + (ciorsq/350. + 1./280.)*std::pow(xx,4))   \
             / ( ciorsq + 2.0 + ( 1.0 - 0.7 * ciorsq ) * std::pow(xx,2) - tmp1 + tmp2);
    b[0] = (std::pow(xx,2)/30.) * ctmp * ( 1.0 + (ciorsq/35. - 1./14.) *std::pow(xx,2) )   \
             / ( 1.0 - ( ciorsq/15. - 1./6. ) * std::pow(xx,2) );
    a[1] = ( 0.1 * std::pow(xx,2) ) * ctmp * ( 1.0 - std::pow(xx,2) / 14. )   \
             / ( 2. * ciorsq + 3. - ( ciorsq/7. - 0.5 ) * std::pow(xx,2) );
    qsca[n_] = 6. * std::pow(xx,4) * ( std::norm(a[0]) + std::norm(b[0]) + fivthr * std::norm(a[1]) );
    gqsc[n_] = 6. * std::pow(xx,4) * std::real( a[0] * std::conj( a[1] + b[0] ) );
    qext[n_] = qsca[n_];
    if ( calcqe ) qext[n_] = 6. * xx * std::real( a[0] + b[0] + complex<double>(fivthr,0.0) * a[1] );
    rtmp = 1.5 * std::pow(xx,3);
    sforw      = complex<double>(rtmp,0.0) * ( a[0] + b[0] + complex<double>(fivthr,0.0) * a[1] );
    sback      = complex<double>(rtmp,0.0) * ( a[0]- b[0] - complex<double>(fivthr,0.0) * a[1] );
    tforw[0]   = complex<double>(rtmp,0.0) * ( b[0] - complex<double>(fivthr,0.0) * a[1] );
    tforw[1]   = complex<double>(rtmp,0.0) * ( a[0] + 2. * complex<double>(fivthr,0.0) * a[1] );
    tback[0]   = tforw[0];
    tback[1]   = complex<double>(rtmp,0.0) * ( a[0] - 2. * complex<double>(fivthr,0.0) * a[1] );
    for(j = 0;j < numang; ++j){
        s1[j] = complex<double>(rtmp,0.0) * ( a[0] + ( b[0] + complex<double>(fivthr, 0.) * a[1] ) * complex<double>(xmu[j], 0.) );
        s2[j] = complex<double>(rtmp,0.0) * ( b[0] + a[0] * complex<double>(xmu[j], 0.) + complex<double>(fivthr, 0.) * a[1]   \
                            * complex<double>( 2.*std::pow(xmu[j],2) - 1. ,0.) );
    }
                                                /** recover actual mie coefficients **/
    a[0] = complex<double>(std::pow(xx,3),0.) * a[0];
    a[1] = complex<double>(std::pow(xx,3),0.) * a[1];
    b[0] = complex<double>(std::pow(xx,3),0.) * b[0];
    b[1] = complex<double>(0.,0.);
    return;
}

KOKKOS_INLINE_FUNCTION bool wrong(real calc,real exact,real accur){
    return std::abs( (calc - exact) / exact ) > accur;
}

KOKKOS_INLINE_FUNCTION void testmi(real *qext,real *qsca,real *gqsc,int n_,complex<double> &sforw,  \
            complex<double> &sback,complex<double> *s1,complex<double> *s2,complex<double> *tforw,complex<double> *tback, \
            real (*pmom)[8],int &momdim,bool &ok){
    /**
        compare mie code test case results with correct answers
        and return  ok=false  if even one result is inaccurate.

        the test case is :  mie size parameter = 10
                            refractive index   = 1.5 - 0.1 i
                            scattering angle = 140 degrees
                            1 sekera moment
        results for this case may be found among the test cases
        at the end of reference (1).
        *** note *** when running on some computers, esp. in single
        precision, the 'accur' criterion below may have to be relaxed.
        however, if 'accur' must be set larger than 10**-3 for some
        size parameters, your computer is probably not accurate
        enough to do mie computations for those size parameters.
    **/
   int m,n;
   real accur = 1.e-4;
   real testqe = 2.459791,testqs = 1.235144,testgq = 1.139235;
   real testpm[2] = {227.1975, 183.6898};
   complex<double> testsf = complex<double>(61.49476, -3.177994), \
                   testsb = complex<double>(1.493434, 0.2963657), \
                   tests1 = complex<double>(-0.1548380, -1.128972),\
                   tests2 = complex<double>(0.05669755, 0.5425681);
   complex<double>testtf[2] ={complex<double>(12.95238, -136.6436),complex<double>(48.54238, 133.4656)};
   complex<double>testtb[2] ={complex<double>(41.88414, -15.57833),complex<double>(43.37758, -15.28196)};
   //wrong( calc, exact ) = dabs( (calc - exact) / exact ) .gt. accur
   ok = true;
   if(wrong(qext[n_],testqe,accur)) tstbad( "qext", std::abs((qext[n_] - testqe) / testqe), ok );
   if(wrong(qsca[n_],testqs,accur)) tstbad( "qsca", std::abs((qsca[n_] - testqs) / testqs), ok );
   if(wrong(gqsc[n_],testgq,accur)) tstbad( "gqsc", std::abs((gqsc[n_] - testgq) / testgq), ok );
   if(wrong(std::real(sforw),std::real(testsf),accur) \
        || wrong(std::imag(sforw),std::imag(testsf),accur)) tstbad( "sforw", std::abs((sforw - testsf) / testsf), ok );
   if(wrong(std::real(sback),std::real(testsb),accur) \
        || wrong(std::imag(sback),std::imag(testsb),accur)) tstbad( "sback", std::abs((sback - testsb) / testsb), ok );
   if(wrong(std::real(s1[0]),std::real(tests1),accur) \
        || wrong(std::imag(s1[0]),std::imag(tests1),accur)) tstbad( "s1", std::abs((s1[0] - tests1) / tests1), ok );
   if(wrong(std::real(s2[0]),std::real(tests2),accur) \
        || wrong(std::imag(s2[0]),std::imag(tests2),accur)) tstbad( "s2", std::abs((s2[0] - tests2) / tests2), ok );
   for(n = 0; n < 2; ++n){
      if(wrong(std::real(tforw[n]),std::real(testtf[n]),accur) \
        || wrong(std::imag(tforw[n]),std::imag(testtf[n]),accur)) tstbad( "tforw", std::abs( (tforw[n] - testtf[n]) / testtf[n] ), ok );
      if(wrong(std::real(tback[n]),std::real(testtb[n]),accur) \
        || wrong(std::imag(tback[n]),std::imag(testtb[n]),accur)) tstbad( "tback", std::abs( (tback[n] - testtb[n]) / testtb[n] ), ok );
   }
   for(m = 0; m <= 1;++m){
      if(wrong(pmom[0][m],testpm[m],accur)) tstbad( "pmom", std::abs( (pmom[m][0]-testpm[m]) / testpm[m] ), ok );
   }
   return;
}

KOKKOS_INLINE_FUNCTION void tstbad(string varnam,real relerr,bool &ok){
    /**
        write name (-varnam-) of variable failing self-test and its
        percent error from the correct value.  return  ok = false.
    **/
   ok = false;
   /**
    write( *, '(/,3a,1p,e11.2,a)' )   &
    ' output variable  ', varnam,'  differed by', 100.*relerr,   &
    **/
}

KOKKOS_INLINE_FUNCTION void biga(complex<double> &cior,real &xx,int &ntrm,bool &noabs,bool &yesang,real *rbiga,complex<double> *cbiga){
/**
     calculate logarithmic derivatives of j-bessel-function
     input :  cior, xx, ntrm, noabs, yesang  (defined in 'miev0')
     output :  rbiga or cbiga  (defined in 'miev0')
     internal variables :
     confra     value of lentz continued fraction for -cbiga(ntrm)-,
                used to initialize downward recurrence.
     down       = true, use down-recurrence.  false, do not.
     f1,f2,f3   arithmetic statement functions used in determining
                whether to use up-  or down-recurrence
                ( ref. 2, eqs. 6-8 )
     mre        real refractive index
     mim        imaginary refractive index
     rezinv     1 / ( mre * xx ); temporary variable for recurrence
     zinv       1 / ( cior * xx ); temporary variable for recurrence
**/
    bool down;
    int n;
    real mre,mim,rezinv,rtmp;
    complex<double> zinv,ctmp;
                                    /**
                                     decide whether 'biga' can be
                                     calculated by up-recurrence
                                    **/
    mre = std::real(cior);
    mim = std::abs(std::imag(cior));
    if(mre < 1. || mre > 10. || mim > 10.){
        down = true;
    }else if(yesang){
        down = true;
        if(mim * xx < f2(mre)) down = false;
    }else{
        down = true;
        if(mim * xx < f1(mre)) down = false;
    }
    zinv = 1.0 / (cior * complex<double>(xx,0.));
    rezinv = 1.0 / (mre * xx);
    if(down){
                                    /**
                                     compute initial high-order 'biga' using
                                     lentz method ( ref. 1, pp. 17-20 )
                                    **/
        ctmp = confra(ntrm,zinv,xx);
                                    /**
                                     downward recurrence for 'biga'
                                     ( ref. 1, eq. 22 )
                                    **/
        if(noabs){
                                    /** no-absorption case **/
            rbiga[ntrm - 1] = std::real(ctmp);
            for(n = ntrm - 1; n >= 1;n--){
                rbiga[n - 1] = ((n+1) * rezinv) - 1. / (((n+1) * rezinv) + rbiga[n]);
            }
        }else{
                                    /** absorptive case **/
            cbiga[ntrm - 1] = ctmp;
            for(n = ntrm - 1; n >= 1;n--){
                cbiga[n - 1] = (complex<double>(n+1,0.) * zinv) - 1. / ((complex<double>(n+1,0.) * zinv) + cbiga[n]);
            }
        }
    }else{
                                    /**
                                     upward recurrence for 'biga'
                                     ( ref. 1, eqs. 20-21 )
                                    **/
        if(noabs){
                                    /** no-absorption case **/
            rtmp = std::sin(mre * xx);
            rbiga[0] = -rezinv + rtmp / (rtmp*rezinv - std::cos( mre*xx ));
            for(n = 1; n < ntrm; ++n){
                rbiga[n] = -((n+1) * rezinv) + 1. / (((n+1) * rezinv) - rbiga[n - 1]);
            }
        }else{
                                    /** absorptive case **/
            ctmp = std::exp(-complex<double>(0.,2.) * cior * complex<double>(xx,0.));
            cbiga[0] = -zinv + (1. - ctmp) / ((1. - ctmp)*zinv - complex<double>(0.,1.) * (1. + ctmp));
            for(n = 1; n < ntrm; ++n){
                cbiga[n] = -(complex<double>(n+1,0.) * zinv) + 1. / ((complex<double>(n+1,0.) * zinv) - cbiga[n - 1]);
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION real f1(real mre){
    return - 8.0 + std::pow(mre,2) * ( 26.22 + mre * ( - 0.4474 + std::pow(mre,3) * ( 0.00204 - 0.000175 * mre ) ) );
}
KOKKOS_INLINE_FUNCTION real f2(real mre){
    return 3.9 + mre * ( - 10.8 + 13.78 * mre );
}

KOKKOS_INLINE_FUNCTION complex<double> confra(int n,complex<double> &zinv,real &xx){
/**
            compute bessel function ratio capital-a-sub-n from its
            continued fraction using lentz method ( ref. 1, pp. 17-20 )
            zinv = reciprocal of argument of capital-a
      i n t e r n a l    v a r i a b l e s
      ------------------------------------
      cak      term in continued fraction expansion of capital-a
                ( ref. 1, eq. 25 )
      capt     factor used in lentz iteration for capital-a
                ( ref. 1, eq. 27 )
      cdenom   denominator in -capt-  ( ref. 1, eq. 28b )
      cnumer   numerator   in -capt-  ( ref. 1, eq. 28a )
      cdtd     product of two successive denominators of -capt-
                factors  ( ref. 1, eq. 34c )
      cntn     product of two successive numerators of -capt-
                factors  ( ref. 1, eq. 34b )
      eps1     ill-conditioning criterion
      eps2     convergence criterion
      kk       subscript k of -cak-  ( ref. 1, eq. 25b )
      kount    iteration counter ( used only to prevent runaway )
      maxit    max. allowed no. of iterations
      mm       + 1  and - 1, alternately
**/
    int mm,kk,kount;
    static int maxit = 10000;
    static real eps1 = 1.e-2,eps2 = 1.e-8;
    complex<double> cak,capt,cdenom,cdtd,cnumer,cntn;
                                            /** ref. 1, eqs. 25a, 27 **/
    complex<double> confra_ = complex<double>(( n + 1 ),0.) * zinv;
    mm     =  - 1;
    kk     = 2 * n + 3;
    cak    = complex<double>(( mm * kk ),0.) * zinv;
    cdenom = cak;
    cnumer = cdenom + 1.0 / confra_;
    kount  = 1;

  do {
    kount += 1;
    if (kount > maxit) {
      errmsg("confra--iteration failed to converge$", true);
    }
    /** ref. 2, eq. 25b **/
    mm = -mm;
    kk = kk + 2;
    cak = complex<double>((mm * kk), 0.) * zinv;
    /** ref. 2, eq. 32 **/
    if (std::abs(cnumer / cak) <= eps1 || std::abs(cdenom / cak) <= eps1) {
      /**
     ill-conditioned case -- stride
     two terms instead of one
    **/
      /** ref. 2, eqs. 34 **/
      cntn = cak * cnumer + 1.0;
      cdtd = cak * cdenom + 1.0;
      confra_ = (cntn / cdtd) * confra_;
      /** ref. 2, eq. 25b **/
      mm = -mm;
      kk = kk + 2;
      cak = complex<double>((mm * kk), 0.) * zinv;
      /** ref. 2, eqs. 35 **/
      cnumer = cak + cnumer / cntn;
      cdenom = cak + cdenom / cdtd;
      kount = kount + 1;
      //   goto label20;
    } else {
      /** well-conditioned case **/
      /** ref. 2, eqs. 26, 27 **/
      capt = cnumer / cdenom;
      confra_ = capt * confra_;
      /** check for convergence **/
      /** ( ref. 2, eq. 31 ) **/
      if (std::abs(std::real(capt) - 1.0) >= eps2 ||
          std::abs(std::imag(capt)) >= eps2) {
        /** ref. 2, eqs. 30a-b **/
        cnumer = cak + 1.0 / cnumer;
        cdenom = cak + 1.0 / cdenom;
        // goto label20;
      }
    }
  } while (std::abs(cnumer / cak) <= eps1 || std::abs(cdenom / cak) <= eps1 ||
           std::abs(std::real(capt) - 1.0) >= eps2 ||
           std::abs(std::imag(capt)) >= eps2);
    return confra_;
}

KOKKOS_INLINE_FUNCTION void lpco1t(int &nmom,int &ipolzn,int &momdim,bool *calcmo,complex<double> *a,complex<double> *b,real (*pmom)[8]){
/**
  calculate legendre polynomial expansion coefficients (also
  called moments) for phase quantities in special case where
  no. terms in mie series = 1
  input:  nmom, ipolzn, momdim     'miev0' arguments
          calcmo                   flags calculated from -ipolzn-
          a(1), b(1)               mie series coefficients
  output: pmom                     legendre moments
**/
    int nummom,l;
    real a1sq,b1sq;
    complex<double> ctmp,a1b1c;
    a1sq = std::norm(a[0]);
    b1sq = std::norm(b[0]);
    a1b1c = a[0] * std::conj( b[0] );
    if(ipolzn < 0){
        if( calcmo[0] )  pmom[0][0] = 2.25 * b1sq;
        if( calcmo[1] )  pmom[1][0] = 2.25 * a1sq;
        if( calcmo[2] )  pmom[2][0] = 2.25 * std::real( a1b1c );
        if( calcmo[3] )  pmom[3][0] = 2.25 *std::imag( a1b1c );
    }else{
        nummom = std::min( nmom, 2 );
                                    /** loop over moments **/
        for(l = 0; l <= nummom; ++l){
            if(ipolzn == 0){
                if(l == 0) pmom[0][l] = 1.5 * ( a1sq + b1sq );
                if(l == 1) pmom[0][l] = 1.5 * std::real( a1b1c );
                if(l == 2) pmom[0][l] = 0.15 * ( a1sq + b1sq );
		continue;
            }
            if(calcmo[0]){
                if(l == 0) pmom[0][l] = 2.25 * ( a1sq + b1sq / 3. );
                if(l == 1) pmom[0][l] = 1.5 * std::real( a1b1c );
                if(l == 2) pmom[0][l] = 0.3 * b1sq;
            }
            if(calcmo[1]){
                if(l == 0) pmom[1][l] = 2.25 * ( b1sq + a1sq / 3. );
                if(l == 1) pmom[1][l] = 1.5 * std::real( a1b1c );
                if(l == 2) pmom[1][l] = 0.3 * a1sq;
            }
            if(calcmo[2]){
                if(l == 0) pmom[2][l] = 3.0 * std::real( a1b1c );
                if(l == 1) pmom[2][l] = 0.75 * ( a1sq + b1sq );
                if(l == 2) pmom[2][l] = 0.3 * std::real( a1b1c );
            }
            if(calcmo[3]){
                if(l == 0) pmom[3][l] = - 1.5 * std::imag( a1b1c );
                if(l == 1) pmom[3][l] = 0.;
                if(l == 2) pmom[3][l] = 0.3 * std::imag( a1b1c );
            }
        }

    }
    return;
}

KOKKOS_INLINE_FUNCTION void lpco2t(int &nmom,int &ipolzn,int &momdim,bool *calcmo,complex<double> *a,complex<double> *b,real (*pmom)[8]){
/**
  calculate legendre polynomial expansion coefficients (also
  called moments) for phase quantities in special case where
  no. terms in mie series = 2
  input:  nmom, ipolzn, momdim     'miev0' arguments
          calcmo                   flags calculated from -ipolzn-
          a(1-2), b(1-2)           mie series coefficients
  output: pmom                     legendre moments
**/
    int l,nummom;
    real pm1,pm2,a2sq,b2sq;
    complex<double> a2c, b2c, ctmp, ca, cac, cat, cb, cbc, cbt, cg, ch;
    ca = 3. * a[0] - 5. * b[1];
    cat= 3. * b[0] - 5. * a[1];
    cac = std::conj( ca );
    a2sq = std::norm(a[1]);
    b2sq = std::norm(b[1]);
    a2c = std::conj( a[1] );
    b2c = std::conj( b[1] );
    if(ipolzn < 0){
                                    /** loop over sekera moments **/
        nummom = std::min(nmom,2);
        for(l = 0; l <= nummom; ++l){
            if(calcmo[0]){
                if( l == 0 ) pmom[0][l] = 0.25 * ( std::norm(cat) +   \
                                                    (100./3.) * b2sq );
                if( l == 1 ) pmom[0][l] = (5./3.) * std::real( cat * b2c );
                if( l == 2 ) pmom[0][l] = (10./3.) * b2sq;
            }
            if(calcmo[1]){
                if( l == 0 ) pmom[1][l] = 0.25 * ( std::norm(ca) +   \
                                                    (100./3.) * a2sq );
                if( l == 1 ) pmom[1][l] = (5./3.) * std::real( ca * a2c );
                if( l == 2 ) pmom[1][l] = (10./3.) * a2sq;
            }
            if(calcmo[2]){
                if( l == 0 ) pmom[2][l] = 0.25 * std::real(cat * cac + \
                                                    (100. / 3.) * b[1] * a2c);
                if( l == 1 ) pmom[2][l] = (5./6.) * std::real( b[1] *cac + \
                                                    cat * a2c );
                if( l == 2 ) pmom[2][l] = (10./3.) * std::real(b[1] * a2c);
            }
            if(calcmo[3]){
                if( l == 0 ) pmom[3][l] = -0.25 * std::imag(cat * cac + \
                                                    (100. / 3.) * b[1] * a2c);
                if( l == 1 ) pmom[3][l] = -(5./6.) * std::imag( b[1] *cac + \
                                                    cat * a2c );
                if( l == 2 ) pmom[3][l] = -(10./3.) * std::imag(b[1] * a2c);
            }
        }
    }else{
        cb = 3. * b[0] + 5. * a[1];
        cbt= 3. * a[0] + 5. * b[1];
        cbc = std::conj( cb );
        cg = ( cbc*cbt + 10.*( cac*a[1] + b2c*cat) ) / 3.;
        ch = 2.*( cbc*a[1] + b2c*cbt );
                                             /** loop over mueller moments **/
        nummom = std::min( nmom, 4 );
        for(l = 0; l <= nummom; ++l){
            if(ipolzn == 0 || calcmo[0]){
                if( l == 0 ) pm1 = 0.25 * std::norm(ca) + std::norm(cb) / 12.   \
                                    + (5./3.) * std::real(ca*b2c) + 5.*b2sq;
                if( l == 1 ) pm1 = std::real( cb * ( cac/6. + b2c ) );
                if( l == 2 ) pm1 = std::norm(cb)/30. + (20./7.) * b2sq   \
                                    + (2./3.) * std::real( ca * b2c );
                if( l == 3 ) pm1 = (2./7.) * std::real( cb * b2c );
                if( l == 4 ) pm1 = (40./63.) * b2sq;
                if ( calcmo[0] )  pmom[0][l] = pm1;
            }
            if(ipolzn == 0 || calcmo[1]){
                if( l == 0 ) pm2 = 0.25 * std::norm(cat) + std::norm(cbt) / 12.   \
                                    + (5./3.) * std::real(cat*a2c) + 5.*a2sq;
                if( l == 1 ) pm2 = std::real( cbt * ( std::conj(cat)/6. + a2c ) );
                if( l == 2 ) pm2 = std::norm(cbt)/30. + (20./7.) * a2sq   \
                                    + (2./3.) * std::real( cat * a2c );
                if( l == 3 ) pm2 = (2./7.) * std::real( cbt * a2c );
                if( l == 4 ) pm2 = (40./63.) * a2sq;
                if ( calcmo[1] )  pmom[1][l] = pm2;
            }
            if(ipolzn == 0){
                pmom[0][l] = 0.5 * ( pm1 + pm2 );
                continue;
            }
            if( calcmo[2] ){
               if( l == 0 ) pmom[2][l] = 0.25 * std::real( cac*cat + cg +   \
                                                       20.*b2c*a[1] );
               if( l == 1 ) pmom[2][l] = std::real( cac*cbt + cbc*cat +   \
                                                3.*ch ) / 12.;
               if( l == 2 ) pmom[2][l] = 0.1 * std::real( cg + (200./7.) *   \
                                                      b2c * a[1] );
               if( l == 3 ) pmom[2][l] = std::real( ch ) / 14.;
               if( l == 4 ) pmom[2][l] = 40./63. * std::real( b2c * a[1] );
            }
            if( calcmo[3] ){
               if( l == 0 ) pmom[3][l] = 0.25 * std::imag( cac*cat + cg +   \
                                                       20.*b2c*a[1] );
               if( l == 1 ) pmom[3][l] = std::imag( cac*cbt + cbc*cat +   \
                                                3.*ch ) / 12.;
               if( l == 2 ) pmom[3][l] = 0.1 * std::imag( cg + (200./7.) *   \
                                                      b2c * a[1] );
               if( l == 3 ) pmom[3][l] = std::imag( ch ) / 14.;
               if( l == 4 ) pmom[3][l] = 40./63. * std::imag( b2c * a[1] );
            }
        }
    }
    return;
}

KOKKOS_INLINE_FUNCTION void lpcoef(int &ntrm,int &nmom,int &ipolzn,int &momdim,bool *calcmo,int &npquan,complex<double> *a,complex<double> *b, \
            real (*pmom)[8]){
/**
        calculate legendre polynomial expansion coefficients (also
        called moments) for phase quantities ( ref. 5 formulation )
    input:  ntrm                    number terms in mie series
            nmom, ipolzn, momdim    'miev0' arguments
            calcmo                  flags calculated from -ipolzn-
            npquan                  defined in 'miev0'
            a, b                    mie series coefficients
    output: pmom                   legendre moments ('miev0' argument)
    *** notes ***

        (1)  eqs. 2-5 are in error in dave, appl. opt. 9,
        1888 (1970).  eq. 2 refers to m1, not m2;  eq. 3 refers to
        m2, not m1.  in eqs. 4 and 5, the subscripts on the second
        term in square brackets should be interchanged.

        (2)  the general-case logic in this subroutine works correctly
        in the two-term mie series case, but subroutine  'lpco2t'
        is called instead, for speed.

        (3)  subroutine  'lpco1t', to do the one-term case, is never
        called within the context of 'miev0', but is included for
        complete generality.

        (4)  some improvement in speed is obtainable by combining the
        310- and 410-loops, if moments for both the third and fourth
        phase quantities are desired, because the third phase quantity
        is the real part of a complex series, while the fourth phase
        quantity is the imaginary part of that very same series.  but
        most users are not interested in the fourth phase quantity,
        which is related to circular polarization, so the present
        scheme is usually more efficient.
        ** specification of local variables
    am(m)       numerical coefficients  a-sub-m-super-l
                    in dave, eqs. 1-15, as simplified in ref. 5.
    bi(i)       numerical coefficients  b-sub-i-super-l
                    in dave, eqs. 1-15, as simplified in ref. 5.
    bidel(i)    1/2 bi(i) times factor capital-del in dave
    cm,      ()     arrays c and d in dave, eqs. 16-17 (mueller form),
                    calculated using recurrence derived in ref. 5
    cs,ds()     arrays c and d in ref. 4, eqs. a5-a6 (sekera form),
                    calculated using recurrence derived in ref. 5
    c,d()       either -cm, - or -cs,ds-, depending on -ipolzn-
    evenl       true for even-numbered moments;  false otherwise
    idel        1 + little-del  in dave
    maxtrm      max. no. of terms in mie series
    maxmom      max. no. of non-zero moments
    nummom      number of non-zero moments
    recip(k)    1 / k
**/
    const int maxtrm = 1502;
    const int maxmom = 2 * maxtrm;
    const int mxmom2 = maxmom / 2;
    const int maxrcp = 4 * maxtrm + 2;
    real am[maxtrm + 1],bi[mxmom2 + 1],bidel[mxmom2 + 1];
    static real recip[maxrcp + 1];
    complex<double> cs[maxtrm],ds[maxtrm];
    complex<double> c[maxtrm];
    complex<double> cm[maxtrm];
    complex<double> d[maxtrm];
    complex<double> dm[maxtrm];
    int k,j,l,nummom,ld2,idel,m,i,mmax,imax;

    bool evenl;
    static bool pass1 = true;

    if(pass1){
        for(k = 1; k < maxrcp + 1; ++k){
            recip[k] = 1. / k;
        }
        pass1 = false;
    }
    for (j = 0; j < std::max(1,npquan); ++j){
        for(l = 0; l <= nmom; ++l){
            pmom[j][l] = 0.;
        }
    }
    if(ntrm == 1){
        lpco1t ( nmom, ipolzn, momdim, calcmo, a, b, pmom );
	return;
    }else if(ntrm == 2){
        lpco2t ( nmom, ipolzn, momdim, calcmo, a, b, pmom );
	return;
    }
    if(ntrm + 2 > maxtrm){

        errmsg( "lpcoef--parameter maxtrm too small", true );
    }
                                                /** calculate mueller c, d arrays **/
    cm[ ntrm+1 ] = complex<double>( 0., 0. );
    dm[ ntrm+1 ] = complex<double>( 0., 0. );
    cm[ ntrm ] = complex<double>( 1. - recip[ ntrm + 1] ,0) * b[ntrm - 1];
    dm[ ntrm ] = complex<double>( 1. - recip[ ntrm + 1 ] ,0) * a[ntrm - 1];
    cm[ ntrm - 1 ] = complex<double>( recip[ ntrm ] + recip[ntrm + 1] ,0.) * a[ntrm - 1]   \
                + complex<double>( 1. - recip[ ntrm] ,0) * b[ ntrm-2 ];
    dm[ ntrm - 1 ] = complex<double>( recip[ ntrm ] + recip[ntrm + 1] ,0.) * b[ntrm - 1]   \
                + complex<double>( 1. - recip[ ntrm] ,0) * a[ ntrm-2 ];
    for(k = ntrm - 2; k >= 1; --k){
         cm[k] = cm[k + 2] - complex<double>( 1. + recip[k + 2] ,0) * b[k + 1]   \
                             + complex<double>( recip[k + 1] + recip[k + 2] ,0) * a[k]   \
                             + complex<double>( 1. - recip[k + 1] ,0) * b[k - 1];
         dm[k] = dm[k + 2] - complex<double>( 1. + recip[k + 2] ,0) * a[k + 1]   \
                             + complex<double>( recip[k + 1] + recip[k + 2] ,0) * b[k]   \
                             + complex<double>( 1. - recip[k + 1] ,0) * a[k - 1];
    }
    cm[0] = cm[2] + 1.5 * ( a[0] - b[1] );
    dm[0] = dm[2] + 1.5 * ( b[0] - a[1] );
    if(ipolzn >= 0){
        for(k = 0; k < ntrm + 2; ++k){
            c[k] = complex<double>( 2* (k + 1) - 1  ,0) * cm[k];
            d[k] = complex<double>( 2 * (k + 1) - 1 ,0) * dm[k];
        }
    }else{
                                                /** compute sekera c and d arrays **/
         cs[ntrm + 1] = complex<double>( 0., 0. );
         ds[ntrm + 1] = complex<double>( 0., 0. );
         cs[ntrm]     = complex<double>( 0., 0. );
         ds[ntrm]     = complex<double>( 0., 0. );
         for(k = ntrm - 1; k >= 0; --k){
            cs[k] = cs[k + 2] + complex<double>( 2*(k+1) + 1 ,0) * ( cm[k + 1] - b[k] );
            ds[k] = ds[k + 2] + complex<double>( 2*(k+1) + 1 ,0) * ( dm[k + 1] - a[k] );
         }
         for(k = 0; k < ntrm + 2; ++k){
            c[k] = complex<double>( 2*(k+1) -1 ,0) * cs[k];
            d[k] = complex<double>( 2*(k+1) - 1 ,0) * ds[k];
         }
    }
    if(ipolzn < 0) {
        nummom = std::min( nmom, 2*ntrm - 2 );
    }
    if(ipolzn >= 0) {
        nummom = std::min( nmom, 2*ntrm );
    }

    if(nummom > maxmom){
        // while(true){}
        // errmsg( "lpcoef--parameter maxtrm too small", true );
    }
    /** loop over moments **/
        //parallel_for("lpcoef",Bounds<1>(LBnd(0,nummom)),YAKL_LAMBDA(int l){
    for(l = 0; l <= nummom; ++l){
        ld2 = l / 2;
        evenl = (l % 2 == 0);
        /** calculate numerical coefficients
            a-sub-m and b-sub-i in dave
            double-sums for moments
        **/
        if(l == 0){
            idel = 1;
            for(m = 0;m <= ntrm;++m){
                am[m] = 2.0 * recip[2 * m + 1] ;
            }
            bi[0] = 1.0;
        }else if(evenl){
            idel = 1;
            for(m = ld2;m <= ntrm; ++m){
                am[m] = ( 1. + recip[2*m-l + 1] ) * am[m];
            }
            for(i = 0; i <= ld2 - 1; ++i){
                bi[i] = ( 1. - recip[ l-2*i] ) * bi[i];
            }
            bi[ld2] = ( 2. - recip[l] ) * bi[ld2 - 1];
        }else{
            idel = 2;
            for(m = ld2;m <= ntrm; ++m){
                am[m] = ( 1. - recip[2*m + l + 2] ) * am[m];
            }
            for(i = 0; i <= ld2; ++i){
                bi[i] = ( 1. - recip[ l + 2*i +   1] ) * bi[i];
            }
        }
                                                /** establish upper limits for sums
                                                    and incorporate factor capital-
                                                    del into b-sub-i
                                                **/
        mmax = ntrm - idel;
        if( ipolzn >= 0 )  mmax = mmax + 1;
        imax = std::min( ld2, mmax - ld2 );
        if( imax < 0 )  return;
        for(i = 0; i <= imax; ++i){
            bidel[i] = bi [i];
        }

        if( evenl )  bidel[0] = 0.5 * bidel[0];
                                                /** perform double sums just for
                                                    phase quantities desired by user
                                                **/
        if(ipolzn == 0){
            for(i = 0; i <= imax; ++i){
        // //                                         /** vectorizable loop (cray) **/
                real thesum = 0.0;
                for(m = ld2; m <= mmax-i; ++m){
                    thesum = thesum + am[m] * ( std::real( c[m - i] * std::conj( c[m + i + idel - 1] ) ) \
                                                +  std::real( d[m - i]  * std::conj( d[m + i + idel - 1] ) ) );
                }
                pmom[0][l] += bidel[i] * thesum;
             }
            pmom[0][l] = 0.5 * pmom[0][l];
            continue;
        }
        if(calcmo[0]){
            for(i = 0; i<= imax; ++i){
                                                /** vectorizable loop (cray) **/
                real thesum = 0.0;
                for(m = ld2; m <= mmax - i; ++m){
                    thesum = thesum + am[m]*   \
                              std::real( c[m - i] * std::conj( c[m + i + idel - 1] ) );
                }
                pmom[0][l] += bidel[i] * thesum;
            }
        }
        if(calcmo[1]){
            for(i = 0; i<= imax; ++i){
                                                /** vectorizable loop (cray) **/
                real thesum = 0.0;
                for(m = ld2; m <= mmax - i; ++m){
                    thesum = thesum + am[m]*   \
                              std::real( d[m - i] * std::conj( d[m + i + idel - 1] ) );
                }
                pmom[1][l] += bidel[i] * thesum;
            }
        }
        if(calcmo[2]){
            for(i = 0; i<= imax; ++i){
                                                /** vectorizable loop (cray) **/
                real thesum = 0.0;
                for(m = ld2; m <= mmax - i; ++m){
                    thesum = thesum + am[m]*   \
                              (std::real( c[m - i] * std::conj( d[m + i + idel - 1] )) +\
                              std::real(  c[m + i + idel - 1] * std::conj( d[m - i] )));
                }
                pmom[2][l] += bidel[i] * thesum;
            }
            pmom[2][l] = 0.5 * pmom[2][l];
        }
        if(calcmo[3]){
            for(i = 0; i<= imax; ++i){
                                                /** vectorizable loop (cray) **/
                real thesum = 0.0;
                for(m = ld2; m <= mmax - i; ++m){
                    thesum = thesum + am[m]*   \
                              (std::imag( c[m - i] * std::conj( d[m + i + idel - 1] )) + \
                              std::imag(c[m + i + idel - 1] * std::conj( d[m - i] )));
                }
                pmom[3][l] += bidel[i] * thesum;
            }
            pmom[3][l] = -0.5 * pmom[3][l];
        }
    }
    return;
}

KOKKOS_INLINE_FUNCTION void miprnt(bool *prnt,real &xx,bool &perfct,complex<double> &crefin,int &numang,real *xmu,real *qext,real *qsca,real *gqsc,int n_,\
            int &nmom,int &ipolzn,int &momdim,bool *calcmo,real (*pmom)[8],complex<double> &sforw,complex<double> &sback,\
            complex<double> *tforw,complex<double> *tback,complex<double> *s1,complex<double> *s2){
    int i,m,j;
    real fi1,fi2,fnorm;
    char fmt[22];
    /**
    if ( perfct )  write ( *, '(''1'',10x,a,1p,e11.4)' )   &
        'perfectly conducting case, size parameter =', xx
    if ( .not.perfct )  write ( *, '(''1'',10x,3(a,1p,e11.4))' )   &
                    'refractive index:  real ', dble(crefin),   &
                    '  imag ', dimag(crefin), ',   size parameter =', xx
     **/
    if(prnt[0] && numang > 0){
    /**
    write ( *, '(/,a)' )   &
        '    cos(angle)  ------- s1 ---------  ------- s2 ---------'//   &
        '  --- s1*conjg(s2) ---   i1=s1**2   i2=s2**2  (i1+i2)/2'//   &
        '  deg polzn'
     **/
        for(i = 0; i < numang; ++i){
            fi1 = std::pow(std::real( s1[i] ),2) + std::pow(std::imag( s1[i] ),2);
            fi2 = std::pow(std::real( s2[i] ),2) + std::pow(std::imag( s2[i] ),2);
            /**
            write( *, '( i4, f10.6, 1p,10e11.3 )'   )   &
                i, xmu(i), s1(i), s2(i), s1(i)*dconjg(s2(i)),   &
                fi1, fi2, 0.5*(fi1+fi2), (fi2-fi1)/(fi2+fi1)
             **/
        }

    }
    if(prnt[1]){
        /**
        write ( *, '(/,a,9x,a,17x,a,17x,a,/,(0p,f7.2, 1p,6e12.3) )' )   &
            '  angle', 's-sub-1', 't-sub-1', 't-sub-2',   &
            0.0,     sforw,    tforw(1),  tforw(2),   &
            180.,     sback,    tback(1),  tback(2)
        write ( *, '(/,4(a,1p,e11.4))' )   &
                    ' efficiency factors,  extinction:', qext,   &
                    '   scattering:', qsca,   &
                    '   absorption:', qext-qsca,   &
                    '   rad. pressure:', qext-gqsc
        **/
        if(nmom > 0.){
            /** write( *, '(/,a)' )  ' normalized moments of :'
                if ( ipolzn.eq.0 ) write ( *, '(''+'',27x,a)' ) 'phase fcn'
                if ( ipolzn.gt.0 )  write ( *, '(''+'',33x,a)' )   &
                'm1           m2          s21          d21'
                if ( ipolzn.lt.0 )  write ( *, '(''+'',33x,a)' )   &
                'r1           r2           r3           r4'
            **/
            fnorm = 4. / (std::pow(xx,2) * qsca[n_]);
            for(m = 0; m <= nmom; ++m){
                //write ( *, '(a,i4)' )  '      moment no.', m
                for(j = 0; j < 4; ++j){
                    if(calcmo[j]){
                    /**
                        write( fmt, 98 )  24 + (j-1)*13
                        write ( *,fmt )  fnorm * pmom(m,j)
                    **/
                    }
                }
            }
        }
    }
    return;
}
// clang-format on
