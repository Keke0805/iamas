#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>
#include <flcl-cxx.hpp>
#include <cstring>
#include <iostream>
#include "rrtmg_lw.hpp"
using std::string;
using real1d_u = Kokkos::Experimental::OffsetView<real*, Kokkos::LayoutLeft, flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real2d_u = Kokkos::Experimental::OffsetView<real**, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using real3d_u = Kokkos::Experimental::OffsetView<real***, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using real4d_u = Kokkos::Experimental::OffsetView<real****, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using int1d_u = Kokkos::Experimental::OffsetView<int*, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using int3d_u = Kokkos::Experimental::OffsetView<int***, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using real1d = Kokkos::Experimental::OffsetView<real*, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using real2d = Kokkos::Experimental::OffsetView<real**, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using real3d = Kokkos::Experimental::OffsetView<real***, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using real4d = Kokkos::Experimental::OffsetView<real****, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using real5d = Kokkos::Experimental::OffsetView<real*****, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using int2d = Kokkos::Experimental::OffsetView<int**, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using int5d = Kokkos::Experimental::OffsetView<int*****, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using int1d = Kokkos::Experimental::OffsetView<int*, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using int3d = Kokkos::Experimental::OffsetView<int***, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
using bool3d = Kokkos::Experimental::OffsetView<bool***, Kokkos::LayoutLeft,flcl::HostMemorySpace>;
struct rrlw lw_rad;
// struct lw_workspace ws;
static lw_workspace *lw_ws = nullptr;
static bool is_lw_ws_inited = false;
static int lw_ws_kts = 0;
static int lw_ws_kte = 0;
static int lw_ws_its = 0;
static int lw_ws_ite = 0;
static int lw_ws_jts = 0;
static int lw_ws_jte = 0;
static int lw_ws_nbndlw = 0;
static int lw_ws_ngptlw = 0;
static int lw_ws_mxmol = 0;
static int lw_ws_maxxsec = 0;
static int lw_ws_noznlevels = 0;
const int ngptlw = lw_rad.ngptlw;
const int mg     = lw_rad.mg;
const int nbndlw = lw_rad.nbndlw;
const int maxinpx= lw_rad.maxinpx;
const int mxmol= lw_rad.mxmol;
const int maxxsec= lw_rad.maxxsec;
const int nmol= lw_rad.nmol;
const int ntbl= lw_rad.ntbl;
const int SIZE_NLAY = 110;
const int SIZE_NCOL = 3;
const int SIZE_NSUBCOL = ngptlw + 2;
const int SIZE_LEVSIZ = 110;
const int RAD_LW_CHUNK_SIZE = 64;

static bool lw_workspace_matches(int kts, int kte, int its, int ite, int jts, int jte,
                                 int nbndlw_in, int ngptlw_in, int mxmol_in,
                                 int maxxsec_in, int noznlevels_in) {
    return is_lw_ws_inited && lw_ws != nullptr &&
           lw_ws_kts == kts &&
           lw_ws_kte == kte &&
           lw_ws_its == its &&
           lw_ws_ite == ite &&
           lw_ws_jts == jts &&
           lw_ws_jte == jte &&
           lw_ws_nbndlw == nbndlw_in &&
           lw_ws_ngptlw == ngptlw_in &&
           lw_ws_mxmol == mxmol_in &&
           lw_ws_maxxsec == maxxsec_in &&
           lw_ws_noznlevels == noznlevels_in;
}

static lw_workspace& get_lw_ws(int kts, int kte, int its, int ite, int jts, int jte,
                               int nbndlw_in, int ngptlw_in, int mxmol_in,
                               int maxxsec_in, int noznlevels_in) {
    if (!lw_workspace_matches(kts, kte, its, ite, jts, jte,
                              nbndlw_in, ngptlw_in, mxmol_in,
                              maxxsec_in, noznlevels_in)) {
        if (lw_ws != nullptr) {
            Kokkos::fence();
            delete lw_ws;
            lw_ws = nullptr;
        }

        lw_ws = new lw_workspace();
        lw_ws->init(kts, kte, its, ite, jts, jte,
                    nbndlw_in, ngptlw_in, mxmol_in, maxxsec_in, noznlevels_in);

        lw_ws_kts = kts;
        lw_ws_kte = kte;
        lw_ws_its = its;
        lw_ws_ite = ite;
        lw_ws_jts = jts;
        lw_ws_jte = jte;
        lw_ws_nbndlw = nbndlw_in;
        lw_ws_ngptlw = ngptlw_in;
        lw_ws_mxmol = mxmol_in;
        lw_ws_maxxsec = maxxsec_in;
        lw_ws_noznlevels = noznlevels_in;
        is_lw_ws_inited = true;
    }

    return *lw_ws;
}

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
// Tempering parameters
#define TMASKB 0x9d2c5680UL 
#define TMASKC 0xefc60000UL 
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))
static unsigned long state[N]; /* the array for the state vector  */
static int left = 1;
static int initf = 0;
static unsigned long *next;
KOKKOS_FUNCTION void inirad(const real3d &o3prof,const real3d &plev,int kts,int kte,int i,int j);
KOKKOS_FUNCTION void o3data(const real3d &o3prof,const real3d &plev,int kts,int kte,int i, int j);
KOKKOS_FUNCTION string trim(char *config_chem_opt,string &schemes);
KOKKOS_FUNCTION void vinterp_ozn(int lchnk, int ncol, int pcols, int pver, const real4d &pmid_in, const real1d_u &pin_in, int levsiz, const real3d &ozmix_in, const real3d &o3vmr_out,int m,int j);
KOKKOS_FUNCTION void relcalc(int ncol,int pcols,int pver, const real4d &t, const real3d &landfrac, const real3d &landm, const real3d &icefrac,const real4d &rel, const real3d &snowh,int m, int n);
KOKKOS_FUNCTION void reicalc(int ncol,int pcols,int pver,const real4d &t, const real4d &re,const real1d_u &retab,int m, int n);
KOKKOS_FUNCTION void mcica_subcol_lw(int& iplon,int& ncol,int& nlay,int& icld,int& permuteseed,int& irng, const real4d &play,
const real4d &cldfrac, const real4d &ciwp, const real4d &clwp, const real4d &cswp, const real4d &rei, const real4d &rel, const real4d &res, const real5d &tauc,
const real5d &cldfmcl, const real5d &ciwpmcl, const real5d &clwpmcl, const real5d &cswpmcl, const real4d &reicmcl, const real4d &relqmcl, const real4d &resnmcl, const real5d &taucmcl,
int m , int n);
KOKKOS_FUNCTION void rrtmg_lw(int ncol, int nlay,int &icld, const real4d &play,const real4d &plev,const real4d &tlay, const real4d &tlev, const real3d &tsfc, const real4d &h2ovmr, const real4d &o3vmr,
const real4d &co2vmr, const real4d &ch4vmr, const real4d &n2ovmr, const real4d &o2vmr, const real4d &cfc11vmr, const real4d &cfc12vmr, const real4d &cfc22vmr, const real4d &ccl4vmr,
const real4d &emis,int &inflglw,int &iceflglw,int &liqflglw, const real5d &cldfmcl, const real5d &taucmcl, const real5d &ciwpmcl, const real5d &clwpmcl,
const real5d &cswpmcl, const real4d &reicmcl, const real4d &relqmcl, const real4d &resnmcl, const real5d &tauaer, const real4d &uflx, const real4d &dflx, const real4d &hr,
const real4d &uflxc, const real4d &dflxc, const real4d &hrc,const real1d_u &tau_tbl,const real1d_u &tfn_tbl,const real1d_u &exp_tbl,const real1d_u &fracrefa_01,const real1d_u &fracrefb_01,const real2d_u &absa_01,const real3d_u &ka_01,const real2d_u &absb_01,const real3d_u &kb_01,const real2d_u &ka_mn2_01,const real2d_u &kb_mn2_01,const real2d_u &selfref_01,const real2d_u &forref_01,                 
const real1d_u &fracrefa_02,const real1d_u &fracrefb_02,const real2d_u &absa_02,const real3d_u &ka_02,const real2d_u &absb_02,const real3d_u &kb_02,const real2d_u &selfref_02,const real2d_u &forref_02,                                     
const real2d_u &fracrefa_03,const real2d_u &fracrefb_03,const real2d_u &absa_03,const real4d_u &ka_03,const real2d_u &absb_03,const real4d_u &kb_03,const real3d_u &ka_mn2o_03,const real3d_u &kb_mn2o_03,const real2d_u &selfref_03,const real2d_u &forref_03,               
const real2d_u &fracrefa_04,const real2d_u &fracrefb_04,const real2d_u &absa_04,const real4d_u &ka_04,const real2d_u &absb_04,const real4d_u &kb_04,const real2d_u &selfref_04,const real2d_u &forref_04,                                     
const real2d_u &fracrefa_05,const real2d_u &fracrefb_05,const real2d_u &absa_05,const real4d_u &ka_05,const real2d_u &absb_05,const real4d_u &kb_05,const real3d_u &ka_mo3_05,const real2d_u &selfref_05,const real2d_u &forref_05,const real1d_u &ccl4_05,                   
const real1d_u &fracrefa_06,const real2d_u &absa_06,const real3d_u &ka_06,const real2d_u &ka_mco2_06,const real2d_u &selfref_06,const real2d_u &forref_06,const real1d_u &cfc11adj_06,const real1d_u &cfc12_06,                           
const real2d_u &fracrefa_07,const real1d_u &fracrefb_07,const real2d_u &absa_07,const real4d_u &ka_07,const real2d_u &absb_07,const real3d_u &kb_07,const real3d_u &ka_mco2_07,const real2d_u &kb_mco2_07,const real2d_u &selfref_07,const real2d_u &forref_07,               
const real1d_u &fracrefa_08,const real1d_u &fracrefb_08,const real2d_u &absa_08,const real3d_u &ka_08,const real2d_u &absb_08,const real3d_u &kb_08,const real2d_u &ka_mco2_08,const real2d_u &ka_mn2o_08,const real2d_u &ka_mo3_08,const real2d_u &kb_mco2_08,const real2d_u &kb_mn2o_08,    
const real2d_u &selfref_08,const real2d_u &forref_08,const real1d_u &cfc12_08,const real1d_u &cfc22adj_08,const real2d_u &fracrefa_09,const real1d_u &fracrefb_09,const real2d_u &absa_09,const real4d_u &ka_09,const real2d_u &absb_09,const real3d_u &kb_09,      
const real3d_u &ka_mn2o_09,const real2d_u &kb_mn2o_09,const real2d_u &selfref_09,const real2d_u &forref_09,const real1d_u &fracrefa_10,const real1d_u &fracrefb_10,const real2d_u &absa_10,const real3d_u &ka_10,const real2d_u &absb_10,const real3d_u &kb_10,      
const real2d_u &selfref_10,const real2d_u &forref_10,const real1d_u &fracrefa_11,const real1d_u &fracrefb_11,const real2d_u &absa_11,const real3d_u &ka_11,const real2d_u &absb_11,const real3d_u &kb_11,      
const real2d_u &ka_mo2_11,const real2d_u &kb_mo2_11,const real2d_u &selfref_11,const real2d_u &forref_11,const real2d_u &fracrefa_12,const real2d_u &absa_12,const real4d_u &ka_12,const real2d_u &selfref_12,const real2d_u &forref_12,           
const real2d_u &fracrefa_13,const real1d_u &fracrefb_13,const real2d_u &absa_13,const real4d_u &ka_13,const real3d_u &ka_mco2_13,const real3d_u &ka_mco_13,const real2d_u &kb_mo3_13,const real2d_u &selfref_13,const real2d_u &forref_13,                 
const real1d_u &fracrefa_14,const real1d_u &fracrefb_14,const real2d_u &absa_14,const real3d_u &ka_14,const real2d_u &absb_14,const real3d_u &kb_14,const real2d_u &selfref_14,const real2d_u &forref_14,                                     
const real2d_u &fracrefa_15,const real2d_u &absa_15,const real4d_u &ka_15,const real3d_u &ka_mn2_15,const real2d_u &selfref_15,const real2d_u &forref_15,const real2d_u &fracrefa_16,const real1d_u &fracrefb_16,const real2d_u &absa_16,const real4d_u &ka_16,const real2d_u &absb_16,const real3d_u &kb_16,      
const real2d_u &selfref_16,const real2d_u &forref_16,const real4d &fracs,const real4d &taug,const real4d &taut,const real4d &cldfmc,const real4d &ciwpmc,const real4d &clwpmc,const real4d &cswpmc,const real4d &taucmc,
const real4d &wkl,const real4d &wx,const real4d &taua,const real4d &planklay,const real4d &planklev,const real4d &efclfrac,const real4d &odcld,int m, int n);
// const real3d &bbugas,const real3d &bbutot,const real3d &secdiff,const int3d &icldlyr,const real3d &uclfl,const real3d &dclfl,const real3d &clrurad,const real3d &clrdrad,
// const real3d &atot,const real3d &atrans,int m, int n);

/* initializes state[N] with a seed */
KOKKOS_FUNCTION
void init_genrand(unsigned long s)
{
    int j;
    state[0]= s & 0xffffffffUL;
    for (j=1; j<N; j++) {
        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array state[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        state[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    left = 1; initf = 1;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
KOKKOS_FUNCTION
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
    }

    state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    left = 1; initf = 1;
}
KOKKOS_FUNCTION
static void next_state(void)
{
    unsigned long *p=state;
    int j;

    /* if init_genrand() has not been called, */
    /* a default initial seed is used         */
    if (initf==0) init_genrand(5489UL);

    left = N;
    next = state;
    
    for (j=N-M+1; --j; p++) 
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (j=M; --j; p++) 
        *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
KOKKOS_FUNCTION
unsigned long genrand_int32(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
KOKKOS_FUNCTION
long genrand_int31(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (long)(y>>1);
}

/* generates a random number on [0,1]-real-interval */
KOKKOS_FUNCTION double genrand_real1(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
KOKKOS_FUNCTION
double genrand_real2(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
KOKKOS_FUNCTION
double genrand_real3(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y + 0.5) * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
KOKKOS_FUNCTION
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return (a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

extern "C" void rrtmg_lw_workspace_finalize() {
    if (lw_ws != nullptr) {
        if (Kokkos::is_initialized()) {
            Kokkos::fence();
        }
        delete lw_ws;
        lw_ws = nullptr;
    }

    is_lw_ws_inited = false;
    lw_ws_kts = 0;
    lw_ws_kte = 0;
    lw_ws_its = 0;
    lw_ws_ite = 0;
    lw_ws_jts = 0;
    lw_ws_jte = 0;
    lw_ws_nbndlw = 0;
    lw_ws_ngptlw = 0;
    lw_ws_mxmol = 0;
    lw_ws_maxxsec = 0;
    lw_ws_noznlevels = 0;
}

extern "C" void rrtmg_lwrad_kokkos(int &ids,int &ide,int &jds,int &jde,int &kds,int &kde,int &ims,int &ime,int &jms,int &jme,int &kms,int &kme,int &its,int &ite,
    int &jts,int &jte,int &kts,int &kte,int &icloud,int &has_reqc,int &has_reqi,int &has_reqs,int &o3input,flcl_ndarray_t *emiss_p,flcl_ndarray_t *tsk_p,
    flcl_ndarray_t *snow_p,flcl_ndarray_t *xice_p,flcl_ndarray_t *xland_p,flcl_ndarray_t *t3d_p,flcl_ndarray_t *p3d_p,flcl_ndarray_t *pi3d_p,flcl_ndarray_t *dz8w_p,flcl_ndarray_t *p8w_p,flcl_ndarray_t *t8w_p,flcl_ndarray_t *re_cloud_p,flcl_ndarray_t *re_ice_p,flcl_ndarray_t *re_snow_p,
    flcl_ndarray_t *cldfra3d_p,flcl_ndarray_t *qv3d_p,flcl_ndarray_t *qc3d_p,flcl_ndarray_t *qr3d_p,flcl_ndarray_t *qi3d_p,flcl_ndarray_t *qs3d_p,flcl_ndarray_t *qg3d_p,flcl_ndarray_t *o33d_p,int &noznlevels,flcl_ndarray_t *pin_p,flcl_ndarray_t *o3clim_p,
    char *config_chem_opt_p, bool &config_aer_ra_feedback_p,flcl_ndarray_t *tauaerlw_p,flcl_ndarray_t *glw_p,flcl_ndarray_t *olr_p,flcl_ndarray_t *lwcf_p,flcl_ndarray_t *lwupt_p,flcl_ndarray_t *lwuptc_p,
    flcl_ndarray_t *lwdnt_p,flcl_ndarray_t *lwdntc_p,flcl_ndarray_t *lwupb_p,flcl_ndarray_t *lwupbc_p,flcl_ndarray_t *lwdnb_p,flcl_ndarray_t *lwdnbc_p,flcl_ndarray_t *rthratenlw_p,flcl_ndarray_t *rre_cloud_p,flcl_ndarray_t *rre_ice_p,
    flcl_ndarray_t *rre_snow_p,flcl_ndarray_t *lwupflx_p,flcl_ndarray_t *lwupflxc_p,flcl_ndarray_t *lwdnflx_p,flcl_ndarray_t *lwdnflxc_p,flcl_ndarray_t *pprof_p,flcl_ndarray_t *tprof_p,int &nproflevs,bool &has_mpas,bool &has_chemistry,
    flcl_ndarray_t *retab_p,flcl_ndarray_t *tau_tbl_p,flcl_ndarray_t *tfn_tbl_p,flcl_ndarray_t *exp_tbl_p,flcl_ndarray_t *fracrefa_01_p,flcl_ndarray_t *fracrefb_01_p,flcl_ndarray_t *absa_01_p,flcl_ndarray_t *ka_01_p,flcl_ndarray_t *absb_01_p,flcl_ndarray_t *kb_01_p,flcl_ndarray_t *ka_mn2_01_p,flcl_ndarray_t *kb_mn2_01_p,flcl_ndarray_t *selfref_01_p,flcl_ndarray_t *forref_01_p,                 
    flcl_ndarray_t *fracrefa_02_p,flcl_ndarray_t *fracrefb_02_p,flcl_ndarray_t *absa_02_p,flcl_ndarray_t *ka_02_p,flcl_ndarray_t *absb_02_p,flcl_ndarray_t *kb_02_p,flcl_ndarray_t *selfref_02_p,flcl_ndarray_t *forref_02_p,                                     
    flcl_ndarray_t *fracrefa_03_p,flcl_ndarray_t *fracrefb_03_p,flcl_ndarray_t *absa_03_p,flcl_ndarray_t *ka_03_p,flcl_ndarray_t *absb_03_p,flcl_ndarray_t *kb_03_p,flcl_ndarray_t *ka_mn2o_03_p,flcl_ndarray_t *kb_mn2o_03_p,flcl_ndarray_t *selfref_03_p,flcl_ndarray_t *forref_03_p,               
    flcl_ndarray_t *fracrefa_04_p,flcl_ndarray_t *fracrefb_04_p,flcl_ndarray_t *absa_04_p,flcl_ndarray_t *ka_04_p,flcl_ndarray_t *absb_04_p,flcl_ndarray_t *kb_04_p,flcl_ndarray_t *selfref_04_p,flcl_ndarray_t *forref_04_p,                                     
    flcl_ndarray_t *fracrefa_05_p,flcl_ndarray_t *fracrefb_05_p,flcl_ndarray_t *absa_05_p,flcl_ndarray_t *ka_05_p,flcl_ndarray_t *absb_05_p,flcl_ndarray_t *kb_05_p,flcl_ndarray_t *ka_mo3_05_p,flcl_ndarray_t *selfref_05_p,flcl_ndarray_t *forref_05_p,flcl_ndarray_t *ccl4_05_p,                   
    flcl_ndarray_t *fracrefa_06_p,flcl_ndarray_t *absa_06_p,flcl_ndarray_t *ka_06_p,flcl_ndarray_t *ka_mco2_06_p,flcl_ndarray_t *selfref_06_p,flcl_ndarray_t *forref_06_p,flcl_ndarray_t *cfc11adj_06_p,flcl_ndarray_t *cfc12_06_p,                           
    flcl_ndarray_t *fracrefa_07_p,flcl_ndarray_t *fracrefb_07_p,flcl_ndarray_t *absa_07_p,flcl_ndarray_t *ka_07_p,flcl_ndarray_t *absb_07_p,flcl_ndarray_t *kb_07_p,flcl_ndarray_t *ka_mco2_07_p,flcl_ndarray_t *kb_mco2_07_p,flcl_ndarray_t *selfref_07_p,flcl_ndarray_t *forref_07_p,               
    flcl_ndarray_t *fracrefa_08_p,flcl_ndarray_t *fracrefb_08_p,flcl_ndarray_t *absa_08_p,flcl_ndarray_t *ka_08_p,flcl_ndarray_t *absb_08_p,flcl_ndarray_t *kb_08_p,flcl_ndarray_t *ka_mco2_08_p,flcl_ndarray_t *ka_mn2o_08_p,flcl_ndarray_t *ka_mo3_08_p,flcl_ndarray_t *kb_mco2_08_p,flcl_ndarray_t *kb_mn2o_08_p,    
    flcl_ndarray_t *selfref_08_p,flcl_ndarray_t *forref_08_p,flcl_ndarray_t *cfc12_08_p,flcl_ndarray_t *cfc22adj_08_p,flcl_ndarray_t *fracrefa_09_p,flcl_ndarray_t *fracrefb_09_p,flcl_ndarray_t *absa_09_p,flcl_ndarray_t *ka_09_p,flcl_ndarray_t *absb_09_p,flcl_ndarray_t *kb_09_p,      
    flcl_ndarray_t *ka_mn2o_09_p,flcl_ndarray_t *kb_mn2o_09_p,flcl_ndarray_t *selfref_09_p,flcl_ndarray_t *forref_09_p,flcl_ndarray_t *fracrefa_10_p,flcl_ndarray_t *fracrefb_10_p,flcl_ndarray_t *absa_10_p,flcl_ndarray_t *ka_10_p,flcl_ndarray_t *absb_10_p,flcl_ndarray_t *kb_10_p,      
    flcl_ndarray_t *selfref_10_p,flcl_ndarray_t *forref_10_p,flcl_ndarray_t *fracrefa_11_p,flcl_ndarray_t *fracrefb_11_p,flcl_ndarray_t *absa_11_p,flcl_ndarray_t *ka_11_p,flcl_ndarray_t *absb_11_p,flcl_ndarray_t *kb_11_p,      
    flcl_ndarray_t *ka_mo2_11_p,flcl_ndarray_t *kb_mo2_11_p,flcl_ndarray_t *selfref_11_p,flcl_ndarray_t *forref_11_p,flcl_ndarray_t *fracrefa_12_p,flcl_ndarray_t *absa_12_p,flcl_ndarray_t *ka_12_p,flcl_ndarray_t *selfref_12_p,flcl_ndarray_t *forref_12_p,           
    flcl_ndarray_t *fracrefa_13_p,flcl_ndarray_t *fracrefb_13_p,flcl_ndarray_t *absa_13_p,flcl_ndarray_t *ka_13_p,flcl_ndarray_t *ka_mco2_13_p,flcl_ndarray_t *ka_mco_13_p,flcl_ndarray_t *kb_mo3_13_p,flcl_ndarray_t *selfref_13_p,flcl_ndarray_t *forref_13_p,                 
    flcl_ndarray_t *fracrefa_14_p,flcl_ndarray_t *fracrefb_14_p,flcl_ndarray_t *absa_14_p,flcl_ndarray_t *ka_14_p,flcl_ndarray_t *absb_14_p,flcl_ndarray_t *kb_14_p,flcl_ndarray_t *selfref_14_p,flcl_ndarray_t *forref_14_p,                                     
    flcl_ndarray_t *fracrefa_15_p,flcl_ndarray_t *absa_15_p,flcl_ndarray_t *ka_15_p,flcl_ndarray_t *ka_mn2_15_p,flcl_ndarray_t *selfref_15_p,flcl_ndarray_t *forref_15_p,flcl_ndarray_t *fracrefa_16_p,flcl_ndarray_t *fracrefb_16_p,flcl_ndarray_t *absa_16_p,flcl_ndarray_t *ka_16_p,flcl_ndarray_t *absb_16_p,flcl_ndarray_t *kb_16_p,      
    flcl_ndarray_t *selfref_16_p,flcl_ndarray_t *forref_16_p                                      ){
    auto emiss_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*emiss_p); auto emiss = real2d_u(emiss_raw,{ims,jms});
    auto tsk_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tsk_p); auto tsk = real2d_u(tsk_raw,{ims,jms});
    auto snow_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*snow_p); auto snow = real2d_u(snow_raw,{ims,jms});
    auto xice_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*xice_p); auto xice = real2d_u(xice_raw,{ims,jms});
    auto xland_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*xland_p); auto xland = real2d_u(xland_raw,{ims,jms});
    auto t3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*t3d_p); auto t3d = real3d_u(t3d_raw,{ims,kms,jms});
    auto p3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*p3d_p); auto p3d = real3d_u(p3d_raw,{ims,kms,jms});
    auto pi3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*pi3d_p); auto pi3d = real3d_u(pi3d_raw,{ims,kms,jms});  
    auto dz8w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*dz8w_p); auto dz8w = real3d_u(dz8w_raw,{ims,kms,jms});
    auto p8w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*p8w_p); auto p8w = real3d_u(p8w_raw,{ims,kms,jms});
    auto t8w_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*t8w_p); auto t8w = real3d_u(t8w_raw,{ims,kms,jms});
    auto re_cloud_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*re_cloud_p); auto re_cloud = real3d_u(re_cloud_raw,{ims,kms,jms});
    auto re_ice_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*re_ice_p); auto re_ice = real3d_u(re_ice_raw,{ims,kms,jms});
    auto re_snow_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*re_snow_p); auto re_snow = real3d_u(re_snow_raw,{ims,kms,jms});
    real3d_u cldfra3d;
    if(cldfra3d_p->data != nullptr){auto cldfra3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*cldfra3d_p); cldfra3d = real3d_u(cldfra3d_raw,{ims,kms,jms});}
    real3d_u qv3d;
    if(qv3d_p->data != nullptr){auto qv3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*qv3d_p); qv3d = real3d_u(qv3d_raw,{ims,kms,jms});}
    real3d_u qc3d;
    if(qc3d_p->data != nullptr){auto qc3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*qc3d_p); qc3d = real3d_u(qc3d_raw,{ims,kms,jms});}
    real3d_u qr3d;
    if(qr3d_p->data != nullptr){auto qr3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*qr3d_p); qr3d = real3d_u(qr3d_raw,{ims,kms,jms});}
    real3d_u qi3d;
    if(qi3d_p->data != nullptr){auto qi3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*qi3d_p); qi3d = real3d_u(qi3d_raw,{ims,kms,jms});}
    real3d_u qs3d;
    if(qs3d_p->data != nullptr){auto qs3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*qs3d_p); qs3d = real3d_u(qs3d_raw,{ims,kms,jms});}
    real3d_u qg3d;
    if(qg3d_p->data != nullptr){auto qg3d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*qg3d_p); qg3d = real3d_u(qg3d_raw,{ims,kms,jms});}
    real3d_u o33d;
    if(o33d_p->data != nullptr){auto o33d_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*o33d_p); o33d = real3d_u(o33d_raw,{ims,kms,jms});}
    real1d_u pin;
    if(pin_p->data != nullptr){auto pin_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*pin_p); pin = real1d_u(pin_raw,{1});}
    auto pprof_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*pprof_p); auto pprof = real1d_u(pprof_raw,{1});
    auto tprof_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*tprof_p); auto tprof = real1d_u(tprof_raw,{1});
    auto retab_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*retab_p); auto retab = real1d_u(retab_raw,{1});
    real3d_u o3clim;
    if(o3clim_p->data != nullptr){auto o3clim_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*o3clim_p); o3clim = real3d_u(o3clim_raw,{ims,1,jms});}
    real4d_u tauaerlw;
    if(tauaerlw_p->data != nullptr){auto tauaerlw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*tauaerlw_p); tauaerlw = real4d_u(tauaerlw_raw,{ims,kms,jms,1});}
    auto glw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*glw_p); auto glw = real2d_u(glw_raw,{ims,jms});
    auto olr_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*olr_p); auto olr = real2d_u(olr_raw,{ims,jms});
    auto lwcf_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwcf_p); auto lwcf = real2d_u(lwcf_raw,{ims,jms});
    real2d_u lwupt;
    if(lwupt_p->data != nullptr){auto lwupt_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwupt_p); lwupt = real2d_u(lwupt_raw,{ims,jms});}
    real2d_u lwuptc;
    if(lwuptc_p->data != nullptr){auto lwuptc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwuptc_p); lwuptc = real2d_u(lwuptc_raw,{ims,jms});}
    real2d_u lwdnt;
    if(lwdnt_p->data != nullptr){auto lwdnt_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwdnt_p); lwdnt = real2d_u(lwdnt_raw,{ims,jms});}
    real2d_u lwdntc;
    if(lwdntc_p->data != nullptr){auto lwdntc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwdntc_p); lwdntc = real2d_u(lwdntc_raw,{ims,jms});}
    real2d_u lwupb;
    if(lwupb_p->data != nullptr){auto lwupb_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwupb_p); lwupb = real2d_u(lwupb_raw,{ims,jms});}
    real2d_u lwupbc;
    if(lwupbc_p->data != nullptr){auto lwupbc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwupbc_p); lwupbc = real2d_u(lwupbc_raw,{ims,jms});}
    real2d_u lwdnb;
    if(lwdnb_p->data != nullptr){auto lwdnb_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwdnb_p); lwdnb = real2d_u(lwdnb_raw,{ims,jms});}
    real2d_u lwdnbc;
    if(lwdnbc_p->data != nullptr){auto lwdnbc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*lwdnbc_p); lwdnbc = real2d_u(lwdnbc_raw,{ims,jms});}
    auto rthratenlw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*rthratenlw_p); auto rthratenlw = real3d_u(rthratenlw_raw,{ims,kms,jms});
    real3d_u rre_cloud;
    if(rre_cloud_p->data != nullptr){auto rre_cloud_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*rre_cloud_p); rre_cloud = real3d_u(rre_cloud_raw,{ims,kms,jms});}
    real3d_u rre_ice;
    if(rre_ice_p->data != nullptr){auto rre_ice_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*rre_ice_p); rre_ice = real3d_u(rre_ice_raw,{ims,kms,jms});}
    real3d_u rre_snow;
    if(rre_snow_p->data != nullptr){auto rre_snow_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*rre_snow_p); rre_snow = real3d_u(rre_snow_raw,{ims,kms,jms});}
    real3d_u lwupflx;
    if(lwupflx_p->data != nullptr){auto lwupflx_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*lwupflx_p); lwupflx = real3d_u(lwupflx_raw,{ims,kms,jms});}
    real3d_u lwupflxc;
    if(lwupflxc_p->data != nullptr){auto lwupflxc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*lwupflxc_p); lwupflxc = real3d_u(lwupflxc_raw,{ims,kms,jms});}
    real3d_u lwdnflx;
    if(lwdnflx_p->data != nullptr){auto lwdnflx_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*lwdnflx_p); lwdnflx = real3d_u(lwdnflx_raw,{ims,kms,jms});}
    real3d_u lwdnflxc;
    if(lwdnflxc_p->data != nullptr){auto lwdnflxc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*lwdnflxc_p); lwdnflxc = real3d_u(lwdnflxc_raw,{ims,kms,jms});}
    auto tau_tbl_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*tau_tbl_p); auto tau_tbl = real1d_u(tau_tbl_raw,{0});
    auto tfn_tbl_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*tfn_tbl_p); auto tfn_tbl = real1d_u(tfn_tbl_raw,{0});
    auto exp_tbl_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*exp_tbl_p); auto exp_tbl = real1d_u(exp_tbl_raw,{0});
    auto fracrefa_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefa_01_p); auto fracrefa_01 = real1d_u(fracrefa_01_raw,{1});
    auto fracrefb_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_01_p); auto fracrefb_01 = real1d_u(fracrefb_01_raw,{1});
    auto absa_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_01_p); auto absa_01 = real2d_u(absa_01_raw,{1,1});
    auto ka_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_01_p); auto ka_01 = real3d_u(ka_01_raw,{1,1,1});
    auto absb_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_01_p); auto absb_01 = real2d_u(absb_01_raw,{1,1});
    auto kb_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_01_p); auto kb_01 = real3d_u(kb_01_raw,{1,1,1});
    auto ka_mn2_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ka_mn2_01_p); auto ka_mn2_01 = real2d_u(ka_mn2_01_raw,{1,1});
    auto kb_mn2_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kb_mn2_01_p); auto kb_mn2_01 = real2d_u(kb_mn2_01_raw,{1,1});
    auto selfref_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_01_p); auto selfref_01 = real2d_u(selfref_01_raw,{1,1});
    auto forref_01_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_01_p); auto forref_01 = real2d_u(forref_01_raw,{1,1});
    auto fracrefa_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefa_02_p); auto fracrefa_02 = real1d_u(fracrefa_02_raw,{1});
    auto fracrefb_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_02_p); auto fracrefb_02 = real1d_u(fracrefb_02_raw,{1});
    auto absa_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_02_p); auto absa_02 = real2d_u(absa_02_raw,{1,1});
    auto ka_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_02_p); auto ka_02 = real3d_u(ka_02_raw,{1,1,1});
    auto absb_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_02_p); auto absb_02 = real2d_u(absb_02_raw,{1,1});
    auto kb_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_02_p); auto kb_02 = real3d_u(kb_02_raw,{1,1,1});
    auto selfref_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_02_p); auto selfref_02 = real2d_u(selfref_02_raw,{1,1});
    auto forref_02_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_02_p); auto forref_02 = real2d_u(forref_02_raw,{1,1});
    auto fracrefa_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_03_p); auto fracrefa_03 = real2d_u(fracrefa_03_raw,{1,1});
    auto fracrefb_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefb_03_p); auto fracrefb_03 = real2d_u(fracrefb_03_raw,{1,1});
    auto absa_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_03_p); auto absa_03 = real2d_u(absa_03_raw,{1,1});
    auto ka_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_03_p); auto ka_03 = real4d_u(ka_03_raw,{1,1,1,1});
    auto absb_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_03_p); auto absb_03 = real2d_u(absb_03_raw,{1,1});
    auto kb_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*kb_03_p); auto kb_03 = real4d_u(kb_03_raw,{1,1,1,1});
    auto ka_mn2o_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_mn2o_03_p); auto ka_mn2o_03 = real3d_u(ka_mn2o_03_raw,{1,1,1});
    auto kb_mn2o_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_mn2o_03_p); auto kb_mn2o_03 = real3d_u(kb_mn2o_03_raw,{1,1,1});
    auto selfref_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_03_p); auto selfref_03 = real2d_u(selfref_03_raw,{1,1});
    auto forref_03_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_03_p); auto forref_03 = real2d_u(forref_03_raw,{1,1});
    auto fracrefa_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_04_p); auto fracrefa_04 = real2d_u(fracrefa_04_raw,{1,1});
    auto fracrefb_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefb_04_p); auto fracrefb_04 = real2d_u(fracrefb_04_raw,{1,1});
    auto absa_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_04_p); auto absa_04 = real2d_u(absa_04_raw,{1,1});
    auto ka_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_04_p); auto ka_04 = real4d_u(ka_04_raw,{1,1,1,1});
    auto absb_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_04_p); auto absb_04 = real2d_u(absb_04_raw,{1,1});
    auto kb_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*kb_04_p); auto kb_04 = real4d_u(kb_04_raw,{1,1,1,1});
    auto selfref_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_04_p); auto selfref_04 = real2d_u(selfref_04_raw,{1,1});
    auto forref_04_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_04_p); auto forref_04 = real2d_u(forref_04_raw,{1,1});
    auto fracrefa_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_05_p); auto fracrefa_05 = real2d_u(fracrefa_05_raw,{1,1});
    auto fracrefb_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefb_05_p); auto fracrefb_05 = real2d_u(fracrefb_05_raw,{1,1});
    auto absa_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_05_p); auto absa_05 = real2d_u(absa_05_raw,{1,1});
    auto ka_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_05_p); auto ka_05 = real4d_u(ka_05_raw,{1,1,1,1});
    auto absb_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_05_p); auto absb_05 = real2d_u(absb_05_raw,{1,1});
    auto kb_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*kb_05_p); auto kb_05 = real4d_u(kb_05_raw,{1,1,1,1});
    auto ka_mo3_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_mo3_05_p); auto ka_mo3_05 = real3d_u(ka_mo3_05_raw,{1,1,1});
    auto selfref_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_05_p); auto selfref_05 = real2d_u(selfref_05_raw,{1,1});
    auto forref_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_05_p); auto forref_05 = real2d_u(forref_05_raw,{1,1});
    auto ccl4_05_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*ccl4_05_p); auto ccl4_05 = real1d_u(ccl4_05_raw,{1});
    auto fracrefa_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefa_06_p); auto fracrefa_06 = real1d_u(fracrefa_06_raw,{1});
    auto absa_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_06_p); auto absa_06 = real2d_u(absa_06_raw,{1,1});
    auto ka_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_06_p); auto ka_06 = real3d_u(ka_06_raw,{1,1,1});
    auto ka_mco2_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ka_mco2_06_p); auto ka_mco2_06 = real2d_u(ka_mco2_06_raw,{1,1});
    auto selfref_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_06_p); auto selfref_06 = real2d_u(selfref_06_raw,{1,1});
    auto forref_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_06_p); auto forref_06 = real2d_u(forref_06_raw,{1,1});
    auto cfc11adj_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*cfc11adj_06_p); auto cfc11adj_06 = real1d_u(cfc11adj_06_raw,{1});
    auto cfc12_06_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*cfc12_06_p); auto cfc12_06 = real1d_u(cfc12_06_raw,{1});
    auto fracrefa_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_07_p); auto fracrefa_07 = real2d_u(fracrefa_07_raw,{1,1});
    auto fracrefb_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_07_p); auto fracrefb_07 = real1d_u(fracrefb_07_raw,{1});
    auto absa_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_07_p); auto absa_07 = real2d_u(absa_07_raw,{1,1});
    auto ka_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_07_p); auto ka_07 = real4d_u(ka_07_raw,{1,1,1,1});
    auto absb_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_07_p); auto absb_07 = real2d_u(absb_07_raw,{1,1});
    auto kb_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_07_p); auto kb_07 = real3d_u(kb_07_raw,{1,1,1});
    auto ka_mco2_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_mco2_07_p); auto ka_mco2_07 = real3d_u(ka_mco2_07_raw,{1,1,1});
    auto kb_mco2_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kb_mco2_07_p); auto kb_mco2_07 = real2d_u(kb_mco2_07_raw,{1,1});
    auto selfref_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_07_p); auto selfref_07 = real2d_u(selfref_07_raw,{1,1});
    auto forref_07_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_07_p); auto forref_07 = real2d_u(forref_07_raw,{1,1});
    auto fracrefa_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefa_08_p); auto fracrefa_08 = real1d_u(fracrefa_08_raw,{1});
    auto fracrefb_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_08_p); auto fracrefb_08 = real1d_u(fracrefb_08_raw,{1});
    auto absa_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_08_p); auto absa_08 = real2d_u(absa_08_raw,{1,1});
    auto ka_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_08_p); auto ka_08 = real3d_u(ka_08_raw,{1,1,1});
    auto absb_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_08_p); auto absb_08 = real2d_u(absb_08_raw,{1,1});
    auto kb_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_08_p); auto kb_08 = real3d_u(kb_08_raw,{1,1,1});
    auto ka_mco2_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ka_mco2_08_p); auto ka_mco2_08 = real2d_u(ka_mco2_08_raw,{1,1});
    auto ka_mn2o_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ka_mn2o_08_p); auto ka_mn2o_08 = real2d_u(ka_mn2o_08_raw,{1,1});
    auto ka_mo3_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ka_mo3_08_p); auto ka_mo3_08 = real2d_u(ka_mo3_08_raw,{1,1});
    auto kb_mco2_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kb_mco2_08_p); auto kb_mco2_08 = real2d_u(kb_mco2_08_raw,{1,1});
    auto kb_mn2o_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kb_mn2o_08_p); auto kb_mn2o_08 = real2d_u(kb_mn2o_08_raw,{1,1});
    auto selfref_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_08_p); auto selfref_08 = real2d_u(selfref_08_raw,{1,1});
    auto forref_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_08_p); auto forref_08 = real2d_u(forref_08_raw,{1,1});
    auto cfc12_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*cfc12_08_p); auto cfc12_08 = real1d_u(cfc12_08_raw,{1});
    auto cfc22adj_08_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*cfc22adj_08_p); auto cfc22adj_08 = real1d_u(cfc22adj_08_raw,{1});
    auto fracrefa_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_09_p); auto fracrefa_09 = real2d_u(fracrefa_09_raw,{1,1});
    auto fracrefb_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_09_p); auto fracrefb_09 = real1d_u(fracrefb_09_raw,{1});
    auto absa_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_09_p); auto absa_09 = real2d_u(absa_09_raw,{1,1});
    auto ka_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_09_p); auto ka_09 = real4d_u(ka_09_raw,{1,1,1,1});
    auto absb_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_09_p); auto absb_09 = real2d_u(absb_09_raw,{1,1});
    auto kb_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_09_p); auto kb_09 = real3d_u(kb_09_raw,{1,1,1});
    auto ka_mn2o_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_mn2o_09_p); auto ka_mn2o_09 = real3d_u(ka_mn2o_09_raw,{1,1,1});
    auto kb_mn2o_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kb_mn2o_09_p); auto kb_mn2o_09 = real2d_u(kb_mn2o_09_raw,{1,1});
    auto selfref_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_09_p); auto selfref_09 = real2d_u(selfref_09_raw,{1,1});
    auto forref_09_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_09_p); auto forref_09 = real2d_u(forref_09_raw,{1,1});
    auto fracrefa_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefa_10_p); auto fracrefa_10 = real1d_u(fracrefa_10_raw,{1});
    auto fracrefb_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_10_p); auto fracrefb_10 = real1d_u(fracrefb_10_raw,{1});
    auto absa_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_10_p); auto absa_10 = real2d_u(absa_10_raw,{1,1});
    auto ka_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_10_p); auto ka_10 = real3d_u(ka_10_raw,{1,1,1});
    auto absb_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_10_p); auto absb_10 = real2d_u(absb_10_raw,{1,1});
    auto kb_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_10_p); auto kb_10 = real3d_u(kb_10_raw,{1,1,1});
    auto selfref_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_10_p); auto selfref_10 = real2d_u(selfref_10_raw,{1,1});
    auto forref_10_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_10_p); auto forref_10 = real2d_u(forref_10_raw,{1,1});
    auto fracrefa_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefa_11_p); auto fracrefa_11 = real1d_u(fracrefa_11_raw,{1});
    auto fracrefb_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_11_p); auto fracrefb_11 = real1d_u(fracrefb_11_raw,{1});
    auto absa_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_11_p); auto absa_11 = real2d_u(absa_11_raw,{1,1});
    auto ka_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_11_p); auto ka_11 = real3d_u(ka_11_raw,{1,1,1});
    auto absb_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_11_p); auto absb_11 = real2d_u(absb_11_raw,{1,1});
    auto kb_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_11_p); auto kb_11 = real3d_u(kb_11_raw,{1,1,1});
    auto ka_mo2_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*ka_mo2_11_p); auto ka_mo2_11 = real2d_u(ka_mo2_11_raw,{1,1});
    auto kb_mo2_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kb_mo2_11_p); auto kb_mo2_11 = real2d_u(kb_mo2_11_raw,{1,1});
    auto selfref_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_11_p); auto selfref_11 = real2d_u(selfref_11_raw,{1,1});
    auto forref_11_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_11_p); auto forref_11 = real2d_u(forref_11_raw,{1,1});
    auto fracrefa_12_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_12_p); auto fracrefa_12 = real2d_u(fracrefa_12_raw,{1,1});
    auto absa_12_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_12_p); auto absa_12 = real2d_u(absa_12_raw,{1,1});
    auto ka_12_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_12_p); auto ka_12 = real4d_u(ka_12_raw,{1,1,1,1});
    auto selfref_12_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_12_p); auto selfref_12 = real2d_u(selfref_12_raw,{1,1});
    auto forref_12_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_12_p); auto forref_12 = real2d_u(forref_12_raw,{1,1});
    auto fracrefa_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_13_p); auto fracrefa_13 = real2d_u(fracrefa_13_raw,{1,1});
    auto fracrefb_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_13_p); auto fracrefb_13 = real1d_u(fracrefb_13_raw,{1});
    auto absa_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_13_p); auto absa_13 = real2d_u(absa_13_raw,{1,1});
    auto ka_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_13_p); auto ka_13 = real4d_u(ka_13_raw,{1,1,1,1});
    auto ka_mco2_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_mco2_13_p); auto ka_mco2_13 = real3d_u(ka_mco2_13_raw,{1,1,1});
    auto ka_mco_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_mco_13_p); auto ka_mco_13 = real3d_u(ka_mco_13_raw,{1,1,1});
    auto kb_mo3_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*kb_mo3_13_p); auto kb_mo3_13 = real2d_u(kb_mo3_13_raw,{1,1});
    auto selfref_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_13_p); auto selfref_13 = real2d_u(selfref_13_raw,{1,1});
    auto forref_13_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_13_p); auto forref_13 = real2d_u(forref_13_raw,{1,1});
    auto fracrefa_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefa_14_p); auto fracrefa_14 = real1d_u(fracrefa_14_raw,{1});
    auto fracrefb_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_14_p); auto fracrefb_14 = real1d_u(fracrefb_14_raw,{1});
    auto absa_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_14_p); auto absa_14 = real2d_u(absa_14_raw,{1,1});
    auto ka_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_14_p); auto ka_14 = real3d_u(ka_14_raw,{1,1,1});
    auto absb_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_14_p); auto absb_14 = real2d_u(absb_14_raw,{1,1});
    auto kb_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_14_p); auto kb_14 = real3d_u(kb_14_raw,{1,1,1});
    auto selfref_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_14_p); auto selfref_14 = real2d_u(selfref_14_raw,{1,1});
    auto forref_14_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_14_p); auto forref_14 = real2d_u(forref_14_raw,{1,1});
    auto fracrefa_15_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_15_p); auto fracrefa_15 = real2d_u(fracrefa_15_raw,{1,1});
    auto absa_15_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_15_p); auto absa_15 = real2d_u(absa_15_raw,{1,1});
    auto ka_15_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_15_p); auto ka_15 = real4d_u(ka_15_raw,{1,1,1,1});
    auto ka_mn2_15_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka_mn2_15_p); auto ka_mn2_15 = real3d_u(ka_mn2_15_raw,{1,1,1});
    auto selfref_15_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_15_p); auto selfref_15 = real2d_u(selfref_15_raw,{1,1});
    auto forref_15_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_15_p); auto forref_15 = real2d_u(forref_15_raw,{1,1});
    auto fracrefa_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*fracrefa_16_p); auto fracrefa_16 = real2d_u(fracrefa_16_raw,{1,1});
    auto fracrefb_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*fracrefb_16_p); auto fracrefb_16 = real1d_u(fracrefb_16_raw,{1});
    auto absa_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa_16_p); auto absa_16 = real2d_u(absa_16_raw,{1,1});
    auto ka_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka_16_p); auto ka_16 = real4d_u(ka_16_raw,{1,1,1,1});
    auto absb_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb_16_p); auto absb_16 = real2d_u(absb_16_raw,{1,1});
    auto kb_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb_16_p); auto kb_16 = real3d_u(kb_16_raw,{1,1,1});
    auto selfref_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref_16_p); auto selfref_16 = real2d_u(selfref_16_raw,{1,1});
    auto forref_16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref_16_p); auto forref_16 = real2d_u(forref_16_raw,{1,1});

    // int ncol,nlay,icld,inflglw,iceflglw,liqflglw;
    // int iplon,irng,permuteseed;
    // int pcols,pver;
    // int idx_rei;
    // int nlayers;
    int nlayers_all = 110; 
    const real thresh1 = 1.e-9;
    const real thresh2 = 1.e9;
    // real slope,ro,dz,corr,gliqwp,gicewp,gsnowp,gravmks,snow_mass_factor;
    
    #ifdef USE_SERIALIZATION
        #warning "Using serial"
        using real1d = Kokkos::Experimental::OffsetView<real*, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using real2d = Kokkos::Experimental::OffsetView<real**, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using real3d = Kokkos::Experimental::OffsetView<real***, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using real4d = Kokkos::Experimental::OffsetView<real****, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using real5d = Kokkos::Experimental::OffsetView<real*****, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using int2d = Kokkos::Experimental::OffsetView<int**, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using int5d = Kokkos::Experimental::OffsetView<int*****, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using int1d = Kokkos::Experimental::OffsetView<int*, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        using int3d = Kokkos::Experimental::OffsetView<int***, Kokkos::LayoutLeft,flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
        real3d tsfc,landfrac,landm,snowh,icefrac,tten1d,cldfra1d,dz1d,p1d,t1d,qv1d,qc1d,qr1d,qi1d,qs1d,qg1d,o31d,pw1d,tw1d,
        o3clim1d,o3mmr,varint,plev_col;
        real4d pdel,cliqwp,cicewp,csnowp,reliq,reice,resnow,emis,plev,tlev,uflx,dflx,uflxc,dflxc,
        play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr,cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,clwpth,ciwpth,
        cswpth,rel,rei,res,cldfrac,relqmcl,reicmcl,resnmcl,hr,hrc,fracs,taug,taut,cldfmc,ciwpmc,cswpmc,clwpmc,taucmc,wkl,wx,taua,planklay,planklev,efclfrac,odcld;
        real5d taucld,cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl,tauaer;
        int2d mpas_nlay;
        real tsfc_s[2 * 2 * 2];
        tsfc = real3d((real*)tsfc_s,{1,1,1},{2,2,2});
        real landfrac_s[2 * 2 * 2];
        landfrac = real3d((real*)landfrac_s,{1,1,1},{2,2,2});
        real landm_s[2 * 2 * 2];
        landm = real3d((real*)landm_s,{1,1,1},{2,2,2});
        real snowh_s[2 * 2 * 2];
        snowh = real3d((real*)snowh_s,{1,1,1},{2,2,2});
        real icefrac_s[2 * 2 * 2];
        icefrac = real3d((real*)icefrac_s,{1,1,1},{2,2,2});
        real tten1d_s[(kte + 1) * 2 * 2];
        tten1d = real3d((real*)tten1d_s,{kts,1,1},{kte+1,2,2});
        real cldfra1d_s[(kte + 1) * 2 * 2];
        cldfra1d = real3d((real*)cldfra1d_s,{kts,1,1},{kte+1,2,2});
        real dz1d_s[(kte + 1) * 2 * 2];
        dz1d = real3d((real*)dz1d_s,{kts,1,1},{kte+1,2,2});
        real p1d_s[(kte + 1) * 2 * 2];
        p1d = real3d((real*)p1d_s,{kts,1,1},{kte+1,2,2});
        real t1d_s[(kte + 1) * 2 * 2];
        t1d = real3d((real*)t1d_s,{kts,1,1},{kte+1,2,2});
        real qv1d_s[(kte + 1) * 2 * 2];
        qv1d = real3d((real*)qv1d_s,{kts,1,1},{kte+1,2,2});
        real qc1d_s[(kte + 1) * 2 * 2];
        qc1d = real3d((real*)qc1d_s,{kts,1,1},{kte+1,2,2});
        real qr1d_s[(kte + 1) * 2 * 2];
        qr1d = real3d((real*)qr1d_s,{kts,1,1},{kte+1,2,2});
        real qi1d_s[(kte + 1) * 2 * 2];
        qi1d = real3d((real*)qi1d_s,{kts,1,1},{kte+1,2,2});
        real qs1d_s[(kte + 1) * 2 * 2];
        qs1d = real3d((real*)qs1d_s,{kts,1,1},{kte+1,2,2});
        real qg1d_s[(kte + 1) * 2 * 2];
        qg1d = real3d((real*)qg1d_s,{kts,1,1},{kte+1,2,2});
        real o31d_s[(kte + 1) * 2 * 2];
        o31d = real3d((real*)o31d_s,{kts,1,1},{kte+1,2,2});
        real pw1d_s[(kte + 2) * 2 * 2];
        pw1d = real3d((real*)pw1d_s,{kts,1,1},{kte+2,2,2});
        real tw1d_s[(kte + 2) * 2 * 2];
        tw1d = real3d((real*)tw1d_s,{kts,1,1},{kte+2,2,2});
        real o3clim1d_s[(noznlevels + 1) * 2 * 2];
        o3clim1d = real3d((real*)o3clim1d_s,{1,1,1},{noznlevels+1,2,2});
        real pdel_s[(kte-kts+2) * 2 * 2 * 2];
        pdel = real4d((real*)pdel_s,{1,1,1,1},{2,kte-kts+2,2,2});
        real cliqwp_s[(kte-kts+2) * 2 * 2 * 2];
        cliqwp = real4d((real*)cliqwp_s,{1,1,1,1},{2,kte-kts+2,2,2});
        real cicewp_s[(kte-kts+2) * 2 * 2 * 2];
        cicewp = real4d((real*)cicewp_s,{1,1,1,1},{2,kte-kts+2,2,2});
        real csnowp_s[(kte-kts+2) * 2 * 2 * 2];
        csnowp = real4d((real*)csnowp_s,{1,1,1,1},{2,kte-kts+2,2,2});
        real reliq_s[(kte-kts+2) * 2 * 2 * 2];
        reliq = real4d((real*)reliq_s,{1,1,1,1},{2,kte-kts+2,2,2});
        real reice_s[(kte-kts+2) * 2 * 2 * 2];
        reice = real4d((real*)reice_s,{1,1,1,1},{2,kte-kts+2,2,2});
        real resnow_s[(kte-kts+2) * 2 * 2 * 2];
        resnow = real4d((real*)resnow_s,{1,1,1,1},{2,kte-kts+2,2,2});
        real emis_s[(nbndlw + 1) * 2 * 2 * 2];
        emis = real4d((real*)emis_s,{1,1,1,1},{2,nbndlw+1,2,2});
        int mpas_nlay_s[(ite+1) * (jte+1)];
        mpas_nlay = int2d((int*)mpas_nlay_s,{its,jts},{ite+1,jte+1});
        real o3mmr_s[(nlayers_all + 1) * 2 * 2];
        o3mmr = real3d((real*)o3mmr_s,{1,1,1},{nlayers_all+1,2,2});
        real varint_s[(nlayers_all + 2) * 2 * 2];
        varint = real3d((real*)varint_s,{1,1,1},{nlayers_all+2,2,2});
        real plev_s[(nlayers_all + 2) * 2 * 2 * 2];
        plev = real4d((real*)plev_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real tlev_s[(nlayers_all + 2) * 2 * 2 * 2];
        tlev = real4d((real*)tlev_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real uflx_s[(nlayers_all + 2) * 2 * 2 * 2];
        uflx = real4d((real*)uflx_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real uflxc_s[(nlayers_all + 2) * 2 * 2 * 2];
        uflxc = real4d((real*)uflxc_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real dflx_s[(nlayers_all + 2) * 2 * 2 * 2];
        dflx = real4d((real*)dflx_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real dflxc_s[(nlayers_all + 2) * 2 * 2 * 2];
        dflxc = real4d((real*)dflxc_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real play_s[(nlayers_all + 2) * 2 * 2 * 2];
        play = real4d((real*)play_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real tlay_s[(nlayers_all + 2) * 2 * 2 * 2];
        tlay = real4d((real*)tlay_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real h2ovmr_s[(nlayers_all + 2) * 2 * 2 * 2];
        h2ovmr = real4d((real*)h2ovmr_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real o3vmr_s[(nlayers_all + 2) * 2 * 2 * 2];
        o3vmr = real4d((real*)o3vmr_s,{1,1,1,1},{2,nlayers_all+2,2,2});
        real co2vmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        co2vmr = real4d((real*)co2vmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real o2vmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        o2vmr = real4d((real*)o2vmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real ch4vmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        ch4vmr = real4d((real*)ch4vmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real n2ovmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        n2ovmr = real4d((real*)n2ovmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real cfc11vmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        cfc11vmr = real4d((real*)cfc11vmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real cfc12vmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        cfc12vmr = real4d((real*)cfc12vmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real cfc22vmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        cfc22vmr = real4d((real*)cfc22vmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real ccl4vmr_s[(nlayers_all + 1) * 2 * 2 * 2];
        ccl4vmr = real4d((real*)ccl4vmr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real clwpth_s[(nlayers_all + 1) * 2 * 2 * 2];
        clwpth = real4d((real*)clwpth_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real ciwpth_s[(nlayers_all + 1) * 2 * 2 * 2];
        ciwpth = real4d((real*)ciwpth_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real cswpth_s[(nlayers_all + 1) * 2 * 2 * 2];
        cswpth = real4d((real*)cswpth_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real rel_s[(nlayers_all + 1) * 2 * 2 * 2];
        rel = real4d((real*)rel_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real rei_s[(nlayers_all + 1) * 2 * 2 * 2];
        rei = real4d((real*)rei_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real res_s[(nlayers_all + 1) * 2 * 2 * 2];
        res = real4d((real*)res_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real cldfrac_s[(nlayers_all + 1) * 2 * 2 * 2];
        cldfrac = real4d((real*)cldfrac_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real relqmcl_s[(nlayers_all + 1) * 2 * 2 * 2];
        relqmcl = real4d((real*)relqmcl_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real reicmcl_s[(nlayers_all + 1) * 2 * 2 * 2];
        reicmcl = real4d((real*)reicmcl_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real resnmcl_s[(nlayers_all + 1) * 2 * 2 * 2];
        resnmcl = real4d((real*)resnmcl_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real hr_s[(nlayers_all + 1) * 2 * 2 * 2];
        hr = real4d((real*)hr_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real hrc_s[(nlayers_all + 1) * 2 * 2 * 2];
        hrc = real4d((real*)hrc_s,{1,1,1,1},{2,nlayers_all+1,2,2});
        real taucld_s[(nbndlw + 1) * 2 * (nlayers_all + 1) * 2 * 2];
        taucld = real5d((real*)taucld_s,{1,1,1,1,1},{nbndlw+1,2,nlayers_all+1,2,2});
        real cldfmcl_s[(ngptlw + 1) * 2 * (nlayers_all + 1) * 2 * 2];
        cldfmcl = real5d((real*)cldfmcl_s,{1,1,1,1,1},{ngptlw+1,2,nlayers_all+1,2,2});
        real clwpmcl_s[(ngptlw + 1) * 2 * (nlayers_all + 1) * 2 * 2];
        clwpmcl = real5d((real*)clwpmcl_s,{1,1,1,1,1},{ngptlw+1,2,nlayers_all+1,2,2});
        real ciwpmcl_s[(ngptlw + 1) * 2 * (nlayers_all + 1) * 2 * 2];
        ciwpmcl = real5d((real*)ciwpmcl_s,{1,1,1,1,1},{ngptlw+1,2,nlayers_all+1,2,2});
        real cswpmcl_s[(ngptlw + 1) * 2 * (nlayers_all + 1) * 2 * 2];
        cswpmcl = real5d((real*)cswpmcl_s,{1,1,1,1,1},{ngptlw+1,2,nlayers_all+1,2,2});
        real taucmcl_s[(ngptlw + 1) * 2 * (nlayers_all + 1) * 2 * 2];
        taucmcl = real5d((real*)taucmcl_s,{1,1,1,1,1},{ngptlw+1,2,nlayers_all+1,2,2});
        real tauaer_s[(nbndlw + 1) * 2 * (nlayers_all + 1) * 2 * 2];
        tauaer = real5d((real*)tauaer_s,{1,1,1,1,1},{2,nlayers_all+1,nbndlw+1,2,2});
        real plev_col_s[(nlayers_all + 2) * 2 * 2];
        plev_col = real3d((real*)plev_col_s,{1,1,1},{nlayers_all+2,2,2});
        real fracs_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        fracs = real4d((real*)fracs_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real taug_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        taug = real4d((real*)taug_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real taut_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        taut = real4d((real*)taut_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real cldfmc_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        cldfmc = real4d((real*)cldfmc_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real ciwpmc_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        ciwpmc = real4d((real*)ciwpmc_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real clwpmc_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        clwpmc = real4d((real*)clwpmc_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real cswpmc_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        cswpmc = real4d((real*)cswpmc_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real taucmc_s[(ngptlw + 1) * (nlayers_all + 2) * 2 * 2];
        taucmc = real4d((real*)taucmc_s,{1,1,1,1},{ngptlw+1,nlayers_all+2,2,2});
        real wkl_s[(mxmol + 1) * (nlayers_all + 2) * 2 * 2];
        wkl = real4d((real*)wkl_s,{1,1,1,1},{mxmol+1,nlayers_all+2,2,2});
        real wx_s[(maxxsec + 1) * (nlayers_all + 2) * 2 * 2];
        wx = real4d((real*)wx_s,{1,1,1,1},{maxxsec+1,nlayers_all+2,2,2});
        real taua_s[(nbndlw + 1) * (nlayers_all + 2) * 2 * 2];
        taua = real4d((real*)taua_s,{1,1,1,1},{nbndlw+1,nlayers_all+2,2,2});
        real planklay_s[(nlayers_all + 2) * (nbndlw + 1) * 2 * 2];
        planklay = real4d((real*)planklay_s,{1,1,1,1},{nlayers_all+2,nbndlw+1,2,2});
        real planklev_s[(nlayers_all + 2) * (nbndlw + 1) * 2 * 2];
        planklev = real4d((real*)planklev_s,{0,1,1,1},{nlayers_all+2,nbndlw+1,2,2});
        real efclfrac_s[(ngptlw + 1) * (nlayers_all + 1) * 2 * 2];
        efclfrac = real4d((real*)efclfrac_s,{1,1,1,1},{ngptlw+1,nlayers_all+1,2,2});
        real odcld_s[(ngptlw + 1) * (nlayers_all + 1) * 2 * 2];
        odcld = real4d((real*)odcld_s,{1,1,1,1},{ngptlw+1,nlayers_all+1,2,2});
    #endif
    #ifdef USE_SW
        lw_workspace &ws = get_lw_ws(kts, kte, 1, RAD_LW_CHUNK_SIZE, jts, jte,
                                     nbndlw, ngptlw, mxmol, maxxsec, noznlevels);
        auto mpas_nlay = int2d("mpas_nlay",{its,ite+1},{jts,jte+1});
        auto tsfc       = ws.tsfc;
        auto landfrac   = ws.landfrac;
        auto landm      = ws.landm;
        auto snowh      = ws.snowh;
        auto icefrac    = ws.icefrac;

        auto tten1d     = ws.tten1d;
        auto cldfra1d   = ws.cldfra1d;
        auto dz1d       = ws.dz1d;
        auto p1d        = ws.p1d;
        auto t1d        = ws.t1d;
        auto qv1d       = ws.qv1d;
        auto qc1d       = ws.qc1d;
        auto qr1d       = ws.qr1d;
        auto qi1d       = ws.qi1d;
        auto qs1d       = ws.qs1d;
        auto qg1d       = ws.qg1d;
        auto o31d       = ws.o31d;
        auto pw1d       = ws.pw1d;
        auto tw1d       = ws.tw1d;
        auto o3clim1d   = ws.o3clim1d;

        auto pdel       = ws.pdel;
        auto cliqwp     = ws.cliqwp;
        auto cicewp     = ws.cicewp;
        auto csnowp     = ws.csnowp;
        auto reliq      = ws.reliq;
        auto reice      = ws.reice;
        auto resnow     = ws.resnow;
        auto emis       = ws.emis;
       // auto mpas_nlay  = ws.mpas_nlay;
        
        auto o3mmr      = ws.o3mmr;
        auto varint     = ws.varint;
        auto plev       = ws.plev;
        auto tlev       = ws.tlev;
        auto uflx       = ws.uflx;
        auto dflx       = ws.dflx;
        auto uflxc      = ws.uflxc;
        auto dflxc      = ws.dflxc;
        
        auto play       = ws.play;
        auto tlay       = ws.tlay;
        auto h2ovmr     = ws.h2ovmr;
        auto o3vmr      = ws.o3vmr;
        auto co2vmr     = ws.co2vmr;
        auto o2vmr      = ws.o2vmr;
        auto ch4vmr     = ws.ch4vmr;
        auto n2ovmr     = ws.n2ovmr;
        auto cfc11vmr   = ws.cfc11vmr;
        auto cfc12vmr   = ws.cfc12vmr;
        auto cfc22vmr   = ws.cfc22vmr;
        auto ccl4vmr    = ws.ccl4vmr;
        
        auto clwpth     = ws.clwpth;
        auto ciwpth     = ws.ciwpth;
        auto cswpth     = ws.cswpth;
        auto rel        = ws.rel;
        auto rei        = ws.rei;
        auto res        = ws.res;
        auto cldfrac    = ws.cldfrac;
        auto relqmcl    = ws.relqmcl;
        auto reicmcl    = ws.reicmcl;
        auto resnmcl    = ws.resnmcl;
        auto hr         = ws.hr;
        auto hrc        = ws.hrc;
        
        // 5D Arrays (High memory footprint)
        auto taucld     = ws.taucld;
        auto cldfmcl    = ws.cldfmcl;
        auto clwpmcl    = ws.clwpmcl;
        auto ciwpmcl    = ws.ciwpmcl;
        auto cswpmcl    = ws.cswpmcl;
        auto taucmcl    = ws.taucmcl;
        auto tauaer     = ws.tauaer;
        
        auto plev_col   = ws.plev_col;
        auto fracs      = ws.fracs;
        auto taug       = ws.taug;
        auto taut       = ws.taut;
        auto cldfmc     = ws.cldfmc;
        auto ciwpmc     = ws.ciwpmc;
        auto clwpmc     = ws.clwpmc;
        auto cswpmc     = ws.cswpmc;
        auto taucmc     = ws.taucmc;
        auto wkl        = ws.wkl;
        auto wx         = ws.wx;
        auto taua       = ws.taua;
        auto planklay   = ws.planklay;
        auto planklev   = ws.planklev;
        auto efclfrac   = ws.efclfrac;
        auto odcld      = ws.odcld;
    #endif
    // auto bbugas = real3d("bbugas",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto bbutot = real3d("bbutot",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto secdiff = real3d("secdiff",{1,nbndlw+1},{its,ite+1},{jts,jte+1});
    // auto icldlyr = int3d("icldlyr",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto uclfl = real3d("uclfl",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto dclfl = real3d("dclfl",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto clrurad = real3d("clrurad_s",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto clrdrad = real3d("clrdrad",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto atot = real3d("atot",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    // auto atrans = real3d("atrans",{1,nlayers_all+1},{its,ite+1},{jts,jte+1});
    real co2 = 379.e-6;         // carbon dioxide (379 ppmv)
    real ch4 = 1774.e-9;        // methane (1774 ppbv)
    real n2o = 319.e-9;         // nitrous oxide (319 ppbv)
    real cfc11 = 0.251e-9;      // cfc-11 (251 ppt)
    real cfc12 = 0.538e-9;      // cfc-12 (538 ppt)
    real cfc22 = 0.169e-9;      // cfc-22 (169 ppt)
    real ccl4 = 0.093e-9;       // ccl4 (93 ppt)
    real o2 = 0.209488;         // --- set oxygen volume mixing ratio (for o2mmr=0.23143):
    real amdw = 1.607793;       // Molecular weight of dry air / water vapor
    real amdo = 0.603461;       // Molecular weight of dry air / ozone
    real amdo2 = 0.905190;      // Molecular weight of dry air / oxygen
    const real deltap = 4.;     // Pressure interval for buffer layer in mb
    const real rgas = 287.0;    // Constant: Gas constant for dry air [J kg-1 K-1]
    const real cp = 7.*rgas/2.; // Constant: Specific heat of dry air at constant pressure [J kg-1 K-1]
    const real g = 9.80616;     // Constant: Acceleration due to gravity [m s-2]
    // int klev;
    // real wght,vark,vark1;
    //--- all fields are ordered vertically from bottom to top (pressures are in mb):
    //--- initialize option for the calculation of the cloud optical properties:
    for (int chunk_start = its; chunk_start <= ite; chunk_start += RAD_LW_CHUNK_SIZE) {
        int current_chunk_end = std::min(chunk_start + RAD_LW_CHUNK_SIZE - 1, ite);
        int current_chunk_len = current_chunk_end - chunk_start + 1;
        Kokkos::parallel_for(
            "rrtmglw_kokkos", Kokkos::RangePolicy<>(1,current_chunk_len + 1), KOKKOS_LAMBDA(int offset) {
                int m,j,i;
                #ifdef USE_SERIALIZATION
                m = 1; j = 1; i = chunk_start + offset - 1;
                #endif
                #ifdef USE_SW
                m = offset; j =1;
                i = chunk_start + offset - 1;
                #endif
                int ncol = 1; int icld = 2; int inflglw = 2; int iceflglw = 3; int liqflglw = 1;
        // for(int j = jts; j <= jte; j++){
        //   for(int i = its; i<= ite; i++){
                //--- set surface emissivity in each RRTMG longwave band:
                for(int nb =1; nb <= nbndlw; nb++){
                    emis(ncol,nb,m,j) = emiss(i,j);
                }
                //--- INITIALIZE COLUMN SOUNDING (the call to the long wave radiation code is done one column at a time
                //    which is why we set ncol = 1 above):   
                for(int k = kts; k <= kte + 1; k++){
                    pw1d(k,m,j) = p8w(i,k,j) / 100.;
                    tw1d(k,m,j) = t8w(i,k,j);
                }
                for(int k = kts; k <= kte; k++){
                    p1d(k,m,j)      = p3d(i,k,j) / 100.;
                    dz1d(k,m,j)     = dz8w(i,k,j);
                    t1d(k,m,j)      = t3d(i,k,j);
                    qv1d(k,m,j)     = std::max(qv3d(i,k,j),1.e-12);
                    o31d(k,m,j)     = 0.;
                    qc1d(k,m,j)     = 0.;
                    qr1d(k,m,j)     = 0.;
                    qi1d(k,m,j)     = 0.;
                    qs1d(k,m,j)     = 0.;
                    qg1d(k,m,j)     = 0.;
                    cldfra1d(k,m,j) = 0.;
                    if(o33d_p->data != nullptr) o31d(k,m,j) = o33d(i,k,j);
                }
                //--- initialize the local arrays containing the different cloud water and ice condenstates:
                if(icloud > 0) {
                    for(int k = kts; k <= kte; k++){
                        if(qc3d_p->data != nullptr) qc1d(k,m,j) = std::max(qc3d(i,k,j),0.);
                        if(qr3d_p->data != nullptr) qr1d(k,m,j) = std::max(qr3d(i,k,j),0.);
                        if(qi3d_p->data != nullptr) qi1d(k,m,j) = std::max(qi3d(i,k,j),0.);
                        if(qs3d_p->data != nullptr) qs1d(k,m,j) = std::max(qs3d(i,k,j),0.);
                        if(((qc3d_p->data != nullptr) || (qi3d_p->data != nullptr) || (qs3d_p->data != nullptr)) && (cldfra3d_p->data != nullptr)) cldfra1d(k,m,j) = cldfra3d(i,k,j);          
                    }
                }
                
                //--- initialize the local radiative heating rate:
                for(int k = kts; k <= kte; k++){
                    tten1d(k,m,j) = 0.;
                }
                /*--- add extra layers to include absorption and transmission between the top of the model and the top of
                    the atmosphere: in contrast to WRF, MPAS columns have different model-top pressures since MPAS uses
                    a height coordinate system. Therefore, we define nlayers for each individual column:
                    czhao modification +++++
                    nlayers = kte + max(nint(pw1d(kte+1)/deltap), 1)
                */
                int nlayers = kte + std::max(int(std::round(pw1d(kte+1,m,j)/deltap+0.5)),1);
                mpas_nlay(i,j) = nlayers-kte;
                
                int nlay = nlayers; // keep these indices the same.
                //--- initialize local arrays called in the calculation of the cloud optical properties and radiative
                //    fluxes:
                for(int n = 1; n <= ncol; n++){
                    for(int k = kts; k <= kte; k++){
                        reliq(ncol,k,m,j)  = 10.;
                        reice(ncol,k,m,j)  = 10.;
                        resnow(ncol,k,m,j) = 10.;
                        cliqwp(ncol,k,m,j) = 0.;
                        cicewp(ncol,k,m,j) = 0.;
                        csnowp(ncol,k,m,j) = 0.;
                    }
                    for(int k = 1; k <= nlayers; k++){
                        clwpth(n,k,m,j)  = 0.;
                        ciwpth(n,k,m,j)  = 0.;
                        cswpth(n,k,m,j)  = 0.;
                        rel(n,k,m,j)     = 0.;
                        rei(n,k,m,j)     = 0.;
                        res(n,k,m,j)     = 0.;
                        cldfrac(n,k,m,j) = 0.;
                        relqmcl(n,k,m,j) = 0.;
                        reicmcl(n,k,m,j) = 0.;
                        resnmcl(n,k,m,j) = 0.;
                        uflx(n,k,m,j)    = 0.;
                        uflxc(n,k,m,j)   = 0.;
                        dflx(n,k,m,j)    = 0.;
                        dflxc(n,k,m,j)   = 0.;
                        hr(n,k,m,j)      = 0.;
                        hrc(n,k,m,j)     = 0.;
                    }
                }
                for(int k = 1; k <= nlayers; k++){
                    for(int n = 1; n <= ncol; n++){
                        for(int p = 1; p <= nbndlw; p++){
                            taucld(p,n,k,m,j)  = 0.;
                        }
                    }
                }
                for(int p = 1; p <= nbndlw; p++){
                    for(int k = 1; k <= nlayers; k++){
                        for(int n = 1; n <= ncol; n++){
                            tauaer(n,k,p,m,j)  = 0.;
                        }
                    }
                }
                
                    
                for(int k = 1; k <= nlayers; k++){
                    for(int n = 1; n <= ncol; n++){
                        for(int p = 1; p <= ngptlw; p++){
                            cldfmcl(p,n,k,m,j) = 0.;
                            clwpmcl(p,n,k,m,j) = 0.;
                            ciwpmcl(p,n,k,m,j) = 0.;
                            cswpmcl(p,n,k,m,j) = 0.;
                            taucmcl(p,n,k,m,j) = 0.;
                        }
                    }
                } 
                uflx(1,nlayers+1,m,j)  = 0.;
                uflxc(1,nlayers+1,m,j) = 0.;
                dflx(1,nlayers+1,m,j)  = 0.;
                dflxc(1,nlayers+1,m,j) = 0.;
                // //--- fill local arrays with input sounding. convert water vapor mass mixing ratio to volume mixing ratio:
                plev(ncol,1,m,j) = pw1d(1,m,j);
                tlev(ncol,1,m,j) = tw1d(1,m,j);
                tsfc(ncol,m,j)   = tsk(i,j);
                for(int k = kts; k <= kte; k++){
                    plev(ncol,k+1,m,j)   = pw1d(k+1,m,j);
                    tlev(ncol,k+1,m,j)   = tw1d(k+1,m,j);
                    play(ncol,k,m,j)     = p1d(k,m,j);
                    pdel(ncol,k,m,j)     = plev(ncol,k,m,j) - plev(ncol,k+1,m,j);
                    tlay(ncol,k,m,j)     = t1d(k,m,j);
                    h2ovmr(ncol,k,m,j)   = qv1d(k,m,j) * amdw;
                    co2vmr(ncol,k,m,j)   = co2;
                    o2vmr(ncol,k,m,j)    = o2;
                    ch4vmr(ncol,k,m,j)   = ch4;
                    n2ovmr(ncol,k,m,j)   = n2o;
                    cfc11vmr(ncol,k,m,j) = cfc11;
                    cfc12vmr(ncol,k,m,j) = cfc12;
                    cfc22vmr(ncol,k,m,j) = cfc22;
                    ccl4vmr(ncol,k,m,j)  = ccl4;
                }
                //--- the sourcecode below follows Steven Cavallo's method to "fill" the atmospheric layers between the
                //    top of the model and the top of the atmosphere. check if the pressure at the top of the atmosphere
                //    is negative. if negative, set it to zero prior to the calculation of temperatures (tlev and tlay):
                for(int k = kte+1; k <= nlayers; k++){
                    plev(ncol,k+1,m,j) = std::max(plev(ncol,k,m,j) - deltap, 0.00);
                    play(ncol,k,m,j) = 0.5*(plev(ncol,k,m,j) + plev(ncol,k+1,m,j));
                }
                //--- add zero as top level. this gets the temperature max at the stratopause, reducing downward flux
                //    errors in the top levels.  If zero happened to be the top level already, this will add another
                //    level with zero, but will not affect the radiative transfer calculation.
                plev(ncol,nlayers+1,m,j) = 0.00;
                play(ncol,nlayers,m,j) =  0.5*(plev(ncol,nlayers,m,j) + plev(ncol,nlayers+1,m,j));
                for(int k = kts; k <= nlayers + 1; k++){
                    plev_col(k,m,j) = plev(ncol,k,m,j);
                }
                int klev;
                real wght,vark,vark1;
                //--- interpolate the table temperatures to column pressure levels:
                for(int k = 1; k <= nlayers + 1; k++){
                    if(pprof(nproflevs) < plev(ncol,k,m,j)){
                        for(int kk = 2; kk <= nproflevs; kk++){       
                            if(pprof(kk) < plev(ncol,k,m,j)){           
                                klev = kk - 1;
                                break;
                            }
                        }
                    }else{
                        klev = nproflevs;
                    } 

                    if(klev != nproflevs){
                        vark  = tprof(klev); 
                        vark1 = tprof(klev+1);
                        wght  = (plev(ncol,k,m,j)-pprof(klev) )/( pprof(klev+1)-pprof(klev));
                    }else{
                        vark  = tprof(klev); 
                        vark1 = tprof(klev);
                        wght  = 0.0;
                    }
                    varint(k,m,j) = wght*(vark1-vark)+vark; 
                }
                //--- match the interpolated table temperature profile:
                for(int k = kte+1; k <= nlayers + 1; k++){
                    tlev(ncol,k,m,j) = varint(k,m,j) + (tlev(ncol,kte,m,j) - varint(kte,m,j));
                    //if(k .le. nlay) then
                    tlay(ncol,k-1,m,j) = 0.5*(tlev(ncol,k,m,j) + tlev(ncol,k-1,m,j));
                    //endif
                }
                // --- fill water vapor and chemical species volume mixing ratios, except ozone:
                for(int k = kte+1; k <= nlayers; k++){
                    h2ovmr(ncol,k,m,j)   = h2ovmr(ncol,kte,m,j); 
                    co2vmr(ncol,k,m,j)   = co2vmr(ncol,kte,m,j); 
                    o2vmr(ncol,k,m,j)    = o2vmr(ncol,kte,m,j); 
                    ch4vmr(ncol,k,m,j)   = ch4vmr(ncol,kte,m,j); 
                    n2ovmr(ncol,k,m,j)   = n2ovmr(ncol,kte,m,j); 
                    cfc11vmr(ncol,k,m,j) = cfc11vmr(ncol,kte,m,j); 
                    cfc12vmr(ncol,k,m,j) = cfc12vmr(ncol,kte,m,j);
                    cfc22vmr(ncol,k,m,j) = cfc22vmr(ncol,kte,m,j); 
                    ccl4vmr(ncol,k,m,j)  = ccl4vmr(ncol,kte,m,j); 
                }
                // --- initialize the ozone voume mixing ratio:
                inirad(o3mmr,plev_col,kts,nlayers-1,m,j);
                // std::cout<<o3mmr<<std::endl;
                if(o3input == 2){
                    for(int k = 1; k <= noznlevels; k++){
                        o3clim1d(k,m,j) = o3clim(i,k,j);
                    }
                    vinterp_ozn(1,ncol,ncol,nlayers,play,pin,noznlevels,o3clim1d,o3mmr,m,j);
                    // std::cout<<o3mmr<<std::endl;
                    for(int k = kts; k <= nlayers; k++){
                        o3vmr(ncol,k,m,j) = o3mmr(k,m,j);
                    }
                }else{
                    for(int k = kts; k <= nlayers; k++){
                        o3vmr(ncol,k,m,j) = o3mmr(k,m,j)*amdo;
                    }
                }
                // --- CALCULATE CLOUD OPTICAL PROPERTIES:
                if(inflglw > 0){
                    for(int k = kts; k <= kte; k++){
                        cldfrac(ncol,k,m,j) = cldfra1d(k,m,j);
                    }
                    //--- zero out cloud optical properties here (not used when passing physical properties to rrtmg_lw as
                    //  taucld is calculated in radiation:
                    
                    for(int k = kts; k <= kte; k++){
                        for(int nb = 1; nb <= nbndlw; nb++){
                            taucld(nb,ncol,k,m,j) = 0.0;
                        }
                    }
                    int pcols = ncol;
                    int pver  = kte - kts + 1;
                    real gravmks = g;
                    if(has_reqc != 0){
                    // --- fill the effective radius for cloud water with that calculated in the Thompson cloud
                    //     microphysics parameterization:
                        inflglw = 3;
                        for(int k = kts; k <= kte; k++){
                            reliq(ncol,k,m,j) = std::max(2.5,re_cloud(i,k,j)*1.e6);
                            if(reliq(ncol,k,m,j) <= 2.5 && cldfrac(ncol,k,m,j) > 0. && (xland(i,j)-1.5) > 0.){ //--- ocean.
                                reliq(ncol,k,m,j) = 10.5;
                            }else if(reliq(ncol,k,m,j) <= 2.5 && cldfrac(ncol,k,m,j) > 0. && (xland(i,j)-1.5) < 0.){   //--- land.
                                reliq(ncol,k,m,j) = 7.5;
                            }
                        }
                    }else{
                        /*---  calculate the effective radius for cloud water in layers below the model top. we added the
                        !     dimensions pcols, pver so that the calls to subroutines relcalc and reicalc remain the same
                        !     as in WRF although the two dimensions are duplicate and not needed inside the subroutines:
                        */
                        landfrac(ncol,m,j) = 2.-xland(i,j);
                        landm(ncol,m,j)    = landfrac(ncol,m,j);
                        snowh(ncol,m,j)    = 0.001*snow(i,j);
                        icefrac(ncol,m,j)  = xice(i,j);
                        relcalc(ncol,pcols,pver,tlay,landfrac,landm,icefrac,reliq,snowh,m,j);
                        // std::cout<<reliq<<std::endl;
                    }
                    if(has_reqi != 0){
                    //--- fill the effective radius for cloud ice with that calculated in the Thompson cloud
                    //    microphysics parameterization: 
                    inflglw  = 4;
                    iceflglw = 4;
                    for(int k = kts; k <= kte; k++){
                            reice(ncol,k,m,j) = std::max(5.,re_ice(i,k,j)*1.e6);
                            if(reice(ncol,k,m,j) <= 5. && cldfrac(ncol,k,m,j) > 0.){
                                int idx_rei = int(t3d(i,k,j)-179.);
                                idx_rei = std::min(std::max(idx_rei,1),75);
                                real corr = t3d(i,k,j) - int(t3d(i,k,j));
                                reice(ncol,k,m,j) = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr;
                                reice(ncol,k,m,j) = std::max(reice(ncol,k,m,j),5.0);
                            }
                    }
                    }else{
                        /*--- calculate the effective radius for cloud ice in layers below the model top, following
                        !    Kristjansson and Mitchell. limit upper bound of reice for Fu ice parameterization and
                        !    convert effective radius to generalized effective size (*1.0315; Fu, 1996):
                        */
                    reicalc(ncol,pcols,pver,tlay,reice,retab,m,j);
                    //    std::cout<<reice<<std::endl;
                    if(iceflglw == 3){
                            for(int k = kts; k <= kte; k++){
                                reice(ncol,k,m,j) = reice(ncol,k,m,j) * 1.0315;
                                reice(ncol,k,m,j) = std::min(140.0,reice(ncol,k,m,j));
                            }
                        }
                    }
                    if(has_reqs != 0){
                    //--- fill the effective radius for snow with that calculated in the Thompson cloud
                    //  microphysics parameterization:
                        inflglw  = 5;
                        iceflglw = 5;
                        for(int k = kts; k <= kte; k++){
                            resnow(ncol,k,m,j) = std::max(10.,re_snow(i,k,j)*1.e6);
                        }
                    }else{
                        for(int k = kts; k <= kte; k++){
                            resnow(ncol,k,m,j) = 10.;
                        }
                    }
                    //--- calculate the cloud liquid water path in units of g/m2 below the model top:
                    for(int k = kts; k <= kte; k++){
                        real gliqwp = qc1d(k,m,j) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;              //grid box liquid water path.
                        cliqwp(ncol,k,m,j) = gliqwp / std::max(0.01,cldfrac(ncol,k,m,j));              //in-cloud liquid water path.
                    }
                    //--- calculate the cloud ice path in units of g/m2 below the model top:
                    if(iceflglw == 3){
                        for(int k = kts; k <= kte; k++){
                            real gicewp = (qi1d(k,m,j)+qs1d(k,m,j)) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;              //grid box liquid water path.
                            cicewp(ncol,k,m,j) = gicewp / std::max(0.01,cldfrac(ncol,k,m,j));              //in-cloud liquid water path.
                        }
                    }else if(iceflglw >= 4){
                        for(int k = kts; k <= kte; k++){
                            real gicewp = qi1d(k,m,j) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;           //grid box ice water path.
                            cicewp(ncol,k,m,j) = gicewp / std::max(0.01,cldfrac(ncol,k,m,j));                //in-cloud ice water path.
                        }
                    }

                    // --- calculate the snow path in units of g/m2 below the model top:
                    if(iceflglw == 5){
                        for(int k = kts; k <= kte; k++){
                            real snow_mass_factor = 1.0;        
                            if(resnow(ncol,k,m,j) > 130.){
                                snow_mass_factor = (130.0/resnow(ncol,k,m,j))*(130.0/resnow(ncol,k,m,j));
                                resnow(ncol,k,m,j)   = 130.0;
                            }             
                            real gsnowp = qs1d(k,m,j) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;           // grid box snow path.
                            csnowp(ncol,k,m,j) = snow_mass_factor * gsnowp / std::max(0.01,cldfrac(ncol,k,m,j));
                        }
                    }
                    //--- set the cloud physical property between the model top and the top of the atmosphere. do not allow
                    // clouds between the model top and the top of the atmosphere: 
                    for(int k = kts; k <= kte; k++){
                        clwpth(ncol,k,m,j) = cliqwp(ncol,k,m,j);
                        ciwpth(ncol,k,m,j) = cicewp(ncol,k,m,j);
                        cswpth(ncol,k,m,j) = csnowp(ncol,k,m,j);
                        rel(ncol,k,m,j) = reliq(ncol,k,m,j);
                        rei(ncol,k,m,j) = reice(ncol,k,m,j);
                        res(ncol,k,m,j) = resnow(ncol,k,m,j);
                    }
                    for(int k = kte+1; k <= nlayers; k++){
                        clwpth(ncol,k,m,j)  = 0.;
                        ciwpth(ncol,k,m,j)  = 0.;
                        cswpth(ncol,k,m,j)  = 0.;
                        rel(ncol,k,m,j)     = 10.;
                        rei(ncol,k,m,j)     = 10.;
                        res(ncol,k,m,j)     = 10.;
                        cldfrac(ncol,k,m,j) = 0.;
                    }
                    for(int k = kte+1; k <= nlayers; k++){
                        for(int nb = 1; nb <= nbndlw; nb++){
                            taucld(nb,ncol,k,m,j) = 0.;
                        }
                    }
                }
                // --- sub-column generator for McICA:
                int iplon = 1;
                int irng  = 0;
                int permuteseed = 150;
                mcica_subcol_lw(iplon   , ncol    , nlay    , icld    , permuteseed , irng    , play    , 
                    cldfrac , ciwpth  , clwpth  , cswpth  , rei         , rel     , res     , 
                    taucld  , cldfmcl , ciwpmcl , clwpmcl , cswpmcl     , reicmcl , relqmcl , 
                    resnmcl , taucmcl,m,j);
                for(int nb = 1; nb <= nbndlw; nb++){
                    for(int k = kts; k <= nlayers; k++){
                        tauaer(ncol,k,nb,m,j) = 0.;
                    }
                }
                #ifdef DO_CHEMISTRY
                string schemes="";
                if(config_aer_ra_feedback_p && (trim(config_chem_opt_p,schemes) != "no_chemistry")){
                        for(int k = kts; k <= kte; k++){
                            if ( tauaerlw(i,k,j,1) > thresh1 && tauaerlw(i,k,j,16) > thresh1 && 
                                tauaerlw(i,k,j,1) < thresh2 && tauaerlw(i,k,j,16) < thresh2){
                                    for(int l = 1; l <= 16; l++){
                                        tauaer(ncol,k,l,m,j) = tauaerlw(i,k,j,l);
                                    }
                            }
                        }
                        // std::cout<<tauaer<<std::endl;
                        // wig beg
                        for(int nb = 1; nb <= nbndlw; nb++){
                            real slope = 0.0;  //use slope as a sum holder
                            for(int k = kts; k <= kte; k++){
                                slope = slope + tauaer(ncol,k,nb,m,j);
                            }
                            // if( slope < 0.0 ) then
                            //     write(msg,'("ERROR: Negative total lw optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,m,j,nb
                            //     call physics_error_fatal(msg)
                            // else if( slope > 5.0 ) then
                            //     call physics_message("-------------------------")
                            //     write(msg,'("WARNING: Large total lw optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,m,j,nb
                            //     call physics_message(msg)

                            //     call physics_message("Diagnostics 1: k, tauaerlw1, tauaerlw16")
                            //     do k=kts,kte
                            //         write(msg,'(i4,2f8.2)') k, tauaerlw(i,k,j,1), tauaerlw(i,k,j,16) !jwfeng: Fixed a small bug.
                            //         call physics_message(msg)
                            //     end do
                            //     call physics_message("-------------------------")
                            // endif
                        }
                }
                #endif
                //--- CALL TO THE RRTMG LONG WAVE RADIATION MODEL:
                rrtmg_lw(ncol    , nlay    , icld     , play     , plev     , tlay    , 
                    tlev    , tsfc    , h2ovmr   , o3vmr    , co2vmr   , ch4vmr  , 
                    n2ovmr  , o2vmr   , cfc11vmr , cfc12vmr , cfc22vmr , ccl4vmr , 
                    emis    , inflglw , iceflglw , liqflglw , cldfmcl  , taucmcl , 
                    ciwpmcl , clwpmcl , cswpmcl  , reicmcl  , relqmcl  , resnmcl , 
                    tauaer  , uflx    , dflx     , hr       , uflxc    , dflxc   , 
                    hrc, tau_tbl,tfn_tbl,exp_tbl,fracrefa_01, fracrefb_01, absa_01, ka_01, absb_01, kb_01, ka_mn2_01, kb_mn2_01, selfref_01, forref_01,                 
                    fracrefa_02, fracrefb_02, absa_02, ka_02, absb_02, kb_02, selfref_02, forref_02,                                     
                    fracrefa_03, fracrefb_03, absa_03, ka_03, absb_03, kb_03, ka_mn2o_03, kb_mn2o_03, selfref_03, forref_03,               
                    fracrefa_04, fracrefb_04, absa_04, ka_04, absb_04, kb_04, selfref_04, forref_04,                                     
                    fracrefa_05, fracrefb_05, absa_05, ka_05, absb_05, kb_05, ka_mo3_05, selfref_05, forref_05, ccl4_05,                   
                    fracrefa_06, absa_06, ka_06, ka_mco2_06, selfref_06, forref_06, cfc11adj_06, cfc12_06,                           
                    fracrefa_07, fracrefb_07, absa_07, ka_07, absb_07, kb_07, ka_mco2_07, kb_mco2_07, selfref_07, forref_07,               
                    fracrefa_08, fracrefb_08, absa_08, ka_08, absb_08, kb_08, ka_mco2_08, ka_mn2o_08, ka_mo3_08, kb_mco2_08, kb_mn2o_08,    
                    selfref_08, forref_08, cfc12_08, cfc22adj_08, fracrefa_09, fracrefb_09, absa_09, ka_09, absb_09, kb_09,      
                    ka_mn2o_09, kb_mn2o_09, selfref_09, forref_09, fracrefa_10, fracrefb_10, absa_10, ka_10, absb_10, kb_10,      
                    selfref_10, forref_10, fracrefa_11, fracrefb_11, absa_11, ka_11, absb_11, kb_11,      
                    ka_mo2_11, kb_mo2_11, selfref_11, forref_11, fracrefa_12, absa_12, ka_12, selfref_12, forref_12,           
                    fracrefa_13, fracrefb_13, absa_13, ka_13, ka_mco2_13, ka_mco_13, kb_mo3_13, selfref_13, forref_13,                 
                    fracrefa_14, fracrefb_14, absa_14, ka_14, absb_14, kb_14, selfref_14, forref_14,                                     
                    fracrefa_15, absa_15, ka_15, ka_mn2_15, selfref_15, forref_15, fracrefa_16, fracrefb_16, absa_16,ka_16, absb_16,kb_16,      
                    selfref_16, forref_16,fracs,taug,taut,cldfmc,ciwpmc,clwpmc,cswpmc,taucmc,wkl,wx,taua,planklay,planklev,efclfrac,odcld,m,j);
                    //bbugas,bbutot,secdiff,icldlyr,uclfl,dclfl,clrurad,clrdrad,atot,atrans,m,j);
                glw(i,j)  = dflx(1,1,m,j);
                olr(i,j)  = uflx(1,nlayers+1,m,j);
                lwcf(i,j) = uflxc(1,nlayers+1,m,j) - uflx(1,nlayers+1,m,j);
                if(lwupt_p->data != nullptr){
                    //output up and down toa fluxes for total and clear sky:
                    lwupt(i,j)  = uflx(1,nlayers+1,m,j);
                    lwuptc(i,j) = uflxc(1,nlayers+1,m,j);
                    lwdnt(i,j)  = dflx(1,nlayers+1,m,j);
                    lwdntc(i,j) = dflxc(1,nlayers+1,m,j);
                    //output up and down surface fluxes for total and clear sky:
                    lwupb(i,j)  = uflx(1,1,m,j);
                    lwupbc(i,j) = uflxc(1,1,m,j);
                    lwdnb(i,j)  = dflx(1,1,m,j);
                    lwdnbc(i,j) = dflxc(1,1,m,j);
                }
                if(lwupflx_p->data != nullptr){
                    //output up and down fluxes:
                    for(int k = kts; k <= nlayers+1; k++){
                        lwupflx(i,k,j)  = uflx(1,k,m,j);
                        lwupflxc(i,k,j) = uflxc(1,k,m,j);
                        lwdnflx(i,k,j)  = dflx(1,k,m,j);
                        lwdnflxc(i,k,j) = dflxc(1,k,m,j);
                    }
                }
                // --- output heating rate tendency. convert heating rate from K per day to K per second:
                for(int k = kts; k <= kte; k++){
                    tten1d(k,m,j) = hr(ncol,k,m,j)/86400.;
                    rthratenlw(i,k,j) = tten1d(k,m,j)/pi3d(i,k,j);
                }
                // --- output the effective radii for cloud water, cloud ice, and snow:
                if(icloud > 0 && (rre_cloud_p->data != nullptr) && (rre_ice_p->data != nullptr) && (rre_snow_p->data != nullptr) ){
                    for(int k = kts; k <= kte; k++){
                        if(cldfra3d(i,k,j) > 0.){
                            rre_cloud(i,k,j) = rel(ncol,k,m,j);
                            rre_ice(i,k,j)   = rei(ncol,k,m,j);
                            rre_snow(i,k,j)  = res(ncol,k,m,j);
                        }
                    }
                }
        });
    }
    //    }
//    }
    Kokkos::fence();
}
KOKKOS_FUNCTION void inirad(const real3d &o3prof,const real3d &plev,int kts,int kte,int i,int j){
    for(int k = kts; k<= kte+1; k++){
        o3prof(k,i,j) = 0.;
    }
    o3data(o3prof,plev,kts,kte,i,j);
}
KOKKOS_FUNCTION void o3data(const real3d &o3prof,const real3d &plev,int kts,int kte,int i, int j){
    real prlevh[SIZE_NLAY+3];
    real ppwrkh[33];
    real o3wrk[32];
    real ppwrk[32];
    real o3sum[32];
    real ppsum[32];
    real o3win[32];
    real ppwin[32];
    real o3ann[32];
    real ppann[32];
    real pb1,pb2,pt1,pt2;
    o3sum[1] = 5.297E-8;  o3sum[2] = 5.852E-8;  o3sum[3] = 6.579E-8;  o3sum[4] = 7.505E-8;
    o3sum[5] = 8.577E-8;  o3sum[6] = 9.895E-8;  o3sum[7] = 1.175E-7;  o3sum[8] = 1.399E-7;
    o3sum[9] = 1.677E-7;  o3sum[10] = 2.003E-7; o3sum[11] = 2.571E-7; o3sum[12] = 3.325E-7;
    o3sum[13] = 4.438E-7; o3sum[14] = 6.255E-7; o3sum[15] = 8.168E-7; o3sum[16] = 1.036E-6;
    o3sum[17] = 1.366E-6; o3sum[18] = 1.855E-6; o3sum[19] = 2.514E-6; o3sum[20] = 3.240E-6;
    o3sum[21] = 4.033E-6; o3sum[22] = 4.854E-6; o3sum[23] = 5.517E-6; o3sum[24] = 6.089E-6;
    o3sum[25] = 6.689E-6; o3sum[26] = 1.106E-5; o3sum[27] = 1.462E-5; o3sum[28] = 1.321E-5;
    o3sum[29] = 9.856E-6; o3sum[30] = 5.960E-6; o3sum[31] = 5.960E-6; 
    ppsum[1] = 955.890;  ppsum[2] = 850.532;  ppsum[3] = 754.599;  ppsum[4] = 667.742;
    ppsum[5] = 589.841;  ppsum[6] = 519.421;  ppsum[7] = 455.480;  ppsum[8] = 398.085;
    ppsum[9] = 347.171;  ppsum[10] = 301.735; ppsum[11] = 261.310; ppsum[12] = 225.360;
    ppsum[13] = 193.419; ppsum[14] = 165.490; ppsum[15] = 141.032; ppsum[16] = 120.125;
    ppsum[17] = 102.689; ppsum[18] = 87.829;  ppsum[19] = 75.123;  ppsum[20] = 64.306;
    ppsum[21] = 55.086;  ppsum[22] = 47.209;  ppsum[23] = 40.535;  ppsum[24] = 34.795;
    ppsum[25] = 29.865;  ppsum[26] = 19.122;  ppsum[27] = 9.277;   ppsum[28] = 4.660;
    ppsum[29] = 2.421;   ppsum[30] = 1.294;   ppsum[31] = 0.647;
    o3win[1] = 4.629E-8; o3win[2] = 4.686E-8; o3win[3] = 5.017E-8; o3win[4] = 5.613E-8;
    o3win[5] = 6.871E-8; o3win[6] = 8.751E-8; o3win[7] = 1.138E-7; o3win[8] = 1.516E-7;
    o3win[9] = 2.161E-7; o3win[10] = 3.264E-7; o3win[11] = 4.968E-7; o3win[12] = 7.338E-7;
    o3win[13] = 1.017E-6; o3win[14] = 1.308E-6; o3win[15] = 1.625E-6; o3win[16] = 2.011E-6;
    o3win[17] = 2.516E-6; o3win[18] = 3.130E-6; o3win[19] = 3.840E-6; o3win[20] = 4.703E-6;
    o3win[21] = 5.486E-6; o3win[22] = 6.289E-6; o3win[23] = 6.993E-6; o3win[24] = 7.494E-6;
    o3win[25] = 8.197E-6; o3win[26] = 9.632E-6; o3win[27] = 1.113E-5; o3win[28] = 1.146E-5;
    o3win[29] = 9.389E-6; o3win[30] = 6.135E-6; o3win[31] = 6.135E-6;
    ppwin[1] = 955.747; ppwin[2] = 841.783; ppwin[3] = 740.199; ppwin[4] = 649.538;
    ppwin[5] = 568.404; ppwin[6] = 495.815; ppwin[7] = 431.069; ppwin[8] = 373.464;
    ppwin[9] = 322.354; ppwin[10] = 277.190; ppwin[11] = 237.635; ppwin[12] = 203.433;
    ppwin[13] = 174.070; ppwin[14] = 148.949; ppwin[15] = 127.408; ppwin[16] = 108.915;
    ppwin[17] = 93.114;  ppwin[18] = 79.551;  ppwin[19] = 67.940;  ppwin[20] = 58.072;
    ppwin[21] = 49.593;  ppwin[22] = 42.318;  ppwin[23] = 36.138;  ppwin[24] = 30.907;
    ppwin[25] = 26.362;  ppwin[26] = 16.423;  ppwin[27] = 7.583;   ppwin[28] = 3.620;
    ppwin[29] = 1.807;   ppwin[30] = 0.938;   ppwin[31] = 0.469;
    for(int k = 1; k <= 31; k++){
        ppann[k] = ppsum[k];
    }
    o3ann[1] = 0.5 * (o3sum[1] + o3win[1]);
    for(int k = 2; k <= 31; k++){
        o3ann[k] = o3win[k-1] + (o3win[k] - o3win[k-1]) / (ppwin[k] - ppwin[k-1]) * (ppsum[k] - ppwin[k-1]);
    }
    for(int k = 2; k <= 31; k++){
        o3ann[k] = 0.5 * (o3ann[k] + o3sum[k]);
    }
    for(int k = 1; k <= 31; k++){
        o3wrk[k] = o3ann[k];
        ppwrk[k] = ppann[k];
    }
    //  CALCULATE HALF PRESSURE LEVELS FOR MODEL AND DATA LEVELS                     
    //                                                                                

    // Plev is total P at model levels, from bottom to top
    // Plev is in mb
    for(int k = kts; k <= kte+2; k++){
        prlevh[k] = plev(k,i,j);
    }

    ppwrkh[1] = 1100.;
    for(int k = 2; k <= 31; k++){
        ppwrkh[k] = (ppwrk[k] + ppwrk[k-1]) / 2.;
    }
    ppwrkh[32] = 0.;
    for(int k = kts; k <= kte + 1;k++){
        for(int jj = 1; jj <= 31; jj++){
            if((-(prlevh[k] - ppwrkh[jj])) >= 0.){
                pb1 = 0.;
            }else{
                pb1 = prlevh[k] - ppwrkh[jj];
            }
            if((-(prlevh[k] - ppwrkh[jj + 1])) >= 0.){
                pb2 = 0.;
            }else{
                pb2 = prlevh[k] - ppwrkh[jj + 1];
            }
            if((-(prlevh[k + 1] - ppwrkh[jj])) >= 0.){
                pt1 = 0.;
            }else{
                pt1 = prlevh[k + 1] - ppwrkh[jj];
            }
            if((-(prlevh[k + 1] - ppwrkh[jj + 1])) >= 0.){
                pt2 = 0.;
            }else{
                pt2 = prlevh[k + 1] - ppwrkh[jj + 1];
            }
            o3prof(k,i,j) = o3prof(k,i,j) + (pb2 - pb1 - pt2 + pt1) * o3wrk[jj];
        }
        o3prof(k,i,j) = o3prof(k,i,j) / (prlevh[k] - prlevh[k + 1]);
    }
}

KOKKOS_FUNCTION string trim(char *config_chem_opt,string &schemes){
    int length = 80; //-whitespace character
    while (length > 0 && config_chem_opt[length-1] == ' ') {
        length--;
    }
    // Terminate the string at the last non-whitespace character
    int start = 0;
    while (start < length - 1 && config_chem_opt[start] == ' ') {
        start++;
    }
    for(int i = start; i < length; i++){
        schemes += config_chem_opt[i];
    }
    return schemes;
}

KOKKOS_FUNCTION void vinterp_ozn(int lchnk, int ncol, int pcols, int pver, const real4d &pmid_in, const real1d_u &pin_in, int levsiz, const real3d &ozmix_in, const real3d &o3vmr_out,int m,int j){
    /*
        subroutine vinterp_ozn is called from the subroutines rrtmg_lwrad  and rrtmg_swrad 
        for interpolation of climatological ozone onto the MPAS pressure levels.
        Purpose: Interpolate ozone from current time-interpolated values to model levels  
        Method: Use pressure values to determine interpolation levels
    */
    int kupper[SIZE_NCOL+1]; // Level indices for interpolation
    int lats[SIZE_NCOL+1]; // latitude indices
    int lons[SIZE_NCOL+1]; // lontitude indices
    real dpu; // upper level pressure difference
    real dpl; // lower level pressure difference
    real pin[SIZE_LEVSIZ+1]; // ozone data level pressures (mks)
    real ozmix[SIZE_LEVSIZ+1]; // ozone mixing ratio
    real o3vmr[SIZE_NLAY+1];   // ozone volume mixing ratio
    real pmid[SIZE_NLAY+1][SIZE_NCOL+1];      // level pressures (mks)
           
    int kkstart;
    int kcount;
    int kk;
    for(int k = 1; k <= levsiz; k++){
        pin[k] = pin_in(k);
    }
    for(int k = 1; k <= levsiz; k++){
        ozmix[k] = ozmix_in(k,m,j);
    }
    for(int i = 1; i <= pcols; i++){
        for(int k = 1; k <= pver; k++){
            kk = pver - k + 1;
            pmid[kk][i] = pmid_in(i,k,m,j);
        }
    }
    // Initialize index array
    for(int i = 1; i <= ncol; i++){
        kupper[i] = 1;
    }
    bool jump = false;
    for(int k = 1; k <= pver; k++){
    // Top level we need to start looking is the top level for the previous k
    // for all longitude points
        kkstart = levsiz;
        for(int i = 1; i <= ncol; i++){
            kkstart = std::min(kkstart,kupper[i]);
        }
        kcount = 0;
    // Store level indices for interpolation
        for(int kk = kkstart; kk <= levsiz - 1; kk++){
            for(int i = 1; i <= ncol; i++){
                if (pin[kk] < pmid[k][i] && pmid[k][i] <= pin[kk+1]){
                    kupper[i] = kk;
                    kcount = kcount + 1;
                }
            }
            // If all indices for this level have been found, do the interpolation and
            // go to the next level
            if (kcount == ncol){
                for(int i = 1; i <= ncol; i++){
                    dpu = pmid[k][i] - pin[kupper[i]];
                    dpl = pin[kupper[i]+1] - pmid[k][i];
                    o3vmr[k] = (ozmix[kupper[i]]*dpl + ozmix[kupper[i]+1]*dpu)/(dpl + dpu);
                }    
                jump = true;
				break;   
            }   
        }
        if(jump){
            jump = false;
            continue;
        }
    // If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
    // must extrapolate from the bottom or top ozone data level for at least some
    // of the longitude points.
        for(int i = 1; i <= ncol; i++){
            if(pmid[k][i] < pin[1]){
                o3vmr[k] = ozmix[1]*pmid[k][i]/pin[1];
            }else if(pmid[k][i] > pin[levsiz]){
                o3vmr[k] = ozmix[levsiz];
            }else{
                dpu = pmid[k][i] - pin[kupper[i]];
                dpl = pin[kupper[i]+1] - pmid[k][i];
                o3vmr[k] = (ozmix[kupper[i]]*dpl + 
                            ozmix[kupper[i]+1]*dpu)/(dpl + dpu);
            }
        }
        if(kcount > ncol){
            //        call endrun ('VINTERP_OZN: Bad ozone data: non-monotonicity suspected')
        }
    }
    for(int k = 1; k <= pver; k++){
        kk = pver-k+1;
        o3vmr_out(kk,m,j) = o3vmr[k];
    }
}

KOKKOS_FUNCTION void relcalc(int ncol,int pcols,int pver, const real4d &t, const real3d &landfrac, const real3d &landm, const real3d &icefrac,const real4d &rel, const real3d &snowh,int m, int n){
    // Compute cloud water size
    // analytic formula following the formulation originally developed by J. T. Kiehl
    real tmelt;             // freezing temperature of fresh water (K)
    real rliqland;          // liquid drop size if over land
    real rliqocean;        // liquid drop size if over ocean
    real rliqice;          // liquid drop size if over sea ice
    tmelt = 273.16;
    rliqocean = 14.0;
    rliqice   = 14.0;
    rliqland  = 8.0;
    for(int k = 1; k <= pver; k++){
        for(int i = 1; i <= ncol; i++){
            // jrm Reworked effective radius algorithm
            // Start with temperature-dependent value appropriate for continental air
            // Note: findmcnew has a pressure dependence here
            rel(i,k,m,n) = rliqland + (rliqocean-rliqland) * std::min(1.0,std::max(0.0,(tmelt-t(i,k,m,n))*0.05));
            // Modify for snow depth over land
            rel(i,k,m,n) = rel(i,k,m,n) + (rliqocean-rel(i,k,m,n)) * std::min(1.0,std::max(0.0,snowh(i,m,n)*10.));
            // Ramp between polluted value over land to clean value over ocean.
            rel(i,k,m,n) = rel(i,k,m,n) + (rliqocean-rel(i,k,m,n)) * std::min(1.0,std::max(0.0,1.0-landm(i,m,n)));
            // Ramp between the resultant value and a sea ice value in the presence of ice.
            rel(i,k,m,n) = rel(i,k,m,n) + (rliqice-rel(i,k,m,n)) * std::min(1.0,std::max(0.0,icefrac(i,m,n)));
            // end jrm        
        }
    }
}

KOKKOS_FUNCTION void reicalc(int ncol,int pcols,int pver,const real4d &t, const real4d &re,const real1d_u &retab,int m, int n){
    real corr;
    int index;
    // Tabulated values of re(T) in the temperature interval 
    // 180 K -- 274 K; hexagonal columns assumed:
    for(int k = 1; k <= pver; k++){
        for(int i = 1; i <= ncol; i++){
            index = int(t(i,k,m,n)-179.);
            index = std::min(std::max(index,1),94);
            corr = t(i,k,m,n) - int(t(i,k,m,n));
            re(i,k,m,n) = retab(index)*(1.-corr)+retab(index+1)*corr;
        }
    }
}

KOKKOS_FUNCTION int m(int k,int n){
    unsigned int uk = (unsigned int)k;
    if(n < 0){
        return k ^ (uk >> -n);
    }else{
        return k ^ (uk << n);
    }
}

KOKKOS_FUNCTION void kissvec(int1d_u &seed1, int1d_u &seed2, int1d_u &seed3, int1d_u &seed4, real1d_u &ran_arr,int size_ran){
    /** 
      The  KISS (Keep It Simple Stupid) random number generator. Combines:
      (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
      (2) A 3-shift shift-register generator, period 2^32-1,
      (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
       Overall period>2^123;       
    */  
    int k,n,kiss;
    for(int i = 1; i <= size_ran; i++){
        seed1(i) = 69069 * seed1(i) + 1327217885;
        seed2(i) = m(m(m(seed2(i), 13), -17), 5);
        seed3(i) = 18000 * (seed3(i) & 65535) +  ((unsigned int)seed3(i) >> 16);
        seed4(i) = 30903 * (seed4(i) & 65535) +  ((unsigned int)seed4(i) >> 16);
        kiss = seed1(i) + seed2(i) + ((unsigned int)seed3(i) << 16) + seed4(i);
        ran_arr(i) = kiss*2.328306e-10 + 0.5;
    } 
            
}
KOKKOS_FUNCTION void generate_stochastic_clouds(int& ncol,int& nlay,const int nsubcol,int& icld,int& irng,real2d_u& pmid, const real4d& cld, const real4d& clwp, const real4d& ciwp, const real4d& cswp,
    const real5d& tauc, const real5d& cld_stoch, const real5d& clwp_stoch, const real5d& ciwp_stoch, const real5d& cswp_stoch, const real5d& tauc_stoch,int& changeSeed,int m,int n){
/**
   ----------------------------------------------------------------------------------------------------------------
   ---------------------
   Contact: Cecile Hannay (hannay@ucar.edu)
   
   Original code: Based on Raisanen et al., QJRMS, 2004.
   
   Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
     random number generator, which can be changed to the optional kissvec random number generator
     with flag 'irng'. Some extra functionality has been commented or removed.  
     Michael J. Iacono, AER, Inc., February 2007
  
   Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
   Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
   and uniform cloud liquid and cloud ice concentration.
   The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
   and obeys an overlap assumption in the vertical.   
   
   Overlap assumption:
    The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
    The default option is maximum-random (option 3)
    The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
    This is set with the variable "overlap" 
    mji - Exponential overlap option (overlap=4) has been deactivated in this version
    The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
   
   Seed:
    If the stochastic cloud generator is called several times during the same timestep, 
    one should change the seed between the call to insure that the subcolumns are different.
    This is done by changing the argument 'changeSeed'
    For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
    use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
  
   PDF assumption:
    We can use arbitrary complicated PDFS. 
    In the present version, we produce homogeneuous clouds (the simplest case).  
    Future developments include using the PDF scheme of Ben Johnson. 
  
   History file:
    Option to add diagnostics variables in the history file. (using FINCL in the namelist)
    nsubcol = number of subcolumns
    overlap = overlap type (1-3)
    Zo = length scale 
    CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
    CLDLIQ_S = mean of the subcolumn cloud water
    CLDICE_S = mean of the subcolumn cloud ice 
  
   Note:
     Here: we force that the cloud condensate to be consistent with the cloud fraction 
     i.e we only have cloud condensate when the cell is cloudy. 
     In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
     and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
     without cloud condensate or the opposite).
  ---------------------------------------------------------------------------------------------------------------
*/  
    real cldf_s[SIZE_NCOL*SIZE_NLAY];
    real2d_u cldf((real*)cldf_s,{1,1},{ncol+1,nlay+1});
    int overlap;                                    // 1 = random overlap, 2 = maximum/random,
                                                // 3 = maximum overlap, 
    std::vector<real> cldf_spread_s(ncol * nlay * nsubcol);
    real3d_u cldf_spread(cldf_spread_s.data(),{1,1,1},{nsubcol+1,ncol+1,nlay+1});
    std::vector<int> iscloudy_s(ncol * nlay * nsubcol);
    int3d_u iscloudy(iscloudy_s.data(),{1,1,1},{nsubcol+1,ncol+1,nlay+1});
    std::vector<real> CDF_s(ncol * nlay * nsubcol);
    real3d_u CDF(CDF_s.data(),{1,1,1},{nsubcol+1,ncol+1,nlay+1});
    const real cldmin = 1.0e-20;                    //min cloud fraction

    int seed1_s[SIZE_NCOL+1];
    int1d_u seed1((int*)seed1_s,{1},{ncol+1});                //seed to create random number
    int seed2_s[SIZE_NCOL+1];
    int1d_u seed2((int*)seed2_s,{1},{ncol+1});
    int seed3_s[SIZE_NCOL+1];
    int1d_u seed3((int*)seed3_s,{1},{ncol+1});
    int seed4_s[SIZE_NCOL+1];
    int1d_u seed4((int*)seed4_s,{1},{ncol+1});

    real rand_num_s[SIZE_NCOL+1];
    real1d_u rand_num((real*)rand_num_s,{1},{ncol+1});         //random number (kissvec)
    int iseed;                                              //seed to create random number (Mersenne Twister)

    real rand_num_mt;                                       //random number (Mersenne Twister)

    int ilev, isubcol, i;       //indices

    // ------------------------------------------------------------------------------------------ 

    // Check that irng is in bounds; if not, set to default

    if(irng != 0) irng=1;
    overlap = icld;

    //Ensure that cloud fractions are in bounds
    for(int ilev = 1;ilev <= nlay;ilev++){
        for(int i = 1;i <= ncol;i++){
            cldf(i,ilev) = cld(i,ilev,m,n);
            if(cldf(i,ilev) < cldmin) cldf(i,ilev) = 0.0;
        }
    }
    // ----- Create seed  --------

    // Advance randum number generator by changeseed values
    if(irng == 0){
        // For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
        // Must use pmid from bottom four layers.
        for(int i = 1;i<=ncol;i++){
            if(pmid(i,1) < pmid(i,2)) throw std::runtime_error("MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.");
            seed1(i) = (pmid(i,1) - int(pmid(i,1)))  * 1000000000;
            seed2(i) = (pmid(i,2) - int(pmid(i,2)))  * 1000000000;
            seed3(i) = (pmid(i,3) - int(pmid(i,3)))  * 1000000000;
            seed4(i) = (pmid(i,4) - int(pmid(i,4)))  * 1000000000;
        }
        
        for(int i = 1;i <= changeSeed;i++){
            kissvec(seed1, seed2, seed3, seed4, rand_num, ncol);
        }
    }else if(irng == 1){
        init_genrand(changeSeed);
    }


    // ------ Apply overlap assumption --------

    // generate the random numbers

    switch (overlap){
    case 1:
    /** 
    Random overlap
    i) pick a random value at every level
    */
    
        if(irng == 0){
            for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
                for(int ilev = 1;ilev <= nlay;ilev++){
                    kissvec(seed1, seed2, seed3, seed4, rand_num, ncol);
                    for(int i = 1;i<=ncol;i++) CDF(isubcol,i,ilev) = rand_num(i);
                }
            }
        }else if(irng == 1){
            for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
                for(int i = 1;i <= ncol;i++){
                    for(int ilev = 1;ilev <= nlay;ilev++){
                        rand_num_mt = genrand_real1();
                        CDF(isubcol,i,ilev) = rand_num_mt;
                    }
                }
            }
        }
    
    break;

    case 2:
    /** 
    Maximum-Random overlap
    i) pick a random number for top layer.
    ii) walk down the column: 
    - if the layer above is cloudy, we use the same random number than in the layer above
    - if the layer above is clear, we use a new random number 
    */          
        if(irng == 0){
            for(int ilev = 1;ilev <= nlay;ilev++){
                for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
                   kissvec(seed1, seed2, seed3, seed4, rand_num, ncol);
                   for(int i = 1;i <= ncol;i++) CDF(isubcol,i,ilev) = rand_num(i);
                }
            }
        }else if(irng == 1){
            for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
                for(int i = 1;i <= ncol;i++){
                    for(int ilev = 1;ilev <= nlay;ilev++){
                        rand_num_mt = genrand_real1();
                        CDF(isubcol,i,ilev) = rand_num_mt;
                    }
                }
            }

        }
       for(int ilev = 2;ilev <= nlay;ilev++){
           for(int i = 1;i <= ncol;i++){
               for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
                   if (CDF(isubcol, i, ilev-1) > 1.0 - cldf(i,ilev-1) ){
                       CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev-1); 
                   }else{
                       CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev) * (1. - cldf(i,ilev-1)); 
                   }
               }
           }
       }
        
    break;

    case 3:
    /** 
    Maximum overlap
    i) pick the same random numebr at every leve
    */
    if (irng == 0){
        for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
            kissvec(seed1, seed2, seed3, seed4, rand_num, ncol);
            for(int ilev = 1;ilev <= nlay;ilev++){
                    for(int i = 1;i <= ncol;i++) CDF(isubcol,i,ilev) = rand_num(i);
            }
        }
    }else if (irng == 1){
        for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
            for(int i = 1;i <= ncol;i++){
                rand_num_mt = genrand_real1();
                for(int ilev = 1;ilev <= nlay;ilev++){
                    CDF(isubcol,i,ilev) = rand_num_mt;
                }
            }
        }
    }                                    
    break;
   }
// /** 
// where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
// where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
// where there is a cloud, define the subcolumn cloud properties, 
// otherwise set these to zero
// */
        // auto &ngb = lw_rad.ngb;
        // memset(cld_stoch.data(),0,cld_stoch.size()*sizeof(real));
        // memset(clwp_stoch.data(),0,clwp_stoch.size()*sizeof(real));
        // memset(ciwp_stoch.data(),0,ciwp_stoch.size()*sizeof(real));
        // memset(cswp_stoch.data(),0,cswp_stoch.size()*sizeof(real));
        // memset(tauc_stoch.data(),0,tauc_stoch.size()*sizeof(real));
        for(int ilev = 1; ilev <= nlay; ilev++){
	        for(int i = 1;i <= ncol;i++){
		        for(int isubcol = 1;isubcol <= nsubcol;isubcol++){
                    cldf_spread(isubcol,i,ilev) = cldf(i,ilev);
                    iscloudy(isubcol,i,ilev) = CDF(isubcol,i,ilev) >= (1.0 - cldf_spread(isubcol,i,ilev)) ? 1 : 0;
                    if ( iscloudy(isubcol,i,ilev) == 1 ){
                        cld_stoch(isubcol,i,ilev,m,n) = 1.0;
                        clwp_stoch(isubcol,i,ilev,m,n) = clwp(i,ilev,m,n);
                        ciwp_stoch(isubcol,i,ilev,m,n) = ciwp(i,ilev,m,n);
                        cswp_stoch(isubcol,i,ilev,m,n) = cswp(i,ilev,m,n);
                        int nb = lw_rad.ngb[isubcol - 1];
                        tauc_stoch(isubcol,i,ilev,m,n) = tauc(nb,i,ilev,m,n);
                    }else{
                        cld_stoch(isubcol,i,ilev,m,n) = 0;
                        clwp_stoch(isubcol,i,ilev,m,n) = 0.;
                        ciwp_stoch(isubcol,i,ilev,m,n) = 0.;
                        cswp_stoch(isubcol,i,ilev,m,n) = 0.;
                        tauc_stoch(isubcol,i,ilev,m,n) = 0.;
                    }
                }
            }
        }

}
KOKKOS_FUNCTION void mcica_subcol_lw(int& iplon,int& ncol,int& nlay,int& icld,int& permuteseed,int& irng, const real4d &play,
    const real4d &cldfrac, const real4d &ciwp, const real4d &clwp, const real4d &cswp, const real4d &rei, const real4d &rel, const real4d &res, const real5d &tauc,
    const real5d &cldfmcl, const real5d &ciwpmcl, const real5d &clwpmcl, const real5d &cswpmcl, const real4d &reicmcl, const real4d &relqmcl, const real4d &resnmcl, const real5d &taucmcl,
    int m , int n){
    const int nsubclw = ngptlw;
    int ilev;
    real pmid_s[SIZE_NCOL*SIZE_NLAY];
    real2d_u pmid((real*)pmid_s,{1,1},{ncol+1,nlay+1});
    if(icld == 0) return;
    if(icld < 0 || icld >3){
            throw std::runtime_error("MCICA_SUBCOL: INVALID ICLD");
    }
    for(int lay = 1;lay<=nlay;lay++){
        for(int col = 1;col<=ncol;col++){
            reicmcl(col,lay,m,n) = rei(col,lay,m,n);
            relqmcl(col,lay,m,n) = rel(col,lay,m,n);
            resnmcl(col,lay,m,n) = res(col,lay,m,n);
            pmid(col,lay) = play(col,lay,m,n)*1.e2;
        }
    }
    generate_stochastic_clouds(ncol, nlay, nsubclw, icld, irng, pmid, cldfrac, clwp, ciwp, cswp,
		tauc, cldfmcl, clwpmcl, ciwpmcl, cswpmcl,taucmcl, permuteseed,m,n);
}
KOKKOS_FUNCTION void inatm(int &iplon, int &nlay,int &icld,int &iaer, const real4d &play, const real4d &plev, const real4d &tlay, const real4d &tlev, const real3d &tsfc, const real4d &h2ovmr, const real4d &o3vmr,
    const real4d &co2vmr, const real4d &ch4vmr, const real4d &n2ovmr, const real4d &o2vmr, const real4d &cfc11vmr, const real4d &cfc12vmr, const real4d &cfc22vmr, const real4d &ccl4vmr,
    const real4d &emis,int &inflglw,int &iceflglw,int &liqflglw, const real5d &cldfmcl, const real5d &taucmcl, const real5d &ciwpmcl, const real5d &clwpmcl,
    const real5d &cswpmcl, const real4d &reicmcl, const real4d &relqmcl, const real4d &resnmcl, const real5d &tauaer,int &nlayers, real1d_u &pavel, real1d_u &pz, real1d_u &tavel, real1d_u &tz,
    real &tbound, real1d_u &semiss, real1d_u &coldry, const real4d &wkl, real1d_u &wbrodl, const real4d &wx,real &pwvcm,int &inflag,int &iceflag,int &liqflag,
    const real4d &cldfmc, const real4d &taucmc, const real4d &ciwpmc, const real4d &clwpmc, const real4d &cswpmc, real1d_u &reicmc, real1d_u &relqmc, real1d_u &resnmc, const real4d &taua,int m,int n){
    /*
        !***************************************************************************
        !
        !  Input atmospheric profile from GCM, and prepare it for use in RRTMG_LW.
        !  Set other RRTMG_LW input parameters.  
        !
        !***************************************************************************
    */
    const real amd = 28.9660; // Effective molecular weight of dry air (g/mol)
    const real amw = 18.0160; // Molecular weight of water vapor (g/mol)
    // Set molecular weight ratios (for converting mmr to vmr)
    const real amdw = 1.607793;  // Molecular weight of dry air / water vapor
    const real amdc = 0.658114;  // Molecular weight of dry air / carbon dioxide
    const real amdo = 0.603428;  // Molecular weight of dry air / ozone
    const real amdm = 1.805423;  // Molecular weight of dry air / methane
    const real amdn = 0.658090;  // Molecular weight of dry air / nitrous oxide
    const real amdo2 = 0.905140; // Molecular weight of dry air / oxygen
    const real amdc1 = 0.210852; // Molecular weight of dry air / CFC11
    const real amdc2 = 0.239546; // Molecular weight of dry air / CFC12
    const real grav = 9.8066;
    const real avogad = 6.02214199e+23;
    real amm,amttl,wvttl,wvsh,summol;
    // Add one to nlayers here to include extra model layer at top of atmosphere
    nlayers = nlay;
    // Initialize all molecular amounts and cloud properties to zero here, then pass input amounts
    // into RRTM arrays below.
    for(int l = 1; l <= nlayers + 1; l++){
        for(int i = 1; i <= mxmol; i++){
            wkl(i,l,m,n) = 0.0;
        }
        for(int ix = 1; ix <= maxxsec; ix++){
            wx(ix,l,m,n) = 0.0;
        }
        for(int ig = 1; ig <= ngptlw; ig++){
            cldfmc(ig,l,m,n) = 0.0;
            taucmc(ig,l,m,n) = 0.0;
            ciwpmc(ig,l,m,n) = 0.0;
            clwpmc(ig,l,m,n) = 0.0;
            cswpmc(ig,l,m,n) = 0.0;
        }
        for(int ib = 1; ib <= nbndlw; ib++){
            taua(ib,l,m,n) = 0.0;
        }
    }

    amttl = 0.0;
    wvttl = 0.0;
    // Set surface temperature.
    tbound = tsfc(iplon,m,n);
    /*
        !  Install input GCM arrays into RRTMG_LW arrays for pressure, temperature,
        !  and molecular amounts.  
        !  Pressures are input in mb, or are converted to mb here.
        !  Molecular amounts are input in volume mixing ratio, or are converted from 
        !  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
        !  here. These are then converted to molecular amount (molec/cm2) below.  
        !  The dry air column COLDRY (in molec/cm2) is calculated from the level 
        !  pressures, pz (in mb), based on the hydrostatic equation and includes a 
        !  correction to account for h2o in the layer.  The molecular weight of moist 
        !  air (amm) is calculated for each layer.  
        !  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
        !  assumes GCM input fields are also bottom to top. Input layer indexing
        !  from GCM fields should be reversed here if necessary.
    */
    pz(0) = plev(iplon,1,m,n);
    tz(0) = tlev(iplon,1,m,n);
    for(int l = 1; l <= nlayers; l++){
        pavel(l) = play(iplon,l,m,n);
        tavel(l) = tlay(iplon,l,m,n);
        pz(l) = plev(iplon,l+1,m,n);
        tz(l) = tlev(iplon,l+1,m,n);
        // For h2o input in vmr:
        wkl(1,l,m,n) = h2ovmr(iplon,l,m,n);
        // For h2o input in mmr:
        // wkl(1,l) = h2o(iplon,l)*amdw
        // For h2o input in specific humidity;
        // wkl(1,l) = (h2o(iplon,l)/(1._rb - h2o(iplon,l)))*amdw
        wkl(2,l,m,n) = co2vmr(iplon,l,m,n);
        wkl(3,l,m,n) = o3vmr(iplon,l,m,n);
        wkl(4,l,m,n) = n2ovmr(iplon,l,m,n);
        wkl(6,l,m,n) = ch4vmr(iplon,l,m,n);
        wkl(7,l,m,n) = o2vmr(iplon,l,m,n);
        amm = (1. - wkl(1,l,m,n)) * amd + wkl(1,l,m,n) * amw;        
        coldry(l) = (pz(l-1)-pz(l)) * 1.e3 * avogad / (1.e2 * grav * amm * (1. + wkl(1,l,m,n)));
    }
    // Set cross section molecule amounts from input; convert to vmr if necessary
    for(int l = 1; l <= nlayers; l++){
        wx(1,l,m,n) = ccl4vmr(iplon,l,m,n);
        wx(2,l,m,n) = cfc11vmr(iplon,l,m,n);
        wx(3,l,m,n) = cfc12vmr(iplon,l,m,n);
        wx(4,l,m,n) = cfc22vmr(iplon,l,m,n);
    }
    /*
        ! The following section can be used to set values for an additional layer (from
        ! the GCM top level to 1.e-4 mb) for improved calculation of TOA fluxes. 
        ! Temperature and molecular amounts in the extra model layer are set to 
        ! their values in the top GCM model layer, though these can be modified
        ! here if necessary. 
        ! If this feature is utilized, increase nlayers by one above, limit the two
        ! loops above to (nlayers-1), and set the top most (extra) layer values here. 
    */
    /*
        !      pavel(nlayers) = 0.5_rb * pz(nlayers-1)
        !      tavel(nlayers) = tavel(nlayers-1)
        !      pz(nlayers) = 1.e-4_rb
        !      tz(nlayers-1) = 0.5_rb * (tavel(nlayers)+tavel(nlayers-1))
        !      tz(nlayers) = tz(nlayers-1)
        !      wkl(1,nlayers) = wkl(1,nlayers-1)
        !      wkl(2,nlayers) = wkl(2,nlayers-1)
        !      wkl(3,nlayers) = wkl(3,nlayers-1)
        !      wkl(4,nlayers) = wkl(4,nlayers-1)
        !      wkl(6,nlayers) = wkl(6,nlayers-1)
        !      wkl(7,nlayers) = wkl(7,nlayers-1)
        !      amm = (1._rb - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
        !      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_rb * avogad / &
        !                        (1.e2_rb * grav * amm * (1._rb + wkl(1,nlayers-1)))
        !      wx(1,nlayers) = wx(1,nlayers-1)
        !      wx(2,nlayers) = wx(2,nlayers-1)
        !      wx(3,nlayers) = wx(3,nlayers-1)
        !      wx(4,nlayers) = wx(4,nlayers-1)
    */
    /*
        ! At this point all molecular amounts in wkl and wx are in volume mixing ratio; 
        ! convert to molec/cm2 based on coldry for use in rrtm.  also, compute precipitable
        ! water vapor for diffusivity angle adjustments in rtrn and rtrnmr.
    */
    for(int l = 1; l <= nlayers; l++){
        summol = 0.0;
        for(int imol = 2; imol <= nmol; imol++){
            summol = summol + wkl(imol,l,m,n);
        }
        wbrodl(l) = coldry(l) * (1. - summol);
        for(int imol = 1; imol <= nmol; imol++){
            wkl(imol,l,m,n) = coldry(l) * wkl(imol,l,m,n);
        }
        amttl = amttl + coldry(l)+wkl(1,l,m,n);
        wvttl = wvttl + wkl(1,l,m,n);
        for(int ix = 1; ix <= maxxsec; ix++){
            if (lw_rad.ixindx[ix - 1] != 0){
               wx(lw_rad.ixindx[ix - 1],l,m,n) = coldry(l) * wx(ix,l,m,n) * 1.e-20;
            }
        }
    }
    wvsh = (amw * wvttl) / (amd * amttl);
    pwvcm = wvsh * (1.e3 * pz(0)) / (1.e2 * grav);
    // Set spectral surface emissivity for each longwave band.  
    for(int nb = 1; nb <= nbndlw; nb++){
        semiss(nb) = emis(iplon,nb,m,n);
    }
    // Transfer aerosol optical properties to RRTM variable;
    // modify to reverse layer indexing here if necessary.
    if(iaer >= 1){
        for(int l = 1; l <= nlayers; l++){
            for(int ib = 1; ib <= nbndlw; ib++){
                taua(ib,l,m,n) = tauaer(iplon,l,ib,m,n);
            }
        }
    }
    // Transfer cloud fraction and cloud optical properties to RRTM variables,
    // modify to reverse layer indexing here if necessary.
    if(icld >= 1){
        inflag = inflglw;
        iceflag = iceflglw;
        liqflag = liqflglw;
        // Move incoming GCM cloud arrays to RRTMG cloud arrays.
        // For GCM input, incoming reicmcl is defined based on selected ice parameterization (inflglw)
        for(int l = 1; l <= nlayers; l++){
            for(int ig = 1; ig <= ngptlw; ig++){
                cldfmc(ig,l,m,n) = cldfmcl(ig,iplon,l,m,n);
                taucmc(ig,l,m,n) = taucmcl(ig,iplon,l,m,n);
                ciwpmc(ig,l,m,n) = ciwpmcl(ig,iplon,l,m,n);
                clwpmc(ig,l,m,n) = clwpmcl(ig,iplon,l,m,n);
                cswpmc(ig,l,m,n) = cswpmcl(ig,iplon,l,m,n);  
            }
            reicmc(l) = reicmcl(iplon,l,m,n);
            relqmc(l) = relqmcl(iplon,l,m,n);
            resnmc(l) = resnmcl(iplon,l,m,n);
        }
        // If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.
        // !         cldfmc(:,nlayers) = 0.0_rb
        // !         taucmc(:,nlayers) = 0.0_rb
        // !         ciwpmc(:,nlayers) = 0.0_rb
        // !         clwpmc(:,nlayers) = 0.0_rb
        // !         reicmc(nlayers) = 0.0_rb
        // !         relqmc(nlayers) = 0.0_rb
        // !         taua(nlayers,:) = 0.0_rb
    }
}
KOKKOS_FUNCTION void cldprmc(int &nlayers,int &inflag,int &iceflag, int &liqflag, const real4d &cldfmc,const real4d &ciwpmc,const real4d &clwpmc,const real4d &cswpmc,real1d_u &reicmc,real1d_u &relqmc,
    real1d_u &resnmc,int &ncbands,const real4d &taucmc,int m, int n){
    // Purpose:  Compute the cloud optical depth(s) for each cloudy layer.
    int index;
    int ib;
    int icb_s[nbndlw+1];
    int1d_u icb((int*)icb_s,{1},{nbndlw+1});
    real abscoice_s[ngptlw+1];
    real1d_u abscoice((real*)abscoice_s,{1},{ngptlw+1});    // ice absorption coefficients
    real abscoliq_s[ngptlw+1];
    real1d_u abscoliq((real*)abscoliq_s,{1},{ngptlw+1});    // liquid absorption coefficients
    real abscosno_s[ngptlw+1];
    real1d_u abscosno((real*)abscosno_s,{1},{ngptlw+1});    // snow absorption coefficients
    real cwp;                                           // cloud water path
    real radice;                                        // cloud ice effective size (microns)
    real factor;                                        //  
    real fint;                                          //
    real radliq;                                        // cloud liquid droplet radius (microns)
    real radsno;                                        // cloud snow effective size (microns)
    const real eps = 1.e-6;                             // epsilon
    const real cldmin = 1.e-20;                         // minimum value for cloud quantities
    string errmess;
    /*
        ! ------- Definitions -------

        !     Explanation of the method for each value of INFLAG.  Values of
        !     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
        !     INFLAG = 2 does distinguish between liquid and ice clouds, and
        !     requires further user input to specify the method to be used to 
        !     compute the aborption due to each.
        !     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
        !                  optical depth are input.  
        !     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
        !                  water path (g/m2) are input.  The (gray) cloud optical 
        !                  depth is computed as in CCM2.
        !     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
        !                  water path (g/m2), and cloud ice fraction are input.
        !       ICEFLAG = 0:  The ice effective radius (microns) is input and the
        !                     optical depths due to ice clouds are computed as in CCM3.
        !       ICEFLAG = 1:  The ice effective radius (microns) is input and the
        !                     optical depths due to ice clouds are computed as in 
        !                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The 
        !                     spectral regions in this work have been matched with
        !                     the spectral bands in RRTM to as great an extent 
        !                     as possible:  
        !                     E&C 1      IB = 5      RRTM bands 9-16
        !                     E&C 2      IB = 4      RRTM bands 6-8
        !                     E&C 3      IB = 3      RRTM bands 3-5
        !                     E&C 4      IB = 2      RRTM band 2
        !                     E&C 5      IB = 1      RRTM band 1
        !       ICEFLAG = 2:  The ice effective radius (microns) is input and the
        !                     optical properties due to ice clouds are computed from
        !                     the optical properties stored in the RT code,
        !                     STREAMER v3.0 (Reference: Key. J., Streamer 
        !                     User's Guide, Cooperative Institute for
        !                     Meteorological Satellite Studies, 2001, 96 pp.).
        !                     Valid range of values for re are between 5.0 and
        !                     131.0 micron.
        !       ICEFLAG = 3: The ice generalized effective size (dge) is input
        !                    and the optical properties, are calculated as in
        !                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
        !                    tables which were appropriately averaged for the
        !                    bands in RRTM_LW.  Linear interpolation is used to
        !                    get the coefficients from the stored tables.
        !                    Valid range of values for dge are between 5.0 and
        !                    140.0 micron.
        !       LIQFLAG = 0:  The optical depths due to water clouds are computed as
        !                     in CCM3.
        !       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
        !                     and the optical depths due to water clouds are computed 
        !                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
        !                     The values for absorption coefficients appropriate for
        !                     the spectral bands in RRTM have been obtained for a 
        !                     range of effective radii by an averaging procedure 
        !                     based on the work of J. Pinto (private communication).
        !                     Linear interpolation is used to get the absorption 
        !                     coefficients for the input effective radius.
    */
   icb(1) = 1; icb(2) = 2; icb(3) = 3; icb(4) = 3; 
   icb(5) = 3; icb(6) = 4; icb(7) = 4; icb(8) = 4;
   icb(9) = 5; icb(10) = 5; icb(11) = 5; icb(12) = 5;
   icb(13) = 5; icb(14) = 5; icb(15) = 5; icb(16) = 5;
   // jm not thread safe      hvrclc = '$Revision: 1.8 $'
   ncbands = 1;
   /*
        ! This initialization is done in rrtmg_lw_subcol.F90.
        !      do lay = 1, nlayers
        !         do ig = 1, ngptlw
        !            taucmc(ig,lay) = 0.0_rb
        !         enddo
        !      enddo

        ! Main layer loop
    */
    for(int lay = 1; lay <= nlayers; lay++){
        for(int ig = 1; ig <= ngptlw; ig++){
            cwp = ciwpmc(ig,lay,m,n) + clwpmc(ig,lay,m,n) + cswpmc(ig,lay,m,n);
            if ((cldfmc(ig,lay,m,n) >= cldmin) && (cwp >= cldmin || taucmc(ig,lay,m,n) >= cldmin)){
               // Ice clouds and water clouds combined.
                if (inflag == 0) {
                    // Cloud optical depth already defined in taucmc, return to main program
                    return;
                }else if(inflag == 1){
                    throw std::runtime_error("INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA");
                }else if(inflag >= 2){
                    radice = reicmc(lay);
                    if ((ciwpmc(ig,lay,m,n)+cswpmc(ig,lay,m,n)) == 0.0){
                        abscoice(ig) = 0.0;
                        abscosno(ig) = 0.0;
                    }else if(iceflag == 0){
                        if (radice < 10.0) throw std::runtime_error("ICE RADIUS TOO SMALL"); 
                        abscoice(ig) = lw_rad.absice0[0] + lw_rad.absice0[1]/radice;
                        abscosno(ig) = 0.0;
                    }else if(iceflag == 1){
                        if (radice < 13.0 || radice > 130.) throw std::runtime_error("ICE RADIUS OUT OF BOUNDS"); 
                        ncbands = 5;
                        ib = icb(lw_rad.ngb[ig - 1]);
                        abscoice(ig) = lw_rad.absice1[0][ib - 1] + lw_rad.absice1[1][ib - 1]/radice;
                        abscosno(ig) = 0.0;
                    // For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns
                    }else if(iceflag == 2){
                        if (radice < 5.0 || radice > 131.0) throw std::runtime_error("ICE RADIUS OUT OF BOUNDS");
                        ncbands = 16;
                        factor = (radice - 2.)/3.;
                        index = int(factor);
                        if (index == 43) index = 42;
                        fint = factor - float(index);
                        ib = lw_rad.ngb[ig - 1];
                        abscoice(ig) = lw_rad.absice2[ib - 1][index - 1] + fint * (lw_rad.absice2[ib - 1][index] - (lw_rad.absice2[ib - 1][index - 1])); 
                        abscosno(ig) = 0.0;
                // For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

                    }else if (iceflag >= 3){
                         //if (radice .lt. 5.0_rb .or. radice .gt. 140.0_rb) then
                         //write(errmess,'(A,i5,i5,f8.2,f8.2)' )  'ERROR: ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS',ig, lay, ciwpmc(ig,lay), radice
                         //FATAL_ERROR(errmess)
                        ncbands = 16;
                        factor = (radice - 2.)/3.;
                        index = int(factor);
                        if (index == 46) index = 45;
                        fint = factor - float(index);
                        ib = lw_rad.ngb[ig - 1];
                        abscoice(ig) = lw_rad.absice3[ib - 1][index - 1] + fint * (lw_rad.absice3[ib - 1][index] - (lw_rad.absice3[ib - 1][index - 1]));
                        abscosno(ig) = 0.0;
                    }
                    // ..Incorporate additional effects due to snow.
                    if (cswpmc(ig,lay,m,n) > 0.0 && iceflag == 5) {
                        radsno = resnmc(lay);
                        //       if (radsno < 5.0 || radsno > 140.0) then
                        //              write(errmess,'(A,i5,i5,f8.2,f8.2)' )         &
                        //    'ERROR: SNOW GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'   &
                        //    ,ig, lay, cswpmc(ig,lay), radsno
                        //              FATAL_ERROR(errmess)
                        ncbands = 16;
                        factor = (radsno - 2.)/3.;
                        index = int(factor);
                        if (index == 46) index = 45;
                        fint = factor - float(index);
                        ib = lw_rad.ngb[ig - 1];
                        abscosno(ig) = lw_rad.absice3[ib - 1][index - 1] + fint *(lw_rad.absice3[ib - 1][index] - (lw_rad.absice3[ib - 1][index - 1]));
                    }
                    // Calculation of absorption coefficients due to water clouds.
                    if (clwpmc(ig,lay,m,n) == 0.0){
                        abscoliq(ig) = 0.0;
                    }else if(liqflag == 0){
                        abscoliq(ig) = lw_rad.absliq0;
                    }else if(liqflag == 1){
                        radliq = relqmc(lay);
                        if (radliq < 2.5 || radliq > 60.) throw std::runtime_error("LIQUID EFFECTIVE RADIUS OUT OF BOUNDS");
                        index = int(radliq - 1.5);
                        if (index == 0) index = 1;
                        if (index == 58) index = 57;
                        fint = radliq - 1.5 - float(index);
                        ib = lw_rad.ngb[ig - 1];
                        abscoliq(ig) = lw_rad.absliq1[ib - 1][index - 1] + fint * (lw_rad.absliq1[ib - 1][index] - (lw_rad.absliq1[ib - 1][index - 1]));
                    }
                    taucmc(ig,lay,m,n) = ciwpmc(ig,lay,m,n) * abscoice(ig) + 
                                    clwpmc(ig,lay,m,n) * abscoliq(ig) + 
                                    cswpmc(ig,lay,m,n) * abscosno(ig);
                }
            }
        }
    }
}
KOKKOS_FUNCTION void setcoef(int &nlayers, int &istart,real1d_u &pavel,real1d_u &tavel,real1d_u &tz,real &tbound,real1d_u &semiss,real1d_u &coldry,const real4d &wkl,real1d_u &wbroad,int &laytrop,
    int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,
    real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,
    real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,
    real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,int m ,int n){
    int indbound, indlev0,indlay, indlev;
    int jp1;
    real stpfac, tbndfrac, t0frac, tlayfrac, tlevfrac;
    real dbdtlev, dbdtlay;
    real plog, fp, ft, ft1, water, scalefac, factor, compfp;
    // jm not thread safe      hvrset = '$Revision: 1.3 $'
    stpfac = 296./1013.;
    indbound = tbound - 159.;
    if (indbound < 1){
        indbound = 1;
    }else if(indbound > 180){
        indbound = 180;
    }
    tbndfrac = tbound - 159. - float(indbound);
    indlev0 = tz(0) - 159.;
    if (indlev0 < 1){
        indlev0 = 1;
    }else if(indlev0 > 180){
        indlev0 = 180;
    }
    t0frac = tz(0) - 159. - float(indlev0);
    laytrop = 0;
    /*
        ! Begin layer loop 
        !  Calculate the integrated Planck functions for each band at the
        !  surface, level, and layer temperatures.
    */
    for(int lay = 1; lay <= nlayers; lay++){
        indlay = tavel(lay) - 159.;
        if (indlay < 1){
            indlay = 1;
        }else if(indlay > 180){
            indlay = 180;
        }
        tlayfrac = tavel(lay) - 159. - float(indlay);
        indlev = tz(lay) - 159.;
        if (indlev < 1){
            indlev = 1;
        }else if(indlev > 180){
            indlev = 180;
        }
        tlevfrac = tz(lay) - 159. - float(indlev);
        // Begin spectral band loop 
        for(int iband = 1; iband <= 15; iband++){
            if (lay == 1){
                dbdtlev = lw_rad.totplnk[iband - 1][indbound] - lw_rad.totplnk[iband - 1][indbound - 1];
                plankbnd(iband) = semiss(iband) * (lw_rad.totplnk[iband - 1][indbound - 1] + tbndfrac * dbdtlev);
                dbdtlev = lw_rad.totplnk[iband - 1][indlev0]-lw_rad.totplnk[iband - 1][indlev0 - 1];
                planklev(0,iband,m,n) = lw_rad.totplnk[iband - 1][indlev0 - 1] + t0frac * dbdtlev;
            }
            dbdtlev = lw_rad.totplnk[iband - 1][indlev] - lw_rad.totplnk[iband - 1][indlev - 1];
            dbdtlay = lw_rad.totplnk[iband - 1][indlay] - lw_rad.totplnk[iband - 1][indlay - 1];
            planklay(lay,iband,m,n) = lw_rad.totplnk[iband - 1][indlay - 1] + tlayfrac * dbdtlay;
            planklev(lay,iband,m,n) = lw_rad.totplnk[iband - 1][indlev - 1] + tlevfrac * dbdtlev;
        }
        /*
            !  For band 16, if radiative transfer will be performed on just
            !  this band, use integrated Planck values up to 3250 cm-1.  
            !  If radiative transfer will be performed across all 16 bands,
            !  then include in the integrated Planck values for this band
            !  contributions from 2600 cm-1 to infinity.
        */
        int iband = 16;
        if (istart == 16){
            if (lay == 1){
                dbdtlev = lw_rad.totplk16[indbound] - lw_rad.totplk16[indbound - 1];
                plankbnd(iband) = semiss(iband) * (lw_rad.totplk16[indbound - 1] + tbndfrac * dbdtlev);
                dbdtlev = lw_rad.totplnk[iband - 1][indlev0]-lw_rad.totplnk[iband - 1][indlev0 - 1];
                planklev(0,iband,m,n) = lw_rad.totplk16[indlev0 - 1] + t0frac * dbdtlev;
            }
            dbdtlev = lw_rad.totplk16[indlev] - lw_rad.totplk16[indlev - 1];
            dbdtlay = lw_rad.totplk16[indlay] - lw_rad.totplk16[indlay - 1];
            planklay(lay,iband,m,n) = lw_rad.totplk16[indlay - 1] + tlayfrac * dbdtlay;
            planklev(lay,iband,m,n) = lw_rad.totplk16[indlev - 1] + tlevfrac * dbdtlev;
        }else{
            if (lay == 1){
                dbdtlev = lw_rad.totplnk[iband - 1][indbound] - lw_rad.totplnk[iband - 1][indbound - 1];
                plankbnd(iband) = semiss(iband) * (lw_rad.totplnk[iband - 1][indbound - 1] + tbndfrac * dbdtlev);
                dbdtlev = lw_rad.totplnk[iband - 1][indlev0]-lw_rad.totplnk[iband - 1][indlev0 - 1];
                planklev(0,iband,m,n) = lw_rad.totplnk[iband - 1][indlev0 - 1] + t0frac * dbdtlev;
            }
            dbdtlev = lw_rad.totplnk[iband - 1][indlev] - lw_rad.totplnk[iband - 1][indlev - 1];
            dbdtlay = lw_rad.totplnk[iband - 1][indlay] - lw_rad.totplnk[iband - 1][indlay - 1];
            planklay(lay,iband,m,n) = lw_rad.totplnk[iband - 1][indlay - 1] + tlayfrac * dbdtlay;
            planklev(lay,iband,m,n) = lw_rad.totplnk[iband - 1][indlev - 1] + tlevfrac * dbdtlev;
        }
        /*
            !  Find the two reference pressures on either side of the
            !  layer pressure.  Store them in JP and JP1.  Store in FP the
            !  fraction of the difference (in ln(pressure)) between these
            !  two values that the layer pressure lies.
        */
        plog = std::log(pavel(lay));
        //plog = dlog(pavel(lay))
        jp(lay) = int(36.- 5*(plog+0.04));
        if (jp(lay) < 1){
            jp(lay) = 1;
        }else if(jp(lay) > 58){
            jp(lay) = 58;
        }
        jp1 = jp(lay) + 1;
        fp = 5. *(lw_rad.preflog[jp(lay) - 1] - plog);
        /*
            !  Determine, for each reference pressure (JP and JP1), which
            !  reference temperature (these are different for each  
            !  reference pressure) is nearest the layer temperature but does
            !  not exceed it.  Store these indices in JT and JT1, resp.
            !  Store in FT (resp. FT1) the fraction of the way between JT
            !  (JT1) and the next highest reference temperature that the 
            !  layer temperature falls.
        */
        jt(lay) = int(3. + (tavel(lay)-lw_rad.tref[jp(lay) - 1])/15.);
        if (jt(lay) < 1){
            jt(lay) = 1;
        }else if(jt(lay) > 4){
            jt(lay) = 4;
        }
        ft = ((tavel(lay)-lw_rad.tref[jp(lay) - 1])/15.) - float(jt(lay)-3);
        jt1(lay) = int(3. + (tavel(lay)-lw_rad.tref[jp1 - 1])/15.);
        if (jt1(lay) < 1){
            jt1(lay) = 1;
        }else if(jt1(lay) > 4){
            jt1(lay) = 4;
        }
        ft1 = ((tavel(lay)-lw_rad.tref[jp1 - 1])/15.) - float(jt1(lay)-3);
        water = wkl(1,lay,m,n)/coldry(lay);
        scalefac = pavel(lay) * stpfac / tavel(lay);
        // If the pressure is less than ~100mb, perform a different
        // set of species interpolations.
        if (plog > 4.56){
            laytrop =  laytrop + 1;

            forfac(lay) = scalefac / (1.+water);
            factor = (332.0-tavel(lay))/36.0;
            indfor(lay) = std::min(2, std::max(1, int(factor)));
            forfrac(lay) = factor - float(indfor(lay));
            // Set up factors needed to separately include the water vapor
            // self-continuum in the calculation of absorption coefficient.
            selffac(lay) = water * forfac(lay);
            factor = (tavel(lay)-188.0)/7.2;
            indself(lay) = std::min(9, std::max(1, int(factor)-7));
            selffrac(lay) = factor - float(indself(lay) + 7);
            // Set up factors needed to separately include the minor gases
            // in the calculation of absorption coefficient
            scaleminor(lay) = pavel(lay)/tavel(lay);
            scaleminorn2(lay) = (pavel(lay)/tavel(lay)) *(wbroad(lay)/(coldry(lay)+wkl(1,lay,m,n)));
            factor = (tavel(lay)-180.8)/7.2;
            indminor(lay) = std::min(18, std::max(1, int(factor)));
            minorfrac(lay) = factor - float(indminor(lay));
            // Setup reference ratio to be used in calculation of binary
            // species parameter in lower atmosphere.
            rat_h2oco2(lay)=lw_rad.chi_mls[0][jp(lay) - 1]/lw_rad.chi_mls[1][jp(lay) - 1];
            rat_h2oco2_1(lay)=lw_rad.chi_mls[0][jp(lay)]/lw_rad.chi_mls[1][jp(lay)];

            rat_h2oo3(lay)=lw_rad.chi_mls[0][jp(lay) - 1]/lw_rad.chi_mls[2][jp(lay) - 1];
            rat_h2oo3_1(lay)=lw_rad.chi_mls[0][jp(lay)]/lw_rad.chi_mls[2][jp(lay)];

            rat_h2on2o(lay)=lw_rad.chi_mls[0][jp(lay) - 1]/lw_rad.chi_mls[3][jp(lay) - 1];
            rat_h2on2o_1(lay)=lw_rad.chi_mls[0][jp(lay)]/lw_rad.chi_mls[3][jp(lay)];

            rat_h2och4(lay)=lw_rad.chi_mls[0][jp(lay) - 1]/lw_rad.chi_mls[5][jp(lay) - 1];
            rat_h2och4_1(lay)=lw_rad.chi_mls[0][jp(lay)]/lw_rad.chi_mls[5][jp(lay)];

            rat_n2oco2(lay)=lw_rad.chi_mls[3][jp(lay) - 1]/lw_rad.chi_mls[1][jp(lay) - 1];
            rat_n2oco2_1(lay)=lw_rad.chi_mls[3][jp(lay)]/lw_rad.chi_mls[1][jp(lay)];

            // Calculate needed column amounts.
            colh2o(lay) = 1.e-20 * wkl(1,lay,m,n);
            colco2(lay) = 1.e-20 * wkl(2,lay,m,n);
            colo3(lay) = 1.e-20 * wkl(3,lay,m,n);
            coln2o(lay) = 1.e-20 * wkl(4,lay,m,n);
            colco(lay) = 1.e-20 * wkl(5,lay,m,n);
            colch4(lay) = 1.e-20 * wkl(6,lay,m,n);
            colo2(lay) = 1.e-20 * wkl(7,lay,m,n);
            if (colco2(lay) == 0.) colco2(lay) = 1.e-32 * coldry(lay);
            if (colo3(lay) == 0.) colo3(lay) = 1.e-32 * coldry(lay);
            if (coln2o(lay) == 0.) coln2o(lay) = 1.e-32 * coldry(lay);
            if (colco(lay) == 0.) colco(lay) = 1.e-32 * coldry(lay);
            if (colch4(lay) == 0.) colch4(lay) = 1.e-32 * coldry(lay);
            colbrd(lay) = 1.e-20 * wbroad(lay);
        }else{
            forfac(lay) = scalefac / (1.+water);
            factor = (tavel(lay)-188.0)/36.0;
            indfor(lay) = 3;
            forfrac(lay) = factor - 1.0;
            // Set up factors needed to separately include the water vapor
            // self-continuum in the calculation of absorption coefficient.
            selffac(lay) = water * forfac(lay);
            // Set up factors needed to separately include the minor gases
            // in the calculation of absorption coefficient
            scaleminor(lay) = pavel(lay)/tavel(lay);         
            scaleminorn2(lay) = (pavel(lay)/tavel(lay)) * (wbroad(lay)/(coldry(lay)+wkl(1,lay,m,n)));
            factor = (tavel(lay)-180.8)/7.2;
            indminor(lay) = std::min(18, std::max(1, int(factor)));
            minorfrac(lay) = factor - float(indminor(lay));
            // Setup reference ratio to be used in calculation of binary
            // species parameter in upper atmosphere.
            rat_h2oco2(lay)=lw_rad.chi_mls[0][jp(lay) - 1]/lw_rad.chi_mls[1][jp(lay) - 1];
            rat_h2oco2_1(lay)=lw_rad.chi_mls[0][jp(lay)]/lw_rad.chi_mls[1][jp(lay)];         

            rat_o3co2(lay)=lw_rad.chi_mls[2][jp(lay) - 1]/lw_rad.chi_mls[1][jp(lay) - 1];
            rat_o3co2_1(lay)=lw_rad.chi_mls[2][jp(lay)]/lw_rad.chi_mls[1][jp(lay)];
            // Calculate needed column amounts.
            colh2o(lay) = 1.e-20 * wkl(1,lay,m,n);
            colco2(lay) = 1.e-20 * wkl(2,lay,m,n);
            colo3(lay) = 1.e-20 * wkl(3,lay,m,n);
            coln2o(lay) = 1.e-20 * wkl(4,lay,m,n);
            colco(lay) = 1.e-20 * wkl(5,lay,m,n);
            colch4(lay) = 1.e-20 * wkl(6,lay,m,n);
            colo2(lay) = 1.e-20 * wkl(7,lay,m,n);
            if (colco2(lay) == 0.) colco2(lay) = 1.e-32 * coldry(lay);
            if (colo3(lay) == 0.) colo3(lay) = 1.e-32 * coldry(lay);
            if (coln2o(lay) == 0.) coln2o(lay) = 1.e-32 * coldry(lay);
            if (colco(lay)  == 0.) colco(lay) = 1.e-32 * coldry(lay);
            if (colch4(lay) == 0.) colch4(lay) = 1.e-32 * coldry(lay);
            colbrd(lay) = 1.e-20 * wbroad(lay);   
        }
        /*
            !  We have now isolated the layer ln pressure and temperature,
            !  between two reference pressures and two reference temperatures 
            !  (for each reference pressure).  We multiply the pressure 
            !  fraction FP with the appropriate temperature fractions to get 
            !  the factors that will be needed for the interpolation that yields
            !  the optical depths (performed in routines TAUGBn for band n).`     
        */
        compfp = 1. - fp;
        fac10(lay) = compfp * ft;
        fac00(lay) = compfp * (1. - ft);
        fac11(lay) = fp * ft1;
        fac01(lay) = fp * (1. - ft1);
        // Rescale selffac and forfac for use in taumol
        selffac(lay) = colh2o(lay)*selffac(lay);
        forfac(lay) = colh2o(lay)*forfac(lay);
    
    }
}
KOKKOS_FUNCTION void rtrnmc(int &nlayers,int &istart,int &iend,int &iout,real1d_u &pz,real1d_u &semiss,int &ncbands,const real4d &cldfmc,const real4d &taucmc,const real4d &planklay,const real4d &planklev,
real1d_u &plankbnd,real &pwvcm,const real4d &fracs,const real4d &taut,real1d_u &totuflux,real1d_u &totdflux,real1d_u &fnet,real1d_u &htr,real1d_u &totuclfl,real1d_u &totdclfl,
real1d_u &fnetc,real1d_u &htrc,const real1d_u &tau_tbl,const real1d_u &tfn_tbl,const real1d_u &exp_tbl,const real4d &efclfrac,const real4d &odcld,int m,int n){
//const real3d &bbugas,const real3d &bbutot,const real3d &secdiff,const int3d &icldlyr,const real3d &uclfl,const real3d &dclfl,const real3d &clrurad,const real3d &clrdrad,const real3d &atot,const real3d &atrans,int m,int n){
    /*
        !-----------------------------------------------------------------------------
        !
        !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
        !  Revision for GCMs:  Michael J. Iacono; October, 2002
        !  Revision for F90:  Michael J. Iacono; June, 2006
        !
        !  This program calculates the upward fluxes, downward fluxes, and
        !  heating rates for an arbitrary clear or cloudy atmosphere.  The input
        !  to this program is the atmospheric profile, all Planck function
        !  information, and the cloud fraction by layer.  A variable diffusivity 
        !  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
        !  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
        !  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
        !  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
        !  use of the emissivity angle for the flux integration can cause errors of 
        !  1 to 4 W/m2 within cloudy layers.  
        !  Clouds are treated with the McICA stochastic approach and maximum-random
        !  cloud overlap. 
        !***************************************************************************
    */
    const int size = 110;
    real atot[size];
    real atrans[size];
    real bbugas[size];
    real bbutot[size];
    real uflux[size];
    real dflux[size];
    real uclfl[size];
    real dclfl[size];
    real secdiff[nbndlw + 1];
    real clrurad[size];
    real clrdrad[size];
    real urad[size];
    real drad[size];
    int icldlyr[size]; 
    real wtdiff = 0.5;
    real rec_6 = 0.166667;
    real transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn;
    real odepth, odtot, odepth_rec, odtot_rec, gassrc;
    real tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac;
    real rad0, reflect, radlu, radclru;
    int l,ib;
    int igc;                                        // g-point interval counter
    int iclddn;                                     // flag for cloud in down path
    int ittot, itgas, itr;                          // lookup table indices
    /*
        ! ------- Definitions -------
        ! input
        !    nlayers                      ! number of model layers
        !    ngptlw                       ! total number of g-point subintervals
        !    nbndlw                       ! number of longwave spectral bands
        !    ncbands                      ! number of spectral bands for clouds
        !    secdiff                      ! diffusivity angle
        !    wtdiff                       ! weight for radiance to flux conversion
        !    pavel                        ! layer pressures (mb)
        !    pz                           ! level (interface) pressures (mb)
        !    tavel                        ! layer temperatures (k)
        !    tz                           ! level (interface) temperatures(mb)
        !    tbound                       ! surface temperature (k)
        !    cldfrac                      ! layer cloud fraction
        !    taucloud                     ! layer cloud optical depth
        !    itr                          ! integer look-up table index
        !    icldlyr                      ! flag for cloudy layers
        !    iclddn                       ! flag for cloud in column at any layer
        !    semiss                       ! surface emissivities for each band
        !    reflect                      ! surface reflectance
        !    bpade                        ! 1/(pade constant)
        !    tau_tbl                      ! clear sky optical depth look-up table
        !    exp_tbl                      ! exponential look-up table for transmittance
        !    tfn_tbl                      ! tau transition function look-up table

        ! local
        !    atrans                       ! gaseous absorptivity
        !    abscld                       ! cloud absorptivity
        !    atot                         ! combined gaseous and cloud absorptivity
        !    odclr                        ! clear sky (gaseous) optical depth
        !    odcld                        ! cloud optical depth
        !    odtot                        ! optical depth of gas and cloud
        !    tfacgas                      ! gas-only pade factor, used for planck fn
        !    tfactot                      ! gas and cloud pade factor, used for planck fn
        !    bbdgas                       ! gas-only planck function for downward rt
        !    bbugas                       ! gas-only planck function for upward rt
        !    bbdtot                       ! gas and cloud planck function for downward rt
        !    bbutot                       ! gas and cloud planck function for upward calc.
        !    gassrc                       ! source radiance due to gas only
        !    efclfrac                     ! effective cloud fraction
        !    radlu                        ! spectrally summed upward radiance 
        !    radclru                      ! spectrally summed clear sky upward radiance 
        !    urad                         ! upward radiance by layer
        !    clrurad                      ! clear sky upward radiance by layer
        !    radld                        ! spectrally summed downward radiance 
        !    radclrd                      ! spectrally summed clear sky downward radiance 
        !    drad                         ! downward radiance by layer
        !    clrdrad                      ! clear sky downward radiance by layer

        ! output
        !    totuflux                     ! upward longwave flux (w/m2)
        !    totdflux                     ! downward longwave flux (w/m2)
        !    fnet                         ! net longwave flux (w/m2)
        !    htr                          ! longwave heating rate (k/day)
        !    totuclfl                     ! clear sky upward longwave flux (w/m2)
        !    totdclfl                     ! clear sky downward longwave flux (w/m2)
        !    fnetc                        ! clear sky net longwave flux (w/m2)
        !    htrc                         ! clear sky longwave heating rate (k/day)
    */
   for(int ibnd = 1; ibnd <= nbndlw; ibnd++){
        if (ibnd == 1 || ibnd == 4 || ibnd >= 10){
            secdiff[ibnd] = 1.66;
        }else{
            secdiff[ibnd] = lw_rad.a0[ibnd - 1] + lw_rad.a1[ibnd - 1]*std::exp(lw_rad.a2[ibnd - 1]*pwvcm);
            if (secdiff[ibnd] > 1.80) secdiff[ibnd] = 1.80;
            if (secdiff[ibnd] < 1.50) secdiff[ibnd] = 1.50;
        }
    }
    for(int lay = 0; lay <= nlayers; lay++){
        urad[lay] = 0.0;
        drad[lay] = 0.0;
        totuflux(lay) = 0.0;
        totdflux(lay) = 0.0;
        clrurad[lay] = 0.0;
        clrdrad[lay] = 0.0;
        totuclfl(lay) = 0.0;
        totdclfl(lay) = 0.0;
        icldlyr[lay] = 0;
    }
    for(int lay = 1; lay <= nlayers; lay++){
        for(int ig = 1; ig <= ngptlw; ig++){
            if(cldfmc(ig, lay,m,n) == 1.){
                double tauc = taucmc(ig, lay,m,n);
                int ib = lw_rad.ngb[ig - 1];
                double secd = secdiff[ib];
                double od = secd * tauc;
                double transcld = exp(-od);
                double abs_cld = cldfmc(ig, lay,m,n) * (1.0 - transcld);
                odcld(ig, lay,m,n) = od;
                efclfrac(ig, lay,m,n) = abs_cld;
                icldlyr[lay] = 1;
            }else{
                efclfrac(ig, lay,m,n) = 0.;
                odcld(ig, lay,m,n) = 0.;
            }
        }
    }

    igc = 1;
    // Loop over frequency bands.
    for(int iband = istart; iband <= iend; iband++){
        // Reinitialize g-point counter for each band if output for each band is requested.
        if (iout > 0 && iband >= 2) igc = lw_rad.ngs[iband - 2] + 1;
        do{
            // Radiative transfer starts here.
            radld = 0.;
            radclrd = 0.;
            iclddn = 0;
            // Downward radiative transfer loop.  
            for(int lev = nlayers; lev >= 1; lev--){
                plfrac = fracs(igc,lev,m,n);
                blay = planklay(lev,iband,m,n);
                dplankup = planklev(lev,iband,m,n) - blay;
                dplankdn = planklev(lev-1,iband,m,n) - blay;
                odepth = secdiff[iband] * taut(igc,lev,m,n);
                if (odepth < 0.0) odepth = 0.0;
                // Cloudy layer/
                if (icldlyr[lev] == 1){
                    iclddn = 1;
                    odtot = odepth + odcld(igc,lev,m,n);
                    if (odtot < 0.06){
                        atrans[lev] = odepth - 0.5*odepth*odepth;
                        odepth_rec = rec_6*odepth;
                        gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans[lev];

                        atot[lev] =  odtot - 0.5*odtot*odtot;
                        odtot_rec = rec_6*odtot;
                        bbdtot =  plfrac * (blay+dplankdn*odtot_rec);
                        bbd = plfrac*(blay+dplankdn*odepth_rec);
                        radld = radld - radld * (atrans[lev] + efclfrac(igc,lev,m,n) * (1. - atrans[lev])) + gassrc + cldfmc(igc,lev,m,n) * (bbdtot * atot[lev] - gassrc);
                        drad[lev-1] = drad[lev-1] + radld;
                        bbugas[lev] =  plfrac * (blay+dplankup*odepth_rec);
                        bbutot[lev] =  plfrac * (blay+dplankup*odtot_rec);
                    }else if(odepth <= 0.06){
                        atrans[lev] = odepth - 0.5*odepth*odepth;
                        odepth_rec = rec_6*odepth;
                        gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans[lev];
                        odtot = odepth + odcld(igc,lev,m,n);
                        tblind = odtot/(lw_rad.bpade+odtot);
                        ittot = lw_rad.tblint*tblind + 0.5;
                        tfactot = tfn_tbl(ittot);
                        bbdtot = plfrac * (blay + tfactot*dplankdn);
                        bbd = plfrac*(blay+dplankdn*odepth_rec);
                        atot[lev] = 1. - exp_tbl(ittot);
                        radld = radld - radld * (atrans[lev] + efclfrac(igc,lev,m,n) * (1. - atrans[lev])) + gassrc + cldfmc(igc,lev,m,n) * (bbdtot * atot[lev] - gassrc);
                        drad[lev-1] = drad[lev-1] + radld;
                        bbugas[lev] = plfrac * (blay + dplankup*odepth_rec);
                        bbutot[lev] = plfrac * (blay + tfactot * dplankup);
                    }else{
                        tblind = odepth/(lw_rad.bpade+odepth);
                        itgas = lw_rad.tblint*tblind+0.5;
                        odepth = tau_tbl(itgas);
                        atrans[lev] = 1. - exp_tbl(itgas);
                        tfacgas = tfn_tbl(itgas);
                        gassrc = atrans[lev] * plfrac * (blay + tfacgas*dplankdn);
                        odtot = odepth + odcld(igc,lev,m,n);
                        tblind = odtot/(lw_rad.bpade+odtot);
                        ittot = lw_rad.tblint*tblind + 0.5;
                        tfactot = tfn_tbl(ittot);
                        bbdtot = plfrac * (blay + tfactot*dplankdn);
                        bbd = plfrac*(blay+tfacgas*dplankdn);
                        atot[lev] = 1. - exp_tbl(ittot);
                        radld = radld - radld * (atrans[lev] + efclfrac(igc,lev,m,n) * (1. - atrans[lev])) + gassrc + cldfmc(igc,lev,m,n) * (bbdtot * atot[lev] - gassrc);
                        drad[lev-1] = drad[lev-1] + radld;
                        bbugas[lev] = plfrac * (blay + tfacgas * dplankup);
                        bbutot[lev] = plfrac * (blay + tfactot * dplankup);
                    }
						// std::cout<<bbutot<<std::endl;
                }else{
                    if (odepth <= 0.06){
                        atrans[lev] = odepth-0.5*odepth*odepth;
                        odepth = rec_6*odepth;
                        bbd = plfrac*(blay+dplankdn*odepth);
                        bbugas[lev] = plfrac*(blay+dplankup*odepth);
                    }else{
                        tblind = odepth/(lw_rad.bpade+odepth);
                        itr = int(lw_rad.tblint*tblind+0.5);
                        transc = exp_tbl(itr);
                        atrans[lev] = 1.-transc;
                        tausfac = tfn_tbl(itr);
                        bbd = plfrac*(blay+tausfac*dplankdn);
                        bbugas[lev] = plfrac * (blay + tausfac * dplankup);
                    }   
                    radld = radld + (bbd-radld)*atrans[lev];
                    drad[lev-1] = drad[lev-1] + radld;
                }
                /*
                    !  Set clear sky stream to total sky stream as long as layers
                    !  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
                    !  and clear sky stream must be computed separately from that point.
                */
                if(iclddn == 1){
                    radclrd = radclrd + (bbd-radclrd) * atrans[lev];
                    clrdrad[lev-1] = clrdrad[lev-1] + radclrd;
                }else{
                    radclrd = radld;
                    clrdrad[lev-1] = drad[lev-1];
                }
            }
            /*
                ! Spectral emissivity & reflectance
                !  Include the contribution of spectrally varying longwave emissivity
                !  and reflection from the surface to the upward radiative transfer.
                !  Note: Spectral and Lambertian reflection are identical for the
                !  diffusivity angle flux integration used here.
            */
            rad0 = fracs(igc,1,m,n) * plankbnd(iband);
            // Add in specular reflection of surface downward radiance.
            reflect = 1. - semiss(iband);
            radlu = rad0 + reflect * radld;
            radclru = rad0 + reflect * radclrd;
            // Upward radiative transfer loop.
            urad[0] = urad[0] + radlu;
            clrurad[0] = clrurad[0] + radclru;
            for(int lev = 1; lev <= nlayers; lev++){
                // Cloudy layer
                if (icldlyr[lev] == 1){
                    gassrc = bbugas[lev] * atrans[lev];
                    radlu = radlu - radlu * (atrans[lev] + efclfrac(igc,lev,m,n) * (1. - atrans[lev])) + gassrc + cldfmc(igc,lev,m,n) * (bbutot[lev] * atot[lev] - gassrc);
                    urad[lev] = urad[lev] + radlu;
                }else{
                    radlu = radlu + (bbugas[lev]-radlu)*atrans[lev];
                    urad[lev] = urad[lev] + radlu;
                }
                /*
                    !  Set clear sky stream to total sky stream as long as all layers
                    !  are clear (iclddn=0).  Streams must be calculated separately at 
                    !  all layers when a cloud is present (ICLDDN=1), because surface 
                    !  reflectance is different for each stream.
                */
                if (iclddn == 1){
                    radclru = radclru + (bbugas[lev]-radclru)*atrans[lev]; 
                    clrurad[lev] = clrurad[lev] + radclru;
                }else{
                    radclru = radlu;
                    clrurad[lev] = urad[lev];
                }
            }
            igc++;
        }while(igc <= lw_rad.ngs[iband -1] );
        // Process longwave output from band for total and clear streams.
        // Calculate upward, downward, and net flux.
        for(int lev = nlayers; lev >= 0; lev--){
            uflux[lev] = urad[lev]*wtdiff;
            dflux[lev] = drad[lev]*wtdiff;
            urad[lev] = 0.0;
            drad[lev] = 0.0;
            totuflux(lev) = totuflux(lev) + uflux[lev] * lw_rad.delwave[iband - 1];
            totdflux(lev) = totdflux(lev) + dflux[lev] * lw_rad.delwave[iband - 1];
            uclfl[lev] = clrurad[lev]*wtdiff;
            dclfl[lev] = clrdrad[lev]*wtdiff;
            clrurad[lev] = 0.0;
            clrdrad[lev] = 0.0;
            totuclfl(lev) = totuclfl(lev) + uclfl[lev] * lw_rad.delwave[iband - 1];
            totdclfl(lev) = totdclfl(lev) + dclfl[lev] * lw_rad.delwave[iband - 1];
        }
    }
    // Calculate fluxes at surface
    totuflux(0) = totuflux(0) * lw_rad.fluxfac;
    totdflux(0) = totdflux(0) * lw_rad.fluxfac;
    fnet(0) = totuflux(0) - totdflux(0);
    totuclfl(0) = totuclfl(0) * lw_rad.fluxfac;
    totdclfl(0) = totdclfl(0) * lw_rad.fluxfac;
    fnetc(0) = totuclfl(0) - totdclfl(0);
    // Calculate fluxes at model levels
    for(int lev = 1; lev <= nlayers; lev++){
        totuflux(lev) = totuflux(lev) * lw_rad.fluxfac;
        totdflux(lev) = totdflux(lev) * lw_rad.fluxfac;
        fnet(lev) = totuflux(lev) - totdflux(lev);
        totuclfl(lev) = totuclfl(lev) * lw_rad.fluxfac;
        totdclfl(lev) = totdclfl(lev) * lw_rad.fluxfac;
        fnetc(lev) = totuclfl(lev) - totdclfl(lev);
        int l = lev - 1;
        // Calculate heating rates at model layers
        htr(l)=lw_rad.heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev));
        htrc(l)=lw_rad.heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev));
    }
    // Set heating rate to zero in top layer
    htr(nlayers) = 0.0;
    htrc(nlayers) = 0.0;
}
KOKKOS_FUNCTION void taugb1(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
            real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
            real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
            real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
            real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
            const real1d_u &fracrefa_01,const real1d_u &fracrefb_01,const real2d_u &absa_01,const real3d_u &ka_01,const real2d_u &absb_01,const real3d_u &kb_01,const real2d_u &ka_mn2_01,
            const real2d_u &kb_mn2_01,const real2d_u &selfref_01,const real2d_u &forref_01,int m,int n){
/*
    ! ------- Modifications -------
    !  Written by Eli J. Mlawer, Atmospheric & Environmental Research.
    !  Revised by Michael J. Iacono, Atmospheric & Environmental Research.
    !
    !     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
    !                          (high key - h2o; high minor - n2)
    !
    !     note: previous versions of rrtm band 1: 
    !           10-250 cm-1 (low - h2o; high - h2o)
    !----------------------------------------------------------------------------
*/
    // int lay,ind0,ind1,inds,indf,indm,ig;
    // real pp,corradj,scalen2,tauself,taufor,taun2;
/*
    ! Minor gas mapping levels:
    !     lower - n2, p = 142.5490 mbar, t = 215.70 k
    !     upper - n2, p = 142.5490 mbar, t = 215.70 k

    ! Compute the optical depth by interpolating in ln(pressure) and 
    ! temperature.  Below laytrop, the water vapor self-continuum and
    ! foreign continuum is interpolated (in temperature) separately.

    ! Lower atmosphere loop
*/

	int lay, inds, indf, indm, ind0, ind1;
    double pp, corradj,taufor,taun2,scalen2;

	for(int lay = 1; lay <= laytrop; lay++){
	    int ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[0] + 1; 
        int ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[0] + 1;
        int inds = indself(lay);
        int indf = indfor(lay);
        int indm = indminor(lay);
        real corradj = 1.;
        real pp = pavel(lay);
        if(pp < 250.)
        corradj = 1. - 0.15 * (250. - pp) / 154.4;
        for(int ig = 1; ig <= lw_rad.ng1; ig++){
            real scalen2 = colbrd(lay) * scaleminorn2(lay);
            real tauself = selffac(lay) * (selfref_01(ig,inds) + selffrac(lay) * (selfref_01(ig,inds+1) - selfref_01(ig,inds)));
            real taun2_lower = scalen2*(ka_mn2_01(ig,indm) + minorfrac(lay) * (ka_mn2_01(ig,indm+1) - ka_mn2_01(ig,indm)));
            real taufor =  forfac(lay) * (forref_01(ig,indf) + forfrac(lay) * (forref_01(ig,indf+1) -  forref_01(ig,indf)));
            taug(ig,lay,m,n) = corradj * (colh2o(lay) * (fac00(lay) * absa_01(ig,ind0) + fac10(lay) * absa_01(ig,ind0+1) + fac01(lay) * absa_01(ig,ind1) + fac11(lay) * absa_01(ig,ind1+1)) + tauself + taufor + taun2_lower);
            fracs(ig,lay,m,n) = fracrefa_01(ig);
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
		ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[0] + 1;
		ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[0] + 1;
		indf = indfor(lay);
		indm = indminor(lay);
		corradj = 1. - 0.15 * (pavel(lay) / 95.6);
		scalen2 = colbrd(lay) * scaleminorn2(lay);
        for(int ig = 1; ig <= lw_rad.ng1; ig++){
            taufor = forfac(lay) * (forref_01(ig,indf) + forfrac(lay) * (forref_01(ig,indf+1) - forref_01(ig,indf))); 
            taun2 = scalen2*(kb_mn2_01(ig,indm) + minorfrac(lay) * (kb_mn2_01(ig,indm+1) - kb_mn2_01(ig,indm)));
            taug(ig,lay,m,n) = corradj * (colh2o(lay) * (fac00(lay) * absb_01(ig,ind0) + fac10(lay) * absb_01(ig,ind0+1) + fac01(lay) * absb_01(ig,ind1) + fac11(lay) * absb_01(ig,ind1+1)) + taufor + taun2);
            fracs(ig,lay,m,n) = fracrefb_01(ig);
        }
    }
}
KOKKOS_FUNCTION void taugb2(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real1d_u &fracrefa_02,const real1d_u &fracrefb_02,const real2d_u &absa_02,const real3d_u &ka_02,const real2d_u &absb_02,const real3d_u &kb_02,const real2d_u &selfref_02,const real2d_u &forref_02,int m,int n){
    /*
        band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
        note: previous version of rrtm band 2: 
                250 - 500 cm-1 (low - h2o; high - h2o)
    */
    int lay, ind0, ind1, inds, indf, ig;
    real pp, corradj, tauself, taufor;
    // Compute the optical depth by interpolating in ln(pressure) and 
    // temperature.  Below laytrop, the water vapor self-continuum and
    // foreign continuum is interpolated (in temperature) separately.
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[1] + 1;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[1] + 1;
        inds = indself(lay);
        indf = indfor(lay);
        pp = pavel(lay);
        corradj = 1. - .05 * (pp - 100.) / 900.;
        for(int ig = 1; ig <= lw_rad.ng2; ig++){
            tauself = selffac(lay) * (selfref_02(ig,inds) + selffrac(lay) * (selfref_02(ig,inds+1) - selfref_02(ig,inds)));
            taufor =  forfac(lay) * (forref_02(ig,indf) + forfrac(lay) * (forref_02(ig,indf+1) - forref_02(ig,indf)));
            taug(lw_rad.ngs1+ig,lay,m,n) = corradj * (colh2o(lay) * (fac00(lay) * absa_02(ig,ind0) + fac10(lay) * absa_02(ig,ind0+1) + fac01(lay) * absa_02(ig,ind1) + fac11(lay) * absa_02(ig,ind1+1)) + tauself + taufor);
            fracs(lw_rad.ngs1+ig,lay,m,n) = fracrefa_02(ig);
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[1] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[1] + 1;
        indf = indfor(lay);
        for(int ig = 1; ig <= lw_rad.ng2; ig++){
            taufor =  forfac(lay) * (forref_02(ig,indf) + forfrac(lay) * (forref_02(ig,indf+1) - forref_02(ig,indf))); 
            taug(lw_rad.ngs1+ig,lay,m,n) = colh2o(lay) * (fac00(lay) * absb_02(ig,ind0) + fac10(lay) * absb_02(ig,ind0+1) + fac01(lay) * absb_02(ig,ind1) + fac11(lay) * absb_02(ig,ind1+1)) + taufor;
            fracs(lw_rad.ngs1+ig,lay,m,n) = fracrefb_02(ig);
        }
    }

}
KOKKOS_FUNCTION
void taugb3(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_03,const real2d_u &fracrefb_03,const real2d_u &absa_03,const real4d_u &ka_03,const real2d_u &absb_03,const real4d_u &kb_03,const real3d_u &ka_mn2o_03,const real3d_u &kb_mn2o_03,const real2d_u &selfref_03,const real2d_u &forref_03,int m,int n){
    /*
        band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
                            (high key - h2o,co2; high minor - n2o)
    */
    int lay, ind0, ind1, inds, indf, indm, ig, js, js1, jmn2o, jpl;
    real speccomb, specparm, specmult, fs, speccomb1, specparm1, specmult1, fs1, speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o, fmn2omf, chi_n2o, ratn2o, adjfac, adjcoln2o;
    real speccomb_planck, specparm_planck, specmult_planck, fpl, p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor, n2om1, n2om2, absn2o, refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b, tau_major, tau_major1;
    // Minor gas mapping levels:
    // lower - n2o, p = 706.272 mbar, t = 278.94 k
    // upper - n2o, p = 95.58 mbar, t = 215.7 k

    // P = 212.725 mb
    refrat_planck_a = lw_rad.chi_mls[0][8]/lw_rad.chi_mls[1][8];

    // P = 95.58 mb
    refrat_planck_b = lw_rad.chi_mls[0][12]/lw_rad.chi_mls[1][12];

    // P = 706.270mb
    refrat_m_a = lw_rad.chi_mls[0][2]/lw_rad.chi_mls[1][2];

    // P = 95.58 mb 
    refrat_m_b = lw_rad.chi_mls[0][12]/lw_rad.chi_mls[1][12];

    /*
        ! Compute the optical depth by interpolating in ln(pressure) and 
        ! temperature, and appropriate species.  Below laytrop, the water vapor 
        ! self-continuum and foreign continuum is interpolated (in temperature) 
        ! separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult); 
        fs = fmod(specmult,1.0);     
        speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);
        speccomb_mn2o = colh2o(lay) + refrat_m_a*colco2(lay);
        specparm_mn2o = colh2o(lay)/speccomb_mn2o;
        if (specparm_mn2o >= lw_rad.oneminus) specparm_mn2o = lw_rad.oneminus;
        specmult_mn2o = 8.*specparm_mn2o;
        jmn2o = 1 + int(specmult_mn2o);
        fmn2o = fmod(specmult_mn2o,1.0);
        fmn2omf = minorfrac(lay)*fmn2o;
        /*
            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_n2o = coln2o(lay)/coldry(lay);
        ratn2o = 1.e20*chi_n2o/lw_rad.chi_mls[3][jp(lay)];
        if (ratn2o > 1.5){
            adjfac = 0.5+std::pow((ratn2o-0.5),0.65);
            adjcoln2o = adjfac*lw_rad.chi_mls[3][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcoln2o = coln2o(lay);
        }

        speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);
        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[2] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[2] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);

        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }
        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng3; ig++){
            tauself = selffac(lay)* (selfref_03(ig,inds) + selffrac(lay) * (selfref_03(ig,inds+1) - selfref_03(ig,inds)));
            taufor = forfac(lay) * (forref_03(ig,indf) + forfrac(lay) * (forref_03(ig,indf+1) - forref_03(ig,indf)));
            n2om1 = ka_mn2o_03(ig,jmn2o,indm) + fmn2o * (ka_mn2o_03(ig,jmn2o+1,indm) - ka_mn2o_03(ig,jmn2o,indm));
            n2om2 = ka_mn2o_03(ig,jmn2o,indm+1) + fmn2o * (ka_mn2o_03(ig,jmn2o+1,indm+1) - ka_mn2o_03(ig,jmn2o,indm+1));
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1);
            if (specparm < 0.125){
               tau_major = speccomb * (fac000 * absa_03(ig,ind0) + fac100 * absa_03(ig,ind0+1) + fac200 * absa_03(ig,ind0+2) + fac010 * absa_03(ig,ind0+9) + fac110 * absa_03(ig,ind0+10) + fac210 * absa_03(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * (fac200 * absa_03(ig,ind0-1) + fac100 * absa_03(ig,ind0) + fac000 * absa_03(ig,ind0+1) + fac210 * absa_03(ig,ind0+8) + fac110 * absa_03(ig,ind0+9) + fac010 * absa_03(ig,ind0+10));
            }else{
               tau_major = speccomb * (fac000 * absa_03(ig,ind0) + fac100 * absa_03(ig,ind0+1) + fac010 * absa_03(ig,ind0+9) + fac110 * absa_03(ig,ind0+10));
            }
            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * (fac001 * absa_03(ig,ind1) + fac101 * absa_03(ig,ind1+1) + fac201 * absa_03(ig,ind1+2) + fac011 * absa_03(ig,ind1+9) + fac111 * absa_03(ig,ind1+10) + fac211 * absa_03(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * (fac201 * absa_03(ig,ind1-1) + fac101 * absa_03(ig,ind1) + fac001 * absa_03(ig,ind1+1) + fac211 * absa_03(ig,ind1+8) + fac111 * absa_03(ig,ind1+9) + fac011 * absa_03(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * (fac001 * absa_03(ig,ind1) +  fac101 * absa_03(ig,ind1+1) + fac011 * absa_03(ig,ind1+9) + fac111 * absa_03(ig,ind1+10));
            }
            taug(lw_rad.ngs2+ig,lay,m,n) = tau_major + tau_major1 + tauself + taufor + adjcoln2o*absn2o;
            fracs(lw_rad.ngs2+ig,lay,m,n) = fracrefa_03(ig,jpl) + fpl * (fracrefa_03(ig,jpl+1)-fracrefa_03(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 4.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);
        speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 4.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        fac000 = (1. - fs) * fac00(lay);
        fac010 = (1. - fs) * fac10(lay);
        fac100 = fs * fac00(lay);
        fac110 = fs * fac10(lay);
        fac001 = (1. - fs1) * fac01(lay);
        fac011 = (1. - fs1) * fac11(lay);
        fac101 = fs1 * fac01(lay);
        fac111 = fs1 * fac11(lay);

        speccomb_mn2o = colh2o(lay) + refrat_m_b*colco2(lay);
        specparm_mn2o = colh2o(lay)/speccomb_mn2o;

        if (specparm_mn2o >= lw_rad.oneminus) specparm_mn2o = lw_rad.oneminus;
        specmult_mn2o = 4.*specparm_mn2o;
        jmn2o = 1 + int(specmult_mn2o);
        fmn2o = fmod(specmult_mn2o,1.0);
        fmn2omf = minorfrac(lay)*fmn2o;
        /*
            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_n2o = coln2o(lay)/coldry(lay);
        ratn2o = 1.e20*chi_n2o/lw_rad.chi_mls[3][jp(lay)];
        if (ratn2o > 1.5){
            adjfac = 0.5+std::pow((ratn2o-0.5),0.65);
            adjcoln2o = adjfac*lw_rad.chi_mls[3][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcoln2o = coln2o(lay);
        }

        speccomb_planck = colh2o(lay)+refrat_planck_b*colco2(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 4.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[2] + js;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[2] + js1;
        indf = indfor(lay);
        indm = indminor(lay);
        for(int ig = 1; ig <= lw_rad.ng3; ig++){
            taufor = forfac(lay) * (forref_03(ig,indf) + forfrac(lay) * (forref_03(ig,indf+1) - forref_03(ig,indf)));
            n2om1 = kb_mn2o_03(ig,jmn2o,indm) + fmn2o * (kb_mn2o_03(ig,jmn2o+1,indm)-kb_mn2o_03(ig,jmn2o,indm));
            n2om2 = kb_mn2o_03(ig,jmn2o,indm+1) + fmn2o * (kb_mn2o_03(ig,jmn2o+1,indm+1)-kb_mn2o_03(ig,jmn2o,indm+1));
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1);
            taug(lw_rad.ngs2+ig,lay,m,n) = speccomb * (fac000 * absb_03(ig,ind0) + fac100 * absb_03(ig,ind0+1) + fac010 * absb_03(ig,ind0+5) + fac110 * absb_03(ig,ind0+6)) + speccomb1 * (fac001 * absb_03(ig,ind1) 
            +  fac101 * absb_03(ig,ind1+1) + fac011 * absb_03(ig,ind1+5) + fac111 * absb_03(ig,ind1+6))  + taufor + adjcoln2o*absn2o;
            fracs(lw_rad.ngs2+ig,lay,m,n) = fracrefb_03(ig,jpl) + fpl * (fracrefb_03(ig,jpl+1)-fracrefb_03(ig,jpl));
        }
    }

}
KOKKOS_FUNCTION
void taugb4(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_04,const real2d_u &fracrefb_04,const real2d_u &absa_04,const real4d_u &ka_04,const real2d_u &absb_04,const real4d_u &kb_04,const real2d_u &selfref_04,const real2d_u &forref_04,int m,int n){
    /* band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2) */
    int lay, ind0, ind1, inds, indf, ig;
    int js, js1, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor;
    real refrat_planck_a, refrat_planck_b;
    real tau_major, tau_major1;
    // P =   142.5940 mb
    refrat_planck_a = lw_rad.chi_mls[0][10]/lw_rad.chi_mls[1][10];
    // P = 95.58350 mb
    refrat_planck_b = lw_rad.chi_mls[2][12]/lw_rad.chi_mls[1][12];
    /*
        ! Compute the optical depth by interpolating in ln(pressure) and 
        ! temperature, and appropriate species.  Below laytrop, the water 
        ! vapor self-continuum and foreign continuum is interpolated (in temperature) 
        ! separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[3] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[3] + js1;
        inds = indself(lay);
        indf = indfor(lay);

        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }

        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng4; ig++){
            tauself = selffac(lay)* (selfref_04(ig,inds) + selffrac(lay) * (selfref_04(ig,inds+1) - selfref_04(ig,inds)));
            taufor =  forfac(lay) * (forref_04(ig,indf) + forfrac(lay) * (forref_04(ig,indf+1) - forref_04(ig,indf)));
            if (specparm < 0.125){
               tau_major = speccomb * 
                    (fac000 * absa_04(ig,ind0) + 
                    fac100 * absa_04(ig,ind0+1) + 
                    fac200 * absa_04(ig,ind0+2) + 
                    fac010 * absa_04(ig,ind0+9) + 
                    fac110 * absa_04(ig,ind0+10) + 
                    fac210 * absa_04(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * 
                    (fac200 * absa_04(ig,ind0-1) + 
                    fac100 * absa_04(ig,ind0) + 
                    fac000 * absa_04(ig,ind0+1) + 
                    fac210 * absa_04(ig,ind0+8) + 
                    fac110 * absa_04(ig,ind0+9) + 
                    fac010 * absa_04(ig,ind0+10));
            }else{
               tau_major = speccomb * 
                    (fac000 * absa_04(ig,ind0) + 
                    fac100 * absa_04(ig,ind0+1) + 
                    fac010 * absa_04(ig,ind0+9) + 
                    fac110 * absa_04(ig,ind0+10));
            }

            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * 
                    (fac001 * absa_04(ig,ind1) +  
                    fac101 * absa_04(ig,ind1+1) + 
                    fac201 * absa_04(ig,ind1+2) + 
                    fac011 * absa_04(ig,ind1+9) + 
                    fac111 * absa_04(ig,ind1+10) + 
                    fac211 * absa_04(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * 
                    (fac201 * absa_04(ig,ind1-1) + 
                    fac101 * absa_04(ig,ind1) + 
                    fac001 * absa_04(ig,ind1+1) + 
                    fac211 * absa_04(ig,ind1+8) + 
                    fac111 * absa_04(ig,ind1+9) + 
                    fac011 * absa_04(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * 
                    (fac001 * absa_04(ig,ind1) + 
                    fac101 * absa_04(ig,ind1+1) + 
                    fac011 * absa_04(ig,ind1+9) + 
                    fac111 * absa_04(ig,ind1+10));
            }
            taug(lw_rad.ngs3+ig,lay,m,n) = tau_major + tau_major1 + tauself + taufor;
            fracs(lw_rad.ngs3+ig,lay,m,n) = fracrefa_04(ig,jpl) + fpl * (fracrefa_04(ig,jpl+1)-fracrefa_04(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay);
        specparm = colo3(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 4.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay);
        specparm1 = colo3(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 4.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        fac000 = (1. - fs) * fac00(lay);
        fac010 = (1. - fs) * fac10(lay);
        fac100 = fs * fac00(lay);
        fac110 = fs * fac10(lay);
        fac001 = (1. - fs1) * fac01(lay);
        fac011 = (1. - fs1) * fac11(lay);
        fac101 = fs1 * fac01(lay);
        fac111 = fs1 * fac11(lay);

        speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay);
        specparm_planck = colo3(lay)/speccomb_planck;
        if (specparm_planck > lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 4.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[3] + js;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[3] + js1;
        for(int ig = 1; ig <= lw_rad.ng4; ig++){
            taug(lw_rad.ngs3+ig,lay,m,n) =  speccomb * 
                (fac000 * absb_04(ig,ind0) + 
                fac100 * absb_04(ig,ind0+1) + 
                fac010 * absb_04(ig,ind0+5) + 
                fac110 * absb_04(ig,ind0+6)) 
                + speccomb1 * 
                (fac001 * absb_04(ig,ind1) +  
                fac101 * absb_04(ig,ind1+1) + 
                fac011 * absb_04(ig,ind1+5) + 
                fac111 * absb_04(ig,ind1+6));
            fracs(lw_rad.ngs3+ig,lay,m,n) = fracrefb_04(ig,jpl) + fpl * 
                (fracrefb_04(ig,jpl+1)-fracrefb_04(ig,jpl));
        }
        /*
            ! Empirical modification to code to improve stratospheric cooling rates
            ! for co2.  Revised to apply weighting for g-point reduction in this band.
        */
        taug(lw_rad.ngs3+8,lay,m,n)*=0.92;
        taug(lw_rad.ngs3+9,lay,m,n)*=0.88;
        taug(lw_rad.ngs3+10,lay,m,n)*=1.07;
        taug(lw_rad.ngs3+11,lay,m,n)*=1.1;
        taug(lw_rad.ngs3+12,lay,m,n)*=0.99;
        taug(lw_rad.ngs3+13,lay,m,n)*=0.88;
        taug(lw_rad.ngs3+14,lay,m,n)*=0.943;
    }
}
KOKKOS_FUNCTION
void taugb5(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_05,const real2d_u &fracrefb_05,const real2d_u &absa_05,const real4d_u &ka_05,const real2d_u &absb_05,const real4d_u &kb_05,const real3d_u &ka_mo3_05,const real2d_u &selfref_05,const real2d_u &forref_05,const real1d_u &ccl4_05,int m,int n){
    /*
        band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
                            (high key - o3,co2)
    */
    int lay, ind0, ind1, inds, indf, indm, ig;
    int js, js1, jmo3, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_mo3, specparm_mo3, specmult_mo3, fmo3;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor, o3m1, o3m2, abso3;
    real refrat_planck_a, refrat_planck_b, refrat_m_a;
    real tau_major, tau_major1;
    /*
        ! Minor gas mapping level :
        !     lower - o3, p = 317.34 mbar, t = 240.77 k
        !     lower - ccl4

        ! Calculate reference ratio to be used in calculation of Planck
        ! fraction in lower/upper atmosphere.
    */
    // P = 473.420 mb
    refrat_planck_a = lw_rad.chi_mls[0][4]/lw_rad.chi_mls[1][4];
    // P = 0.2369 mb
    refrat_planck_b = lw_rad.chi_mls[2][42]/lw_rad.chi_mls[1][42];
    // P = 317.3480
    refrat_m_a = lw_rad.chi_mls[0][6]/lw_rad.chi_mls[1][6];
    /*
        ! Compute the optical depth by interpolating in ln(pressure) and 
        ! temperature, and appropriate species.  Below laytrop, the 
        ! water vapor self-continuum and foreign continuum is 
        ! interpolated (in temperature) separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_mo3 = colh2o(lay) + refrat_m_a*colco2(lay);
        specparm_mo3 = colh2o(lay)/speccomb_mo3;
        if (specparm_mo3 >= lw_rad.oneminus) specparm_mo3 = lw_rad.oneminus;
        specmult_mo3 = 8.*specparm_mo3;
        jmo3 = 1 + int(specmult_mo3);
        fmo3 = fmod(specmult_mo3,1.0);

        speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[4] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[4] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);

        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }

        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1 ;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng5; ig++){
            tauself = selffac(lay) * (selfref_05(ig,inds) + selffrac(lay) * 
                 (selfref_05(ig,inds+1) - selfref_05(ig,inds)));
            taufor =  forfac(lay) * (forref_05(ig,indf) + forfrac(lay) * 
                 (forref_05(ig,indf+1) - forref_05(ig,indf)));
            o3m1 = ka_mo3_05(ig,jmo3,indm) + fmo3 * 
                 (ka_mo3_05(ig,jmo3+1,indm)-ka_mo3_05(ig,jmo3,indm));
            o3m2 = ka_mo3_05(ig,jmo3,indm+1) + fmo3 * 
                 (ka_mo3_05(ig,jmo3+1,indm+1)-ka_mo3_05(ig,jmo3,indm+1));
            abso3 = o3m1 + minorfrac(lay)*(o3m2-o3m1);
            if (specparm < 0.125){
               tau_major = speccomb * 
                    (fac000 * absa_05(ig,ind0) + 
                    fac100 * absa_05(ig,ind0+1) + 
                    fac200 * absa_05(ig,ind0+2) + 
                    fac010 * absa_05(ig,ind0+9) + 
                    fac110 * absa_05(ig,ind0+10) + 
                    fac210 * absa_05(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * 
                    (fac200 * absa_05(ig,ind0-1) + 
                    fac100 * absa_05(ig,ind0) + 
                    fac000 * absa_05(ig,ind0+1) + 
                    fac210 * absa_05(ig,ind0+8) + 
                    fac110 * absa_05(ig,ind0+9) + 
                    fac010 * absa_05(ig,ind0+10));
            }else{
               tau_major = speccomb * 
                    (fac000 * absa_05(ig,ind0) + 
                    fac100 * absa_05(ig,ind0+1) + 
                    fac010 * absa_05(ig,ind0+9) + 
                    fac110 * absa_05(ig,ind0+10));
            }

            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * 
                    (fac001 * absa_05(ig,ind1) + 
                    fac101 * absa_05(ig,ind1+1) + 
                    fac201 * absa_05(ig,ind1+2) + 
                    fac011 * absa_05(ig,ind1+9) + 
                    fac111 * absa_05(ig,ind1+10) + 
                    fac211 * absa_05(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * 
                    (fac201 * absa_05(ig,ind1-1) + 
                    fac101 * absa_05(ig,ind1) + 
                    fac001 * absa_05(ig,ind1+1) + 
                    fac211 * absa_05(ig,ind1+8) + 
                    fac111 * absa_05(ig,ind1+9) + 
                    fac011 * absa_05(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * 
                    (fac001 * absa_05(ig,ind1) +
                    fac101 * absa_05(ig,ind1+1) + 
                    fac011 * absa_05(ig,ind1+9) + 
                    fac111 * absa_05(ig,ind1+10));
            }
            taug(lw_rad.ngs4+ig,lay,m,n) = tau_major + tau_major1 
                 + tauself + taufor 
                 + abso3*colo3(lay) 
                 + wx(1,lay,m,n) * ccl4_05(ig);
            fracs(lw_rad.ngs4+ig,lay,m,n) = fracrefa_05(ig,jpl) + fpl * 
                 (fracrefa_05(ig,jpl+1)-fracrefa_05(ig,jpl));
        }
    }

    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay);
        specparm = colo3(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 4.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay);
        specparm1 = colo3(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 4.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        fac000 = (1. - fs) * fac00(lay);
        fac010 = (1. - fs) * fac10(lay);
        fac100 = fs * fac00(lay);
        fac110 = fs * fac10(lay);
        fac001 = (1. - fs1) * fac01(lay);
        fac011 = (1. - fs1) * fac11(lay);
        fac101 = fs1 * fac01(lay);
        fac111 = fs1 * fac11(lay);

        speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay);
        specparm_planck = colo3(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 4.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[4] + js;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[4] + js1;
        for(int ig = 1; ig <= lw_rad.ng5; ig++){
            taug(lw_rad.ngs4+ig,lay,m,n) = speccomb * 
                (fac000 * absb_05(ig,ind0) + 
                fac100 * absb_05(ig,ind0+1) + 
                fac010 * absb_05(ig,ind0+5) + 
                fac110 * absb_05(ig,ind0+6)) 
                + speccomb1 * 
                (fac001 * absb_05(ig,ind1) + 
                fac101 * absb_05(ig,ind1+1) + 
                fac011 * absb_05(ig,ind1+5) + 
                fac111 * absb_05(ig,ind1+6))  
                + wx(1,lay,m,n) * ccl4_05(ig);
            fracs(lw_rad.ngs4+ig,lay,m,n) = fracrefb_05(ig,jpl) + fpl * 
                (fracrefb_05(ig,jpl+1)-fracrefb_05(ig,jpl));
        }
    }
}
KOKKOS_FUNCTION
void taugb6(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real1d_u &fracrefa_06,const real2d_u &absa_06,const real3d_u &ka_06,const real2d_u &ka_mco2_06,const real2d_u &selfref_06,const real2d_u &forref_06,const real1d_u &cfc11adj_06,const real1d_u &cfc12_06,int m,int n){
    /*
        band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
                            (high key - nothing; high minor - cfc11, cfc12)            
    */
    int lay, ind0, ind1, inds, indf, indm, ig;
    real chi_co2, ratco2, adjfac, adjcolco2;
    real tauself, taufor, absco2;
    /*
        ! Minor gas mapping level:
        !     lower - co2, p = 706.2720 mb, t = 294.2 k
        !     upper - cfc11, cfc12

        ! Compute the optical depth by interpolating in ln(pressure) and
        ! temperature. The water vapor self-continuum and foreign continuum
        ! is interpolated (in temperature) separately.  

        ! Lower atmosphere loop
    */
    for(int lay = 1; lay <= laytrop; lay++){
        /*
            ! In atmospheres where the amount of CO2 is too great to be considered
            ! a minor species, adjust the column amount of CO2 by an empirical factor 
            ! to obtain the proper contribution.
        */
        chi_co2 = colco2(lay)/(coldry(lay));
        ratco2 = 1.e20*chi_co2/lw_rad.chi_mls[1][jp(lay)];
        if (ratco2 > 3.0){
            adjfac = 2.0+std::pow((ratco2-2.0),0.77);
            adjcolco2 = adjfac*lw_rad.chi_mls[1][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcolco2 = colco2(lay);
        }

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[5] + 1;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[5] + 1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);
        for(int ig = 1; ig <= lw_rad.ng6; ig++){
            tauself = selffac(lay) * (selfref_06(ig,inds) + selffrac(lay) * 
                 (selfref_06(ig,inds+1) - selfref_06(ig,inds)));
            taufor =  forfac(lay) * (forref_06(ig,indf) + forfrac(lay) * 
                 (forref_06(ig,indf+1) - forref_06(ig,indf)));
            absco2 =  (ka_mco2_06(ig,indm) + minorfrac(lay) * 
                 (ka_mco2_06(ig,indm+1) - ka_mco2_06(ig,indm)));
            taug(lw_rad.ngs5+ig,lay,m,n) = colh2o(lay) * 
                (fac00(lay) * absa_06(ig,ind0) + 
                 fac10(lay) * absa_06(ig,ind0+1) + 
                 fac01(lay) * absa_06(ig,ind1) +  
                 fac11(lay) * absa_06(ig,ind1+1))  
                 + tauself + taufor 
                 + adjcolco2 * absco2 
                 + wx(2,lay,m,n) * cfc11adj_06(ig) 
                 + wx(3,lay,m,n) * cfc12_06(ig);
            fracs(lw_rad.ngs5+ig,lay,m,n) = fracrefa_06(ig);
        }
    }
    // Upper atmosphere loop
    // Nothing important goes on above laytrop in this band.
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        for(int ig = 1; ig <= lw_rad.ng6; ig++){
            taug(lw_rad.ngs5+ig,lay,m,n) = 0.0
                 + wx(2,lay,m,n) * cfc11adj_06(ig) 
                 + wx(3,lay,m,n) * cfc12_06(ig);
            fracs(lw_rad.ngs5+ig,lay,m,n) = fracrefa_06(ig);
        }
    }
}
KOKKOS_FUNCTION
void taugb7(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_07,const real1d_u &fracrefb_07,const real2d_u &absa_07,const real4d_u &ka_07,const real2d_u &absb_07,const real3d_u &kb_07,const real3d_u &ka_mco2_07,const real2d_u &kb_mco2_07,const real2d_u &selfref_07,const real2d_u &forref_07,int m,int n){
    /*
        band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
                                (high key - o3; high minor - co2)
    */
    int  lay, ind0, ind1, inds, indf, indm, ig;
    int js, js1, jmco2, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_mco2, specparm_mco2, specmult_mco2, fmco2;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor, co2m1, co2m2, absco2;
    real chi_co2, ratco2, adjfac, adjcolco2;
    real refrat_planck_a, refrat_m_a;
    real tau_major, tau_major1;
    /*
        ! Minor gas mapping level :
        !     lower - co2, p = 706.2620 mbar, t= 278.94 k
        !     upper - co2, p = 12.9350 mbar, t = 234.01 k

        ! Calculate reference ratio to be used in calculation of Planck
        ! fraction in lower atmosphere
    */
    // P = 706.2620 mb
    refrat_planck_a = lw_rad.chi_mls[0][2]/lw_rad.chi_mls[2][2];
    // P = 706.2720 mb
    refrat_m_a = lw_rad.chi_mls[0][2]/lw_rad.chi_mls[2][2];
    /*
        ! Compute the optical depth by interpolating in ln(pressure), 
        ! temperature, and appropriate species.  Below laytrop, the water
        ! vapor self-continuum and foreign continuum is interpolated 
        ! (in temperature) separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2oo3(lay)*colo3(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colh2o(lay) + rat_h2oo3_1(lay)*colo3(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_mco2 = colh2o(lay) + refrat_m_a*colo3(lay);
        specparm_mco2 = colh2o(lay)/speccomb_mco2;
        if (specparm_mco2 >= lw_rad.oneminus) specparm_mco2 = lw_rad.oneminus;
        specmult_mco2 = 8.*specparm_mco2;

        jmco2 = 1 + int(specmult_mco2);
        fmco2 = fmod(specmult_mco2,1.0);
        /*
            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_co2 = colco2(lay)/(coldry(lay));
        ratco2 = 1.e20*chi_co2/lw_rad.chi_mls[1][jp(lay)];
        if (ratco2 > 3.0){
            adjfac = 3.0+std::pow((ratco2-3.0),0.79);
            adjcolco2 = adjfac*lw_rad.chi_mls[1][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcolco2 = colco2(lay);
        }

        speccomb_planck = colh2o(lay)+refrat_planck_a*colo3(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[6] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[6] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);

        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }
        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng7; ig++){
            tauself = selffac(lay)* (selfref_07(ig,inds) + selffrac(lay) * 
                 (selfref_07(ig,inds+1) - selfref_07(ig,inds)));
            taufor = forfac(lay) * (forref_07(ig,indf) + forfrac(lay) * 
                 (forref_07(ig,indf+1) - forref_07(ig,indf)));
            co2m1 = ka_mco2_07(ig,jmco2,indm) + fmco2 * 
                 (ka_mco2_07(ig,jmco2+1,indm) - ka_mco2_07(ig,jmco2,indm));
            co2m2 = ka_mco2_07(ig,jmco2,indm+1) + fmco2 * 
                 (ka_mco2_07(ig,jmco2+1,indm+1) - ka_mco2_07(ig,jmco2,indm+1));
            absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1);
            if (specparm < 0.125){
               tau_major = speccomb * 
                    (fac000 * absa_07(ig,ind0) + 
                    fac100 * absa_07(ig,ind0+1) + 
                    fac200 * absa_07(ig,ind0+2) + 
                    fac010 * absa_07(ig,ind0+9) + 
                    fac110 * absa_07(ig,ind0+10) + 
                    fac210 * absa_07(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * 
                    (fac200 * absa_07(ig,ind0-1) + 
                    fac100 * absa_07(ig,ind0) + 
                    fac000 * absa_07(ig,ind0+1) +
                    fac210 * absa_07(ig,ind0+8) + 
                    fac110 * absa_07(ig,ind0+9) + 
                    fac010 * absa_07(ig,ind0+10));
            }else{
               tau_major = speccomb * 
                    (fac000 * absa_07(ig,ind0) + 
                    fac100 * absa_07(ig,ind0+1) + 
                    fac010 * absa_07(ig,ind0+9) + 
                    fac110 * absa_07(ig,ind0+10));
            }

            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * 
                    (fac001 * absa_07(ig,ind1) + 
                    fac101 * absa_07(ig,ind1+1) + 
                    fac201 * absa_07(ig,ind1+2) + 
                    fac011 * absa_07(ig,ind1+9) + 
                    fac111 * absa_07(ig,ind1+10) + 
                    fac211 * absa_07(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * 
                    (fac201 * absa_07(ig,ind1-1) + 
                    fac101 * absa_07(ig,ind1) + 
                    fac001 * absa_07(ig,ind1+1) + 
                    fac211 * absa_07(ig,ind1+8) + 
                    fac111 * absa_07(ig,ind1+9) + 
                    fac011 * absa_07(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * 
                    (fac001 * absa_07(ig,ind1) +  
                    fac101 * absa_07(ig,ind1+1) + 
                    fac011 * absa_07(ig,ind1+9) + 
                    fac111 * absa_07(ig,ind1+10));
            }
            taug(lw_rad.ngs6+ig,lay,m,n) = tau_major + tau_major1 
                 + tauself + taufor 
                 + adjcolco2*absco2;
            fracs(lw_rad.ngs6+ig,lay,m,n) = fracrefa_07(ig,jpl) + fpl * 
                 (fracrefa_07(ig,jpl+1)-fracrefa_07(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        /*
            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_co2 = colco2(lay)/(coldry(lay));
        ratco2 = 1.e20*chi_co2/lw_rad.chi_mls[1][jp(lay)];
        if (ratco2 > 3.0){
            adjfac = 2.0+std::pow((ratco2-2.0),0.79);
            adjcolco2 = adjfac*lw_rad.chi_mls[1][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcolco2 = colco2(lay);
        }

        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[6] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[6] + 1;
        indm = indminor(lay);
        for(int ig = 1; ig <= lw_rad.ng7; ig++){
            absco2 = kb_mco2_07(ig,indm) + minorfrac(lay) * 
                 (kb_mco2_07(ig,indm+1) - kb_mco2_07(ig,indm));
            taug(lw_rad.ngs6+ig,lay,m,n) = colo3(lay) * 
                 (fac00(lay) * absb_07(ig,ind0) +
                 fac10(lay) * absb_07(ig,ind0+1) + 
                 fac01(lay) * absb_07(ig,ind1) + 
                 fac11(lay) * absb_07(ig,ind1+1))
                 + adjcolco2 * absco2;
            fracs(lw_rad.ngs6+ig,lay,m,n) = fracrefb_07(ig);
        }
        /*
            ! Empirical modification to code to improve stratospheric cooling rates
            ! for o3.  Revised to apply weighting for g-point reduction in this band.
        */
        taug(lw_rad.ngs6+6,lay,m,n)=taug(lw_rad.ngs6+6,lay,m,n)*0.92;
        taug(lw_rad.ngs6+7,lay,m,n)=taug(lw_rad.ngs6+7,lay,m,n)*0.88;
        taug(lw_rad.ngs6+8,lay,m,n)=taug(lw_rad.ngs6+8,lay,m,n)*1.07;
        taug(lw_rad.ngs6+9,lay,m,n)=taug(lw_rad.ngs6+9,lay,m,n)*1.1;
        taug(lw_rad.ngs6+10,lay,m,n)=taug(lw_rad.ngs6+10,lay,m,n)*0.99;
        taug(lw_rad.ngs6+11,lay,m,n)=taug(lw_rad.ngs6+11,lay,m,n)*0.855;
    }
}
KOKKOS_FUNCTION
void taugb8(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real1d_u &fracrefa_08,const real1d_u &fracrefb_08,const real2d_u &absa_08,const real3d_u &ka_08,const real2d_u &absb_08,const real3d_u &kb_08,const real2d_u &ka_mco2_08,const real2d_u &ka_mn2o_08,const real2d_u &ka_mo3_08,const real2d_u &kb_mco2_08,const real2d_u &kb_mn2o_08,    
    const real2d_u &selfref_08,const real2d_u &forref_08,const real1d_u &cfc12_08,const real1d_u &cfc22adj_08,int m,int n){
    /*
        band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
                                (high key - o3; high minor - co2, n2o)
    */
    int lay, ind0, ind1, inds, indf, indm, ig;
    real tauself, taufor, absco2, abso3, absn2o;
    real chi_co2, ratco2, adjfac, adjcolco2;
    /*
        ! Minor gas mapping level:
        !     lower - co2, p = 1053.63 mb, t = 294.2 k
        !     lower - o3,  p = 317.348 mb, t = 240.77 k
        !     lower - n2o, p = 706.2720 mb, t= 278.94 k
        !     lower - cfc12,cfc11
        !     upper - co2, p = 35.1632 mb, t = 223.28 k
        !     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

        ! Compute the optical depth by interpolating in ln(pressure) and 
        ! temperature, and appropriate species.  Below laytrop, the water vapor 
        ! self-continuum and foreign continuum is interpolated (in temperature) 
        ! separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        /*
            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_co2 = colco2(lay)/(coldry(lay));
        ratco2 = 1.e20*chi_co2/lw_rad.chi_mls[1][jp(lay)];
        if (ratco2 > 3.0){
            adjfac = 2.0+std::pow((ratco2-2.0),0.65);
            adjcolco2 = adjfac*lw_rad.chi_mls[1][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcolco2 = colco2(lay);
        }

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[7] + 1;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[7] + 1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);
        for(int ig = 1; ig <= lw_rad.ng8; ig++){
            tauself = selffac(lay) * (selfref_08(ig,inds) + selffrac(lay) * 
                 (selfref_08(ig,inds+1) - selfref_08(ig,inds)));
            taufor = forfac(lay) * (forref_08(ig,indf) + forfrac(lay) * 
                 (forref_08(ig,indf+1) - forref_08(ig,indf)));
            absco2 =  (ka_mco2_08(ig,indm) + minorfrac(lay) * 
                 (ka_mco2_08(ig,indm+1) - ka_mco2_08(ig,indm)));
            abso3 =  (ka_mo3_08(ig,indm) + minorfrac(lay) * 
                 (ka_mo3_08(ig,indm+1) - ka_mo3_08(ig,indm)));
            absn2o =  (ka_mn2o_08(ig,indm) + minorfrac(lay) * 
                 (ka_mn2o_08(ig,indm+1) - ka_mn2o_08(ig,indm)));
            taug(lw_rad.ngs7+ig,lay,m,n) = colh2o(lay) * 
                 (fac00(lay) * absa_08(ig,ind0) + 
                 fac10(lay) * absa_08(ig,ind0+1) + 
                 fac01(lay) * absa_08(ig,ind1) +  
                 fac11(lay) * absa_08(ig,ind1+1)) 
                 + tauself + taufor 
                 + adjcolco2*absco2 
                 + colo3(lay) * abso3 
                 + coln2o(lay) * absn2o 
                 + wx(3,lay,m,n) * cfc12_08(ig) 
                 + wx(4,lay,m,n) * cfc22adj_08(ig);
            fracs(lw_rad.ngs7+ig,lay,m,n) = fracrefa_08(ig);
        }
    }

    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        /*
            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_co2 = colco2(lay)/coldry(lay);
        ratco2 = 1.e20*chi_co2/lw_rad.chi_mls[1][jp(lay)];
        if (ratco2 > 3.0){
            adjfac = 2.0+std::pow((ratco2-2.0),0.65);
            adjcolco2 = adjfac*lw_rad.chi_mls[1][jp(lay)] * coldry(lay)*1.e-20;
        }else{
            adjcolco2 = colco2(lay);
        }

        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[7] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[7] + 1;
        indm = indminor(lay);
        for(int ig = 1; ig <= lw_rad.ng8; ig++){
            absco2 =  (kb_mco2_08(ig,indm) + minorfrac(lay) * 
                 (kb_mco2_08(ig,indm+1) - kb_mco2_08(ig,indm)));
            absn2o =  (kb_mn2o_08(ig,indm) + minorfrac(lay) * 
                 (kb_mn2o_08(ig,indm+1) - kb_mn2o_08(ig,indm)));
            taug(lw_rad.ngs7+ig,lay,m,n) = colo3(lay) * 
                 (fac00(lay) * absb_08(ig,ind0) + 
                 fac10(lay) * absb_08(ig,ind0+1) + 
                 fac01(lay) * absb_08(ig,ind1) + 
                 fac11(lay) * absb_08(ig,ind1+1)) 
                 + adjcolco2*absco2 
                 + coln2o(lay)*absn2o 
                 + wx(3,lay,m,n) * cfc12_08(ig) 
                 + wx(4,lay,m,n) * cfc22adj_08(ig);
            fracs(lw_rad.ngs7+ig,lay,m,n) = fracrefb_08(ig);
        }
    }
}
KOKKOS_FUNCTION
void taugb9(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_09,const real1d_u &fracrefb_09,const real2d_u &absa_09,const real4d_u &ka_09,const real2d_u &absb_09,const real3d_u &kb_09,      
    const real3d_u &ka_mn2o_09,const real2d_u &kb_mn2o_09,const real2d_u &selfref_09,const real2d_u &forref_09,int m,int n){
    /*
        band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
                                (high key - ch4; high minor - n2o)
    */
    int lay, ind0, ind1, inds, indf, indm, ig;
    int js, js1, jmn2o, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor, n2om1, n2om2, absn2o;
    real chi_n2o, ratn2o, adjfac, adjcoln2o;
    real refrat_planck_a, refrat_m_a;
    real tau_major, tau_major1;
    /*
        ! Minor gas mapping level :
        !     lower - n2o, p = 706.272 mbar, t = 278.94 k
        !     upper - n2o, p = 95.58 mbar, t = 215.7 k

        ! Calculate reference ratio to be used in calculation of Planck
        ! fraction in lower/upper atmosphere.
    */
    // P = 212 mb
    refrat_planck_a = lw_rad.chi_mls[0][8]/lw_rad.chi_mls[5][8];
    // P = 706.272 mb 
    refrat_m_a = lw_rad.chi_mls[0][2]/lw_rad.chi_mls[5][2];
    /*
        ! Compute the optical depth by interpolating in ln(pressure), 
        ! temperature, and appropriate species.  Below laytrop, the water
        ! vapor self-continuum and foreign continuum is interpolated 
        ! (in temperature) separately.  
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_mn2o = colh2o(lay) + refrat_m_a*colch4(lay);
        specparm_mn2o = colh2o(lay)/speccomb_mn2o;
        if (specparm_mn2o >= lw_rad.oneminus) specparm_mn2o = lw_rad.oneminus;
        specmult_mn2o = 8.*specparm_mn2o;
        jmn2o = 1 + int(specmult_mn2o);
        fmn2o = fmod(specmult_mn2o,1.0);
        /*
            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_n2o = coln2o(lay)/(coldry(lay));
        ratn2o = 1.e20*chi_n2o/lw_rad.chi_mls[3][jp(lay)];
        if (ratn2o > 1.5){
            adjfac = 0.5+std::pow((ratn2o-0.5),0.65);
            adjcoln2o = adjfac*lw_rad.chi_mls[3][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcoln2o = coln2o(lay);
        }

        speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[8] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[8] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);
        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }

        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng9; ig++){
            tauself = selffac(lay)* (selfref_09(ig,inds) + selffrac(lay) * 
                 (selfref_09(ig,inds+1) - selfref_09(ig,inds)));
            taufor = forfac(lay) * (forref_09(ig,indf) + forfrac(lay) * 
                 (forref_09(ig,indf+1) - forref_09(ig,indf)));
            n2om1 = ka_mn2o_09(ig,jmn2o,indm) + fmn2o * 
                 (ka_mn2o_09(ig,jmn2o+1,indm) - ka_mn2o_09(ig,jmn2o,indm));
            n2om2 = ka_mn2o_09(ig,jmn2o,indm+1) + fmn2o * 
                 (ka_mn2o_09(ig,jmn2o+1,indm+1) - ka_mn2o_09(ig,jmn2o,indm+1));
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1);
            if (specparm < 0.125){
               tau_major = speccomb * 
                    (fac000 * absa_09(ig,ind0) + 
                    fac100 * absa_09(ig,ind0+1) + 
                    fac200 * absa_09(ig,ind0+2) + 
                    fac010 * absa_09(ig,ind0+9) + 
                    fac110 * absa_09(ig,ind0+10) + 
                    fac210 * absa_09(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * 
                    (fac200 * absa_09(ig,ind0-1) + 
                    fac100 * absa_09(ig,ind0) + 
                    fac000 * absa_09(ig,ind0+1) + 
                    fac210 * absa_09(ig,ind0+8) + 
                    fac110 * absa_09(ig,ind0+9) + 
                    fac010 * absa_09(ig,ind0+10));
            }else{
               tau_major = speccomb * 
                    (fac000 * absa_09(ig,ind0) + 
                    fac100 * absa_09(ig,ind0+1) + 
                    fac010 * absa_09(ig,ind0+9) + 
                    fac110 * absa_09(ig,ind0+10));
            }

            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * 
                    (fac001 * absa_09(ig,ind1) + 
                    fac101 * absa_09(ig,ind1+1) + 
                    fac201 * absa_09(ig,ind1+2) + 
                    fac011 * absa_09(ig,ind1+9) + 
                    fac111 * absa_09(ig,ind1+10) + 
                    fac211 * absa_09(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * 
                    (fac201 * absa_09(ig,ind1-1) + 
                    fac101 * absa_09(ig,ind1) + 
                    fac001 * absa_09(ig,ind1+1) + 
                    fac211 * absa_09(ig,ind1+8) + 
                    fac111 * absa_09(ig,ind1+9) + 
                    fac011 * absa_09(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * 
                    (fac001 * absa_09(ig,ind1) + 
                    fac101 * absa_09(ig,ind1+1) + 
                    fac011 * absa_09(ig,ind1+9) + 
                    fac111 * absa_09(ig,ind1+10));
            }
            taug(lw_rad.ngs8+ig,lay,m,n) = tau_major + tau_major1 
                 + tauself + taufor 
                 + adjcoln2o*absn2o;
            fracs(lw_rad.ngs8+ig,lay,m,n) = fracrefa_09(ig,jpl) + fpl * 
                 (fracrefa_09(ig,jpl+1)-fracrefa_09(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        /*
            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_n2o = coln2o(lay)/(coldry(lay));
        ratn2o = 1.e20*chi_n2o/lw_rad.chi_mls[3][jp(lay)];
        if (ratn2o > 1.5){
            adjfac = 0.5+std::pow((ratn2o-0.5),0.65);
            adjcoln2o = adjfac*lw_rad.chi_mls[3][jp(lay)]*coldry(lay)*1.e-20;
        }else{
            adjcoln2o = coln2o(lay);
        }

        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[8] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[8] + 1;
        indm = indminor(lay);

        for(int ig = 1; ig <= lw_rad.ng9; ig++){
            absn2o = kb_mn2o_09(ig,indm) + minorfrac(lay) * 
                (kb_mn2o_09(ig,indm+1) - kb_mn2o_09(ig,indm));
            taug(lw_rad.ngs8+ig,lay,m,n) = colch4(lay) * 
                 (fac00(lay) * absb_09(ig,ind0) + 
                 fac10(lay) * absb_09(ig,ind0+1) + 
                 fac01(lay) * absb_09(ig,ind1) +  
                 fac11(lay) * absb_09(ig,ind1+1))
                 + adjcoln2o*absn2o;
            fracs(lw_rad.ngs8+ig,lay,m,n) = fracrefb_09(ig);
        }
    }
}
KOKKOS_FUNCTION
void taugb10(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real1d_u &fracrefa_10,const real1d_u &fracrefb_10,const real2d_u &absa_10,const real3d_u &ka_10,const real2d_u &absb_10,const real3d_u &kb_10,      
    const real2d_u &selfref_10,const real2d_u &forref_10,int m,int n){
    // band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
    int lay, ind0, ind1, inds, indf, ig;
    real tauself, taufor;
    /*
        ! Compute the optical depth by interpolating in ln(pressure) and 
        ! temperature.  Below laytrop, the water vapor self-continuum and
        ! foreign continuum is interpolated (in temperature) separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[9] + 1;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[9] + 1;
        inds = indself(lay);
        indf = indfor(lay);
        for(int ig = 1; ig <= lw_rad.ng10; ig++){
            tauself = selffac(lay) * (selfref_10(ig,inds) + selffrac(lay) * 
                 (selfref_10(ig,inds+1) - selfref_10(ig,inds)));
            taufor = forfac(lay) * (forref_10(ig,indf) + forfrac(lay) * 
                 (forref_10(ig,indf+1) - forref_10(ig,indf)));
            taug(lw_rad.ngs9+ig,lay,m,n) = colh2o(lay) * 
                 (fac00(lay) * absa_10(ig,ind0) + 
                 fac10(lay) * absa_10(ig,ind0+1) + 
                 fac01(lay) * absa_10(ig,ind1) + 
                 fac11(lay) * absa_10(ig,ind1+1))
                 + tauself + taufor;
            fracs(lw_rad.ngs9+ig,lay,m,n) = fracrefa_10(ig);
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[9] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[9] + 1;
        indf = indfor(lay);
        for(int ig = 1; ig <= lw_rad.ng10; ig++){
            taufor = forfac(lay) * (forref_10(ig,indf) + forfrac(lay) * 
                 (forref_10(ig,indf+1) - forref_10(ig,indf))); 
            taug(lw_rad.ngs9+ig,lay,m,n) = colh2o(lay) * 
                 (fac00(lay) * absb_10(ig,ind0) + 
                 fac10(lay) * absb_10(ig,ind0+1) + 
                 fac01(lay) * absb_10(ig,ind1) +  
                 fac11(lay) * absb_10(ig,ind1+1))
                 + taufor;
            fracs(lw_rad.ngs9+ig,lay,m,n) = fracrefb_10(ig);
        }
    }
}
KOKKOS_FUNCTION
void taugb11(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real1d_u &fracrefa_11,const real1d_u &fracrefb_11,const real2d_u &absa_11,const real3d_u &ka_11,const real2d_u &absb_11,const real3d_u &kb_11,      
    const real2d_u &ka_mo2_11,const real2d_u &kb_mo2_11,const real2d_u &selfref_11,const real2d_u &forref_11,int m,int n){
    /*
        band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
                                 (high key - h2o; high minor - o2)
    */
    int lay, ind0, ind1, inds, indf, indm, ig;
    real scaleo2, tauself, taufor, tauo2;
    /*
        ! Minor gas mapping level :
        !     lower - o2, p = 706.2720 mbar, t = 278.94 k
        !     upper - o2, p = 4.758820 mbarm t = 250.85 k

        ! Compute the optical depth by interpolating in ln(pressure) and 
        ! temperature.  Below laytrop, the water vapor self-continuum and
        ! foreign continuum is interpolated (in temperature) separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[10] + 1;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[10] + 1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);
        scaleo2 = colo2(lay)*scaleminor(lay);
        for(int ig = 1; ig <= lw_rad.ng11; ig++){
            tauself = selffac(lay) * (selfref_11(ig,inds) + selffrac(lay) * 
                 (selfref_11(ig,inds+1) - selfref_11(ig,inds)));
            taufor = forfac(lay) * (forref_11(ig,indf) + forfrac(lay) * 
                 (forref_11(ig,indf+1) - forref_11(ig,indf)));
            tauo2 =  scaleo2 * (ka_mo2_11(ig,indm) + minorfrac(lay) * 
                 (ka_mo2_11(ig,indm+1) - ka_mo2_11(ig,indm)));
            taug(lw_rad.ngs10+ig,lay,m,n) = colh2o(lay) * 
                 (fac00(lay) * absa_11(ig,ind0) + 
                 fac10(lay) * absa_11(ig,ind0+1) + 
                 fac01(lay) * absa_11(ig,ind1) + 
                 fac11(lay) * absa_11(ig,ind1+1)) 
                 + tauself + taufor 
                 + tauo2;
            fracs(lw_rad.ngs10+ig,lay,m,n) = fracrefa_11(ig);
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[10] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[10] + 1;
        indf = indfor(lay);
        indm = indminor(lay);
        scaleo2 = colo2(lay)*scaleminor(lay);
        for(int ig = 1; ig <= lw_rad.ng11; ig++){
            taufor = forfac(lay) * (forref_11(ig,indf) + forfrac(lay) * 
                 (forref_11(ig,indf+1) - forref_11(ig,indf)));
            tauo2 =  scaleo2 * (kb_mo2_11(ig,indm) + minorfrac(lay) * 
                 (kb_mo2_11(ig,indm+1) - kb_mo2_11(ig,indm)));
            taug(lw_rad.ngs10+ig,lay,m,n) = colh2o(lay) * 
                 (fac00(lay) * absb_11(ig,ind0) + 
                 fac10(lay) * absb_11(ig,ind0+1) + 
                 fac01(lay) * absb_11(ig,ind1) + 
                 fac11(lay) * absb_11(ig,ind1+1))
                 + taufor 
                 + tauo2;
            fracs(lw_rad.ngs10+ig,lay,m,n) = fracrefb_11(ig);
        }
    }
}
KOKKOS_FUNCTION
void taugb12(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_12,const real2d_u &absa_12,const real4d_u &ka_12,const real2d_u &selfref_12,const real2d_u &forref_12,int m,int n){
    // band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
    int lay, ind0, ind1, inds, indf, ig;
    int js, js1, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor;
    real refrat_planck_a;
    real tau_major, tau_major1;
    /*
        ! Calculate reference ratio to be used in calculation of Planck
        ! fraction in lower/upper atmosphere.
    */
    // P =   174.164 mb 
    refrat_planck_a = lw_rad.chi_mls[0][9]/lw_rad.chi_mls[1][9];
    /*
        ! Compute the optical depth by interpolating in ln(pressure), 
        ! temperature, and appropriate species.  Below laytrop, the water
        ! vapor self-continuum adn foreign continuum is interpolated 
        ! (in temperature) separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[11] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[11] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }

        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1 ;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng12; ig++){
            tauself = selffac(lay)* (selfref_12(ig,inds) + selffrac(lay) * 
                 (selfref_12(ig,inds+1) - selfref_12(ig,inds)));
            taufor = forfac(lay) * (forref_12(ig,indf) + forfrac(lay) * 
                 (forref_12(ig,indf+1) - forref_12(ig,indf)));
            if (specparm < 0.125){
                tau_major = speccomb * 
                    (fac000 * absa_12(ig,ind0) + 
                    fac100 * absa_12(ig,ind0+1) + 
                    fac200 * absa_12(ig,ind0+2) + 
                    fac010 * absa_12(ig,ind0+9) + 
                    fac110 * absa_12(ig,ind0+10) + 
                    fac210 * absa_12(ig,ind0+11));
            }else if (specparm > 0.875){
                tau_major = speccomb * 
                    (fac200 * absa_12(ig,ind0-1) + 
                    fac100 * absa_12(ig,ind0) + 
                    fac000 * absa_12(ig,ind0+1) + 
                    fac210 * absa_12(ig,ind0+8) + 
                    fac110 * absa_12(ig,ind0+9) + 
                    fac010 * absa_12(ig,ind0+10));
            }else{
                tau_major = speccomb * 
                    (fac000 * absa_12(ig,ind0) + 
                    fac100 * absa_12(ig,ind0+1) + 
                    fac010 * absa_12(ig,ind0+9) + 
                    fac110 * absa_12(ig,ind0+10));
            }

            if (specparm1 < 0.125){
                tau_major1 = speccomb1 * 
                    (fac001 * absa_12(ig,ind1) + 
                    fac101 * absa_12(ig,ind1+1) + 
                    fac201 * absa_12(ig,ind1+2) + 
                    fac011 * absa_12(ig,ind1+9) + 
                    fac111 * absa_12(ig,ind1+10) + 
                    fac211 * absa_12(ig,ind1+11));
            }else if (specparm1 > 0.875){
                tau_major1 = speccomb1 * 
                    (fac201 * absa_12(ig,ind1-1) + 
                    fac101 * absa_12(ig,ind1) + 
                    fac001 * absa_12(ig,ind1+1) + 
                    fac211 * absa_12(ig,ind1+8) + 
                    fac111 * absa_12(ig,ind1+9) + 
                    fac011 * absa_12(ig,ind1+10));
            }else{
                tau_major1 = speccomb1 * 
                    (fac001 * absa_12(ig,ind1) + 
                    fac101 * absa_12(ig,ind1+1) + 
                    fac011 * absa_12(ig,ind1+9) + 
                    fac111 * absa_12(ig,ind1+10));
            }
            taug(lw_rad.ngs11+ig,lay,m,n) = tau_major + tau_major1 
                 + tauself + taufor;
            fracs(lw_rad.ngs11+ig,lay,m,n) = fracrefa_12(ig,jpl) + fpl * 
                 (fracrefa_12(ig,jpl+1)-fracrefa_12(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        for(int ig = 1; ig <= lw_rad.ng12; ig++){
            taug(lw_rad.ngs11+ig,lay,m,n) = 0.0;
            fracs(lw_rad.ngs11+ig,lay,m,n) = 0.0;
        }
    }
}
KOKKOS_FUNCTION
void taugb13(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_13,const real1d_u &fracrefb_13,const real2d_u &absa_13,const real4d_u &ka_13,const real3d_u &ka_mco2_13,const real3d_u &ka_mco_13,const real2d_u &kb_mo3_13,const real2d_u &selfref_13,const real2d_u &forref_13,int m,int n){
    // band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
    int lay, ind0, ind1, inds, indf, indm, ig;
    int js, js1, jmco2, jmco, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_mco2, specparm_mco2, specmult_mco2, fmco2;
    real speccomb_mco, specparm_mco, specmult_mco, fmco;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor, co2m1, co2m2, absco2; 
    real com1, com2, absco, abso3;
    real chi_co2, ratco2, adjfac, adjcolco2;
    real refrat_planck_a, refrat_m_a, refrat_m_a3;
    real tau_major, tau_major1;
    /*
        ! Minor gas mapping levels :
        !     lower - co2, p = 1053.63 mb, t = 294.2 k
        !     lower - co, p = 706 mb, t = 278.94 k
        !     upper - o3, p = 95.5835 mb, t = 215.7 k

        ! Calculate reference ratio to be used in calculation of Planck
        ! fraction in lower/upper atmosphere.
    */
    // P = 473.420 mb (Level 5)
    refrat_planck_a = lw_rad.chi_mls[0][4]/lw_rad.chi_mls[3][4];
    // P = 1053. (Level 1)
    refrat_m_a = lw_rad.chi_mls[0][0]/lw_rad.chi_mls[3][0];
    // P = 706. (Level 3)
    refrat_m_a3 = lw_rad.chi_mls[0][2]/lw_rad.chi_mls[3][2];
    /*
        ! Compute the optical depth by interpolating in ln(pressure), 
        ! temperature, and appropriate species.  Below laytrop, the water
        ! vapor self-continuum and foreign continuum is interpolated 
        ! (in temperature) separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2on2o(lay)*coln2o(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colh2o(lay) + rat_h2on2o_1(lay)*coln2o(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_mco2 = colh2o(lay) + refrat_m_a*coln2o(lay);
        specparm_mco2 = colh2o(lay)/speccomb_mco2;
        if (specparm_mco2 >= lw_rad.oneminus) specparm_mco2 = lw_rad.oneminus;
        specmult_mco2 = 8.*specparm_mco2;
        jmco2 = 1 + int(specmult_mco2);
        fmco2 = fmod(specmult_mco2,1.0);
        /*
            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor 
            !  to obtain the proper contribution.
        */
        chi_co2 = colco2(lay)/(coldry(lay));
        ratco2 = 1.e20*chi_co2/3.55e-4;
        if (ratco2 > 3.0){
            adjfac = 2.0+std::pow((ratco2-2.0),0.68);
            adjcolco2 = adjfac*3.55e-4*coldry(lay)*1.e-20;
        }else{
            adjcolco2 = colco2(lay);
        }

        speccomb_mco = colh2o(lay) + refrat_m_a3*coln2o(lay);
        specparm_mco = colh2o(lay)/speccomb_mco;
        if (specparm_mco >= lw_rad.oneminus) specparm_mco = lw_rad.oneminus;
        specmult_mco = 8.*specparm_mco;
        jmco = 1 + int(specmult_mco);
        fmco = fmod(specmult_mco,1.0);

        speccomb_planck = colh2o(lay)+refrat_planck_a*coln2o(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[12] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[12] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);
        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }

        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng13; ig++){
            tauself = selffac(lay)* (selfref_13(ig,inds) + selffrac(lay) * 
                 (selfref_13(ig,inds+1) - selfref_13(ig,inds)));
            taufor = forfac(lay) * (forref_13(ig,indf) + forfrac(lay) * 
                 (forref_13(ig,indf+1) - forref_13(ig,indf)));
            co2m1 = ka_mco2_13(ig,jmco2,indm) + fmco2 * 
                 (ka_mco2_13(ig,jmco2+1,indm) - ka_mco2_13(ig,jmco2,indm));
            co2m2 = ka_mco2_13(ig,jmco2,indm+1) + fmco2 * 
                 (ka_mco2_13(ig,jmco2+1,indm+1) - ka_mco2_13(ig,jmco2,indm+1));
            absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1);
            com1 = ka_mco_13(ig,jmco,indm) + fmco * 
                 (ka_mco_13(ig,jmco+1,indm) - ka_mco_13(ig,jmco,indm));
            com2 = ka_mco_13(ig,jmco,indm+1) + fmco * 
                 (ka_mco_13(ig,jmco+1,indm+1) - ka_mco_13(ig,jmco,indm+1));
            absco = com1 + minorfrac(lay) * (com2 - com1);
            if (specparm < 0.125){
               tau_major = speccomb * 
                    (fac000 * absa_13(ig,ind0) + 
                    fac100 * absa_13(ig,ind0+1) + 
                    fac200 * absa_13(ig,ind0+2) + 
                    fac010 * absa_13(ig,ind0+9) + 
                    fac110 * absa_13(ig,ind0+10) + 
                    fac210 * absa_13(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * 
                    (fac200 * absa_13(ig,ind0-1) + 
                    fac100 * absa_13(ig,ind0) + 
                    fac000 * absa_13(ig,ind0+1) + 
                    fac210 * absa_13(ig,ind0+8) + 
                    fac110 * absa_13(ig,ind0+9) + 
                    fac010 * absa_13(ig,ind0+10));
            }else{
               tau_major = speccomb * 
                    (fac000 * absa_13(ig,ind0) + 
                    fac100 * absa_13(ig,ind0+1) + 
                    fac010 * absa_13(ig,ind0+9) + 
                    fac110 * absa_13(ig,ind0+10));
            }
            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * 
                    (fac001 * absa_13(ig,ind1) + 
                    fac101 * absa_13(ig,ind1+1) + 
                    fac201 * absa_13(ig,ind1+2) + 
                    fac011 * absa_13(ig,ind1+9) + 
                    fac111 * absa_13(ig,ind1+10) + 
                    fac211 * absa_13(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * 
                    (fac201 * absa_13(ig,ind1-1) + 
                    fac101 * absa_13(ig,ind1) + 
                    fac001 * absa_13(ig,ind1+1) + 
                    fac211 * absa_13(ig,ind1+8) + 
                    fac111 * absa_13(ig,ind1+9) + 
                    fac011 * absa_13(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * 
                    (fac001 * absa_13(ig,ind1) + 
                    fac101 * absa_13(ig,ind1+1) + 
                    fac011 * absa_13(ig,ind1+9) + 
                    fac111 * absa_13(ig,ind1+10));
            }
            taug(lw_rad.ngs12+ig,lay,m,n) = tau_major + tau_major1 
                 + tauself + taufor 
                 + adjcolco2*absco2 
                 + colco(lay)*absco;
            fracs(lw_rad.ngs12+ig,lay,m,n) = fracrefa_13(ig,jpl) + fpl * 
                 (fracrefa_13(ig,jpl+1)-fracrefa_13(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        indm = indminor(lay);
        for(int ig = 1; ig <= lw_rad.ng13;ig++){
            abso3 = kb_mo3_13(ig,indm) + minorfrac(lay) *   
                 (kb_mo3_13(ig,indm+1) - kb_mo3_13(ig,indm));
            taug(lw_rad.ngs12+ig,lay,m,n) = colo3(lay)*abso3;
            fracs(lw_rad.ngs12+ig,lay,m,n) =  fracrefb_13(ig);
        }
    }
}
KOKKOS_FUNCTION
void taugb14(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real1d_u &fracrefa_14,const real1d_u &fracrefb_14,const real2d_u &absa_14,const real3d_u &ka_14,const real2d_u &absb_14,const real3d_u &kb_14,const real2d_u &selfref_14,const real2d_u &forref_14,int m,int n){
    // band 14:  2250-2380 cm-1 (low - co2; high - co2)
    int lay, ind0, ind1, inds, indf, ig;
    real tauself, taufor;
    for(int lay = 1; lay <= laytrop; lay++){
        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[13] + 1;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[13] + 1;
        inds = indself(lay);
        indf = indfor(lay);
        for(int ig = 1; ig <= lw_rad.ng14; ig++){
            tauself = selffac(lay) * (selfref_14(ig,inds) + selffrac(lay) * 
                 (selfref_14(ig,inds+1) - selfref_14(ig,inds)));
            taufor =  forfac(lay) * (forref_14(ig,indf) + forfrac(lay) * 
                 (forref_14(ig,indf+1) - forref_14(ig,indf)));
            taug(lw_rad.ngs13+ig,lay,m,n) = colco2(lay) * 
                 (fac00(lay) * absa_14(ig,ind0) + 
                 fac10(lay) * absa_14(ig,ind0+1) + 
                 fac01(lay) * absa_14(ig,ind1) + 
                 fac11(lay) * absa_14(ig,ind1+1)) 
                 + tauself + taufor;
            fracs(lw_rad.ngs13+ig,lay,m,n) = fracrefa_14(ig);
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[13] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[13] + 1;
        for(int ig = 1; ig <= lw_rad.ng14; ig++){
            taug(lw_rad.ngs13+ig,lay,m,n) = colco2(lay) * 
                 (fac00(lay) * absb_14(ig,ind0) + 
                 fac10(lay) * absb_14(ig,ind0+1) + 
                 fac01(lay) * absb_14(ig,ind1) + 
                 fac11(lay) * absb_14(ig,ind1+1));
            fracs(lw_rad.ngs13+ig,lay,m,n) = fracrefb_14(ig);
        }
    }
}
KOKKOS_FUNCTION
void taugb15(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_15,const real2d_u &absa_15,const real4d_u &ka_15,const real3d_u &ka_mn2_15,const real2d_u &selfref_15,const real2d_u &forref_15,int m,int n){
    /*
        band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
                                (high - nothing)
    */
    int lay, ind0, ind1, inds, indf, indm, ig;
    int js, js1, jmn2, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_mn2, specparm_mn2, specmult_mn2, fmn2;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real scalen2, tauself, taufor, n2m1, n2m2, taun2;
    real refrat_planck_a, refrat_m_a;
    real tau_major, tau_major1;
    /*
        ! Minor gas mapping level : 
        !     Lower - Nitrogen Continuum, P = 1053., T = 294.

        ! Calculate reference ratio to be used in calculation of Planck
        ! fraction in lower atmosphere.
    */
    // P = 1053. mb (Level 1)
    refrat_planck_a = lw_rad.chi_mls[3][0]/lw_rad.chi_mls[1][0];
    // P = 1053. 
    refrat_m_a = lw_rad.chi_mls[3][0]/lw_rad.chi_mls[1][0];
    /*
        ! Compute the optical depth by interpolating in ln(pressure), 
        ! temperature, and appropriate species.  Below laytrop, the water
        ! vapor self-continuum and foreign continuum is interpolated 
        ! (in temperature) separately.
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = coln2o(lay) + rat_n2oco2(lay)*colco2(lay);
        specparm = coln2o(lay)/speccomb;
        if (specparm >= lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = coln2o(lay) + rat_n2oco2_1(lay)*colco2(lay);
        specparm1 = coln2o(lay)/speccomb1;
        if (specparm1 >= lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_mn2 = coln2o(lay) + refrat_m_a*colco2(lay);
        specparm_mn2 = coln2o(lay)/speccomb_mn2;
        if (specparm_mn2 >= lw_rad.oneminus) specparm_mn2 = lw_rad.oneminus;
        specmult_mn2 = 8.*specparm_mn2;
        jmn2 = 1 + int(specmult_mn2);
        fmn2 = fmod(specmult_mn2,1.0);

        speccomb_planck = coln2o(lay)+refrat_planck_a*colco2(lay);
        specparm_planck = coln2o(lay)/speccomb_planck;
        if (specparm_planck >= lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[14] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[14] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        indm = indminor(lay);
        
        scalen2 = colbrd(lay)*scaleminor(lay);
        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }
        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1 ;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng15; ig++){
            tauself = selffac(lay)* (selfref_15(ig,inds) + selffrac(lay) * 
                 (selfref_15(ig,inds+1) - selfref_15(ig,inds)));
            taufor =  forfac(lay) * (forref_15(ig,indf) + forfrac(lay) * 
                 (forref_15(ig,indf+1) - forref_15(ig,indf)));
            n2m1 = ka_mn2_15(ig,jmn2,indm) + fmn2 * 
                 (ka_mn2_15(ig,jmn2+1,indm) - ka_mn2_15(ig,jmn2,indm));
            n2m2 = ka_mn2_15(ig,jmn2,indm+1) + fmn2 * 
                 (ka_mn2_15(ig,jmn2+1,indm+1) - ka_mn2_15(ig,jmn2,indm+1));
            taun2 = scalen2 * (n2m1 + minorfrac(lay) * (n2m2 - n2m1));
            if (specparm < 0.125){
               tau_major = speccomb * 
                    (fac000 * absa_15(ig,ind0) + 
                    fac100 * absa_15(ig,ind0+1) + 
                    fac200 * absa_15(ig,ind0+2) + 
                    fac010 * absa_15(ig,ind0+9) + 
                    fac110 * absa_15(ig,ind0+10) + 
                    fac210 * absa_15(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * 
                    (fac200 * absa_15(ig,ind0-1) + 
                    fac100 * absa_15(ig,ind0) + 
                    fac000 * absa_15(ig,ind0+1) + 
                    fac210 * absa_15(ig,ind0+8) + 
                    fac110 * absa_15(ig,ind0+9) + 
                    fac010 * absa_15(ig,ind0+10));
            }else{
               tau_major = speccomb * 
                    (fac000 * absa_15(ig,ind0) + 
                    fac100 * absa_15(ig,ind0+1) +
                    fac010 * absa_15(ig,ind0+9) + 
                    fac110 * absa_15(ig,ind0+10));
            } 

            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * 
                    (fac001 * absa_15(ig,ind1) + 
                    fac101 * absa_15(ig,ind1+1) + 
                    fac201 * absa_15(ig,ind1+2) + 
                    fac011 * absa_15(ig,ind1+9) + 
                    fac111 * absa_15(ig,ind1+10) + 
                    fac211 * absa_15(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * 
                    (fac201 * absa_15(ig,ind1-1) + 
                    fac101 * absa_15(ig,ind1) + 
                    fac001 * absa_15(ig,ind1+1) + 
                    fac211 * absa_15(ig,ind1+8) + 
                    fac111 * absa_15(ig,ind1+9) + 
                    fac011 * absa_15(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * 
                    (fac001 * absa_15(ig,ind1) + 
                    fac101 * absa_15(ig,ind1+1) + 
                    fac011 * absa_15(ig,ind1+9) + 
                    fac111 * absa_15(ig,ind1+10));
            }
            taug(lw_rad.ngs14+ig,lay,m,n) = tau_major + tau_major1 
                 + tauself + taufor 
                 + taun2;
            fracs(lw_rad.ngs14+ig,lay,m,n) = fracrefa_15(ig,jpl) + fpl * 
                 (fracrefa_15(ig,jpl+1)-fracrefa_15(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        for(int ig = 1; ig <= lw_rad.ng15;ig++){
            taug(lw_rad.ngs14+ig,lay,m,n) = 0.0;
            fracs(lw_rad.ngs14+ig,lay,m,n) = 0.0;
        }
    }
}
KOKKOS_FUNCTION
void taugb16(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
    real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
    real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
    real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
    real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
    const real2d_u &fracrefa_16,const real1d_u &fracrefb_16,const real2d_u &absa_16,const real4d_u &ka_16,const real2d_u &absb_16,const real3d_u &kb_16,      
    const real2d_u &selfref_16,const real2d_u &forref_16,int m,int n){
    // band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
    int lay, ind0, ind1, inds, indf, ig;
    int js, js1, jpl;
    real speccomb, specparm, specmult, fs;
    real speccomb1, specparm1, specmult1, fs1;
    real speccomb_planck, specparm_planck, specmult_planck, fpl;
    real p, p4, fk0, fk1, fk2;
    real fac000, fac100, fac200, fac010, fac110, fac210;
    real fac001, fac101, fac201, fac011, fac111, fac211;
    real tauself, taufor;
    real refrat_planck_a;
    real tau_major, tau_major1;
    /*
        ! Calculate reference ratio to be used in calculation of Planck
        ! fraction in lower atmosphere.
    */
    // P = 387. mb (Level 6)
    refrat_planck_a = lw_rad.chi_mls[0][5]/lw_rad.chi_mls[5][5];
    /*
        ! Compute the optical depth by interpolating in ln(pressure), 
        ! temperature,and appropriate species.  Below laytrop, the water
        ! vapor self-continuum and foreign continuum is interpolated 
        ! (in temperature) separately. 
    */
    // Lower atmosphere loop
    for(int lay = 1; lay <= laytrop; lay++){
        speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay);
        specparm = colh2o(lay)/speccomb;
        if (specparm > lw_rad.oneminus) specparm = lw_rad.oneminus;
        specmult = 8.*(specparm);
        js = 1 + int(specmult);
        fs = fmod(specmult,1.0);

        speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay);
        specparm1 = colh2o(lay)/speccomb1;
        if (specparm1 > lw_rad.oneminus) specparm1 = lw_rad.oneminus;
        specmult1 = 8.*(specparm1);
        js1 = 1 + int(specmult1);
        fs1 = fmod(specmult1,1.0);

        speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay);
        specparm_planck = colh2o(lay)/speccomb_planck;
        if (specparm_planck > lw_rad.oneminus) specparm_planck=lw_rad.oneminus;
        specmult_planck = 8.*specparm_planck;
        jpl= 1 + int(specmult_planck);
        fpl = fmod(specmult_planck,1.0);

        ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*lw_rad.nspa[15] + js;
        ind1 = (jp(lay)*5+(jt1(lay)-1))*lw_rad.nspa[15] + js1;
        inds = indself(lay);
        indf = indfor(lay);
        if (specparm < 0.125){
            p = fs - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else if (specparm > 0.875){
            p = -fs; 
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac000 = fk0*fac00(lay);
            fac100 = fk1*fac00(lay);
            fac200 = fk2*fac00(lay);
            fac010 = fk0*fac10(lay);
            fac110 = fk1*fac10(lay);
            fac210 = fk2*fac10(lay);
        }else{
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
        }

        if (specparm1 < 0.125){
            p = fs1 - 1;
            p4 = std::pow(p,4);
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else if (specparm1 > 0.875){
            p = -fs1;
            p4 = std::pow(p,4);;
            fk0 = p4;
            fk1 = 1 - p - 2.0*p4;
            fk2 = p + p4;
            fac001 = fk0*fac01(lay);
            fac101 = fk1*fac01(lay);
            fac201 = fk2*fac01(lay);
            fac011 = fk0*fac11(lay);
            fac111 = fk1*fac11(lay);
            fac211 = fk2*fac11(lay);
        }else{
            fac001 = (1. - fs1) * fac01(lay);
            fac011 = (1. - fs1) * fac11(lay);
            fac101 = fs1 * fac01(lay);
            fac111 = fs1 * fac11(lay);
        }
        for(int ig = 1; ig <= lw_rad.ng16; ig++){
            tauself = selffac(lay)* (selfref_16(ig,inds) + selffrac(lay) * 
                 (selfref_16(ig,inds+1) - selfref_16(ig,inds)));
            taufor =  forfac(lay) * (forref_16(ig,indf) + forfrac(lay) * 
                 (forref_16(ig,indf+1) - forref_16(ig,indf)));
            if (specparm < 0.125){
               tau_major = speccomb * 
                    (fac000 * absa_16(ig,ind0) + 
                    fac100 * absa_16(ig,ind0+1) + 
                    fac200 * absa_16(ig,ind0+2) + 
                    fac010 * absa_16(ig,ind0+9) + 
                    fac110 * absa_16(ig,ind0+10) + 
                    fac210 * absa_16(ig,ind0+11));
            }else if (specparm > 0.875){
               tau_major = speccomb * 
                    (fac200 * absa_16(ig,ind0-1) + 
                    fac100 * absa_16(ig,ind0) +
                    fac000 * absa_16(ig,ind0+1) + 
                    fac210 * absa_16(ig,ind0+8) + 
                    fac110 * absa_16(ig,ind0+9) + 
                    fac010 * absa_16(ig,ind0+10));
            }else{
               tau_major = speccomb * 
                    (fac000 * absa_16(ig,ind0) + 
                    fac100 * absa_16(ig,ind0+1) + 
                    fac010 * absa_16(ig,ind0+9) + 
                    fac110 * absa_16(ig,ind0+10));
            }

            if (specparm1 < 0.125){
               tau_major1 = speccomb1 * 
                    (fac001 * absa_16(ig,ind1) + 
                    fac101 * absa_16(ig,ind1+1) + 
                    fac201 * absa_16(ig,ind1+2) + 
                    fac011 * absa_16(ig,ind1+9) + 
                    fac111 * absa_16(ig,ind1+10) + 
                    fac211 * absa_16(ig,ind1+11));
            }else if (specparm1 > 0.875){
               tau_major1 = speccomb1 * 
                    (fac201 * absa_16(ig,ind1-1) + 
                    fac101 * absa_16(ig,ind1) + 
                    fac001 * absa_16(ig,ind1+1) + 
                    fac211 * absa_16(ig,ind1+8) + 
                    fac111 * absa_16(ig,ind1+9) + 
                    fac011 * absa_16(ig,ind1+10));
            }else{
               tau_major1 = speccomb1 * 
                    (fac001 * absa_16(ig,ind1) + 
                    fac101 * absa_16(ig,ind1+1) + 
                    fac011 * absa_16(ig,ind1+9) + 
                    fac111 * absa_16(ig,ind1+10));
            }
            taug(lw_rad.ngs15+ig,lay,m,n) = tau_major + tau_major1 
                 + tauself + taufor;
            fracs(lw_rad.ngs15+ig,lay,m,n) = fracrefa_16(ig,jpl) + fpl * 
                 (fracrefa_16(ig,jpl+1)-fracrefa_16(ig,jpl));
        }
    }
    // Upper atmosphere loop
    for(int lay = laytrop+1; lay <= nlayers; lay++){
        ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*lw_rad.nspb[15] + 1;
        ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*lw_rad.nspb[15] + 1;
        for(int ig = 1; ig <= lw_rad.ng16; ig++){
            taug(lw_rad.ngs15+ig,lay,m,n) = colch4(lay) * 
                (fac00(lay) * absb_16(ig,ind0) + 
                fac10(lay) * absb_16(ig,ind0+1) + 
                fac01(lay) * absb_16(ig,ind1) + 
                fac11(lay) * absb_16(ig,ind1+1));
            fracs(lw_rad.ngs15+ig,lay,m,n) = fracrefb_16(ig);
        }
    }
}
KOKKOS_FUNCTION void taumol(int &nlayers,real1d_u &pavel,const real4d &wx,real1d_u &coldry,int &laytrop,int1d_u &jp,int1d_u &jt,int1d_u &jt1,const real4d &planklay,const real4d &planklev,real1d_u &plankbnd,
real1d_u &colh2o,real1d_u &colco2,real1d_u &colo3,real1d_u &coln2o,real1d_u &colco,real1d_u &colch4,real1d_u &colo2,real1d_u &colbrd,real1d_u &fac00,real1d_u &fac01,
real1d_u &fac10,real1d_u &fac11,real1d_u &rat_h2oco2,real1d_u &rat_h2oco2_1,real1d_u &rat_h2oo3,real1d_u &rat_h2oo3_1,real1d_u &rat_h2on2o,real1d_u &rat_h2on2o_1,real1d_u &rat_h2och4,
real1d_u &rat_h2och4_1,real1d_u &rat_n2oco2,real1d_u &rat_n2oco2_1,real1d_u &rat_o3co2,real1d_u &rat_o3co2_1,real1d_u &selffac,real1d_u &selffrac,int1d_u &indself,
real1d_u &forfac,real1d_u &forfrac,int1d_u &indfor,real1d_u &minorfrac,real1d_u &scaleminor,real1d_u &scaleminorn2,int1d_u &indminor,const real4d &fracs,const real4d &taug,
const real1d_u &fracrefa_01,const real1d_u &fracrefb_01,const real2d_u &absa_01,const real3d_u &ka_01,const real2d_u &absb_01,const real3d_u &kb_01,const real2d_u &ka_mn2_01,const real2d_u &kb_mn2_01,const real2d_u &selfref_01,const real2d_u &forref_01,                 
const real1d_u &fracrefa_02,const real1d_u &fracrefb_02,const real2d_u &absa_02,const real3d_u &ka_02,const real2d_u &absb_02,const real3d_u &kb_02,const real2d_u &selfref_02,const real2d_u &forref_02,                                     
const real2d_u &fracrefa_03,const real2d_u &fracrefb_03,const real2d_u &absa_03,const real4d_u &ka_03,const real2d_u &absb_03,const real4d_u &kb_03,const real3d_u &ka_mn2o_03,const real3d_u &kb_mn2o_03,const real2d_u &selfref_03,const real2d_u &forref_03,               
const real2d_u &fracrefa_04,const real2d_u &fracrefb_04,const real2d_u &absa_04,const real4d_u &ka_04,const real2d_u &absb_04,const real4d_u &kb_04,const real2d_u &selfref_04,const real2d_u &forref_04,                                     
const real2d_u &fracrefa_05,const real2d_u &fracrefb_05,const real2d_u &absa_05,const real4d_u &ka_05,const real2d_u &absb_05,const real4d_u &kb_05,const real3d_u &ka_mo3_05,const real2d_u &selfref_05,const real2d_u &forref_05,const real1d_u &ccl4_05,                   
const real1d_u &fracrefa_06,const real2d_u &absa_06,const real3d_u &ka_06,const real2d_u &ka_mco2_06,const real2d_u &selfref_06,const real2d_u &forref_06,const real1d_u &cfc11adj_06,const real1d_u &cfc12_06,                           
const real2d_u &fracrefa_07,const real1d_u &fracrefb_07,const real2d_u &absa_07,const real4d_u &ka_07,const real2d_u &absb_07,const real3d_u &kb_07,const real3d_u &ka_mco2_07,const real2d_u &kb_mco2_07,const real2d_u &selfref_07,const real2d_u &forref_07,               
const real1d_u &fracrefa_08,const real1d_u &fracrefb_08,const real2d_u &absa_08,const real3d_u &ka_08,const real2d_u &absb_08,const real3d_u &kb_08,const real2d_u &ka_mco2_08,const real2d_u &ka_mn2o_08,const real2d_u &ka_mo3_08,const real2d_u &kb_mco2_08,const real2d_u &kb_mn2o_08,    
const real2d_u &selfref_08,const real2d_u &forref_08,const real1d_u &cfc12_08,const real1d_u &cfc22adj_08,const real2d_u &fracrefa_09,const real1d_u &fracrefb_09,const real2d_u &absa_09,const real4d_u &ka_09,const real2d_u &absb_09,const real3d_u &kb_09,      
const real3d_u &ka_mn2o_09,const real2d_u &kb_mn2o_09,const real2d_u &selfref_09,const real2d_u &forref_09,const real1d_u &fracrefa_10,const real1d_u &fracrefb_10,const real2d_u &absa_10,const real3d_u &ka_10,const real2d_u &absb_10,const real3d_u &kb_10,      
const real2d_u &selfref_10,const real2d_u &forref_10,const real1d_u &fracrefa_11,const real1d_u &fracrefb_11,const real2d_u &absa_11,const real3d_u &ka_11,const real2d_u &absb_11,const real3d_u &kb_11,      
const real2d_u &ka_mo2_11,const real2d_u &kb_mo2_11,const real2d_u &selfref_11,const real2d_u &forref_11,const real2d_u &fracrefa_12,const real2d_u &absa_12,const real4d_u &ka_12,const real2d_u &selfref_12,const real2d_u &forref_12,           
const real2d_u &fracrefa_13,const real1d_u &fracrefb_13,const real2d_u &absa_13,const real4d_u &ka_13,const real3d_u &ka_mco2_13,const real3d_u &ka_mco_13,const real2d_u &kb_mo3_13,const real2d_u &selfref_13,const real2d_u &forref_13,                 
const real1d_u &fracrefa_14,const real1d_u &fracrefb_14,const real2d_u &absa_14,const real3d_u &ka_14,const real2d_u &absb_14,const real3d_u &kb_14,const real2d_u &selfref_14,const real2d_u &forref_14,                                     
const real2d_u &fracrefa_15,const real2d_u &absa_15,const real4d_u &ka_15,const real3d_u &ka_mn2_15,const real2d_u &selfref_15,const real2d_u &forref_15,const real2d_u &fracrefa_16,const real1d_u &fracrefb_16,const real2d_u &absa_16,const real4d_u &ka_16,const real2d_u &absb_16,const real3d_u &kb_16,      
const real2d_u &selfref_16,const real2d_u &forref_16,int m,int n){
/*
    !----------------------------------------------------------------------------

    ! *******************************************************************************
    ! *                                                                             *
    ! *                  Optical depths developed for the                           *
    ! *                                                                             *
    ! *                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
    ! *                                                                             *
    ! *                                                                             *
    ! *            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
    ! *                        131 HARTWELL AVENUE                                  *
    ! *                        LEXINGTON, MA 02421                                  *
    ! *                                                                             *
    ! *                                                                             *
    ! *                           ELI J. MLAWER                                     * 
    ! *                         JENNIFER DELAMERE                                   * 
    ! *                         STEVEN J. TAUBMAN                                   *
    ! *                         SHEPARD A. CLOUGH                                   *
    ! *                                                                             *
    ! *                                                                             *
    ! *                                                                             *
    ! *                                                                             *
    ! *                       email:  mlawer@aer.com                                *
    ! *                       email:  jdelamer@aer.com                              *
    ! *                                                                             *
    ! *        The authors wish to acknowledge the contributions of the             *
    ! *        following people:  Karen Cady-Pereira, Patrick D. Brown,             *  
    ! *        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
    ! *                                                                             *
    ! *******************************************************************************
    ! *                                                                             *
    ! *  Revision for g-point reduction: Michael J. Iacono, AER, Inc.               *
    ! *                                                                             *
    ! *******************************************************************************
    ! *     TAUMOL                                                                  *
    ! *                                                                             *
    ! *     This file contains the subroutines TAUGBn (where n goes from            *
    ! *     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
    ! *     per g-value and layer for band n.                                       *
    ! *                                                                             *
    ! *  Output:  optical depths (unitless)                                         *
    ! *           fractions needed to compute Planck functions at every layer       *
    ! *               and g-value                                                   *
    ! *                                                                             *
    ! *     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
    ! *     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
    ! *                                                                             *
    ! *  Input                                                                      *
    ! *                                                                             *
    ! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
    ! *     COMMON /PRECISE/  ONEMINUS                                              *
    ! *     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
    ! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
    ! *     COMMON /PROFDATA/ LAYTROP,                                              *
    ! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
    ! *    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
    ! *    &                  COLO2(MXLAY)
    ! *     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
    ! *    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
    ! *     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
    ! *     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
    ! *                                                                             *
    ! *     Description:                                                            *
    ! *     NG(IBAND) - number of g-values in band IBAND                            *
    ! *     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
    ! *                   atmospheres that are stored for band IBAND per            *
    ! *                   pressure level and temperature.  Each of these            *
    ! *                   atmospheres has different relative amounts of the         *
    ! *                   key species for the band (i.e. different binary           *
    ! *                   species parameters).                                      *
    ! *     NSPB(IBAND) - same for upper atmosphere                                 *
    ! *     ONEMINUS - since problems are caused in some cases by interpolation     *
    ! *                parameters equal to or greater than 1, for these cases       *
    ! *                these parameters are set to this value, slightly < 1.        *
    ! *     PAVEL - layer pressures (mb)                                            *
    ! *     TAVEL - layer temperatures (degrees K)                                  *
    ! *     PZ - level pressures (mb)                                               *
    ! *     TZ - level temperatures (degrees K)                                     *
    ! *     LAYTROP - layer at which switch is made from one combination of         *
    ! *               key species to another                                        *
    ! *     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
    ! *               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
    ! *               respectively (molecules/cm**2)                                *
    ! *     FACij(LAY) - for layer LAY, these are factors that are needed to        *
    ! *                  compute the interpolation factors that multiply the        *
    ! *                  appropriate reference k-values.  A value of 0 (1) for      *
    ! *                  i,j indicates that the corresponding factor multiplies     *
    ! *                  reference k-value for the lower (higher) of the two        *
    ! *                  appropriate temperatures, and altitudes, respectively.     *
    ! *     JP - the index of the lower (in altitude) of the two appropriate        *
    ! *          reference pressure levels needed for interpolation                 *
    ! *     JT, JT1 - the indices of the lower of the two appropriate reference     *
    ! *               temperatures needed for interpolation (for pressure           *
    ! *               levels JP and JP+1, respectively)                             *
    ! *     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
    ! *               (water vapor density)/(atmospheric density at 296K and        *
    ! *               1013 mb)                                                      *
    ! *     SELFFRAC - factor needed for temperature interpolation of reference     *
    ! *                water vapor self-continuum data                              *
    ! *     INDSELF - index of the lower of the two appropriate reference           *
    ! *               temperatures needed for the self-continuum interpolation      *
    ! *     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
    ! *     FORFRAC - factor needed for temperature interpolation of reference      *
    ! *                water vapor foreign-continuum data                           *
    ! *     INDFOR  - index of the lower of the two appropriate reference           *
    ! *               temperatures needed for the foreign-continuum interpolation   *
    ! *                                                                             *
    ! *  Data input                                                                 *
    ! *     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
    ! *                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
    ! *        (note:  n is the band number,'MGAS' is the species name of the minor *
    ! *         gas)                                                                *
    ! *                                                                             *
    ! *     Description:                                                            *
    ! *     KA - k-values for low reference atmospheres (key-species only)          *
    ! *          (units: cm**2/molecule)                                            *
    ! *     KB - k-values for high reference atmospheres (key-species only)         *
    ! *          (units: cm**2/molecule)                                            *
    ! *     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
    ! *          (units: cm**2/molecule)                                            *
    ! *     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
    ! *          (units: cm**2/molecule)                                            *
    ! *     SELFREF - k-values for water vapor self-continuum for reference         *
    ! *               atmospheres (used below LAYTROP)                              *
    ! *               (units: cm**2/molecule)                                       *
    ! *     FORREF  - k-values for water vapor foreign-continuum for reference      *
    ! *               atmospheres (used below/above LAYTROP)                        *
    ! *               (units: cm**2/molecule)                                       *
    ! *                                                                             *
    ! *     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
    ! *     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
    ! *                                                                             *
    !*******************************************************************************
*/
    taugb1(nlayers,pavel,wx,coldry,laytrop,jp,jt,jt1,planklay,planklev,plankbnd,
        colh2o,colco2,colo3,coln2o,colco,colch4,colo2,colbrd,fac00,fac01,
        fac10,fac11,rat_h2oco2,rat_h2oco2_1,rat_h2oo3,rat_h2oo3_1,rat_h2on2o,rat_h2on2o_1,rat_h2och4,
        rat_h2och4_1,rat_n2oco2,rat_n2oco2_1,rat_o3co2,rat_o3co2_1,selffac,selffrac,indself,
        forfac,forfrac,indfor,minorfrac,scaleminor,scaleminorn2,indminor,fracs,taug,
        fracrefa_01,fracrefb_01,absa_01,ka_01,absb_01,kb_01,ka_mn2_01,kb_mn2_01,selfref_01,forref_01,m,n);
    taugb2(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,             
        fracrefa_02, fracrefb_02, absa_02, ka_02, absb_02, kb_02, selfref_02, forref_02,m,n);
    taugb3(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_03, fracrefb_03, absa_03, ka_03, absb_03, kb_03, ka_mn2o_03, kb_mn2o_03, selfref_03, forref_03,m,n);
    taugb4(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_04, fracrefb_04, absa_04, ka_04, absb_04, kb_04, selfref_04, forref_04,m,n);
    taugb5(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_05, fracrefb_05, absa_05, ka_05, absb_05, kb_05, ka_mo3_05, selfref_05, forref_05, ccl4_05,m,n);
    taugb6(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_06, absa_06, ka_06, ka_mco2_06, selfref_06, forref_06, cfc11adj_06, cfc12_06,m,n);
    taugb7(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_07, fracrefb_07, absa_07, ka_07, absb_07, kb_07, ka_mco2_07, kb_mco2_07, selfref_07, forref_07,m,n);
    taugb8(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_08, fracrefb_08, absa_08, ka_08, absb_08, kb_08, ka_mco2_08, ka_mn2o_08, ka_mo3_08, kb_mco2_08, kb_mn2o_08,selfref_08, forref_08, cfc12_08, cfc22adj_08,m,n);
    taugb9(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_09, fracrefb_09, absa_09, ka_09, absb_09, kb_09,ka_mn2o_09, kb_mn2o_09, selfref_09, forref_09,m,n);
    taugb10(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_10, fracrefb_10, absa_10, ka_10, absb_10, kb_10,selfref_10, forref_10,m,n);
    taugb11(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_11, fracrefb_11, absa_11, ka_11, absb_11, kb_11,ka_mo2_11, kb_mo2_11, selfref_11, forref_11,m,n);
    taugb12(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_12, absa_12, ka_12, selfref_12, forref_12,m,n);
    taugb13(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_13, fracrefb_13, absa_13, ka_13, ka_mco2_13, ka_mco_13, kb_mo3_13, selfref_13, forref_13,m,n);
    taugb14(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_14, fracrefb_14, absa_14, ka_14, absb_14, kb_14, selfref_14, forref_14,m,n);
    taugb15(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_15, absa_15, ka_15, ka_mn2_15, selfref_15, forref_15,m,n);
    taugb16(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, 
        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, 
        rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug,
        fracrefa_16, fracrefb_16, absa_16,ka_16, absb_16,kb_16,selfref_16, forref_16,m,n);
}
KOKKOS_FUNCTION void rrtmg_lw(int ncol, int nlay,int &icld, const real4d &play,const real4d &plev,const real4d &tlay, const real4d &tlev, const real3d &tsfc, const real4d &h2ovmr, const real4d &o3vmr,
    const real4d &co2vmr, const real4d &ch4vmr, const real4d &n2ovmr, const real4d &o2vmr, const real4d &cfc11vmr, const real4d &cfc12vmr, const real4d &cfc22vmr, const real4d &ccl4vmr,
    const real4d &emis,int &inflglw,int &iceflglw,int &liqflglw, const real5d &cldfmcl, const real5d &taucmcl, const real5d &ciwpmcl, const real5d &clwpmcl,
    const real5d &cswpmcl, const real4d &reicmcl, const real4d &relqmcl, const real4d &resnmcl, const real5d &tauaer, const real4d &uflx, const real4d &dflx, const real4d &hr,
    const real4d &uflxc, const real4d &dflxc, const real4d &hrc,const real1d_u &tau_tbl,const real1d_u &tfn_tbl,const real1d_u &exp_tbl,const real1d_u &fracrefa_01,const real1d_u &fracrefb_01,const real2d_u &absa_01,const real3d_u &ka_01,const real2d_u &absb_01,const real3d_u &kb_01,const real2d_u &ka_mn2_01,const real2d_u &kb_mn2_01,const real2d_u &selfref_01,const real2d_u &forref_01,                 
    const real1d_u &fracrefa_02,const real1d_u &fracrefb_02,const real2d_u &absa_02,const real3d_u &ka_02,const real2d_u &absb_02,const real3d_u &kb_02,const real2d_u &selfref_02,const real2d_u &forref_02,                                     
    const real2d_u &fracrefa_03,const real2d_u &fracrefb_03,const real2d_u &absa_03,const real4d_u &ka_03,const real2d_u &absb_03,const real4d_u &kb_03,const real3d_u &ka_mn2o_03,const real3d_u &kb_mn2o_03,const real2d_u &selfref_03,const real2d_u &forref_03,               
    const real2d_u &fracrefa_04,const real2d_u &fracrefb_04,const real2d_u &absa_04,const real4d_u &ka_04,const real2d_u &absb_04,const real4d_u &kb_04,const real2d_u &selfref_04,const real2d_u &forref_04,                                     
    const real2d_u &fracrefa_05,const real2d_u &fracrefb_05,const real2d_u &absa_05,const real4d_u &ka_05,const real2d_u &absb_05,const real4d_u &kb_05,const real3d_u &ka_mo3_05,const real2d_u &selfref_05,const real2d_u &forref_05,const real1d_u &ccl4_05,                   
    const real1d_u &fracrefa_06,const real2d_u &absa_06,const real3d_u &ka_06,const real2d_u &ka_mco2_06,const real2d_u &selfref_06,const real2d_u &forref_06,const real1d_u &cfc11adj_06,const real1d_u &cfc12_06,                           
    const real2d_u &fracrefa_07,const real1d_u &fracrefb_07,const real2d_u &absa_07,const real4d_u &ka_07,const real2d_u &absb_07,const real3d_u &kb_07,const real3d_u &ka_mco2_07,const real2d_u &kb_mco2_07,const real2d_u &selfref_07,const real2d_u &forref_07,               
    const real1d_u &fracrefa_08,const real1d_u &fracrefb_08,const real2d_u &absa_08,const real3d_u &ka_08,const real2d_u &absb_08,const real3d_u &kb_08,const real2d_u &ka_mco2_08,const real2d_u &ka_mn2o_08,const real2d_u &ka_mo3_08,const real2d_u &kb_mco2_08,const real2d_u &kb_mn2o_08,    
    const real2d_u &selfref_08,const real2d_u &forref_08,const real1d_u &cfc12_08,const real1d_u &cfc22adj_08,const real2d_u &fracrefa_09,const real1d_u &fracrefb_09,const real2d_u &absa_09,const real4d_u &ka_09,const real2d_u &absb_09,const real3d_u &kb_09,      
    const real3d_u &ka_mn2o_09,const real2d_u &kb_mn2o_09,const real2d_u &selfref_09,const real2d_u &forref_09,const real1d_u &fracrefa_10,const real1d_u &fracrefb_10,const real2d_u &absa_10,const real3d_u &ka_10,const real2d_u &absb_10,const real3d_u &kb_10,      
    const real2d_u &selfref_10,const real2d_u &forref_10,const real1d_u &fracrefa_11,const real1d_u &fracrefb_11,const real2d_u &absa_11,const real3d_u &ka_11,const real2d_u &absb_11,const real3d_u &kb_11,      
    const real2d_u &ka_mo2_11,const real2d_u &kb_mo2_11,const real2d_u &selfref_11,const real2d_u &forref_11,const real2d_u &fracrefa_12,const real2d_u &absa_12,const real4d_u &ka_12,const real2d_u &selfref_12,const real2d_u &forref_12,           
    const real2d_u &fracrefa_13,const real1d_u &fracrefb_13,const real2d_u &absa_13,const real4d_u &ka_13,const real3d_u &ka_mco2_13,const real3d_u &ka_mco_13,const real2d_u &kb_mo3_13,const real2d_u &selfref_13,const real2d_u &forref_13,                 
    const real1d_u &fracrefa_14,const real1d_u &fracrefb_14,const real2d_u &absa_14,const real3d_u &ka_14,const real2d_u &absb_14,const real3d_u &kb_14,const real2d_u &selfref_14,const real2d_u &forref_14,                                     
    const real2d_u &fracrefa_15,const real2d_u &absa_15,const real4d_u &ka_15,const real3d_u &ka_mn2_15,const real2d_u &selfref_15,const real2d_u &forref_15,const real2d_u &fracrefa_16,const real1d_u &fracrefb_16,const real2d_u &absa_16,const real4d_u &ka_16,const real2d_u &absb_16,const real3d_u &kb_16,      
    const real2d_u &selfref_16,const real2d_u &forref_16,const real4d &fracs,const real4d &taug,const real4d &taut,const real4d &cldfmc,const real4d &ciwpmc,const real4d &clwpmc,const real4d &cswpmc,const real4d &taucmc,
    const real4d &wkl,const real4d &wx,const real4d &taua,const real4d &planklay,const real4d &planklev,const real4d &efclfrac,const real4d &odcld,int m, int n){
    // const real3d &bbugas,const real3d &bbutot,const real3d &secdiff,const int3d &icldlyr,const real3d &uclfl,const real3d &dclfl,const real3d &clrurad,const real3d &clrdrad,
    // const real3d &atot,const real3d &atrans,int m, int n){
/*
    ! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
    ! model for application to GCMs, that has been adapted from RRTM_LW for
    ! improved efficiency.
    !
    ! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
    !  area, since this has to be called only once. 
    !
    ! This routine:
    !    a) calls INATM to read in the atmospheric profile from GCM;
    !       all layering in RRTMG is ordered from surface to toa. 
    !    b) calls CLDPRMC to set cloud optical depth for McICA based 
    !       on input cloud properties 
    !    c) calls SETCOEF to calculate various quantities needed for 
    !       the radiative transfer algorithm
    !    d) calls TAUMOL to calculate gaseous optical depths for each 
    !       of the 16 spectral bands
    !    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
    !       radiative transfer calculation using McICA, the Monte-Carlo 
    !       Independent Column Approximation, to represent sub-grid scale 
    !       cloud variability
    !    f) passes the necessary fluxes and cooling rates back to GCM
    !
    ! Two modes of operation are possible:
    !     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
    !     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
    !
    !    1) Standard, single forward model calculation (imca = 0)
    !    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
    !       JC, 2003) method is applied to the forward model calculation (imca = 1)
    !
    ! This call to RRTMG_LW must be preceeded by a call to the module
    !     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
    !     which will provide the cloud physical or cloud optical properties
    !     on the RRTMG quadrature point (ngpt) dimension.
    !     Two random number generators are available for use when imca = 1.
    !     This is chosen by setting flag irnd on input to mcica_subcol_gen_lw.
    !     1) KISSVEC (irnd = 0)
    !     2) Mersenne-Twister (irnd = 1)
    !
    ! Two methods of cloud property input are possible:
    !     Cloud properties can be input in one of two ways (controlled by input 
    !     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
    !     and subroutine rrtmg_lw_cldprop.f90 for further details):
    !
    !    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
    !    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);  
    !       cloud optical properties are calculated by cldprop or cldprmc based
    !       on input settings of iceflglw and liqflglw.  Ice particle size provided
    !       must be appropriately defined for the ice parameterization selected. 
    !
    ! One method of aerosol property input is possible:
    !     Aerosol properties can be input in only one way (controlled by input 
    !     flag iaer; see text file rrtmg_lw_instructions for further details):
    !
    !    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
    !       band average optical depth at the mid-point of each spectral band.
    !       RRTMG_LW currently treats only aerosol absorption;
    !       scattering capability is not presently available.
    !
    !
    ! ------- Modifications -------
    !
    ! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
    ! set of g-points for application to GCMs.  
    !
    !-- Original version (derived from RRTM_LW), reduction of g-points, other
    !   revisions for use with GCMs.  
    !     1999: M. J. Iacono, AER, Inc.
    !-- Adapted for use with NCAR/CAM.
    !     May 2004: M. J. Iacono, AER, Inc.
    !-- Revised to add McICA capability. 
    !     Nov 2005: M. J. Iacono, AER, Inc.
    !-- Conversion to F90 formatting for consistency with rrtmg_sw.
    !     Feb 2007: M. J. Iacono, AER, Inc.
    !-- Modifications to formatting to use assumed-shape arrays.
    !     Aug 2007: M. J. Iacono, AER, Inc.
    !-- Modified to add longwave aerosol absorption.
    !     Apr 2008: M. J. Iacono, AER, Inc.
*/
    int nlayers;             // total number of layers
    int istart;              // beginning band of calculation
    int iend;                // ending band of calculation
    int iout;                // output option flag (inactive)
    int iaer;                // aerosol option flag
    int iplon;               // column loop index
    int imca;                // flag for mcica [0=off, 1=on]
    int ims;                 // value for changing mcica permute seed
    int k;                   // layer loop index
    int ig;                  // g-point loop index
    real pwvcm;
    real tbound;                                                        // surface temperature (K)
// Atmosphere - setcoef
    int laytrop;                                                      // tropopause layer index
    
//Atmosphere/clouds - cldprop
    int ncbands;                                                        // number of cloud spectral bands
    int inflag;                                                         // flag for cloud property method
    int iceflag;                                                        // flag for ice cloud properties
    int liqflag;                                                        // flag for liquid cloud properties
    const int size = 110;
    // Atmosphere
    real pavel_s[size + 2];
    real1d_u pavel((real*)pavel_s,{1},{size+2});                          // layer pressures (mb) 
    real tavel_s[size + 2];
    real1d_u tavel((real*)tavel_s,{1},{size+2});                          // layer temperatures (K)
    real pz_s[size + 2];
    real1d_u pz((real*)pz_s,{0},{size+2});                                // level (interface) pressures (hPa, mb)
    real tz_s[size + 2];
    real1d_u tz((real*)tz_s,{0},{size+2});                                // level (interface) temperatures (K)
    real coldry_s[size + 2];
    real1d_u coldry((real*)coldry_s,{1},{size+2});                        // dry air column density (mol/cm2)
    real semiss_s[nbndlw + 1];
    real1d_u semiss((real*)semiss_s,{1},{nbndlw+1});                        // lw surface emissivity
    real wbrodl_s[size + 2];
    real1d_u wbrodl((real*)wbrodl_s,{1},{size+2});                        // broadening gas column density (mol/cm2)
    
// Atmosphere - setcoef
    int jp_s[size + 2];
    int jt_s[size + 2];
    int jt1_s[size + 2];
    int indself_s[size + 2];
    int indfor_s[size + 2];
    int indminor_s[size + 2];
    int1d_u jp((int*)jp_s,{1},{size+2});                               // lookup table index 
    int1d_u jt((int*)jt_s,{1},{size+2});                               // lookup table index 
    int1d_u jt1((int*)jt1_s,{1},{size+2});                             // lookup table index 
    int1d_u indself((int*)indself_s,{1},{size+2});        
    int1d_u indfor((int*)indfor_s,{1},{size+2});           
    int1d_u indminor((int*)indminor_s,{1},{size+2});   
    real plankbnd_s[nbndlw + 1];
    real1d_u plankbnd((real*)plankbnd_s,{1},{nbndlw+1});
    real colh2o_s[size + 2];
    real1d_u colh2o((real*)colh2o_s,{1},{size+2});                     // column amount (h2o)
    real colco2_s[size + 2];
    real1d_u colco2((real*)colco2_s,{1},{size+2});                    // column amount (co2)
    real colo3_s[size + 2];
    real1d_u colo3((real*)colo3_s,{1},{size+2});                      // column amount (o3)
    real coln2o_s[size + 2];
    real1d_u coln2o((real*)coln2o_s,{1},{size+2});                    // column amount (n2o)
    real colco_s[size + 2];
    real1d_u colco((real*)colco_s,{1},{size+2});                      // column amount (co)
    real colch4_s[size + 2];
    real1d_u colch4((real*)colch4_s,{1},{size+2});                    // column amount (ch4)
    real colo2_s[size + 2];
    real1d_u colo2((real*)colo2_s,{1},{size+2});                      // column amount (o2)
    real colbrd_s[size + 2];
    real1d_u colbrd((real*)colbrd_s,{1},{size+2});                    // column amount (broadening gases)
    real selffac_s[size + 2];
    real1d_u selffac((real*)selffac_s,{1},{size+2});
    real selffrac_s[size + 2];
    real1d_u selffrac((real*)selffrac_s,{1},{size+2});
    real forfac_s[size + 2];
    real1d_u forfac((real*)forfac_s,{1},{size+2});
    real forfrac_s[size + 2];
    real1d_u forfrac((real*)forfrac_s,{1},{size+2});
    real minorfrac_s[size + 2];
    real1d_u minorfrac((real*)minorfrac_s,{1},{size+2});
    real scaleminor_s[size + 2];
    real1d_u scaleminor((real*)scaleminor_s,{1},{size+2});
    real scaleminorn2_s[size + 2];
    real1d_u scaleminorn2((real*)scaleminorn2_s,{1},{size+2});
    real fac00_s[size + 2];
    real1d_u fac00((real*)fac00_s,{1},{size+2});
    real fac01_s[size + 2];
    real1d_u fac01((real*)fac01_s,{1},{size+2});
    real fac10_s[size + 2];
    real1d_u fac10((real*)fac10_s,{1},{size+2});
    real fac11_s[size + 2];
    real1d_u fac11((real*)fac11_s,{1},{size+2});
    real rat_h2oco2_s[size + 2];
    real1d_u rat_h2oco2((real*)rat_h2oco2_s,{1},{size+2});
    real rat_h2oco2_1_s[size + 2];
    real1d_u rat_h2oco2_1((real*)rat_h2oco2_1_s,{1},{size+2});
    real rat_h2oo3_s[size + 2];
    real1d_u rat_h2oo3((real*)rat_h2oo3_s,{1},{size+2});
    real rat_h2on2o_1_s[size + 2];
    real1d_u rat_h2on2o_1((real*)rat_h2on2o_1_s,{1},{size+2});
    real rat_h2oo3_1_s[size + 2];
    real1d_u rat_h2oo3_1((real*)rat_h2oo3_1_s,{1},{size+2});
    real rat_h2on2o_s[size + 2];
    real1d_u rat_h2on2o((real*)rat_h2on2o_s,{1},{size+2});
    real rat_h2och4_s[size + 2];
    real1d_u rat_h2och4((real*)rat_h2och4_s,{1},{size+2});
    real rat_h2och4_1_s[size + 2];
    real1d_u rat_h2och4_1((real*)rat_h2och4_1_s,{1},{size+2});
    real rat_n2oco2_s[size + 2];
    real1d_u rat_n2oco2((real*)rat_n2oco2_s,{1},{size+2});
    real rat_n2oco2_1_s[size + 2];
    real1d_u rat_n2oco2_1((real*)rat_n2oco2_1_s,{1},{size+2});
    real rat_o3co2_s[size + 2];
    real1d_u rat_o3co2((real*)rat_o3co2_s,{1},{size+2});
    real rat_o3co2_1_s[size + 2];
    real1d_u rat_o3co2_1((real*)rat_o3co2_1_s,{1},{size+2});
    
// Atmosphere/clouds - cldprmc [mcica]
    real relqmc_s[size + 2];
    real1d_u relqmc((real*)relqmc_s,{1},{size+2});                        // liquid particle effective radius (microns)
    real reicmc_s[size + 2];
    real1d_u reicmc((real*)reicmc_s,{1},{size+2});                        // ice particle effective size (microns)
    real resnmc_s[size + 2];
    real1d_u resnmc((real*)resnmc_s,{1},{size+2});                        // snow particle effective size (microns)
    
// Output
    real totuflux_s[size + 2];
    real1d_u totuflux((real*)totuflux_s,{0},{size+2});                    // upward longwave flux (w/m2)
    real totdflux_s[size + 2];
    real1d_u totdflux((real*)totdflux_s,{0},{size+2});                    // downward longwave flux (w/m2)
    real fnet_s[size + 2];
    real1d_u fnet((real*)fnet_s,{0},{size+2});                            // net longwave flux (w/m2)
    real htr_s[size + 2];
    real1d_u htr((real*)htr_s,{0},{size+2});                              // longwave heating rate (k/day)
    real totuclfl_s[size + 2];
    real1d_u totuclfl((real*)totuclfl_s,{0},{size+2});                    // clear sky upward longwave flux (w/m2)
    real totdclfl_s[size + 2];
    real1d_u totdclfl((real*)totdclfl_s,{0},{size+2});                    // clear sky downward longwave flux (w/m2)
    real fnetc_s[size + 2];
    real1d_u fnetc((real*)fnetc_s,{0},{size+2});                          // clear sky net longwave flux (w/m2)
    real htrc_s[size + 2];
    real1d_u htrc((real*)htrc_s,{0},{size+2});                            //  clear sky longwave heating rate (k/day)
    // Initializations
    // jm not thread safe      oneminus = 1. - 1.e-6_rb
    // jm not thread safe      pi = 2._rb * asin(1._rb)
    // jm not thread safe      fluxfac = pi * 2.e4_rb                  ! orig:   fluxfac = pi * 2.d4  
    istart = 1;
    iend = 16;
    iout = 0;
    ims = 1;
    /*
        ! Set imca to select calculation type:
        !  imca = 0, use standard forward model calculation
        !  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

        ! *** This version uses McICA (imca = 1) ***

        ! Set icld to select of clear or cloud calculation and cloud overlap method  
        ! icld = 0, clear only
        ! icld = 1, with clouds using random cloud overlap
        ! icld = 2, with clouds using maximum/random cloud overlap
        ! icld = 3, with clouds using maximum cloud overlap (McICA only)
    */
    if (icld < 0 || icld > 3) icld = 2;
    /*
        ! Set iaer to select aerosol option
        ! iaer = 0, no aerosols
        ! icld = 10, input total aerosol optical depth (tauaer) directly
    */
    iaer = 10;
    /*
        ! Call model and data initialization, compute lookup tables, perform
        ! reduction of g-points from 256 to 140 for input absorption coefficient 
        ! data and other arrays.
        !
        ! In a GCM this call should be placed in the model initialization
        ! area, since this has to be called only once.  
        !      call rrtmg_lw_ini(cpdair)

        !  This is the main longitude/column loop within RRTMG.
    */
    for(int iplon = 1; iplon <= ncol; iplon++){
        //  Prepare atmospheric profile from GCM for use in RRTMG, and define
        //  other input parameters. 
        inatm (iplon, nlay, icld, iaer, 
                play, plev, tlay, tlev, tsfc, h2ovmr, 
                o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, 
                cfc22vmr, ccl4vmr, emis, inflglw, iceflglw, liqflglw, 
                cldfmcl, taucmcl, ciwpmcl, clwpmcl, cswpmcl, reicmcl, relqmcl, resnmcl, tauaer, 
                nlayers, pavel, pz, tavel, tz, tbound, semiss, coldry, 
                wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, 
                cldfmc, taucmc, ciwpmc, clwpmc, cswpmc, reicmc, relqmc, resnmc, taua,m,n);
        // /*
        //     !  For cloudy atmosphere, use cldprop to set cloud optical properties based on
        //     !  input cloud physical properties.  Select method based on choices described
        //     !  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
        //     !  effective radius must be passed into cldprop.  Cloud fraction and cloud
        //     !  optical depth are transferred to rrtmg_lw arrays in cldprop. 
        // */
        cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, ciwpmc, 
                clwpmc, cswpmc, reicmc, relqmc, resnmc, ncbands, taucmc,m,n);
        // /*
        //     ! Calculate information needed by the radiative transfer routine
        //     ! that is specific to this atmosphere, especially some of the 
        //     ! coefficients and indices needed to compute the optical depths
        //     ! by interpolating data from stored reference atmospheres. 
        // */
        setcoef(nlayers, istart, pavel, tavel, tz, tbound, semiss, 
                coldry, wkl, wbrodl, 
                laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
                colh2o, colco2, colo3, coln2o, colco, colch4, colo2, 
                colbrd, fac00, fac01, fac10, fac11, 
                rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, 
                rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, 
                rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, 
                selffac, selffrac, indself, forfac, forfrac, indfor, 
                minorfrac, scaleminor, scaleminorn2, indminor,m,n);
        //  Calculate the gaseous optical depths and Planck fractions for 
        //  each longwave spectral band.
        taumol(nlayers, pavel, wx, coldry, 
                laytrop, jp, jt, jt1, planklay, planklev, plankbnd, 
                colh2o, colco2, colo3, coln2o, colco, colch4, colo2, 
                colbrd, fac00, fac01, fac10, fac11, 
                rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, 
                rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, 
                rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, 
                selffac, selffrac, indself, forfac, forfrac, indfor, 
                minorfrac, scaleminor, scaleminorn2, indminor, 
                fracs, taug,fracrefa_01, fracrefb_01, absa_01, ka_01, absb_01, kb_01, ka_mn2_01, kb_mn2_01, selfref_01, forref_01,                 
                fracrefa_02, fracrefb_02, absa_02, ka_02, absb_02, kb_02, selfref_02, forref_02,                                     
                fracrefa_03, fracrefb_03, absa_03, ka_03, absb_03, kb_03, ka_mn2o_03, kb_mn2o_03, selfref_03, forref_03,               
                fracrefa_04, fracrefb_04, absa_04, ka_04, absb_04, kb_04, selfref_04, forref_04,                                     
                fracrefa_05, fracrefb_05, absa_05, ka_05, absb_05, kb_05, ka_mo3_05, selfref_05, forref_05, ccl4_05,                   
                fracrefa_06, absa_06, ka_06, ka_mco2_06, selfref_06, forref_06, cfc11adj_06, cfc12_06,                           
                fracrefa_07, fracrefb_07, absa_07, ka_07, absb_07, kb_07, ka_mco2_07, kb_mco2_07, selfref_07, forref_07,               
                fracrefa_08, fracrefb_08, absa_08, ka_08, absb_08, kb_08, ka_mco2_08, ka_mn2o_08, ka_mo3_08, kb_mco2_08, kb_mn2o_08,    
                selfref_08, forref_08, cfc12_08, cfc22adj_08, fracrefa_09, fracrefb_09, absa_09, ka_09, absb_09, kb_09,      
                ka_mn2o_09, kb_mn2o_09, selfref_09, forref_09, fracrefa_10, fracrefb_10, absa_10, ka_10, absb_10, kb_10,      
                selfref_10, forref_10, fracrefa_11, fracrefb_11, absa_11, ka_11, absb_11, kb_11,      
                ka_mo2_11, kb_mo2_11, selfref_11, forref_11, fracrefa_12, absa_12, ka_12, selfref_12, forref_12,           
                fracrefa_13, fracrefb_13, absa_13, ka_13, ka_mco2_13, ka_mco_13, kb_mo3_13, selfref_13, forref_13,                 
                fracrefa_14, fracrefb_14, absa_14, ka_14, absb_14, kb_14, selfref_14, forref_14,                                     
                fracrefa_15, absa_15, ka_15, ka_mn2_15, selfref_15, forref_15, fracrefa_16, fracrefb_16, absa_16,ka_16, absb_16,kb_16,      
                selfref_16, forref_16,m,n);
        
           
        // Combine gaseous and aerosol optical depths, if aerosol active
            if(iaer == 0){
                for(int k = 1; k <= nlayers; k++){
                    for(int ig = 1; ig <= ngptlw; ig++){
                        taut(ig,k,m,n) = taug(ig,k,m,n);
                    }
                }
            }else if(iaer == 10){
                auto &ngb = lw_rad.ngb;
                for(int k = 1; k <= nlayers; k++){
                    for(int ig = 1; ig <= ngptlw; ig++){
                        taut(ig,k,m,n) = taug(ig,k,m,n) + taua(ngb[ig - 1],k,m,n);
                    }
                }
            }
            /*
                ! Call the radiative transfer routine.
                ! Either routine can be called to do clear sky calculation.  If clouds
                ! are present, then select routine based on cloud overlap assumption
                ! to be used.  Clear sky calculation is done simultaneously.
                ! For McICA, RTRNMC is called for clear and cloudy calculations.
            */
            rtrnmc(nlayers, istart, iend, iout, pz, semiss, ncbands, 
                    cldfmc, taucmc, planklay, planklev, plankbnd, 
                    pwvcm, fracs, taut, 
                    totuflux, totdflux, fnet, htr, 
                    totuclfl, totdclfl, fnetc, htrc,tau_tbl,tfn_tbl,exp_tbl,efclfrac,odcld,m,n);
            for(int k = 0; k <= nlayers; k++){
                uflx(iplon,k+1,m,n) = totuflux(k);
                dflx(iplon,k+1,m,n) = totdflux(k);
                uflxc(iplon,k+1,m,n) = totuclfl(k);
                dflxc(iplon,k+1,m,n) = totdclfl(k);
            }
            for(int k = 0; k <= nlayers - 1; k++){
                hr(iplon,k+1,m,n) = htr(k);
                hrc(iplon,k+1,m,n) = htrc(k);
            }
    }
}
