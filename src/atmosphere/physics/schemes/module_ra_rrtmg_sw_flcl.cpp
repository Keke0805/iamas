#include <float.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>
#include <flcl-cxx.hpp>
#include <cstring>
#include <iostream>
#include "rrtmg_sw.hpp"
typedef double real;
using real1d_u = Kokkos::Experimental::OffsetView<real*, Kokkos::LayoutLeft, flcl::HostMemorySpace, Kokkos::MemoryUnmanaged>;
using real2d_u = Kokkos::Experimental::OffsetView<real**, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using real3d_u = Kokkos::Experimental::OffsetView<real***, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using real4d_u = Kokkos::Experimental::OffsetView<real****, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
using real5d_u = Kokkos::Experimental::OffsetView<real*****, Kokkos::LayoutLeft, flcl::HostMemorySpace,Kokkos::MemoryUnmanaged>;
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

using std::string;

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
struct rrtmg_sw rrtmgsw;

static sw_workspace *sw_ws = nullptr;
static bool is_sw_ws_inited = false;
static int sw_ws_kts = 0;
static int sw_ws_kte = 0;
static int sw_ws_its = 0;
static int sw_ws_ite = 0;
static int sw_ws_jts = 0;
static int sw_ws_jte = 0;
static int sw_ws_nbndsw = 0;
static int sw_ws_ngptsw = 0;
static int sw_ws_naerec = 0;
static int sw_ws_mxmol = 0;
static int sw_ws_noznlevels = 0;
static int sw_ws_size_nlay = 0;
const int RAD_SW_CHUNK_SIZE = 64;
const int ngptsw = rrtmgsw.ngptsw;
const int mg     = rrtmgsw.mg;
const int nbndsw = rrtmgsw.nbndsw;
const int maxinpx= rrtmgsw.maxinpx;
const int naerec = rrtmgsw.naerec;
const int mxmol = rrtmgsw.mxmol;
const int jpband = rrtmgsw.jpband;
const int jpb1  = rrtmgsw.jpb1;
const int jpb2  = rrtmgsw.jpb2;
const real heatfac = rrtmgsw.heatfac;
const real pi = rrtmgsw.pi;
const real rrsw_scon = rrtmgsw.rrsw_scon;
const real grav = rrtmgsw.grav;
const int nmol = rrtmgsw.nmol;
const real avogad = rrtmgsw.avogad;
const real od_lo = rrtmgsw.od_lo;
const real bpade = rrtmgsw.bpade;
const real tblint = rrtmgsw.tblint;
auto &rsrtaua = rrtmgsw.rsrtaua;
auto &rsrpiza = rrtmgsw.rsrpiza;
auto &rsrasya = rrtmgsw.rsrasya;
const int SIZE_NLAY = 110;
const int SIZE_NCOL = 2;
const real g = 9.80616;   //gravity
static unsigned long state[N]; /* the array for the state vector  */
static int left = 1;
static int initf = 0;
static unsigned long *next;

static bool sw_workspace_matches(int kts, int kte, int its, int ite, int jts, int jte,
                                 int nbndsw_in, int ngptsw_in, int naerec_in,
                                 int mxmol_in, int noznlevels_in, int size_nlay_in) {
    return is_sw_ws_inited && sw_ws != nullptr &&
           sw_ws_kts == kts &&
           sw_ws_kte == kte &&
           sw_ws_its == its &&
           sw_ws_ite == ite &&
           sw_ws_jts == jts &&
           sw_ws_jte == jte &&
           sw_ws_nbndsw == nbndsw_in &&
           sw_ws_ngptsw == ngptsw_in &&
           sw_ws_naerec == naerec_in &&
           sw_ws_mxmol == mxmol_in &&
           sw_ws_noznlevels == noznlevels_in &&
           sw_ws_size_nlay == size_nlay_in;
}

static sw_workspace& get_sw_ws(int kts, int kte, int its, int ite, int jts, int jte,
                               int nbndsw_in, int ngptsw_in, int naerec_in,
                               int mxmol_in, int noznlevels_in, int size_nlay_in) {
    if (!sw_workspace_matches(kts, kte, its, ite, jts, jte,
                              nbndsw_in, ngptsw_in, naerec_in,
                              mxmol_in, noznlevels_in, size_nlay_in)) {
        if (sw_ws != nullptr) {
            Kokkos::fence();
            delete sw_ws;
            sw_ws = nullptr;
        }

        sw_ws = new sw_workspace();
        sw_ws->init(kts, kte, its, ite, jts, jte,
                    nbndsw_in, ngptsw_in, naerec_in,
                    mxmol_in, noznlevels_in, size_nlay_in);

        sw_ws_kts = kts;
        sw_ws_kte = kte;
        sw_ws_its = its;
        sw_ws_ite = ite;
        sw_ws_jts = jts;
        sw_ws_jte = jte;
        sw_ws_nbndsw = nbndsw_in;
        sw_ws_ngptsw = ngptsw_in;
        sw_ws_naerec = naerec_in;
        sw_ws_mxmol = mxmol_in;
        sw_ws_noznlevels = noznlevels_in;
        sw_ws_size_nlay = size_nlay_in;
        is_sw_ws_inited = true;
    }

    return *sw_ws;
}

extern "C" void rrtmg_sw_workspace_finalize() {
    if (sw_ws != nullptr) {
        if (Kokkos::is_initialized()) {
            Kokkos::fence();
        }
        delete sw_ws;
        sw_ws = nullptr;
    }

    is_sw_ws_inited = false;
    sw_ws_kts = 0;
    sw_ws_kte = 0;
    sw_ws_its = 0;
    sw_ws_ite = 0;
    sw_ws_jts = 0;
    sw_ws_jte = 0;
    sw_ws_nbndsw = 0;
    sw_ws_ngptsw = 0;
    sw_ws_naerec = 0;
    sw_ws_mxmol = 0;
    sw_ws_noznlevels = 0;
    sw_ws_size_nlay = 0;
}

/* initializes state[N] with a seed */
void init_genrand_sw(unsigned long s)
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
void init_by_array_sw(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    init_genrand_sw(19650218UL);
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

void next_state_sw(void)
{
    unsigned long *p=state;
    int j;

    /* if init_genrand() has not been called, */
    /* a default initial seed is used         */
    if (initf==0) init_genrand_sw(5489UL);

    left = N;
    next = state;
    
    for (j=N-M+1; --j; p++) 
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (j=M; --j; p++) 
        *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32_sw(void)
{
    unsigned long y;

    if (--left == 0) next_state_sw();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31_sw(void)
{
    unsigned long y;

    if (--left == 0) next_state_sw();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (long)(y>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1_sw(void)
{
    unsigned long y;

    if (--left == 0) next_state_sw();
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
double genrand_real2_sw(void)
{
    unsigned long y;

    if (--left == 0) next_state_sw();
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
double genrand_real3_sw(void)
{
    unsigned long y;

    if (--left == 0) next_state_sw();
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
double genrand_res53_sw(void) 
{ 
    unsigned long a=genrand_int32_sw()>>5, b=genrand_int32_sw()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
}

string trim_sw(char *config_chem_opt,string &schemes){
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

KOKKOS_FUNCTION void inirad_sw(const real3d &o3prof,const real4d &plev,int kts,int kte,int m,int j);

KOKKOS_FUNCTION void vinterp_ozn_sw(int lchnk,int ncol,int pcols,int pver,const real4d& pmid_in,const real1d_u& pin_in,int levsiz,const real3d& ozmix_in,const real3d& o3vmr_out,int m,int j);

KOKKOS_FUNCTION void relcalc_sw(int ncol,int pcols,int pver,const real4d &t, const real3d &landfrac, const real3d &landm, const real3d &icefrac,const real4d &rel, const real3d &snowh,int m,int j);

KOKKOS_FUNCTION void reicalc_sw(int ncol,int pcols,int pver,const real4d &t, const real4d &re,const real1d_u &retab,int m,int j);

KOKKOS_FUNCTION void mcica_subcol_sw_kokkos(int iplon,int ncol,int nlay,int icld,int permuteseed,int irng,const real4d& play,
    const real4d& cldfrac,const real4d& ciwp,const real4d& clwp,const real4d& cswp,const real4d& rei,const real4d& rel,const real4d& res,const real5d& tauc,const real5d& ssac,const real5d& asmc,const real5d& fsfc,
    const real5d& cldfmcl,const real5d& ciwpmcl,const real5d& clwpmcl,const real5d& cswpmcl,const real4d& reicmcl,const real4d& relqmcl,const real4d& resnmcl,
    const real5d& taucmcl,const real5d& ssacmcl,const real5d& asmcmcl,const real5d& fsfcmcl,int m,int j);

KOKKOS_FUNCTION void rrtmg_sw (int ncol, int nlay, int icld,
    const real4d& play, const real4d& plev, const real4d& tlay, const real4d& tlev, const real3d& tsfc, 
    const real4d& h2ovmr, const real4d& o3vmr, const real4d& co2vmr, const real4d& ch4vmr, const real4d&n2ovmr, const real4d& o2vmr,
    const real3d& asdir,const real3d& asdif,const real3d& aldir,const real3d& aldif,
    const real3d& coszen,real& adjes,int& dyofyr,real& scon, 
    int& inflgsw ,int& iceflgsw,int& liqflgsw,const real5d& cldfmcl ,
    const real5d& taucmcl ,const real5d& ssacmcl ,const real5d& asmcmcl ,const real5d& fsfcmcl ,
    const real5d& ciwpmcl ,const real5d& clwpmcl ,const real5d& cswpmcl ,const real4d& reicmcl ,const real4d& relqmcl ,const real4d& resnmcl,
    const real5d& tauaer  ,const real5d& ssaaer  ,const real5d& asmaer  ,const real5d& ecaer   ,
    const real4d& swuflx  ,const real4d& swdflx  ,const real4d& swhr    ,const real4d& swuflxc ,const real4d& swdflxc ,const real4d& swhrc,int& aer_opt,
    const real4d& sibvisdir,const real4d& sibvisdif,const real4d& sibnirdir,const real4d& sibnirdif,
    const real4d& swdkdir,const real4d& swdkdif,
    const real2d_u& absa16,const real4d_u& ka16,const real2d_u& absb16,const real3d_u& kb16,const real2d_u& forref16,const real2d_u& selfref16,const real1d_u& sfluxref16, real rayl16,int layreffr16, real strrat16,
    const real2d_u& absa17,const real4d_u& ka17,const real2d_u& absb17,const real4d_u& kb17,const real2d_u& forref17,const real2d_u& selfref17,const real2d_u& sfluxref17, real rayl17,int layreffr17, real strrat17,
    const real2d_u& absa18,const real4d_u& ka18,const real2d_u& absb18,const real3d_u& kb18,const real2d_u& forref18,const real2d_u& selfref18,const real2d_u& sfluxref18, real rayl18,int layreffr18, real strrat18,
    const real2d_u& absa19,const real4d_u& ka19,const real2d_u& absb19,const real3d_u& kb19,const real2d_u& forref19,const real2d_u& selfref19,const real2d_u& sfluxref19, real rayl19,int layreffr19, real strrat19,
    const real2d_u& absa20,const real3d_u& ka20,const real2d_u& absb20,const real3d_u& kb20,const real2d_u& forref20,const real2d_u& selfref20,const real1d_u& sfluxref20,const real1d_u& absch420, real rayl20,int layreffr20,
    const real2d_u& absa21,const real4d_u& ka21,const real2d_u& absb21,const real4d_u& kb21,const real2d_u& forref21,const real2d_u& selfref21,const real2d_u& sfluxref21, real rayl21,int layreffr21, real strrat21,
    const real2d_u& absa22,const real4d_u& ka22,const real2d_u& absb22,const real3d_u& kb22,const real2d_u& forref22,const real2d_u& selfref22,const real2d_u& sfluxref22, real rayl22,int layreffr22, real strrat22,
    const real2d_u& absa23,const real3d_u& ka23,const real2d_u& forref23,const real2d_u& selfref23,const real1d_u& sfluxref23,const real1d_u& rayl23,int layreffr23, real givfac23,
    const real2d_u& absa24,const real4d_u& ka24,const real2d_u& absb24,const real3d_u& kb24,const real2d_u& forref24,const real2d_u& selfref24,const real2d_u& sfluxref24,const real1d_u& abso3a24,const real1d_u& abso3b24,const real2d_u& rayla24,const real1d_u& raylb24,int layreffr24, real strrat24,
    const real2d_u& absa25,const real3d_u& ka25,const real1d_u& sfluxref25,const real1d_u& abso3a25,const real1d_u& abso3b25,const real1d_u& rayl25,int layreffr25,
    const real1d_u& sfluxref26,const real1d_u& rayl26,
    const real2d_u& absa27,const real3d_u& ka27,const real2d_u& absb27,const real3d_u& kb27,const real1d_u& sfluxref27,const real1d_u& rayl27,int layreffr27,real scalekur27,
    const real2d_u& absa28,const real4d_u& ka28,const real2d_u& absb28,const real4d_u& kb28,const real2d_u& sfluxref28,real rayl28,int layreffr28,real strrat28,
    const real2d_u& absa29,const real3d_u& ka29,const real2d_u& absb29,const real3d_u& kb29,const real2d_u& forref29,const real2d_u& selfref29,const real1d_u& sfluxref29,const real1d_u& absh2o29,const real1d_u& absco229,
    real rayl29,int layreffr29,const real4d& cldfmc,const real4d& ciwpmc,const real4d& clwpmc,const real4d& cswpmc,const real4d& taucmc,
    const real4d& taormc,const real4d& ssacmc,const real4d& asmcmc,const real4d& fsfcmc,const real4d& zcldfmc,const real4d& ztaucmc,const real4d& ztaormc,const real4d& zasycmc,const real4d& zomgcmc,
    const real4d& taua,const real4d& ssaa,const real4d& asma,const real4d& wkl,const real4d& ztaug,const real4d& ztaur,const real4d& zcd,const real4d& zcu,const real4d& zfd,const real4d& zfu,
    const real3d &zdbt,const real3d &zdbt_nodel,const real3d &zgc,const real3d &zgcc,const real3d &zgco,const real3d &zomc,const real3d &zomcc,const real3d &zomco,const real3d &zrdnd,const real3d &zrdndc,const real3d &zref,
    const real3d &zrefc,const real3d &zrefo,const real3d &zrefd,const real3d &zrefdc,const real3d &zrefdo,const real3d &zrup,const real3d &zrupd,const real3d &zrupc,const real3d &zrupdc,const real3d &zs1,const real3d &ztauc,const real3d &ztauo,const real3d &ztdn,const real3d &ztdnd,const real3d &ztdbt,const real3d &ztoc,const real3d &ztor,const real3d &ztra,const real3d &ztrac,const real3d &ztrao,
    const real3d &ztrad,const real3d &ztradc,const real3d &ztrado,const real3d &zdbtc,const real3d &ztdbtc,const real3d &zincflx,const real3d &zdbtc_nodel,const real3d &ztdbt_nodel,const real3d &ztdbtc_nodel,const real3d &zsflxzen,int m,int j);
   

extern "C" void rrtmg_swrad_kokkos(flcl_ndarray_t *p3d_p,flcl_ndarray_t *p8w_p,flcl_ndarray_t *pi3d_p,flcl_ndarray_t *t3d_p,flcl_ndarray_t *t8w_p,flcl_ndarray_t *dz8w_p,flcl_ndarray_t *qv3d_p,flcl_ndarray_t *qc3d_p,flcl_ndarray_t *qr3d_p,       
    flcl_ndarray_t *qi3d_p,flcl_ndarray_t *qs3d_p,flcl_ndarray_t *qg3d_p,flcl_ndarray_t *cldfra3d_p,flcl_ndarray_t *o33d_p,flcl_ndarray_t *tsk_p,flcl_ndarray_t *albedo_p,        
    flcl_ndarray_t *xland_p,flcl_ndarray_t *xice_p,flcl_ndarray_t *snow_p,flcl_ndarray_t *coszr_p,flcl_ndarray_t *hrang_p,real& xtime,real& gmt,int& julday,real& radt,    
    real& degrad,real& declin,real& solcon,flcl_ndarray_t *xlat_p,flcl_ndarray_t *xlong_p,int& icloud,int& o3input, 
    int& noznlevels,flcl_ndarray_t *pin_p,flcl_ndarray_t *o3clim_p,flcl_ndarray_t *gsw_p,flcl_ndarray_t *swcf_p,                 
    int& config_aer_opt_p,bool& config_aer_ra_feedback_p,          
    char* config_chem_opt_p,                                
    flcl_ndarray_t *rthratensw_p,int& has_reqc,int& has_reqi,int& has_reqs,flcl_ndarray_t *re_cloud_p, 
    flcl_ndarray_t *re_ice_p,flcl_ndarray_t *re_snow_p,                                 
    flcl_ndarray_t *swupt_p,flcl_ndarray_t *swuptc_p,flcl_ndarray_t *swdnt_p,flcl_ndarray_t *swdntc_p,                   
    flcl_ndarray_t *swupb_p,flcl_ndarray_t *swupbc_p,flcl_ndarray_t *swdnb_p,flcl_ndarray_t *swdnbc_p,                   
    flcl_ndarray_t *swupflx_p,flcl_ndarray_t *swupflxc_p,flcl_ndarray_t *swdnflx_p,flcl_ndarray_t *swdnflxc_p,           
    flcl_ndarray_t *tauaersw_p,flcl_ndarray_t *gaersw_p,flcl_ndarray_t *waersw_p,                         
    int& ids,int& ide,int& jds,int& jde,int& kds,int& kde,                       
    int& ims,int& ime,int& jms,int& jme,int& kms,int& kme,                      
    int& its,int& ite,int& jts,int& jte,int& kts,int& kte,                       
    flcl_ndarray_t *wavemin_p,flcl_ndarray_t *wavemax_p,flcl_ndarray_t *retab_p,
    flcl_ndarray_t *absa16_p,flcl_ndarray_t *ka16_p,flcl_ndarray_t *absb16_p,flcl_ndarray_t *kb16_p,flcl_ndarray_t *forref16_p,flcl_ndarray_t *selfref16_p,flcl_ndarray_t *sfluxref16_p,real& rayl16,int& layreffr16,real& strrat16,
    flcl_ndarray_t *absa17_p,flcl_ndarray_t *ka17_p,flcl_ndarray_t *absb17_p,flcl_ndarray_t *kb17_p,flcl_ndarray_t *forref17_p,flcl_ndarray_t *selfref17_p,flcl_ndarray_t *sfluxref17_p,real& rayl17,int& layreffr17,real& strrat17,
    flcl_ndarray_t *absa18_p,flcl_ndarray_t *ka18_p,flcl_ndarray_t *absb18_p,flcl_ndarray_t *kb18_p,flcl_ndarray_t *forref18_p,flcl_ndarray_t *selfref18_p,flcl_ndarray_t *sfluxref18_p,real& rayl18,int& layreffr18,real& strrat18,
    flcl_ndarray_t *absa19_p,flcl_ndarray_t *ka19_p,flcl_ndarray_t *absb19_p,flcl_ndarray_t *kb19_p,flcl_ndarray_t *forref19_p,flcl_ndarray_t *selfref19_p,flcl_ndarray_t *sfluxref19_p,real& rayl19,int& layreffr19,real& strrat19,
    flcl_ndarray_t *absa20_p,flcl_ndarray_t *ka20_p,flcl_ndarray_t *absb20_p,flcl_ndarray_t *kb20_p,flcl_ndarray_t *forref20_p,flcl_ndarray_t *selfref20_p,flcl_ndarray_t *sfluxref20_p,flcl_ndarray_t *absch420_p,real& rayl20,int& layreffr20,
    flcl_ndarray_t *absa21_p,flcl_ndarray_t *ka21_p,flcl_ndarray_t *absb21_p,flcl_ndarray_t *kb21_p,flcl_ndarray_t *forref21_p,flcl_ndarray_t *selfref21_p,flcl_ndarray_t *sfluxref21_p,real& rayl21,int& layreffr21,real& strrat21,
    flcl_ndarray_t *absa22_p,flcl_ndarray_t *ka22_p,flcl_ndarray_t *absb22_p,flcl_ndarray_t *kb22_p,flcl_ndarray_t *forref22_p,flcl_ndarray_t *selfref22_p,flcl_ndarray_t *sfluxref22_p,real& rayl22,int& layreffr22,real& strrat22,
    flcl_ndarray_t *absa23_p,flcl_ndarray_t *ka23_p,flcl_ndarray_t *forref23_p,flcl_ndarray_t *selfref23_p,flcl_ndarray_t *sfluxref23_p,flcl_ndarray_t *rayl23_p,int& layreffr23,real& givfac23,
    flcl_ndarray_t *absa24_p,flcl_ndarray_t *ka24_p,flcl_ndarray_t *absb24_p,flcl_ndarray_t *kb24_p,flcl_ndarray_t *forref24_p,flcl_ndarray_t *selfref24_p,flcl_ndarray_t *sfluxref24_p,flcl_ndarray_t *abso3a24_p,flcl_ndarray_t *abso3b24_p,flcl_ndarray_t *rayla24_p,flcl_ndarray_t *raylb24_p,int& layreffr24,real& strrat24,
    flcl_ndarray_t *absa25_p,flcl_ndarray_t *ka25_p,flcl_ndarray_t *sfluxref25_p,flcl_ndarray_t *abso3a25_p,flcl_ndarray_t *abso3b25_p,flcl_ndarray_t *rayl25_p,int& layreffr25,
    flcl_ndarray_t *sfluxref26_p,flcl_ndarray_t *rayl26_p,
    flcl_ndarray_t *absa27_p,flcl_ndarray_t *ka27_p,flcl_ndarray_t *absb27_p,flcl_ndarray_t *kb27_p,flcl_ndarray_t *sfluxref27_p,flcl_ndarray_t *rayl27_p,int& layreffr27,real& scalekur27,
    flcl_ndarray_t *absa28_p,flcl_ndarray_t *ka28_p,flcl_ndarray_t *absb28_p,flcl_ndarray_t *kb28_p,flcl_ndarray_t *sfluxref28_p,real& rayl28,int& layreffr28,real& strrat28,
    flcl_ndarray_t *absa29_p,flcl_ndarray_t *ka29_p,flcl_ndarray_t *absb29_p,flcl_ndarray_t *kb29_p,flcl_ndarray_t *forref29_p,flcl_ndarray_t *selfref29_p,flcl_ndarray_t *sfluxref29_p,flcl_ndarray_t *absh2o29_p,flcl_ndarray_t *absco229_p,real& rayl29,int& layreffr29){
    real1d_u pin;
    if(pin_p->data != nullptr){auto pin_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*pin_p); pin = real1d_u(pin_raw,{1});}
    auto wavemin_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*wavemin_p); auto wavemin = real1d_u(wavemin_raw,{1});
    auto wavemax_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*wavemax_p); auto wavemax = real1d_u(wavemax_raw,{1});
    auto retab_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*retab_p); auto retab = real1d_u(retab_raw,{1});
    auto hrang_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*hrang_p); auto hrang = real1d_u(hrang_raw,{its});
    auto xlat_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*xlat_p); auto xlat = real2d_u(xlat_raw,{ims,jms});
    auto xlong_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*xlong_p); auto xlong = real2d_u(xlong_raw,{ims,jms});
    auto albedo_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*albedo_p); auto albedo = real2d_u(albedo_raw,{ims,jms});
    auto tsk_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*tsk_p); auto tsk = real2d_u(tsk_raw,{ims,jms});
    auto snow_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*snow_p); auto snow = real2d_u(snow_raw,{ims,jms});
    auto xice_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*xice_p); auto xice = real2d_u(xice_raw,{ims,jms});
    auto xland_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*xland_p); auto xland = real2d_u(xland_raw,{ims,jms});
    auto coszr_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*coszr_p); auto coszr = real2d_u(coszr_raw,{ims,jms});
    auto gsw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*gsw_p); auto gsw = real2d_u(gsw_raw,{ims,jms});
    auto swcf_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swcf_p); auto swcf = real2d_u(swcf_raw,{ims,jms});
    real2d_u swupt;
    if(swupt_p->data != nullptr){auto swupt_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swupt_p); swupt = real2d_u(swupt_raw,{ims,jms});}
    real2d_u swuptc;
    if(swuptc_p->data != nullptr){auto swuptc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swuptc_p); swuptc = real2d_u(swuptc_raw,{ims,jms});}
    real2d_u swdnt;
    if(swdnt_p->data != nullptr){auto swdnt_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swdnt_p); swdnt = real2d_u(swdnt_raw,{ims,jms});}
    real2d_u swdntc;
    if(swdntc_p->data != nullptr){auto swdntc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swdntc_p); swdntc = real2d_u(swdntc_raw,{ims,jms});}
    real2d_u swupb;
    if(swupb_p->data != nullptr){auto swupb_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swupb_p);  swupb = real2d_u(swupb_raw,{ims,jms});}
    real2d_u swupbc;
    if(swupbc_p->data != nullptr){auto swupbc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swupbc_p); swupbc = real2d_u(swupbc_raw,{ims,jms});}
    real2d_u swdnb;
    if(swdnb_p->data != nullptr){auto swdnb_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swdnb_p); swdnb = real2d_u(swdnb_raw,{ims,jms});}
    real2d_u swdnbc;
    if(swdnbc_p->data != nullptr){auto swdnbc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*swdnbc_p); swdnbc = real2d_u(swdnbc_raw,{ims,jms});}
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
    real3d_u o3clim;
    if(o3clim_p->data != nullptr){auto o3clim_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*o3clim_p); o3clim = real3d_u(o3clim_raw,{ims,1,jms});}
    auto rthratensw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*rthratensw_p); auto rthratensw = real3d_u(rthratensw_raw,{ims,kms,jms});
    real3d_u swupflx;
    if(swupflx_p->data != nullptr){auto swupflx_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*swupflx_p); swupflx = real3d_u(swupflx_raw,{ims,kms,jms});}
    real3d_u swupflxc;
    if(swupflxc_p->data != nullptr){auto swupflxc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*swupflxc_p); swupflxc = real3d_u(swupflxc_raw,{ims,kms,jms});}
    real3d_u swdnflx;
    if(swdnflx_p->data != nullptr){auto swdnflx_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*swdnflx_p); swdnflx = real3d_u(swdnflx_raw,{ims,kms,jms});}
    real3d_u swdnflxc;
    if(swdnflxc_p->data != nullptr){auto swdnflxc_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*swdnflxc_p); swdnflxc = real3d_u(swdnflxc_raw,{ims,kms,jms});}
    real4d_u tauaersw;
    if(tauaersw_p->data != nullptr){auto tauaersw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*tauaersw_p); tauaersw = real4d_u(tauaersw_raw,{ims,kms,jms,1});}
    real4d_u gaersw;
    if(gaersw_p->data != nullptr){auto gaersw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*gaersw_p); gaersw = real4d_u(gaersw_raw,{ims,kms,jms,1});}
    real4d_u waersw;
    if(waersw_p->data != nullptr){auto waersw_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*waersw_p); waersw = real4d_u(waersw_raw,{ims,kms,jms,1});}
    //local variables and arrays:
    
    #ifdef USE_SERIALIZATION
    #warning "Using serial"
    auto tsfc = real3d("tsfc",{1,2},{1,2},{1,2});
    auto landfrac = real3d("landfrac",{1,2},{1,2},{1,2});
    auto landm = real3d("landm",{1,2},{1,2},{1,2});
    auto snowh = real3d("snowh",{1,2},{1,2},{1,2});
    auto icefrac = real3d("icefrac",{1,2},{1,2},{1,2});
    auto asdir = real3d("asdir",{1,2},{1,2},{1,2});
    auto asdif = real3d("asdif",{1,2},{1,2},{1,2});
    auto aldir = real3d("aldir",{1,2},{1,2},{1,2});
    auto aldif = real3d("aldif",{1,2},{1,2},{1,2});
    auto coszen = real3d("coszen",{1,2},{1,2},{1,2});
    auto pdel = real4d("pdel",{1,2},{1,kte-kts+2},{1,2},{1,2});
    auto cicewp = real4d("cicewp",{1,2},{1,kte-kts+2},{1,2},{1,2});
    auto cliqwp = real4d("cliqwp",{1,2},{1,kte-kts+2},{1,2},{1,2});
    auto csnowp = real4d("csnowp",{1,2},{1,kte-kts+2},{1,2},{1,2});
    auto reliq = real4d("reliq",{1,2},{1,kte-kts+2},{1,2},{1,2});
    auto reice = real4d("reice",{1,2},{1,kte-kts+2},{1,2},{1,2});
    auto resnow = real4d("resnow",{1,2},{1,kte-kts+2},{1,2},{1,2});
    auto tten1d = real3d("tten1d",{kts,kte+1},{1,2},{1,2});
    auto cldfra1d = real3d("cldfra1d",{kts,kte+1},{1,2},{1,2});
    auto dz1d = real3d("dz1d",{kts,kte+1},{1,2},{1,2});
    auto p1d = real3d("p1d",{kts,kte+1},{1,2},{1,2});
    auto t1d = real3d("t1d",{kts,kte+1},{1,2},{1,2});
    auto qv1d = real3d("qv1d",{kts,kte+1},{1,2},{1,2});
    auto qc1d = real3d("qc1d",{kts,kte+1},{1,2},{1,2});
    auto qr1d = real3d("qr1d",{kts,kte+1},{1,2},{1,2});
    auto qi1d = real3d("qi1d",{kts,kte+1},{1,2},{1,2});
    auto qs1d = real3d("qs1d",{kts,kte+1},{1,2},{1,2});
    auto qg1d = real3d("qg1d",{kts,kte+1},{1,2},{1,2});
    auto o31d = real3d("o31d",{kts,kte+1},{1,2},{1,2});
    auto pw1d = real3d("pw1d",{kts,kte+2},{1,2},{1,2});
    auto tw1d = real3d("tw1d",{kts,kte+2},{1,2},{1,2});
    auto o3mmr = real3d("o3mmr",{kts,kte+2},{1,2},{1,2});
    //--- additional local variables and arrays needed to include additional layers between the model top
    //    and the top of the atmosphere:
    auto play = real4d("play",{1,2},{kts,kte+2},{1,2},{1,2});
    auto tlay = real4d("tlay",{1,2},{kts,kte+2},{1,2},{1,2});
    auto h2ovmr = real4d("h2ovmr",{1,2},{kts,kte+2},{1,2},{1,2});
    auto o3vmr = real4d("o3vmr",{1,2},{kts,kte+2},{1,2},{1,2});
    auto co2vmr = real4d("co2vmr",{1,2},{kts,kte+2},{1,2},{1,2});
    auto o2vmr = real4d("o2vmr",{1,2},{kts,kte+2},{1,2},{1,2});
    auto ch4vmr = real4d("ch4vmr",{1,2},{kts,kte+2},{1,2},{1,2});
    auto n2ovmr = real4d("n2ovmr",{1,2},{kts,kte+2},{1,2},{1,2});
    auto clwpth = real4d("clwpth",{1,2},{kts,kte+2},{1,2},{1,2});
    auto ciwpth = real4d("ciwpth",{1,2},{kts,kte+2},{1,2},{1,2});
    auto cswpth = real4d("cswpth",{1,2},{kts,kte+2},{1,2},{1,2});
    auto rel = real4d("rel",{1,2},{kts,kte+2},{1,2},{1,2});
    auto rei = real4d("rei",{1,2},{kts,kte+2},{1,2},{1,2});
    auto res = real4d("res",{1,2},{kts,kte+2},{1,2},{1,2});
    auto cldfrac = real4d("cldfrac",{1,2},{kts,kte+2},{1,2},{1,2});
    auto relqmcl = real4d("relqmcl",{1,2},{kts,kte+2},{1,2},{1,2});
    auto reicmcl = real4d("reicmcl",{1,2},{kts,kte+2},{1,2},{1,2});
    auto resnmcl = real4d("resnmcl",{1,2},{kts,kte+2},{1,2},{1,2});
    auto swhr = real4d("swhr",{1,2},{kts,kte+2},{1,2},{1,2});
    auto swhrc = real4d("swhrc",{1,2},{kts,kte+2},{1,2},{1,2});
    auto plev = real4d("plev",{1,2},{kts,kte+3},{1,2},{1,2});
    auto tlev = real4d("tlev",{1,2},{kts,kte+3},{1,2},{1,2});
    auto swuflx = real4d("swuflx",{1,2},{kts,kte+3},{1,2},{1,2});
    auto swdflx = real4d("swdflx",{1,2},{kts,kte+3},{1,2},{1,2});
    auto swuflxc = real4d("swuflxc",{1,2},{kts,kte+3},{1,2},{1,2});
    auto swdflxc = real4d("swdflxc",{1,2},{kts,kte+3},{1,2},{1,2});
    auto sibvisdir = real4d("sibvisdir",{1,2},{kts,kte+3},{1,2},{1,2});
    auto sibvisdif = real4d("sibvisdif",{1,2},{kts,kte+3},{1,2},{1,2});
    auto sibnirdir = real4d("sibnirdir",{1,2},{kts,kte+3},{1,2},{1,2});
    auto sibnirdif = real4d("sibnirdif",{1,2},{kts,kte+3},{1,2},{1,2});
    auto swdkdir = real4d("swdkdir",{1,2},{kts,kte+3},{1,2},{1,2});
    auto swdkdif = real4d("swdkdif",{1,2},{kts,kte+3},{1,2},{1,2});
    auto tauaer = real5d("tauaer",{1,2},{kts,kte+2},{1,nbndsw+1},{1,2},{1,2});
    auto ssaaer = real5d("ssaaer",{1,2},{kts,kte+2},{1,nbndsw+1},{1,2},{1,2});
    auto asmaer = real5d("asmaer",{1,2},{kts,kte+2},{1,nbndsw+1},{1,2},{1,2});
    auto taucld = real5d("taucld",{1,nbndsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto ssacld = real5d("ssacld",{1,nbndsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto asmcld = real5d("asmcld",{1,nbndsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto fsfcld = real5d("fsfcld",{1,nbndsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto cldfmcl = real5d("cldfmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto clwpmcl = real5d("clwpmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto ciwpmcl = real5d("ciwpmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto cswpmcl = real5d("cswpmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto taucmcl = real5d("taucmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto ssacmcl = real5d("ssacmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto asmcmcl = real5d("asmcmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto fsfcmcl = real5d("fsfcmcl",{1,ngptsw+1},{1,2},{kts,kte+2},{1,2},{1,2});
    auto cldfmc = real4d("cldfmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! cloud fraction [mcica]
    auto ciwpmc = real4d("ciwpmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! in-cloud ice water path [mcica]
    auto clwpmc = real4d("clwpmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! in-cloud liquid water path [mcica]
    auto cswpmc = real4d("cswpmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! in-cloud snow water path [mcica]
    auto taucmc = real4d("taucmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! in-cloud optical depth [mcica]
    auto taormc = real4d("taormc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! unscaled in-cloud optical depth [mcica]s
    auto ssacmc = real4d("ssacmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! in-cloud single scattering albedo [mcica]s
    auto asmcmc = real4d("asmcmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! in-cloud asymmetry parameter [mcica]
    auto fsfcmc = real4d("fsfcmc",{1,ngptsw+1},{1,kte+2},{1,2},{1,2});              //! in-cloud forward scattering fraction [mcica]  
    auto zcldfmc = real4d("zcldfmc",{1,kte+2},{1,ngptsw+1},{1,2},{1,2});            //! cloud fraction [mcica]
    auto ztaucmc = real4d("ztaucmc",{1,kte+2},{1,ngptsw+1},{1,2},{1,2});            //! cloud optical depth [mcica]
    auto ztaormc = real4d("ztaormc",{1,kte+2},{1,ngptsw+1},{1,2},{1,2});            //! unscaled cloud optical depth [mcica]
    auto zasycmc = real4d("zasycmc",{1,kte+2},{1,ngptsw+1},{1,2},{1,2});            //! cloud asymmetry parameter [mcica]
    auto zomgcmc = real4d("zomgcmc",{1,kte+2},{1,ngptsw+1},{1,2},{1,2});            //! cloud single scattering albedo [mcica]
    auto taua = real4d("taua",{1,kte+2},{1,nbndsw+1},{1,2},{1,2});                   //! Aerosol optical depth
    auto ssaa = real4d("ssaa",{1,kte+2},{1,nbndsw+1},{1,2},{1,2});                   //! Aerosol single scattering albedo                
    auto asma = real4d("asma",{1,kte+2},{1,nbndsw+1},{1,2},{1,2});                   //! Aerosol asymmetry parameter  
    // auto ztauc = real4d("ztauc",{1,kte+2},{1,nbndsw+1},{1,2},{1,2});                //! cloud optical depth
    // auto ztaucorig = real4d("ztaucorig",{1,kte+2},{1,nbndsw+1},{1,2},{1,2});        //! unscaled cloud optical depth
    // auto zasyc = real4d("zasyc",{1,kte+2},{1,nbndsw+1},{1,2},{1,2});               //! cloud asymmetry parameter  (first moment of phase function)
    // auto zomgc = real4d("zomgc",{1,kte+2},{1,nbndsw+1},{1,2},{1,2});                //! cloud single scattering albedo
    auto wkl = real4d("wkl",{1,mxmol+1},{1,kte+2},{1,2},{1,2});                 //! molecular amounts (mol/cm-2)
    auto ztaug = real4d("ztaug",{1,SIZE_NLAY+1},{1,ngptsw + 1},{1,2},{1,2});
    auto ztaur = real4d("ztaur",{1,SIZE_NLAY+1},{1,ngptsw + 1},{1,2},{1,2});
    auto zcd = real4d("zcd",{1,SIZE_NLAY+1},{1,ngptsw + 2},{1,2},{1,2});
    auto zcu = real4d("zcu",{1,SIZE_NLAY+1},{1,ngptsw + 2},{1,2},{1,2});
    auto zfd = real4d("zfd",{1,SIZE_NLAY+1},{1,ngptsw + 2},{1,2},{1,2});
    auto zfu = real4d("zfu",{1,SIZE_NLAY+1},{1,ngptsw + 2},{1,2},{1,2});
    auto zdbt = real3d("zdbt",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zdbt_nodel = real3d("zdbt_nodel",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zgc = real3d("zgc",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto zgcc = real3d("zgcc",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto zgco = real3d("zgco",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto zomc = real3d("zomc",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto zomcc = real3d("zomcc",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto zomco = real3d("zomco",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto zrdnd = real3d("zrdnd",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrdndc = real3d("zrdndc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zref = real3d("zref",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrefc = real3d("zrefc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrefo = real3d("zrefo",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrefd = real3d("zrefd",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrefdc = real3d("zrefdc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrefdo = real3d("zrefdo",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrup = real3d("zrup",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrupd = real3d("zrupd",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrupc = real3d("zrupc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zrupdc = real3d("zrupdc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zs1 = real3d("zs1",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztauc = real3d("ztauc",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto ztauo = real3d("ztauo",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto ztdn = real3d("ztdn",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztdnd = real3d("ztdnd",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztdbt = real3d("ztdbt",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztoc = real3d("ztoc",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto ztor = real3d("ztor",{1,SIZE_NLAY+1},{1,2},{1,2});
    auto ztra = real3d("ztra",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztrac = real3d("ztrac",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztrao = real3d("ztrao",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztrad = real3d("ztrad",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztradc = real3d("ztradc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztrado = real3d("ztrado",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zdbtc = real3d("zdbtc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztdbtc = real3d("ztdbtc",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zincflx = real3d("zincflx",{1,ngptsw+1},{1,2},{1,2});
    auto zdbtc_nodel = real3d("zdbtc_nodel",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztdbt_nodel = real3d("ztdbt_nodel",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto ztdbtc_nodel = real3d("ztdbtc_nodel",{1,SIZE_NLAY+2},{1,2},{1,2});
    auto zsflxzen = real3d("zsflxzen",{1,ngptsw+1},{1,2},{1,2});
    
    //--- additional local variables and arrays needed for the CAM ozone climatologyL
    auto o3clim1d = real3d("o3clim1d",{1,noznlevels + 1},{1,2},{1,2});
    //--- additional local variables related to the implementation of aerosols in rrtmg_swrad in WRF 3.8.
    //    In WRF 3.8, these variables are in the argument list of subroutine rrtmg_swrad, but are made
    //    local here:
    
    auto ecaer = real5d("ecaer",{1,2},{kts,kte+2},{1,naerec + 1},{1,2},{1,2});
    // jwfeng from AEMOL USTC: 2020/10/10
    // wavelength corresponding to wavenum1 and wavenum2 (cm-1)
    auto wavemid = real3d("wavemid",{1,nbndsw + 1},{1,2},{1,2});
    #endif
    #ifdef USE_SW
    #warning "Using SW"
        sw_workspace &ws = get_sw_ws(kts, kte, its, RAD_SW_CHUNK_SIZE, jts, jte,
                                     nbndsw, ngptsw, naerec, mxmol, noznlevels, SIZE_NLAY);
        auto tsfc       = ws.tsfc;
        auto landfrac   = ws.landfrac;
        auto landm      = ws.landm;
        auto snowh      = ws.snowh;
        auto icefrac    = ws.icefrac;
        auto asdir      = ws.asdir;
        auto asdif      = ws.asdif;
        auto aldir      = ws.aldir;
        auto aldif      = ws.aldif;
        auto coszen     = ws.coszen;

        // --- 2. 1D Profiles (expanded to 3D in Kokkos View) ---
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
        auto o3mmr      = ws.o3mmr;

        // --- 3. Atmosphere & Physics (4D) ---
        auto pdel       = ws.pdel;
        auto cicewp     = ws.cicewp;
        auto cliqwp     = ws.cliqwp;
        auto csnowp     = ws.csnowp;
        auto reliq      = ws.reliq;
        auto reice      = ws.reice;
        auto resnow     = ws.resnow;
        
        auto play       = ws.play;
        auto tlay       = ws.tlay;
        auto h2ovmr     = ws.h2ovmr;
        auto o3vmr      = ws.o3vmr;
        auto co2vmr     = ws.co2vmr;
        auto o2vmr      = ws.o2vmr;
        auto ch4vmr     = ws.ch4vmr;
        auto n2ovmr     = ws.n2ovmr;
        
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
        auto swhr       = ws.swhr;
        auto swhrc      = ws.swhrc;

        // --- 4. Fluxes (4D) ---
        auto plev       = ws.plev;
        auto tlev       = ws.tlev;
        auto swuflx     = ws.swuflx;
        auto swdflx     = ws.swdflx;
        auto swuflxc    = ws.swuflxc;
        auto swdflxc    = ws.swdflxc;
        auto sibvisdir  = ws.sibvisdir;
        auto sibvisdif  = ws.sibvisdif;
        auto sibnirdir  = ws.sibnirdir;
        auto sibnirdif  = ws.sibnirdif;
        auto swdkdir    = ws.swdkdir;
        auto swdkdif    = ws.swdkdif;

        // --- 5. High-Dimension Optical Properties (5D - The Memory Killers) ---
        auto tauaer     = ws.tauaer;
        auto ssaaer     = ws.ssaaer;
        auto asmaer     = ws.asmaer;
        auto taucld     = ws.taucld;
        auto ssacld     = ws.ssacld;
        auto asmcld     = ws.asmcld;
        auto fsfcld     = ws.fsfcld;
        
        auto cldfmcl    = ws.cldfmcl;
        auto clwpmcl    = ws.clwpmcl;
        auto ciwpmcl    = ws.ciwpmcl;
        auto cswpmcl    = ws.cswpmcl;
        auto taucmcl    = ws.taucmcl;
        auto ssacmcl    = ws.ssacmcl;
        auto asmcmcl    = ws.asmcmcl;
        auto fsfcmcl    = ws.fsfcmcl;
        
        auto ecaer      = ws.ecaer; // This was the one crashing

        // --- 6. McICA Intermediates (4D) ---
        auto cldfmc     = ws.cldfmc;
        auto ciwpmc     = ws.ciwpmc;
        auto clwpmc     = ws.clwpmc;
        auto cswpmc     = ws.cswpmc;
        auto taucmc     = ws.taucmc;
        auto taormc     = ws.taormc;
        auto ssacmc     = ws.ssacmc;
        auto asmcmc     = ws.asmcmc;
        auto fsfcmc     = ws.fsfcmc;

        auto zcldfmc    = ws.zcldfmc;
        auto ztaucmc    = ws.ztaucmc;
        auto ztaormc    = ws.ztaormc;
        auto zasycmc    = ws.zasycmc;
        auto zomgcmc    = ws.zomgcmc;

        auto taua       = ws.taua;
        auto ssaa       = ws.ssaa;
        auto asma       = ws.asma;
        auto wkl        = ws.wkl;

        // --- 7. Two-Stream Solver Variables (Z-Variables) ---
        // 4D Z-vars
        auto ztaug      = ws.ztaug;
        auto ztaur      = ws.ztaur;
        auto zcd        = ws.zcd;
        auto zcu        = ws.zcu;
        auto zfd        = ws.zfd;
        auto zfu        = ws.zfu;

        // 3D Z-vars
        auto zdbt       = ws.zdbt;
        auto zdbt_nodel = ws.zdbt_nodel;
        auto zgc        = ws.zgc;
        auto zgcc       = ws.zgcc;
        auto zgco       = ws.zgco;
        auto zomc       = ws.zomc;
        auto zomcc      = ws.zomcc;
        auto zomco      = ws.zomco;
        auto zrdnd      = ws.zrdnd;
        auto zrdndc     = ws.zrdndc;
        auto zref       = ws.zref;
        auto zrefc      = ws.zrefc;
        auto zrefo      = ws.zrefo;
        auto zrefd      = ws.zrefd;
        auto zrefdc     = ws.zrefdc;
        auto zrefdo     = ws.zrefdo;
        auto zrup       = ws.zrup;
        auto zrupd      = ws.zrupd;
        auto zrupc      = ws.zrupc;
        auto zrupdc     = ws.zrupdc;
        auto zs1        = ws.zs1;
        auto ztauc      = ws.ztauc;
        auto ztauo      = ws.ztauo;
        auto ztdn       = ws.ztdn;
        auto ztdnd      = ws.ztdnd;
        auto ztdbt      = ws.ztdbt;
        auto ztoc       = ws.ztoc;
        auto ztor       = ws.ztor;
        auto ztra       = ws.ztra;
        auto ztrac      = ws.ztrac;
        auto ztrao      = ws.ztrao;
        auto ztrad      = ws.ztrad;
        auto ztradc     = ws.ztradc;
        auto ztrado     = ws.ztrado;
        auto zdbtc      = ws.zdbtc;
        auto ztdbtc     = ws.ztdbtc;
        
        auto zincflx    = ws.zincflx;
        auto zdbtc_nodel= ws.zdbtc_nodel;
        auto ztdbt_nodel= ws.ztdbt_nodel;
        auto ztdbtc_nodel= ws.ztdbtc_nodel;
        auto zsflxzen   = ws.zsflxzen;

        // --- 8. Others ---
        auto o3clim1d   = ws.o3clim1d;
        auto wavemid    = ws.wavemid;
        
    #endif
    //define kokkos_array in rrsw_kg
    auto absa16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa16_p); auto absa16 = real2d_u(absa16_raw,{1,1});
    auto absa28_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa28_p); auto absa28 = real2d_u(absa28_raw,{1,1});
    auto absa17_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa17_p); auto absa17 = real2d_u(absa17_raw,{1,1});
    auto absa18_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa18_p); auto absa18 = real2d_u(absa18_raw,{1,1});
    auto absa19_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa19_p); auto absa19 = real2d_u(absa19_raw,{1,1});
    auto absa24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa24_p); auto absa24 = real2d_u(absa24_raw,{1,1});
    auto absa21_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa21_p); auto absa21 = real2d_u(absa21_raw,{1,1});
    auto absa22_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa22_p); auto absa22 = real2d_u(absa22_raw,{1,1});
    auto absa20_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa20_p); auto absa20 = real2d_u(absa20_raw,{1,1});
    auto absa23_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa23_p); auto absa23 = real2d_u(absa23_raw,{1,1});
    auto absa25_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa25_p); auto absa25 = real2d_u(absa25_raw,{1,1});
    auto absa27_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa27_p); auto absa27 = real2d_u(absa27_raw,{1,1});
    auto absa29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absa29_p); auto absa29 = real2d_u(absa29_raw,{1,1});
    auto ka16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka16_p); auto ka16 = real4d_u(ka16_raw,{1,1,1,1});
    auto ka28_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka28_p); auto ka28 = real4d_u(ka28_raw,{1,1,1,1});
    auto ka17_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka17_p); auto ka17 = real4d_u(ka17_raw,{1,1,1,1});
    auto ka18_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka18_p); auto ka18 = real4d_u(ka18_raw,{1,1,1,1});
    auto ka19_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka19_p); auto ka19 = real4d_u(ka19_raw,{1,1,1,1});
    auto ka24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka24_p); auto ka24 = real4d_u(ka24_raw,{1,1,1,1});
    auto ka21_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka21_p); auto ka21 = real4d_u(ka21_raw,{1,1,1,1});
    auto ka22_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*ka22_p); auto ka22 = real4d_u(ka22_raw,{1,1,1,1});
    auto ka20_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka20_p); auto ka20 = real3d_u(ka20_raw,{1,1,1});
    auto ka23_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka23_p); auto ka23 = real3d_u(ka23_raw,{1,1,1});
    auto ka25_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka25_p); auto ka25 = real3d_u(ka25_raw,{1,1,1});
    auto ka27_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka27_p); auto ka27 = real3d_u(ka27_raw,{1,1,1});
    auto ka29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*ka29_p); auto ka29 = real3d_u(ka29_raw,{1,1,1});
    auto absb16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb16_p); auto absb16 = real2d_u(absb16_raw,{1,1});
    auto absb18_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb18_p); auto absb18 = real2d_u(absb18_raw,{1,1});
    auto absb19_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb19_p); auto absb19 = real2d_u(absb19_raw,{1,1});
    auto absb24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb24_p); auto absb24 = real2d_u(absb24_raw,{1,1});
    auto absb27_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb27_p); auto absb27 = real2d_u(absb27_raw,{1,1});
    auto absb20_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb20_p); auto absb20 = real2d_u(absb20_raw,{1,1});
    auto absb22_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb22_p); auto absb22 = real2d_u(absb22_raw,{1,1});
    auto absb29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb29_p); auto absb29 = real2d_u(absb29_raw,{1,1});
    auto absb17_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb17_p); auto absb17 = real2d_u(absb17_raw,{1,1});
    auto absb21_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb21_p); auto absb21 = real2d_u(absb21_raw,{1,1});
    auto absb28_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*absb28_p); auto absb28 = real2d_u(absb28_raw,{1,1});
    auto kb16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb16_p); auto kb16 = real3d_u(kb16_raw,{1,13,1});
    auto kb18_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb18_p); auto kb18 = real3d_u(kb18_raw,{1,13,1});
    auto kb19_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb19_p); auto kb19 = real3d_u(kb19_raw,{1,13,1});
    auto kb24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb24_p); auto kb24 = real3d_u(kb24_raw,{1,13,1});
    auto kb27_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb27_p); auto kb27 = real3d_u(kb27_raw,{1,13,1});
    auto kb20_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb20_p); auto kb20 = real3d_u(kb20_raw,{1,13,1});
    auto kb22_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb22_p); auto kb22 = real3d_u(kb22_raw,{1,13,1});
    auto kb29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ***>(*kb29_p); auto kb29 = real3d_u(kb29_raw,{1,13,1});
    auto kb17_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*kb17_p); auto kb17 = real4d_u(kb17_raw,{1,1,13,1});
    auto kb21_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*kb21_p); auto kb21 = real4d_u(kb21_raw,{1,1,13,1});
    auto kb28_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t ****>(*kb28_p); auto kb28 = real4d_u(kb28_raw,{1,1,13,1});
    auto forref16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref16_p); auto forref16 = real2d_u(forref16_raw,{1,1});
    auto forref18_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref18_p); auto forref18 = real2d_u(forref18_raw,{1,1});
    auto forref19_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref19_p); auto forref19 = real2d_u(forref19_raw,{1,1});
    auto forref24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref24_p); auto forref24 = real2d_u(forref24_raw,{1,1});
    auto forref22_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref22_p); auto forref22 = real2d_u(forref22_raw,{1,1});
    auto forref23_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref23_p); auto forref23 = real2d_u(forref23_raw,{1,1});
    auto forref17_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref17_p); auto forref17 = real2d_u(forref17_raw,{1,1});
    auto forref21_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref21_p); auto forref21 = real2d_u(forref21_raw,{1,1});
    auto forref29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref29_p); auto forref29 = real2d_u(forref29_raw,{1,1});
    auto forref20_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*forref20_p); auto forref20 = real2d_u(forref20_raw,{1,1});
    auto selfref16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref16_p); auto selfref16 = real2d_u(selfref16_raw,{1,1});
    auto selfref17_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref17_p); auto selfref17 = real2d_u(selfref17_raw,{1,1});
    auto selfref18_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref18_p); auto selfref18 = real2d_u(selfref18_raw,{1,1});
    auto selfref19_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref19_p); auto selfref19 = real2d_u(selfref19_raw,{1,1});
    auto selfref24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref24_p); auto selfref24 = real2d_u(selfref24_raw,{1,1});
    auto selfref20_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref20_p); auto selfref20 = real2d_u(selfref20_raw,{1,1});
    auto selfref21_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref21_p); auto selfref21 = real2d_u(selfref21_raw,{1,1});
    auto selfref23_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref23_p); auto selfref23 = real2d_u(selfref23_raw,{1,1});
    auto selfref22_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref22_p); auto selfref22 = real2d_u(selfref22_raw,{1,1});
    auto selfref29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*selfref29_p); auto selfref29 = real2d_u(selfref29_raw,{1,1});
    auto sfluxref16_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*sfluxref16_p); auto sfluxref16 = real1d_u(sfluxref16_raw,{1});
    auto sfluxref25_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*sfluxref25_p); auto sfluxref25 = real1d_u(sfluxref25_raw,{1});
    auto sfluxref26_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*sfluxref26_p); auto sfluxref26 = real1d_u(sfluxref26_raw,{1});
    auto sfluxref20_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*sfluxref20_p); auto sfluxref20 = real1d_u(sfluxref20_raw,{1});
    auto sfluxref23_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*sfluxref23_p); auto sfluxref23 = real1d_u(sfluxref23_raw,{1});
    auto sfluxref27_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*sfluxref27_p); auto sfluxref27 = real1d_u(sfluxref27_raw,{1});
    auto sfluxref29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*sfluxref29_p); auto sfluxref29 = real1d_u(sfluxref29_raw,{1});
    auto sfluxref17_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*sfluxref17_p); auto sfluxref17 = real2d_u(sfluxref17_raw,{1,1});
    auto sfluxref28_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*sfluxref28_p); auto sfluxref28 = real2d_u(sfluxref28_raw,{1,1});
    auto sfluxref18_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*sfluxref18_p); auto sfluxref18 = real2d_u(sfluxref18_raw,{1,1});
    auto sfluxref19_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*sfluxref19_p); auto sfluxref19 = real2d_u(sfluxref19_raw,{1,1});
    auto sfluxref24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*sfluxref24_p); auto sfluxref24 = real2d_u(sfluxref24_raw,{1,1});
    auto sfluxref21_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*sfluxref21_p); auto sfluxref21 = real2d_u(sfluxref21_raw,{1,1});
    auto sfluxref22_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*sfluxref22_p); auto sfluxref22 = real2d_u(sfluxref22_raw,{1,1});
    auto absch420_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*absch420_p); auto absch420 = real1d_u(absch420_raw,{1});
    auto rayl23_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rayl23_p); auto rayl23 = real1d_u(rayl23_raw,{1});
    auto abso3a24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*abso3a24_p); auto abso3a24 = real1d_u(abso3a24_raw,{1});
    auto abso3b24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*abso3b24_p); auto abso3b24 = real1d_u(abso3b24_raw,{1});
    auto rayl27_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rayl27_p); auto rayl27 = real1d_u(rayl27_raw,{1});
    auto raylb24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*raylb24_p); auto raylb24 = real1d_u(raylb24_raw,{1});
    auto abso3a25_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*abso3a25_p); auto abso3a25 = real1d_u(abso3a25_raw,{1});
    auto abso3b25_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*abso3b25_p); auto abso3b25 = real1d_u(abso3b25_raw,{1});
    auto rayl25_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rayl25_p); auto rayl25 = real1d_u(rayl25_raw,{1});
    auto rayl26_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*rayl26_p); auto rayl26 = real1d_u(rayl26_raw,{1});
    auto absh2o29_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*absh2o29_p); auto absh2o29 = real1d_u(absh2o29_raw,{1});
    auto absco229_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t *>(*absco229_p); auto absco229 = real1d_u(absco229_raw,{1});
    auto rayla24_raw = flcl::view_from_ndarray<flcl::flcl_view_r64_c_t **>(*rayla24_p); auto rayla24 = real2d_u(rayla24_raw,{1,1});
    const real thresh1 = 1.e-9;
    const real thresh2 = 1.e9;

    //jwfeng: 2020/10/14
    //is_init is a var used to debug, can be deleted later.
    

    // --- set trace gas volume mixing ratios, 2005 values, IPCC (2007):
    // carbon dioxide (379 ppmv)
    real co2 = 379.e-6;
    //methane (1774 ppbv)
    real ch4 = 1774.e-9;
    //nitrous oxide (319 ppbv)
    real n2o = 319.e-9;

    //set oxygen volume mixing ratio (for o2mmr=0.23143):
    real o2 = 0.209488;

    real amdw = 1.607793;   //Molecular weight of dry air / water vapor
    real amdo = 0.603461;   //Molecular weight of dry air / ozone
    real amdo2 = 0.905190;  //Molecular weight of dry air / oxygen
    //--- all fields are ordered vertically from bottom to top (pressures are in mb):
    

    //--- initialize option for the calculation of the cloud optical properties:
    
    // for(int j = jts; j<=jte; j++){
    //     //--- longitude loop:
    //     for(int i = its; i<=ite; i++){
    for (int chunk_start = its; chunk_start <= ite; chunk_start += RAD_SW_CHUNK_SIZE) {
        int current_chunk_end = std::min(chunk_start + RAD_SW_CHUNK_SIZE - 1, ite);
        int current_chunk_len = current_chunk_end - chunk_start + 1;  
        Kokkos::parallel_for(
            "rrtmgsw_kokkos", Kokkos::RangePolicy<>(its,current_chunk_len + 1), KOKKOS_LAMBDA(int offset) {
                //--- calculate the cosine of the solar zenith angle at the current time step to determine if the sun is
                //    above or below the horizon (xt24 is the fractional part of simulation days plus half of radt in 
                //    units of minutes, julian is in days, and radt is in minutes). do not call rrtmg_sw is night-time:
                int m,j,i;
                #ifdef USE_SERIALIZATION
                m = 1; j = 1; i = chunk_start + offset - 1;
                #endif
                #ifdef USE_SW
                m = offset; j =1;
                i = chunk_start + offset - 1;
                #endif
                bool is_init = true;
                int ncol = 1,icld = 2,inflgsw = 2,iceflgsw = 3,liqflgsw = 1;
                bool dorrsw = true;
                real xt24 = fmod(xtime+radt*0.5,1440.0);
                real tloctm = gmt + xt24/60. + xlong(i,j)/15.;                       
                hrang(i) = 15. * (tloctm-12.) * degrad;                                         //qkyang replace hrang with hrang(i)
                real xxlat = xlat(i,j) * degrad;
                real coszrs = sin(xxlat) * sin(declin) + cos(xxlat) * cos(declin) * cos(hrang(i));      //qkyang replace hrang with hrang(i)
                coszr(i,j) = coszrs;
                    
                if(coszrs<=0.0) dorrsw = false;
                    
                if(dorrsw){
                    //--- INITIALIZE COLUMN SOUNDING (the call to the short wave radiation code is done one column
                    //    at a time):
                    for( int k = kts; k<=kte+1; k++){
                        pw1d(k,m,j) = p8w(i,k,j) / 100.;
                        tw1d(k,m,j) = t8w(i,k,j);
                    }

                    for( int k = kts; k<=kte; k++){
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

                    //initialize the local arrays containing the different cloud water and ice condenstates:
                    if(icloud > 0){
                        for(int k = kts; k<=kte; k++){
                            if(qc3d_p->data != nullptr) qc1d(k,m,j) = std::max(qc3d(i,k,j),0.0);
                            if(qr3d_p->data != nullptr) qr1d(k,m,j) = std::max(qr3d(i,k,j),0.0);
                            if(qi3d_p->data != nullptr) qi1d(k,m,j) = std::max(qi3d(i,k,j),0.0);
                            if(qs3d_p->data != nullptr) qs1d(k,m,j) = std::max(qs3d(i,k,j),0.0);
                            if((qc3d_p->data != nullptr || qi3d_p->data != nullptr || qs3d_p->data != nullptr) && cldfra3d_p->data != nullptr) cldfra1d(k,m,j) = cldfra3d(i,k,j);
                        }
                    }

                    //--- initialize the local radiative heating rate:
                    for(int k = kts; k<=kte; k++) tten1d(k,m,j) = 0.;

                    //--- add extra layer to include absorption between the top of the model and the top of the atmosphere:
                    int nlay = (kte-kts+1)+1;
                    
                    //--- initialize local arrays called in the calculation of the cloud optical properties and radiative
                    //    fluxes:
                    for(int n = 1; n<=ncol; n++){
                        for(int k = kts; k<=kte; k++){
                            reliq(ncol,k,m,j)  = 10.;
                            reice(ncol,k,m,j)  = 10.;
                            resnow(ncol,k,m,j) = 10.;
                            cliqwp(ncol,k,m,j) = 0.;
                            cicewp(ncol,k,m,j) = 0.;
                            csnowp(ncol,k,m,j) = 0.;
                        }

                        for(int k = 1; k<=nlay; k++){
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
                            swuflx(n,k,m,j)  = 0.;
                            swuflxc(n,k,m,j) = 0.;
                            swdflx(n,k,m,j)  = 0.;
                            swdflxc(n,k,m,j) = 0.;
                            swhr(n,k,m,j)    = 0.;
                            swhrc(n,k,m,j)   = 0.;
                            for(int nbnd = 1; nbnd<=nbndsw; nbnd++){
                                taucld(nbnd,n,k,m,j)  = 0.;
                                tauaer(n,k,nbnd,m,j)  = 0.;
                                ssaaer(n,k,nbnd,m,j)  = 0.;
                                asmaer(n,k,nbnd,m,j)  = 0.;
                            }
                            for(int ngpt = 1; ngpt<=ngptsw; ngpt++){
                                cldfmcl(ngpt,n,k,m,j) = 0.;
                                clwpmcl(ngpt,n,k,m,j) = 0.;
                                ciwpmcl(ngpt,n,k,m,j) = 0.;
                                cswpmcl(ngpt,n,k,m,j) = 0.;
                                taucmcl(ngpt,n,k,m,j) = 0.;
                            }
                        }

                        for(int k = 1; k<=nlay+1; k++){
                            sibvisdir(ncol,k,m,j) = 0.;
                            sibvisdif(ncol,k,m,j) = 0.;
                            sibnirdir(ncol,k,m,j) = 0.;
                            sibnirdif(ncol,k,m,j) = 0.;
                            swdkdir(n,k,m,j)      = 0.;
                            swdkdif(n,k,m,j)      = 0.;
                        }

                        swuflx(n,nlay+1,m,j)  = 0.;
                        swuflxc(n,nlay+1,m,j) = 0.;
                        swdflx(n,nlay+1,m,j)  = 0.;
                        swdflxc(n,nlay+1,m,j) = 0.;
                    }

                    // !--- initialization of aerosol optical properties:
                    // !czhao add namelist
                    // !aer_opt = 0
                    int aer_opt = config_aer_opt_p;
                    for(int n = 1; n<=ncol; n++){
                        for(int k = 1; k<=nlay; k++){
                            for(int na = 1; na<=naerec; na++) ecaer(n,k,na,m,j) = 0;
                        }
                    }

                    // !--- fill local arrays with input sounding. convert water vapor mass mixing ratio to
                    // !    volume mixing ratio:
                    plev(ncol,1,m,j) = pw1d(1,m,j);
                    tlev(ncol,1,m,j) = tw1d(1,m,j);
                    tsfc(ncol,m,j)   = tsk(i,j);
                    for( int k = kts; k<=kte; k++){
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
                    }
                    plev(ncol,kte+2,m,j)   = 1.0e-5;
                    tlev(ncol,kte+2,m,j)   = tlev(ncol,kte+1,m,j) + 0.0;
                    play(ncol,kte+1,m,j)   = 0.5 * plev(ncol,kte+1,m,j);
                    tlay(ncol,kte+1,m,j)   = tlev(ncol,kte+1,m,j) + 0.0;
                    h2ovmr(ncol,kte+1,m,j) = h2ovmr(ncol,kte,m,j);
                    co2vmr(ncol,kte+1,m,j) = co2vmr(ncol,kte,m,j); 
                    o2vmr(ncol,kte+1,m,j)  = o2vmr(ncol,kte,m,j); 
                    ch4vmr(ncol,kte+1,m,j) = ch4vmr(ncol,kte,m,j); 
                    n2ovmr(ncol,kte+1,m,j) = n2ovmr(ncol,kte,m,j);

                    //!--- initialize the ozone voume mixing ratio:
                    inirad_sw(o3mmr,plev,kts,kte,m,j);
                    if(o3input == 2){
                        for( int k = 1; k<=noznlevels; k++) o3clim1d(k,m,j) = o3clim(i,k,j);
                        vinterp_ozn_sw(1,ncol,ncol,kte+1,play,pin,noznlevels,o3clim1d,o3mmr,m,j);
                        for( int k = kts; k<=kte+1; k++) o3vmr(ncol,k,m,j) = o3mmr(k,m,j);
                    }else{
                        for(int k = kts; k<=kte+1; k++) o3vmr(ncol,k,m,j) = o3mmr(k,m,j)*amdo;
                    }

                    // //!--- initialize the surface albedo:
                    asdir(ncol,m,j) = albedo(i,j);
                    asdif(ncol,m,j) = albedo(i,j);
                    aldir(ncol,m,j) = albedo(i,j);
                    aldif(ncol,m,j) = albedo(i,j);

                    //!--- set solar constant:
                    real scon = solcon;
                    //!--- set cosine of solar zenith angle:
                    coszen(ncol,m,j) = coszrs;
                    // !--- in MPAS, the solar constant is already provided with eccentricity adjustment, so do not do
                    // !    this here:
                    int dyofyr = 0;
                    real adjes  = 1.0;

                    //--- CALCULATE CLOUD OPTICAL PROPERTIES:
                    if(inflgsw > 0){
                        for(int k = kts; k<=kte; k++){
                            cldfrac(ncol,k,m,j) = cldfra1d(k,m,j);
                        }

                        //--- zero out cloud optical properties here (not used when passing physical properties to rrtmg_sw):
                        for(int k = kts; k<=kte; k++){
                            for(int nb = 1; nb<=nbndsw; nb++){
                                taucld(nb,ncol,k,m,j) = 0.0;
                                ssacld(nb,ncol,k,m,j) = 1.0;
                                asmcld(nb,ncol,k,m,j) = 0.0;
                                fsfcld(nb,ncol,k,m,j) = 0.0;  
                            }
                        }

                        int pcols = ncol;
                        int pver = kte - kts + 1;
                        real gravmks = g;

                        if(has_reqc != 0 ){
                            // !--- fill the effective radius for cloud water with that calculated in the Thompson cloud
                            // !    microphysics parameterization:
                            inflgsw = 3;
                            for(int k = kts;k<=kte;k++){
                                reliq(ncol,k,m,j) = std::max(2.5,re_cloud(i,k,j)*1.e6);
                                if(reliq(ncol,k,m,j) <= 2.5 && cldfrac(ncol,k,m,j) > 0. && xland(i,j)-1.5 > 0.)         reliq(ncol,k,m,j) = 10.5;          //--- ocean
                                else if (reliq(ncol,k,m,j) <= 2.5 && cldfrac(ncol,k,m,j) > 0. && xland(i,j)-1.5 < 0.)    reliq(ncol,k,m,j) = 7.5;           //--- land
                            }
                        }else{
                            // !---  calculate the effective radius for cloud water in layers below the model top. we added the
                            // !     dimensions pcols, pver so that the calls to subroutines relcalc and reicalc remain the same
                            // !     as in WRF although the two dimensions are duplicate and not needed inside the subroutines:
                            landfrac(ncol,m,j) = 2.-xland(i,j);
                            landm(ncol,m,j)    = landfrac(ncol,m,j);
                            snowh(ncol,m,j)    = 0.001*snow(i,j);
                            icefrac(ncol,m,j)  = xice(i,j);

                            relcalc_sw(ncol,pcols,pver,tlay,landfrac,landm,icefrac,reliq,snowh,m,j);
                        }

                        if(has_reqi != 0){
                            // !--- fill the effective radius for cloud ice with that calculated in the Thompson cloud
                            // !    microphysics parameterization:
                            inflgsw  = 4;
                            iceflgsw = 4;
                            for(int k = kts; k<=kte; k++){
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
                            // !--- calculate the effective radius for cloud ice in layers below the model top, following
                            // !    Kristjansson and Mitchell. limit upper bound of reice for Fu ice parameterization and
                            // !    convert effective radius to generalized effective size (*1.0315; Fu, 1996):

                            reicalc_sw(ncol,pcols,pver,tlay,reice,retab,m,j);

                            if(iceflgsw == 3){
                                for(int k = kts; k<=kte; k++){
                                    reice(ncol,k,m,j) = reice(ncol,k,m,j) * 1.0315;
                                    reice(ncol,k,m,j) = std::min(140.0,reice(ncol,k,m,j));
                                }
                            }
                        }

                        if(has_reqs != 0.){
                            // !--- fill the effective radius for snow with that calculated in the Thompson cloud
                            // !    microphysics parameterization:
                            inflgsw  = 5;
                            iceflgsw = 5;
                            for(int k = kts; k<= kte; k++){
                                resnow(ncol,k,m,j) = std::max(10.,re_snow(i,k,j)*1.e6);
                            }
                        }else{
                            for(int k = kts; k<= kte; k++){
                                resnow(ncol,k,m,j) = 10.;
                            }
                        }

                        //--- calculate the cloud liquid water path in units of g/m2 below the model top:
                        for(int k = kts; k<= kte; k++){
                            real gliqwp = qc1d(k,m,j) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;              //!grid box liquid water path.
                            cliqwp(ncol,k,m,j) = gliqwp / std::max(0.01,cldfrac(ncol,k,m,j));                   //!in-cloud liquid water path.
                        }

                        //--- calculate the cloud ice path in units of g/m2 below the model top:
                        if(iceflgsw == 3){
                            for(int k = kts; k<=kte; k++){
                                real gicewp = (qi1d(k,m,j)+qs1d(k,m,j)) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;     //!grid box ice water path.
                                cicewp(ncol,k,m,j) = gicewp / std::max(0.01,cldfrac(ncol,k,m,j));                     //!in-cloud ice water path.
                            }
                        }
                        else if(iceflgsw == 4){
                            for(int k = kts; k<=kte; k++){
                                real gicewp = qi1d(k,m,j) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;               //!grid box ice water path.
                                cicewp(ncol,k,m,j) = gicewp / std::max(0.01,cldfrac(ncol,k,m,j));                     //!in-cloud ice water path.
                            }
                        }

                        //--- calculate the snow path in units of g/m2 below the model top:
                        if(iceflgsw == 5){
                            for(int k = kts; k<=kte; k++){
                                real snow_mass_factor = 1.0;
                                if(resnow(ncol,k,m,j) > 130.){
                                    snow_mass_factor = (130.0/resnow(ncol,k,m,j))*(130.0/resnow(ncol,k,m,j));
                                    resnow(ncol,k,m,j)   = 130.0;
                                }
                                real gsnowp = qs1d(k,m,j) * pdel(ncol,k,m,j)*100.0 / gravmks * 1000.0;                   //!grid box snow path.
                                csnowp(ncol,k,m,j) = snow_mass_factor * gsnowp / std::max(0.01,cldfrac(ncol,k,m,j));
                            }
                        }

                        //--- set the cloud physical property between the model top and the top of the atmosphere.
                        //--- do not allow clouds between the model top and the top of the atmosphere:
                        for(int k = kts; k<=kte; k++){
                            clwpth(ncol,k,m,j) = cliqwp(ncol,k,m,j);
                            ciwpth(ncol,k,m,j) = cicewp(ncol,k,m,j);
                            cswpth(ncol,k,m,j) = csnowp(ncol,k,m,j);
                            rel(ncol,k,m,j) = reliq(ncol,k,m,j);
                            rei(ncol,k,m,j) = reice(ncol,k,m,j);
                            res(ncol,k,m,j) = resnow(ncol,k,m,j); 
                        }
                        //---- no clouds are allowed in the extra layer from model top to TOA:
                        clwpth(ncol,kte+1,m,j)  = 0.;
                        ciwpth(ncol,kte+1,m,j)  = 0.;
                        cswpth(ncol,kte+1,m,j)  = 0.;
                        rel(ncol,kte+1,m,j)     = 10.;
                        rei(ncol,kte+1,m,j)     = 10.;
                        res(ncol,kte+1,m,j)     = 10.;
                        cldfrac(ncol,kte+1,m,j) = 0.;
                        for(int nb = 1; nb<= nbndsw; nb++){
                            taucld(nb,ncol,kte+1,m,j) = 0.;
                            ssacld(nb,ncol,kte+1,m,j) = 1.;
                            asmcld(nb,ncol,kte+1,m,j) = 0.;
                            fsfcld(nb,ncol,kte+1,m,j) = 0.;
                        }                                        
                    }

                    //--- sub-column generator for McICA:
                    int iplon = 1;
                    int irng  = 0;
                    int permuteseed = 1;
                    mcica_subcol_sw_kokkos(iplon   , ncol    , nlay    , icld    , permuteseed , irng    , play    , 
                        cldfrac , ciwpth  , clwpth  , cswpth  , rei         , rel     , res     ,
                        taucld  , ssacld  , asmcld  , fsfcld  , cldfmcl     , ciwpmcl , clwpmcl ,
                        cswpmcl , reicmcl , relqmcl , resnmcl , taucmcl     , ssacmcl , asmcmcl ,
                        fsfcmcl,m,j);
                    // --- initialization of aerosol optical properties:
                    for(int nb = 1; nb<= nbndsw; nb++){
                        for(int k = kts; k<=kte+1; k++){
                            tauaer(ncol,k,nb,m,j) = 0.;
                            ssaaer(ncol,k,nb,m,j) = 1.;
                            asmaer(ncol,k,nb,m,j) = 0.;                                                       
                        }
                    }

                    //!jwfeng from AEMOL USTC: 2020/10/10
                    // !Aerosol optical feedback from MPAS-Chem
                    // !jwfeng modified: 2022/03/21
                    // !Directly use all 14 sw bands from tauaersw, gaersw, waersw                   
                    #ifdef DO_CHEMISTRY
                    string schemes="";
                    if( config_aer_ra_feedback_p && trim_sw(config_chem_opt_p,schemes) != "no_chemistry"){
                        for(int nb = 1; nb<= nbndsw; nb++){
                            wavemid(nb,m,j)=0.5*(wavemin(nb)+wavemax(nb));
                            for(int k = kts; k<=kte; k++){
                                if ( tauaersw(i,k,j,1) > thresh1 && tauaersw(i,k,j,nbndsw) > thresh1 && tauaersw(i,k,j,1) < thresh2 && tauaersw(i,k,j,nbndsw) < thresh2){
                                    tauaer(ncol,k,nb,m,j) = tauaersw(i,k,j,nb);
                                }

                                if(is_init){
                                    //call physics_message('*** tauaer is changed by tauaersw!')
                                    is_init = false;
                                }

                                ssaaer(ncol,k,nb,m,j) = waersw(i,k,j,nb);
                                if ( ssaaer(ncol,k,nb,m,j) < 0.4 ) ssaaer(ncol,k,nb,m,j) = 0.4;
                                if ( ssaaer(ncol,k,nb,m,j) >= 1.0) ssaaer(ncol,k,nb,m,j) = 1.0;

                                asmaer(ncol,k,nb,m,j) = gaersw(i,k,j,nb);
                                if( asmaer(ncol,k,nb,m,j) < 0.5) asmaer(ncol,k,nb,m,j) = 0.5;
                                if( asmaer(ncol,k,nb,m,j) >= 1.0) asmaer(ncol,k,nb,m,j) = 1.0;
                            } 
                        }

                        //Diagnostics
                        for(int nb = 1; nb<= nbndsw; nb++){
                            real slope = 0.0;  //use slope as a sum holder
                    
                            for(int k = kts; k<=kte; k++){
                                slope = slope + tauaer(ncol,k,nb,m,j);        
                            }
                            if ( slope < -1.0e-10 ){
                                // write(msg,'("ERROR: Negative total optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,i,j,nb
                                // call physics_error_fatal(msg)
                            }
                            else if(slope > 5.0 ){
                                // call physics_message("-------------------------")
                                // write(msg,'("WARNING: Large total sw optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,i,j,nb
                                // call physics_message(msg)
                                // call physics_message("Diagnostics 1: k, tauaer(1:14)")

                                // call physics_message("Limit tauaer to tauaer*5.0/sum(tauaer)")
                                for(int k = kts; k<=kte; k++){
                                    // write(msg,'(i4,15f8.2)') k, tauaersw(i,k,j,1:14), tauaer(ncol,k,nb)
                                    // call physics_message(msg)
                                    //czhao set an up-limit here to avoid segmentation fault from extreme AOD
                                    tauaer(ncol,k,nb,m,j) = tauaer(ncol,k,nb,m,j) * 5.0 / slope;           
                                }

                                //  call physics_message("Diagnostics 2: k, gaer(1:14)")
                                for(int k = kts; k<=kte; k++){
                                    // write(msg,'(i4,14f8.2)') k, gaersw(i,k,j,1:14)
                                    // call physics_message(msg)
                                }
                                
                                // call physics_message("Diagnostics 3: k, waer(1:14)")
                                for(int k = kts; k<=kte; k++){
                                    // write(msg,'(i4,14f8.2)') k, waersw(i,k,j,1:14)
                                    // call physics_message(msg)
                                }

                                //call physics_message("Diagnostics 4: k, ssaal, asyal, taual")
                                for(int k = kts; k<=kte; k++){
                                    // write(msg,'(i4,3f8.2)') k, ssaaer(ncol,k,nb), asmaer(ncol,k,nb), tauaer(ncol,k,nb)
                                    // call physics_message(msg)
                                }
                                // call physics_message("-------------------------")
                            }
                        }
                    }//end config_aer_opt
                    #endif
                    //--- CALL TO THE RRTMG SHORT WAVE RADIATION MODEL:
                    rrtmg_sw(ncol      , nlay      , icld      , play      , plev    , tlay     ,
                            tlev      , tsfc      , h2ovmr    , o3vmr     , co2vmr  , ch4vmr   ,
                            n2ovmr    , o2vmr     , asdir     , asdif     , aldir   , aldif    ,
                            coszen    , adjes     , dyofyr    , scon      , inflgsw , iceflgsw ,
                            liqflgsw  , cldfmcl   , taucmcl   , ssacmcl   , asmcmcl , fsfcmcl  ,
                            ciwpmcl   , clwpmcl   , cswpmcl   , reicmcl   , relqmcl , resnmcl  ,
                            tauaer    , ssaaer    , asmaer    , ecaer     , swuflx  , swdflx   ,
                            swhr      , swuflxc   , swdflxc   , swhrc     ,               
                            aer_opt   , 
                            sibvisdir , sibvisdif , sibnirdir , sibnirdif , // !added for ssib coupling.
                            swdkdir   , swdkdif,
                            absa16, ka16, absb16, kb16, forref16, selfref16, sfluxref16, rayl16, layreffr16, strrat16,
                            absa17, ka17, absb17, kb17, forref17, selfref17, sfluxref17, rayl17, layreffr17, strrat17,
                            absa18, ka18, absb18, kb18, forref18, selfref18, sfluxref18, rayl18, layreffr18, strrat18,
                            absa19, ka19, absb19, kb19, forref19, selfref19, sfluxref19, rayl19, layreffr19, strrat19,
                            absa20, ka20, absb20, kb20, forref20, selfref20, sfluxref20, absch420, rayl20, layreffr20,
                            absa21, ka21, absb21, kb21, forref21, selfref21, sfluxref21, rayl21, layreffr21, strrat21,
                            absa22, ka22, absb22, kb22, forref22, selfref22, sfluxref22, rayl22, layreffr22, strrat22,
                            absa23, ka23, forref23, selfref23, sfluxref23, rayl23, layreffr23, givfac23,
                            absa24, ka24, absb24, kb24, forref24, selfref24, sfluxref24, abso3a24, abso3b24, rayla24, raylb24, layreffr24, strrat24,
                            absa25, ka25, sfluxref25, abso3a25, abso3b25, rayl25, layreffr25,
                            sfluxref26, rayl26,
                            absa27, ka27, absb27, kb27, sfluxref27, rayl27, layreffr27, scalekur27,
                            absa28, ka28, absb28, kb28, sfluxref28, rayl28, layreffr28, strrat28,
                            absa29, ka29, absb29, kb29, forref29, selfref29, sfluxref29, absh2o29, absco229, rayl29, layreffr29,
                            cldfmc,ciwpmc,clwpmc,cswpmc,taucmc,taormc,ssacmc,asmcmc,fsfcmc,zcldfmc,ztaucmc,ztaormc,zasycmc,zomgcmc,
                            taua,ssaa,asma,wkl,ztaug,ztaur,zcd,zcu,zfd,zfu,zdbt,zdbt_nodel,zgc,zgcc,zgco,zomc,zomcc,zomco,zrdnd,zrdndc,zref,
                            zrefc,zrefo,zrefd,zrefdc,zrefdo,zrup,zrupd,zrupc,zrupdc,zs1,ztauc,ztauo,ztdn,ztdnd,ztdbt,ztoc,ztor,ztra,ztrac,ztrao,
                            ztrad,ztradc,ztrado,zdbtc,ztdbtc,zincflx,zdbtc_nodel,ztdbt_nodel,ztdbtc_nodel,zsflxzen,m,j);
                            

                    //--- OUTPUTS:
                    gsw(i,j)  = swdflx(1,1,m,j) - swuflx(1,1,m,j);
                    swcf(i,j) = (swdflx(1,kte+2,m,j)-swuflx(1,kte+2,m,j)) - (swdflxc(1,kte+2,m,j)-swuflxc(1,kte+2,m,j));
            
                    if(swupt_p->data != nullptr){
                        //output up and down toa fluxes for total and clear sky:
                        swupt(i,j)  = swuflx(1,kte+2,m,j);
                        swuptc(i,j) = swuflxc(1,kte+2,m,j);
                        swdnt(i,j)  = swdflx(1,kte+2,m,j);
                        swdntc(i,j) = swdflxc(1,kte+2,m,j);
                        // output up and down surface fluxes for total and clear sky:
                        swupb(i,j)  = swuflx(1,1,m,j);
                        swupbc(i,j) = swuflxc(1,1,m,j);
                        swdnb(i,j)  = swdflx(1,1,m,j);
                        swdnbc(i,j) = swdflxc(1,1,m,j);
                    }

                    if(swupflx_p->data != nullptr){
                        for(int k = kts; k<=kte+2; k++){
                            swupflx(i,k,j)  = swuflx(1,k,m,j);
                            swupflxc(i,k,j) = swuflxc(1,k,m,j);
                            swdnflx(i,k,j)  = swdflx(1,k,m,j);
                            swdnflxc(i,k,j) = swdflxc(1,k,m,j);  
                        }
                    }

                    //--- output heating rate tendency; convert heating rate from K/d to K/s:
                    for(int k = kts; k<=kte; k++){
                        tten1d(k,m,j) = swhr(ncol,k,m,j)/86400.;
                        rthratensw(i,k,j) = tten1d(k,m,j)/pi3d(i,k,j);       
                    }                       
                }else{
                    if(swupt_p->data != nullptr){
                        //output up and down toa fluxes for total and clear sky:
                        swupt(i,j)     = 0.;
                        swuptc(i,j)    = 0.;
                        swdnt(i,j)     = 0.;
                        swdntc(i,j)    = 0.;
                        //output up and down surface fluxes for total and clear sky:
                        swupb(i,j)     = 0.;
                        swupbc(i,j)    = 0.;
                        swdnb(i,j)     = 0.;
                        swdnbc(i,j)    = 0.;
                    }                                        
                }
            });
    }
        Kokkos::fence();
        // }
    // }
    
}

KOKKOS_FUNCTION void o3data_sw(const real3d &o3prof,const real4d &plev,int kts,int kte,int i, int j){
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
        prlevh[k] = plev(1,k,i,j);
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

KOKKOS_FUNCTION void inirad_sw(const real3d &o3prof,const real4d &plev,int kts,int kte,int m,int j){
    for(int k = kts; k<= kte+1; k++){
        o3prof(k,m,j) = 0.;
    }
    o3data_sw(o3prof,plev,kts,kte,m,j);
}


KOKKOS_FUNCTION void vinterp_ozn_sw(int lchnk, int ncol, int pcols, int pver, const real4d &pmid_in, const real1d_u &pin_in, int levsiz, const real3d &ozmix_in, const real3d &o3vmr_out,int m,int j){
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
    real pin[SIZE_NLAY+1]; // ozone data level pressures (mks)
    real ozmix[SIZE_NLAY+1]; // ozone mixing ratio
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

KOKKOS_FUNCTION void relcalc_sw(int ncol,int pcols,int pver,const real4d &t, const real3d &landfrac, const real3d &landm, const real3d &icefrac,const real4d &rel, const real3d &snowh,int m,int j){
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
          rel(i,k,m,j) = rliqland + (rliqocean-rliqland) * std::min(1.0,std::max(0.0,(tmelt-t(i,k,m,j))*0.05));
          // Modify for snow depth over land
          rel(i,k,m,j) = rel(i,k,m,j) + (rliqocean-rel(i,k,m,j)) * std::min(1.0,std::max(0.0,snowh(i,m,j)*10.));
          // Ramp between polluted value over land to clean value over ocean.
          rel(i,k,m,j) = rel(i,k,m,j) + (rliqocean-rel(i,k,m,j)) * std::min(1.0,std::max(0.0,1.0-landm(i,m,j)));
          // Ramp between the resultant value and a sea ice value in the presence of ice.
          rel(i,k,m,j) = rel(i,k,m,j) + (rliqice-rel(i,k,m,j)) * std::min(1.0,std::max(0.0,icefrac(i,m,j)));
    // end jrm        
        }
    }
}

KOKKOS_FUNCTION void reicalc_sw(int ncol,int pcols,int pver,const real4d &t, const real4d &re,const real1d_u &retab,int m,int j){
    real corr;
    int index;
    // Tabulated values of re(T) in the temperature interval 
    // 180 K -- 274 K; hexagonal columns assumed:
    for(int k = 1; k <= pver; k++){
        for(int i = 1; i <= ncol; i++){
            index = int(t(i,k,m,j)-179.);
            index = std::min(std::max(index,1),94);
            corr = t(i,k,m,j) - int(t(i,k,m,j));
            re(i,k,m,j) = retab(index)*(1.-corr)+retab(index+1)*corr;
        }
    }
}

KOKKOS_INLINE_FUNCTION int m_sw(int k,int n){
    unsigned int uk = (unsigned int)k;
    if(n < 0){
        return k ^ (uk >> -n);
    }else{
        return k ^ (uk << n);
    }
}

KOKKOS_FUNCTION void kissvec_sw(const int1d_u &seed1,const int1d_u &seed2,const int1d_u &seed3,const int1d_u &seed4,const real1d_u &ran_arr,int size_ran){
    /**     
    ! public domain code
    ! made available from http://www.fortran.com/
    ! downloaded by pjr on 03/16/04 for NCAR CAM
    ! converted to vector form, functions inlined by pjr,mvr on 05/10/2004
    
    ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
    ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
    ! (2) A 3-shift shift-register generator, period 2^32-1,
    ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
    !  Overall period>2^123; 
    */  
    int k,n,kiss;
    for(int i = 1; i <= size_ran; i++){
        seed1(i) = 69069 * seed1(i) + 1327217885;
        seed2(i) = m_sw(m_sw(m_sw(seed2(i), 13), -17), 5);
        seed3(i) = 18000 * (seed3(i) & 65535) +  ((unsigned int)seed3(i) >> 16);
        seed4(i) = 30903 * (seed4(i) & 65535) +  ((unsigned int)seed4(i) >> 16);
        kiss = seed1(i) + seed2(i) + ((unsigned int)seed3(i) << 16) + seed4(i);
        ran_arr(i) = kiss*2.328306e-10 + 0.5;
    } 
}

KOKKOS_FUNCTION void generate_stochastic_clouds_sw (int ncol,int nlay,const int nsubcol,int icld,int irng,const real2d_u& pmid,const real4d& cld,const real4d& clwp,const real4d& ciwp,const real4d& cswp,
    const real5d& tauc,const real5d& ssac,const real5d& asmc,const real5d& fsfc,const real5d& cld_stoch,const real5d& clwp_stoch,const real5d& ciwp_stoch,const real5d& cswp_stoch,
    const real5d& tauc_stoch,const real5d& ssac_stoch,const real5d& asmc_stoch,const real5d& fsfc_stoch,int changeSeed,int m,int j){

    std::vector<real> cldf_spread_s(ncol * nlay * nsubcol);
    real3d_u cldf_spread(cldf_spread_s.data(),{1,1,1},{nsubcol+1,ncol+1,nlay+1});

    real cldf_s[SIZE_NCOL * SIZE_NLAY];
    real2d_u cldf((real *)cldf_s,{1,1},{ncol+1,nlay+1});

    int overlap;                                    // 1 = random overlap, 2 = maximum/random,
                                                     // 3 = maximum overlap, 

    const real cldmin = 1.0e-20;                    //min cloud fraction

    std::vector<real> CDF_s(ncol * nlay * nsubcol);
    real3d_u CDF(CDF_s.data(),{1,1,1},{nsubcol+1,ncol+1,nlay+1});          //random numbers
    int seed1_s[SIZE_NCOL];
    int1d_u seed1((int*)seed1_s,{1},{ncol+1});                //seed to create random number
    int seed2_s[SIZE_NCOL];
    int1d_u seed2((int*)seed2_s,{1},{ncol+1});
    int seed3_s[SIZE_NCOL];
    int1d_u seed3((int*)seed3_s,{1},{ncol+1});
    int seed4_s[SIZE_NCOL];
    int1d_u seed4((int*)seed4_s,{1},{ncol+1});

    real rand_num_s[SIZE_NCOL];
    real1d_u rand_num((real*)rand_num_s,{1},{ncol+1});         //random number (kissvec)
    int iseed;                                              //seed to create random number (Mersenne Twister)

    real rand_num_mt;                                       //random number (Mersenne Twister)

    std::vector<int> isCloudy_s(ncol * nlay * nsubcol);
    int3d_u isCloudy(isCloudy_s.data(),{1,1,1},{nsubcol+1,ncol+1,nlay+1});        // flag that says whether a gridbox is cloudy

    int ilev, isubcol, i, n, ngbm;       //indices

     // ------------------------------------------------------------------------------------------ 

     // Check that irng is in bounds; if not, set to default

     if(irng != 0) irng=1;

     //Pass input cloud overlap setting to local variable
     overlap = icld;

     //Ensure that cloud fractions are in bounds
     for(int ilev = 1;ilev<=nlay;ilev++){
        for(int i = 1;i<=ncol;i++){
            cldf(i,ilev) = cld(i,ilev,m,j);
            if(cldf(i,ilev) < cldmin) cldf(i,ilev) = 0.0;
        }
     }

     // ----- Create seed  --------
     
     // Advance randum number generator by changeseed values
     if(irng == 0){
        for(int i = 1;i<=ncol;i++){
            if(pmid(i,1) < pmid(i,2)) throw std::runtime_error("MCICA_SUBCOL: INVALID ICLD");
            seed1(i) = (pmid(i,1) - int(pmid(i,1)))  * 1000000000;
            seed2(i) = (pmid(i,2) - int(pmid(i,2)))  * 1000000000;
            seed3(i) = (pmid(i,3) - int(pmid(i,3)))  * 1000000000;
            seed4(i) = (pmid(i,4) - int(pmid(i,4)))  * 1000000000;
        }
        for(int i = 1;i<=changeSeed;i++){
            kissvec_sw(seed1, seed2, seed3, seed4, rand_num,ncol);
        }
     }
     else if(irng == 1){
        init_genrand_sw(changeSeed);
     }

     // ------ Apply overlap assumption --------

     // generate the random numbers
     
     switch (overlap){
     case 1:
        if(irng == 0){
            for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                for(int ilev = 1;ilev<=nlay;ilev++){
                    kissvec_sw(seed1, seed2, seed3, seed4, rand_num,ncol);
                    for(int i = 1;i<=ncol;i++) CDF(isubcol,i,ilev) = rand_num(i);
                }
            }
        }
        else if(irng == 1){
            for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                for(int i = 1;i<=ncol;i++){
                    for(int ilev = 1;ilev<=nlay;ilev++){
                        rand_num_mt = genrand_real1_sw();
                        CDF(isubcol,i,ilev) = rand_num_mt;
                    }
                }
            }

        }
     break;
     
     case 2:
        if(irng == 0){
            for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                for(int ilev = 1;ilev<=nlay;ilev++){
                    kissvec_sw(seed1, seed2, seed3, seed4, rand_num,ncol);
                    for(int i = 1;i<=ncol;i++) CDF(isubcol,i,ilev) = rand_num(i);
                }
            }
        }
        else if(irng == 1){
            for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                for(int i = 1;i<=ncol;i++){
                    for(int ilev = 1;ilev<=nlay;ilev++){
                        rand_num_mt = genrand_real1_sw();
                        CDF(isubcol,i,ilev) = rand_num_mt;
                    }
                }
            }

        }

        for(int ilev = 2;ilev<=nlay;ilev++){
            for(int i = 1;i<=ncol;i++){
                for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                    if (CDF(isubcol, i, ilev-1) > 1.0 - cldf(i,ilev-1) ){
                        CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev-1);
                    }
                    else{
                        CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev) * (1. - cldf(i,ilev-1)); 
                    }
                }
            }
        }
     break;
     

     case 3:
        if (irng == 0){
            for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                kissvec_sw(seed1, seed2, seed3, seed4, rand_num,ncol);
                for(int ilev = 1;ilev<=nlay;ilev++){
                    for(int i = 1;i<=ncol;i++) CDF(isubcol,i,ilev) = rand_num(i);
                }
            }
        }        
        else if (irng == 1){
            for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                for(int i = 1;i<=ncol;i++){
                    rand_num_mt = genrand_real1_sw();
                    for(int ilev = 1;ilev<=nlay;ilev++){
                        CDF(isubcol,i,ilev) = rand_num_mt;
                    }
                }
            }
        }                                      
     break;
     }
     
     ngbm = rrtmgsw.ngb[0] - 1;
     for(int ilev = 1;ilev<=nlay;ilev++){
        for(int i = 1;i<=ncol;i++){
            for(int isubcol = 1;isubcol<=nsubcol;isubcol++){
                cldf_spread(isubcol,i,ilev) = cldf(i,ilev);
                isCloudy(isubcol,i,ilev) = CDF(isubcol,i,ilev) >= (1.0 - cldf_spread(isubcol,i,ilev)) ? 1 : 0;
                if ( isCloudy(isubcol,i,ilev) == 1){
                    cld_stoch(isubcol,i,ilev,m,j) = 1.0;
                    clwp_stoch(isubcol,i,ilev,m,j) = clwp(i,ilev,m,j);
                    ciwp_stoch(isubcol,i,ilev,m,j) = ciwp(i,ilev,m,j);
                    cswp_stoch(isubcol,i,ilev,m,j) = cswp(i,ilev,m,j);
                    n = rrtmgsw.ngb[isubcol-1] - ngbm;
                    tauc_stoch(isubcol,i,ilev,m,j) = tauc(n,i,ilev,m,j);
                    ssac_stoch(isubcol,i,ilev,m,j) = ssac(n,i,ilev,m,j);
                    asmc_stoch(isubcol,i,ilev,m,j) = asmc(n,i,ilev,m,j);
                    fsfc_stoch(isubcol,i,ilev,m,j) = fsfc(n,i,ilev,m,j);
                }
                else{
                    cld_stoch(isubcol,i,ilev,m,j) = 0.0;
                    clwp_stoch(isubcol,i,ilev,m ,j) = 0.0;
                    ciwp_stoch(isubcol,i,ilev,m,j) = 0.0;
                    cswp_stoch(isubcol,i,ilev,m,j) = 0.0;
                    tauc_stoch(isubcol,i,ilev,m,j) = 0.0;
                    ssac_stoch(isubcol,i,ilev,m,j) = 1.0;
                    asmc_stoch(isubcol,i,ilev,m,j) = 0.0;
                    fsfc_stoch(isubcol,i,ilev,m,j) = 0.0;
                }
            }
        }
     }
}
KOKKOS_FUNCTION void mcica_subcol_sw_kokkos(int iplon,int ncol,int nlay,int icld,int permuteseed,int irng,const real4d& play,
    const real4d& cldfrac,const real4d& ciwp,const real4d& clwp,const real4d& cswp,const real4d& rei,const real4d& rel,const real4d& res,const real5d& tauc,const real5d& ssac,const real5d& asmc,const real5d& fsfc,
    const real5d& cldfmcl,const real5d& ciwpmcl,const real5d& clwpmcl,const real5d& cswpmcl,const real4d& reicmcl,const real4d& relqmcl,const real4d& resnmcl,
    const real5d& taucmcl,const real5d& ssacmcl,const real5d& asmcmcl,const real5d& fsfcmcl,int m,int j){
    real pmid_s[SIZE_NCOL * SIZE_NLAY];
    real2d_u pmid((real*)pmid_s,{1,1},{ncol+1,nlay+1});

    const int nsubcsw = ngptsw;
    int ilev;

    if(icld == 0) return;
    if(icld < 0 || icld >3){
        throw std::runtime_error("MCICA_SUBCOL: INVALID ICLD");
    }

    for(int col = 1;col<=ncol;col++){
        for(int lay = 1;lay<=nlay;lay++){
            reicmcl(col,lay,m,j) = rei(col,lay,m,j);
            relqmcl(col,lay,m,j) = rel(col,lay,m,j);
            resnmcl(col,lay,m,j) = res(col,lay,m,j);
            pmid(col,lay) = play(col,lay,m,j)*1.e2;
        }
    }
    
    generate_stochastic_clouds_sw (ncol, nlay, nsubcsw, icld, irng, pmid, cldfrac, clwp, ciwp, cswp,
            tauc, ssac, asmc, fsfc, cldfmcl, clwpmcl, ciwpmcl, cswpmcl,
            taucmcl, ssacmcl, asmcmcl, fsfcmcl, permuteseed,m,j);

}

KOKKOS_INLINE_FUNCTION real earth_sun(int idn){
    // !  Purpose: Function to calculate the correction factor of Earth's orbit
    // !  for current day of the year

    // !  idn        : Day of the year
    // !  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

    real gamma;

    gamma = 2.*pi*(idn-1)/365.;

    // ! Use Iqbal's equation 1.2.1
    return 1.000110 + .034221 * cos(gamma) + .001289 * sin(gamma) + .000719 * cos(2.*gamma) + .000077 * sin(2.*gamma);
}

KOKKOS_FUNCTION void inatm_sw (int iplon,int nlay,int icld,int iaer,
    const real4d& play,const real4d& plev,const real4d& tlay,const real4d& tlev,const real3d& tsfc,const real4d& h2ovmr,
    const real4d& o3vmr,const real4d& co2vmr,const real4d& ch4vmr,const real4d& n2ovmr,const real4d& o2vmr,
    real& adjes,int& dyofyr,real& scon,int& inflgsw,int& iceflgsw,int& liqflgsw,
    const real5d& cldfmcl,const real5d& taucmcl,const real5d& ssacmcl,const real5d& asmcmcl,const real5d& fsfcmcl,const real5d& ciwpmcl,const real5d& clwpmcl,const real5d& cswpmcl,
    const real4d& reicmcl,const real4d& relqmcl,const real4d& resnmcl,const real5d& tauaer,const real5d& ssaaer,const real5d& asmaer,
    int& nlayers,const real1d_u& pavel,const real1d_u& pz,const real1d_u& pdp,const real1d_u& tavel,const real1d_u& tz,real tbound,const real1d_u& coldry,const real4d& wkl,
    const real1d_u& adjflux,const real1d_u& solvar,int& inflag,int& iceflag,int& liqflag,const real4d& cldfmc,const real4d& taucmc,
    const real4d& ssacmc,const real4d& asmcmc,const real4d& fsfcmc,const real4d& ciwpmc,const real4d& clwpmc,const real4d& cswpmc,const real1d_u& reicmc,const real1d_u& relqmc,const real1d_u& resnmc,
    const real4d& taua,const real4d& ssaa,const real4d& asma,int m,int j){
    // ! ----- Local -----
    real amd = 28.9660;    //! Effective molecular weight of dry air (g/mol)
    real amw = 18.0160;    //! Molecular weight of water vapor (g/mol)]

    // ! Set molecular weight ratios (for converting mmr to vmr)
    // !  e.g. h2ovmr = h2ommr * amdw)
    real amdw = 1.607793;  //! Molecular weight of dry air / water vapor
    real amdc = 0.658114;  //! Molecular weight of dry air / carbon dioxide
    real amdo = 0.603428;  //! Molecular weight of dry air / ozone
    real amdm = 1.805423;  //! Molecular weight of dry air / methane
    real amdn = 0.658090;  //! Molecular weight of dry air / nitrous oxide
    real amdo2 = 0.905140; //! Molecular weight of dry air / oxygen

    real sbc = 5.67e-08;   //! Stefan-Boltzmann constant (W/m2K4)

    int isp, l, ix, n, imol, ib, ig;        //Loop indices
    real amm, summol;
    real adjflx;                            //flux adjustment for Earth/Sun distance

    nlayers = nlay;

    // !  Initialize all molecular amounts to zero here, then pass input amounts
    // !  into RRTM array WKL below.

    for(int mx = 1; mx<=mxmol; mx++){
        for(int nl = 1; nl<=nlay; nl++){
            wkl(mx,nl,m,j) = 0.0;
        }
    }
    for(int ngp = 1; ngp<=ngptsw; ngp++){
        for(int nl = 1; nl<=nlay; nl++){
            cldfmc(ngp,nl,m,j) = 0.0;
            taucmc(ngp,nl,m,j) = 0.0;
            ssacmc(ngp,nl,m,j) = 1.0;
            asmcmc(ngp,nl,m,j) = 0.0;
            fsfcmc(ngp,nl,m,j) = 0.0;
            ciwpmc(ngp,nl,m,j) = 0.0;
            clwpmc(ngp,nl,m,j) = 0.0;
            cswpmc(ngp,nl,m,j) = 0.0;
        }
    }
    for(int nl = 1; nl<=nlay; nl++){
        reicmc(nl) = 0.0;
        relqmc(nl) = 0.0;
        resnmc(nl) = 0.0;
    }
    for(int nl = 1; nl<=nlay; nl++){
        for(int nbnd = 1;nbnd <= nbndsw; nbnd++){
            taua(nl,nbnd,m,j) = 0.0;
            ssaa(nl,nbnd,m,j) = 1.0;
            asma(nl,nbnd,m,j) = 0.0;
        }
    }

    // ! Set flux adjustment for current Earth/Sun distance (two options).
    // ! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
    adjflx = adjes;

    // ! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
    // !    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU).
    if(dyofyr > 0) adjflx = earth_sun(dyofyr);

    //! Set incoming solar flux adjustment to include adjustment for
    // ! current Earth/Sun distance (ADJFLX) and scaling of default internal
    // ! solar constant (rrsw_scon = 1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
    // ! to a single scaling factor as needed, or to a different value in each 
    // ! band, which may be necessary for paleoclimate simulations. 

    for(ib = jpb1; ib <= jpb2; ib++){
        solvar(ib) = scon / rrsw_scon;
        adjflux(ib) = adjflx * solvar(ib);   
    }
    tbound = tsfc(iplon,m,j);

    // !  Set surface temperature.

    // !  Install input GCM arrays into RRTMG_SW arrays for pressure, temperature,
    // !  and molecular amounts.  
    // !  Pressures are input in mb, or are converted to mb here.
    // !  Molecular amounts are input in volume mixing ratio, or are converted from 
    // !  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
    // !  here. These are then converted to molecular amount (molec/cm2) below.  
    // !  The dry air column COLDRY (in molec/cm2) is calculated from the level 
    // !  pressures, pz (in mb), based on the hydrostatic equation and includes a 
    // !  correction to account for h2o in the layer.  The molecular weight of moist 
    // !  air (amm) is calculated for each layer.  
    // !  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
    // !  assumes GCM input fields are also bottom to top. Input layer indexing
    // !  from GCM fields should be reversed here if necessary.

    pz(0) = plev(iplon,1,m,j);
    tz(0) = tlev(iplon,1,m,j);
    for(l = 1; l <= nlayers; l++){
        pavel(l) = play(iplon,l,m,j);
        tavel(l) = tlay(iplon,l,m,j);
        pz(l) = plev(iplon,l+1,m,j);
        tz(l) = tlev(iplon,l+1,m,j);
        pdp(l) = pz(l-1) - pz(l);
        //! For h2o input in vmr:
        wkl(1,l,m,j) = h2ovmr(iplon,l,m,j);
        wkl(2,l,m,j) = co2vmr(iplon,l,m,j);
        wkl(3,l,m,j) = o3vmr(iplon,l,m,j);
        wkl(4,l,m,j) = n2ovmr(iplon,l,m,j);
        wkl(6,l,m,j) = ch4vmr(iplon,l,m,j);
        wkl(7,l,m,j) = o2vmr(iplon,l,m,j);
        amm = (1. - wkl(1,l,m,j)) * amd + wkl(1,l,m,j) * amw;         
        coldry(l) = (pz(l-1)-pz(l)) * 1.e3 * avogad / (1.e2 * grav * amm * (1. + wkl(1,l,m,j))); 
    }

    // ! At this point all molecular amounts in wkl are in volume mixing ratio; 
    // ! convert to molec/cm2 based on coldry for use in rrtm. 

    for(l = 1; l <= nlayers; l++){
        for( imol = 1; imol <= nmol; imol++){
            wkl(imol,l,m,j) = coldry(l) * wkl(imol,l,m,j);
        }
    }

    // ! Transfer aerosol optical properties to RRTM variables;
    // ! modify to reverse layer indexing here if necessary.

    if (iaer >= 1){
        for(l = 1; l <= nlayers; l++){
            for(ib = 1;ib <= nbndsw; ib++){
                taua(l,ib,m,j) = tauaer(iplon,l,ib,m,j);
                ssaa(l,ib,m,j) = ssaaer(iplon,l,ib,m,j);
                asma(l,ib,m,j) = asmaer(iplon,l,ib,m,j);
            }
        }

    }

    // ! Transfer cloud fraction and cloud optical properties to RRTM variables;
    // ! modify to reverse layer indexing here if necessary.
    if (icld >= 1){
        inflag = inflgsw;
        iceflag = iceflgsw;
        liqflag = liqflgsw;

        // ! Move incoming GCM cloud arrays to RRTMG cloud arrays.
        // ! For GCM input, incoming reicmcl is defined based on selected ice parameterization (inflgsw)

        for(l = 1; l <= nlayers; l++){
            for(ig = 1; ig <= ngptsw; ig++){
                cldfmc(ig,l,m,j) = cldfmcl(ig,iplon,l,m,j);
                taucmc(ig,l,m,j) = taucmcl(ig,iplon,l,m,j);
                ssacmc(ig,l,m,j) = ssacmcl(ig,iplon,l,m,j);
                asmcmc(ig,l,m,j) = asmcmcl(ig,iplon,l,m,j);
                fsfcmc(ig,l,m,j) = fsfcmcl(ig,iplon,l,m,j);
                ciwpmc(ig,l,m,j) = ciwpmcl(ig,iplon,l,m,j);
                clwpmc(ig,l,m,j) = clwpmcl(ig,iplon,l,m,j);
                if(iceflag == 5) cswpmc(ig,l,m,j)=cswpmcl(ig,iplon,l,m,j);     
            }
            reicmc(l) = reicmcl(iplon,l,m,j);
            relqmc(l) = relqmcl(iplon,l,m,j);
            if(iceflag == 5) resnmc(l) = resnmcl(iplon,l,m,j);
        }
    }
}

KOKKOS_FUNCTION void cldprmc_sw(int nlayers,int& inflag,int& iceflag,int& liqflag,const real4d& cldfmc,const real4d& ciwpmc,const real4d& clwpmc,const real4d& cswpmc,const real1d_u& reicmc,
    const real1d_u& relqmc,const real1d_u& resnmc,const real4d& taormc,const real4d& taucmc,const real4d& ssacmc,const real4d& asmcmc,const real4d& fsfcmc,int m,int j){

    //! ------- Local -------
    int ib, lay, istr, index, icx, ig;
    real const eps = 1.e-06;
    real const cldmin = 1.e-20;
    real cwp;                       //! total cloud water path
    real radliq;                    //! cloud liquid droplet radius (microns)
    real radice;                    //! cloud ice effective size (microns)
    real radsno;                    //! cloud snow effective size (microns)
    real factor, fint;

    real taucldorig_a, taucloud_a, ssacloud_a, ffp, ffp1, ffpssa;
    real tauiceorig, scatice, ssaice, tauice, tauliqorig, scatliq, ssaliq, tauliq;
    real tausnoorig, scatsno, ssasno, tausno;
    real fdelta_s[ngptsw + 1];
    real extcoice_s[ngptsw + 1];
    real gice_s[ngptsw + 1];
    real ssacoice_s[ngptsw + 1];
    real forwice_s[ngptsw + 1];
    real extcoliq_s[ngptsw + 1];
    real gliq_s[ngptsw + 1];
    real ssacoliq_s[ngptsw + 1];
    real forwliq_s[ngptsw + 1];
    real extcosno_s[ngptsw + 1];
    real gsno_s[ngptsw + 1];
    real ssacosno_s[ngptsw + 1];
    real forwsno_s[ngptsw + 1];
    real1d_u fdelta((real*)fdelta_s,{1},{ngptsw+1});
    real1d_u extcoice((real*)extcoice_s,{1},{ngptsw+1});
    real1d_u gice((real*)gice_s,{1},{ngptsw+1});
    real1d_u ssacoice((real*)ssacoice_s,{1},{ngptsw+1});
    real1d_u forwice((real*)forwice_s,{1},{ngptsw+1});
    real1d_u gliq((real*)gliq_s,{1},{ngptsw+1});
    real1d_u ssacoliq((real*)ssacoliq_s,{1},{ngptsw+1});
    real1d_u forwliq((real*)forwliq_s,{1},{ngptsw+1});
    real1d_u extcoliq((real*)extcoliq_s,{1},{ngptsw+1});
    real1d_u extcosno((real*)extcosno_s,{1},{ngptsw+1});
    real1d_u gsno((real*)gsno_s,{1},{ngptsw+1});
    real1d_u ssacosno((real*)ssacosno_s,{1},{ngptsw+1});
    real1d_u forwsno((real*)forwsno_s,{1},{ngptsw+1});

    // ! Initialize

    // !jm not thread safe      hvrclc = '$Revision: 1.3 $'

    // ! Some of these initializations are done elsewhere
    for(lay = 1; lay <= nlayers; lay++){
        for(ig = 1; ig<= ngptsw; ig++){
            taormc(ig,lay,m,j) = taucmc(ig,lay,m,j);
        }
    }

    //! Main layer loop
    for(lay = 1; lay <= nlayers; lay++){
        //! Main g-point interval loop  
        for(ig = 1; ig<= ngptsw; ig++){
            cwp = ciwpmc(ig,lay,m,j) + clwpmc(ig,lay,m,j) + cswpmc(ig,lay,m,j);

            if (cldfmc(ig,lay,m,j) >= cldmin && (cwp >= cldmin || taucmc(ig,lay,m,j) >= cldmin)){
                if (inflag == 0){
                    // ! Cloud optical properties already defined in taucmc, ssacmc, asmcmc are unscaled;
                    // ! Apply delta-M scaling here (using Henyey-Greenstein approximation)
                    taucldorig_a = taucmc(ig,lay,m,j);
                    ffp = fsfcmc(ig,lay,m,j);
                    ffp1 = 1.0 - ffp;
                    ffpssa = 1.0 - ffp * ssacmc(ig,lay,m,j);
                    ssacloud_a = ffp1 * ssacmc(ig,lay,m,j) / ffpssa;
                    taucloud_a = ffpssa * taucldorig_a;

                    taormc(ig,lay,m,j) = taucldorig_a;
                    ssacmc(ig,lay,m,j) = ssacloud_a;
                    taucmc(ig,lay,m,j) = taucloud_a;
                    asmcmc(ig,lay,m,j) = (asmcmc(ig,lay,m,j) - ffp) / (ffp1);
                }else if(inflag == 1) {
                    throw std::runtime_error("INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA");//stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'
                    //! (inflag=2): Separate treatement of ice clouds and water clouds.
                }else if(inflag >= 2){
                    radice = reicmc(lay);
                    //! Calculation of absorption coefficients due to ice clouds.
                    if ((ciwpmc(ig,lay,m,j)+cswpmc(ig,lay,m,j)) == 0.0){
                        extcoice(ig) = 0.0;
                        ssacoice(ig) = 0.0;
                        gice(ig)     = 0.0;
                        forwice(ig)  = 0.0;

                        extcosno(ig) = 0.0;
                        ssacosno(ig) = 0.0;
                        gsno(ig)     = 0.0;
                        forwsno(ig)  = 0.0; 
                    }else if(iceflag == 1){
                        // ! (iceflag = 1): 
                        // ! Note: This option uses Ebert and Curry approach for all particle sizes similar to
                        // ! CAM3 implementation, though this is somewhat unjustified for large ice particles
                        ib = rrtmgsw.ngb[ig-1];
                        if (rrtmgsw.wavenum2[ib-16] > 1.43e04) icx = 1;
                        else if (rrtmgsw.wavenum2[ib-16] > 7.7e03) icx =2;
                        else if (rrtmgsw.wavenum2[ib-16] > 5.3e03) icx =3;
                        else if (rrtmgsw.wavenum2[ib-16] > 4.0e03) icx =4;
                        else if (rrtmgsw.wavenum2[ib-16] >= 2.5e03) icx =5;

                        extcoice(ig) = (rrtmgsw.abari[icx-1] + rrtmgsw.bbari[icx-1]/radice);
                        ssacoice(ig) = 1. - rrtmgsw.cbari[icx-1] - rrtmgsw.dbari[icx-1] * radice;
                        gice(ig) = rrtmgsw.ebari[icx-1] + rrtmgsw.fbari[icx-1] * radice;

                        //! Check to ensure upper limit of gice is within physical limits for large particles
                        if (gice(ig)>=1.) gice(ig) = 1. - eps;
                        forwice(ig) = gice(ig)*gice(ig);

                        // ! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ig) < 0.0) throw std::runtime_error ("ICE EXTINCTION LESS THAN 0.0");
                        if (ssacoice(ig) > 1.0) throw std::runtime_error ("ICE SSA GRTR THAN 1.0");
                        if (ssacoice(ig) < 0.0) throw std::runtime_error ("ICE SSA LESS THAN 0.0");
                        if (gice(ig) > 1.0) throw std::runtime_error ("ICE ASYM GRTR THAN 1.0");
                        if (gice(ig) < 0.0) throw std::runtime_error ("ICE ASYM LESS THAN 0.0");
                    }else if(iceflag == 2){
                        // ! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns
                        if (radice < 5.0 || radice > 131.0) throw std::runtime_error ("ICE RADIUS OUT OF BOUNDS");
                        factor = (radice - 2.)/3.;
                        index = int(factor);
                        if (index == 43) index = 42;
                        fint = factor - float(index);
                        ib = rrtmgsw.ngb[ig-1];
                        extcoice(ig) = rrtmgsw.extice2[ib-16][index-1] + fint * (rrtmgsw.extice2[ib-16][index] -  rrtmgsw.extice2[ib-16][index-1]);
                        ssacoice(ig) = rrtmgsw.ssaice2[ib-16][index-1] + fint * (rrtmgsw.ssaice2[ib-16][index] -  rrtmgsw.ssaice2[ib-16][index-1]);
                        gice(ig) = rrtmgsw.asyice2[ib-16][index-1] + fint * (rrtmgsw.asyice2[ib-16][index] -  rrtmgsw.asyice2[ib-16][index-1]);
                        forwice(ig) = gice(ig)*gice(ig);
                        //! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ig) < 0.0) throw std::runtime_error ("ICE EXTINCTION LESS THAN 0.0");
                        if (ssacoice(ig) > 1.0) throw std::runtime_error ("ICE SSA GRTR THAN 1.0");
                        if (ssacoice(ig) < 0.0) throw std::runtime_error ("ICE SSA LESS THAN 0.0");
                        if (gice(ig) > 1.0) throw std::runtime_error ("ICE ASYM GRTR THAN 1.0");
                        if (gice(ig) < 0.0) throw std::runtime_error ("ICE ASYM LESS THAN 0.0");
                    }else if(iceflag >= 3){
                        // ! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns
                        // C++ exceptions cannot be raised safely from a Sunway slave-core
                        // Kokkos kernel.  Keep the lookup index inside the documented
                        // parameterization range instead.
                        if (radice < 5.0) radice = 5.0;
                        if (radice > 140.0) radice = 140.0;
                        factor = (radice - 2.)/3.;
                        index = int(factor);
                        if (index == 46) index = 45;
                        fint = factor - float(index);
                        ib = rrtmgsw.ngb[ig-1];
                        extcoice(ig) = rrtmgsw.extice3[ib-16][index-1] + fint * (rrtmgsw.extice3[ib-16][index] - rrtmgsw.extice3[ib-16][index-1]);
                        ssacoice(ig) = rrtmgsw.ssaice3[ib-16][index-1] + fint * (rrtmgsw.ssaice3[ib-16][index] - rrtmgsw.ssaice3[ib-16][index-1]);
                        gice(ig) = rrtmgsw.asyice3[ib-16][index-1] + fint * (rrtmgsw.asyice3[ib-16][index] - rrtmgsw.asyice3[ib-16][index-1]);
                        fdelta(ig) = rrtmgsw.fdlice3[ib-16][index-1] + fint * (rrtmgsw.fdlice3[ib-16][index] - rrtmgsw.fdlice3[ib-16][index-1]);
                        // Roundoff at the table edge can put fdelta slightly outside
                        // its physical [0, 1] interval.  Clamp it on device rather than
                        // invoking the unsupported slave-core exception runtime.
                        if (fdelta(ig) < 0.0) fdelta(ig) = 0.0;
                        if (fdelta(ig) > 1.0) fdelta(ig) = 1.0;
                        forwice(ig) = fdelta(ig) + 0.5 / ssacoice(ig);
                        //! See Fu 1996 p. 2067 
                        if (forwice(ig) > gice(ig)) forwice(ig) = gice(ig);

                        //! Check to ensure all calculated quantities are within physical limits.  
                        if (extcoice(ig) < 0.0) throw std::runtime_error ("ICE EXTINCTION LESS THAN 0.0");
                        if (ssacoice(ig) > 1.0) throw std::runtime_error ("ICE SSA GRTR THAN 1.0");
                        if (ssacoice(ig) < 0.0) throw std::runtime_error ("ICE SSA LESS THAN 0.0");
                        if (gice(ig) > 1.0) throw std::runtime_error ("ICE ASYM GRTR THAN 1.0");
                        if (gice(ig) < 0.0) throw std::runtime_error ("ICE ASYM LESS THAN 0.0");
                    }

                    // !!!!!!!!!!!!!!!!!! Mukul !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // !!!!  INSERT THE EQUIVALENT SNOW VARIABLE CODE HERE
                    // !!!! Although far from perfect, the snow will utilize the
                    // !!!! same lookup table constants as cloud ice.  Changes
                    // !!!! to those constants for larger particle snow would be
                    // !!!! an improvement.
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    if (cswpmc(ig,lay,m,j)>0.0 && iceflag == 5){
                        radsno = resnmc(lay);
                        if (radsno < 5.0 || radsno > 140.0){
                            // write(errmess,'(A,i5,i5,f8.2,f8.2)' )         &
                            // 'ERROR: SNOW GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'   &
                            // ,ig, lay, cswpmc(ig,lay), radsno
                            // FATAL_ERROR(errmess)
                        }
                        factor = (radsno - 2.)/3.;
                        index = int(factor);
                        if (index == 46) index = 45;
                        fint = factor - float(index);
                        ib = rrtmgsw.ngb[ig-1];
                        extcosno(ig) = rrtmgsw.extice3[ib-16][index-1] + fint * (rrtmgsw.extice3[ib-16][index] - rrtmgsw.extice3[ib-16][index-1]);
                        ssacosno(ig) = rrtmgsw.ssaice3[ib-16][index-1] + fint * (rrtmgsw.ssaice3[ib-16][index] - rrtmgsw.ssaice3[ib-16][index-1]);
                        gsno(ig) = rrtmgsw.asyice3[ib-16][index-1] + fint * (rrtmgsw.asyice3[ib-16][index] - rrtmgsw.asyice3[ib-16][index-1]);
                        fdelta(ig) = rrtmgsw.fdlice3[ib-16][index-1] + fint * (rrtmgsw.fdlice3[ib-16][index] - rrtmgsw.fdlice3[ib-16][index-1]);
                        if (fdelta(ig) < 0.0){
                            // write(errmess, *) 'FDELTA LESS THAN 0.0'
                            // FATAL_ERROR(errmess)
                        }
                        if (fdelta(ig) > 1.0){
                            // write(errmess, *) 'FDELTA GT THAN 1.0'
                            // FATAL_ERROR(errmess)
                        }
                        forwsno(ig) = fdelta(ig) + 0.5 / ssacosno(ig);
                        //! See Fu 1996 p. 2067
                        if (forwsno(ig) > gsno(ig)) forwsno(ig) = gsno(ig);
                        //! Check to ensure all calculated quantities are within physical limits.
                        if (extcosno(ig) < 0.0) throw std::runtime_error ("FATAL ERROR : SNOW EXTINCTION LESS THAN 0.0");
                        if (ssacosno(ig) > 1.0) throw std::runtime_error ("FATAL ERROR : SNOW SSA GRTR THAN 1.0");
                        if (ssacosno(ig) < 0.0) throw std::runtime_error ("FATAL ERROR : SNOW SSA LESS THAN 0.0");
                        if (gsno(ig) > 1.0) throw std::runtime_error ("FATAL ERROR : SNOW ASYM GRTR THAN 1.0");
                        if (gsno(ig) < 0.0) throw std::runtime_error ("FATAL ERROR : SNOW ASYM LESS THAN 0.0");
                    }
                    else{
                        extcosno(ig) = 0.0;
                        ssacosno(ig) = 0.0;
                        gsno(ig)     = 0.0;
                        forwsno(ig)  = 0.0;
                    }

                    // ! Calculation of absorption coefficients due to water clouds.
                    if (clwpmc(ig,lay,m,j) == 0.0){
                        extcoliq(ig) = 0.0;
                        ssacoliq(ig) = 0.0;
                        gliq(ig) = 0.0;
                        forwliq(ig) = 0.0;
                    }
                    else if(liqflag == 1){
                        radliq = relqmc(lay);
                        if (radliq < 1.5 || radliq > 60.) throw std::runtime_error("liquid effective radius out of bounds");
                        index = int(radliq - 1.5);
                        if (index == 0) index = 1;
                        if (index == 58) index = 57;
                        fint = radliq - 1.5 - float(index);
                        ib = rrtmgsw.ngb[ig-1];
                        extcoliq(ig) = rrtmgsw.extliq1[ib-16][index-1] + fint * (rrtmgsw.extliq1[ib-16][index] - rrtmgsw.extliq1[ib-16][index-1]);
                        ssacoliq(ig) = rrtmgsw.ssaliq1[ib-16][index-1] + fint * (rrtmgsw.ssaliq1[ib-16][index] - rrtmgsw.ssaliq1[ib-16][index-1]);
                        if (fint < 0. && ssacoliq(ig) > 1.) ssacoliq(ig) = rrtmgsw.ssaliq1[ib-16][index-1];
                        gliq(ig) = rrtmgsw.asyliq1[ib-16][index-1] + fint * (rrtmgsw.asyliq1[ib-16][index] - rrtmgsw.asyliq1[ib-16][index-1]);
                        forwliq(ig) = gliq(ig)*gliq(ig);
                        //! Check to ensure all calculated quantities are within physical limits.
                        if (extcoliq(ig) < 0.0) throw std::runtime_error ("LIQUID EXTINCTION LESS THAN 0.0");
                        if (ssacoliq(ig) > 1.0) throw std::runtime_error ("LIQUID SSA GRTR THAN 1.0");
                        if (ssacoliq(ig) < 0.0) throw std::runtime_error ("LIQUID SSA LESS THAN 0.0");
                        if (gliq(ig) > 1.0) throw std::runtime_error ("LIQUID ASYM GRTR THAN 1.0");
                        if (gliq(ig) < 0.0) throw std::runtime_error ("LIQUID ASYM LESS THAN 0.0");
                    }

                    if(iceflag < 5){
                        tauliqorig = clwpmc(ig,lay,m,j) * extcoliq(ig);
                        tauiceorig = ciwpmc(ig,lay,m,j) * extcoice(ig);
                        taormc(ig,lay,m,j) = tauliqorig + tauiceorig;
        
                        ssaliq = ssacoliq(ig) * (1. - forwliq(ig)) / (1. - forwliq(ig) * ssacoliq(ig));
                        tauliq = (1. - forwliq(ig) * ssacoliq(ig)) * tauliqorig;
                        ssaice = ssacoice(ig) * (1. - forwice(ig)) / (1. - forwice(ig) * ssacoice(ig));
                        tauice = (1. - forwice(ig) * ssacoice(ig)) * tauiceorig;
                        scatliq = ssaliq * tauliq;
                        scatice = ssaice * tauice;
                        scatsno = 0.0; 
                        taucmc(ig,lay,m,j) = tauliq + tauice;
                    }
                    else{
                        tauliqorig = clwpmc(ig,lay,m,j) * extcoliq(ig);
                        tauiceorig = ciwpmc(ig,lay,m,j) * extcoice(ig);
                        tausnoorig = cswpmc(ig,lay,m,j) * extcosno(ig);
                        taormc(ig,lay,m,j) = tauliqorig + tauiceorig + tausnoorig;

                        ssaliq = ssacoliq(ig) * (1. - forwliq(ig)) / (1. - forwliq(ig) * ssacoliq(ig));
                        tauliq = (1. - forwliq(ig) * ssacoliq(ig)) * tauliqorig;
                        ssaice = ssacoice(ig) * (1. - forwice(ig)) / (1. - forwice(ig) * ssacoice(ig));
                        tauice = (1. - forwice(ig) * ssacoice(ig)) * tauiceorig;
                        ssasno = ssacosno(ig) * (1. - forwsno(ig)) / (1. - forwsno(ig) * ssacosno(ig));
                        tausno = (1. - forwsno(ig) * ssacosno(ig)) * tausnoorig;
                        scatliq = ssaliq * tauliq;
                        scatice = ssaice * tauice;
                        scatsno = ssasno * tausno;
                        taucmc(ig,lay,m,j) = tauliq + tauice + tausno;
                    }

                    // ! Ensure non-zero taucmc and scatice
                    if(taucmc(ig,lay,m,j) == 0.) taucmc(ig,lay,m,j) = cldmin;
                    if(scatice == 0.) scatice = cldmin;
                    if(scatsno == 0.) scatsno = cldmin;

                    if (iceflag < 5) ssacmc(ig,lay,m,j) = (scatliq + scatice) / taucmc(ig,lay,m,j);
                    else ssacmc(ig,lay,m,j) = (scatliq + scatice + scatsno) / taucmc(ig,lay,m,j);

                    if (iceflag == 3 || iceflag == 4){
                        // ! In accordance with the 1996 Fu paper, equation A.3, 
                        // ! the moments for ice were calculated depending on whether using spheres
                        // ! or hexagonal ice crystals.
                        // ! Set asymetry parameter to first moment (istr=1)
                        istr = 1;
                        asmcmc(ig,lay,m,j) = (1.0/(scatliq+scatice))*
                        (scatliq*(std::pow(gliq(ig) , istr) - forwliq(ig)) / 
                        (1.0 - forwliq(ig)) + scatice * 
                        std::pow(((gice(ig)-forwice(ig)) /(1.0 - forwice(ig))) , istr));
                    }
                    else if(iceflag == 5){
                        istr = 1;
                        asmcmc(ig,lay,m,j) = (1.0/(scatliq+scatice+scatsno))
                                        *  (scatliq*(std::pow(gliq(ig) , istr) - forwliq(ig))/(1.0 - forwliq(ig))
                                        + scatice * ((gice(ig)-forwice(ig))/(1.0 - forwice(ig)))
                                        + scatsno * std::pow(((gsno(ig)-forwsno(ig))/(1.0 - forwsno(ig))) , istr));
                    }
                    else{
                        // ! This code is the standard method for delta-m scaling. 
                        // ! Set asymetry parameter to first moment (istr=1)
                        istr = 1;
                        asmcmc(ig,lay,m,j) = (scatliq *(std::pow(gliq(ig),istr) - forwliq(ig)) /
                                (1.0 - forwliq(ig)) + scatice * (std::pow(gice(ig) , istr) - forwice(ig)) /
                                (1.0 - forwice(ig)))/(scatliq + scatice);
                    }
                }//endif inflag == 2
            }
        }//! End g-point interval loop
    }//! End layer loop
}

KOKKOS_FUNCTION void setcoef_sw(int nlayers,const real1d_u& pavel,const real1d_u& tavel,const real1d_u& pz,const real1d_u& tz,real tbound,const real1d_u& coldry,const real4d& wkl,int& laytrop,int& layswtch,int& laylow,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,
	const real1d_u& co2mult,const real1d_u& colch4,real1d_u& colco2,real1d_u& colh2o,real1d_u& colmol,real1d_u& coln2o,
	const real1d_u& colo2, const real1d_u& colo3,const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,
	const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,int m,int j){

	//! ------- Local -------
	int indbound;
	int indlev0;
	int lay;
	int jp1;

	real stpfac;
	real tbndfrac;
	real t0frac;
	real plog;
	real fp;
	real ft;
	real ft1;
	real water;
	real scalefac;
	real factor;
	real co2reg;
	real compfp;

	// Initializations
	stpfac = 296./1013.;

	indbound = tbound - 159.;
	tbndfrac = tbound - int(tbound);
	indlev0  = tz(0) - 159.;
	t0frac   = tz(0) - int(tz(0));

	laytrop  = 0;
	layswtch = 0;
	laylow   = 0;

	//! Begin layer loop
	for(lay = 1; lay <= nlayers; lay++){
        // ! Find the two reference pressures on either side of the
        // ! layer pressure.  Store them in JP and JP1.  Store in FP the
        // ! fraction of the difference (in ln(pressure)) between these
        // ! two values that the layer pressure lies.

        plog = std::log(pavel(lay));
        jp(lay) = int(36. - 5*(plog+0.04));
        if (jp(lay) < 1){
            jp(lay) = 1;

        }
        else if (jp(lay) > 58){
            jp(lay) = 58;
        }
        jp1 = jp(lay) + 1;
        fp = 5. * (rrtmgsw.preflog[jp(lay)-1] - plog);

        // ! Determine, for each reference pressure (JP and JP1), which
        // ! reference temperature (these are different for each  
        // ! reference pressure) is nearest the layer temperature but does
        // ! not exceed it.  Store these indices in JT and JT1, resp.
        // ! Store in FT (resp. FT1) the fraction of the way between JT
        // ! (JT1) and the next highest reference temperature that the 
        // ! layer temperature falls.

        jt(lay) = int(3. + (tavel(lay)-rrtmgsw.tref[jp(lay)-1])/15.);
        if (jt(lay) < 1) {
            jt(lay) = 1;
        }
        else if (jt(lay) > 4){
            jt(lay) = 4;
        }
        ft = ((tavel(lay)-rrtmgsw.tref[jp(lay)-1])/15.) - float(jt(lay)-3);
        jt1(lay) = int(3. + (tavel(lay)-rrtmgsw.tref[jp1-1])/15.);
        if (jt1(lay) < 1){
            jt1(lay) = 1;
        }
        else if (jt1(lay) > 4){
            jt1(lay) = 4;
        }
        ft1 = ((tavel(lay)-rrtmgsw.tref[jp1-1])/15.) - float(jt1(lay)-3);

        water = wkl(1,lay,m,j)/coldry(lay);
        scalefac = pavel(lay) * stpfac / tavel(lay);

        // ! If the pressure is less than ~100mb, perform a different
        // ! set of species interpolations.

        if (plog > 4.56){
            laytrop =  laytrop + 1;
            if (plog >= 6.62) laylow = laylow + 1;

            // ! Set up factors needed to separately include the water vapor
            // ! foreign-continuum in the calculation of absorption coefficient.

            forfac(lay) = scalefac / (1.+water);
            factor = (332.0-tavel(lay))/36.0;
            indfor(lay) = std::min(2, std::max(1, int(factor)));
            forfrac(lay) = factor - float(indfor(lay));

            // ! Set up factors needed to separately include the water vapor
            // ! self-continuum in the calculation of absorption coefficient.

            selffac(lay) = water * forfac(lay);
            factor = (tavel(lay)-188.0)/7.2;
            indself(lay) = std::min(9, std::max(1, int(factor)-7));
            selffrac(lay) = factor - float(indself(lay) + 7);

            // ! Calculate needed column amounts.

            colh2o(lay) = 1.e-20 * wkl(1,lay,m,j);
            colco2(lay) = 1.e-20 * wkl(2,lay,m,j);
            colo3(lay) = 1.e-20 * wkl(3,lay,m,j);
            coln2o(lay) = 1.e-20 * wkl(4,lay,m,j);
            colch4(lay) = 1.e-20 * wkl(6,lay,m,j);
            colo2(lay) = 1.e-20 * wkl(7,lay,m,j);
            colmol(lay) = 1.e-20 * coldry(lay) + colh2o(lay);
            if (colco2(lay) == 0.) colco2(lay) = 1.e-32 * coldry(lay);
            if (coln2o(lay) == 0.) coln2o(lay) = 1.e-32 * coldry(lay);
            if (colch4(lay) == 0.) colch4(lay) = 1.e-32 * coldry(lay);
            if (colo2(lay) == 0.) colo2(lay) = 1.e-32 * coldry(lay);
            // ! Using E = 1334.2 cm-1
            co2reg = 3.55e-24 * coldry(lay);
            co2mult(lay)= (colco2(lay) - co2reg) * 272.63*std::exp(-1919.4/tavel(lay))/(8.7604e-4*tavel(lay));
        }
        else{
            // ! Set up factors needed to separately include the water vapor
            // ! foreign-continuum in the calculation of absorption coefficient.

            forfac(lay) = scalefac / (1.+water);
            factor = (tavel(lay)-188.0)/36.0;
            indfor(lay) = 3;
            forfrac(lay) = factor - 1.0;

            //! Calculate needed column amounts.

            colh2o(lay) = 1.e-20 * wkl(1,lay,m,j);
            colco2(lay) = 1.e-20 * wkl(2,lay,m,j);
            colo3(lay)  = 1.e-20 * wkl(3,lay,m,j);
            coln2o(lay) = 1.e-20 * wkl(4,lay,m,j);
            colch4(lay) = 1.e-20 * wkl(6,lay,m,j);
            colo2(lay)  = 1.e-20 * wkl(7,lay,m,j);
            colmol(lay) = 1.e-20 * coldry(lay) + colh2o(lay);
            if (colco2(lay) == 0.) colco2(lay) = 1.e-32 * coldry(lay);
            if (coln2o(lay) == 0.) coln2o(lay) = 1.e-32 * coldry(lay);
            if (colch4(lay) == 0.) colch4(lay) = 1.e-32 * coldry(lay);
            if (colo2(lay)  == 0.) colo2(lay)  = 1.e-32 * coldry(lay);
            co2reg = 3.55e-24 * coldry(lay);
            co2mult(lay)= (colco2(lay) - co2reg) * 272.63*std::exp(-1919.4/tavel(lay))/(8.7604e-4*tavel(lay));

            selffac(lay) = 0.;
            selffrac(lay)= 0.;
            indself(lay) = 0;
        }

        // ! We have now isolated the layer ln pressure and temperature,
        // ! between two reference pressures and two reference temperatures 
        // ! (for each reference pressure).  We multiply the pressure 
        // ! fraction FP with the appropriate temperature fractions to get 
        // ! the factors that will be needed for the interpolation that yields
        // ! the optical depths (performed in routines TAUGBn for band n).
        compfp = 1. - fp;
        fac10(lay) = compfp * ft;
        fac00(lay) = compfp * (1. - ft);
        fac11(lay) = fp * ft1;
        fac01(lay) = fp * (1. - ft1);
	}//end layer
}

KOKKOS_FUNCTION void reftra_sw(int nlayers,const int1d_u& lrtchk,const real3d& pgg,real prmuz,const real3d& ptau,const real3d& pw,
                const real3d& pref,const real3d& prefd,const real3d& ptra,const real3d& ptrad,int m,int j){
        
    //!Local

    int jk, jl, kmodts;
    int itind;

    real tblind;
    real za, za1, za2;
    real zbeta, zdend, zdenr, zdent;
    real ze1, ze2, zem1, zem2, zemm, zep1, zep2;
    real zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt;
    real zr1, zr2, zr3, zr4, zr5;
    real zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp;
    real zsr3, zt1, zt2, zt3, zt4, zt5, zto1;
    real zw, zwcrit, zwo;
    real denom;

    real eps = 1.e-08;

    //!     ------------------------------------------------------------------
    // ! Initialize

    zsr3=std::sqrt(3.);
    zwcrit=0.9999995;
    kmodts=2;

    for (jk = 1; jk <= nlayers; jk++){
        if(lrtchk(jk) == 0){
            pref(jk,m,j) =0.;
            ptra(jk,m,j) =1.;
            prefd(jk,m,j)=0.;
            ptrad(jk,m,j)=1.;
        }
        else{
            zto1=ptau(jk,m,j);
            zw  =pw(jk,m,j);
            zg  =pgg(jk,m,j);
            // ! General two-stream expressions

            zg3= 3. * zg;
            if(kmodts == 1){
                zgamma1= (7. - zw * (4. + zg3)) * 0.25;
                zgamma2=-(1. - zw * (4. - zg3)) * 0.25;
                zgamma3= (2. - zg3 * prmuz ) * 0.25;
            }
            else if(kmodts == 2){
                zgamma1= (8. - zw * (5. + zg3)) * 0.25;
                zgamma2=  3. *(zw * (1. - zg )) * 0.25;
                zgamma3= (2. - zg3 * prmuz ) * 0.25;
            }
            else if (kmodts == 3){
                zgamma1= zsr3 * (2. - zw * (1. + zg)) * 0.5;
                zgamma2= zsr3 * zw * (1. - zg ) * 0.5;
                zgamma3= (1. - zsr3 * zg * prmuz ) * 0.5;
            }
            zgamma4= 1. - zgamma3;
            // ! Recompute original s.s.a. to test for conservative solution
            // !Balwinder.Singh@pnnl.gov: Code added to avoid 'divide by zero' error in zwo computation

            denom = std::max((1. - (1. - zw) * std::pow((zg /(1. - zg)),2)),1.0E-30);
            zwo= zw / denom;

            if(zwo >= zwcrit){
                // ! Conservative scattering

                za  = zgamma1 * prmuz;
                za1 = za - zgamma3;
                zgt = zgamma1 * zto1;

                // ! Homogeneous reflectance and transmittance,
                // ! collimated beam

                ze1 = std::min ( zto1 / prmuz , 500.);

                // ! Use exponential lookup table for transmittance, or expansion of 
                // ! exponential for low tau
                if(ze1 <= od_lo){
                    ze2 = 1. - ze1 + 0.5 * ze1 * ze1;
                }
                else{
                    tblind = ze1 / (bpade + ze1);
                    itind = tblint * tblind + 0.5;
                    ze2 = rrtmgsw.exp_tbl[itind];
                }
                pref(jk,m,j) = (zgt - za1 * (1. - ze2)) / (1. + zgt);
                ptra(jk,m,j) = 1. - pref(jk,m,j);

                // ! isotropic incidence

                prefd(jk,m,j) = zgt / (1. + zgt);
                ptrad(jk,m,j) = 1. - prefd(jk,m,j) ;

                // ! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
                // ! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
                // ! table returns a transmittance of 1.0.
                if (ze2 == 1.0){
                    pref(jk,m,j) = 0.0;
                    ptra(jk,m,j) = 1.0;
                    prefd(jk,m,j) = 0.0;
                    ptrad(jk,m,j) = 1.0;
                }
            }
            else{
                za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3;
                za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4;
                zrk = std::sqrt ( std::pow(zgamma1,2) - std::pow(zgamma2,2));
                zrp = zrk * prmuz;     
                zrp1 = 1. + zrp;
                zrm1 = 1. - zrp;
                zrk2 = 2. * zrk;
                zrpp = 1. - zrp*zrp;
                zrkg = zrk + zgamma1;
                zr1  = zrm1 * (za2 + zrk * zgamma3);
                zr2  = zrp1 * (za2 - zrk * zgamma3);
                zr3  = zrk2 * (zgamma3 - za2 * prmuz );
                zr4  = zrpp * zrkg;
                zr5  = zrpp * (zrk - zgamma1);
                zt1  = zrp1 * (za1 + zrk * zgamma4);
                zt2  = zrm1 * (za1 - zrk * zgamma4);
                zt3  = zrk2 * (zgamma4 + za1 * prmuz );
                zt4  = zr4;
                zt5  = zr5;

                // ! mji - reformulated code to avoid potential floating point exceptions
                // !               zbeta = - zr5 / zr4
                zbeta = (zgamma1 - zrk) / zrkg;

                // ! Homogeneous reflectance and transmittance
                ze1 = std::min ( zrk * zto1, 500.);
                ze2 = std::min ( zto1 / prmuz , 500.);

                // ! Use exponential lookup table for transmittance, or expansion of 
                // ! exponential for low tau
                if(ze1 <= od_lo){
                    zem1 = 1. - ze1 + 0.5 * ze1 * ze1;
                    zep1 = 1. / zem1;
                }
                else{
                    tblind = ze1 / (bpade + ze1);
                    itind = tblint * tblind + 0.5;
                    zem1 = rrtmgsw.exp_tbl[itind];
                    zep1 = 1. / zem1 ;
                }

                if(ze2 <= od_lo){
                    zem2 = 1. - ze2 + 0.5 * ze2 * ze2;
                    zep2 = 1. / zem2;
                }
                else{
                    tblind = ze2 / (bpade + ze2);
                    itind = tblint * tblind + 0.5;
                    zem2 = rrtmgsw.exp_tbl[itind];
                    zep2 = 1. / zem2;
                }

                zdenr = zr4*zep1 + zr5*zem1;
                zdent = zt4*zep1 + zt5*zem1;
                if (zdenr >= -eps && zdenr <= eps){
                    pref(jk,m,j) = eps;
                    ptra(jk,m,j) = zem2;
                }
                else{
                    pref(jk,m,j) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr;
                    ptra(jk,m,j) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent;
                }

                //! diffuse beam
                zemm = zem1*zem1;
                zdend = 1. / ( (1. - zbeta*zemm ) * zrkg);
                prefd(jk,m,j) =  zgamma2 * (1. - zemm) * zdend;
                ptrad(jk,m,j) =  zrk2*zem1*zdend;
            }
        }
    }
}

KOKKOS_FUNCTION void vrtqdr_sw(int klev,int kw, 
    const real3d& pref,const real3d& prefd,const real3d& ptra,const real3d& ptrad,
    const real3d& pdbt,const real3d& prdnd,const real3d& prup,const real3d& prupd,const real3d& ptdbt,
    const real4d& pfd,const real4d& pfu,int m, int j){
                
    //!Local
    int ikp, ikx, jk;
    real zreflect;
    real ztdn_s[SIZE_NLAY];
    real1d_u ztdn((real*)ztdn_s,{1},{klev+2});

    //! Link lowest layer with surface
    zreflect = 1. / (1. - prefd(klev+1,m,j) * prefd(klev,m,j));
    prup(klev,m,j) = pref(klev,m,j) + (ptrad(klev,m,j) * 
            ((ptra(klev,m,j) - pdbt(klev,m,j)) * prefd(klev+1,m,j) +
            pdbt(klev,m,j) * pref(klev+1,m,j))) * zreflect;
    prupd(klev,m,j) = prefd(klev,m,j) + ptrad(klev,m,j) * ptrad(klev,m,j) *
            prefd(klev+1,m,j) * zreflect;

    //! Pass from bottom to top 
    for(jk = 1; jk <= klev-1; jk++){
        ikp = klev+1-jk;                       
        ikx = ikp-1;
        zreflect = 1. / (1. -prupd(ikp,m,j) * prefd(ikx,m,j));
        prup(ikx,m,j) = pref(ikx,m,j) + (ptrad(ikx,m,j) * 
            ((ptra(ikx,m,j) - pdbt(ikx,m,j)) * prupd(ikp,m,j) +
                pdbt(ikx,m,j) * prup(ikp,m,j))) * zreflect;
        prupd(ikx,m,j) = prefd(ikx,m,j) + ptrad(ikx,m,j) * ptrad(ikx,m,j) *
                prupd(ikp,m,j) * zreflect;
    }

    //! Upper boundary conditions

    ztdn(1) = 1.;
    prdnd(1,m,j) = 0.;
    ztdn(2) = ptra(1,m,j);
    prdnd(2,m,j) = prefd(1,m,j);
    //! Pass from top to bottom

    for(jk = 2; jk <= klev; jk++){
        ikp = jk+1;
        zreflect = 1. / (1. - prefd(jk,m,j) * prdnd(jk,m,j));
        ztdn(ikp) = ptdbt(jk,m,j) * ptra(jk,m,j) +
            (ptrad(jk,m,j) * ((ztdn(jk) - ptdbt(jk,m,j)) +
                ptdbt(jk,m,j) * pref(jk,m,j) * prdnd(jk,m,j))) * zreflect;
        prdnd(ikp,m,j) = prefd(jk,m,j) + ptrad(jk,m,j) * ptrad(jk,m,j) *
                prdnd(jk,m,j) * zreflect;
    }

    //! Up and down-welling fluxes at levels
    for(jk = 1; jk <= klev+1; jk++){
        zreflect = 1. / (1. - prdnd(jk,m,j) * prupd(jk,m,j));
        pfu(jk,kw,m,j) = (ptdbt(jk,m,j) * prup(jk,m,j) +
                (ztdn(jk) - ptdbt(jk,m,j)) * prupd(jk,m,j)) * zreflect;
        pfd(jk,kw,m,j) = ptdbt(jk,m,j) + (ztdn(jk) - ptdbt(jk,m,j)+
                ptdbt(jk,m,j) * prup(jk,m,j) * prdnd(jk,m,j)) * zreflect;
    }
} 

KOKKOS_FUNCTION void taumol16(int nlayers,const real1d_u& colh2o,const real1d_u& colch4,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa16,const real2d_u& absb16,const real1d_u& sfluxref16,const real2d_u& selfref16,const real2d_u& forref16, real strrat16,
     int layreffr16, real rayl16,int m,int j){
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;

    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure), 
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.  

    // ! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            speccomb = colh2o(lay) + strrat16*colch4(lay);
            specparm = colh2o(lay)/speccomb;
            if (specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[0] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[0] + js;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl16;

            for(ig = 1; ig <= rrtmgsw.ng16; ig++){
                    taug(lay,ig,m,j) = speccomb *
                    (fac000 * absa16(ind0   ,ig) +
                    fac100 * absa16(ind0 +1,ig) +
                    fac010 * absa16(ind0 +9,ig) +
                    fac110 * absa16(ind0+10,ig) +
                    fac001 * absa16(ind1   ,ig) +
                    fac101 * absa16(ind1 +1,ig) +
                    fac011 * absa16(ind1 +9,ig) +
                    fac111 * absa16(ind1+10,ig)) +
                    colh2o(lay) *
                    (selffac(lay) * (selfref16(inds,ig) +
                    selffrac(lay) *
                    (selfref16(inds+1,ig) - selfref16(inds,ig))) +
                    forfac(lay) * (forref16(indf,ig) +
                    forfrac(lay) *
                    (forref16(indf+1,ig) - forref16(indf,ig))));
                    taur(lay,ig,m,j) = tauray;
            }

    }

    laysolfr = nlayers;
    
    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            if (jp(lay-1) < layreffr16 && jp(lay) >= layreffr16) laysolfr = lay;
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[0] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[0] + 1;
            tauray = colmol(lay) * rayl16;

            for(ig = 1; ig<=rrtmgsw.ng16; ig++){
                    taug(lay,ig,m,j) = colch4(lay) *
                    (fac00(lay) * absb16(ind0  ,ig) +
                    fac10(lay) * absb16(ind0+1,ig) +
                    fac01(lay) * absb16(ind1  ,ig) +
                    fac11(lay) * absb16(ind1+1,ig)); 
                    if (lay == laysolfr) sfluxzen(ig,m,j) = sfluxref16(ig);
                    taur(lay,ig,m,j) = tauray;
            }
    }

}

KOKKOS_FUNCTION void taumol17(int nlayers,const real1d_u& colh2o,const real1d_u& colco2,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa17,const real4d_u& ka17,const real2d_u& absb17,const real4d_u& kb17,const real2d_u& forref17,const real2d_u& selfref17,
    const real2d_u& sfluxref17,real rayl17,int layreffr17,real strrat17,int m, int j){

    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure), 
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.  

    // ! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            speccomb = colh2o(lay) + strrat17*colco2(lay);
            specparm = colh2o(lay)/speccomb;
            if (specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[1] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[1] + js;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl17;

            for(ig = 1; ig <= rrtmgsw.ng17; ig++){
                    taug(lay,rrtmgsw.ngs16+ig,m,j) = speccomb *
                    (fac000 * absa17(ind0   ,ig) +
                    fac100 * absa17(ind0 +1,ig) +
                    fac010 * absa17(ind0 +9,ig) +
                    fac110 * absa17(ind0+10,ig) +
                    fac001 * absa17(ind1   ,ig) +
                    fac101 * absa17(ind1 +1,ig) +
                    fac011 * absa17(ind1 +9,ig) +
                    fac111 * absa17(ind1+10,ig)) +
                    colh2o(lay) *
                    (selffac(lay) * (selfref17(inds,ig) +
                    selffrac(lay) *
                    (selfref17(inds+1,ig) - selfref17(inds,ig))) +
                    forfac(lay) * (forref17(indf,ig) +
                    forfrac(lay) *
                    (forref17(indf+1,ig) - forref17(indf,ig))));
                    taur(lay,rrtmgsw.ngs16+ig,m,j) = tauray;
            }

    }
    laysolfr = nlayers;
    
    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            if (jp(lay-1) < layreffr17 && jp(lay) >= layreffr17) laysolfr = lay;
            speccomb = colh2o(lay) + strrat17*colco2(lay);
            specparm = colh2o(lay)/speccomb;
            if (specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 4.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[1] + js;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[1] + js;
            indf = indfor(lay);
            tauray = colmol(lay) * rayl17;

            for(ig = 1; ig<= rrtmgsw.ng17; ig++){
                    taug(lay,rrtmgsw.ngs16+ig,m,j) = speccomb *
                    (fac000 * absb17(ind0,ig) + 
                    fac100 * absb17(ind0+1,ig) + 
                    fac010 * absb17(ind0+5,ig) + 
                    fac110 * absb17(ind0+6,ig) + 
                    fac001 * absb17(ind1,ig) + 
                    fac101 * absb17(ind1+1,ig) + 
                    fac011 * absb17(ind1+5,ig) + 
                    fac111 * absb17(ind1+6,ig)) + 
                    colh2o(lay) * 
                    forfac(lay) * (forref17(indf,ig) + 
                    forfrac(lay) * 
                    (forref17(indf+1,ig) - forref17(indf,ig)));
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs16+ig,m,j) = sfluxref17(ig,js) + fs * (sfluxref17(ig,js+1) - sfluxref17(ig,js));
                    taur(lay,rrtmgsw.ngs16+ig,m,j) = tauray;
            }
    }
}

KOKKOS_FUNCTION void taumol18(int nlayers,const real1d_u& colh2o,const real1d_u& colch4,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,const real2d_u& absa18,
    const real4d_u& ka18,const real2d_u& absb18,const real3d_u& kb18,const real2d_u& forref18,const real2d_u& selfref18,const real2d_u& sfluxref18,real rayl18,int layreffr18,real strrat18,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure), 
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    // ! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr18 && jp(lay+1) >= layreffr18) laysolfr = std::min(lay+1,laytrop);
            speccomb = colh2o(lay) + strrat18*colch4(lay);
            specparm = colh2o(lay)/speccomb;
            if(specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[2] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[2] + js;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl18;

            for(ig = 1; ig <= rrtmgsw.ng18; ig++){
                    taug(lay,rrtmgsw.ngs17+ig,m,j) = speccomb *
                    (fac000 * absa18(ind0,ig) +
                    fac100 * absa18(ind0+1,ig) +
                    fac010 * absa18(ind0+9,ig) +
                    fac110 * absa18(ind0+10,ig) +
                    fac001 * absa18(ind1,ig) +
                    fac101 * absa18(ind1+1,ig) +
                    fac011 * absa18(ind1+9,ig) +
                    fac111 * absa18(ind1+10,ig)) +
                    colh2o(lay) *
                    (selffac(lay) * (selfref18(inds,ig) +
                    selffrac(lay) * 
                    (selfref18(inds+1,ig) - selfref18(inds,ig))) + 
                    forfac(lay) * (forref18(indf,ig) + 
                    forfrac(lay) * 
                    (forref18(indf+1,ig) - forref18(indf,ig)))); 
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs17+ig,m,j) = sfluxref18(ig,js)+ fs * (sfluxref18(ig,js+1) - sfluxref18(ig,js));
                    taur(lay,rrtmgsw.ngs17+ig,m,j) = tauray;
            }
    }

        //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[2] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[2] + 1;
            tauray = colmol(lay) * rayl18;

            for(ig = 1; ig<= rrtmgsw.ng18; ig++){
                    taug(lay,rrtmgsw.ngs17+ig,m,j) = colch4(lay) *
                    (fac00(lay) * absb18(ind0  ,ig) +
                    fac10(lay) * absb18(ind0+1,ig) +
                    fac01(lay) * absb18(ind1  ,ig) +
                    fac11(lay) * absb18(ind1+1,ig)); 
                    taur(lay,rrtmgsw.ngs17+ig,m,j) = tauray;
            }
    }       
}

KOKKOS_FUNCTION void taumol19(int nlayers,const real1d_u& colh2o,const real1d_u& colco2,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa19,const real4d_u& ka19,const real2d_u& absb19,const real3d_u& kb19,const real2d_u& forref19,const real2d_u& selfref19,
    const real2d_u& sfluxref19,real rayl19,int layreffr19,real strrat19,int m,int j){

    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure), 
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.  

    laysolfr = laytrop;

    //! Lower atmosphere loop 
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr19 && jp(lay+1) >= layreffr19) laysolfr = std::min(lay+1,laytrop);
            speccomb = colh2o(lay) + strrat19*colco2(lay);
            specparm = colh2o(lay)/speccomb;
            if(specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[3] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[3] + js;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl19;

            for(ig = 1; ig <= rrtmgsw.ng19; ig++){
                    taug(lay,rrtmgsw.ngs18+ig,m,j) = speccomb *
                    (fac000 * absa19(ind0,ig) +
                    fac100 * absa19(ind0+1,ig) +
                    fac010 * absa19(ind0+9,ig) +
                    fac110 * absa19(ind0+10,ig) +
                    fac001 * absa19(ind1,ig) +
                    fac101 * absa19(ind1+1,ig) +
                    fac011 * absa19(ind1+9,ig) +
                    fac111 * absa19(ind1+10,ig)) +
                    colh2o(lay) *
                    (selffac(lay) * (selfref19(inds,ig) +
                    selffrac(lay) * 
                    (selfref19(inds+1,ig) - selfref19(inds,ig))) + 
                    forfac(lay) * (forref19(indf,ig) + 
                    forfrac(lay) * 
                    (forref19(indf+1,ig) - forref19(indf,ig)))); 

            if (lay == laysolfr) sfluxzen(rrtmgsw.ngs18+ig,m,j) = sfluxref19(ig,js)+ fs * (sfluxref19(ig,js+1) - sfluxref19(ig,js));
            taur(lay,rrtmgsw.ngs18+ig,m,j) = tauray;
            }
    }

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[3] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[3] + 1;
            tauray = colmol(lay) * rayl19;

            for(ig = 1; ig<= rrtmgsw.ng19; ig++){
                    taug(lay,rrtmgsw.ngs18+ig,m,j) = colco2(lay) *
                    (fac00(lay) * absb19(ind0  ,ig) +
                    fac10(lay) * absb19(ind0+1,ig) +
                    fac01(lay) * absb19(ind1  ,ig) +
                    fac11(lay) * absb19(ind1+1,ig)); 
                    taur(lay,rrtmgsw.ngs18+ig,m,j) = tauray;
            }
    } 
}

KOKKOS_FUNCTION void taumol20(int nlayers,const real1d_u& colh2o,const real1d_u& colch4,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur, 
    const real2d_u& absa20,const real3d_u& ka20,const real2d_u& absb20,const real3d_u& kb20,const real2d_u& forref20,const real2d_u& selfref20,const real1d_u& sfluxref20,const real1d_u& absch420,
    real rayl20,int layreffr20,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure), 
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.  

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr20 && jp(lay+1) >= layreffr20) laysolfr = std::min(lay+1,laytrop);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[4] + 1;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[4] + 1;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl20;

            for(ig = 1; ig <= rrtmgsw.ng20; ig++){
                    taug(lay,rrtmgsw.ngs19+ig,m,j) = colh2o(lay) * 
                    ((fac00(lay) * absa20(ind0,ig) + 
                            fac10(lay) * absa20(ind0+1,ig) + 
                            fac01(lay) * absa20(ind1,ig) + 
                            fac11(lay) * absa20(ind1+1,ig)) + 
                            selffac(lay) * (selfref20(inds,ig) + 
                            selffrac(lay) * 
                            (selfref20(inds+1,ig) - selfref20(inds,ig))) + 
                            forfac(lay) * (forref20(indf,ig) + 
                            forfrac(lay) * 
                            (forref20(indf+1,ig) - forref20(indf,ig)))) 
                            + colch4(lay) * absch420(ig);

                    taur(lay,rrtmgsw.ngs19+ig,m,j) = tauray; 
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs19+ig,m,j) = sfluxref20(ig) ;
            }
    }

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[4] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[4] + 1;
            indf = indfor(lay);
            tauray = colmol(lay) * rayl20;

            for(ig = 1; ig<= rrtmgsw.ng20; ig++){
                    taug(lay,rrtmgsw.ngs19+ig,m,j) = colh2o(lay) *
                    (fac00(lay) * absb20(ind0  ,ig) +
                    fac10(lay) * absb20(ind0+1,ig) +
                    fac01(lay) * absb20(ind1  ,ig) +
                    fac11(lay) * absb20(ind1+1,ig) +
                    forfac(lay) * (forref20(indf,ig) + 
                    forfrac(lay) * 
                    (forref20(indf+1,ig) - forref20(indf,ig)))) + 
                    colch4(lay) * absch420(ig);
                    taur(lay,rrtmgsw.ngs19+ig,m,j) = tauray;
            }
    }       
}

KOKKOS_FUNCTION void taumol21(int nlayers,const real1d_u& colh2o,const real1d_u& colco2,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur, 
    const real2d_u& absa21,const real4d_u& ka21,const real2d_u& absb21,const real4d_u& kb21,const real2d_u& forref21,const real2d_u& selfref21,const real2d_u& sfluxref21,real rayl21,
    int layreffr21,real strrat21,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure), 
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr21 && jp(lay+1) >= layreffr21) laysolfr = std::min(lay+1,laytrop);
            speccomb = colh2o(lay) + strrat21*colco2(lay);
            specparm = colh2o(lay)/speccomb;
            if(specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[5] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[5] + js;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl21;

            for(ig = 1; ig <= rrtmgsw.ng21; ig++){
                    taug(lay,rrtmgsw.ngs20+ig,m,j) = speccomb *
                    (fac000 * absa21(ind0,ig) +
                    fac100 * absa21(ind0+1,ig) +
                    fac010 * absa21(ind0+9,ig) +
                    fac110 * absa21(ind0+10,ig) +
                    fac001 * absa21(ind1,ig) +
                    fac101 * absa21(ind1+1,ig) +
                    fac011 * absa21(ind1+9,ig) +
                    fac111 * absa21(ind1+10,ig)) +
                    colh2o(lay) *
                    (selffac(lay) * (selfref21(inds,ig) +
                    selffrac(lay) * 
                    (selfref21(inds+1,ig) - selfref21(inds,ig))) + 
                    forfac(lay) * (forref21(indf,ig) + 
                    forfrac(lay) * 
                    (forref21(indf+1,ig) - forref21(indf,ig)))); 
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs20+ig,m,j) = sfluxref21(ig,js) + fs * (sfluxref21(ig,js+1) - sfluxref21(ig,js));
                    taur(lay,rrtmgsw.ngs20+ig,m,j) = tauray;
            }
    }

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            speccomb = colh2o(lay) + strrat21*colco2(lay);
            specparm = colh2o(lay)/speccomb;
            if (specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 4.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[5] + js;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[5] + js;
            indf = indfor(lay);
            tauray = colmol(lay) * rayl21;

            for(ig = 1; ig<= rrtmgsw.ng21; ig++){
                    taug(lay,rrtmgsw.ngs20+ig,m,j) = speccomb *
                    (fac000 * absb21(ind0,ig) + 
                    fac100 * absb21(ind0+1,ig) + 
                    fac010 * absb21(ind0+5,ig) + 
                    fac110 * absb21(ind0+6,ig) + 
                    fac001 * absb21(ind1,ig) + 
                    fac101 * absb21(ind1+1,ig) + 
                    fac011 * absb21(ind1+5,ig) + 
                    fac111 * absb21(ind1+6,ig)) + 
                    colh2o(lay) * 
                    forfac(lay) * (forref21(indf,ig) + 
                    forfrac(lay) * 
                    (forref21(indf+1,ig) - forref21(indf,ig)));
                    taur(lay,rrtmgsw.ngs20+ig,m,j) = tauray;
            }
    }                  
}

KOKKOS_FUNCTION void taumol22(int nlayers,const real1d_u& colh2o,const real1d_u& colo2,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa22,const real4d_u& ka22,const real2d_u& absb22,const real3d_u& kb22,const real2d_u& forref22,const real2d_u& selfref22,const real2d_u& sfluxref22,real rayl22,
    int layreffr22,real strrat22,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;
    real o2adj, o2cont;

    // ! The following factor is the ratio of total O2 band intensity (lines 
    // ! and Mate continuum) to O2 band intensity (line only).  It is needed
    // ! to adjust the optical depths since the k's include only lines.
    o2adj = 1.6;

    // ! Compute the optical depth by interpolating in ln(pressure), 
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr22 && jp(lay+1) >= layreffr22) laysolfr = std::min(lay+1,laytrop);
            o2cont = (4.35e-4)*colo2(lay)/(350.0*2.0);
            speccomb = colh2o(lay) + o2adj*strrat22*colo2(lay);
            specparm = colh2o(lay)/speccomb;
            if(specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);

            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[6] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[6] + js;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl22;

            for(ig = 1; ig <= rrtmgsw.ng22; ig++){
                    taug(lay,rrtmgsw.ngs21+ig,m,j) = speccomb *
                    (fac000 * absa22(ind0,ig) +
                    fac100 * absa22(ind0+1,ig) +
                    fac010 * absa22(ind0+9,ig) +
                    fac110 * absa22(ind0+10,ig) +
                    fac001 * absa22(ind1,ig) +
                    fac101 * absa22(ind1+1,ig) +
                    fac011 * absa22(ind1+9,ig) +
                    fac111 * absa22(ind1+10,ig)) +
                    colh2o(lay) *
                    (selffac(lay) * (selfref22(inds,ig) +
                    selffrac(lay) * 
                    (selfref22(inds+1,ig) - selfref22(inds,ig))) + 
                    forfac(lay) * (forref22(indf,ig) + 
                    forfrac(lay) * 
                    (forref22(indf+1,ig) - forref22(indf,ig))))+ o2cont; 
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs21+ig,m,j) = sfluxref22(ig,js)+ fs * (sfluxref22(ig,js+1) - sfluxref22(ig,js));
                    taur(lay,rrtmgsw.ngs21+ig,m,j) = tauray;
            }
    }
    
    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            o2cont = 4.35e-4*colo2(lay)/(350.0*2.0);
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[6] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[6] + 1;
            tauray = colmol(lay) * rayl22;

            for(ig = 1; ig<= rrtmgsw.ng22; ig++){
                    taug(lay,rrtmgsw.ngs21+ig,m,j) = colo2(lay) * o2adj *
                    (fac00(lay) * absb22(ind0  ,ig) +
                    fac10(lay) * absb22(ind0+1,ig) +
                    fac01(lay) * absb22(ind1  ,ig) +
                    fac11(lay) * absb22(ind1+1,ig)) + o2cont;
                    taur(lay,rrtmgsw.ngs21+ig,m,j) = tauray;
            }
    }              
}

KOKKOS_FUNCTION void taumol23(int nlayers,const real1d_u& colh2o,const real1d_u& colco2,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa23,const real3d_u& ka23,const real2d_u& forref23,const real2d_u& selfref23,const real1d_u& sfluxref23,const real1d_u& rayl23,int layreffr23,real givfac23,int m,int j){

    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure),
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr23 && jp(lay+1) >= layreffr23) laysolfr = std::min(lay+1,laytrop);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[7] + 1;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[7] + 1;
            inds = indself(lay);
            indf = indfor(lay);

            for(ig = 1; ig <= rrtmgsw.ng23; ig++){
                    tauray = colmol(lay) * rayl23(ig);
                    taug(lay,rrtmgsw.ngs22+ig,m,j) = colh2o(lay) * 
                    (givfac23*(fac00(lay) * absa23(ind0,ig) + 
                            fac10(lay) * absa23(ind0+1,ig) + 
                            fac01(lay) * absa23(ind1,ig) + 
                            fac11(lay) * absa23(ind1+1,ig)) + 
                            selffac(lay) * (selfref23(inds,ig) + 
                            selffrac(lay) * 
                            (selfref23(inds+1,ig) - selfref23(inds,ig))) + 
                            forfac(lay) * (forref23(indf,ig) + 
                            forfrac(lay) * 
                            (forref23(indf+1,ig) - forref23(indf,ig))));
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs22+ig,m,j) = sfluxref23(ig);
                    taur(lay,rrtmgsw.ngs22+ig,m,j) = tauray;
            }
    }

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            for(ig = 1; ig<= rrtmgsw.ng23; ig++){
                    taug(lay,rrtmgsw.ngs22+ig,m,j) = 0;
                    taur(lay,rrtmgsw.ngs22+ig,m,j) = colmol(lay) * rayl23(ig);
            }
    }        
}

KOKKOS_FUNCTION void taumol24(int nlayers,const real1d_u& colh2o,const real1d_u& colo2,const real1d_u& colo3,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa24,const real4d_u& ka24,const real2d_u& absb24,const real3d_u& kb24,const real2d_u& forref24,const real2d_u& selfref24,const real2d_u& sfluxref24,const real1d_u& abso3a24,
    const real1d_u& abso3b24,const real2d_u& rayla24,const real1d_u& raylb24,int layreffr24,real strrat24,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure),
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr24 && jp(lay+1) >= layreffr24) laysolfr = std::min(lay+1,laytrop);
            speccomb = colh2o(lay) + strrat24*colo2(lay);
            specparm = colh2o(lay)/speccomb;
            if(specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[8] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[8] + js;
            inds = indself(lay);
            indf = indfor(lay);

            for(ig = 1; ig <= rrtmgsw.ng24; ig++){
                    tauray = colmol(lay) * (rayla24(ig,js) +
                    fs * (rayla24(ig,js+1) - rayla24(ig,js)));
                    taug(lay,rrtmgsw.ngs23+ig,m,j) = speccomb *
                    (fac000 * absa24(ind0,ig) +
                    fac100 * absa24(ind0+1,ig) +
                    fac010 * absa24(ind0+9,ig) +
                    fac110 * absa24(ind0+10,ig) +
                    fac001 * absa24(ind1,ig) +
                    fac101 * absa24(ind1+1,ig) +
                    fac011 * absa24(ind1+9,ig) +
                    fac111 * absa24(ind1+10,ig)) +
                    colo3(lay) * abso3a24(ig) +
                    colh2o(lay) *
                    (selffac(lay) * (selfref24(inds,ig) +
                    selffrac(lay) * 
                    (selfref24(inds+1,ig) - selfref24(inds,ig))) + 
                    forfac(lay) * (forref24(indf,ig) + 
                    forfrac(lay) * 
                    (forref24(indf+1,ig) - forref24(indf,ig)))); 
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs23+ig,m,j) = sfluxref24(ig,js)+ fs * (sfluxref24(ig,js+1) - sfluxref24(ig,js));
                    taur(lay,rrtmgsw.ngs23+ig,m,j) = tauray;
            }
    }

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[8] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[8] + 1;

            for(ig = 1; ig<= rrtmgsw.ng24; ig++){
                    tauray = colmol(lay) * raylb24(ig);
                    taug(lay,rrtmgsw.ngs23+ig,m,j) = colo2(lay) *
                    (fac00(lay) * absb24(ind0  ,ig) +
                    fac10(lay) * absb24(ind0+1,ig) +
                    fac01(lay) * absb24(ind1  ,ig) +
                    fac11(lay) * absb24(ind1+1,ig)) + 
                    colo3(lay) * abso3b24(ig);
                    taur(lay,rrtmgsw.ngs23+ig,m,j) = tauray;
            }
    }    
}

KOKKOS_FUNCTION void taumol25(int nlayers,const real1d_u& colh2o,const real1d_u& colo3,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa25,const real3d_u& ka25,const real1d_u& sfluxref25,const real1d_u& abso3a25,const real1d_u& abso3b25,const real1d_u& rayl25,int layreffr25,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure),
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            if (jp(lay) < layreffr25 && jp(lay+1) >= layreffr25) laysolfr = std::min(lay+1,laytrop);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[9] + 1;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[9] + 1;

            for(ig = 1; ig <= rrtmgsw.ng25; ig++){
                    tauray = colmol(lay) * rayl25(ig);
                    taug(lay,rrtmgsw.ngs24+ig,m,j) = colh2o(lay) * 
                    (fac00(lay) * absa25(ind0,ig) + 
                            fac10(lay) * absa25(ind0+1,ig) + 
                            fac01(lay) * absa25(ind1,ig) + 
                            fac11(lay) * absa25(ind1+1,ig)) + 
                            colo3(lay) * abso3a25(ig);

                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs24+ig,m,j) = sfluxref25(ig);
                    taur(lay,rrtmgsw.ngs24+ig,m,j) = tauray;
            }
    }

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            for(ig = 1; ig<= rrtmgsw.ng25; ig++){
                    tauray = colmol(lay) * rayl25(ig);
                    taug(lay,rrtmgsw.ngs24+ig,m,j) = colo3(lay) * abso3b25(ig);
                    taur(lay,rrtmgsw.ngs24+ig,m,j) = tauray;
            }
    }
}

KOKKOS_FUNCTION void taumol26(int nlayers,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real1d_u& sfluxref26,const real1d_u& rayl26,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure),
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            for(ig = 1; ig <= rrtmgsw.ng26; ig++){
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs25+ig,m,j) = sfluxref26(ig);
                    taug(lay,rrtmgsw.ngs25+ig,m,j) = 0.;
                    taur(lay,rrtmgsw.ngs25+ig,m,j) = colmol(lay) * rayl26(ig); 
            }
    }

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            for(ig = 1; ig<= rrtmgsw.ng26; ig++){
                    taug(lay,rrtmgsw.ngs25+ig,m,j) = 0;
                    taur(lay,rrtmgsw.ngs25+ig,m,j) = colmol(lay) * rayl26(ig);
            }
    }
}

KOKKOS_FUNCTION void taumol27(int nlayers,const real1d_u& colo3,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa27,const real3d_u& ka27,const real2d_u& absb27,const real3d_u& kb27,const real1d_u& sfluxref27,const real1d_u& rayl27,int layreffr27,real scalekur27,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure),
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[11] + 1;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[11] + 1;

            for(ig = 1; ig <= rrtmgsw.ng27; ig++){
                    tauray = colmol(lay) * rayl27(ig);
                    taug(lay,rrtmgsw.ngs26+ig,m,j) = colo3(lay) * 
                    (fac00(lay) * absa27(ind0,ig) + 
                            fac10(lay) * absa27(ind0+1,ig) + 
                            fac01(lay) * absa27(ind1,ig) + 
                            fac11(lay) * absa27(ind1+1,ig));
                    taur(lay,rrtmgsw.ngs26+ig,m,j) = tauray;
            }
    }

    laysolfr = nlayers;

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            if (jp(lay-1) < layreffr27 && jp(lay) >= layreffr27) laysolfr = lay;
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[11] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[11] + 1;

            for(ig = 1; ig<= rrtmgsw.ng27; ig++){
                    tauray = colmol(lay) * rayl27(ig);
                    taug(lay,rrtmgsw.ngs26+ig,m,j) = colo3(lay) *
                    (fac00(lay) * absb27(ind0  ,ig) +
                    fac10(lay) * absb27(ind0+1,ig) +
                    fac01(lay) * absb27(ind1  ,ig) +
                    fac11(lay) * absb27(ind1+1,ig));
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs26+ig,m,j) = scalekur27 * sfluxref27(ig); 
                    taur(lay,rrtmgsw.ngs26+ig,m,j) = tauray;
            }
    }  
}

KOKKOS_FUNCTION void taumol28(int nlayers ,const real1d_u& colo2,const real1d_u& colo3,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa28,const real4d_u& ka28,const real2d_u& absb28,const real4d_u& kb28,const real2d_u& sfluxref28,real rayl28,int layreffr28,real strrat28,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure),
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            speccomb = colo3(lay) + strrat28*colo2(lay);
            specparm = colo3(lay)/speccomb;
            if(specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 8.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[12] + js;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[12] + js;
            tauray = colmol(lay) * rayl28;

            for(ig = 1; ig <= rrtmgsw.ng28; ig++){
                    tauray = colmol(lay) * rayl28;
                    taug(lay,rrtmgsw.ngs27+ig,m,j) = speccomb *
                    (fac000 * absa28(ind0,ig) +
                    fac100 * absa28(ind0+1,ig) +
                    fac010 * absa28(ind0+9,ig) +
                    fac110 * absa28(ind0+10,ig) +
                    fac001 * absa28(ind1,ig) +
                    fac101 * absa28(ind1+1,ig) +
                    fac011 * absa28(ind1+9,ig) +
                    fac111 * absa28(ind1+10,ig));
                    taur(lay,rrtmgsw.ngs27+ig,m,j) = tauray;
            }

    }

    laysolfr = nlayers;

    // ! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            if (jp(lay-1) < layreffr28 && jp(lay) >= layreffr28) laysolfr = lay;
            speccomb = colo3(lay) + strrat28*colo2(lay);
            specparm = colo3(lay)/speccomb ;
            if (specparm >= rrtmgsw.oneminus) specparm = rrtmgsw.oneminus;
            specmult = 4.*(specparm);
            js = 1 + int(specmult);
            fs = std::fmod(specmult, 1. );
            fac000 = (1. - fs) * fac00(lay);
            fac010 = (1. - fs) * fac10(lay);
            fac100 = fs * fac00(lay);
            fac110 = fs * fac10(lay);
            fac001 = (1. - fs) * fac01(lay);
            fac011 = (1. - fs) * fac11(lay);
            fac101 = fs * fac01(lay);
            fac111 = fs * fac11(lay);
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[12] + js;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[12] + js;
            tauray = colmol(lay) * rayl28;

            for(ig = 1; ig<= rrtmgsw.ng28; ig++){
                    taug(lay,rrtmgsw.ngs27+ig,m,j) = speccomb *
                    (fac000 * absb28(ind0  ,ig) +
                    fac100 * absb28(ind0+1,ig) +
                    fac010 * absb28(ind0+5,ig) +
                    fac110 * absb28(ind0+6,ig) +
                    fac001 * absb28(ind1,ig) + 
                    fac101 * absb28(ind1+1,ig) + 
                    fac011 * absb28(ind1+5,ig) + 
                    fac111 * absb28(ind1+6,ig));
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs27+ig,m,j) = sfluxref28(ig,js) + fs * (sfluxref28(ig,js+1) - sfluxref28(ig,js));
                    taur(lay,rrtmgsw.ngs27+ig,m,j) = tauray;
            }
    }
}

KOKKOS_FUNCTION void taumol29(int nlayers,const real1d_u& colh2o,const real1d_u& colco2,int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,const real1d_u& colmol,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,
    const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa29,const real3d_u& ka29,const real2d_u& absb29,const real3d_u& kb29,const real2d_u& forref29,const real2d_u& selfref29,const real1d_u& sfluxref29,const real1d_u& absh2o29,const real1d_u& absco229,real rayl29,int layreffr29,int m,int j){
    
    //!Local
    int ig, ind0, ind1, inds, indf, js, lay, laysolfr;
    real fac000, fac001, fac010, fac011, fac100, fac101;
    real fac110, fac111, fs, speccomb, specmult, specparm, tauray;

    // ! Compute the optical depth by interpolating in ln(pressure),
    // ! temperature, and appropriate species.  Below LAYTROP, the water
    // ! vapor self-continuum is interpolated (in temperature) separately.

    laysolfr = laytrop;

    //! Lower atmosphere loop
    for(lay = 1;lay <= laytrop; lay++){
            ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*rrtmgsw.nspa[13] + 1;
            ind1 = (jp(lay)*5+(jt1(lay)-1))*rrtmgsw.nspa[13] + 1;
            inds = indself(lay);
            indf = indfor(lay);
            tauray = colmol(lay) * rayl29;

            for(ig = 1; ig <= rrtmgsw.ng29; ig++){
                    taug(lay,rrtmgsw.ngs28+ig,m,j) = colh2o(lay) * 
                    ((fac00(lay) * absa29(ind0,ig) + 
                    fac10(lay) * absa29(ind0+1,ig) + 
                    fac01(lay) * absa29(ind1,ig) + 
                    fac11(lay) * absa29(ind1+1,ig)) + 
                    selffac(lay) * (selfref29(inds,ig) + 
                    selffrac(lay) * 
                    (selfref29(inds+1,ig) - selfref29(inds,ig))) + 
                    forfac(lay) * (forref29(indf,ig) + 
                    forfrac(lay) * 
                    (forref29(indf+1,ig) - forref29(indf,ig)))) 
                    + colco2(lay) * absco229(ig) ;
                    taur(lay,rrtmgsw.ngs28+ig,m,j) = tauray;
            }
    }

    laysolfr = nlayers;

    //! Upper atmosphere loop
    for(lay = laytrop+1; lay <= nlayers; lay++){
            if (jp(lay-1) < layreffr29 && jp(lay) >= layreffr29) laysolfr = lay;
            ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*rrtmgsw.nspb[13] + 1;
            ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*rrtmgsw.nspb[13] + 1;
            tauray = colmol(lay) * rayl29;

            for(ig = 1; ig<= rrtmgsw.ng29; ig++){
                    taug(lay,rrtmgsw.ngs28+ig,m,j) = colco2(lay) * 
                    (fac00(lay) * absb29(ind0  ,ig) +
                    fac10(lay) * absb29(ind0+1,ig) +
                    fac01(lay) * absb29(ind1  ,ig) +
                    fac11(lay) * absb29(ind1+1,ig))
                    + colh2o(lay) * absh2o29(ig);
                    if (lay == laysolfr) sfluxzen(rrtmgsw.ngs28+ig,m,j) = sfluxref29(ig) ;
                    taur(lay,rrtmgsw.ngs28+ig,m,j) = tauray;
            }
    }
}

KOKKOS_FUNCTION void  taumol_sw(int nlayers,const real1d_u& colh2o,const real1d_u& colco2,const real1d_u& colch4,const real1d_u& colo2,const real1d_u& colo3,const real1d_u& colmol,
    int laytrop,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,
    const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,
    const real3d& sfluxzen,const real4d& taug,const real4d& taur,
    const real2d_u& absa16,const real4d_u& ka16,const real2d_u& absb16,const real3d_u& kb16,const real2d_u& forref16,const real2d_u& selfref16,const real1d_u& sfluxref16,real rayl16,int layreffr16,real strrat16,
    const real2d_u& absa17,const real4d_u& ka17,const real2d_u& absb17,const real4d_u& kb17,const real2d_u& forref17,const real2d_u& selfref17,const real2d_u& sfluxref17,real rayl17,int layreffr17,real strrat17,
    const real2d_u& absa18,const real4d_u& ka18,const real2d_u& absb18,const real3d_u& kb18,const real2d_u& forref18,const real2d_u& selfref18,const real2d_u& sfluxref18,real rayl18,int layreffr18,real strrat18,
    const real2d_u& absa19,const real4d_u& ka19,const real2d_u& absb19,const real3d_u& kb19,const real2d_u& forref19,const real2d_u& selfref19,const real2d_u& sfluxref19,real rayl19,int layreffr19,real strrat19,
    const real2d_u& absa20,const real3d_u& ka20,const real2d_u& absb20,const real3d_u& kb20,const real2d_u& forref20,const real2d_u& selfref20,const real1d_u& sfluxref20,const real1d_u& absch420,real rayl20,int layreffr20,
    const real2d_u& absa21,const real4d_u& ka21,const real2d_u& absb21,const real4d_u& kb21,const real2d_u& forref21,const real2d_u& selfref21,const real2d_u& sfluxref21,real rayl21,int layreffr21,real strrat21,
    const real2d_u& absa22,const real4d_u& ka22,const real2d_u& absb22,const real3d_u& kb22,const real2d_u& forref22,const real2d_u& selfref22,const real2d_u& sfluxref22,real rayl22,int layreffr22,real strrat22,
    const real2d_u& absa23,const real3d_u& ka23,const real2d_u& forref23,const real2d_u& selfref23,const real1d_u& sfluxref23,const real1d_u& rayl23,int layreffr23,real givfac23,
    const real2d_u& absa24,const real4d_u& ka24,const real2d_u& absb24,const real3d_u& kb24,const real2d_u& forref24,const real2d_u& selfref24,const real2d_u& sfluxref24,const real1d_u& abso3a24,const real1d_u& abso3b24,const real2d_u& rayla24,const real1d_u& raylb24,int layreffr24,real strrat24,
    const real2d_u& absa25,const real3d_u& ka25,const real1d_u& sfluxref25,const real1d_u& abso3a25,const real1d_u& abso3b25,const real1d_u& rayl25,int layreffr25,
    const real1d_u& sfluxref26,const real1d_u& rayl26,
    const real2d_u& absa27,const real3d_u& ka27,const real2d_u& absb27,const real3d_u& kb27,const real1d_u& sfluxref27,const real1d_u& rayl27,int layreffr27,real scalekur27,
    const real2d_u& absa28,const real4d_u& ka28,const real2d_u& absb28,const real4d_u& kb28,const real2d_u& sfluxref28,real rayl28,int layreffr28,real strrat28,
    const real2d_u& absa29,const real3d_u& ka29,const real2d_u& absb29,const real3d_u& kb29,const real2d_u& forref29,const real2d_u& selfref29,const real1d_u& sfluxref29,const real1d_u& absh2o29,const real1d_u& absco229,real rayl29,int layreffr29,int m,int j){

    // ! Calculate gaseous optical depth and planck fractions for each spectral band.
    taumol16(nlayers,colh2o, colch4, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa16, absb16, sfluxref16, selfref16, forref16,strrat16,
        layreffr16,rayl16,m,j);
    
    taumol17( nlayers, colh2o, colco2, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa17, ka17, absb17, kb17, forref17, selfref17,
        sfluxref17, rayl17, layreffr17, strrat17,m,j);

    taumol18( nlayers, colh2o, colch4, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur, absa18,
        ka18, absb18, kb18, forref18, selfref18, sfluxref18, rayl18, layreffr18, strrat18,m,j);

    taumol19( nlayers, colh2o, colco2, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa19, ka19, absb19, kb19, forref19, selfref19,
        sfluxref19, rayl19, layreffr19, strrat19,m,j);
        
    taumol20( nlayers, colh2o, colch4, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur, 
        absa20, ka20, absb20, kb20, forref20, selfref20, sfluxref20, absch420,
        rayl20, layreffr20,m,j);
    
    taumol21(nlayers, colh2o, colco2, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa21, ka21, absb21, kb21, forref21, selfref21,
        sfluxref21, rayl21, layreffr21, strrat21,m,j);
    
    taumol22(nlayers, colh2o, colo2, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa22, ka22, absb22, kb22, forref22, selfref22, sfluxref22, rayl22,
        layreffr22, strrat22,m,j);

    taumol23(nlayers, colh2o, colco2, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa23, ka23, forref23, selfref23, sfluxref23, rayl23,
        layreffr23, givfac23,m,j);

    taumol24(nlayers, colh2o, colo2, colo3, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa24, ka24, absb24, kb24, forref24, selfref24, sfluxref24, abso3a24,
        abso3b24, rayla24, raylb24, layreffr24, strrat24,m,j);

    taumol25(nlayers, colh2o , colo3, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa25, ka25, sfluxref25,abso3a25, abso3b25, rayl25, layreffr25,m,j);

    taumol26(nlayers, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        sfluxref26, rayl26,m,j);
    
    taumol27(nlayers , colo3, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
        absa27, ka27, absb27, kb27, sfluxref27, rayl27, layreffr27, scalekur27,m,j);

    taumol28(nlayers , colo2, colo3, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
    absa28, ka28, absb28, kb28, sfluxref28, rayl28, layreffr28, strrat28,m,j);
    
    taumol29(nlayers, colh2o, colco2, laytrop, jp, jt, jt1, colmol,
        fac00, fac01, fac10, fac11, selffac, selffrac, indself,
        forfac, forfrac, indfor, sfluxzen, taug, taur,
    absa29, ka29, absb29, kb29, forref29, selfref29, sfluxref29, absh2o29, absco229, rayl29, layreffr29,m,j);
}

void spcvmc_sw(int nlayers,int istart,int iend,int icpr,int iout,
    const real1d_u& pavel,const real1d_u& tavel,const real1d_u& pz,const real1d_u& tz,real tbound,const real1d_u& palbd,const real1d_u& palbp,
    const real4d& pcldfmc,const real4d& ptaucmc,const real4d& pasycmc,const real4d& pomgcmc,const real4d& ptaormc,
    const real2d_u& ptaua,const real2d_u& pasya,const real2d_u& pomga,real prmu0,const real1d_u& coldry,const real4d& wkl,const real1d_u& adjflux,
    int& laytrop,int layswtch,int& laylow,const int1d_u& jp,const int1d_u& jt,const int1d_u& jt1,
    const real1d_u& co2mult,const real1d_u& colch4,const real1d_u& colco2,const real1d_u& colh2o,const real1d_u& colmol,const real1d_u& coln2o,const real1d_u& colo2,const real1d_u& colo3,
    const real1d_u& fac00,const real1d_u& fac01,const real1d_u& fac10,const real1d_u& fac11,
    const real1d_u& selffac,const real1d_u& selffrac,const int1d_u& indself,const real1d_u& forfac,const real1d_u& forfrac,const int1d_u& indfor,
    const real1d_u& pbbfd,const real1d_u& pbbfu,const real1d_u& pbbcd,const real1d_u& pbbcu,const real1d_u& puvfd,const real1d_u& puvcd,const real1d_u& pnifd,const real1d_u& pnicd,
    const real1d_u& pbbfddir,const real1d_u& pbbcddir,const real1d_u& puvfddir,const real1d_u& puvcddir,const real1d_u& pnifddir,const real1d_u& pnicddir,
    const real4d& ztaug,const real4d& ztaur,const real4d& zcd,const real4d& zcu,const real4d& zfd,const real4d& zfu,
    const real3d &zdbt,const real3d &zdbt_nodel,const real3d &zgc,const real3d &zgcc,const real3d &zgco,const real3d &zomc,const real3d &zomcc,const real3d &zomco,const real3d &zrdnd,const real3d &zrdndc,const real3d &zref,
    const real3d &zrefc,const real3d &zrefo,const real3d &zrefd,const real3d &zrefdc,const real3d &zrefdo,const real3d &zrup,const real3d &zrupd,const real3d &zrupc,const real3d &zrupdc,const real3d &zs1,const real3d &ztauc,const real3d &ztauo,const real3d &ztdn,const real3d &ztdnd,const real3d &ztdbt,const real3d &ztoc,const real3d &ztor,const real3d &ztra,const real3d &ztrac,const real3d &ztrao,
    const real3d &ztrad,const real3d &ztradc,const real3d &ztrado,const real3d &zdbtc,const real3d &ztdbtc,const real3d &zincflx,const real3d &zdbtc_nodel,const real3d &ztdbt_nodel,const real3d &ztdbtc_nodel,const real3d &zsflxzen,
    const real2d_u& absa16,const real4d_u& ka16,const real2d_u& absb16,const real3d_u& kb16,const real2d_u& forref16,const real2d_u& selfref16,const real1d_u& sfluxref16, real rayl16,int layreffr16, real strrat16,
    const real2d_u& absa17,const real4d_u& ka17,const real2d_u& absb17,const real4d_u& kb17,const real2d_u& forref17,const real2d_u& selfref17,const real2d_u& sfluxref17, real rayl17,int layreffr17, real strrat17,
    const real2d_u& absa18,const real4d_u& ka18,const real2d_u& absb18,const real3d_u& kb18,const real2d_u& forref18,const real2d_u& selfref18,const real2d_u& sfluxref18, real rayl18,int layreffr18, real strrat18,
    const real2d_u& absa19,const real4d_u& ka19,const real2d_u& absb19,const real3d_u& kb19,const real2d_u& forref19,const real2d_u& selfref19,const real2d_u& sfluxref19, real rayl19,int layreffr19, real strrat19,
    const real2d_u& absa20,const real3d_u& ka20,const real2d_u& absb20,const real3d_u& kb20,const real2d_u& forref20,const real2d_u& selfref20,const real1d_u& sfluxref20,const real1d_u& absch420, real rayl20,int layreffr20,
    const real2d_u& absa21,const real4d_u& ka21,const real2d_u& absb21,const real4d_u& kb21,const real2d_u& forref21,const real2d_u& selfref21,const real2d_u& sfluxref21, real rayl21,int layreffr21, real strrat21,
    const real2d_u& absa22,const real4d_u& ka22,const real2d_u& absb22,const real3d_u& kb22,const real2d_u& forref22,const real2d_u& selfref22,const real2d_u& sfluxref22, real rayl22,int layreffr22, real strrat22,
    const real2d_u& absa23,const real3d_u& ka23,const real2d_u& forref23,const real2d_u& selfref23,const real1d_u& sfluxref23,const real1d_u& rayl23,int layreffr23, real givfac23,
    const real2d_u& absa24,const real4d_u& ka24,const real2d_u& absb24,const real3d_u& kb24,const real2d_u& forref24,const real2d_u& selfref24,const real2d_u& sfluxref24,const real1d_u& abso3a24,const real1d_u& abso3b24,const real2d_u& rayla24,const real1d_u& raylb24,int layreffr24, real strrat24,
    const real2d_u& absa25,const real3d_u& ka25,const real1d_u& sfluxref25,const real1d_u& abso3a25,const real1d_u& abso3b25,const real1d_u& rayl25,int layreffr25,
    const real1d_u& sfluxref26,const real1d_u& rayl26,
    const real2d_u& absa27,const real3d_u& ka27,const real2d_u& absb27,const real3d_u& kb27,const real1d_u& sfluxref27,const real1d_u& rayl27,int layreffr27,real scalekur27,
    const real2d_u& absa28,const real4d_u& ka28,const real2d_u& absb28,const real4d_u& kb28,const real2d_u& sfluxref28,real rayl28,int layreffr28,real strrat28,
    const real2d_u& absa29,const real3d_u& ka29,const real2d_u& absb29,const real3d_u& kb29,const real2d_u& forref29,const real2d_u& selfref29,const real1d_u& sfluxref29,const real1d_u& absh2o29,const real1d_u& absco229,
    real rayl29,int layreffr29,int m,int j){

    // ! ------- Local -------
    int lrtchkclr_s[SIZE_NLAY];
    int lrtchkcld_s[SIZE_NLAY];
    int1d_u lrtchkclr((int*)lrtchkclr_s,{1},{nlayers + 1});
    int1d_u lrtchkcld((int*)lrtchkcld_s,{1},{nlayers + 1});
    int klev;
    int ib1, ib2, ibm, igt, ikl, ikp, ikx;
    int iw, jb, jg, jl, jk;
    int itind;

    real tblind, ze1;
    real zclear, zcloud;
    

    real zdbtmc, zdbtmo, zf, zgw, zreflect;
    real zwf, tauorig, repclc;

    
    

    // ! Arrays from rrtmg_sw_vrtqdr routine

    

    //! Initializations
    ib1 = istart;
    ib2 = iend;
    klev = nlayers;
    iw = 0;
    repclc = 1.e-12;

    for(jk = 1; jk<= klev+1; jk++){
        pbbcd(jk)=0.;
        pbbcu(jk)=0.;
        pbbfd(jk)=0.;
        pbbfu(jk)=0.;
        pbbcddir(jk)=0.;
        pbbfddir(jk)=0.;
        puvcd(jk)=0.;
        puvfd(jk)=0.;
        puvcddir(jk)=0.;
        puvfddir(jk)=0.;
        pnicd(jk)=0.;
        pnifd(jk)=0.;
        pnicddir(jk)=0.;
        pnifddir(jk)=0.;
    }

    // ! Calculate the optical depths for gaseous absorption and Rayleigh scattering
    taumol_sw(klev,
        colh2o, colco2, colch4, colo2, colo3, colmol,
        laytrop, jp, jt, jt1,
        fac00, fac01, fac10, fac11,
        selffac, selffrac, indself, forfac, forfrac, indfor,
        zsflxzen, ztaug, ztaur,
        absa16, ka16, absb16, kb16, forref16, selfref16, sfluxref16, rayl16, layreffr16, strrat16,
        absa17, ka17, absb17, kb17, forref17, selfref17, sfluxref17, rayl17, layreffr17, strrat17,
        absa18, ka18, absb18, kb18, forref18, selfref18, sfluxref18, rayl18, layreffr18, strrat18,
        absa19, ka19, absb19, kb19, forref19, selfref19, sfluxref19, rayl19, layreffr19, strrat19,
        absa20, ka20, absb20, kb20, forref20, selfref20, sfluxref20, absch420, rayl20, layreffr20,
        absa21, ka21, absb21, kb21, forref21, selfref21, sfluxref21, rayl21, layreffr21, strrat21,
        absa22, ka22, absb22, kb22, forref22, selfref22, sfluxref22, rayl22, layreffr22, strrat22,
        absa23, ka23, forref23, selfref23, sfluxref23, rayl23, layreffr23, givfac23,
        absa24, ka24, absb24, kb24, forref24, selfref24, sfluxref24, abso3a24, abso3b24, rayla24, raylb24, layreffr24, strrat24,
        absa25, ka25, sfluxref25, abso3a25, abso3b25, rayl25, layreffr25,
        sfluxref26, rayl26,
        absa27, ka27, absb27, kb27, sfluxref27, rayl27, layreffr27, scalekur27,
        absa28, ka28, absb28, kb28, sfluxref28, rayl28, layreffr28, strrat28,
        absa29, ka29, absb29, kb29, forref29, selfref29, sfluxref29, absh2o29, absco229, rayl29, layreffr29,m,j);
    // ! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14
    for(jb = ib1; jb <= ib2; jb++){
        ibm = jb-15;
        igt = rrtmgsw.ngc[ibm-1];

        // ! Reinitialize g-point counter for each band if output for each band is requested.
        if (iout > 0 && ibm >= 2) iw = rrtmgsw.ngs[ibm-2];
        for(jg = 1; jg<= igt; jg++){
            iw = iw + 1;
            zincflx(iw,m,j) = adjflux(jb) * zsflxzen(iw,m,j) * prmu0;
            // ! Clear-sky    
            // !   TOA direct beam    
            ztdbtc(1,m,j)=1.0;
            ztdbtc_nodel(1,m,j)=1.0;
            // !   Surface values
            zdbtc(klev+1,m,j) =0.0;
            ztrac(klev+1,m,j) =0.0;
            ztradc(klev+1,m,j)=0.0;
            zrefc(klev+1,m,j) =palbp(ibm);
            zrefdc(klev+1,m,j)=palbd(ibm);
            zrupc(klev+1,m,j) =palbp(ibm);
            zrupdc(klev+1,m,j)=palbd(ibm);

            // ! Total sky    
            // !   TOA direct beam    
            ztdbt(1,m,j)=1.0;
            ztdbt_nodel(1,m,j)=1.0;
            // !   Surface values
            zdbt(klev+1,m,j) =0.0;
            ztra(klev+1,m,j) =0.0;
            ztrad(klev+1,m,j)=0.0;
            zref(klev+1,m,j) =palbp(ibm);
            zrefd(klev+1,m,j)=palbd(ibm);
            zrup(klev+1,m,j) =palbp(ibm);
            zrupd(klev+1,m,j)=palbd(ibm);
                
            for(jk = 1; jk <= klev; jk++){
                // ! Note: two-stream calculations proceed from top to bottom; 
                // !   RRTMG_SW quantities are given bottom to top and are reversed here

                ikl=klev+1-jk;

                // ! Set logical flag to do REFTRA calculation
                // !   Do REFTRA for all clear layers
                lrtchkclr(jk)=1;

                // !   Do REFTRA only for cloudy layers in profile, since already done for clear layers
                lrtchkcld(jk)=0;
                lrtchkcld(jk)=(pcldfmc(ikl,iw,m,j) > repclc ? 1 : 0);
                
                // ! Clear-sky optical parameters including aerosols
                ztauc(jk,m,j) = ztaur(ikl,iw,m,j) + ztaug(ikl,iw,m,j) + ptaua(ikl,ibm);
                zomcc(jk,m,j) = ztaur(ikl,iw,m,j) * 1.0 + ptaua(ikl,ibm) * pomga(ikl,ibm);
                zgcc(jk,m,j) = pasya(ikl,ibm) * pomga(ikl,ibm) * ptaua(ikl,ibm) / zomcc(jk,m,j);
                zomcc(jk,m,j) = zomcc(jk,m,j) / ztauc(jk,m,j);

                // ! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
                // !   \/\/\/ This block of code is only needed for direct beam calculation

                zclear = 1.0 - pcldfmc(ikl,iw,m,j);
                zcloud = pcldfmc(ikl,iw,m,j);

                // ! Use exponential lookup table for transmittance, or expansion of 
                // ! exponential for low tau
                ze1 = ztauc(jk,m,j) / prmu0;
                if (ze1 <= od_lo)
                    zdbtmc = 1. - ze1 + 0.5 * ze1 * ze1;
                else{
                    tblind = ze1 / (bpade + ze1);
                    itind = tblint * tblind + 0.5;
                    zdbtmc = rrtmgsw.exp_tbl[itind];
                }

                zdbtc_nodel(jk,m,j) = zdbtmc;
                ztdbtc_nodel(jk+1,m,j) = zdbtc_nodel(jk,m,j) * ztdbtc_nodel(jk,m,j);

                // ! Clear + Cloud
                tauorig = ztauc(jk,m,j) + ptaormc(ikl,iw,m,j);

                // ! Use exponential lookup table for transmittance, or expansion of 
                // ! exponential for low tau
                ze1 = tauorig / prmu0;
                if (ze1 <= od_lo)
                    zdbtmo = 1. - ze1 + 0.5 * ze1 * ze1;
                else{
                    tblind = ze1 / (bpade + ze1);
                    itind = tblint * tblind + 0.5;
                    zdbtmo = rrtmgsw.exp_tbl[itind];
                }
                zdbt_nodel(jk,m,j) = zclear*zdbtmc + zcloud*zdbtmo;
                ztdbt_nodel(jk+1,m,j) = zdbt_nodel(jk,m,j) * ztdbt_nodel(jk,m,j);

                // ! Delta scaling - clear   
                zf = zgcc(jk,m,j) * zgcc(jk,m,j);
                zwf = zomcc(jk,m,j) * zf;
                ztauc(jk,m,j) = (1.0 - zwf) * ztauc(jk,m,j);
                zomcc(jk,m,j) = (zomcc(jk,m,j) - zwf) / (1.0 - zwf);
                zgcc (jk,m,j) = (zgcc(jk,m,j) - zf) / (1.0 - zf);

                // ! Total sky optical parameters (cloud properties already delta-scaled)
                // !   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
                if (icpr >= 1){
                    ztauo(jk,m,j) = ztauc(jk,m,j) + ptaucmc(ikl,iw,m,j);
                    zomco(jk,m,j) = ztauc(jk,m,j) * zomcc(jk,m,j) + ptaucmc(ikl,iw,m,j) * pomgcmc(ikl,iw,m,j);
                    zgco (jk,m,j) = (ptaucmc(ikl,iw,m,j) * pomgcmc(ikl,iw,m,j) * pasycmc(ikl,iw,m,j) + ztauc(jk,m,j) * zomcc(jk,m,j) * zgcc(jk,m,j)) / zomco(jk,m,j);
                    zomco(jk,m,j) = zomco(jk,m,j) / ztauo(jk,m,j);
                }
                // ! Total sky optical parameters (if cloud properties not delta scaled)
                // !   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
                else if (icpr == 0){
                    ztauo(jk,m,j) = ztaur(ikl,iw,m,j) + ztaug(ikl,iw,m,j) + ptaua(ikl,ibm) + ptaucmc(ikl,iw,m,j);
                    zomco(jk,m,j) = ptaua(ikl,ibm) * pomga(ikl,ibm) + ptaucmc(ikl,iw,m,j) * pomgcmc(ikl,iw,m,j) + ztaur(ikl,iw,m,j) * 1.0;
                    zgco (jk,m,j) = (ptaucmc(ikl,iw,m,j) * pomgcmc(ikl,iw,m,j) * pasycmc(ikl,iw,m,j) + ptaua(ikl,ibm)*pomga(ikl,ibm)*pasya(ikl,ibm)) / zomco(jk,m,j);
                    zomco(jk,m,j) = zomco(jk,m,j) / ztauo(jk,m,j);
                    // ! Delta scaling - clouds 
                    // !   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
                    zf = zgco(jk,m,j) * zgco(jk,m,j);
                    zwf = zomco(jk,m,j) * zf;
                    ztauo(jk,m,j) = (1. - zwf) * ztauo(jk,m,j);
                    zomco(jk,m,j) = (zomco(jk,m,j) - zwf) / (1.0 - zwf);
                    zgco (jk,m,j) = (zgco(jk,m,j) - zf) / (1.0 - zf);
                }
            }
            // ! Clear sky reflectivities
            reftra_sw (klev,
                lrtchkclr, zgcc, prmu0, ztauc, zomcc,
                zrefc, zrefdc, ztrac, ztradc,m,j);

            //! Total sky reflectivities      
            reftra_sw (klev,
                lrtchkcld, zgco, prmu0, ztauo, zomco,
                zrefo, zrefdo, ztrao, ztrado,m,j);

            for(jk = 1; jk <= klev; jk++){
                // ! Combine clear and cloudy contributions for total sky
                ikl = klev+1-jk; 
                zclear = 1.0 - pcldfmc(ikl,iw,m,j);
                zcloud = pcldfmc(ikl,iw,m,j);

                zref(jk,m,j) = zclear*zrefc(jk,m,j) + zcloud*zrefo(jk,m,j);
                zrefd(jk,m,j)= zclear*zrefdc(jk,m,j) + zcloud*zrefdo(jk,m,j);
                ztra(jk,m,j) = zclear*ztrac(jk,m,j) + zcloud*ztrao(jk,m,j);
                ztrad(jk,m,j)= zclear*ztradc(jk,m,j) + zcloud*ztrado(jk,m,j);

                // ! Use exponential lookup table for transmittance, or expansion of 
                // ! exponential for low tau
                ze1 = ztauc(jk,m,j) / prmu0;
                if (ze1 <= od_lo)
                    zdbtmc = 1. - ze1 + 0.5 * ze1 * ze1;
                else{
                    tblind = ze1 / (bpade + ze1);
                    itind = tblint * tblind + 0.5;
                    zdbtmc = rrtmgsw.exp_tbl[itind];
                }

                zdbtc(jk,m,j) = zdbtmc;
                ztdbtc(jk+1,m,j) = zdbtc(jk,m,j)*ztdbtc(jk,m,j);

                // ! Use exponential lookup table for transmittance, or expansion of 
                // ! exponential for low tau
                ze1 = ztauo(jk,m,j) / prmu0;
                if (ze1 <= od_lo)
                    zdbtmo = 1. - ze1 + 0.5 * ze1 * ze1;
                else{
                    tblind = ze1 / (bpade + ze1);
                    itind = tblint * tblind + 0.5;
                    zdbtmo = rrtmgsw.exp_tbl[itind];
                }

                zdbt(jk,m,j) = zclear*zdbtmc + zcloud*zdbtmo;
                ztdbt(jk+1,m,j) = zdbt(jk,m,j)*ztdbt(jk,m,j);
            }
            // ! Vertical quadrature for clear-sky fluxes                
            vrtqdr_sw(klev, iw,
                    zrefc, zrefdc, ztrac, ztradc,
                    zdbtc, zrdndc, zrupc, zrupdc, ztdbtc,
                    zcd, zcu,m,j);
        
            //! Vertical quadrature for cloudy fluxes
            vrtqdr_sw(klev, iw,
                    zref, zrefd, ztra, ztrad,
                    zdbt, zrdnd, zrup, zrupd, ztdbt,
                    zfd, zfu,m,j);

            // ! Upwelling and downwelling fluxes at levels
            // !   Two-stream calculations go from top to bottom; 
            // !   layer indexing is reversed to go bottom to top for output arrays
            for(jk =1; jk<=klev+1; jk++){
                ikl = klev+2-jk;
                // ! Accumulate spectral fluxes over whole spectrum  
                pbbfu(ikl) = pbbfu(ikl) + zincflx(iw,m,j)*zfu(jk,iw,m,j);
                pbbfd(ikl) = pbbfd(ikl) + zincflx(iw,m,j)*zfd(jk,iw,m,j);
                pbbcu(ikl) = pbbcu(ikl) + zincflx(iw,m,j)*zcu(jk,iw,m,j);
                pbbcd(ikl) = pbbcd(ikl) + zincflx(iw,m,j)*zcd(jk,iw,m,j);
                pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw,m,j)*ztdbt_nodel(jk,m,j);
                pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw,m,j)*ztdbtc_nodel(jk,m,j);
                // ! Accumulate direct fluxes for UV/visible bands
                if (ibm >= 10 && ibm <= 13){
                    puvcd(ikl) = puvcd(ikl) + zincflx(iw,m,j)*zcd(jk,iw,m,j);
                    puvfd(ikl) = puvfd(ikl) + zincflx(iw,m,j)*zfd(jk,iw,m,j);
                    puvcddir(ikl) = puvcddir(ikl) + zincflx(iw,m,j)*ztdbtc_nodel(jk,m,j);
                    puvfddir(ikl) = puvfddir(ikl) + zincflx(iw,m,j)*ztdbt_nodel(jk,m,j);
                }
                // ! Accumulate direct fluxes for near-IR bands
                else if (ibm == 14 || ibm <= 9){
                    pnicd(ikl) = pnicd(ikl) + zincflx(iw,m,j)*zcd(jk,iw,m,j);
                    pnifd(ikl) = pnifd(ikl) + zincflx(iw,m,j)*zfd(jk,iw,m,j);
                    pnicddir(ikl) = pnicddir(ikl) + zincflx(iw,m,j)*ztdbtc_nodel(jk,m,j);
                    pnifddir(ikl) = pnifddir(ikl) + zincflx(iw,m,j)*ztdbt_nodel(jk,m,j);
                }

            }

        }//end loop on jg
    }// end loop on jg
}

KOKKOS_FUNCTION void rrtmg_sw (int ncol, int nlay, int icld,
    const real4d& play, const real4d& plev, const real4d& tlay, const real4d& tlev, const real3d& tsfc, 
    const real4d& h2ovmr, const real4d& o3vmr, const real4d& co2vmr, const real4d& ch4vmr, const real4d&n2ovmr, const real4d& o2vmr,
    const real3d& asdir,const real3d& asdif,const real3d& aldir,const real3d& aldif,
    const real3d& coszen,real& adjes,int& dyofyr,real& scon, 
    int& inflgsw ,int& iceflgsw,int& liqflgsw,const real5d& cldfmcl ,
    const real5d& taucmcl ,const real5d& ssacmcl ,const real5d& asmcmcl ,const real5d& fsfcmcl ,
    const real5d& ciwpmcl ,const real5d& clwpmcl ,const real5d& cswpmcl ,const real4d& reicmcl ,const real4d& relqmcl ,const real4d& resnmcl,
    const real5d& tauaer  ,const real5d& ssaaer  ,const real5d& asmaer  ,const real5d& ecaer   ,
    const real4d& swuflx  ,const real4d& swdflx  ,const real4d& swhr    ,const real4d& swuflxc ,const real4d& swdflxc ,const real4d& swhrc,int& aer_opt,
    const real4d& sibvisdir,const real4d& sibvisdif,const real4d& sibnirdir,const real4d& sibnirdif,
    const real4d& swdkdir,const real4d& swdkdif,
    const real2d_u& absa16,const real4d_u& ka16,const real2d_u& absb16,const real3d_u& kb16,const real2d_u& forref16,const real2d_u& selfref16,const real1d_u& sfluxref16, real rayl16,int layreffr16, real strrat16,
    const real2d_u& absa17,const real4d_u& ka17,const real2d_u& absb17,const real4d_u& kb17,const real2d_u& forref17,const real2d_u& selfref17,const real2d_u& sfluxref17, real rayl17,int layreffr17, real strrat17,
    const real2d_u& absa18,const real4d_u& ka18,const real2d_u& absb18,const real3d_u& kb18,const real2d_u& forref18,const real2d_u& selfref18,const real2d_u& sfluxref18, real rayl18,int layreffr18, real strrat18,
    const real2d_u& absa19,const real4d_u& ka19,const real2d_u& absb19,const real3d_u& kb19,const real2d_u& forref19,const real2d_u& selfref19,const real2d_u& sfluxref19, real rayl19,int layreffr19, real strrat19,
    const real2d_u& absa20,const real3d_u& ka20,const real2d_u& absb20,const real3d_u& kb20,const real2d_u& forref20,const real2d_u& selfref20,const real1d_u& sfluxref20,const real1d_u& absch420, real rayl20,int layreffr20,
    const real2d_u& absa21,const real4d_u& ka21,const real2d_u& absb21,const real4d_u& kb21,const real2d_u& forref21,const real2d_u& selfref21,const real2d_u& sfluxref21, real rayl21,int layreffr21, real strrat21,
    const real2d_u& absa22,const real4d_u& ka22,const real2d_u& absb22,const real3d_u& kb22,const real2d_u& forref22,const real2d_u& selfref22,const real2d_u& sfluxref22, real rayl22,int layreffr22, real strrat22,
    const real2d_u& absa23,const real3d_u& ka23,const real2d_u& forref23,const real2d_u& selfref23,const real1d_u& sfluxref23,const real1d_u& rayl23,int layreffr23, real givfac23,
    const real2d_u& absa24,const real4d_u& ka24,const real2d_u& absb24,const real3d_u& kb24,const real2d_u& forref24,const real2d_u& selfref24,const real2d_u& sfluxref24,const real1d_u& abso3a24,const real1d_u& abso3b24,const real2d_u& rayla24,const real1d_u& raylb24,int layreffr24, real strrat24,
    const real2d_u& absa25,const real3d_u& ka25,const real1d_u& sfluxref25,const real1d_u& abso3a25,const real1d_u& abso3b25,const real1d_u& rayl25,int layreffr25,
    const real1d_u& sfluxref26,const real1d_u& rayl26,
    const real2d_u& absa27,const real3d_u& ka27,const real2d_u& absb27,const real3d_u& kb27,const real1d_u& sfluxref27,const real1d_u& rayl27,int layreffr27,real scalekur27,
    const real2d_u& absa28,const real4d_u& ka28,const real2d_u& absb28,const real4d_u& kb28,const real2d_u& sfluxref28,real rayl28,int layreffr28,real strrat28,
    const real2d_u& absa29,const real3d_u& ka29,const real2d_u& absb29,const real3d_u& kb29,const real2d_u& forref29,const real2d_u& selfref29,const real1d_u& sfluxref29,const real1d_u& absh2o29,const real1d_u& absco229,
    real rayl29,int layreffr29,const real4d& cldfmc,const real4d& ciwpmc,const real4d& clwpmc,const real4d& cswpmc,const real4d& taucmc,
    const real4d& taormc,const real4d& ssacmc,const real4d& asmcmc,const real4d& fsfcmc,const real4d& zcldfmc,const real4d& ztaucmc,const real4d& ztaormc,const real4d& zasycmc,const real4d& zomgcmc,
    const real4d& taua,const real4d& ssaa,const real4d& asma,const real4d& wkl,const real4d& ztaug,const real4d& ztaur,const real4d& zcd,const real4d& zcu,const real4d& zfd,const real4d& zfu,
    const real3d &zdbt,const real3d &zdbt_nodel,const real3d &zgc,const real3d &zgcc,const real3d &zgco,const real3d &zomc,const real3d &zomcc,const real3d &zomco,const real3d &zrdnd,const real3d &zrdndc,const real3d &zref,
    const real3d &zrefc,const real3d &zrefo,const real3d &zrefd,const real3d &zrefdc,const real3d &zrefdo,const real3d &zrup,const real3d &zrupd,const real3d &zrupc,const real3d &zrupdc,const real3d &zs1,const real3d &ztauc,const real3d &ztauo,const real3d &ztdn,const real3d &ztdnd,const real3d &ztdbt,const real3d &ztoc,const real3d &ztor,const real3d &ztra,const real3d &ztrac,const real3d &ztrao,
    const real3d &ztrad,const real3d &ztradc,const real3d &ztrado,const real3d &zdbtc,const real3d &ztdbtc,const real3d &zincflx,const real3d &zdbtc_nodel,const real3d &ztdbt_nodel,const real3d &ztdbtc_nodel,const real3d &zsflxzen,int m,int j){
    

    //! ----- Local -----
    //! Control

    int nlayers;             //! total number of layers
    int istart;              //! beginning band of calculation
    int iend;                //! ending band of calculation
    int icpr;                //! cldprop/cldprmc use flag
    int iout;                //! output option flag (inactive)
    int iaer;                //! aerosol option flag
    int idelm;               //! delta-m scaling flag (inactive)
    int isccos;              //! instrumental cosine response flag (inactive)
    int ims;                 //! value for changing mcica permute seed
    int imca;                //! flag for mcica [0=off, 1=on]

    real zepsec, zepzen;         //! epsilon
    real zdpgcp;                 //! flux to heating conversion ratio
    const int size = 110;
    //! Atmosphere
    real pavel_s[size + 2];
    real1d_u pavel((real*)pavel_s,{1},{size+2});                     //! layer pressures (mb) 
    real tavel_s[size + 2];
    real1d_u tavel((real*)tavel_s,{1},{size+2});                     //! layer temperatures (K)s
    real pz_s[size + 2];
    real1d_u pz((real*)pz_s,{0},{size+2});                           //! level (interface) pressures (hPa, mb)
    real tz_s[size + 2];
    real1d_u tz((real*)tz_s,{0},{size+2});                           //! level (interface) temperatures (K)     
    real tbound;                                                    //! surface temperature (K)
    real pdp_s[size + 2];
    real1d_u pdp((real*)pdp_s,{1},{size+2});                         //! layer pressure thickness (hPa, mb)
    real coldry_s[size + 2];
    real1d_u coldry((real*)coldry_s,{1},{size+2});                   //! dry air column amount
    real cossza;                                                     //! Cosine of solar zenith angle
    real adjflux_s[jpband + 1];
    real1d_u adjflux((real*)adjflux_s,{1},{jpband+1});                 //!adjustment for current Earth/Sun distance
    real solvar_s[jpband + 1];
    real1d_u solvar((real*)solvar_s,{1},{jpband+1});                   //! solar constant scaling factor from rrtmg_sw
    real albdir_s[nbndsw + 1];
    real1d_u albdir((real*)albdir_s,{1},{nbndsw+1});                   //! surface albedo, direct          ! zalbp
    real albdif_s[nbndsw + 1];
    real1d_u albdif((real*)albdif_s,{1},{nbndsw+1});                   //! surface albedo, diffuse         ! zalbd
    

    //! Atmosphere - setcoef
    int laytrop;             //! tropopause layer index
    int layswtch;            //! tropopause layer index
    int laylow;              //! tropopause layer index
    int jp_s[size + 2];
    int jt_s[size + 2];
    int jt1_s[size + 2];
    real colh2o_s[size + 2];
    real colco2_s[size + 2];
    real colo3_s[size + 2];
    real coln2o_s[size + 2];
    real colch4_s[size + 2];
    real colo2_s[size + 2];
    real colmol_s[size + 2];
    real co2mult_s[size + 2];
    int1d_u jp((int*)jp_s,{1},{size+2});                    
    int1d_u jt((int*)jt_s,{1},{size+2});
    int1d_u jt1((int*)jt1_s,{1},{size+2});
    real1d_u colh2o((real*)colh2o_s,{1},{size+2});           //! column amount (h2o)
    real1d_u colco2((real*)colco2_s,{1},{size+2});           //! column amount (co2)
    real1d_u colo3((real*)colo3_s,{1},{size+2});             //! column amount (o3)
    real1d_u coln2o((real*)coln2o_s,{1},{size+2});           //! column amount (n2o)
    real1d_u colch4((real*)colch4_s,{1},{size+2});           //! column amount (ch4)
    real1d_u colo2((real*)colo2_s,{1},{size+2});             //! column amount (o2)
    real1d_u colmol((real*)colmol_s,{1},{size+2});           //! column amount
    real1d_u co2mult((real*)co2mult_s,{1},{size+2});         //! column amount 
    int indself_s[size + 2];
    int indfor_s[size + 2];
    real selffac_s[size + 2];
    real selffrac_s[size + 2];
    real forfac_s[size + 2];
    real forfrac_s[size + 2];
    real fac00_s[size + 2];
    real fac01_s[size + 2];
    real fac10_s[size + 2];
    real fac11_s[size + 2];
    int1d_u indself((int*)indself_s,{1},{size+2}); 
    int1d_u indfor((int*)indfor_s,{1},{size+2}); 
    real1d_u selffac((real*)selffac_s,{1},{size+2}); 
    real1d_u selffrac((real*)selffrac_s,{1},{size+2}); 
    real1d_u forfac((real*)forfac_s,{1},{size+2}); 
    real1d_u forfrac((real*)forfrac_s,{1},{size+2}); 
    real1d_u fac00((real*)fac00_s,{1},{size+2}); 
    real1d_u fac01((real*)fac01_s,{1},{size+2}); 
    real1d_u fac10((real*)fac10_s,{1},{size+2}); 
    real1d_u fac11((real*)fac11_s,{1},{size+2}); 
    //! Atmosphere/clouds - cldprop
    int ncbands;             //! number of cloud spectral bands
    int inflag;              //! flag for cloud property method
    int iceflag;             //! flag for ice cloud properties
    int liqflag;             //! flag for liquid cloud properties

    //! Atmosphere/clouds - cldprmc [mcica]
    real relqmc_s[size + 2];
    real reicmc_s[size + 2];
    real resnmc_s[size + 2];
    real1d_u relqmc((real*)relqmc_s,{1},{size+2});                 //! liquid particle effective radius (microns)
    real1d_u reicmc((real*)reicmc_s,{1},{size+2});                 //! ice particle effective size (microns)
    real1d_u resnmc((real*)resnmc_s,{1},{size+2});                 //! snow particle effective size (microns)            
    //! Atmosphere/clouds/aerosol - spcvrt,spcvmc
    real zbbfu_s[size + 3];
    real zbbfd_s[size + 3];
    real zbbcu_s[size + 3];
    real zbbcd_s[size + 3];
    real zbbfddir_s[size + 3];
    real zbbcddir_s[size + 3];
    real zuvfd_s[size + 3];
    real zuvcd_s[size + 3];
    real zuvfddir_s[size + 3];
    real zuvcddir_s[size + 3];
    real znifd_s[size + 3];
    real znicd_s[size + 3];
    real znifddir_s[size + 3];
    real znicddir_s[size + 3];
    real swnflx_s[size + 3];
    real swnflxc_s[size + 3];
    real dirdflux_s[size + 3];
    real difdflux_s[size + 3];
    real uvdflx_s[size + 3];
    real nidflx_s[size + 3];
    real dirdnuv_s[size + 3];
    real difdnuv_s[size + 3];
    real dirdnir_s[size + 3];
    real difdnir_s[size + 3];
    

    real1d_u zbbfu((real*)zbbfu_s,{1},{size+3});                                     //! temporary upward shortwave flux (w/m2)
    real1d_u zbbfd((real*)zbbfd_s,{1},{size+3});                                     //! temporary downward shortwave flux (w/m2)
    real1d_u zbbcu((real*)zbbcu_s,{1},{size+3});                                     //! temporary clear sky upward shortwave flux (w/m2)
    real1d_u zbbcd((real*)zbbcd_s,{1},{size+3});                                     //! temporary clear sky downward shortwave flux (w/m2)
    real1d_u zbbfddir((real*)zbbfddir_s,{1},{size+3});                                  //! temporary downward direct shortwave flux (w/m2)s
    real1d_u zbbcddir((real*)zbbcddir_s,{1},{size+3});                                  //! temporary clear sky downward direct shortwave flux (w/m2)
    real1d_u zuvfd((real*)zuvfd_s,{1},{size+3});                                     //! temporary UV downward shortwave flux (w/m2)
    real1d_u zuvcd((real*)zuvcd_s,{1},{size+3});                                     //! temporary clear sky UV downward shortwave flux (w/m2)
    real1d_u zuvfddir((real*)zuvfddir_s,{1},{size+3});                                  //! temporary UV downward direct shortwave flux (w/m2)
    real1d_u zuvcddir((real*)zuvcddir_s,{1},{size+3});                                  //! temporary clear sky UV downward direct shortwave flux (w/m2)s
    real1d_u znifd((real*)znifd_s,{1},{size+3});                                     //! temporary near-IR downward shortwave flux (w/m2)
    real1d_u znicd((real*)znicd_s,{1},{size+3});                                     //! temporary clear sky near-IR downward shortwave flux (w/m2)
    real1d_u znifddir((real*)znifddir_s,{1},{size+3});                                  //! temporary near-IR downward direct shortwave flux (w/m2)
    real1d_u znicddir((real*)znicddir_s,{1},{size+3});                                  //! temporary clear sky near-IR downward direct shortwave flux (w/m2)
    
    //! Optional output fields      
    real1d_u swnflx((real*)swnflx_s,{1},{size+3});                                   //! Total sky shortwave net flux (W/m2)
    real1d_u swnflxc((real*)swnflxc_s,{1},{size+3});                                 //! Clear sky shortwave net flux (W/m2)
    real1d_u dirdflux((real*)dirdflux_s,{1},{size+3});                               //! Direct downward shortwave surface flux
    real1d_u difdflux((real*)difdflux_s,{1},{size+3});                               //! Diffuse downward shortwave surface flux
    real1d_u uvdflx((real*)uvdflx_s,{1},{size+3});                                   //! Total sky downward shortwave flux, UV/vis  
    real1d_u nidflx((real*)nidflx_s,{1},{size+3});                                   //! Total sky downward shortwave flux, near-IR s
    real1d_u dirdnuv((real*)dirdnuv_s,{1},{size+3});                                 //! Direct downward shortwave flux, UV/viss
    real1d_u difdnuv((real*)difdnuv_s,{1},{size+3});                                 //! Diffuse downward shortwave flux, UV/viss
    real1d_u dirdnir((real*)dirdnir_s,{1},{size+3});                                 //! Direct downward shortwave flux, near-IRs
    real1d_u difdnir((real*)difdnir_s,{1},{size+3});                                 //! Diffuse downward shortwave flux, near-IR
    real ztaua_s[(size + 1) * nbndsw];
    real zasya_s[(size + 1) * nbndsw];  
    real zomga_s[(size + 1) * nbndsw];      
    real2d_u ztaua((real*)ztaua_s,{1,1},{size+2,nbndsw+1});                //! total aerosol optical depths
    real2d_u zasya((real*)zasya_s,{1,1},{size+2,nbndsw+1});                //! total aerosol asymmetry parameter 
    real2d_u zomga((real*)zomga_s,{1,1},{size+2,nbndsw+1});                //! total aerosol single scattering albedo  

    // ! Initializations
    iout = 0;               //!BSINGH(PNNL) initializing iout to zero(Might be wrong!) as this variable is never initialized but used in spcvmc_sw
    zepsec = 1.e-06;
    zepzen = 1.e-10;
    // !jm not thread safe      oneminus = 1.0 - zepsec
    // !jm not thread safe      pi = 2. * asin(1.)

    istart = jpb1;
    iend = jpb2;
    icpr = 0;
    ims = 2;

    // ! In a GCM with or without McICA, set nlon to the longitude dimension
    // !
    // ! Set imca to select calculation type:
    // !  imca = 0, use standard forward model calculation (clear and overcast only)
    // !  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
    // !            (clear, overcast or partial cloud conditions)

    // ! *** This version uses McICA (imca = 1) ***

    // ! Set icld to select of clear or cloud calculation and cloud 
    // ! overlap method (read by subroutine readprof from input file INPUT_RRTM):  
    // ! icld = 0, clear only
    // ! icld = 1, with clouds using random cloud overlap (McICA only)
    // ! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
    // ! icld = 3, with clouds using maximum cloud overlap (McICA only)
    if (icld < 0 || icld > 3) icld = 2;

    // ! Set iaer to select aerosol option
    // ! iaer = 0, no aerosols
    // ! iaer = 6, use six ECMWF aerosol types
    // !           input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
    // ! iaer = 10, input total aerosol optical depth, single scattering albedo 
    // !            and asymmetry parameter (tauaer, ssaaer, asmaer) directly
    if ( aer_opt == 0 || aer_opt == 2 || aer_opt == 3)     iaer = 10;
    else if( aer_opt == 1 ) iaer = 6;
    // ! Call model and data initialization, compute lookup tables, perform
    // ! reduction of g-points from 224 to 112 for input absorption
    // ! coefficient data and other arrays.
    // !
    // ! In a GCM this call should be placed in the model initialization
    // ! area, since this has to be called only once.  
    // !      call rrtmg_sw_ini(cpdair)

    // ! This is the main longitude/column loop in RRTMG.
    // ! Modify to loop over all columns (nlon) or over daylight columns

    for(int iplon = 1; iplon <= ncol; iplon++){
        // ! Prepare atmosphere profile from GCM for use in RRTMG, and define
        // ! other input parameters

        inatm_sw (iplon, nlay, icld, iaer,
            play, plev, tlay, tlev, tsfc, h2ovmr,
            o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
            adjes, dyofyr, scon, inflgsw, iceflgsw, liqflgsw,
            cldfmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl, ciwpmcl, clwpmcl, cswpmcl,
            reicmcl, relqmcl, resnmcl, tauaer, ssaaer, asmaer,
            nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl,
            adjflux, solvar, inflag, iceflag, liqflag, cldfmc, taucmc,
            ssacmc, asmcmc, fsfcmc, ciwpmc, clwpmc, cswpmc, reicmc, relqmc, resnmc,
            taua, ssaa, asma,m,j);
        // !  For cloudy atmosphere, use cldprop to set cloud optical properties based on
        // !  input cloud physical properties.  Select method based on choices described
        // !  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
        // !  effective radius must be passed in cldprop.  Cloud fraction and cloud
        // !  optical properties are transferred to rrtmg_sw arrays in cldprop.
        cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc,
            ciwpmc, clwpmc, cswpmc, reicmc, relqmc, resnmc,
            taormc, taucmc, ssacmc, asmcmc, fsfcmc,m,j);
        icpr = 1;
        // ! Calculate coefficients for the temperature and pressure dependence of the 
        // ! molecular absorption coefficients by interpolating data from stored
        // ! reference atmospheres.
        setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl,
                laytrop, layswtch, laylow, jp, jt, jt1,
                co2mult, colch4, colco2, colh2o, colmol, coln2o,
                colo2, colo3, fac00, fac01, fac10, fac11,
                selffac, selffrac, indself, forfac, forfrac, indfor,m,j);
        // ! Cosine of the solar zenith angle 
        // !  Prevent using value of zero; ideally, SW model is not called from host model when sun 
        // !  is below horizon
        cossza = coszen(iplon,m,j);
        if (cossza <= zepzen) cossza = zepzen;

        // ! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 

        // ! Surface albedo
        // !  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
        for(int ib = 1;ib <= 9; ib++){
            albdir(ib) = aldir(iplon,m,j);
            albdif(ib) = aldif(iplon,m,j);
        }
        albdir(nbndsw) = aldir(iplon,m,j);
        albdif(nbndsw) = aldif(iplon,m,j);
        //!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
        for(int ib = 10; ib <= 13; ib++){
            albdir(ib) = asdir(iplon,m,j);
            albdif(ib) = asdif(iplon,m,j);
        }

        //! Clouds
        if (icld == 0){
            for(int ingp = 1; ingp <= ngptsw; ingp++){
                for(int ilay = 1; ilay <= nlay+1; ilay++){
                    zcldfmc(ilay,ingp,m,j) = 0.;
                    ztaucmc(ilay,ingp,m,j) = 0.;
                    ztaormc(ilay,ingp,m,j) = 0.;
                    zasycmc(ilay,ingp,m,j) = 0.;
                    zomgcmc(ilay,ingp,m,j) = 1.;
                }
            }
        }else if(icld >= 1){
            for(int ig = 1; ig<= ngptsw; ig++){
                for(int i = 1; i <= nlayers; i++){
                    zcldfmc(i,ig,m,j) = cldfmc(ig,i,m,j);
                    ztaucmc(i,ig,m,j) = taucmc(ig,i,m,j);
                    ztaormc(i,ig,m,j) = taormc(ig,i,m,j);
                    zasycmc(i,ig,m,j) = asmcmc(ig,i,m,j);
                    zomgcmc(i,ig,m,j) = ssacmc(ig,i,m,j);
                }
            }
        }

        // ! Aerosol
        // ! IAER = 0: no aerosols
        if (iaer == 0){
            for(int inbnd = 1; inbnd <= nbndsw; inbnd++){
                for(int ilay = 1; ilay <= nlay+1; ilay++){
                    ztaua(ilay,inbnd) = 0.;
                    zasya(ilay,inbnd) = 0.;
                    zomga(ilay,inbnd) = 1.;
                }
            }
        }
        // ! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
        // ! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer), 
        // ! or set manually here for each aerosol and layer.
        else if(iaer == 6){
            for(int ib = 1; ib<= nbndsw; ib++){
                for(int i = 1; i <= nlayers; i++){
                    ztaua(i,ib) = 0.;
                    zasya(i,ib) = 0.;
                    zomga(i,ib) = 0.;
                    for(int ia = 1; ia<=naerec; ia++){
                        ztaua(i,ib) = ztaua(i,ib) + rsrtaua[ia-1][ib-1] * ecaer(iplon,i,ia,m,j);
                        zomga(i,ib) = zomga(i,ib) + rsrtaua[ia-1][ib-1] * ecaer(iplon,i,ia,m,j) * rsrpiza[ia-1][ib-1];
                        zasya(i,ib) = zasya(i,ib) + rsrtaua[ia-1][ib-1] * ecaer(iplon,i,ia,m,j) * rsrpiza[ia-1][ib-1] * rsrasya[ia-1][ib-1];
                    }
                    if(zomga(i,ib) != 0.){
                        zasya(i,ib) = zasya(i,ib) / zomga(i,ib);
                    }
                    if (ztaua(i,ib) != 0.){
                        zomga(i,ib) = zomga(i,ib) / ztaua(i,ib);
                    }
                }
            }
        }

        //! IAER=10: Direct specification of aerosol optical properties from GCM
        else if(iaer==10){
            for(int ib = 1; ib<= nbndsw; ib++){
                for(int i = 1; i <= nlayers; i++){   
                    ztaua(i,ib) = taua(i,ib,m,j);
                    zasya(i,ib) = asma(i,ib,m,j);
                    zomga(i,ib) = ssaa(i,ib,m,j);
                }
            }
        }

        // ! Call the 2-stream radiation transfer model

        for(int i = 1; i <= nlayers+1; i++){
            zbbcu(i) = 0.;
            zbbcd(i) = 0.;
            zbbfu(i) = 0.;
            zbbfd(i) = 0.;
            zbbcddir(i) = 0.;
            zbbfddir(i) = 0.;
            zuvcd(i) = 0.;
            zuvfd(i) = 0.;
            zuvcddir(i) = 0.;
            zuvfddir(i) = 0.;
            znicd(i) = 0.;
            znifd(i) = 0.;
            znicddir(i) = 0.;
            znifddir(i) = 0.;
        }
        spcvmc_sw(nlayers, istart, iend, icpr, iout,
            pavel, tavel, pz, tz, tbound, albdif, albdir,
            zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc,
            ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, 
            laytrop, layswtch, laylow, jp, jt, jt1,
            co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3,
            fac00, fac01, fac10, fac11,
            selffac, selffrac, indself, forfac, forfrac, indfor,
            zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd,
            zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir,
            ztaug,ztaur,zcd,zcu,zfd,zfu,zdbt,zdbt_nodel,zgc,zgcc,zgco,zomc,zomcc,zomco,zrdnd,zrdndc,zref,
            zrefc,zrefo,zrefd,zrefdc,zrefdo,zrup,zrupd,zrupc,zrupdc,zs1,ztauc,ztauo,ztdn,ztdnd,ztdbt,ztoc,ztor,ztra,ztrac,ztrao,
            ztrad,ztradc,ztrado,zdbtc,ztdbtc,zincflx,zdbtc_nodel,ztdbt_nodel,ztdbtc_nodel,zsflxzen,
            absa16, ka16, absb16, kb16, forref16, selfref16, sfluxref16, rayl16, layreffr16, strrat16,
            absa17, ka17, absb17, kb17, forref17, selfref17, sfluxref17, rayl17, layreffr17, strrat17,
            absa18, ka18, absb18, kb18, forref18, selfref18, sfluxref18, rayl18, layreffr18, strrat18,
            absa19, ka19, absb19, kb19, forref19, selfref19, sfluxref19, rayl19, layreffr19, strrat19,
            absa20, ka20, absb20, kb20, forref20, selfref20, sfluxref20, absch420, rayl20, layreffr20,
            absa21, ka21, absb21, kb21, forref21, selfref21, sfluxref21, rayl21, layreffr21, strrat21,
            absa22, ka22, absb22, kb22, forref22, selfref22, sfluxref22, rayl22, layreffr22, strrat22,
            absa23, ka23, forref23, selfref23, sfluxref23, rayl23, layreffr23, givfac23,
            absa24, ka24, absb24, kb24, forref24, selfref24, sfluxref24, abso3a24, abso3b24, rayla24, raylb24, layreffr24, strrat24,
            absa25, ka25, sfluxref25, abso3a25, abso3b25, rayl25, layreffr25,
            sfluxref26, rayl26,
            absa27, ka27, absb27, kb27, sfluxref27, rayl27, layreffr27, scalekur27,
            absa28, ka28, absb28, kb28, sfluxref28, rayl28, layreffr28, strrat28,
            absa29, ka29, absb29, kb29, forref29, selfref29, sfluxref29, absh2o29, absco229, rayl29, layreffr29,m,j);
        // ! Transfer up and down, clear and total sky fluxes to output arrays.
        // ! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.
        for(int i = 1; i <= nlayers+1; i++){
            swuflxc(iplon,i,m,j) = zbbcu(i);
            swdflxc(iplon,i,m,j) = zbbcd(i);
            swuflx(iplon,i,m,j) = zbbfu(i);
            swdflx(iplon,i,m,j) = zbbfd(i);
            uvdflx(i) = zuvfd(i);
            nidflx(i) = znifd(i);

            // !  Direct/diffuse fluxes
            dirdflux(i) = zbbfddir(i);
            difdflux(i) = swdflx(iplon,i,m,j) - dirdflux(i);
            swdkdir(iplon,i,m,j) = dirdflux(i);                  //! all-sky direct flux   jararias, 2013/08/10
            swdkdif(iplon,i,m,j) = difdflux(i);                 //! all-sky diffuse flux  jararias, 2013/08/10

            //!  UV/visible direct/diffuse fluxes
            dirdnuv(i) = zuvfddir(i);
            difdnuv(i) = zuvfd(i) - dirdnuv(i);
            // !  ------- Zhenxin add vis/uv downwards dir or dif here --!
            sibvisdir(iplon,i,m,j) = dirdnuv(i);
            sibvisdif(iplon,i,m,j) = difdnuv(i);
            // !  ----- End of Zhenxin addition  ------------!
            // !  Near-IR direct/diffuse fluxes
            dirdnir(i) = znifddir(i);
            difdnir(i) = znifd(i) - dirdnir(i);
            // !  ---------Zhenxin add nir downwards dir and dif here --!
            sibnirdir(iplon,i,m,j) = dirdnir(i);
            sibnirdif(iplon,i,m,j) = difdnir(i);
            // !  --------    End of Zhenxin addition 2011-05  ---------!
        }

        //!  Total and clear sky net fluxes
        for(int i = 1; i <= nlayers+1; i++){
            swnflxc(i) = swdflxc(iplon,i,m,j) - swuflxc(iplon,i,m,j);
            swnflx(i) = swdflx(iplon,i,m,j) - swuflx(iplon,i,m,j);
        }

        //!  Total and clear sky heating rates
        for(int i = 1; i <= nlayers; i++){
            zdpgcp = heatfac / pdp(i);
            swhrc(iplon,i,m,j) = (swnflxc(i+1) - swnflxc(i)) * zdpgcp;
            swhr(iplon,i,m,j) = (swnflx(i+1) - swnflx(i)) * zdpgcp;

        }
        swhrc(iplon,nlayers,m,j) = 0.;
        swhr(iplon,nlayers,m,j) = 0.;
    }//End longitude loops
}
