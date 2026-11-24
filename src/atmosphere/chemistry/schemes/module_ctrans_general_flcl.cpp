#include <Kokkos_Core.hpp>
#include "flcl-cxx.hpp"
#include <math.h>
#include <string>
#include <cstring>

//3.10 DS

typedef double real;

typedef Kokkos::View<int*, Kokkos::LayoutLeft> int1d;

typedef Kokkos::View<int****, Kokkos::LayoutLeft> int4d;

typedef Kokkos::View<real*, Kokkos::LayoutLeft> real1d;

typedef Kokkos::View<real**, Kokkos::LayoutLeft> real2d;

typedef Kokkos::View<real***, Kokkos::LayoutLeft> real3d;

constexpr int CHEM_WORKSPACE_CHUNK_SIZE = 64;

struct ctrans_general_workspace {
    real2d aq_aer_ratio_comp;

    ctrans_general_workspace(int ncomp, int ncol)
        : aq_aer_ratio_comp("ctrans_ws_aq_aer_ratio_comp", ncomp, ncol) {}
};

static ctrans_general_workspace *ctrans_general_ws = nullptr;
static int ctrans_general_ws_ncomp = -1;
static int ctrans_general_ws_ncol = -1;

static bool ctrans_general_workspace_matches(int ncomp, int ncol) {
    return ctrans_general_ws != nullptr &&
           ctrans_general_ws_ncomp == ncomp &&
           ctrans_general_ws_ncol == ncol;
}

static ctrans_general_workspace &get_ctrans_general_workspace(
        int ncomp, int ncol) {
    if (!ctrans_general_workspace_matches(ncomp, ncol)) {
        if (ctrans_general_ws != nullptr) {
            Kokkos::fence();
            delete ctrans_general_ws;
        }
        ctrans_general_ws = new ctrans_general_workspace(ncomp, ncol);
        ctrans_general_ws_ncomp = ncomp;
        ctrans_general_ws_ncol = ncol;
    }
    return *ctrans_general_ws;
}

extern "C" void ctrans_general_workspace_finalize() {
    if (ctrans_general_ws != nullptr) {
        if (Kokkos::is_initialized()) {
            Kokkos::fence();
        }
        delete ctrans_general_ws;
        ctrans_general_ws = nullptr;
    }
    ctrans_general_ws_ncomp = -1;
    ctrans_general_ws_ncol = -1;
}

real static qcldwtr_cutoff = 1.0e-6;  //! kg / kg
real static mwdry = 28.966;           //! molecular mass of dry air (g/mol)
real static mwso4 = 96.00;            //! molecular mass of so4-- (g/mol)
real static mwno3 = 62.0;             //! molecular mass of no3- (g/mol)
real static mwnh4 = 18.0985;          //! molecular mass of nh4+ (g/mol)

extern "C" void c_ctrans_aqchem_bridge(
        real &temp, real &pres_pa, real &taucld, real &prcrate,
        real &wcavg, real &wtavg, real &airm, real &alfa0, real &alfa2,
        real &alfa3, real *gas, real *aerosol, real *liquid,
        real *gaswdep, real *aerwdep, real &hpwdep);

KOKKOS_INLINE_FUNCTION void cup_env(real *z,real *qes,real *he,real *hes,const real2d &t,const real2d &q,const real2d &p,const real1d & z1,
        const real1d &psur,int &ierr,real &tcrit,int itest,real &xl,real &cp,
        int &i,int &kts,int &kte );

KOKKOS_INLINE_FUNCTION void cup_env_clev(const real2d &t,real *qes,const real2d &q,real *he,real *hes,real *z,const real2d &p,real *qes_cup,real *q_cup,
        real *he_cup,real *hes_cup,real *z_cup,real *p_cup,real *gamma_cup,real *t_cup,const real1d &psur, 
        int &ierr,const real1d &z1,real &xl,real &rv,real &cp,                                
        int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_env_clev_tr(const real3d &tracer,real *tre_cup,int &num_chem,int &ierr,int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_maximi(real *array,int ks,int &ke,int &maxx,int &ierr,
              int &i, int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_kbcon(real &cap_inc,int iloop,int &k22,const int1d &kbcon,real *he_cup,real *hes_cup,
        int &ierr,int &kbmax,real *p_cup, real &cap_max,
        int &i, int &kts, int &kte);

KOKKOS_INLINE_FUNCTION void cup_minimi(real *array,const int1d &ks,int &kend,int &kt,int &ierr, 
        int &i,int & kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_minimi_l(real *array,int &ks,int &kend,int &kt,int &ierr, 
int &i,int & kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_jmin_fix_single_column(int &i,int &kte,int& jmin,int& kdet,const int1d &ktop,
    int& ierr,real *hes_cup,real *z_cup,real *hcd
);

KOKKOS_INLINE_FUNCTION void cup_up_he(int &k22,real &hkb,real *z_cup,real *cd,real &entr,real *he_cup,real *hc, 
        const int1d &kbcon,int &ierr,real *dby,real *he,real *hes_cup,             
        int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_ktop(int ilo,real *dby,const int1d &kbcon,const int1d &ktop,int &ierr, 
        int &i, int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_up_nms(real *zu,real *z_cup,real &entr,real *cd,const int1d &kbcon,const int1d &ktop,int &ierr,const int1d &k22,
        int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_up_nms_l(real *zu,real *z_cup,real &entr,real *cd,const int1d &kbcon,const int1d &ktop,int &ierr,int &k22,
        int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_dd_nms(real *zd,real *z_cup,real *cdd,real &entr,int &jmin,int &ierr,
        int itest,int &kdet,const real1d &z1,
        int &i, int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_dd_he(real *hes_cup,real *zd,real *hcd,real *z_cup,real *cdd,real &entr,
        int &jmin,int &ierr,real *he,real *dby,real *he_cup,
        int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_dd_moisture(real *zd,real *hcd,real *hes_cup,real *qcd,real *qes_cup,
        real *pwd,real *q_cup,real *z_cup,real *cdd,real &entr,int &jmin,int &ierr,
        real *gamma_cup,real &pwev,real &bu,real *qrcd,
        const real2d &q,real *he,real *t_cup,int iloop,real &xl,
        int &i, int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_up_moisture(int &ierr,real *z_cup,real *qc,real *qrc,real *pw,real &pwav,
        const int1d &kbcon,const int1d &ktop,real *cd,real *dby,real &mentr_rate,real *clw_all,
        const real2d &q,real *gamma_cup,real *zu,real *qes_cup,int &k22,real *qe_cup,real &xl,
        int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_dd_edt(int &ierr,const real2d &us,const real2d &vs,real *z,const int1d &ktop,const int1d &kbcon,real &edt,const real2d &p,real &pwav,
        real &pwev,real &edtmax,real &edtmin,int maxens2,real *edtc,
        int &i,int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_up_tracer(int &ierr, real &tcrit,const real2d &t, const real1d &pre,real *z_cup,const real2d &p,const real3d &tracer,real *tre_cup,real *tr_up,
        real *tr_pw,real *tr_c,const real2d &hstary,real *cupclw,real *clw_all,const int1d &kbcon,const int1d &ktop,real *cd,real &mentr_rate,real *zu,
        const int1d &k22,int &num_gas,char *chem_opt,bool &cu_scav_on,int &conv_tr_aqchem,const real3d &chem,const real2d &aq_aer_ratio,
        const real2d &aq_aer_ratio_comp_workspace, int workspace_col,
        int &i, int &kts, int &kte, int &ipr, int &jpr, int &j, int &npr, int &num_chem,const char* name,const int1d &ncomp_aer,
        int &mastercompindx_cl_aer, int &mastercompindx_na_aer, int &mastercompindx_oc_aer, int &mastercompindx_bc_aer, int &mastercompindx_so4_aer, 
        int &mastercompindx_oin_aer, int &mastercompindx_dust_aer, int &mastercompindx_nh4_aer, int &mastercompindx_msa_aer,
        int &mastercompindx_no3_aer, int &mastercompindx_ca_aer, int &mastercompindx_co3_aer,
        const int1d &nsize_aer,const int4d &massptr_aer,const real2d &dens_aer,
        int &p_so2,int &p_hno3,int &p_hno4,int &p_dms,int &p_n2o5,int &p_nh3,int &p_h2o2,int &p_o3,int &p_ora1,
        int &p_op1,int &p_paa,int &p_h2so4,int &p_so4aj,int &p_nh4aj ,int &p_no3aj);

KOKKOS_INLINE_FUNCTION void cup_up_tracer_l(int &ierr, real &tcrit,const real2d &t, const real1d &pre,real *z_cup,const real2d &p,const real3d &tracer,real *tre_cup,real *tr_up,
        real *tr_pw,real *tr_c,const real2d &hstary,real *cupclw,real *clw_all,const int1d &kbcon,const int1d &ktop,real *cd,real &mentr_rate,real *zu,
        int &k22,int &num_gas,char *chem_opt,bool &cu_scav_on,int &conv_tr_aqchem,const real3d &chem,const real2d &aq_aer_ratio,
        const real2d &aq_aer_ratio_comp_workspace, int workspace_col,
        int &i, int &kts, int &kte, int &ipr, int &jpr, int &j, int &npr, int &num_chem,const char* name,const int1d &ncomp_aer,
        int &mastercompindx_cl_aer, int &mastercompindx_na_aer, int &mastercompindx_oc_aer, int &mastercompindx_bc_aer, int &mastercompindx_so4_aer, 
        int &mastercompindx_oin_aer, int &mastercompindx_dust_aer, int &mastercompindx_nh4_aer, int &mastercompindx_msa_aer,
        int &mastercompindx_no3_aer, int &mastercompindx_ca_aer, int &mastercompindx_co3_aer,
        const int1d &nsize_aer,const int4d &massptr_aer,const real2d &dens_aer,
        int &p_so2,int &p_hno3,int &p_hno4,int &p_dms,int &p_n2o5,int &p_nh3,int &p_h2o2,int &p_o3,int &p_ora1,
        int &p_op1,int &p_paa,int &p_h2so4,int &p_so4aj,int &p_nh4aj ,int &p_no3aj);

KOKKOS_INLINE_FUNCTION void cup_dellas_tr(int &ierr,real *z_cup,real *p_cup,real *tr_dd,real &edt,real *zd,real *cdd,
                    const real3d &tracer, const real3d &tracert,int &j,real &mentrd_rate,real *zu,real &g,const real1d &xmb,
                    real *cd,real *tr_up,const int1d &ktop,const int1d &k22, const int1d &kbcon, real &mentr_rate,int &jmin,
                    real *tre_cup,int &kdet, const int1d &kpbl,
                    int &ipr, int &jpr, int &npr,const char* name, int &num_chem,
                    int &i, int &kts, int &kte);

KOKKOS_INLINE_FUNCTION void cup_dellas_tr_l(int &ierr,real *z_cup,real *p_cup,real *tr_dd,real &edt,real *zd,real *cdd,
                           const real3d &tracer, const real3d &tracert,int &j,real &mentrd_rate,real *zu,real &g,real &xmb,
                           real *cd,real *tr_up,const int1d &ktop,int &k22, const int1d &kbcon, real &mentr_rate,int &jmin,
                           real *tre_cup,int &kdet,int &kpbl,
                           int &ipr, int &jpr, int &npr,const char* name, int &num_chem,
                           int &i, int &kts, int &kte);

KOKKOS_INLINE_FUNCTION void cup_dd_tracer(int &ierr,real *z_cup,real *qrcd, const real3d &tracer,real *tre_cup,real *tr_up,real *tr_dd,
         real *tr_pw,real *tr_pwd,int &jmin,real *cdd,real &entr,real *zd,real *pwdper,const real3d &wetdep,real &xmb,int &k22,const real2d &rho_phy,const real2d &dz8w,
         int &num_gas,int &num_chem,
         int &i, int &kts,int &kte);

KOKKOS_INLINE_FUNCTION void cup_dellabot_tr(int &ipr, int &jpr, real *tre_cup,int &ierr,real *z_cup,const real2d &p_cup,
                              real *tr_dd,real &edt,real *zd,real *cdd,const real3d &tracer, const real3d &tracert,int &j,
                              real &mentrd_rate,real *z, real &g,real &xmb,
                              int &num_chem,
                              int &i, int &kts, int &kte);

// KOKKOS_INLINE_FUNCTION void aqchem(const real2d &t,const real2d &q,const real2d &p,const real2d &z,const real3d &chem,const real2d &aq_aer_ratio,
//         const int1d &ierr,int &num_chem, int &num_gas, char *chem_opt, bool &cu_scav_on, int &conv_tr_aqchem,
//         int &its,int &ite,int &kts,int &kte);

extern "C" void c_cup_ct(flcl_ndarray_t* ktop_p,flcl_ndarray_t* kbcon_p,flcl_ndarray_t* k23_p,flcl_ndarray_t* kbcon3_p,flcl_ndarray_t* ktop3_p,
                        flcl_ndarray_t*xmb3_p,flcl_ndarray_t* tracer_p,int &j,flcl_ndarray_t* aaeq_p,flcl_ndarray_t* t_p,flcl_ndarray_t* q_p,flcl_ndarray_t* z1_p,
                        flcl_ndarray_t* pre_p,flcl_ndarray_t* p_p,flcl_ndarray_t* tracert_p,flcl_ndarray_t* tracert3_p,flcl_ndarray_t* hstary_p,
                        real &dtime,flcl_ndarray_t* psur_p,flcl_ndarray_t* us_p,flcl_ndarray_t* vs_p,real &tcrit,    
                        flcl_ndarray_t* wetdep_p,flcl_ndarray_t* rho_phy_p,flcl_ndarray_t* dz8w_p,real &xl,real &rv,real &cp,real &g,int &ipr,int &jpr,int &npr,int &num_chem,char *chem_opt,         
                        bool cu_scav_on,int &conv_tr_aqchem,flcl_ndarray_t* chem_p,flcl_ndarray_t* aq_aer_ratio_p,           
                        int &ishallow,int &num_gas,int &its,int &ite,int &kts,int &kte,
                        flcl_ndarray_t* nsize_aer_p,flcl_ndarray_t* ncomp_aer_p,flcl_ndarray_t* massptr_aer_p,flcl_ndarray_t* dens_aer_p,
                        int &mastercompindx_cl_aer, int &mastercompindx_na_aer, int &mastercompindx_oc_aer, int &mastercompindx_bc_aer, int &mastercompindx_so4_aer, 
                        int &mastercompindx_oin_aer, int &mastercompindx_dust_aer, int &mastercompindx_nh4_aer, int &mastercompindx_msa_aer,
                        int &mastercompindx_no3_aer, int &mastercompindx_ca_aer, int &mastercompindx_co3_aer, 
                        int &p_co,int &p_so2,int &p_hno3,int &p_hno4,int &p_dms,int &p_n2o5,int &p_nh3,int &p_h2o2,int &p_o3,int &p_ora1 ,
                        int &p_op1,int &p_paa,int &p_h2so4,int &p_so4aj,int &p_nh4aj ,int &p_no3aj){

    auto ktop      = flcl::view_from_ndarray<int*>(*ktop_p);
    auto kbcon     = flcl::view_from_ndarray<int*>(*kbcon_p);
    auto k23       = flcl::view_from_ndarray<int*>(*k23_p);
    auto kbcon3    = flcl::view_from_ndarray<int*>(*kbcon3_p);
    auto ktop3     = flcl::view_from_ndarray<int*>(*ktop3_p);

    auto xmb3      = flcl::view_from_ndarray<real*>(*xmb3_p);
    auto tracer    = flcl::view_from_ndarray<real***>(*tracer_p);
    auto aaeq      = flcl::view_from_ndarray<real*>(*aaeq_p);
    auto t         = flcl::view_from_ndarray<real**>(*t_p);
    auto q         = flcl::view_from_ndarray<real**>(*q_p);
    auto z1        = flcl::view_from_ndarray<real*>(*z1_p);
    auto pre       = flcl::view_from_ndarray<real*>(*pre_p);
    auto p         = flcl::view_from_ndarray<real**>(*p_p);
    auto tracert    = flcl::view_from_ndarray<real***>(*tracert_p);
    auto tracert3   = flcl::view_from_ndarray<real***>(*tracert3_p);
    auto hstary     = flcl::view_from_ndarray<real**>(*hstary_p);
    auto psur      = flcl::view_from_ndarray<real*>(*psur_p);
    auto us        = flcl::view_from_ndarray<real**>(*us_p);
    auto vs        = flcl::view_from_ndarray<real**>(*vs_p);
    auto wetdep    = flcl::view_from_ndarray<real***>(*wetdep_p);
    auto rho_phy   = flcl::view_from_ndarray<real**>(*rho_phy_p);
    auto dz8w   = flcl::view_from_ndarray<real**>(*dz8w_p);
    auto chem      = flcl::view_from_ndarray<real***>(*chem_p);
    auto aq_aer_ratio = flcl::view_from_ndarray<real**>(*aq_aer_ratio_p);

    //variable form module_data_binaer_asect
    auto nsize_aer = flcl::view_from_ndarray<int*>(*nsize_aer_p);
    auto ncomp_aer = flcl::view_from_ndarray<int*>(*ncomp_aer_p);
    auto massptr_aer = flcl::view_from_ndarray<int****>(*massptr_aer_p);

    auto dens_aer = flcl::view_from_ndarray<real**>(*dens_aer_p);

#if defined(KOKKOS_ENABLE_SERIAL)
    const int workspace_ncol = 1;
#else
    const int workspace_ncol = CHEM_WORKSPACE_CHUNK_SIZE;
#endif
    auto &workspace = get_ctrans_general_workspace(
        ncomp_aer(0), workspace_ncol);
    auto aq_aer_ratio_comp_workspace = workspace.aq_aer_ratio_comp;

    // loacl variables
    real edtc[3];
    real he[kte];
    real hes[kte];
    real qes[kte];
    real z[kte];
    real pwdper[kte];
    real qes_cup[kte];
    real q_cup[kte];
    real he_cup[kte];
    real hes_cup[kte];
    real z_cup[kte];
    real p_cup[kte];
    real gamma_cup[kte];
    real t_cup[kte];
    real dby[kte];
    real qc[kte];
    real qrcd[kte];
    real pwd[kte];
    real pw[kte];
    real hcd[kte];
    real qcd[kte];
    real dbyd[kte];
    real hc[kte];
    real qrc[kte];
    real zu[kte];
    real zd[kte];
    real clw_all[kte];
    real zd3[kte];
    real clw_all3[kte];
    real cd3[kte];
    real pw3[kte];
    real zu3[kte];
    real cd[kte];
    real cdd[kte];
    real cdd3[kte];

    real edt ,hkb ,edt3 ,xmb ,pwav ,pwev ,bu ,cap_max ,cap_max_increment;

    real tr3_c[num_chem*kte];
    real tr3_up[num_chem*kte];
    real tr3_pw[num_chem*kte];
    real trkb3[num_chem];
    real tr_c[num_chem*kte];
    real tr_up[num_chem*kte];
    real tr_dd[num_chem*kte];
    real tr_dd3[num_chem*kte];
    real tre_cup[num_chem*kte];
    real tr_pw[num_chem*kte];
    real tr_pwd[num_chem*kte];
    int kzdown ,kdet ,k22 ,jmin ,kstabi ,kstabm ,ierr ,kbmax ,ierr5;

    int ki, i, k, kk;
    int itf, jtf, ktf;
    real day, dz, mbdt, entr_rate, radius;
    real entrd_rate, mentr_rate, mentrd_rate;
    real zcutdown, edtmax, edtmin, depth_min;
    real zkbmax, z_detr, zktop;
    real dh, cap_maxs, entr_rate3, mentr_rate3;

    char* message;
    //!sms$distribute end
    day = 86400.;
    //!--- specify entrainmentrate and detrainmentrate
    radius = 12000.;
    // !--- gross entrainment rate (these may be changed later on in the
    // !--- program, depending what your detrainment is!!)
    entr_rate = .2 / radius;
    entr_rate3 = .2 / 200.;
    //!--- entrainment of mass
    mentrd_rate = 0.;
    mentr_rate = entr_rate;
    mentr_rate3 = entr_rate3;
    //!--- initial detrainmentrates

    // !--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
    // !    base mass flux
    edtmax = .8;
    edtmin = .2;
    // !--- minimum depth (m), clouds must have
    depth_min = 500.;
    //!--- maximum depth (mb) of capping
    //!--- inversion (larger cap = no convection)
    cap_maxs = 175.;

    //!--- max height(m) above ground where updraft air can originate
    zkbmax = 4000.;
    //!--- height(m) above which no downdrafts are allowed to originate
    zcutdown = 3000.;
    //!--- depth(m) over which downdraft detrains all its mass
    z_detr = 1250.;
    mbdt=dtime*4.E-03;

    //opt
    for(int i = 0;i < ite;i++){         // i-loop start
        // Columns are processed serially here, so slots can be reused safely
        // when a local block contains more than the fixed workspace width.
        const int workspace_col = i % workspace_ncol;
        // These variables are column-indexed arrays in the Fortran source.
        // Reset the scalar workspace for every fused single-column pass.
        kbmax = 1;
        cap_max_increment = 0.;
        edt = 0.;
        edt3 = 0.;
        kstabm = kte - 1;
        ierr = 0;
        ierr5 = 0;

        for(int k = 0;k < kte;k++){
            cd[k] = 0.1 * entr_rate;
            cd3[k] = entr_rate3;
            cdd[k] = 0.;
            cdd3[k] = 0.;
            clw_all[k] = 0.;
        }

        for(int k = 0;k < kte;k++){
            for(int ki = 0;ki < num_chem;ki++){
                tr_dd3[ki*kte+k] = 0.;
            }
        }

        //!sms$to_local(grid_dh: <1, mix :size>, <2, mjx :size>) begin

        if (ktop3(i) < 2)ierr5 = 20;
        if (xmb3(i) == 0.)ierr5 = 21;

        cap_max = cap_maxs;
        for(int k = 0;k < kte;k++){
            zd3[k] = 0.;
        }

        //!--- calculate moist static energy, heights, qes

        cup_env(z,qes,he,hes,t,q,p,z1,
        psur,ierr,tcrit,0,xl,cp,
        i,kts,kte);

        // !--- environmental values on cloud levels

        cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup,
            hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,
            ierr,z1,xl,rv,cp,
            i, kts,kte);

        cup_env_clev_tr(tracer,tre_cup,num_chem,ierr,i,kts,kte);

        if (ierr == 0){
            for(int k = 0;k < kte-2;k++){
                if (z_cup[k] > zkbmax + z1(i)){
                  kbmax = k + 1;
                  break;//exit ! jump out loop k
                }
            }
            for(int k = 0;k < kte;k++){
                if (z_cup[k] > z_detr + z1(i)){
                  kdet = k + 1;
                  break;//exit ! jump out loop k
                }
            }
        }

        //!--- determine level with highest moist static energy content - k22

        cup_maximi(he_cup,3,kbmax,k22,ierr,
            i, kts,kte);

        if (ierr == 0.) if (k22 >= kbmax) ierr = 2;

        if (aaeq(i) != 0.) ierr = 20;

        //!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
        cup_kbcon(cap_max_increment,1,k22,kbcon,he_cup,hes_cup,
            ierr,kbmax,p_cup,cap_max,
            i, kts,kte);

        //!--- increase detrainment in stable layers
        cup_minimi(hes_cup,kbcon,kstabm,kstabi,ierr,
            i, kts,kte);

        if (ierr == 0.){
            if (kstabm - 1 > kstabi){
                for(int k = kstabi-1; k < kstabm - 1;k++){
                    cd[k] = cd[k - 1] + 1.5 * entr_rate;
                    if (cd[k] > 10.0 * entr_rate){
                        cd[k] = 10.0 * entr_rate;
                    }
                }
            }
        }

        cup_up_he(k22,hkb,z_cup,cd,mentr_rate,he_cup,hc,
           kbcon,ierr,dby,he,hes_cup,
           i, kts,kte);

        cup_ktop(1,dby,kbcon,ktop,ierr,i, kts,kte);

        kzdown = 0;
        if (ierr == 0){
            zktop = (z_cup[ktop(i)-1] - z1(i)) * .6;
            zktop = std::min(zktop + z1(i), zcutdown + z1(i));
            for(int k = 0;k < kte;k++){
                if (z_cup[k] > zktop){
                    kzdown = k + 1;
                    break;
                }
            }
        }

        //   !--- doWNDRAFT ORIGINATING LEVEL - JMIN
        cup_minimi_l(hes_cup,k22,kzdown,jmin,ierr,i, kts,kte);

        cup_jmin_fix_single_column(i,kte,jmin,kdet,ktop,ierr,
            hes_cup,z_cup,hcd);

        // - Must have at least depth_min m between cloud convective base
        //     and cloud top.
        if (ierr == 0.0f) {
            if (jmin <= 3) {
                ierr = 9;
            }
            if (z_cup[ktop(i)-1] - z_cup[kbcon(i)-1] < depth_min) {
                ierr = 6;
            }
        }

        //   !--- normalized updraft mass flux profile
        cup_up_nms_l(zu,z_cup,mentr_rate,cd,kbcon,ktop,ierr,k22,
            i, kts,kte);

        //   !--- normalized downdraft mass flux profile,also work on
        //   !--- bottom detrainment in this routine

        cup_dd_nms(zd,z_cup,cdd,mentrd_rate,jmin,ierr,
            0,kdet,z1,
            i, kts,kte);

        //   !--- downdraft moist static energy

        cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,mentrd_rate,
            jmin,ierr,he,dbyd,he_cup,
            i, kts,kte);

        //   !--- calculate moisture properties of downdraft

        cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup,
            pwd,q_cup,z_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup,
            pwev,bu,qrcd,q,he,t_cup,1,xl,
            i, kts,kte);

        //   !--- calculate moisture properties of updraft

        cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav,
            kbcon,ktop,cd,dby,mentr_rate,clw_all,
            q,gamma_cup,zu,qes_cup,k22,q_cup,xl,
            i, kts,kte);

        //   !--- DETERMINE doWNDRAFT STRENGTH IN TERMS OF WINDSHEAR

        cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav,
            pwev,edtmax,edtmin,3,edtc,
            i, kts,kte);

        if (ierr == 0) edt = edtc[1];

        //! massflux from precip and normalized cloud properties
        for(int k = 0; k < kte; k++){
            pwdper[k]= 0.;
        }

        if (ierr > 0) pre(i) = 0.;
        if (ierr == 0){
            xmb = pre(i) / (pwav + edt * pwev);
            if (i == ipr){
            //    write(message, "(a,x,3F9.3)") 'xmb, edt, pwav = ', xmb(i), edt(i), pwav(i)
            //    call chemistry_message(trim(message))
            //    write(message, "(a,x,F9.3)") 'k, pwdper(k, i), pw, pwd(k, i)', z1(i)
            //    call chemistry_message(trim(message))
            }
            for(int k =0; k < ktop(i); k++){
                pwdper[k] = - edt * pwd[k] / pwav;
                if (i == ipr){
                //   write(message, "(I3,x,3F9.3)") , k, pwdper(k, i), pw(k, i), pwd(k, i)
                //   call chemistry_message(trim(message))
                }
            }
        }
    
        // !!! NOW WE HAVE EVREYTHING TO CALCULATE TRACER TRANSPORT AND WET DEPOSITION !!!
        //   !--- calculate incloud tracer distribution

        cup_up_tracer_l(ierr,tcrit,t,pre,z_cup,p,tracer,tre_cup,
        tr_up,tr_pw,tr_c,hstary,pw,
        clw_all,kbcon,ktop,cd,mentr_rate,zu,k22,
        num_gas,chem_opt,cu_scav_on,conv_tr_aqchem, chem, aq_aer_ratio,
        aq_aer_ratio_comp_workspace,workspace_col,
        i, kts,kte,
        ipr,jpr,j,npr,num_chem,"deep",ncomp_aer,
        mastercompindx_cl_aer, mastercompindx_na_aer, mastercompindx_oc_aer, mastercompindx_bc_aer, mastercompindx_so4_aer, 
        mastercompindx_oin_aer, mastercompindx_dust_aer, mastercompindx_nh4_aer, mastercompindx_msa_aer,
        mastercompindx_no3_aer, mastercompindx_ca_aer, mastercompindx_co3_aer,
        nsize_aer,massptr_aer,dens_aer,
        p_so2,p_hno3,p_hno4,p_dms,p_n2o5,p_nh3,p_h2o2,p_o3,p_ora1,
        p_op1,p_paa,p_h2so4,p_so4aj,p_nh4aj ,p_no3aj);

        if (ishallow == 1){
            cup_up_nms(zu3,z_cup,mentr_rate3,cd3,kbcon3,ktop3,ierr5,k23,
                i, kts,kte);

            cup_up_tracer(ierr5,tcrit,t,pre,z_cup,p,tracer,
                tre_cup,tr3_up,tr3_pw,tr3_c,hstary,pw3,
                clw_all3,kbcon3,ktop3,cd3,mentr_rate3,zu3,k23,
                num_gas,chem_opt,cu_scav_on,conv_tr_aqchem, chem, aq_aer_ratio, 
                aq_aer_ratio_comp_workspace,workspace_col,
                i, kts,kte,
                ipr,jpr,j,npr,num_chem,"shallow",ncomp_aer,
                mastercompindx_cl_aer, mastercompindx_na_aer, mastercompindx_oc_aer, mastercompindx_bc_aer, mastercompindx_so4_aer, 
                mastercompindx_oin_aer, mastercompindx_dust_aer, mastercompindx_nh4_aer, mastercompindx_msa_aer,
                mastercompindx_no3_aer, mastercompindx_ca_aer, mastercompindx_co3_aer,
                nsize_aer,massptr_aer,dens_aer,
                p_so2,p_hno3,p_hno4,p_dms,p_n2o5,p_nh3,p_h2o2,p_o3,p_ora1,
                p_op1,p_paa,p_h2so4,p_so4aj,p_nh4aj ,p_no3aj);

            cup_dellas_tr(ierr5,z_cup,p_cup,tr_dd3,edt3,zd3,cdd3,
                tracer,tracert3,j,mentrd_rate,zu3,g,xmb3,
                cd3,tr3_up,ktop3,k23,kbcon3,mentr_rate3,jmin,tre_cup,kdet,
                k23,ipr,jpr,npr,"shallow",num_chem,
                i, kts,kte);     
        }

        cup_dd_tracer(ierr,z_cup,qrcd,tracer,tre_cup,tr_up,tr_dd,
        tr_pw,tr_pwd,jmin,cdd,mentrd_rate,zd,pwdper,wetdep,xmb,k22,rho_phy,dz8w,
        num_gas,num_chem,
        i, kts,kte);

        cup_dellabot_tr(ipr,jpr,tre_cup,ierr,z_cup,p,tr_dd,edt,
        zd,cdd,tracer,tracert,j,mentrd_rate,z,g,xmb,
        num_chem,
        i, kts,kte);

        cup_dellas_tr_l(ierr,z_cup,p_cup,tr_dd,edt,zd,cdd,
        tracer,tracert,j,mentrd_rate,zu,g,xmb,
        cd,tr_up,ktop,k22,kbcon,mentr_rate,jmin,tre_cup,kdet,
        k22,ipr,jpr,npr,"deep",num_chem,
        i, kts,kte);
    }   //i -loop end
}

KOKKOS_INLINE_FUNCTION void cup_env(real *z,real *qes,real *he,real *hes,const real2d &t,const real2d &q,const real2d &p,const real1d & z1,
        const real1d &psur,int &ierr,real &tcrit,int itest,real &xl,real &cp,
        int &i,int &kts,int &kte ){
    //cup_env local variables
    //auto tv = real2d("tv",kte,ite);
    real tv[kte];
    real e, tvbar;
    real ae[2], be[2], ht[2];
    int iph;
    ht[0] = xl / cp;
    ht[1] = 2.834e6 / cp;
    be[0] = .622 * ht[0] / .286;
    ae[0] = be[0] / 273. + std::log(610.71);
    be[1] = .622 * ht[1] / .286;
    ae[1] = be[1] / 273. + std::log(610.71);

    for(int k = 0;k < kte; k++){
        if (ierr == 0){
            iph = 0;
            if (t(k,i) <= tcrit) iph = 1;
            e = std::exp(ae[iph] - be[iph] / t(k, i));
            qes[k] = .622 * e / (100. * p(k, i) - e);
            if (qes[k] <= 1.e-08) qes[k] = 1.e-08;
            if (q(k, i) > qes[k]) q(k, i) = qes[k];
            tv[k] = t(k, i) + .608 * q(k, i) * t(k, i);              
        }
    }

    if (itest != 2){
        if (ierr == 0){
            z[0] = std::max(0., z1(i)) - (std::log(p(0, i)) - std::log(psur(i))) * 287. * tv[0] / 9.81;
        }

        for(int k = 1;k < kte; k++){
            if(ierr == 0){
                tvbar = .5 * tv[k] + .5 * tv[k-1];
                z[k] = z[k-1] - (std::log(p(k,i)) - std::log(p(k-1,i)))
                        *287.*tvbar/9.81;
            }
        }
    }
    else{
        for(int k = 0;k < kte; k++){
            if (ierr == 0){
                z[k]=(he[k]-1004.*t(k,i)-2.5e6*q(k,i))/9.81;
                z[k] = std::max(1.e-3, z[k]);
            }
        }      
    }
    for(int k = 0;k < kte; k++){
        if (ierr == 0){
            if(itest == 0) he[k]=9.81*z[k]+1004.*t(k,i)+2.5e06*q(k,i);
            hes[k]=9.81*z[k]+1004.*t(k,i)+2.5e06*qes[k];
            if (he[k] >= hes[k]) he[k] = hes[k];
        }
    }
}
KOKKOS_INLINE_FUNCTION void cup_env_clev(const real2d &t,real *qes,const real2d &q,real *he,real *hes,real *z,const real2d &p,real *qes_cup,real *q_cup,
              real *he_cup,real *hes_cup,real *z_cup,real *p_cup,real *gamma_cup,real *t_cup,const real1d &psur, 
              int &ierr,const real1d &z1,real &xl,real &rv,real &cp,                                
              int &i,int &kts,int &kte                        ){

    for(int k = 1;k < kte; k++){
        if (ierr == 0){
            qes_cup[k]=.5*(qes[k-1]+qes[k]);
            q_cup[k]=.5*(q(k-1,i)+q(k,i));
            hes_cup[k]=.5*(hes[k-1]+hes[k]);
            he_cup[k]=.5*(he[k-1]+he[k]);
            if (he_cup[k] > hes_cup[k]) he_cup[k] = hes_cup[k];
            z_cup[k]=.5*(z[k-1]+z[k]);
            p_cup[k]=.5*(p(k-1,i)+p(k,i));
            t_cup[k]=.5*(t(k-1,i)+t(k,i));
            gamma_cup[k] = (xl/cp)*(xl/(rv*t_cup[k]*t_cup[k]))*qes_cup[k];
        }
    }
    if(ierr == 0){
        qes_cup[0]=qes[0];
        q_cup[0]=q(0,i);
        hes_cup[0]=hes[0];
        he_cup[0]=he[0];
        z_cup[0]=.5*(z[0]+z1(i));
        p_cup[0]=.5*(p(0,i)+psur(i));
        t_cup[0]=t(0,i);
        gamma_cup[0] = (xl/cp)*(xl/(rv*t_cup[0]*t_cup[0]))*qes_cup[0];
    }
}

KOKKOS_INLINE_FUNCTION void cup_env_clev_tr(const real3d &tracer,real *tre_cup,int &num_chem,int &ierr,int &i,int &kts,int &kte){

    for(int k = 1;k < kte; k++){
        for(int nv = 0;nv < num_chem;nv++){
            if(ierr == 0) tre_cup[nv*kte+k]=.5*(tracer(nv,k-1,i)+tracer(nv,k,i));
        }
    }
    for(int nv = 0;nv < num_chem;nv++){
        if(ierr == 0)
        tre_cup[nv*kte+ 0] = tracer(nv, 0, i);
    }
}

KOKKOS_INLINE_FUNCTION void cup_maximi(real *array,int ks,int &ke,int &maxx,int &ierr,
              int &i, int &kts,int &kte){
    real xar;
    real x;

    maxx = ks;
    if(ierr == 0.){
        x = array[ks-1];

        for(int k = ks-1;k < ke;k++){
            xar = array[k];
            if (xar >= x){
                x = xar;
                maxx = k+1;
                
            }
        }
    }
}
KOKKOS_INLINE_FUNCTION void cup_kbcon(real &cap_inc,int iloop,int &k22,const int1d &kbcon,real *he_cup,real *hes_cup,
                        int &ierr,int &kbmax,real *p_cup,real &cap_max,
                        int &i, int &kts, int &kte){
    //strat cup_kbcon
    real pbcdif, plus;

    kbcon(i)=1;
    if(ierr != 0) return;

    kbcon(i) = k22;

    bool found = false;
    while (!found) {
        while (he_cup[k22-1] < hes_cup[kbcon(i)-1]) {
            kbcon(i)++;
            if (kbcon(i) > (kbmax + 2)) {
                if (iloop < 4) {
                    ierr = 3;
                }
                return;
            }
        }

        if (kbcon(i) - k22 == 1) {
            found = true;
        } else {
            pbcdif = -p_cup[kbcon(i)-1] + p_cup[k22-1];
            plus = std::max(25.0, cap_max - static_cast<real>(iloop - 1) * cap_inc);
            if(iloop == 4) plus = cap_max;

            if (pbcdif > plus) {
                k22++;
                kbcon(i) = k22;
            } else {
                found = true;
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_minimi(real *array,const int1d &ks,int &kend,int &kt,int &ierr, 
              int &i,int & kts,int &kte         ){

    real x;
    int kstop;
    kt = ks(i);
    if(ierr == 0){
        x = array[ks(i)-1];
        kstop = std::max(ks(i) + 1, kend);

        for(int k = ks(i); k < kstop; k++){
            if (array[k] < x){
                x = array[k];
                kt = k+1;
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_minimi_l(real *array,int &ks,int &kend,int &kt,int &ierr, 
              int &i,int & kts,int &kte         ){

    real x;
    int kstop;
    kt = ks;
    if(ierr == 0){
        x = array[ks-1];
        kstop = std::max(ks + 1, kend);

        for(int k = ks; k < kstop; k++){
            if (array[k] < x){
                x = array[k];
                kt = k + 1;
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_jmin_fix_single_column(int &i,int &kte,int& jmin,int& kdet,const int1d &ktop,
    int& ierr,real *hes_cup,real *z_cup,real *hcd){
    if (ierr != 0) return;

    if (jmin <= 3) {
        ierr = 9;
        return;
    }

    bool redo_search = true;
    while (redo_search) {
        redo_search = false;

        if (jmin - 1 < kdet) kdet = jmin - 1;
        if (jmin >= ktop(i) - 1) jmin = ktop(i) - 2;

        int ki = jmin - 1;
        hcd[ki] = hes_cup[ki];

        real dh = 0.0;

        for (int k = ki - 1; k >= 0; --k) {
            hcd[k] = hes_cup[jmin - 1];
            real dz = z_cup[k + 1] - z_cup[k];
            dh += dz * (hcd[k] - hes_cup[k]);
           if (dh > 0.0f) {
                jmin--;

                if (jmin > 3) {
                    redo_search = true;
                    break;
                } else {
                    ierr = 9;
                    return;
                }
            }
        }
    }

    if (jmin <= 3) {
        ierr = 4;
    }
}

KOKKOS_INLINE_FUNCTION void cup_up_he(int &k22,real &hkb,real *z_cup,real *cd,real &entr,real *he_cup,real *hc, 
              const int1d &kbcon,int &ierr,real *dby,real *he,real *hes_cup,             
              int &i,int &kts,int &kte                  ){
    real dz;

    if (ierr == 0.){
        hkb = he_cup[k22-1];
        for(int k =0; k < k22; k++){
            hc[k] = he_cup[k];
            dby[k] = 0.;
        }
        for(int k = k22 - 1; k < kbcon(i) -1 ; k++){
            hc[k] = hkb;
            dby[k] = 0.;
        }
        int k = kbcon(i)-1;
        hc[k] = hkb;
        dby[kbcon(i)-1] = hkb - hes_cup[k];
    }

    for(int k = kts; k < kte; k++){
        if (k >= kbcon(i) && ierr == 0.){
            dz = z_cup[k] - z_cup[k - 1];
            hc[k]=(hc[k-1]*(1.-.5*cd[k]*dz)+entr * dz*he[k-1])/(1.+entr*dz-.5*cd[k]*dz);
            dby[k]= hc[k] - hes_cup[k];
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_ktop(int ilo,real *dby,const int1d &kbcon,const int1d &ktop,int &ierr, 
         int &i, int &kts,int &kte){

    ktop(i) = 1;

    if (ierr != 0) return;

    bool found_top = false;

    for (int k = kbcon(i); k < kte - 1; ++k) {
        if (dby[k] <= 0.0) {
            ktop(i) = k;
            found_top = true;
            break;
        }
    }
    if (!found_top) {
        if (ilo == 1) ierr = 5;
        return; 
    }
    for (int k = ktop(i); k < kte; ++k) {
        dby[k] = 0.0;
    }
}

KOKKOS_INLINE_FUNCTION void cup_up_nms(real *zu,real *z_cup,real &entr,real *cd,const int1d &kbcon,const int1d &ktop,int &ierr,const int1d &k22,
              int &i,int &kts,int &kte){
    real dz;

    for(int k = 0;k < kte;k++){
        zu[k] = 0.;
    }          

    if (ierr == 0){
        for(int k = k22(i)-1; k < kbcon(i); k++){
            zu[k] = 1.;
        }
        for(int k =kbcon(i); k < ktop(i); k++){
            dz = z_cup[k] - z_cup[k - 1];
            zu[k] = zu[k-1] * (1. + (entr  - cd[k]) * dz);
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_up_nms_l(real *zu,real *z_cup,real &entr,real *cd,const int1d &kbcon,const int1d &ktop,int &ierr,int &k22,
              int &i,int &kts,int &kte){
    real dz;

    for(int k = 0;k < kte;k++){
        zu[k] = 0.;
    }          

    if (ierr == 0){
        for(int k = k22-1; k < kbcon(i); k++){
            zu[k] = 1.;
        }
        for(int k =kbcon(i); k < ktop(i); k++){
            dz = z_cup[k] - z_cup[k - 1];
            zu[k] = zu[k-1] * (1. + (entr  - cd[k]) * dz);
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_dd_nms(real *zd,real *z_cup,real *cdd,real &entr,int &jmin,int &ierr,
              int itest,int &kdet,const real1d &z1,
              int &i,int &kts,int &kte){
    real a;           
    real perc =0.03;
    real dz;

    for(int k = 0;k < kte;k++){
        zd[k] = 0.;
        if (itest == 0) cdd[k] = 0.;
    }
    a = 1. - perc;

    if (ierr == 0){
        zd[jmin-1] = 1.;

        for(int ki = jmin-2; ki >= 0; ki--){
            dz = z_cup[ki + 1] - z_cup[ki];
            if (ki <= kdet-1 && itest == 0){
                cdd[ki]=entr+(1.-(a*(z_cup[ki]-z1(i))
                        +perc*(z_cup[kdet-1]-z1(i)))
                        /(a*(z_cup[ki+1]-z1(i))
                        +perc*(z_cup[kdet-1]-z1(i))))/dz;
            }
            zd[ki]=zd[ki+1]*(1.+(entr-cdd[ki])*dz);
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_dd_he(real *hes_cup,real *zd,real *hcd,real *z_cup,real *cdd,real &entr,
        int &jmin,int &ierr,real *he,real *dby,real *he_cup,
        int &i,int &kts,int &kte){
    int k;
    real dz;

    for(int k = 1;k < kte;k++){
        dby[k] = 0.;
        if (ierr == 0) hcd[k] = hes_cup[k];
    }

    if (ierr == 0){
        k = jmin-1;
        hcd[k] = hes_cup[k];
        dby[k] = hcd[jmin-1] - hes_cup[k];

        for(int ki = jmin - 2; ki >= 0; ki--){
            dz = z_cup[ki + 1] - z_cup[ki];
            hcd[ki]=(hcd[ki+1]*(1.-.5*cdd[ki]*dz)
                    +entr*dz*he[ki])
                    /(1.+entr*dz-.5*cdd[ki]*dz);
            dby[ki] = hcd[ki] - hes_cup[ki];
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_dd_moisture(real *zd,real *hcd,real *hes_cup,real *qcd,real *qes_cup,
              real *pwd,real *q_cup,real *z_cup,real *cdd,real &entr,int &jmin,int &ierr,
              real *gamma_cup,real &pwev,real &bu,real *qrcd,
              const real2d &q,real *he,real *t_cup,int iloop,real &xl,
              int &i, int &kts,int &kte){
    int k;
    real dh, dz, dqeva;
    bu = 0.;
    pwev = 0.;   

    for(k = 0;k < kte;k++){
        qcd[k] = 0.;
        qrcd[k] = 0.;
        pwd[k] = 0.;        
    }

    if (ierr == 0){
        k = jmin - 1;
        dz = z_cup[k + 1] - z_cup[k];
        qcd[k] = q_cup[k];
        qrcd[k] = qes_cup[k];
        pwd[jmin-1] = std::min(0., qcd[k] - qrcd[k]);
        pwev = pwev + pwd[jmin-1];
        qcd[k]= qes_cup[k];
        dh = hcd[k] - hes_cup[k];
        bu = dz * dh;
        for(int ki = jmin - 2; ki >= 0; ki--){
            dz = z_cup[ki + 1] - z_cup[ki];
            qcd[ki]=(qcd[ki+1]*(1.-.5*cdd[ki]*dz)
                    +entr*dz*q(ki,i))
                    /(1.+entr*dz-.5*cdd[ki]*dz);
            dh = hcd[ki] - hes_cup[ki];
            bu = bu + dz * dh;
            qrcd[ki]=qes_cup[ki]+(1./xl)*(gamma_cup[ki]
                    /(1.+gamma_cup[ki]))*dh;
            dqeva = qcd[ki] - qrcd[ki];
            if (dqeva > 0.) dqeva = 0.;
            pwd[ki] = zd[ki] * dqeva;
            qcd[ki] = qrcd[ki];
            pwev = pwev + pwd[ki];
        }

        if (pwev == 0 && iloop == 1) ierr = 7;
        if (bu >= 0 && iloop == 1) ierr = 7;
    }
}

KOKKOS_INLINE_FUNCTION void cup_up_moisture(int &ierr,real *z_cup,real *qc,real *qrc,real *pw,real &pwav,
              const int1d &kbcon,const int1d &ktop,real *cd,real *dby,real &mentr_rate,real *clw_all,
              const real2d &q,real *gamma_cup,real *zu,real *qes_cup,int &k22,real *qe_cup,real &xl,
              int &i,int &kts,int &kte){
    
    real  dh, qrch, dz, radius;
    int iall = 0;
    real c0 = 0.002;
    
    if (mentr_rate > 0.){
        radius = .2 / mentr_rate;
        if (radius < 900.) c0 = 0.;
    }

    pwav = 0.;

    for(int k = 0;k < kte;k++){
        pw[k] = 0.;
        if (ierr == 0) qc[k] = qes_cup[k];
        clw_all[k] = 0.;
        qrc[k] = 0.;
    }
    if (ierr == 0.){
        for(int k = k22-1; k < kbcon(i)-1;k++) qc[k] = qe_cup[k22-1];
    }

    for(int k = 1;k < kte;k++){
        if (ierr != 0) continue;
        if (k < kbcon(i)-1) continue;
        if (k > ktop(i)-1) continue;
        dz = z_cup[k] - z_cup[k - 1];
        qc[k]=(qc[k-1]*(1.-.5*cd[k]*dz)+mentr_rate*
                dz*q(k-1,i))/(1.+mentr_rate*dz-.5*cd[k]*dz);
        qrch=qes_cup[k]+(1./xl)*(gamma_cup[k]
                /(1.+gamma_cup[k]))*dby[k];
        clw_all[k] = qc[k] - qrch;
        qrc[k]=(qc[k]-qrch)/(1.+c0*dz*zu[k]);
        if (qrc[k] < 0.)  qrc[k] = 0.;

        pw[k] = c0 * dz * qrc[k] * zu[k];
        if (iall == 1){
            qrc[k] = 0.;
            pw[k] = (qc[k] - qrch) * zu[k];
            if (pw[k] < 0.) pw[k] = 0.;
        }

        qc[k] = qrc[k] + qrch;

        pwav = pwav + pw[k];
    }
}

KOKKOS_INLINE_FUNCTION void cup_dd_edt(int &ierr,const real2d &us,const real2d &vs,real *z,const int1d &ktop,const int1d &kbcon,real &edt,const real2d &p,real &pwav,
              real &pwev,real &edtmax,real &edtmin,int maxens2,real *edtc,
              int &i,int &kts,int &kte){
    
    real    einc, pef, pefb, prezk, zkbc;
    real  vshear,sdp,vws;

    edt = 0.;
    vws = 0.;
    sdp = 0.;
    vshear = 0.;

    for(int kk = 0;kk < kte-1;kk++){
        if (ierr != 0) continue;
        if (kk <= ktop(i)-1 && kk >= kbcon(i)-1){
            vws = vws + (abs((us(kk+1,i)-us(kk,i))/(z[kk+1]-z[kk]))+   abs((vs(kk+1,i)-vs(kk,i))/(z[kk+1]-z[kk]))) *(p(kk,i) - p(kk+1,i));
            sdp = sdp + p(kk,i) - p(kk+1,i);
        }
        if (kk == kte-2) vshear = 1.e3 * vws / sdp;
    }

    if (ierr == 0){

        pef=(1.591-.639*vshear+.0953*(vshear*vshear)
                -.00496*(vshear*vshear*vshear));
        if (pef > edtmax) pef = edtmax;
        if (pef < edtmin) pef = edtmin;

        zkbc=z[kbcon(i)-1]*3.281e-3;
        prezk = 0.02;
        if (zkbc > 3.) prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc *(- 1.2569798e-2+zkbc*(4.2772e-4-zkbc*5.44e-6))));
        if (zkbc > 25) prezk = 2.4;
        pefb = 1. / (1. + prezk);
        if (pefb > edtmax) pefb = edtmax;
        if (pefb < edtmin) pefb = edtmin;
        edt = 1. - .5 * (pefb + pef);
        einc = .2 * edt;
        for(int k = 0; k < maxens2; k++) edtc[k] = edt + static_cast<real>(k - 1) * einc;
    }

    if (ierr == 0){
        for(int k = 0; k < maxens2; k++){
            edtc[k] = -edtc[k] * pwav / pwev;
            if (edtc[k] > edtmax)edtc[k] = edtmax;
            if (edtc[k] < edtmin)edtc[k] = edtmin;

        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_up_tracer(int &ierr, real &tcrit,const real2d &t, const real1d &pre,real *z_cup,const real2d &p,const real3d &tracer,real *tre_cup,real *tr_up,
            real *tr_pw,real *tr_c,const real2d &hstary,real *cupclw,real *clw_all,const int1d &kbcon,const int1d &ktop,real *cd,real &mentr_rate,real *zu,
        const int1d &k22,int &num_gas,char *chem_opt,bool &cu_scav_on,int &conv_tr_aqchem,const real3d &chem,const real2d &aq_aer_ratio,
            const real2d &aq_aer_ratio_comp_workspace, int workspace_col,
            int &i, int &kts, int &kte, int &ipr, int &jpr, int &j, int &npr, int &num_chem,const char* name,const int1d &ncomp_aer,
            int &mastercompindx_cl_aer, int &mastercompindx_na_aer, int &mastercompindx_oc_aer, int &mastercompindx_bc_aer, int &mastercompindx_so4_aer, 
            int &mastercompindx_oin_aer, int &mastercompindx_dust_aer, int &mastercompindx_nh4_aer, int &mastercompindx_msa_aer,
            int &mastercompindx_no3_aer, int &mastercompindx_ca_aer, int &mastercompindx_co3_aer,
            const int1d &nsize_aer,const int4d &massptr_aer,const real2d &dens_aer,
            int &p_so2,int &p_hno3,int &p_hno4,int &p_dms,int &p_n2o5,int &p_nh3,int &p_h2o2,int &p_o3,int &p_ora1,
            int &p_op1,int &p_paa,int &p_h2so4,int &p_so4aj,int &p_nh4aj ,int &p_no3aj){

    int static ngas = 12;     //! number of gas-phase species for aqchem
    int static naer = 36;     //! number of aerosol species for aqchem
    int static nliqs = 41;    //! number of liquid-phase species in aqchem
        //!...pointers for the aqchem array gas
    int static lso2    =  1;  //! sulfur dioxide
    int static lhno3   =  2;  //! nitric acid
    int static ln2o5   =  3;  //! dinitrogen pentoxide
    int static lco2    =  4;  //! carbon dioxide
    int static lnh3    =  5;  //! ammonia
    int static lh2o2   =  6;  //! hydrogen perioxide
    int static lo3     =  7;  //! ozone
    int static lfoa    =  8;  //! formic acid
    int static lmhp    =  9;  //! methyl hydrogen peroxide
    int static lpaa    = 10;  //! peroxyacidic acid
    int static lh2so4  = 11;  //! sulfuric acid
    int static lhcl    = 12;  //! hydrogen chloride
        //!...pointers for the aqchem array aerosol
    int static lso4akn  =  1;  //! aitken-mode sulfate
    int static lso4acc  =  2;  //! accumulation-mode sulfate
    int static lso4cor  =  3;  //! coarse-mode sulfate
    int static lnh4akn  =  4;  //! aitken-mode ammonium
    int static lnh4acc  =  5;  //! accumulation-mode ammonium
    int static lno3akn  =  6;  //! aitken-mode nitrate
    int static lno3acc  =  7;  //! accumulation-mode nitrate
    int static lno3cor  =  8;  //! coarse-mode nitrate
    int static lorgaakn =  9;  //! aitken-mode anthropogenic soa
    int static lorgaacc = 10;  //! accumulation-mode anthropogenic soa
    int static lorgpakn = 11;  //! aitken-mode primary organic aerosol
    int static lorgpacc = 12;  //! accumulation-mode primary organic aerosol
    int static lorgbakn = 13;  //! aitken-mode biogenic soa
    int static lorgbacc = 14;  //! accumulation-mode biogenic soa
    int static lecakn   = 15;  //! aitken-mode elemental carbon
    int static lecacc   = 16;  //! accumulation-mode elemental carbon
    int static lpriakn  = 17;  //! aitken-mode primary aerosol
    int static lpriacc  = 18;  //! accumulation-mode primary aerosol
    int static lpricor  = 19;  //! coarse-mode primary aerosol
    int static lnaakn   = 20;  //! aitken-mode sodium
    int static lnaacc   = 21;  //! accumulation-mode sodium
    int static lnacor   = 22;  //! coarse-mode sodium
    int static lclakn   = 23;  //! aitken-mode chloride ion
    int static lclacc   = 24;  //! accumulation-mode chloride ion
    int static lclcor   = 25;  //! coarse-mode chloride ion
    int static lnumakn  = 26;  //! aitken-mode number
    int static lnumacc  = 27;  //! accumulation-mode number
    int static lnumcor  = 28;  //! coarse-mode number
    int static lsrfakn  = 29;  //! aitken-mode surface area
    int static lsrfacc  = 30;  //! accumulation-mode surface area
    int static lnacl    = 31;  //! sodium chloride aerosol for ae3 only {depreciated in ae4}
    int static lcaco3   = 32;  //! calcium carbonate aerosol (place holder)
    int static lmgco3   = 33;  //! magnesium carbonate aerosol (place holder)
    int static la3fe    = 34;  //! iron aerosol (place holder)
    int static lb2mn    = 35;  //! manganese aerosol (place holder)
    int static lk       = 36;  //! potassium aerosol (cl- tracked separately) (place holder)

    //!  local variables in this routine
    real conc_equi, conc_mxr, partialp, taucld;
    real trcc, trch, dh, c0, dz, radius, airm, dens;
    int iall, k, iwd, nv;
    int iaer, igas;
    int icomp, ichem, itype, isize, iphase;
    //! I/O for AQCHEM:
    real precip; //! Precipitation rate (mm/h)
    real gas[ngas];
    real gaswdep[ngas];
    real aerosol[naer];
    real aerwdep[naer];
    real liquid[nliqs];
     
    real aq_gas_ratio;
    real hpwdep;
    real alfa0, alfa2, alfa3; //! Aerosol scavenging coefficients for Aitken mode
    //! jwfeng: Total aerosol mass inside grid box.
    real total_ratio, total_vol, temp_vol;
    //character(len=300) message
    
    iall = 0;
    c0 = 0.8;
    for(int k = 0;k < kte;k++) aq_aer_ratio(k,i) = 0.0;
    iwd = 0;

    if(std::strcmp(name, "shallow") == 0) c0 = 0;
    for(int k = 0;k < kte;k++){
        for(int nv = 0; nv <num_chem; nv++){
            tr_pw[nv*kte+k]=0.;
            if(ierr == 0) tr_up[nv*kte+k]=tre_cup[nv*kte+k];
            tr_c[nv*kte+k]=0.;
        }
    }

    for(int nv = 0; nv <num_chem; nv++){
        if(ierr==0.){
            for(int k = k22(i)-1; k < kbcon(i)-1; k++){
                tr_up[nv*kte+k]=tre_cup[nv*kte+k22(i)-1];
            }
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!              Diff from wrfchem, maybe cause error. zhxia 10/2024                !!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (cu_scav_on){
        itype = 0;
        iphase = 0;
        for(int nc = 0; nc < ncomp_aer(itype); nc++){
            aq_aer_ratio_comp_workspace(nc,workspace_col) = 0.0;
        }
        if (mastercompindx_cl_aer   >= 1 && mastercompindx_cl_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_cl_aer-1,workspace_col)    = 1.0;
        if (mastercompindx_na_aer   >= 1 && mastercompindx_na_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_na_aer-1,workspace_col)    = 1.0;
        if (mastercompindx_oc_aer   >= 1 && mastercompindx_oc_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_oc_aer-1,workspace_col)    = 0.5;
        if (mastercompindx_bc_aer   >= 1 && mastercompindx_bc_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_bc_aer-1,workspace_col)    = 0.1;
        if (mastercompindx_so4_aer  >= 1 && mastercompindx_so4_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_so4_aer-1,workspace_col)   = 1.0;
        if (mastercompindx_no3_aer  >= 1 && mastercompindx_no3_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_no3_aer-1,workspace_col)   = 1.0;
        if (mastercompindx_oin_aer  >= 1 && mastercompindx_oin_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_oin_aer-1,workspace_col)   = 0.1;
        if (mastercompindx_dust_aer >= 1 && mastercompindx_dust_aer <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_dust_aer-1,workspace_col)  = 0.1;
        if (mastercompindx_ca_aer   >= 1 && mastercompindx_ca_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_ca_aer-1,workspace_col)    = 0.1;
        if (mastercompindx_co3_aer  >= 1 && mastercompindx_co3_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_co3_aer-1,workspace_col)   = 0.1;
        if (mastercompindx_nh4_aer  >= 1 && mastercompindx_nh4_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_nh4_aer-1,workspace_col)   = 1.0;
        if (mastercompindx_msa_aer  >= 1 && mastercompindx_msa_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_msa_aer-1,workspace_col)   = 0.1;
        for(int k = 0;k < kte;k++){
            total_vol = 0.0;
            total_ratio = 0.0;
            for(int isize = 0; isize < nsize_aer(itype); isize++){
                for(int icomp = 0; icomp < ncomp_aer(itype); icomp++){
                    ichem = massptr_aer(icomp,isize,itype,iphase);
                    temp_vol = std::max(chem(ichem - 1, k, i), 0.0) / dens_aer(icomp, itype);
                    total_vol = total_vol + temp_vol;
                    total_ratio = total_ratio + temp_vol * aq_aer_ratio_comp_workspace(icomp,workspace_col);
                }
            }
            aq_aer_ratio(k,i) = total_ratio / total_vol;
            // if (aq_aer_ratio(k,i)<0.0 || aq_aer_ratio(k,i)>1.0){
            //     call chemistry_error_fatal("Wrong aq_aer_ratio"); }
        }
    }

    for(int k = 1;k < kte-1;k++){
        if (ierr != 0) continue;
        if (k < kbcon(i)-1) continue;
        if (k > ktop(i)) continue;
        dz = z_cup[k] - z_cup[k - 1];
        if (cupclw[k] <= 0.) cupclw[k] = 0.;
        if (clw_all[k] <= 0.) clw_all[k] = 0.;

        for(int nv = 0; nv < num_gas; nv++){
            tr_up[nv*kte+k]=(tr_up[nv*kte+k-1]*(1.-.5*cd[k]*dz)+mentr_rate*
                    dz*tracer(nv,k-1,i))/(1.+mentr_rate*dz-.5*cd[k]*dz);
            tr_up[nv*kte+k]=std::max(1.e-16,tr_up[nv*kte+k]);
        }

        if (conv_tr_aqchem == 1){
            dens = 0.1*p(k,i)/t(k,i)*mwdry/8.314472;
            airm = 1000.0*dens*dz/mwdry;
            for(int ng = 0; ng < ngas; ng++){
                gaswdep[ng] = 0.0;
            }
            for(int na = 0; na < naer; na++){
                aerwdep[na] = 0.0;
            }
            hpwdep = 0.0;

            precip = 0.0; //! mm/hr
            alfa0 = 0.0;
            alfa2 = 0.0;
            alfa3 = 0.0;

            for(int g = 0; g < ngas; g++){
                gas[g] = 0.0;
            }
            gas[lco2-1]   = 380.0e-6;
            gas[lso2-1]   = tr_up[(p_so2-1)*kte+k]*1.0e-6;
            gas[lhno3-1]  = tr_up[(p_hno3-1)*kte+k]*1.0e-6;
            gas[ln2o5-1]  = tr_up[(p_n2o5-1)*kte+k]*1.0e-6;
            gas[lnh3-1]   = tr_up[(p_nh3-1)*kte+k]*1.0e-6;
            gas[lh2o2-1]  = tr_up[(p_h2o2-1)*kte+k]*1.0e-6;
            gas[lo3-1]    = tr_up[(p_o3-1)*kte+k]*1.0e-6;
            gas[lfoa-1]   = tr_up[(p_ora1-1)*kte+k]*1.0e-6;
            gas[lmhp-1]   = tr_up[(p_op1-1)*kte+k]*1.0e-6;
            gas[lpaa-1]   = tr_up[(p_paa-1)*kte+k]*1.0e-6;
            gas[lh2so4-1] = tr_up[(p_h2so4-1)*kte+k]*1.0e-6;

            for(int na = 0; na < naer; na++){
                aerosol[na] = 0.0;
            }

            aerosol[lso4acc-1] = tr_up[(p_so4aj-1)*kte+k]*1.0e-9*mwdry/mwso4;
            aerosol[lnh4acc-1] = tr_up[(p_nh4aj-1)*kte+k]*1.0e-9*mwdry/mwnh4;
            aerosol[lno3acc-1] = tr_up[(p_no3aj-1)*kte+k]*1.0e-9*mwdry/mwno3;
            taucld = 1800.0;

            if (clw_all[k] * dens > qcldwtr_cutoff){
                real temp = t(k,i);
                real pres_pa = p(k,i)*100.;
                real wcavg = clw_all[k]*dens;
                real wtavg = clw_all[k]*dens;
                c_ctrans_aqchem_bridge(temp, pres_pa, taucld, precip,
                        wcavg, wtavg, airm, alfa0, alfa2, alfa3, gas,
                        aerosol, liquid, gaswdep, aerwdep, hpwdep);
            }

            tr_up[(p_so2-1)*kte+k]  =  gas[lso2-1]*1.0e6;
            tr_up[(p_hno3-1)*kte+k] =  gas[lhno3-1]*1.0e6;
            tr_up[(p_n2o5-1)*kte+k] =  gas[ln2o5-1]*1.0e6;
            tr_up[(p_nh3-1)*kte+k]  =  gas[lnh3-1]*1.0e6;
            tr_up[(p_h2o2-1)*kte+k] =  gas[lh2o2-1]*1.0e6;
            tr_up[(p_o3-1)*kte+k]   =  gas[lo3-1]*1.0e6;
            tr_up[(p_ora1-1)*kte+k] =  gas[lfoa-1]*1.0e6;
            tr_up[(p_op1-1)*kte+k]  =  gas[lmhp-1]*1.0e6;
            tr_up[(p_paa-1)*kte+k]  =  gas[lpaa-1]*1.0e6;
            tr_up[(p_h2so4-1)*kte+k] =  gas[lh2so4-1]*1.0e6;

            tr_up[(p_so4aj-1)*kte+k] = aerosol[lso4acc-1]*1.0e9*mwso4/mwdry;
            tr_up[(p_nh4aj-1)*kte+k] = aerosol[lnh4acc-1]*1.0e9*mwnh4/mwdry;
            tr_up[(p_no3aj-1)*kte+k] = aerosol[lno3acc-1]*1.0e9*mwno3/mwdry;
        }

        if (cu_scav_on){
            for(int nv = 0; nv < num_gas; nv++){
                aq_gas_ratio = 0.0;
                if (nv == p_so2-1) aq_gas_ratio = 1.0;
                if (nv == p_h2so4-1) aq_gas_ratio = 1.0;
                if (nv == p_nh3-1) aq_gas_ratio = 1.0;
                if (nv == p_hno3-1) aq_gas_ratio = 1.0;
                tr_c[nv*kte+k] = aq_gas_ratio*tr_up[nv*kte+k];
                trch = tr_up[nv*kte+k]-tr_c[nv*kte+k];
                tr_pw[nv*kte+k] = c0*dz*tr_c[nv*kte+k]*zu[k];
                trcc = tr_c[nv*kte+k]*(1-c0);
                tr_up[nv*kte+k] = trcc+trch;
            }
            for(int nv = num_gas; nv < num_chem; nv++){
                if (aq_aer_ratio(k,i) > 0.0){
                    tr_c[nv*kte+k] = aq_aer_ratio(k,i)*tr_up[nv*kte+k];
                    trch = tr_up[nv*kte+k]-tr_c[nv*kte+k];
                    tr_pw[nv*kte+k] = c0*dz*tr_c[nv*kte+k]*zu[k];
                    trcc = tr_c[nv*kte+k]*(1-c0);
                    tr_up[nv*kte+k] = trcc+trch;
                }
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_up_tracer_l(int &ierr, real &tcrit,const real2d &t, const real1d &pre,real *z_cup,const real2d &p,const real3d &tracer,real *tre_cup,real *tr_up,
            real *tr_pw,real *tr_c,const real2d &hstary,real *cupclw,real *clw_all,const int1d &kbcon,const int1d &ktop,real *cd,real &mentr_rate,real *zu,
            int &k22,int &num_gas,char *chem_opt,bool &cu_scav_on,int &conv_tr_aqchem,const real3d &chem,const real2d &aq_aer_ratio,
            const real2d &aq_aer_ratio_comp_workspace, int workspace_col,
            int &i, int &kts, int &kte, int &ipr, int &jpr, int &j, int &npr, int &num_chem,const char* name,const int1d &ncomp_aer,
            int &mastercompindx_cl_aer, int &mastercompindx_na_aer, int &mastercompindx_oc_aer, int &mastercompindx_bc_aer, int &mastercompindx_so4_aer, 
            int &mastercompindx_oin_aer, int &mastercompindx_dust_aer, int &mastercompindx_nh4_aer, int &mastercompindx_msa_aer,
            int &mastercompindx_no3_aer, int &mastercompindx_ca_aer, int &mastercompindx_co3_aer,
            const int1d &nsize_aer,const int4d &massptr_aer,const real2d &dens_aer,
            int &p_so2,int &p_hno3,int &p_hno4,int &p_dms,int &p_n2o5,int &p_nh3,int &p_h2o2,int &p_o3,int &p_ora1,
            int &p_op1,int &p_paa,int &p_h2so4,int &p_so4aj,int &p_nh4aj ,int &p_no3aj){

    int static ngas = 12;     //! number of gas-phase species for aqchem
    int static naer = 36;     //! number of aerosol species for aqchem
    int static nliqs = 41;    //! number of liquid-phase species in aqchem
        //!...pointers for the aqchem array gas
    int static lso2    =  1;  //! sulfur dioxide
    int static lhno3   =  2;  //! nitric acid
    int static ln2o5   =  3;  //! dinitrogen pentoxide
    int static lco2    =  4;  //! carbon dioxide
    int static lnh3    =  5;  //! ammonia
    int static lh2o2   =  6;  //! hydrogen perioxide
    int static lo3     =  7;  //! ozone
    int static lfoa    =  8;  //! formic acid
    int static lmhp    =  9;  //! methyl hydrogen peroxide
    int static lpaa    = 10;  //! peroxyacidic acid
    int static lh2so4  = 11;  //! sulfuric acid
    int static lhcl    = 12;  //! hydrogen chloride
        //!...pointers for the aqchem array aerosol
    int static lso4akn  =  1;  //! aitken-mode sulfate
    int static lso4acc  =  2;  //! accumulation-mode sulfate
    int static lso4cor  =  3;  //! coarse-mode sulfate
    int static lnh4akn  =  4;  //! aitken-mode ammonium
    int static lnh4acc  =  5;  //! accumulation-mode ammonium
    int static lno3akn  =  6;  //! aitken-mode nitrate
    int static lno3acc  =  7;  //! accumulation-mode nitrate
    int static lno3cor  =  8;  //! coarse-mode nitrate
    int static lorgaakn =  9;  //! aitken-mode anthropogenic soa
    int static lorgaacc = 10;  //! accumulation-mode anthropogenic soa
    int static lorgpakn = 11;  //! aitken-mode primary organic aerosol
    int static lorgpacc = 12;  //! accumulation-mode primary organic aerosol
    int static lorgbakn = 13;  //! aitken-mode biogenic soa
    int static lorgbacc = 14;  //! accumulation-mode biogenic soa
    int static lecakn   = 15;  //! aitken-mode elemental carbon
    int static lecacc   = 16;  //! accumulation-mode elemental carbon
    int static lpriakn  = 17;  //! aitken-mode primary aerosol
    int static lpriacc  = 18;  //! accumulation-mode primary aerosol
    int static lpricor  = 19;  //! coarse-mode primary aerosol
    int static lnaakn   = 20;  //! aitken-mode sodium
    int static lnaacc   = 21;  //! accumulation-mode sodium
    int static lnacor   = 22;  //! coarse-mode sodium
    int static lclakn   = 23;  //! aitken-mode chloride ion
    int static lclacc   = 24;  //! accumulation-mode chloride ion
    int static lclcor   = 25;  //! coarse-mode chloride ion
    int static lnumakn  = 26;  //! aitken-mode number
    int static lnumacc  = 27;  //! accumulation-mode number
    int static lnumcor  = 28;  //! coarse-mode number
    int static lsrfakn  = 29;  //! aitken-mode surface area
    int static lsrfacc  = 30;  //! accumulation-mode surface area
    int static lnacl    = 31;  //! sodium chloride aerosol for ae3 only {depreciated in ae4}
    int static lcaco3   = 32;  //! calcium carbonate aerosol (place holder)
    int static lmgco3   = 33;  //! magnesium carbonate aerosol (place holder)
    int static la3fe    = 34;  //! iron aerosol (place holder)
    int static lb2mn    = 35;  //! manganese aerosol (place holder)
    int static lk       = 36;  //! potassium aerosol (cl- tracked separately) (place holder)

    //!  local variables in this routine
    real conc_equi, conc_mxr, partialp, taucld;
    real trcc, trch, dh, c0, dz, radius, airm, dens;
    int iall, k, iwd, nv;
    int iaer, igas;
    int icomp, ichem, itype, isize, iphase;
    //! I/O for AQCHEM:
    real precip; //! Precipitation rate (mm/h)
    real gas[ngas];
    real gaswdep[ngas];
    real aerosol[naer];
    real aerwdep[naer];
    real liquid[nliqs];
     
    real aq_gas_ratio;
    real hpwdep;
    real alfa0, alfa2, alfa3; //! Aerosol scavenging coefficients for Aitken mode
    //! jwfeng: Total aerosol mass inside grid box.
    real total_ratio, total_vol, temp_vol;
    //character(len=300) message
    
    iall = 0;
    c0 = 0.8;
    for(int k = 0;k < kte;k++) aq_aer_ratio(k,i) = 0.0;
    iwd = 0;

    if(std::strcmp(name, "shallow") == 0) c0 = 0;
    for(int k = 0;k < kte;k++){
        for(int nv = 0; nv <num_chem; nv++){
            tr_pw[nv*kte+k]=0.;
            if(ierr == 0) tr_up[nv*kte+k]=tre_cup[nv*kte+k];
            tr_c[nv*kte+k]=0.;
        }
    }

    for(int nv = 0; nv <num_chem; nv++){
        if(ierr==0.){
            for(int k = k22-1; k < kbcon(i)-1; k++){
                tr_up[nv*kte+k]=tre_cup[nv*kte+k22-1];
            }
        }
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!              Diff from wrfchem, maybe cause error. zhxia 10/2024                !!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (cu_scav_on){
        itype = 0;
        iphase = 0;
        for(int nc = 0; nc < ncomp_aer(itype); nc++){
            aq_aer_ratio_comp_workspace(nc,workspace_col) = 0.0;
        }
        if (mastercompindx_cl_aer   >= 1 && mastercompindx_cl_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_cl_aer-1,workspace_col)    = 1.0;
        if (mastercompindx_na_aer   >= 1 && mastercompindx_na_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_na_aer-1,workspace_col)    = 1.0;
        if (mastercompindx_oc_aer   >= 1 && mastercompindx_oc_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_oc_aer-1,workspace_col)    = 0.5;
        if (mastercompindx_bc_aer   >= 1 && mastercompindx_bc_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_bc_aer-1,workspace_col)    = 0.1;
        if (mastercompindx_so4_aer  >= 1 && mastercompindx_so4_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_so4_aer-1,workspace_col)   = 1.0;
        if (mastercompindx_no3_aer  >= 1 && mastercompindx_no3_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_no3_aer-1,workspace_col)   = 1.0;
        if (mastercompindx_oin_aer  >= 1 && mastercompindx_oin_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_oin_aer-1,workspace_col)   = 0.1;
        if (mastercompindx_dust_aer >= 1 && mastercompindx_dust_aer <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_dust_aer-1,workspace_col)  = 0.1;
        if (mastercompindx_ca_aer   >= 1 && mastercompindx_ca_aer   <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_ca_aer-1,workspace_col)    = 0.1;
        if (mastercompindx_co3_aer  >= 1 && mastercompindx_co3_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_co3_aer-1,workspace_col)   = 0.1;
        if (mastercompindx_nh4_aer  >= 1 && mastercompindx_nh4_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_nh4_aer-1,workspace_col)   = 1.0;
        if (mastercompindx_msa_aer  >= 1 && mastercompindx_msa_aer  <= ncomp_aer(itype)) aq_aer_ratio_comp_workspace(mastercompindx_msa_aer-1,workspace_col)   = 0.1;
        for(int k = 0;k < kte;k++){
            total_vol = 0.0;
            total_ratio = 0.0;
            for(int isize = 0; isize < nsize_aer(itype); isize++){
                for(int icomp = 0; icomp < ncomp_aer(itype); icomp++){
                    ichem = massptr_aer(icomp,isize,itype,iphase);
                    temp_vol = std::max(chem(ichem - 1, k, i), 0.0) / dens_aer(icomp, itype);
                    total_vol = total_vol + temp_vol;
                    total_ratio = total_ratio + temp_vol * aq_aer_ratio_comp_workspace(icomp,workspace_col);
                }
            }
            aq_aer_ratio(k,i) = total_ratio / total_vol;
            // if (aq_aer_ratio(k,i)<0.0 || aq_aer_ratio(k,i)>1.0){
            //     call chemistry_error_fatal("Wrong aq_aer_ratio"); }
        }
    }

    for(int k = 1;k < kte-1;k++){
        if (ierr != 0) continue;
        if (k < kbcon(i)-1) continue;
        if (k > ktop(i)) continue;
        dz = z_cup[k] - z_cup[k - 1];
        if (cupclw[k] <= 0.) cupclw[k] = 0.;
        if (clw_all[k] <= 0.) clw_all[k] = 0.;

        for(int nv = 0; nv < num_gas; nv++){
            tr_up[nv*kte+k]=(tr_up[nv*kte+k-1]*(1.-.5*cd[k]*dz)+mentr_rate*
                    dz*tracer(nv,k-1,i))/(1.+mentr_rate*dz-.5*cd[k]*dz);
            tr_up[nv*kte+k]=std::max(1.e-16,tr_up[nv*kte+k]);
        }

        if (conv_tr_aqchem == 1){
            dens = 0.1*p(k,i)/t(k,i)*mwdry/8.314472;
            airm = 1000.0*dens*dz/mwdry;
            for(int ng = 0; ng < ngas; ng++){
                gaswdep[ng] = 0.0;
            }
            for(int na = 0; na < naer; na++){
                aerwdep[na] = 0.0;
            }
            hpwdep = 0.0;

            precip = 0.0; //! mm/hr
            alfa0 = 0.0;
            alfa2 = 0.0;
            alfa3 = 0.0;

            for(int g = 0; g < ngas; g++){
                gas[g] = 0.0;
            }
            gas[lco2-1]   = 380.0e-6;
            gas[lso2-1]   = tr_up[(p_so2-1)*kte+k]*1.0e-6;
            gas[lhno3-1]  = tr_up[(p_hno3-1)*kte+k]*1.0e-6;
            gas[ln2o5-1]  = tr_up[(p_n2o5-1)*kte+k]*1.0e-6;
            gas[lnh3-1]   = tr_up[(p_nh3-1)*kte+k]*1.0e-6;
            gas[lh2o2-1]  = tr_up[(p_h2o2-1)*kte+k]*1.0e-6;
            gas[lo3-1]    = tr_up[(p_o3-1)*kte+k]*1.0e-6;
            gas[lfoa-1]   = tr_up[(p_ora1-1)*kte+k]*1.0e-6;
            gas[lmhp-1]   = tr_up[(p_op1-1)*kte+k]*1.0e-6;
            gas[lpaa-1]   = tr_up[(p_paa-1)*kte+k]*1.0e-6;
            gas[lh2so4-1] = tr_up[(p_h2so4-1)*kte+k]*1.0e-6;

            for(int na = 0; na < naer; na++){
                aerosol[na] = 0.0;
            }

            aerosol[lso4acc-1] = tr_up[(p_so4aj-1)*kte+k]*1.0e-9*mwdry/mwso4;
            aerosol[lnh4acc-1] = tr_up[(p_nh4aj-1)*kte+k]*1.0e-9*mwdry/mwnh4;
            aerosol[lno3acc-1] = tr_up[(p_no3aj-1)*kte+k]*1.0e-9*mwdry/mwno3;
            taucld = 1800.0;

            if (clw_all[k] * dens > qcldwtr_cutoff){
                real temp = t(k,i);
                real pres_pa = p(k,i)*100.;
                real wcavg = clw_all[k]*dens;
                real wtavg = clw_all[k]*dens;
                c_ctrans_aqchem_bridge(temp, pres_pa, taucld, precip,
                        wcavg, wtavg, airm, alfa0, alfa2, alfa3, gas,
                        aerosol, liquid, gaswdep, aerwdep, hpwdep);
            }

            tr_up[(p_so2-1)*kte+k]  =  gas[lso2-1]*1.0e6;
            tr_up[(p_hno3-1)*kte+k] =  gas[lhno3-1]*1.0e6;
            tr_up[(p_n2o5-1)*kte+k] =  gas[ln2o5-1]*1.0e6;
            tr_up[(p_nh3-1)*kte+k]  =  gas[lnh3-1]*1.0e6;
            tr_up[(p_h2o2-1)*kte+k] =  gas[lh2o2-1]*1.0e6;
            tr_up[(p_o3-1)*kte+k]   =  gas[lo3-1]*1.0e6;
            tr_up[(p_ora1-1)*kte+k] =  gas[lfoa-1]*1.0e6;
            tr_up[(p_op1-1)*kte+k]  =  gas[lmhp-1]*1.0e6;
            tr_up[(p_paa-1)*kte+k]  =  gas[lpaa-1]*1.0e6;
            tr_up[(p_h2so4-1)*kte+k] =  gas[lh2so4-1]*1.0e6;

            tr_up[(p_so4aj-1)*kte+k] = aerosol[lso4acc-1]*1.0e9*mwso4/mwdry;
            tr_up[(p_nh4aj-1)*kte+k] = aerosol[lnh4acc-1]*1.0e9*mwnh4/mwdry;
            tr_up[(p_no3aj-1)*kte+k] = aerosol[lno3acc-1]*1.0e9*mwno3/mwdry;
        }

        if (cu_scav_on){
            for(int nv = 0; nv < num_gas; nv++){
                aq_gas_ratio = 0.0;
                if (nv == p_so2-1) aq_gas_ratio = 1.0;
                if (nv == p_h2so4-1) aq_gas_ratio = 1.0;
                if (nv == p_nh3-1) aq_gas_ratio = 1.0;
                if (nv == p_hno3-1) aq_gas_ratio = 1.0;
                tr_c[nv*kte+k] = aq_gas_ratio*tr_up[nv*kte+k];
                trch = tr_up[nv*kte+k]-tr_c[nv*kte+k];
                tr_pw[nv*kte+k] = c0*dz*tr_c[nv*kte+k]*zu[k];
                trcc = tr_c[nv*kte+k]*(1-c0);
                tr_up[nv*kte+k] = trcc+trch;
            }
            for(int nv = num_gas; nv < num_chem; nv++){
                if (aq_aer_ratio(k,i) > 0.0){
                    tr_c[nv*kte+k] = aq_aer_ratio(k,i)*tr_up[nv*kte+k];
                    trch = tr_up[nv*kte+k]-tr_c[nv*kte+k];
                    tr_pw[nv*kte+k] = c0*dz*tr_c[nv*kte+k]*zu[k];
                    trcc = tr_c[nv*kte+k]*(1-c0);
                    tr_up[nv*kte+k] = trcc+trch;
                }
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_dellas_tr(int &ierr,real *z_cup,real *p_cup,real *tr_dd,real &edt,real *zd,real *cdd,
                    const real3d &tracer, const real3d &tracert,int &j,real &mentrd_rate,real *zu,real &g,const real1d &xmb,
                           real *cd,real *tr_up,const int1d &ktop,const int1d &k22, const int1d &kbcon, real &mentr_rate,int &jmin,
                           real *tre_cup, int &kdet, const int1d &kpbl,
                           int &ipr, int &jpr, int &npr,const char* name, int &num_chem,
                           int &i, int &kts, int &kte){

    //loacal varibales
    int k, nv;
    real detdo1, detdo2, entdo, dp, dz, subin, detdo, entup, detup, subdown, entdoj, entupk, detupk, totmas;
    int itf, ktf, kstart;
    //character(len=300) message
    
    kstart = kts;
    if (std::strcmp(name, "shallow") == 0) kstart = 0;
    for(int k = kstart;k < kte;k++){
        for(int nv = 0; nv <num_chem; nv++){
            tracert(nv, k, i) = 0.;
        }
    }

    for(int k = 1;k < kte-1;k++){
        if (ierr != 0) continue;
        if (k > ktop(i)-1) continue;

        //!--- specify detrainment of downdraft, has to be consistent
        //!--- with zd calculations in soundd.
        dz = z_cup[k + 1] - z_cup[k];
        detdo = edt * cdd[k] * dz * zd[k + 1];
        entdo = edt * mentrd_rate * dz * zd[k + 1];
        subin = zu[k + 1] - zd[k + 1] * edt;
        entup = 0.;
        detup = 0.;
        if (k >= kbcon(i)-1 && k < ktop(i)-1){
            entup = mentr_rate * dz * zu[k];
            detup = cd[k + 1] * dz * zu[k];
        }
        subdown = (zu[k] - zd[k] * edt);
        entdoj = 0.;
        entupk = 0.;
        detupk = 0.;

        if (k == jmin - 1) entdoj = edt * zd[k];

        if (k == k22(i) - 2) entupk = zu[kpbl(i)-1];
        if (k > kdet - 1) detdo = 0.;

        if (k == ktop(i) - 1){
            detupk = zu[ktop(i)-1];
            subin = 0.;
        }

        if (k < kbcon(i)-1) detup = 0.;

        //!--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
        totmas = subin - subdown + detup - entup - entdo + detdo - entupk - entdoj + detupk;
        // if ( i .eq. ipr) then
        //    write(message, "(a,x,I3,x,3F9.3)") 'k, totmas, sui, sud = ', k, totmas, subin, subdown
        //    call chemistry_message(trim(message))
        // end if

        // if (abs(totmas) .gt. 1.e-6) then
        //    write(message, "(a,3I3,x,F9.3,a)") '*********************',i,j,k,totmas,name
        //    call chemistry_message(trim(message))
        //    write(message, "(4I5)") kpbl(i), k22(i), kbcon(i), ktop(i)
        //    call chemistry_message(trim(message))
        //    call chemistry_error_fatal( 'cup_dellas_tr: TOTMAS > CRITICAL VALUE')
        // endif

        dp = 100. * (p_cup[k - 1] - p_cup[k]);
        for(int nv = 0; nv <num_chem; nv++){
            tracert(nv,k,i)=(subin*tracer(nv,k+1,i)
                -subdown*tracer(nv,k,i)
                +detup*.5*(tr_up[nv*kte+k+1]+tr_up[nv*kte+k])
                +detdo*.5*(tr_dd[nv*kte+k+1]+tr_dd[nv*kte+k])
                -entup*tracer(nv,k,i)
                -entdo*tracer(nv,k,i)
                -entupk*tre_cup[nv*kte+k22(i)-1]
                -entdoj*tre_cup[nv*kte+jmin-1]
                +detupk*tr_up[nv*kte+ktop(i)-1]
                )*g/dp*xmb(i);
        }
        // if(i.eq.ipr)then
        //    write(message, "(I3,x,4F9.3)") k, tracert(npr, k, i), subin * tre_cup(npr, k + 1, i), subdown * tre_cup(npr, k, i), &
        //                      detdo * .5 * (tr_dd(npr, k + 1, i) + tr_dd(npr, k, i))
        //    call chemistry_message(trim(message))
        //    write(message, "(I3,x,4F9.3)") k, detup * .5 * (tr_up(npr, k + 1, i) + tr_up(npr, k, i)), detupk * tr_up(npr, ktop(i), i), &
        //                      entup * tracer(npr, k, i), entdo * tracer(npr, k, i)
        //    call chemistry_message(trim(message))
        //    write(message, "(I3,x,3F9.3)") k, entupk * tre_cup(npr, k, i), detupk, tr_up(npr, ktop(i), i)
        //    call chemistry_message(trim(message))
        // endif
    }
}

KOKKOS_INLINE_FUNCTION void cup_dellas_tr_l(int &ierr,real *z_cup,real *p_cup,real *tr_dd,real &edt,real *zd,real *cdd,
                           const real3d &tracer, const real3d &tracert,int &j,real &mentrd_rate,real *zu,real &g,real &xmb,
                           real *cd,real *tr_up,const int1d &ktop,int &k22, const int1d &kbcon, real &mentr_rate,int &jmin,
                           real *tre_cup, int &kdet,int &kpbl,
                           int &ipr, int &jpr, int &npr,const char* name, int &num_chem,
                           int &i, int &kts, int &kte){

    //loacal varibales
    int k, nv;
    real detdo1, detdo2, entdo, dp, dz, subin, detdo, entup, detup, subdown, entdoj, entupk, detupk, totmas;
    int itf, ktf, kstart;
    //character(len=300) message
    
    kstart = kts;
    if (std::strcmp(name, "shallow") == 0) kstart = 0;
    for(int k = kstart;k < kte;k++){
        for(int nv = 0; nv <num_chem; nv++){
            tracert(nv, k, i) = 0.;
        }
    }

    for(int k = 1;k < kte-1;k++){
        if (ierr != 0) continue;
        if (k > ktop(i) -1) continue;

        //!--- specify detrainment of downdraft, has to be consistent
        //!--- with zd calculations in soundd.
        dz = z_cup[k + 1] - z_cup[k];
        detdo = edt * cdd[k] * dz * zd[k + 1];
        entdo = edt * mentrd_rate * dz * zd[k + 1];
        subin = zu[k + 1] - zd[k + 1] * edt;
        entup = 0.;
        detup = 0.;
        if (k >= kbcon(i)-1 && k < ktop(i)-1){
            entup = mentr_rate * dz * zu[k];
            detup = cd[k + 1] * dz * zu[k];
        }
        subdown = (zu[k] - zd[k] * edt);
        entdoj = 0.;
        entupk = 0.;
        detupk = 0.;

        if (k == jmin - 1) entdoj = edt * zd[k];

        if (k == k22 - 2) entupk = zu[kpbl-1];
        if (k > kdet - 1) detdo = 0.;

        if (k == ktop(i) - 1){
            detupk = zu[ktop(i)-1];
            subin = 0.;
        }

        if (k < kbcon(i)-1) detup = 0.;

        //!--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
        totmas = subin - subdown + detup - entup - entdo + detdo - entupk - entdoj + detupk;
        // if ( i .eq. ipr) then
        //    write(message, "(a,x,I3,x,3F9.3)") 'k, totmas, sui, sud = ', k, totmas, subin, subdown
        //    call chemistry_message(trim(message))
        // end if

        // if (abs(totmas) .gt. 1.e-6) then
        //    write(message, "(a,3I3,x,F9.3,a)") '*********************',i,j,k,totmas,name
        //    call chemistry_message(trim(message))
        //    write(message, "(4I5)") kpbl(i), k22(i), kbcon(i), ktop(i)
        //    call chemistry_message(trim(message))
        //    call chemistry_error_fatal( 'cup_dellas_tr: TOTMAS > CRITICAL VALUE')
        // endif

        dp = 100. * (p_cup[k - 1] - p_cup[k]);
        for(int nv = 0; nv <num_chem; nv++){
            tracert(nv,k,i)=(subin*tracer(nv,k+1,i)
                -subdown*tracer(nv,k,i)
                +detup*.5*(tr_up[nv*kte+k+1]+tr_up[nv*kte+k])
                +detdo*.5*(tr_dd[nv*kte+k+1]+tr_dd[nv*kte+k])
                -entup*tracer(nv,k,i)
                -entdo*tracer(nv,k,i)
                -entupk*tre_cup[nv*kte+k22-1]
                -entdoj*tre_cup[nv*kte+jmin-1]
                +detupk*tr_up[nv*kte+ktop(i)-1]
                )*g/dp*xmb;
        }
        // if(i.eq.ipr)then
        //    write(message, "(I3,x,4F9.3)") k, tracert(npr, k, i), subin * tre_cup(npr, k + 1, i), subdown * tre_cup(npr, k, i), &
        //                      detdo * .5 * (tr_dd(npr, k + 1, i) + tr_dd(npr, k, i))
        //    call chemistry_message(trim(message))
        //    write(message, "(I3,x,4F9.3)") k, detup * .5 * (tr_up(npr, k + 1, i) + tr_up(npr, k, i)), detupk * tr_up(npr, ktop(i), i), &
        //                      entup * tracer(npr, k, i), entdo * tracer(npr, k, i)
        //    call chemistry_message(trim(message))
        //    write(message, "(I3,x,3F9.3)") k, entupk * tre_cup(npr, k, i), detupk, tr_up(npr, ktop(i), i)
        //    call chemistry_message(trim(message))
        // endif
    }
}

KOKKOS_INLINE_FUNCTION void cup_dd_tracer(int &ierr,real *z_cup,real *qrcd, const real3d &tracer,real *tre_cup,real *tr_up,real *tr_dd,
         real *tr_pw,real *tr_pwd,int &jmin,real *cdd,real &entr,real *zd,real *pwdper,const real3d &wetdep,real &xmb,int &k22,const real2d &rho_phy,const real2d &dz8w,
         int &num_gas,int &num_chem,
         int &i, int &kts,int &kte){
    real mwdry   = 28.966;

    int iall;
    real dh, c0, dz, radius;
    real evaporate, condensate, net_removed, air_mass;

    for(int num = 0; num < num_chem; num++){
        for(int k = 0;k < kte;k++){
            tr_pwd[num*kte + k] = 0.0;
        }
    }
    for(int num = 0; num < num_chem; num++){
        for(int k = 0;k < kte;k++){
            wetdep(num,k,i) = 0.0;
        }
    }
    // ! Calculate wet deposition with re-evaporation (based on wet scavenging in cup_up_tracer);
    // ! assume no gas takeup by rain during fall for now
    if(ierr == 0){
        for(int nv = 0; nv < num_gas; nv++){
            for(int k = kte-1; k >= 0; k--){
                condensate = tr_pw[nv*kte+k];
                for(int kj = k; kj >= 0; kj--){
                    evaporate = condensate*pwdper[kj];
                    tr_pwd[nv*kte+kj] = tr_pwd[nv*kte+kj] + evaporate;
                    condensate = std::max(0.0,condensate-evaporate);
                }
            }
        }
    }
    // xyruan: calculate model-layer mixing-ratio tendency.
    if(ierr == 0){
        for(int k = kts - 1; k < kte; k++){
            air_mass = rho_phy(k,i) * dz8w(k,i);
            if (air_mass > 0.0){
                for(int nv = 0; nv < num_chem; nv++){
                    if (nv < num_gas){
                        net_removed = tr_pw[nv*kte+k] - tr_pwd[nv*kte + k];
                    }else{
                        net_removed = tr_pw[nv*kte+k];
                    }
                    wetdep(nv,k,i) = - net_removed * xmb / air_mass;
                }
            }
        }
    }
    //! In downdraft, do only transport of tracers
    if(ierr == 0){
        for(int nv = 0; nv < num_chem; nv++) tr_dd[nv*kte + jmin-1] = tre_cup[nv*kte + jmin-1];
        for(int ki = jmin-2; ki >=0; ki--){
            dz=z_cup[ki+1] - z_cup[ki];
            for(int nv = 0; nv < num_chem; nv++){
                tr_dd[nv*kte+ki]=(tr_dd[nv*kte+ki+1]*(1.-.5*cdd[ki]*dz)
                        +entr*dz*tracer(nv,ki,i))
                        /(1.+entr*dz-.5*cdd[ki]*dz);
            }
        }
    }

    if(ierr == 0){
        for(int k = 0;k < kte;k++){
            for(int nv = 0; nv < num_chem; nv++){
                tr_dd[nv*kte+k] = tr_dd[nv*kte+k] + tr_pwd[nv*kte+k];
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION void cup_dellabot_tr(int &ipr, int &jpr,real *tre_cup,int &ierr,real *z_cup,const real2d &p_cup,
                              real *tr_dd,real &edt,real *zd,real *cdd,const real3d &tracer, const real3d &tracert,int &j,
                              real &mentrd_rate,real *z, real &g,real &xmb,
                              int &num_chem,
                              int &i, int &kts, int &kte){
    real detdo, detdo1, detdo2, entdo, dp, dz, subin, totmas;
    int  npr = 24;
    //character(len=300) message
    for(int num = 0; num < num_chem; num++){
        for(int k = 0;k < kte;k++){
            tracert(num, k, i) = 0.;
        }
    }
    if (ierr != 0) return;
    dz = z_cup[1] - z_cup[0];
    dp = 100. * (p_cup(0,i) - p_cup(1,i));
    detdo1 = edt * zd[1] * cdd[0] * dz;
    detdo2 = edt * zd[0];
    entdo = edt * zd[1] * mentrd_rate * dz;
    subin = - edt * zd[1];
    detdo = detdo1 + detdo2 - entdo + subin;
    for(int nv = 0; nv < num_chem; nv++){
        tracert(nv,0,i)=(detdo1*.5*(tr_dd[nv*kte+0]+tr_dd[nv*kte+1])
                +detdo2*tr_dd[nv*kte+0]
                +subin*tre_cup[nv*kte+1]
                -entdo*tracer(nv,0,i))*g/dp*xmb;
    }
    // if (i .eq. ipr) then
    //     write(message, "(a,x,8F9.3)") 'in cup dellabot ', tracert(npr, 1, i) , &
    //     detdo1, detdo2, subin, entdo, tr_dd(npr, 1, i), tr_dd(npr, 2, i), tracer(npr, 1, i)
    //     call chemistry_message(trim(message))
    // end if
}
