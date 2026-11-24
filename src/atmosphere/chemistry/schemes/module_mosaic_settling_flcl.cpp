#include <flcl-cxx.hpp>
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <math.h>
#include <utility>

typedef double real;

typedef Kokkos::View<int**, Kokkos::LayoutLeft> int2d;

typedef Kokkos::View<real*, Kokkos::LayoutLeft> real1d;

typedef Kokkos::View<real**, Kokkos::LayoutLeft> real2d;

typedef Kokkos::View<real***, Kokkos::LayoutLeft> real3d;

typedef Kokkos::View<real****, Kokkos::LayoutLeft> real4d;

constexpr int CHEM_WORKSPACE_CHUNK_SIZE = 64;

//===========================================
// @module module_data_binaer_asect
//===========================================
constexpr int maxd_atype = 1;
constexpr int maxd_asize = 10;
constexpr int maxd_acomp = 120;
constexpr int maxd_aphase = 2;

constexpr real hygro_na_aer = 1.16;
constexpr real hygro_dust_aer = 0.10;

// extern "C" int ncomp_aer[maxd_atype];
// extern "C" int massptr_aer[maxd_aphase][maxd_atype][maxd_asize][maxd_acomp];
// extern "C" real dens_aer[maxd_atype][maxd_acomp];
// extern "C" real hygro_aer[maxd_atype][maxd_acomp];
// extern "C" real dcen_sect[maxd_atype][maxd_asize];

// extern "C" int num_moist;  // @module mpas_atmphys_vars
// extern "C" int num_chem;   // @module mpas_atmchem_vars
// extern "C" int nbin_a;     // @module module_data_mosaic_therm

constexpr real dyn_visc = 1.5E-5;  // @module module_data_gocart_dust

struct mosaic_settling_workspace {
  real2d tmp, airmas, p_mid, delz, rh, rho_air;
  real3d v_grav_1d, avghygro, avgden, avgreff;
  real4d totmass0, totmass, bstl_mass;
  real2d dt_settl;
  int2d ndt_settl;
  real4d tc1;
  real3d vd_delz_dt;

  mosaic_settling_workspace(int nlev, int ncol, int nbin, int ncomp)
      : tmp("settling_ws_tmp", nlev, ncol),
        airmas("settling_ws_airmas", nlev, ncol),
        p_mid("settling_ws_p_mid", nlev, ncol),
        delz("settling_ws_delz", nlev, ncol),
        rh("settling_ws_rh", nlev, ncol),
        rho_air("settling_ws_rho_air", nlev, ncol),
        v_grav_1d("settling_ws_v_grav_1d", nbin, nlev, ncol),
        avghygro("settling_ws_avghygro", nbin, nlev, ncol),
        avgden("settling_ws_avgden", nbin, nlev, ncol),
        avgreff("settling_ws_avgreff", nbin, nlev, ncol),
        totmass0("settling_ws_totmass0", ncomp, nbin, nlev+1, ncol),
        totmass("settling_ws_totmass", ncomp, nbin, nlev+1, ncol),
        bstl_mass("settling_ws_bstl_mass", ncomp, nbin, nlev+1, ncol),
        dt_settl("settling_ws_dt_settl", nbin, ncol),
        ndt_settl("settling_ws_ndt_settl", nbin, ncol),
        tc1("settling_ws_tc1", ncomp, nbin, nlev+1, ncol),
        vd_delz_dt("settling_ws_vd_delz_dt", nbin, nlev+1, ncol) {}
};

static mosaic_settling_workspace *mosaic_settling_ws = nullptr;
static int mosaic_settling_ws_nlev = -1;
static int mosaic_settling_ws_ncol = -1;
static int mosaic_settling_ws_nbin = -1;
static int mosaic_settling_ws_ncomp = -1;

static bool mosaic_settling_workspace_matches(int nlev, int ncol,
                                               int nbin, int ncomp) {
  return mosaic_settling_ws != nullptr &&
         mosaic_settling_ws_nlev == nlev &&
         mosaic_settling_ws_ncol == ncol &&
         mosaic_settling_ws_nbin == nbin &&
         mosaic_settling_ws_ncomp == ncomp;
}

static mosaic_settling_workspace &get_mosaic_settling_workspace(
    int nlev, int ncol, int nbin, int ncomp) {
  if (!mosaic_settling_workspace_matches(nlev, ncol, nbin, ncomp)) {
    if (mosaic_settling_ws != nullptr) {
      Kokkos::fence();
      delete mosaic_settling_ws;
    }
    mosaic_settling_ws = new mosaic_settling_workspace(
        nlev, ncol, nbin, ncomp);
    mosaic_settling_ws_nlev = nlev;
    mosaic_settling_ws_ncol = ncol;
    mosaic_settling_ws_nbin = nbin;
    mosaic_settling_ws_ncomp = ncomp;
  }
  return *mosaic_settling_ws;
}

extern "C" void mosaic_settling_workspace_finalize() {
  if (mosaic_settling_ws != nullptr) {
    if (Kokkos::is_initialized()) {
      Kokkos::fence();
    }
    delete mosaic_settling_ws;
    mosaic_settling_ws = nullptr;
  }
  mosaic_settling_ws_nlev = -1;
  mosaic_settling_ws_ncol = -1;
  mosaic_settling_ws_nbin = -1;
  mosaic_settling_ws_ncomp = -1;
}

void settling(int &its,int &ite,
          int& kts,int& kte,
          int& nmx,real g0,const real& dyn_visc,
          const real4d& tc,const real2d& tmp,const real2d& p_mid,const real2d& delz,const real2d& airmas,
          const real3d& hygro,const real3d& den,const real3d& reff,real& dt,const real4d& bstl,const real2d& rh,const real3d& v_grav_1d,const real2d& rho_air,int& ncomp_aer1,
          mosaic_settling_workspace& workspace);

extern "C" void c_gocart_settling_scheme(real &dt, flcl_ndarray_t* t_phy_p, flcl_ndarray_t* moist_p,
                                         flcl_ndarray_t* chem_p, flcl_ndarray_t* rho_phy_p, flcl_ndarray_t* dz8w_p,
                                         flcl_ndarray_t* p8w_p, flcl_ndarray_t* p_phy_p, flcl_ndarray_t* dx_p,
                                         real &g, flcl_ndarray_t* v_grav_3d_p, int& its, int& ite, int& kts, int& kte,
                                         int& p_qv,flcl_ndarray_t* ncomp_aer_p,flcl_ndarray_t* massptr_aer_p,flcl_ndarray_t* dens_aer_p,flcl_ndarray_t* hygro_aer_p,flcl_ndarray_t* dcen_sect_p,
                                         int& num_moist,int& num_chem,int& nbin_a){
  const int nlev = kte - kts + 1;
  const int ncol = ite - its + 1;
  auto t_phy_all   = flcl::view_from_ndarray<real**>(*t_phy_p);
  auto rho_phy_all = flcl::view_from_ndarray<real**>(*rho_phy_p);
  auto dz8w_all    = flcl::view_from_ndarray<real**>(*dz8w_p);
  auto p8w_all     = flcl::view_from_ndarray<real**>(*p8w_p);
  auto p_phy_all   = flcl::view_from_ndarray<real**>(*p_phy_p);
  auto dx_all      = flcl::view_from_ndarray<real*>(*dx_p);


  auto moist_all = flcl::view_from_ndarray<real***>(*moist_p);
  auto chem_all  = flcl::view_from_ndarray<real***>(*chem_p);
  auto v_grav_3d_all  = flcl::view_from_ndarray<real***>(*v_grav_3d_p);

  auto ncomp_aer  = flcl::view_from_ndarray<int*>(*ncomp_aer_p);
  auto massptr_aer= flcl::view_from_ndarray<int****>(*massptr_aer_p);
  auto dens_aer   = flcl::view_from_ndarray<real**>(*dens_aer_p);
  auto hygro_aer  = flcl::view_from_ndarray<real**>(*hygro_aer_p);
  auto dcen_sect  = flcl::view_from_ndarray<real**>(*dcen_sect_p);

  int ncomp = ncomp_aer(0);
#ifdef USE_SERIALIZATION
  const int workspace_ncol = 1;
#else
  const int workspace_ncol = CHEM_WORKSPACE_CHUNK_SIZE;
#endif
  auto &workspace = get_mosaic_settling_workspace(
      nlev, workspace_ncol, nbin_a, ncomp);
  auto tmp       = workspace.tmp;
  auto airmas    = workspace.airmas;
  auto p_mid     = workspace.p_mid;
  auto delz      = workspace.delz;
  auto rh        = workspace.rh;
  auto rho_air   = workspace.rho_air;
  auto v_grav_1d = workspace.v_grav_1d;
  auto avghygro  = workspace.avghygro;
  auto avgden    = workspace.avgden;
  auto avgreff   = workspace.avgreff;
  auto totmass0  = workspace.totmass0;
  auto totmass   = workspace.totmass;
  auto bstl_mass = workspace.bstl_mass;

  int itype, iphase;

  constexpr real conver=1.e-9;

  itype = 0;
  iphase = 0;

  // Serial execution reuses one scratch column. Parallel execution reuses a
  // fixed 64-column workspace, including for the tail chunk.
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

  auto t_phy = Kokkos::subview(t_phy_all, Kokkos::ALL(), column_range);
  auto rho_phy = Kokkos::subview(rho_phy_all, Kokkos::ALL(), column_range);
  auto dz8w = Kokkos::subview(dz8w_all, Kokkos::ALL(), column_range);
  auto p8w = Kokkos::subview(p8w_all, Kokkos::ALL(), column_range);
  auto p_phy = Kokkos::subview(p_phy_all, Kokkos::ALL(), column_range);
  auto dx = Kokkos::subview(dx_all, column_range);
  auto moist = Kokkos::subview(
      moist_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
  auto chem = Kokkos::subview(
      chem_all, Kokkos::ALL(), Kokkos::ALL(), column_range);
  auto v_grav_3d = Kokkos::subview(
      v_grav_3d_all, Kokkos::ALL(), Kokkos::ALL(), column_range);

  Kokkos::parallel_for("mosaic_settling_prepare_grid",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
          {0, 0}, {nlev+1, current_chunk_len}),
      KOKKOS_LAMBDA(int k, int i) {
        for (int isize = 0; isize < nbin_a; ++isize) {
          for (int l = 0; l < ncomp_aer(itype); ++l) {
            totmass(l,isize,k,i) = 0.0;
          }
        }
        if (k < nlev) {
          airmas(k,i) = (p8w(k,i) - p8w(k+1,i))*dx(i)*dx(i)/g;
          p_mid(k,i) = 0.01*p_phy(k,i);
          delz(k,i)  = dz8w(k,i);
          tmp(k,i)   = t_phy(k,i);
          rh(k,i)    = moist(p_qv-1,k,i) /
              (3.8*Kokkos::exp(17.27*(t_phy(k,i)-273.) /
              (t_phy(k,i)-36.0))/p_mid(k,i));
          rh(k,i) = Kokkos::min(0.95,Kokkos::max(1.0E-1,rh(k,i)));

          for (int isize = 0; isize < nbin_a; ++isize) {
            real sumhygro = 0.0;
            real sumvol = 0.0;
            real total_mass = 0.0;
            for (int l = 0; l < ncomp_aer(itype); ++l) {
              int lmass = massptr_aer(l,isize,itype,iphase)-1;
              real chem_val = chem(lmass,k,i);
              real comp_vol = chem_val/(dens_aer(l,itype)*1.e12);
              real tm = chem_val*conver;
              sumhygro += hygro_aer(l,itype)*comp_vol;
              sumvol += comp_vol;
              total_mass += tm;
              totmass(l,isize,k,i) = tm;
            }
            if(sumvol > 0.0){
                avghygro(isize,k,i) = sumhygro/sumvol;
                avgden(isize,k,i) = total_mass/sumvol;
            }else{
                avghygro(isize,k,i) = 1e-6;
                avgden(isize,k,i) = 1e-6;
            }
            avgreff(isize,k,i) = dcen_sect(isize,itype)*0.01/2.0;
          }
        }
  });

  Kokkos::deep_copy(totmass0, totmass);

  int chunk_its = 0;
  int chunk_ite = current_chunk_len - 1;
  settling(chunk_its,chunk_ite,kts,kte,
            nbin_a,g,dyn_visc,
            totmass, tmp, p_mid, delz, airmas,
            avghygro,avgden,avgreff, dt, bstl_mass, rh,
            v_grav_1d,rho_air,ncomp,workspace);

  Kokkos::parallel_for("mosaic_settling_writeback",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
          {0, 0}, {current_chunk_len, nlev}),
      KOKKOS_LAMBDA(int i, int k) {
        for (int isize = 0; isize < nbin_a; ++isize) {
          for (int l = 0; l < ncomp_aer(itype); ++l) {
            int lmass=massptr_aer(l,isize,itype,iphase)-1;
            chem(lmass,k,i) = chem(lmass,k,i)
                            * totmass(l,isize,k,i)
                            / totmass0(l,isize,k,i);
            v_grav_3d(lmass,k,i) = v_grav_1d(isize,k,i);
          }
        }
  });

  Kokkos::fence();
  }
}


void settling(int &its,int &ite,
          int& kts,int& kte,
          int& nmx,real g0,const real& dyn_visc,
          const real4d& tc,const real2d& tmp,const real2d& p_mid,const real2d& delz,const real2d& airmas,
          const real3d& hygro,const real3d& den,const real3d& reff,real& dt,const real4d& bstl,const real2d& rh,const real3d& v_grav_1d,const real2d& rho_air,int& ncomp_aer1,
          mosaic_settling_workspace& workspace){

            //loacal variables
            const int nlev = kte - kts + 1;
            const int ncol = ite - its + 1;
            auto dt_settl  = workspace.dt_settl;
            auto ndt_settl = workspace.ndt_settl;
            auto tc1       = workspace.tc1;
            auto vd_delz_dt = workspace.vd_delz_dt;

            //for sea-salt:
            constexpr real c1=0.7674;
            constexpr real c2=3.079;
            constexpr real c3=2.573E-11;
            constexpr real c4=-1.424;

            const int ntdt = int(dt);

            Kokkos::parallel_for("settling_copy_grid",
                Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev+1}),
                KOKKOS_LAMBDA(int i, int k) {
                for (int n = 0; n < nmx; ++n) {
                  for (int l = 0; l < ncomp_aer1; ++l) {
                    tc1(l,n,k,i)=tc(l,n,k,i);
                  }
                }
            });

            Kokkos::parallel_for("settling_substep_dt",
                Kokkos::RangePolicy<>(0,ncol),
                KOKKOS_LAMBDA(int i) {
              for (int n = 0; n < nmx; ++n) {
                real min_dt = 1.0e300;
                for (int k = 0; k < nlev; ++k) {
                  real growth_fac = 1.0 +
                      ((hygro(n,k,i)-hygro_dust_aer)/
                      (hygro_na_aer-hygro_dust_aer))*2.0;
                  real vsettl_val = 2.0/9.0*g0*den(n,k,i)*
                      Kokkos::pow(growth_fac*reff(n,k,i),2)/(0.5*dyn_visc);
                  real dt_val = delz(k,i)/vsettl_val;
                  if (dt_val < min_dt) min_dt=dt_val;
                }
                ndt_settl(n,i)=Kokkos::min(Kokkos::max(1,int(ntdt/min_dt)),15);
                dt_settl(n,i)=real(ntdt)/real(ndt_settl(n,i));
              }
            });

            Kokkos::parallel_for("settling_velocity_grid",
                Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev}),
                KOKKOS_LAMBDA(int i, int k) {
              real rhb_ki=Kokkos::min(9.9E-1,rh(k,i));
              real c_stokes=1.458E-6*Kokkos::pow(tmp(k,i),1.5)/(tmp(k,i)+110.4);
              real free_path=1.1E-3/p_mid(k,i)/Kokkos::sqrt(tmp(k,i));
              for (int n = 0; n < nmx; ++n) {
                real rcm=reff(n,k,i)*100.0;
                real rwet_seas=0.01*Kokkos::pow(
                    c1*Kokkos::pow(rcm,c2)/
                    (c3*Kokkos::pow(rcm,c4)-Kokkos::log10(rhb_ki))+
                    Kokkos::pow(rcm,3),0.33);
                real hygro_fac=(hygro(n,k,i)-hygro_dust_aer)/
                    (hygro_na_aer-hygro_dust_aer);
                real rwet=reff(n,k,i)+hygro_fac*(rwet_seas-reff(n,k,i));
                real ratio_r=Kokkos::pow(reff(n,k,i)/rwet,3.0);
                real rho_priv=ratio_r*den(n,k,i)+(1.0-ratio_r)*1000.0;
                real c_cun=1.0+free_path/rwet*
                    (1.257+0.4*Kokkos::exp(-1.1*rwet/free_path));
                real viscosity=c_stokes/c_cun;
                real vdc=2.0/9.0*g0*rho_priv*Kokkos::pow(rwet,2)/viscosity;
                v_grav_1d(n,k,i)=vdc;
                vd_delz_dt(n,k,i)=dt_settl(n,i)*vdc/delz(k,i);
              }
            });

          Kokkos::parallel_for("settling_apply_column",
              Kokkos::RangePolicy<>(0,ncol),
              KOKKOS_LAMBDA(int i) {
              for (int isize = 0; isize < nmx; ++isize) {
                for (int n = 0; n < ndt_settl(isize,i); ++n) {
                  for (int k = nlev-1; k >= 0; --k) {
                    for (int l = 0; l < ncomp_aer1; ++l) {
                      if (k == nlev-1) {
                        tc(l,isize,k,i)=tc(l,isize,k,i)/
                            (1.0+vd_delz_dt(isize,k,i));
                      }
                      else if (k == 0) {
                        tc(l,isize,k,i)=tc(l,isize,k,i)+
                            tc(l,isize,k+1,i)*vd_delz_dt(isize,k+1,i);
                      }
                      else {
                        tc(l,isize,k,i)=1.0/(1.0+vd_delz_dt(isize,k,i))*
                            (tc(l,isize,k,i)+tc(l,isize,k+1,i)*
                            vd_delz_dt(isize,k+1,i));
                      }
                    }
                  }
                }
              }
            });

          Kokkos::parallel_for("settling_tendency",
              Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev+1}),
              KOKKOS_LAMBDA(int i, int k) {
              for (int isize = 0; isize < nmx; ++isize) {
                for (int l = 0; l < ncomp_aer1; ++l) {
                  tc(l,isize,k,i)=Kokkos::max(1.0E-32,tc(l,isize,k,i));
                  bstl(l,isize,k,i)=(tc(l,isize,k,i)-tc1(l,isize,k,i))*1.e9;
                }
              }
          });
          Kokkos::fence();
}
