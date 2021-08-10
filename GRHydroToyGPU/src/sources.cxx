#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace HydroToyGPU {
using namespace std;
using namespace Loop;

enum class reconstruction_t { Godunov, minmod };

namespace {
template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T minmod(const T &x,
                                                                   const T &y) {
  if (signbit(x) != signbit(y))
    return T(0);
  if (fabs(x) < fabs(y))
    return x;
  else
    return y;
}
} // namespace

// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir> void CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  // const array<CCTK_REAL, dim> dx = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
  //                                   CCTK_DELTA_SPACE(2)};
  // const CCTK_REAL dV = dx[0] * dx[1] * dx[2]; // cell volume
  // const CCTK_REAL dA = dV / dx[dir];          // face area

  // Cell centred grid functions
  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  constexpr array<int, dim> vertex_centred = {0, 0, 0};
  const GF3D2layout gf_layout_cell(cctkGH, cell_centred);
  const GF3D2layout gf_layout_vertex(cctkGH, vertex_centred);

  const GF3D2<const CCTK_REAL> gf_gxx(gf_layout_vertex, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(gf_layout_vertex, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(gf_layout_vertex, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(gf_layout_vertex, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(gf_layout_vertex, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(gf_layout_vertex, gzz);

  const GF3D2<const CCTK_REAL> gf_alp(gf_layout_vertex, alp);
  const GF3D2<const CCTK_REAL> gf_betax(gf_layout_vertex, betax);
  const GF3D2<const CCTK_REAL> gf_betay(gf_layout_vertex, betay);
  const GF3D2<const CCTK_REAL> gf_betaz(gf_layout_vertex, betaz);

  const GF3D2<const CCTK_REAL> gf_dens(gf_layout_cell, dens);
  const GF3D2<const CCTK_REAL> gf_momx(gf_layout_cell, momx);
  const GF3D2<const CCTK_REAL> gf_momy(gf_layout_cell, momy);
  const GF3D2<const CCTK_REAL> gf_momz(gf_layout_cell, momz);
  const GF3D2<const CCTK_REAL> gf_tau(gf_layout_cell, tau);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout_cell, rho);
  const GF3D2<const CCTK_REAL> gf_velx(gf_layout_cell, velx);
  const GF3D2<const CCTK_REAL> gf_vely(gf_layout_cell, vely);
  const GF3D2<const CCTK_REAL> gf_velz(gf_layout_cell, velz);
  const GF3D2<const CCTK_REAL> gf_press(gf_layout_cell, press);

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};
  const GF3D2layout gf_fluxlayout(cctkGH, face_centred);

  // Get the grid function pointers for fluxes in direction `dir`
  const array<CCTK_REAL *, dim> fluxdenss = {fxdens, fydens, fzdens};
  const array<CCTK_REAL *, dim> fluxmomxs = {fxmomx, fymomx, fzmomx};
  const array<CCTK_REAL *, dim> fluxmomys = {fxmomy, fymomy, fzmomy};
  const array<CCTK_REAL *, dim> fluxmomzs = {fxmomz, fymomz, fzmomz};
  const array<CCTK_REAL *, dim> fluxtaus = {fxtau, fytau, fztau};
  const GF3D2<CCTK_REAL> gf_fluxdens(gf_fluxlayout, fluxdenss[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomx(gf_fluxlayout, fluxmomxs[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomy(gf_fluxlayout, fluxmomys[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomz(gf_fluxlayout, fluxmomzs[dir]);
  const GF3D2<CCTK_REAL> gf_fluxtau(gf_fluxlayout, fluxtaus[dir]);

  // fdens^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // ftau^i = (tau + press) vel^i

  reconstruction_t reconstruction;
  if (CCTK_EQUALS(reconstruction_method, "Godunov"))
    reconstruction = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(reconstruction_method, "minmod"))
    reconstruction = reconstruction_t::minmod;
  else
    CCTK_ERROR("Unknown value for parameter \"reconstruction_method\"");

  switch (reconstruction) {
  case reconstruction_t::Godunov:
    assert(cctk_nghostzones[dir] >= 1);
  case reconstruction_t::minmod:
    assert(cctk_nghostzones[dir] >= 2);
  }

  constexpr auto DI = PointDesc::DI;
  const auto reconstruct =
      [=] CCTK_DEVICE CCTK_HOST(
          const GF3D2<const CCTK_REAL> &gf_var,
          const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" cell indices
        const auto Imm = p.I - 2 * DI[dir];
        const auto Im = p.I - DI[dir];
        const auto Ip = p.I;
        const auto Ipp = p.I + DI[dir];

        switch (reconstruction) {
        case reconstruction_t::Godunov: {
          CCTK_REAL var_m = gf_var(Im);
          CCTK_REAL var_p = gf_var(Ip);
          return array<CCTK_REAL, 2>{var_m, var_p};
        }
        case reconstruction_t::minmod: {
          CCTK_REAL var_slope_p = gf_var(Ipp) - gf_var(Ip);
          CCTK_REAL var_slope_c = gf_var(Ip) - gf_var(Im);
          CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
          CCTK_REAL var_m = gf_var(Im) + minmod(var_slope_c, var_slope_m) / 2;
          CCTK_REAL var_p = gf_var(Ip) - minmod(var_slope_p, var_slope_c) / 2;
          return array<CCTK_REAL, 2>{var_m, var_p};
        }
        default:
          CCTK_BUILTIN_UNREACHABLE();
        }
      };

  const auto calcflux = [=] CCTK_DEVICE CCTK_HOST(
                            CCTK_REAL var_m, CCTK_REAL var_p, CCTK_REAL flux_m,
                            CCTK_REAL flux_p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    CCTK_REAL lambda_m = +1.0;
    CCTK_REAL lambda_p = -1.0;
    CCTK_REAL llf =
        0.5 * ((flux_m + flux_p) -
               fmax(fabs(lambda_m), fabs(lambda_p)) * (var_p - var_m));
    // return dA * llf;
    return llf;
  };

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /*
     Set up multipatch stuff
  */

  const CCTK_REAL * restrict vup;
  const CCTK_REAL * restrict Bprim;
  const CCTK_REAL * restrict g11;
  const CCTK_REAL * restrict g12;
  const CCTK_REAL * restrict g13;
  const CCTK_REAL * restrict g22;
  const CCTK_REAL * restrict g23;
  const CCTK_REAL * restrict g33;
  const CCTK_REAL * restrict k11;
  const CCTK_REAL * restrict k12;
  const CCTK_REAL * restrict k13;
  const CCTK_REAL * restrict k22;
  const CCTK_REAL * restrict k23;
  const CCTK_REAL * restrict k33;
  const CCTK_REAL * restrict beta1;
  const CCTK_REAL * restrict beta2;
  const CCTK_REAL * restrict beta3;

  const CCTK_REAL one = 1.00;
  const CCTK_REAL two = 2.00;
  const CCTK_REAL half = 0.50;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  
  const CCTK_REAL ih[3] = { 1.0/dx, 1.0/dy, 1.0/dz };

  grid.loop_int_device<1,1,1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // TODO: these are not needed.
        const CCTK_REAL localgxx = g11[ijk];
        const CCTK_REAL localgxy = g12[ijk];
        const CCTK_REAL localgxz = g13[ijk];
        const CCTK_REAL localgyy = g22[ijk];
        const CCTK_REAL localgyz = g23[ijk];
        const CCTK_REAL localgzz = g33[ijk];

        const CCTK_REAL sqrtdet = sdetg[ijk];
        const CCTK_REAL invsqrtdet = 1./sqrtdet;
        
        CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
        UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, sqrtdet*sqrtdet, localgxx,
                    localgxy, localgxz, localgyy, localgyz, localgzz);
        
        const CCTK_REAL shiftx = beta1[ijk];
        const CCTK_REAL shifty = beta2[ijk];
        const CCTK_REAL shiftz = beta3[ijk];

        //  Derivatives of the lapse, metric and shift

        const int nvars = 10;
        const CCTK_REAL* const restrict vars[nvars] = { beta1, beta2, beta3, alp,
                                                     g11, g12, g13,
                                                     g22, g23,
                                                     g33 };
        CCTK_REAL dvars[nvars][3];
        
        alldiff::apply(cctkGH, dvars, vars, i, j, k, ih, nvars);
        
        const CCTK_REAL dx_betax = dvars[0][0];
        const CCTK_REAL dx_betay = dvars[1][0];
        const CCTK_REAL dx_betaz = dvars[2][0];
           
        const CCTK_REAL dy_betax = dvars[0][1];
        const CCTK_REAL dy_betay = dvars[1][1];
        const CCTK_REAL dy_betaz = dvars[2][1];
           
        const CCTK_REAL dz_betax = dvars[0][2];
        const CCTK_REAL dz_betay = dvars[1][2];
        const CCTK_REAL dz_betaz = dvars[2][2];
        
        const CCTK_REAL dx_alp = dvars[3][0];
        const CCTK_REAL dy_alp = dvars[3][1];
        const CCTK_REAL dz_alp = dvars[3][2];

        const CCTK_REAL dx_gxx = dvars[4][0];
        const CCTK_REAL dx_gxy = dvars[5][0];
        const CCTK_REAL dx_gxz = dvars[6][0];
        const CCTK_REAL dx_gyy = dvars[7][0];
        const CCTK_REAL dx_gyz = dvars[8][0];
        const CCTK_REAL dx_gzz = dvars[9][0];
        const CCTK_REAL dy_gxx = dvars[4][1];
        const CCTK_REAL dy_gxy = dvars[5][1];
        const CCTK_REAL dy_gxz = dvars[6][1];
        const CCTK_REAL dy_gyy = dvars[7][1];
        const CCTK_REAL dy_gyz = dvars[8][1];
        const CCTK_REAL dy_gzz = dvars[9][1];
        const CCTK_REAL dz_gxx = dvars[4][2];
        const CCTK_REAL dz_gxy = dvars[5][2];
        const CCTK_REAL dz_gxz = dvars[6][2];
        const CCTK_REAL dz_gyy = dvars[7][2];
        const CCTK_REAL dz_gyz = dvars[8][2];
        const CCTK_REAL dz_gzz = dvars[9][2];


        const CCTK_REAL invalp = 1.0 / alp[ijk];
        const CCTK_REAL invalp2 = SQR(invalp);
        const CCTK_REAL velxshift = velx[ijk] - shiftx*invalp;
        const CCTK_REAL velyshift = vely[ijk] - shifty*invalp;
        const CCTK_REAL velzshift = velz[ijk] - shiftz*invalp;

        // vel_i  = g_ij v^j
        // B_i = g_ij B^i

        const CCTK_REAL vlowx = g11[ijk]*velx[ijk] + g12[ijk]*vely[ijk] + g13[ijk]*velz[ijk];
        const CCTK_REAL vlowy = g12[ijk]*velx[ijk] + g22[ijk]*vely[ijk] + g23[ijk]*velz[ijk];
        const CCTK_REAL vlowz = g13[ijk]*velx[ijk] + g23[ijk]*vely[ijk] + g33[ijk]*velz[ijk];
        const CCTK_REAL Bvecxlow = do_mhd ? g11[ijk]*Bvecx[ijk] + g12[ijk]*Bvecy[ijk] + g13[ijk]*Bvecz[ijk] : 0;
        const CCTK_REAL Bvecylow = do_mhd ? g12[ijk]*Bvecx[ijk] + g22[ijk]*Bvecy[ijk] + g23[ijk]*Bvecz[ijk] : 0;
        const CCTK_REAL Bveczlow = do_mhd ? g13[ijk]*Bvecx[ijk] + g23[ijk]*Bvecy[ijk] + g33[ijk]*Bvecz[ijk] : 0;

        //  B^i v_i (= b^0/u^0)
        const CCTK_REAL Bdotv = do_mhd ? vlowx*Bvecx[ijk]+vlowy*Bvecy[ijk]+vlowz*Bvecz[ijk] : 0;

        // v^2 = v_i v^i; w=(1-v^2)^{-1/2}

        const CCTK_REAL v2 = vlowx*velx[ijk] + vlowy*vely[ijk] + vlowz*velz[ijk];
        const CCTK_REAL invw = sqrt(1.0-v2);
        const CCTK_REAL w = 1./invw;

        // b^2 = B^i B_i / w^2 + (b^0/u^0)^2

        const CCTK_REAL b2 = do_mhd ? (Bvecx[ijk]*Bvecxlow+Bvecy[ijk]*Bvecylow+Bvecz[ijk]*Bveczlow)*SQR(invw)+SQR(Bdotv) : 0;

        // b_i = B_i/w +w*(B dot v)*v_i
        const CCTK_REAL bxlow = do_mhd ? Bvecxlow*invw+w*Bdotv*vlowx : 0;
        const CCTK_REAL bylow = do_mhd ? Bvecylow*invw+w*Bdotv*vlowy : 0;
        const CCTK_REAL bzlow = do_mhd ? Bveczlow*invw+w*Bdotv*vlowz : 0;


        // These are the contravariant components
        const CCTK_REAL bt = do_mhd ? w*invalp*Bdotv : 0;
        const CCTK_REAL bx = do_mhd ? Bvecx[ijk]*invw+w*Bdotv*velxshift : 0;
        const CCTK_REAL by = do_mhd ? Bvecy[ijk]*invw+w*Bdotv*velyshift : 0;
        const CCTK_REAL bz = do_mhd ? Bvecz[ijk]*invw+w*Bdotv*velzshift : 0;

        // TODO: all of these can be expressed much more easily in terms of the
        // conservatives
        const CCTK_REAL rhohstarW2 = (rho[ijk]*(one + eps[ijk]) + (press[ijk] + b2)) *
                                       SQR(w);
        const CCTK_REAL pstar = press[ijk]+0.50*b2;

        //  For a change, these are T^{ij}

        // TODO: strict IEEE compliance does not let the compiler remove "+ 0"
        // terms, so we have to do something else here
        const CCTK_REAL t00 = (rhohstarW2 - pstar)*invalp2-SQR(bt);
        const CCTK_REAL t0x = rhohstarW2*velxshift*invalp +
             pstar*shiftx*invalp2-bt*bx;
        const CCTK_REAL t0y = rhohstarW2*velyshift*invalp +
             pstar*shifty*invalp2-bt*by;
        const CCTK_REAL t0z = rhohstarW2*velzshift*invalp +
             pstar*shiftz*invalp2-bt*bz;
        const CCTK_REAL txx = rhohstarW2*velxshift*velxshift +
             pstar*(uxx - shiftx*shiftx*invalp2)-SQR(bx);
        const CCTK_REAL txy = rhohstarW2*velxshift*velyshift +
             pstar*(uxy - shiftx*shifty*invalp2)-bx*by;
        const CCTK_REAL txz = rhohstarW2*velxshift*velzshift +
             pstar*(uxz - shiftx*shiftz*invalp2)-bx*bz;
        const CCTK_REAL tyy = rhohstarW2*velyshift*velyshift +
             pstar*(uyy - shifty*shifty*invalp2)-SQR(by);
        const CCTK_REAL tyz = rhohstarW2*velyshift*velzshift +
             pstar*(uyz - shifty*shiftz*invalp2)-by*bz;
        const CCTK_REAL tzz = rhohstarW2*velzshift*velzshift +
             pstar*(uzz - shiftz*shiftz*invalp2)-SQR(bz);

//        Contract the shift with the extrinsic curvature

        const CCTK_REAL shiftshiftk = shiftx*shiftx*k11[ijk] +
                                      shifty*shifty*k22[ijk] +
                                      shiftz*shiftz*k33[ijk] +
             two*(shiftx*shifty*k12[ijk] +
                  shiftx*shiftz*k13[ijk] +
                  shifty*shiftz*k23[ijk]);

        const CCTK_REAL shiftkx = shiftx*k11[ijk] + shifty*k12[ijk] + shiftz*k13[ijk];
        const CCTK_REAL shiftky = shiftx*k12[ijk] + shifty*k22[ijk] + shiftz*k23[ijk];
        const CCTK_REAL shiftkz = shiftx*k13[ijk] + shifty*k23[ijk] + shiftz*k33[ijk];

//        Contract the matter terms with the extrinsic curvature

        const CCTK_REAL sumTK = txx*k11[ijk] + tyy*k22[ijk] + tzz*k33[ijk]
                         + two*(txy*k12[ijk] + txz*k13[ijk] + tyz*k23[ijk]);

//        Update term for tau
        
        const CCTK_REAL tau_source = t00*
             (shiftshiftk - (shiftx*dx_alp + shifty*dy_alp + shiftz*dz_alp) )
             + t0x*(-dx_alp + two*shiftkx)
             + t0y*(-dy_alp + two*shiftky)
             + t0z*(-dz_alp + two*shiftkz)
             + sumTK;

//        The following looks very little like the terms in the
//        standard papers. Take a look in the ThornGuide to see why
//        it is really the same thing.

//        Contract the shift with derivatives of the metric

        const CCTK_REAL halfshiftdgx = half*(shiftx*shiftx*dx_gxx +
             shifty*shifty*dx_gyy + shiftz*shiftz*dx_gzz) +
             shiftx*shifty*dx_gxy + shiftx*shiftz*dx_gxz +
             shifty*shiftz*dx_gyz;
        const CCTK_REAL halfshiftdgy = half*(shiftx*shiftx*dy_gxx +
             shifty*shifty*dy_gyy + shiftz*shiftz*dy_gzz) +
             shiftx*shifty*dy_gxy + shiftx*shiftz*dy_gxz +
             shifty*shiftz*dy_gyz;
        const CCTK_REAL halfshiftdgz = half*(shiftx*shiftx*dz_gxx +
             shifty*shifty*dz_gyy + shiftz*shiftz*dz_gzz) +
             shiftx*shifty*dz_gxy + shiftx*shiftz*dz_gxz +
             shifty*shiftz*dz_gyz;

//        Contract the matter with derivatives of the metric

        const CCTK_REAL halfTdgx = half*(txx*dx_gxx + tyy*dx_gyy + tzz*dx_gzz) +
             txy*dx_gxy + txz*dx_gxz + tyz*dx_gyz;
        const CCTK_REAL halfTdgy = half*(txx*dy_gxx + tyy*dy_gyy + tzz*dy_gzz) +
             txy*dy_gxy + txz*dy_gxz + tyz*dy_gyz;
        const CCTK_REAL halfTdgz = half*(txx*dz_gxx + tyy*dz_gyy + tzz*dz_gzz) +
             txy*dz_gxy + txz*dz_gxz + tyz*dz_gyz;

     
       const CCTK_REAL sx_source = t00*
             (halfshiftdgx - alp[ijk]*dx_alp) + halfTdgx +
             t0x*(shiftx*dx_gxx + shifty*dx_gxy + shiftz*dx_gxz) +
             t0y*(shiftx*dx_gxy + shifty*dx_gyy + shiftz*dx_gyz) +
             t0z*(shiftx*dx_gxz + shifty*dx_gyz + shiftz*dx_gzz) +
             rhohstarW2*invalp*(vlowx*dx_betax + vlowy*dx_betay + vlowz*dx_betaz) -
             bt*(bxlow*dx_betax + bylow*dx_betay + bzlow*dx_betaz);
        
       const CCTK_REAL sy_source = t00*
             (halfshiftdgy - alp[ijk]*dy_alp) + halfTdgy +
             t0x*(shiftx*dy_gxx + shifty*dy_gxy + shiftz*dy_gxz) +
             t0y*(shiftx*dy_gxy + shifty*dy_gyy + shiftz*dy_gyz) +
             t0z*(shiftx*dy_gxz + shifty*dy_gyz + shiftz*dy_gzz) +
             rhohstarW2*invalp*(vlowx*dy_betax + vlowy*dy_betay + vlowz*dy_betaz) -
             bt*(bxlow*dy_betax + bylow*dy_betay + bzlow*dy_betaz);

       const CCTK_REAL sz_source = t00*
             (halfshiftdgz - alp[ijk]*dz_alp) + halfTdgz +
             t0x*(shiftx*dz_gxx + shifty*dz_gxy + shiftz*dz_gxz) +
             t0y*(shiftx*dz_gxy + shifty*dz_gyy + shiftz*dz_gyz) +
             t0z*(shiftx*dz_gxz + shifty*dz_gyz + shiftz*dz_gzz) +
             rhohstarW2*invalp*(vlowx*dz_betax + vlowy*dz_betay + vlowz*dz_betaz) -
             bt*(bxlow*dz_betax + bylow*dz_betay + bzlow*dz_betaz);

        densrhs[ijk] = 0.0;
        srhs[ijk]        = alp[ijk]*sqrtdet*sx_source;
        srhs[ijk + N]    = alp[ijk]*sqrtdet*sy_source;
        srhs[ijk + 2*N]  = alp[ijk]*sqrtdet*sz_source;
        taurhs[ijk]      = alp[ijk]*sqrtdet*tau_source;

        if (do_Avec) {

          // B^i and A^i both live in cell centers currently
          const CCTK_REAL Avcx_source = alp[ijk]*sqrtdet*(velyshift*Bvecz[ijk] - velzshift*Bvecy[ijk]);
          const CCTK_REAL Avcy_source = alp[ijk]*sqrtdet*(velzshift*Bvecx[ijk] - velxshift*Bvecz[ijk]);
          const CCTK_REAL Avcz_source = alp[ijk]*sqrtdet*(velxshift*Bvecy[ijk] - velyshift*Bvecx[ijk]);

          Avecrhsx[ijk] = Avcx_source;
          Avecrhsy[ijk] = Avcy_source;
          Avecrhsz[ijk] = Avcz_source;

        }

        if(do_clean_divergence) {
   
           // g^{jk} d_i g_{kj} = d_i (g) / det
           const CCTK_REAL dx_det_bydet = uxx*dx_gxx + uyy*dx_gyy + uzz*dx_gzz +
                two*(uxy*dx_gxy+uxz*dx_gxz+uyz*dx_gyz);
           const CCTK_REAL dy_det_bydet = uxx*dy_gxx + uyy*dy_gyy + uzz*dy_gzz +
                two*(uxy*dy_gxy+uxz*dy_gxz+uyz*dy_gyz);
           const CCTK_REAL dz_det_bydet = uxx*dz_gxx + uyy*dz_gyy + uzz*dz_gzz +
                two*(uxy*dz_gxy+uxz*dz_gxz+uyz*dz_gyz);

           // g^{ik} d_k g_{li}
           const CCTK_REAL gdg_x = uxx*dx_gxx + uxy*dy_gxx + uxz*dz_gxx +
                   uxy*dx_gxy + uyy*dy_gxy + uyz*dz_gxy +
                   uxz*dx_gxz + uyz*dy_gxz + uzz*dz_gxz;

           const CCTK_REAL gdg_y = uxx*dx_gxy + uxy*dy_gxy + uxz*dz_gxy +
                   uxy*dx_gyy + uyy*dy_gyy + uyz*dz_gyy +
                   uxz*dx_gyz + uyz*dy_gyz + uzz*dz_gyz;

           const CCTK_REAL gdg_z = uxx*dx_gxz + uxy*dy_gxz + uxz*dz_gxz +
                   uxy*dx_gyz + uyy*dy_gyz + uyz*dz_gyz +
                   uxz*dx_gzz + uyz*dy_gzz + uzz*dz_gzz;

           CCTK_REAL bvcx_source, bvcy_source, bvcz_source;
           psidcrhs[ijk] = -one * (kap_dc*alp[ijk] + 
                dx_betax + dy_betay + dz_betaz ) * psidc[ijk] + 
                Bconsx[ijk] * (dx_alp - half*alp[ijk] * 
                 ( uxx*dx_gxx + uyy*dx_gyy + uzz*dx_gzz + two*uxy*dx_gxy + 
                   two*uxz*dx_gxz + two*uyz*dx_gyz ) )*invsqrtdet + 
                Bconsy[ijk] * (dy_alp - half*alp[ijk] * 
                 ( uxx*dy_gxx + uyy*dy_gyy + uzz*dy_gzz + two*uxy*dy_gxy + 
                   two*uxz*dy_gxz + two*uyz*dy_gyz ) )*invsqrtdet + 
                Bconsz[ijk] * (dz_alp - half*alp[ijk] * 
                 ( uxx*dz_gxx + uyy*dz_gyy + uzz*dz_gzz + two*uxy*dz_gxy + 
                   two*uxz*dz_gxz + two*uyz*dz_gyz ) )*invsqrtdet;

           bvcx_source = -one * ( Bconsx[ijk]*dx_betax + 
                Bconsy[ijk]*dy_betax + Bconsz[ijk]*dz_betax ) + 
                psidc[ijk]*sqrtdet*(( uxx*dx_alp+uxy*dy_alp+uxz*dz_alp ) + 
                alp[ijk]*(half*( uxx*dx_det_bydet + 
                  uxy*dy_det_bydet + uxz*dz_det_bydet) - 
                ( uxx*gdg_x + uxy*gdg_y + uxz*gdg_z )));

           bvcy_source = -one * ( Bconsx[ijk]*dx_betay + 
                Bconsy[ijk]*dy_betay + Bconsz[ijk]*dz_betay ) + 
                psidc[ijk]*sqrtdet*(( uxy*dx_alp+uyy*dy_alp+uyz*dz_alp ) + 
                alp[ijk]*(half*( uxy*dx_det_bydet + 
                  uyy*dy_det_bydet + uyz*dz_det_bydet ) - 
                ( uxy*gdg_x + uyy*gdg_y + uyz*gdg_z )));

           bvcz_source = -one * ( Bconsx[ijk]*dx_betaz + 
                Bconsy[ijk]*dy_betaz + Bconsz[ijk]*dz_betaz ) + 
                psidc[ijk]*sqrtdet*(( uxz*dx_alp+uyz*dy_alp+uzz*dz_alp ) + 
                alp[ijk]*(half*( uxz*dx_det_bydet + 
                  uyz*dy_det_bydet + uzz*dz_det_bydet ) - 
                ( uxz*gdg_x + uyz*gdg_y + uzz*gdg_z )));

           Bconsrhsx[ijk] = bvcx_source;
           Bconsrhsy[ijk] = bvcy_source;
           Bconsrhsz[ijk] = bvcz_source;
        } // if(do_clean_divergence)



                          /// OLD CODE BELOW
        // Reconstruct values from the cells on left and right side of this face
        array<CCTK_REAL, 2> rho_r = reconstruct(gf_rho, p);
        array<CCTK_REAL, 2> velx_r = reconstruct(gf_velx, p);
        array<CCTK_REAL, 2> vely_r = reconstruct(gf_vely, p);
        array<CCTK_REAL, 2> velz_r = reconstruct(gf_velz, p);
        array<CCTK_REAL, 2> press_r = reconstruct(gf_press, p);
        array<array<CCTK_REAL, 2>, 3> vels_r = {velx_r, vely_r, velz_r};
        array<CCTK_REAL, 2> vel_r = vels_r[dir];

        array<CCTK_REAL, 2> tau_r;
        for (int f = 0; f < 2; ++f) {
          CCTK_REAL ekin =
              0.5 * rho_r[f] *
              (pow2(velx_r[f]) + pow2(vely_r[f]) + pow2(velz_r[f]));
          CCTK_REAL eint = press_r[f] / (gamma - 1);
          tau_r[f] = ekin + eint;
        }

        CCTK_REAL gxx_r = 0;
        // loop over cube corners except for direction in whihc we reconstruct
        for(int di = 0 ; di < 2 - (dir == 0) ; ++di) {
        for(int dj = 0 ; dj < 2 - (dir == 1) ; ++dj) {
        for(int dk = 0 ; dk < 2 - (dir == 2) ; ++dk) {
          gxx_r += gf_gxx(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
        }}}
        gxx_r /= 4;
        // TODO: do the same for other metric variables
        CCTK_REAL sqrt_detg_r = sqrt(gxx_r...);

        array<CCTK_REAL, 2> dens_r = {
                sqrt_detg_r * rho_r[0] * w_lorentz_r[0],
                sqrt_detg_r * rho_r[1] * w_lorentz_r[1]
        };
        const CCTK_REAL flux_dens_m = dens_r[0] * (vel_r[0] - beta_r[0] / alp_r[0]);
        const CCTK_REAL flux_dens_p = dens_r[1] * (vel_r[1] - beta_r[1] / alp_r[1]);

        gf_fluxdens(p.I) = calcflux(dens_r[0], dens_r[1],
                                    flux_dens_m, flux_dens_p);
       
        gf_fluxmomx(p.I) =
            calcflux(rho_r[0] * velx_r[0], rho_r[1] * velx_r[1],
                     rho_r[0] * velx_r[0] * vel_r[0] + (dir == 0) * press_r[0],
                     rho_r[1] * velx_r[1] * vel_r[1] + (dir == 0) * press_r[1]);
        gf_fluxmomy(p.I) =
            calcflux(rho_r[0] * vely_r[0], rho_r[1] * vely_r[1],
                     rho_r[0] * vely_r[0] * vel_r[0] + (dir == 1) * press_r[0],
                     rho_r[1] * vely_r[1] * vel_r[1] + (dir == 1) * press_r[1]);
        gf_fluxmomz(p.I) =
            calcflux(rho_r[0] * velz_r[0], rho_r[1] * velz_r[1],
                     rho_r[0] * velz_r[0] * vel_r[0] + (dir == 2) * press_r[0],
                     rho_r[1] * velz_r[1] * vel_r[1] + (dir == 2) * press_r[1]);
        gf_fluxtau(p.I) =
            calcflux(tau_r[0], tau_r[1], (tau_r[0] + press_r[0]) * vel_r[0],
                     (tau_r[1] + press_r[1]) * vel_r[1]);
      });
}

extern "C" void HydroToyGPU_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  CalcFlux<0>(cctkGH);
  CalcFlux<1>(cctkGH);
  CalcFlux<2>(cctkGH);
}

} // namespace HydroToyGPU
