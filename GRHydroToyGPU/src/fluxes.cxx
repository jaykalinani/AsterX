#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>


namespace GRHydroToyGPU {
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
        DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Fluxes;
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
        const GF3D2<const CCTK_REAL> gf_eps(gf_layout_cell, eps);

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
            const auto Im  = p.I - DI[dir];
            const auto Ip  = p.I;
            const auto Ipp = p.I + DI[dir];

            switch (reconstruction) {
                case reconstruction_t::Godunov: {
                    CCTK_REAL var_m = gf_var(Im);
                    CCTK_REAL var_p = gf_var(Ip);
                    return array<CCTK_REAL, 2> {var_m, var_p};
                }

                case reconstruction_t::minmod: {
                    CCTK_REAL var_slope_p = gf_var(Ipp) - gf_var(Ip);
                    CCTK_REAL var_slope_c = gf_var(Ip)  - gf_var(Im);
                    CCTK_REAL var_slope_m = gf_var(Im)  - gf_var(Imm);
                    CCTK_REAL var_m = gf_var(Im) + minmod(var_slope_c, var_slope_m) / 2;
                    CCTK_REAL var_p = gf_var(Ip) - minmod(var_slope_p, var_slope_c) / 2;
                    return array<CCTK_REAL, 2>{var_m, var_p};
                }

                default:
                    CCTK_BUILTIN_UNREACHABLE();
            }
        };


        const auto calcflux =
            [=] CCTK_DEVICE CCTK_HOST(
                array<CCTK_REAL, 2> var, array<CCTK_REAL, 2> flux)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                array<CCTK_REAL, 2> lambda = {+1.0, -1.0};
            CCTK_REAL llf = 0.5 * ((flux[0] + flux[1])
                          - fmax(fabs(lambda[0]), fabs(lambda[1])) * (var[1] - var[0]));
            // return dA * llf;
            return llf;
        };



        grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
            grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

            // Reconstruct primitives from the cells on left (indice 0) and right (indice 1) side of this face
            // rc = reconstructed variables or computed from reconstructed variables 
            const array<CCTK_REAL, 2> rho_rc   = reconstruct(gf_rho, p);
            const array<CCTK_REAL, 2> velx_rc  = reconstruct(gf_velx, p);
            const array<CCTK_REAL, 2> vely_rc  = reconstruct(gf_vely, p);
            const array<CCTK_REAL, 2> velz_rc  = reconstruct(gf_velz, p);
            const array<CCTK_REAL, 2> press_rc = reconstruct(gf_press, p);

            const array<array<CCTK_REAL, 2>, 3> vels_rc = {velx_rc, vely_rc, velz_rc};
            const array<CCTK_REAL, 2> vel_rc = vels_rc[dir];
        
            //TODO: to reconstruct w_lorentz*vel or 4-velocity u_i

            // Computing metric components
            CCTK_REAL alp_avg = 0;

            CCTK_REAL betax_avg = 0;
            CCTK_REAL betay_avg = 0;
            CCTK_REAL betaz_avg = 0;

            CCTK_REAL gxx_avg = 0;
            CCTK_REAL gxy_avg = 0;
            CCTK_REAL gxz_avg = 0;
            CCTK_REAL gyy_avg = 0;
            CCTK_REAL gyz_avg = 0;
            CCTK_REAL gzz_avg = 0;

            for(int dk = 0 ; dk < (dir == 2 ? 1 : 2) ; ++dk) {
                for(int dj = 0 ; dj < (dir == 1 ? 1 : 2) ; ++dj) {
                    for(int di = 0 ; di < (dir == 0 ? 1 : 2) ; ++di) {
                        alp_avg += gf_alp(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);

                        betax_avg += gf_betax(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        betay_avg += gf_betay(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        betaz_avg += gf_betaz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                     
                        gxx_avg += gf_gxx(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gxy_avg += gf_gxy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gxz_avg += gf_gxz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gyy_avg += gf_gyy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gyz_avg += gf_gyz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gzz_avg += gf_gzz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                    }
                }
            }

            alp_avg   /= 4;
            betax_avg /= 4;
            betay_avg /= 4;
            betaz_avg /= 4;
            gxx_avg   /= 4;
            gxy_avg   /= 4;
            gxz_avg   /= 4;
            gyy_avg   /= 4;
            gyz_avg   /= 4;
            gzz_avg   /= 4;

            const array<CCTK_REAL, 3> betas_avg = {betax_avg, betay_avg, betaz_avg};
            const CCTK_REAL beta_avg = betas_avg[dir];

            //TODO: Compute specific internal energy based on user-specified EOS.
            //      Currently, computing eps for classical ideal gas

            const array<CCTK_REAL, 2> eps_rc = {
                press_rc[0]/(rho_rc[0] * (gamma - 1)),
                press_rc[1]/(rho_rc[1] * (gamma - 1))
            };

            // Determinant of spatial metric
            const CCTK_REAL detg = - gxz_avg*gxz_avg*gyy_avg
                                   + 2.0*gxy_avg*gxz_avg*gyz_avg
                                   - gxx_avg*gyz_avg*gyz_avg
                                   - gxy_avg*gxy_avg*gzz_avg
                                   + gxx_avg*gyy_avg*gzz_avg;
            const CCTK_REAL sqrt_detg = sqrt(detg);



            // v_j
            const array<CCTK_REAL, 2> vlowx_rc = { 
                gxx_avg*velx_rc[0] + gxy_avg*vely_rc[0] + gxz_avg*velz_rc[0],
                gxx_avg*velx_rc[1] + gxy_avg*vely_rc[1] + gxz_avg*velz_rc[1]
            };

            const array<CCTK_REAL, 2> vlowy_rc = {
                gxy_avg*velx_rc[0] + gyy_avg*vely_rc[0] + gyz_avg*velz_rc[0],
                gxy_avg*velx_rc[1] + gyy_avg*vely_rc[1] + gyz_avg*velz_rc[1]
            };

            const array<CCTK_REAL, 2> vlowz_rc = {
                gxz_avg*velx_rc[0] + gyz_avg*vely_rc[0] + gzz_avg*velz_rc[0],
                gxz_avg*velx_rc[1] + gyz_avg*vely_rc[1] + gzz_avg*velz_rc[1]
            };

            // w_lorentz
            const array<CCTK_REAL, 2> w_lorentz_rc = {
                1.0 / sqrt(1 - (vlowx_rc[0]*velx_rc[0] +
                                vlowy_rc[0]*vely_rc[0] +
                                vlowz_rc[0]*velz_rc[0])),
                1.0 / sqrt(1 - (vlowx_rc[1]*velx_rc[1] +
                                vlowy_rc[1]*vely_rc[1] +
                                vlowz_rc[1]*velz_rc[1]))
            };



            // Auxiliary variables to compute conservatives
            const CCTK_REAL aux_mom_0 = 1 + eps_rc[0] + press_rc[0]/rho_rc[0];
            const CCTK_REAL aux_mom_1 = 1 + eps_rc[1] + press_rc[1]/rho_rc[1];

            // Computing conservatives from primitives
            const array<CCTK_REAL, 2> dens_rc = {
                sqrt_detg * rho_rc[0] * w_lorentz_rc[0],
                sqrt_detg * rho_rc[1] * w_lorentz_rc[1]
            };

            const array<CCTK_REAL, 2> momx_rc = {
                dens_rc[0] * aux_mom_0 * vlowx_rc[0],
                dens_rc[1] * aux_mom_1 * vlowx_rc[1]
            };

            const array<CCTK_REAL, 2> momy_rc = {
                dens_rc[0] * aux_mom_0 * vlowy_rc[0],
                dens_rc[1] * aux_mom_1 * vlowy_rc[1]
            };

            const array<CCTK_REAL, 2> momz_rc = {
                dens_rc[0] * aux_mom_0 * vlowz_rc[0],
                dens_rc[1] * aux_mom_1 * vlowz_rc[1]
            };

            const array<CCTK_REAL, 2> tau_rc = {
                dens_rc[0] * (aux_mom_0 * w_lorentz_rc[0] - 1) - press_rc[0],
                dens_rc[1] * (aux_mom_1 * w_lorentz_rc[1] - 1) - press_rc[1]
            }; 



            // Auxiliary variables to compute fluxes of conservatives
            const CCTK_REAL aux_flux_0a = vel_rc[0] - beta_avg/alp_avg;
            const CCTK_REAL aux_flux_0b = sqrt_detg*press_rc[0];

            const CCTK_REAL aux_flux_1a = vel_rc[1] - beta_avg/alp_avg;
            const CCTK_REAL aux_flux_1b = sqrt_detg*press_rc[1];

            // Computing fluxes of conserved variables
            const array<CCTK_REAL, 2> flux_dens = {
                dens_rc[0]*aux_flux_0a, dens_rc[1]*aux_flux_1a
            };
            const array<CCTK_REAL, 2> flux_momx = {
                momx_rc[0]*aux_flux_0a + (dir == 0)*aux_flux_0b,
                momx_rc[1]*aux_flux_1a + (dir == 0)*aux_flux_1b
            };

            const array<CCTK_REAL, 2> flux_momy = {
                momy_rc[0]*aux_flux_0a + (dir == 1)*aux_flux_0b,
                momy_rc[1]*aux_flux_1a + (dir == 1)*aux_flux_1b
            };

            const array<CCTK_REAL, 2> flux_momz = {
                momz_rc[0]*aux_flux_0a + (dir == 2)*aux_flux_0b,
                momz_rc[1]*aux_flux_1a + (dir == 2)*aux_flux_1b
            };

            const array<CCTK_REAL, 2> flux_tau = {
                tau_rc[0]*aux_flux_0a + aux_flux_0b*vel_rc[0],
                tau_rc[1]*aux_flux_1a + aux_flux_1b*vel_rc[1]
            };

            gf_fluxdens(p.I) = calcflux(dens_rc, flux_dens);
            gf_fluxmomx(p.I) = calcflux(momx_rc, flux_momx);
            gf_fluxmomy(p.I) = calcflux(momy_rc, flux_momy);
            gf_fluxmomz(p.I) = calcflux(momz_rc, flux_momz);
            gf_fluxtau(p.I)  = calcflux(tau_rc,  flux_tau);
        });
    }



    extern "C" void GRHydroToyGPU_Fluxes(CCTK_ARGUMENTS) {
      DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Fluxes;
      DECLARE_CCTK_PARAMETERS;

      CalcFlux<0>(cctkGH);
      CalcFlux<1>(cctkGH);
      CalcFlux<2>(cctkGH);
    }

} // namespace GRHydroToyGPU



