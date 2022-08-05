
/* for debug output the location of the points */
void TOV_debug_input_points(CCTK_INT size,
                       const CCTK_REAL *x,
                       const CCTK_REAL *y,
                       const CCTK_REAL *z);
void TOV_debug_input_points(CCTK_INT size,
                       const CCTK_REAL *x,
                       const CCTK_REAL *y,
                       const CCTK_REAL *z)
{
    FILE *f;
    int i;
    f=fopen("test.out", "w");
    for (i=0; i<size; i++)
      fprintf(f, "%.8g %.8g %.8g\n", x[i], y[i], z[i]);
    fclose(f);
}

static inline int get_nearest_star(const CCTK_REAL *x, const CCTK_REAL *y,
                            const CCTK_REAL *z, const int i)
{
    DECLARE_CCTK_PARAMETERS
    int star = 0, star_loop;
    CCTK_REAL r_tmp, r_to_star=1.e40;  /* something very large */
    for (star_loop=0; star_loop<TOV_Num_TOVs; star_loop++)
    {
        r_tmp = sqrt( (x[i]-TOV_Position_x[star_loop]) *
                      (x[i]-TOV_Position_x[star_loop]) +
                      (y[i]-TOV_Position_y[star_loop]) *
                      (y[i]-TOV_Position_y[star_loop]) +
                      (z[i]-TOV_Position_z[star_loop]) *
                      (z[i]-TOV_Position_z[star_loop]) );
        if (r_tmp < r_to_star)
        {
            r_to_star = r_tmp;
            star = star_loop;
        }
    }
    return star;
}
/* This function provides an interface for other thorns which can input a list
 * of points (p_x, p_y, p_z) with a size of 'size' and get back the rho_adm
 * at these points (source) */
/* NOTE: when used with two or more NS, there might be problems if they are
 * very close to each other, because this function does not use the fancy
 * algorithms above to determine from which NS the data for a perticular grid
 * point should be taken, but just uses the closest one */
CCTK_INT TOV_Set_Rho_ADM(CCTK_POINTER_TO_CONST cctkGH,
                         CCTK_INT  size,
                         CCTK_REAL *source,
                         const CCTK_REAL *x,
                         const CCTK_REAL *y,
                         const CCTK_REAL *z )
{
  DECLARE_CCTK_PARAMETERS

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);

  /* CCTK_INFO("Request for Sources from TOVSolver"); */
  /* TOV_debug_input_points(size,x,y,z); */
  /* loop over all points we got */

  #pragma omp parallel for
  for (CCTK_INT i=0; i<size; i++)
  {
    CCTK_REAL r_to_star, press, phi, r;
    CCTK_REAL rho, gxx, w_lorentz_2;
    CCTK_INT  star_i=0, star=0;
    /* calculate the distance to the star */
    star=get_nearest_star(x,y,z,i);
    star_i = star * TOV_Num_Radial;
    r_to_star = sqrt( (x[i]-TOV_Position_x[star]) *
                      (x[i]-TOV_Position_x[star]) +
                      (y[i]-TOV_Position_y[star]) *
                      (y[i]-TOV_Position_y[star]) +
                      (z[i]-TOV_Position_z[star]) *
                      (z[i]-TOV_Position_z[star]) );
    /* call the interpolator */
    TOV_C_interp_tov_isotropic(star,
                               &(TOV_press_1d[star_i]), &(TOV_phi_1d[star_i]),
                               &(TOV_rbar_1d[star_i]),  &(TOV_r_1d[star_i]),
                               &r_to_star, TOV_Surface[star],
                               &press, &phi, &r);
    press = (press > 0.0) ? press : 0.0;
    rho = pow(press/TOV_K, 1.0/TOV_Gamma);
    gxx = r / (r_to_star + 1.0e-30) * r / (r_to_star + 1.0e-30);
    /* velocity components as in gr-qc/9811015 */
    w_lorentz_2 = 1. / (1. - gxx*TOV_Velocity_x[star]*TOV_Velocity_x[star]
                           - gxx*TOV_Velocity_y[star]*TOV_Velocity_y[star]
                           - gxx*TOV_Velocity_z[star]*TOV_Velocity_z[star]);
    if (rho > 0.0)
      source[i]=gxx*gxx* (w_lorentz_2*rho +
                          (w_lorentz_2*TOV_Gamma/(TOV_Gamma-1.)-1.)*
                          press);
    else
      source[i]=0.0;
  }
  return 1;
}

/* This function provides an interface for other thorns which can input a list
 * of points (p_x, p_y, p_z) with a size of 'size' and get back the momentum
 * sources at these points (source) */
/* NOTE: when used with two or more NS, there might be problems if they are
 * very close to each other, because this function does not use the fancy
 * algorithms above to determine from which NS the data for a perticular grid
 * point should be taken, but just uses the closest one */
CCTK_INT TOV_Set_Momentum_Source(
                         CCTK_POINTER_TO_CONST cctkGH,
                         CCTK_INT  dir,
                         CCTK_INT  size,
                         CCTK_REAL *source,
                         const CCTK_REAL *x,
                         const CCTK_REAL *y,
                         const CCTK_REAL *z )
{
  DECLARE_CCTK_PARAMETERS
  CCTK_REAL r_to_star, press, phi, r, gxx, rho, rho_eps, psip, w_lorentz_2;
  CCTK_INT  star_i=0, star=0;
  CCTK_INT  i;

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);

  /* CCTK_INFO("Request for Momentum-Sources from TOVSolver"); */
  /* loop over all points we got */

  for (i=0; i<size; i++)
  {
    /* calculate the distance to the star */
    star=get_nearest_star(x,y,z,i);
    star_i = star * TOV_Num_Radial;
    r_to_star = sqrt( (x[i]-TOV_Position_x[star]) *
                      (x[i]-TOV_Position_x[star]) +
                      (y[i]-TOV_Position_y[star]) *
                      (y[i]-TOV_Position_y[star]) +
                      (z[i]-TOV_Position_z[star]) *
                      (z[i]-TOV_Position_z[star]) );
    /* call the interpolator */
    TOV_C_interp_tov_isotropic(star,
                               &(TOV_press_1d[star_i]), &(TOV_phi_1d[star_i]),
                               &(TOV_rbar_1d[star_i]),  &(TOV_r_1d[star_i]),
                               &r_to_star, TOV_Surface[star],
                               &press, &phi, &r);
    rho = pow(press/TOV_K, 1.0/TOV_Gamma);
    gxx = r / (r_to_star + 1.0e-30) * r / (r_to_star + 1.0e-30);
    psip = pow(gxx, TOV_Momentum_Psi_Power/4.);
    /* velocity components as in gr-qc/9811015 */
    w_lorentz_2 = 1./(1. - gxx*TOV_Velocity_x[star]*TOV_Velocity_x[star]
                         - gxx*TOV_Velocity_y[star]*TOV_Velocity_y[star]
                         - gxx*TOV_Velocity_z[star]*TOV_Velocity_z[star]);
    rho_eps = w_lorentz_2 * (rho + TOV_Gamma/(TOV_Gamma-1.) * press);
    switch(dir)
    {
      case 0: source[i]=psip*rho_eps*gxx*TOV_Velocity_x[star]; break;
      case 1: source[i]=psip*rho_eps*gxx*TOV_Velocity_y[star]; break;
      case 2: source[i]=psip*rho_eps*gxx*TOV_Velocity_z[star]; break;
      default: CCTK_WARN(0,
                         "dir has to be in [0:2] in TOV_Set_Momentum_Source");
    }
    /* Only set velocity where matter is */
    if (press<1.e-15)
      source[i]=0.0;
  }
  return 1;
}

/* NOTE: when used with two or more NS, there might be problems if they are
 * very close to each other, because this function does not use the fancy
 * algorithms above to determine from which NS the data for a perticular grid
 * point should be taken, but just uses the closest one */
CCTK_INT TOV_Set_Initial_Guess_for_u(
                         CCTK_POINTER_TO_CONST cctkGH,
                         CCTK_INT  size,
                         CCTK_REAL *u,
                         const CCTK_REAL *x,
                         const CCTK_REAL *y,
                         const CCTK_REAL *z )
{
  DECLARE_CCTK_PARAMETERS
  CCTK_REAL r_to_star, press, phi, r;
  CCTK_INT  star_i=0, star=0;
  CCTK_INT  i;

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);

  CCTK_INFO("Using initial guess from TOVSolver");
  /* TOV_debug_input_points(size,x,y,z); */

  /* loop over all points we got */
  for (i=0; i<size; i++)
  {
    CCTK_REAL gxx;
    /* calculate the distance to the stars and choose the closest */
    star=get_nearest_star(x,y,z,i);
    star_i = star * TOV_Num_Radial;
    r_to_star = sqrt( (x[i]-TOV_Position_x[star]) *
                      (x[i]-TOV_Position_x[star]) +
                      (y[i]-TOV_Position_y[star]) *
                      (y[i]-TOV_Position_y[star]) +
                      (z[i]-TOV_Position_z[star]) *
                      (z[i]-TOV_Position_z[star]) );
    /* call the interpolator */
    TOV_C_interp_tov_isotropic(star,
                               &(TOV_press_1d[star_i]), &(TOV_phi_1d[star_i]),
                               &(TOV_rbar_1d[star_i]),  &(TOV_r_1d[star_i]),
                               &r_to_star, TOV_Surface[star],
                               &press, &phi, &r);
    gxx = r / (r_to_star + 1.0e-30) * r / (r_to_star + 1.0e-30);
    u[i]=pow(gxx, 0.25) - 1.0;
  }
  return 1;
}


void bisection     (CCTK_REAL *vx, CCTK_REAL *vy, CCTK_REAL *vz,
                    CCTK_REAL *rhoold, CCTK_REAL *rhonew,
                    CCTK_REAL *w_lorentz_2, CCTK_REAL *v_2,
                    CCTK_REAL Gamma, CCTK_REAL K,
                    CCTK_REAL source,
                    CCTK_REAL mom_source_x,
                    CCTK_REAL mom_source_y,
                    CCTK_REAL mom_source_z,
                    CCTK_REAL my_psi4,
                    CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL rho_orig);
void bisection     (CCTK_REAL *vx, CCTK_REAL *vy, CCTK_REAL *vz,
                    CCTK_REAL *rhoold, CCTK_REAL *rhonew,
                    CCTK_REAL *w_lorentz_2, CCTK_REAL *v_2,
                    CCTK_REAL Gamma, CCTK_REAL K,
                    CCTK_REAL source,
                    CCTK_REAL mom_source_x,
                    CCTK_REAL mom_source_y,
                    CCTK_REAL mom_source_z,
                    CCTK_REAL my_psi4,
                    CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL rho_orig)
{
    CCTK_REAL lb, ub, f;
    lb=0.0;
    ub=2*rho_orig;
    do
    {
        *rhonew = (ub-lb)/2.;
        *vx = mom_source_x/my_psi4/
              (source/my_psi4/my_psi4+K * pow(*rhonew, Gamma));
        *vy = mom_source_y/my_psi4/
              (source/my_psi4/my_psi4+K * pow(*rhonew, Gamma));
        *vz = mom_source_z/my_psi4/
              (source/my_psi4/my_psi4+K * pow(*rhonew, Gamma));
        *v_2 = my_psi4*(*vx**vx+*vy**vy+*vz**vz);
        /* velocity components as in gr-qc/9811015 */
        *w_lorentz_2 = 1./(1.- *v_2);
        f = my_psi4*my_psi4*(
            *w_lorentz_2 * *rhonew +
            (*w_lorentz_2*Gamma/(Gamma-1.)-1.) *
            K * pow(*rhonew, Gamma)) - source;
        if (f < 0.0)
          ub = *rhonew;
        else
          lb = *rhonew;
        *rhoold = *rhonew;
    } while (ub-lb < 1.e-10);
}

void newton_raphson(CCTK_REAL *vx, CCTK_REAL *vy, CCTK_REAL *vz,
                    CCTK_REAL *rhoold, CCTK_REAL *rhonew,
                    CCTK_REAL *w_lorentz_2, CCTK_REAL *v_2,
                    CCTK_REAL Gamma, CCTK_REAL K,
                    CCTK_REAL source,
                    CCTK_REAL mom_source_x,
                    CCTK_REAL mom_source_y,
                    CCTK_REAL mom_source_z,
                    CCTK_REAL my_psi4,
                    CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL rho_orig);
void newton_raphson(CCTK_REAL *vx, CCTK_REAL *vy, CCTK_REAL *vz,
                    CCTK_REAL *rhoold, CCTK_REAL *rhonew,
                    CCTK_REAL *w_lorentz_2, CCTK_REAL *v_2,
                    CCTK_REAL Gamma, CCTK_REAL K,
                    CCTK_REAL source,
                    CCTK_REAL mom_source_x,
                    CCTK_REAL mom_source_y,
                    CCTK_REAL mom_source_z,
                    CCTK_REAL my_psi4,
                    CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL rho_orig)
{
    int count=0;
    CCTK_REAL d_w_lorentz_2, f, df;
      do
      {
        count++;
        *vx = mom_source_x/my_psi4/
              (source/my_psi4/my_psi4+K * pow(*rhonew, Gamma));
        *vy = mom_source_y/my_psi4/
              (source/my_psi4/my_psi4+K * pow(*rhonew, Gamma));
        *vz = mom_source_z/my_psi4/
              (source/my_psi4/my_psi4+K * pow(*rhonew, Gamma));
        *v_2 = my_psi4*(*vx**vx+*vy**vy+*vz**vz);
        if (count > 100)
        {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, \
                     "Failed to converge IVP rescaling, trying slow bisection method\n"\
                     "  at (%.3g %.3g %.3g)\n"\
                     "  orig_rho:%g rho:%g, v^2=%g, rel_diff=%g", \
                     x, y, z, rho_orig, *rhonew, *v_2,
                     fabs(*rhonew-*rhoold)/fabs(*rhonew));
          bisection(vx,vy,vz,rhoold,rhonew,w_lorentz_2,v_2,
                    Gamma, K,
                    source, mom_source_x, mom_source_y, mom_source_z,
                    my_psi4, x, y, z, rho_orig);
          return;
        }
        /* velocity components as in gr-qc/9811015 */
        *w_lorentz_2 = 1./(1.- *v_2);
        f = my_psi4*my_psi4*(
            *w_lorentz_2 * *rhonew +
            (*w_lorentz_2*Gamma/(Gamma-1.)-1.) *
            K * pow(*rhonew, Gamma)) - source;
        /* d_w_lorentz_2/drhonew */
        d_w_lorentz_2 =
             -2 * (*w_lorentz_2)*(*w_lorentz_2) * 
             (*v_2) * K*Gamma*pow(*rhonew, Gamma-1) / 
             (source/my_psi4/my_psi4 + K*pow(*rhonew, Gamma));
        df = my_psi4*my_psi4*(
             d_w_lorentz_2**rhonew + *w_lorentz_2 +
             d_w_lorentz_2*Gamma/(Gamma-1.)*
              K * pow(*rhonew, Gamma) +
             (*w_lorentz_2*Gamma/(Gamma-1.)-1.) *
               K * Gamma * pow(*rhonew, Gamma-1.));
        *rhoold = *rhonew;
        *rhonew = *rhoold - f/df;
      }
      while ( (fabs(*rhonew - *rhoold)/fabs(*rhonew) > 1.e-12) &&
              (fabs(*rhonew - *rhoold) > 1.e-12)
            );
}

CCTK_INT TOV_Rescale_Sources(
    CCTK_POINTER_TO_CONST cctkGH,
    CCTK_INT              size,
    const CCTK_REAL       *x,
    const CCTK_REAL       *y,
    const CCTK_REAL       *z,
    const CCTK_REAL       *psi,
    const CCTK_REAL       *gxx,
    const CCTK_REAL       *gyy,
    const CCTK_REAL       *gzz,
    const CCTK_REAL       *gxy,
    const CCTK_REAL       *gxz,
    const CCTK_REAL       *gyz)
{
  DECLARE_CCTK_PARAMETERS
  CCTK_REAL *rho, *press, *eps, *w_lorentz, *vel0, *vel1, *vel2;
  CCTK_REAL *rho_p, *press_p, *eps_p, *w_lorentz_p, *vel0_p, *vel1_p, *vel2_p;
  CCTK_REAL *rho_p_p, *press_p_p, *eps_p_p, *w_lorentz_p_p,
            *vel0_p_p, *vel1_p_p, *vel2_p_p;
  CCTK_REAL *source, *mom_source[3];
  CCTK_REAL psip, my_psi4, vx, vy, vz, v_2,
            w_lorentz_2;
  int       type;
  CCTK_REAL rhoold, rhonew;
  int i3D;
  static int debug = 1;
  FILE *debugfile = NULL;

  CCTK_INFO("Rescaling Sources");

  /* get Hydro variables */
  rho      =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "HydroBase::rho");
  press    =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "HydroBase::press");
  eps      =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "HydroBase::eps");
  w_lorentz=(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "HydroBase::w_lorentz");
  vel0     =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "HydroBase::vel[0]");
  vel1     =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "HydroBase::vel[1]");
  vel2     =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "HydroBase::vel[2]");

  rho_p    =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 1, "HydroBase::rho");
  press_p  =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 1, "HydroBase::press");
  eps_p    =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 1, "HydroBase::eps");
  w_lorentz_p=(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 1, "HydroBase::w_lorentz");
  vel0_p   =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 1, "HydroBase::vel[0]");
  vel1_p   =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 1, "HydroBase::vel[1]");
  vel2_p   =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 1, "HydroBase::vel[2]");

  rho_p_p  =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 2, "HydroBase::rho");
  press_p_p=(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 2, "HydroBase::press");
  eps_p_p  =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 2, "HydroBase::eps");
  w_lorentz_p_p=(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 2, "HydroBase::w_lorentz");
  vel0_p_p =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 2, "HydroBase::vel[0]");
  vel1_p_p =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 2, "HydroBase::vel[1]");
  vel2_p_p =(CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 2, "HydroBase::vel[2]");

  /* get the old source */
  source=calloc(size, sizeof(CCTK_REAL));
  TOV_Set_Rho_ADM(cctkGH, size, source, x, y, z);
  mom_source[0]=calloc(size, sizeof(CCTK_REAL));
  mom_source[1]=calloc(size, sizeof(CCTK_REAL));
  mom_source[2]=calloc(size, sizeof(CCTK_REAL));
  TOV_Set_Momentum_Source(cctkGH, 0, size, mom_source[0], x, y, z);
  TOV_Set_Momentum_Source(cctkGH, 1, size, mom_source[1], x, y, z);
  TOV_Set_Momentum_Source(cctkGH, 2, size, mom_source[2], x, y, z);

  if (debug)
  {
    debugfile=fopen("tovdebug.dat", "w");
    fprintf(debugfile, "# 1:x 2:y 3:psi4 4:source 5-7:mom_source_?, 8,9: ham-source, 10,11: momx-source\n");
  }
  /* loop over all points we got */
  for (i3D=0; i3D<size; i3D++)
  {
    if (source[i3D]>0.0)
    {
      if (psi != NULL)
        my_psi4 = pow(psi[i3D], 4.0)*gxx[i3D];
      else
        my_psi4 = gxx[i3D];
      /* We store the original J^i in mom_source[i] */
      psip = pow(my_psi4, TOV_Momentum_Psi_Power/4.);
      mom_source[0][i3D] /= psip;
      mom_source[1][i3D] /= psip;
      mom_source[2][i3D] /= psip;
      if (debug && (fabs(z[i3D])<1.e-15))
      {
        fprintf(debugfile,
                "%.15g %.15g %.15g %.15g %.15g %.15g %.15g ",
                x[i3D], y[i3D], my_psi4, source[i3D],
                mom_source[0][i3D], mom_source[1][i3D], mom_source[1][i3D]);
      }

      if (rho[i3D] < 1.e-10)
          rho[i3D] = 1.e-10;
      rhoold = rhonew = rho[i3D];
#ifdef OLD_NR
      int count=0;
      do
      {
        CCTK_REAL f, df;
        count++;
        vx = mom_source[0][i3D]/my_psi4/
             (source[i3D]/my_psi4/my_psi4+TOV_K * pow(rhonew, TOV_Gamma));
        vy = mom_source[1][i3D]/my_psi4/
             (source[i3D]/my_psi4/my_psi4+TOV_K * pow(rhonew, TOV_Gamma));
        vz = mom_source[2][i3D]/my_psi4/
             (source[i3D]/my_psi4/my_psi4+TOV_K * pow(rhonew, TOV_Gamma));
        v_2 = my_psi4*(vx*vx+vy*vy+vz*vz);
        if (count > 100)
          CCTK_VWarn(count<110, __LINE__, __FILE__, CCTK_THORNSTRING, \
                     "Failed to converge IVP rescaling\nat (%.3g %.3g %.3g)\n"\
                     " orig_rho:%g rho:%g, v^2=%g, rel_diff=%g", \
                     x[i3D], y[i3D], z[i3D], rho[i3D], rhonew, v_2,
                     fabs(rhonew-rhoold)/fabs(rhonew));
        /* velocity components as in gr-qc/9811015 */
        w_lorentz_2 = 1./(1.- v_2);
        /*
        f = TOV_K / (TOV_Gamma-1.0) * pow(rhonew, TOV_Gamma) +
            rhonew - source[i3D];
        df = TOV_K * TOV_Gamma / (TOV_Gamma-1.0) *
             pow(rhonew, TOV_Gamma-1.0) + 1.0;
        */
        f = my_psi4*my_psi4*(
            w_lorentz_2 * rhonew +
            (w_lorentz_2*TOV_Gamma/(TOV_Gamma-1.)-1.) *
            TOV_K * pow(rhonew, TOV_Gamma)) - source[i3D];
        /*
        df= w_lorentz_2 +
            (w_lorentz_2*TOV_Gamma/(TOV_Gamma-1.)-1.) *
            TOV_K * TOV_Gamma * pow(rhonew, TOV_Gamma-1.);
        */
        /* d_w_lorentz_2/drhonew */
        CCTK_REAL d_w_lorentz_2 =
             -2 * w_lorentz_2*w_lorentz_2 *
             v_2  * TOV_K*TOV_Gamma*pow(rhonew, TOV_Gamma-1) /
             (source[i3D]/my_psi4/my_psi4 + TOV_K*pow(rhonew, TOV_Gamma));
        df = my_psi4*my_psi4*(
             d_w_lorentz_2*rhonew + w_lorentz_2 +
             d_w_lorentz_2*TOV_Gamma/(TOV_Gamma-1.)*
              TOV_K * pow(rhonew, TOV_Gamma) +
             (w_lorentz_2*TOV_Gamma/(TOV_Gamma-1.)-1.) *
               TOV_K * TOV_Gamma * pow(rhonew, TOV_Gamma-1.));
        /*
        df= 1. + my_psi4*v_2/(rhonew*rhonew) +
            ((TOV_Gamma+(2.-TOV_Gamma)*my_psi4*v_2/(rhonew*rhonew))
             / (TOV_Gamma-1.) - 1.)
            * TOV_K * TOV_Gamma * pow(rhonew, TOV_Gamma-1.);
        */

        rhoold = rhonew;
        rhonew = rhoold - f/df;
      }
      while ( (fabs(rhonew-rhoold)/fabs(rhonew) > 1.e-12) &&
              (fabs(rhonew - rhoold) > 1.e-12)
            );
      if (count==1)
          printf("Fast rescale! %g %g %g %g\n", my_psi4, my_psi4,
                 mom_source[0][i3D]/my_psi4/
                  (source[i3D]/my_psi4/my_psi4+TOV_K*
                  pow(rhonew, TOV_Gamma)),
                 mom_source[0][i3D]/my_psi4/
                  (source[i3D]/my_psi4/my_psi4+TOV_K*
                  pow(rhonew, TOV_Gamma)));
#else
      newton_raphson(&vx,&vy,&vz,&rhoold,&rhonew,&w_lorentz_2,&v_2,
                     TOV_Gamma, TOV_K,
                     source[i3D],
                     mom_source[0][i3D],mom_source[1][i3D],mom_source[2][i3D],
                     my_psi4, x[i3D], y[i3D], z[i3D], rho[i3D]);
#endif
      rho[i3D]   = rhonew;
      press[i3D] = TOV_K * pow(rhonew, TOV_Gamma);
      eps[i3D]   = TOV_K * pow(rhonew, TOV_Gamma - 1.0) /
                              (TOV_Gamma - 1.0);
      w_lorentz[i3D] = sqrt(w_lorentz_2);
      vel0[i3D] = vx;
      vel1[i3D] = vy;
      vel2[i3D] = vz;

      if (debug && (fabs(z[i3D])<1.e-15))
      {
        fprintf(debugfile,
                "%.15g %.15g %.15g %.15g\n", /*7*/
                gxx[i3D]*gxx[i3D]*( /* ham - source*/
                  w_lorentz_2*rho[i3D]+
                  (w_lorentz_2*TOV_Gamma/(TOV_Gamma-1.)-1.)*press[i3D]),
                gxx[i3D]*gxx[i3D]*( /* ham - source */
                  w_lorentz_2*(rho[i3D]*(1.+eps[i3D])+press[i3D])-press[i3D]),
                (w_lorentz_2 * rho[i3D] + /* momx - source */
                 w_lorentz_2*TOV_Gamma/(TOV_Gamma-1.) * press[i3D])*
                 gxx[i3D]*vel0[i3D],
                w_lorentz_2 * (rho[i3D]* /* momx - source */
                 (1.+eps[i3D])+ press[i3D])*
                 gxx[i3D]*vel0[i3D]
                );
      }

    }
  }
  if (debug)
    fclose(debugfile);

  switch(*(const CCTK_INT*)CCTK_ParameterGet("TOV_Populate_Timelevels",
                                             "TOVSolver", &type))
  {
    case 3:
        TOV_Copy(size, rho_p_p,  rho);
        TOV_Copy(size, press_p_p,press);
        TOV_Copy(size, eps_p_p,  eps);
        TOV_Copy(size, vel0_p_p, vel0);
        TOV_Copy(size, vel1_p_p, vel1);
        TOV_Copy(size, vel2_p_p, vel2);
        TOV_Copy(size, w_lorentz_p_p, w_lorentz);
    case 2:
        TOV_Copy(size, rho_p,  rho);
        TOV_Copy(size, press_p,press);
        TOV_Copy(size, eps_p,  eps);
        TOV_Copy(size, vel0_p, vel0);
        TOV_Copy(size, vel1_p, vel1);
        TOV_Copy(size, vel2_p, vel2);
        TOV_Copy(size, w_lorentz_p, w_lorentz);
  }
  free(source);
  return 1;
}

void TOV_write_1D_datafile(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TOV_write_1D_datafile
  DECLARE_CCTK_PARAMETERS
  int i;

  if (TOV_Num_TOVs > 1)
    CCTK_WARN(0, "Writing a data file for multiple stars is not (yet) "
                 "supported");
  FILE *file;
  file = fopen(TOV_save_to_datafile, "w");
  fprintf(file, "TOVSolver data file\n");
  fprintf(file, "version 1.0\n");
  fprintf(file, "TOV_Num_Radial %d\n", (int)TOV_Num_Radial);
  fprintf(file, "\n");
  for (i=0; i < TOV_Num_Radial; i++)
  {
      fprintf(file, "%g %g %g %g %g\n",
              TOV_rbar_1d[i],
              TOV_r_1d[i],
              TOV_m_1d[i],
              TOV_phi_1d[i],
              TOV_press_1d[i]);
  }
  fclose(file);
  CCTK_WARN(0, "Simulation stopped as requested after writing data to file");
}

void TOVSolverC_export_local_variables(CCTK_REAL **exported_TOV_press_1d, CCTK_REAL **exported_TOV_m_1d, CCTK_REAL **exported_TOV_phi_1d, CCTK_REAL **exported_TOV_rbar_1d, CCTK_REAL **exported_TOV_r_1d)
{
  *exported_TOV_press_1d = TOV_press_1d;
  *exported_TOV_m_1d     = TOV_m_1d;
  *exported_TOV_phi_1d   = TOV_phi_1d;
  *exported_TOV_rbar_1d  = TOV_rbar_1d;
  *exported_TOV_r_1d     = TOV_r_1d;
}

