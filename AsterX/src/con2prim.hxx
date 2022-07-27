#ifndef CON2PRIM_H
#define CON2PRIM_H
#include "AMReX_GpuQualifiers.H"
#define CCTK_DEVICE AMREX_GPU_DEVICE
#define CCTK_HOST AMREX_GPU_HOST
#include "utils.hxx"

namespace AsterX
{

constexpr int NCONS = 8;
constexpr int NPRIMS = 8;

constexpr int D = 0;
constexpr int S1_COV = 1;
constexpr int S2_COV = 2;
constexpr int S3_COV = 3;
constexpr int TAU = 4;
constexpr int B1 = 5;
constexpr int B2 = 6;
constexpr int B3 = 7;


constexpr int RHO = 0;
constexpr int V1_CON = 1;
constexpr int V2_CON = 2;
constexpr int V3_CON = 3;
constexpr int V1_COV = 1;
constexpr int V2_COV = 2;
constexpr int V3_COV = 3;
constexpr int EPS = 4;

constexpr int TT = 0;
constexpr int TX = 1;
constexpr int TY = 2;
constexpr int TZ = 3;
constexpr int XX = 4;
constexpr int XY = 5;
constexpr int XZ = 6;
constexpr int YY = 7;
constexpr int YZ = 8;
constexpr int ZZ = 9;

/* Abstract class con2primFactory */
class con2primFactory
{
public:
    /* The constructor must initialize the following vectors */
    //std::vector<PrimitiveLabel> PrimitiveLabels;
    CCTK_REAL ConservedVars[NCONS];            // Conserved to solve
    CCTK_REAL PrimitiveVarsSeed[NPRIMS]; // Primitive seeds
    CCTK_REAL PrimitiveVars[NPRIMS];     // Primitive solution

    /* These must be set for 2DNRNoble scheme */
    CCTK_INT Failed_2DNRNoble;
    CCTK_INT Nit_2DNRNoble; 
    CCTK_REAL W_Seed, vsq_Sol, Ssq, Press_Seed, Z_Seed, Z_Sol, vsq_Seed, bsq;
    CCTK_REAL gcov[10], gcon[10];
    CCTK_REAL Bsq;
    CCTK_REAL BiSi;
    CCTK_HOST CCTK_DEVICE void get_Ssq_Exact();          // From cons (exact)
    CCTK_HOST CCTK_DEVICE void get_Press_Seed();         // From seed prims and cons
    CCTK_HOST CCTK_DEVICE void get_Z_Seed();             // From seed prims and cons
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq); 
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE void WZ2Prim();
};


class idealFluid : public con2primFactory
{
public:
    /* Some attributes */
    CCTK_REAL GammaIdealFluid;
    /* Constructor */
    CCTK_HOST CCTK_DEVICE idealFluid(CCTK_REAL gamma, CCTK_REAL (&cons)[NCONS], CCTK_REAL (&prim)[NPRIMS], CCTK_REAL (&gcov)[4][4], CCTK_REAL (&gcon)[4][4]);

    /* Called by 2DNRNoble */
    CCTK_HOST CCTK_DEVICE void get_Ssq_Exact();
    CCTK_HOST CCTK_DEVICE void get_Press_Seed();
    CCTK_HOST CCTK_DEVICE void get_Z_Seed();
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE void WZ2Prim();

    /* Destructor */
    CCTK_HOST CCTK_DEVICE ~idealFluid();
};
}
#endif // #ifndef CON2PRIM_H

