#ifndef CARPET_PROLONGATE_3D_RF2_H_
#define CARPET_PROLONGATE_3D_RF2_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Interpolater.H>

/**
* \brief Vertex centered sampling interpolation.
*
* Vertex centered interpolation copied from CarpetLib.
*/

namespace amrex {
template <int ORDER>
class Prolongate3D_VC_RF2
    :
    public Interpolater
{
public:

    /**
    * \brief The destructor.
    */
    virtual ~Prolongate3D_VC_RF2 () override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            gpu_or_cpu) override;
};

extern Prolongate3D_VC_RF2<5> prolongate3d_vc_rf2_o5;
}

#endif // CARPET_PROLONGATE_3D_RF2_H_
