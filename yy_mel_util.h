/* FILE: yy_MEL_util.h                 -*-Mode: c++-*-
 *
 * OOMMF magnetoelastic coupling extension module.
 * yy_mel_util.* contain YY_MELField class , which is to be used by other
 * YY_*MEL classes
 * 
 */

#define YY_DEBUG
#ifdef YY_DEBUG
#define YY_DEBUGMSG(x) fprintf(stderr,(x))
#else
#define YY_DEBUGMSG(x) ;
#endif

#ifndef _YY_MEL_UTIL
#define _YY_MEL_UTIL

#include <string>

#include "oc.h"
#include "director.h"
#include "energy.h"
#include "mesh.h"
#include "rectangularmesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "vectorfield.h"
#include "util.h"

/* End includes */
class YY_MELField {
private:
  Oxs_MeshValue<OC_REAL8m> MELCoef1, MELCoef2;   // MEL coefficients

  // Displacement u and strain e.
  // Displacement
  mutable Oxs_MeshValue<ThreeVector> u, u_cache;
  mutable Oxs_MeshValue<ThreeVector> du;  // du/dt
  // Diagonal elements of strain
  mutable Oxs_MeshValue<ThreeVector> e_diag, e_diag_cache;
  mutable Oxs_MeshValue<ThreeVector> de_diag, de_diag_cache;
  // Off-diagonal elements
  mutable Oxs_MeshValue<ThreeVector> e_offdiag, e_offdiag_cache;
  mutable Oxs_MeshValue<ThreeVector> de_offdiag, de_offdiag_cache;
  // The *_cache stores the displacement/strain when they are
  // updated at each stage. When TransformDisplacement() or
  // TransformStrain() is called, transformed values are stored
  // in u or e_*diag and are used for MEL field calculation.
  // Strain is expressed with abbreviated suffix notation
  //  Strain      Diagonal   Off-diagonal
  // / 0 5 4 \   / 0     \   /   2 1 \
  // | . 1 3 | = |   1   | + |     0 |
  // \ . . 2 /   \     2 /   \       /
  // e.offdiag.x = e_yz ([1][2])
  // e.offdiag.y = e_xz ([0][2])
  // e.offdiag.z = e_xy ([0][1])

  OC_BOOL displacement_valid; // True if u_cache has been calculated
  OC_BOOL strain_valid;       // True if e_*diag_cache has been calculated
  OC_BOOL MELCoef1_valid, MELCoef2_valid;

  mutable ThreeVector max_field;
  mutable OC_REAL8m pE_pt_buf;

public:
  YY_MELField();
  ~YY_MELField() {};
  void Release();

  // Note: Before calculation of MEL field, the MEL coefficient and the
  // strain must be set with the following member functions. The strain can
  // be set either directly with SetStrain or by specifying the displacement
  // with SetDisplacement, which calculates the strain from the displacement.
  void SetMELCoef(
      const Oxs_SimState& state,
      const Oxs_OwnedPointer<Oxs_ScalarField>& MELCoef1_init,
      const Oxs_OwnedPointer<Oxs_ScalarField>& MELCoef2_init);
  void SetDisplacement(
      const Oxs_SimState& state,
      const Oxs_OwnedPointer<Oxs_VectorField>& u_init);
  void CalculateStrain(const Oxs_SimState& state);
  void SetStrain(
      const Oxs_SimState& state,
      const Oxs_OwnedPointer<Oxs_VectorField>& e_diag_init,
      const Oxs_OwnedPointer<Oxs_VectorField>& e_offdiag_init);

  // Note: The following member functions are for transforming the
  // displacement u or strain e_diag, e_offdiag with the transformation
  // matrix T, which is defined with three row vectors row1, row2, row3
  // and its time derivative drow1, drow2, drow3.
  // Note that T is applied differently for displacement u and strain e.
  // The transformed value will be
  // u = T u_cache
  // e = T e_cache T_t
  // where T_t is the transpose of T.
  void TransformDisplacement(
      const Oxs_SimState& state,
      ThreeVector& row1, ThreeVector& row2, ThreeVector& row3,
      ThreeVector& drow1, ThreeVector& drow2, ThreeVector& drow3);
  void TransformStrain(
      const Oxs_SimState& state,
      ThreeVector& row1, ThreeVector& row2, ThreeVector& row3,
      ThreeVector& drow1, ThreeVector& drow2, ThreeVector& drow3);

  void CalculateMELField(
    const Oxs_SimState& state,
    OC_REAL8m hmult,
    Oxs_MeshValue<ThreeVector>& field_buf,
    Oxs_MeshValue<OC_REAL8m>& energy_buf) const;

  ThreeVector GetMaxField() const
  { return max_field; };
  OC_REAL8m GetPEPT() const
  { return pE_pt_buf; };

#ifdef YY_DEBUG
  // For debug. Display variables in specified range
  void DisplayValues(
      const Oxs_SimState& state,
      OC_INDEX xmin, OC_INDEX xmax,
      OC_INDEX ymin, OC_INDEX ymax,
      OC_INDEX zmin, OC_INDEX zmax) const;
#endif
};

#endif // _YY_MEL_UTIL
