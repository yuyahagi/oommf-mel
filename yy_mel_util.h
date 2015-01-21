/* FILE: yy_MEL_util.h                 -*-Mode: c++-*-
 *
 * Misc
 * 
 */

//#define YY_DEBUG
#ifdef YY_DEBUG
#define YY_DEBUGMSG(x) fprintf(stderr,x)
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
  Oxs_MeshValue<OC_REAL8m> MELCoef;   // Isotropic MEL coefficient

  mutable Oxs_MeshValue<ThreeVector> u;       // Displacement
  mutable Oxs_MeshValue<ThreeVector> diag;    // Diagonal elements of strain
  mutable Oxs_MeshValue<ThreeVector> offdiag; // Off-diagonal elements
  // in abbreviated suffix notation
  //  Strain      Diagonal   Off-diagonal
  // / 0 5 4 \   / 0     \   /   2 1 \
  // | . 1 3 | = |   1   | + |     0 |
  // \ . . 2 /   \     2 /   \       /
  // offdiag.x = e_yz ([1][2])
  // offdiag.y = e_xz ([0][2])
  // offdiag.z = e_xy ([0][1])

  OC_BOOL strain_valid; // True if strain has been calculated
  OC_BOOL MELCoef_valid;

  mutable ThreeVector max_field;

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
      const Oxs_OwnedPointer<Oxs_ScalarField>& MELCoef_init);
  void SetDisplacement(
      const Oxs_SimState& state,
      const Oxs_OwnedPointer<Oxs_VectorField>& u_init);
  void SetStrain(
      const Oxs_SimState& state,
      const Oxs_OwnedPointer<Oxs_VectorField>& diag_init,
      const Oxs_OwnedPointer<Oxs_VectorField>& offdiag_init);

  void CalculateMELField(
    const Oxs_SimState& state,
    OC_REAL8m hmult,
    Oxs_MeshValue<ThreeVector>& field_buf) const;

  ThreeVector GetMaxField() const
  { return max_field; };

  // For debug. Display variables in specified range
  void DisplayValues(
      const Oxs_SimState& state,
      OC_INDEX xmin, OC_INDEX xmax,
      OC_INDEX ymin, OC_INDEX ymax,
      OC_INDEX zmin, OC_INDEX zmax) const;

};

#endif // _YY_MEL_UTIL
