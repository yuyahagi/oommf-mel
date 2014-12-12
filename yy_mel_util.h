/* FILE: yy_MEL_util.h                 -*-Mode: c++-*-
 *
 * Misc
 * 
 */

#ifndef _YY_MEL_UTIL
#define _YY_MEL_UTIL

#include <string>
#include <vector>

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
#include "output.h"

/* End includes */

class YY_Strain {
private:
  Oxs_MeshValue<ThreeVector> u;       // Displacement
  Oxs_MeshValue<OC_REAL8m> Ms;
  OC_INDEX nx, ny, nz;

  OC_BOOL displacement_valid;
  OC_BOOL strain_valid; // True if strain has been calculated

  //mutable Oxs_OwnedPointer<Oxs_VectorField> diag_init;

  Nb_TclCommand cmd;

public:
  //virtual const char* ClassName() const;
  //// ClassName() is automatically generatoed by OXS_EXT_REGISTER macro.
  YY_Strain();
  ~YY_Strain() {};

  Oxs_MeshValue<ThreeVector> diag;    // Diagonal elements of strain
  Oxs_MeshValue<ThreeVector> offdiag; // Off-diagonal elements
  // in abbreviated suffix notation
  //  Strain      Diagonal   Off-diagonal
  // / 0 5 4 \   / 0     \   /   2 1 \
  // | . 1 3 | = |   1   | + |     0 |
  // \ . . 2 /   \     2 /   \       /
  // offdiag.x = e_yz ([1][2])
  // offdiag.y = e_xz ([0][2])
  // offdiag.z = e_xy ([0][1])

  void Init(const Oxs_SimState& state);
  void SetDisplacement(const Oxs_SimState& state, const Oxs_MeshValue<ThreeVector>& u_in);
  //void SetDisplacement(const Oxs_MeshValue<ThreeVector>& u_in)
  //{ u = u_in; displacement_valid = 1; strain_valid = 0;};

  void CalculateStrain(const Oxs_SimState&);

  // For debug. Display variables in specified range
  void DisplayValues(
      const Oxs_SimState& state,
      OC_INDEX xmin, OC_INDEX xmax,
      OC_INDEX ymin, OC_INDEX ymax,
      OC_INDEX zmin, OC_INDEX zmax) const;

};

#endif // _YY_MEL_UTIL
