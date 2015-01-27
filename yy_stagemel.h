/* FILE: yy_stagemel.h                 -*-Mode: c++-*-
 *
 * OOMMF magnetoelastic coupling extension module.
 * YY_StageMEL class.
 * Calculates MEL field/energy for each stage.
 * 
 */

#ifndef _YY_STAGEMEL
#define _YY_STAGEMEL

#include <string>
#include <vector>

#include "nb.h"
#include "director.h"
#include "energy.h"
#include "simstate.h"
#include "threevector.h"
#include "vectorfield.h"

#include "yy_mel_util.h"

OC_USE_STRING;

/* End includes */

class YY_StageMEL:public Oxs_Energy {
private:
  // magnetoelastic coefficient MELCoef1,2 = B1,2 (J/m^3)
  Oxs_OwnedPointer<Oxs_ScalarField> MELCoef1_init, MELCoef2_init;
  mutable Oxs_OwnedPointer<Oxs_VectorField> u_init;
  mutable Oxs_OwnedPointer<Oxs_VectorField> e_diag_init, e_offdiag_init;

  mutable YY_MELField MELField;

  // The following members were cribbed from stagezeeman.h
  OC_REAL8m hmult;  // multiplier
  OC_UINT4m number_of_stages;

  mutable OC_UINT4m mesh_id;
  mutable OC_BOOL stage_valid;
  mutable OC_UINT4m working_stage;  // Stage index
  mutable ThreeVector max_field;

  // Vector field may be specified by *either* a list of
  // files, or else a Tcl command that returns a vector
  // field spec.  This is set up in the Oxs_StageZeeman
  // constructor.  The choice is determined elsewhere by
  // examining the length of filelist; if it is >0, then
  // the list of files method is used, otherwise cmd is
  // called with the stage number as the single appended
  // argument.
  vector<String> u_filelist;
  vector<String> e_diag_filelist, e_offdiag_filelist;
  Nb_TclCommand u_cmd, e_diag_cmd, e_offdiag_cmd;

  // Flags to specify the input source(s) of strain.
  // Displacement u or strain e and with scripts or
  // file lists as a shorhand.
  OC_BOOL use_u, use_u_filelist, use_u_script;
  OC_BOOL use_e, use_e_filelist, use_e_script;

  // Member function to determine the way of specifying
  // strain. Whether directly using strain or indirectly by
  // specifying displacement and calculating strain with it.
  // Also determines whether filelists or scripts are used.
  // It sets use_u and use_e (== !use_u) and use_u_filelist
  // or use_e_filelist depending on use_u.
  void SelectElasticityInputType();

  // Update the initializers for displacement or strain. ChangeInitializer()
  // checks use_u flag, calls either one of ChangeDisplacementInitializer()
  // or ChangeStrainInitializer(). After that, call SetStrain() to sets the
  // MELField strain in a proper way by calling either 
  // MELField.SetDisplacement() or MELField.SetStrain().
  void ChangeInitializer(const Oxs_SimState& state) const;
  void ChangeDisplacementInitializer(OC_UINT4m stage, const Oxs_Mesh* mesh) const;
  void ChangeStrainInitializer(OC_UINT4m stage, const Oxs_Mesh* mesh) const;
  void SetStrain(const Oxs_SimState& state) const;
  void UpdateCache(const Oxs_SimState& state) const;

  // Additional outputs
  Oxs_ScalarOutput<YY_StageMEL> B_MEL_output;
  Oxs_ScalarOutput<YY_StageMEL> B_MELx_output;
  Oxs_ScalarOutput<YY_StageMEL> B_MELy_output;
  Oxs_ScalarOutput<YY_StageMEL> B_MELz_output;
  void Fill__B_MEL_output(const Oxs_SimState& state);

protected:
  virtual void GetEnergy(const Oxs_SimState& state,
      Oxs_EnergyData& oed) const;

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  YY_StageMEL(const char* name, // Child instance id
		    Oxs_Director* newdtr,   // App director
		    const char* argstr);    // MIF input block parameters
  virtual ~YY_StageMEL();
  virtual OC_BOOL Init();
  virtual void StageRequestCount(unsigned int& min,
      unsigned int& max) const;

#ifdef YY_DEBUG
  // For debug. Display variables in specified range
  void DisplayValues(
      const Oxs_SimState& state,
      OC_INDEX xmin, OC_INDEX xmax,
      OC_INDEX ymin, OC_INDEX ymax,
      OC_INDEX zmin, OC_INDEX zmax) const;
#endif
};

#endif // _YY_STAGEMEL
