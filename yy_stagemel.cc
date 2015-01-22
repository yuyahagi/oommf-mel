/* FILE: yy_stagemel.cc                 -*-Mode: c++-*-
 *
 * OOMMF magnetoelastic coupling extension module.
 * YY_StageMEL class.
 * Calculates MEL field/energy for each stage.
 * 
 */

#include "oc.h"
#include "nb.h"
#include "director.h"
#include "mesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "rectangularmesh.h"
#include "energy.h"             // Needed to make MSVC++ 5 happy

#include "yy_stagemel.h"

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(YY_StageMEL);

/* End includes */

// Constructor
YY_StageMEL::YY_StageMEL(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_Energy(name,newdtr,argstr), number_of_stages(0),
  mesh_id(0), stage_valid(0)
{
  YY_DEBUGMSG("YY_StageMEL constructor start.\n");
  working_stage = static_cast<OC_UINT4m>(static_cast<OC_INT4m>(-1));
  // First ChangeFileInitializer() call should be triggered
  // by the stage_valid boolean.

  // Process arguments
  hmult = GetRealInitValue("multiplier", 1.0);

  // Script name
  cmd.SetBaseCommand(InstanceName(),
                     director->GetMifInterp(),
                     GetStringInitValue("script"),1);
  number_of_stages = GetUIntInitValue("stage_count",0);

  // Generate MELCoef initializer
  if(HasInitValue("B1")) {
    OXS_GET_INIT_EXT_OBJECT("B1",Oxs_ScalarField,MELCoef1_init);
  } else {
    MELCoef1_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
        (MakeNew("Oxs_UniformScalarField",director,"value 7.85e6")));
  }
  if(HasInitValue("B2")) {
    OXS_GET_INIT_EXT_OBJECT("B2",Oxs_ScalarField,MELCoef2_init);
  } else {
    MELCoef2_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
        (MakeNew("Oxs_UniformScalarField",director,"value 7.85e6")));
  }

  // Initialize outputs.
  YY_DEBUGMSG("YY_StageMEL constructor initialize outputs.\n");
  B_MEL_output.Setup(this, InstanceName(), "B max", "mT", 1,
    &YY_StageMEL::Fill__B_MEL_output);
  B_MELx_output.Setup(this, InstanceName(), "Bx max", "mT", 1,
    &YY_StageMEL::Fill__B_MEL_output);
  B_MELy_output.Setup(this, InstanceName(), "By max", "mT", 1,
    &YY_StageMEL::Fill__B_MEL_output);
  B_MELz_output.Setup(this, InstanceName(), "Bz max", "mT", 1,
    &YY_StageMEL::Fill__B_MEL_output);

  // Register outputs.
  B_MEL_output.Register(director, 0);
  B_MELx_output.Register(director, 0);
  B_MELy_output.Register(director, 0);
  B_MELz_output.Register(director, 0);

  VerifyAllInitArgsUsed();
}

YY_StageMEL::~YY_StageMEL()
{}

OC_BOOL YY_StageMEL::Init()
{
  YY_DEBUGMSG("YY_StageMEL::Init(): start.\n");
  stage_valid = 0;

  mesh_id = 0;
  MELField.Release();

  YY_DEBUGMSG("YY_StageMEL::Init(): Calling Oxs_Energy::Init()\n");
  return Oxs_Energy::Init();
}

void YY_StageMEL::StageRequestCount(unsigned int& min,
    unsigned int& max) const
{
  if(number_of_stages == 0) {
    min = 0; max = UINT_MAX;  // No restriction on stage count
  } else {
    min = max = number_of_stages;
  }
}

void YY_StageMEL::ChangeDisplacementInitializer(
    OC_UINT4m stage, const Oxs_Mesh* mesh) const
{
  YY_DEBUGMSG("YY_StageMEL::ChangeDisplacementInitializer(): start.\n");
  // Setup displacement
  vector<String> params;
  if(filelist.empty()) {
    // Use cmd to generate field initializer
    YY_DEBUGMSG("YY_StageMEL::ChangeDisplacementInitializer(): filelist.empty.\n");
    cmd.SaveInterpResult();
    cmd.SetCommandArg(0,stage);
    cmd.Eval();
    cmd.GetResultList(params);
    cmd.RestoreInterpResult();
  } else {
    YY_DEBUGMSG("YY_StageMEL::ChangeDisplacementInitializer(): !filelist.empty.\n");
    // Construct field initializer using Oxs_FileVectorField
    // with filename from filelist and range from mesh.
    OC_UINT4m index = stage;
    OC_UINT4m filecount = static_cast<OC_UINT4m>(filelist.size());
    if(index >= filecount) index = filecount - 1;
    vector<String> options;
    options.push_back(String("file"));
    options.push_back(filelist[index]);

    Oxs_Box bbox;    mesh->GetBoundingBox(bbox);
    char buf[64];
    options.push_back(String("xrange"));
    Oc_Snprintf(buf,sizeof(buf),"%.17g %.17g",
                static_cast<double>(bbox.GetMinX()),
                static_cast<double>(bbox.GetMaxX()));
    options.push_back(String(buf));

    options.push_back(String("yrange"));
    Oc_Snprintf(buf,sizeof(buf),"%.17g %.17g",
                static_cast<double>(bbox.GetMinY()),
                static_cast<double>(bbox.GetMaxY()));
    options.push_back(String(buf));

    options.push_back(String("zrange"));
    Oc_Snprintf(buf,sizeof(buf),"%.17g %.17g",
                static_cast<double>(bbox.GetMinZ()),
                static_cast<double>(bbox.GetMaxZ()));
    options.push_back(String(buf));

    params.push_back(String("Oxs_FileVectorField"));
    params.push_back(Nb_MergeList(options));
  }

  //YY_DEBUGMSG("YY_StageMEL::ChangeDisplacementInitializer(): OXS_GET_EXT_OBJECT, displacement_init.\n");
  OXS_GET_EXT_OBJECT(params,Oxs_VectorField,displacement_init);
  working_stage = stage;
  stage_valid = 1;
}

void
YY_StageMEL::FillDisplacementCache(const Oxs_SimState& state) const
{
  const Oxs_Mesh* mesh = state.mesh;
  // Displacement
  YY_DEBUGMSG("YY_StageMEL::FillDisplacementCache(): FillMeshValue.\n");
  MELField.SetDisplacement(state, displacement_init);

}

void YY_StageMEL::UpdateCache(const Oxs_SimState& state) const
{
  // Update cache as necessary
  if(!stage_valid) {
    // The first go.
    YY_DEBUGMSG("YY_StageMEL:UpdateCache(): !stage_valide.\n");
    mesh_id = 0;
    MELField.SetMELCoef(state,MELCoef1_init,MELCoef2_init);
    ChangeDisplacementInitializer(state.stage_number,state.mesh);
    FillDisplacementCache(state);
  } else if(working_stage != state.stage_number) {
    mesh_id = 0;
    YY_DEBUGMSG("YY_StageMEL:UpdateCache(): ChangeDisplacementInitializer().\n");
    ChangeDisplacementInitializer(state.stage_number,state.mesh);
    YY_DEBUGMSG("YY_StageMEL:UpdateCache(): FillDisplacementCache().\n");
    FillDisplacementCache(state);
    mesh_id = state.mesh->Id();
  } else if(mesh_id != state.mesh->Id()) {
    mesh_id = 0;
    YY_DEBUGMSG("YY_StageMEL:UpdateCache(): in else, FillDisplacementCache().\n");
    FillDisplacementCache(state);
    mesh_id = state.mesh->Id();
  }
}


void YY_StageMEL::GetEnergy
(const Oxs_SimState& state, Oxs_EnergyData& oed) const
{
  YY_DEBUGMSG("YY_StageMEL::GetEnergy(): Beginning the method.\n");
  OC_INDEX size = state.mesh->Size();
  if(size<1) return;

  YY_DEBUGMSG("YY_StageMEL::GetEnergy(): UpdateCache().\n");
  UpdateCache(state);

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);

  YY_DEBUGMSG("YY_StageMEL::GetEnergy(): Starting H_MEL calculation.\n");

  // Use supplied buffer space, and reflect that use in oed.
  oed.energy = oed.energy_buffer;
  oed.field = oed.field_buffer;
  Oxs_MeshValue<OC_REAL8m>& energy = *oed.energy_buffer;
  Oxs_MeshValue<ThreeVector>& field = *oed.field_buffer;
  energy.AdjustSize(state.mesh);
  field.AdjustSize(state.mesh);

  MELField.CalculateMELField(state, hmult, field);
  max_field = MELField.GetMaxField();

  // Calculate pointwise energy density -0.5*MU0*<M,Hmel>
  for(OC_INDEX i=0; i<size; i++) {
    energy[i] = -0.5 * MU0 * Ms[i] * (spin[i]*field[i]);
  }

}

void
YY_StageMEL::Fill__B_MEL_output(const Oxs_SimState& state)
{
  UpdateCache(state);

  ThreeVector B = max_field;
  B *= MU0 * 1000;  // Report B_MEL in mT

  if(B_MEL_output.GetCacheRequestCount() > 0) {
    B_MEL_output.cache.state_id = 0;
    B_MEL_output.cache.value = sqrt(B.MagSq());
    B_MEL_output.cache.state_id = state.Id();
  }

  if(B_MELx_output.GetCacheRequestCount() > 0) {
    B_MELx_output.cache.state_id = 0;
    B_MELx_output.cache.value = sqrt(B.MagSq());
    B_MELx_output.cache.state_id = state.Id();
  }

  if(B_MELy_output.GetCacheRequestCount() > 0) {
    B_MELy_output.cache.state_id = 0;
    B_MELy_output.cache.value = sqrt(B.MagSq());
    B_MELy_output.cache.state_id = state.Id();
  }

  if(B_MELz_output.GetCacheRequestCount() > 0) {
    B_MELz_output.cache.state_id = 0;
    B_MELz_output.cache.value = sqrt(B.MagSq());
    B_MELz_output.cache.state_id = state.Id();
  }
}
