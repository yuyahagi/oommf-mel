/* FILE: yystageisotropicmel.cc                 -*-Mode: c++-*-
 *
 * Magnetoelastic energy with isotropic coefficient
 * 
 */

#include "nb.h"
#include "director.h"
#include "mesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "rectangularmesh.h"
#include "energy.h"             // Needed to make MSVC++ 5 happy

#include "yy_stageisotropicmel.h"

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(YY_StageIsotropicMEL);

/* End includes */

// Constructor
YY_StageIsotropicMEL::YY_StageIsotropicMEL(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_Energy(name,newdtr,argstr), number_of_stages(0),
  mesh_id(0), stage_valid(0),
  stagefield(NULL)
{
  fprintf(stderr, "YY_StageIsotropicMEL constructor start.\n");
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
  if(HasInitValue("B")) {
    OXS_GET_INIT_EXT_OBJECT("B",Oxs_ScalarField,MELCoef_init);
  } else {
    MELCoef_init.SetAsOwner(dynamic_cast<Oxs_ScalarField *>
        (MakeNew("Oxs_UniformScalarField",director,"value 7.85e6")));
  }
  
  // Initialize outputs.
  fprintf(stderr, "YY_StageIsotropicMEL constructor initialize outputs.\n");
  B_MEL_output.Setup(this, InstanceName(), "B max", "mT", 1,
    &YY_StageIsotropicMEL::Fill__B_MEL_output);
  B_MELx_output.Setup(this, InstanceName(), "Bx max", "mT", 1,
    &YY_StageIsotropicMEL::Fill__B_MEL_output);
  B_MELy_output.Setup(this, InstanceName(), "By max", "mT", 1,
    &YY_StageIsotropicMEL::Fill__B_MEL_output);
  B_MELz_output.Setup(this, InstanceName(), "Bz max", "mT", 1,
    &YY_StageIsotropicMEL::Fill__B_MEL_output);

  // Register outputs.
  B_MEL_output.Register(director, 0);
  B_MELx_output.Register(director, 0);
  B_MELy_output.Register(director, 0);
  B_MELz_output.Register(director, 0);

  VerifyAllInitArgsUsed();
}

YY_StageIsotropicMEL::~YY_StageIsotropicMEL()
{}

OC_BOOL YY_StageIsotropicMEL::Init()
{
  fprintf(stderr, "YY_StageIsotropicMEL::Init(): start.\n");
  stage_valid = 0;

  mesh_id = 0;
  displacement.Release();
  stagefield.Release();

  fprintf(stderr, "YY_StageIsotropicMEL::Init(): Calling Oxs_Energy::Init()\n");
  return Oxs_Energy::Init();
}

void YY_StageIsotropicMEL::StageRequestCount(unsigned int& min,
    unsigned int& max) const
{
  if(number_of_stages == 0) {
    min = 0; max = UINT_MAX;
  } else {
    min = max = number_of_stages;
  }
}

void YY_StageIsotropicMEL::ChangeDisplacementInitializer(
    OC_UINT4m stage, const Oxs_Mesh* mesh) const
{
  fprintf(stderr,"YY_StageIsotropicMEL::ChangeDisplacementInitializer(): start.\n");
  // Setup displacement
  vector<String> params;
  if(filelist.empty()) {
    // Use cmd to generate field initializer
    fprintf(stderr,"YY_StageIsotropicMEL::ChangeDisplacementInitializer(): filelist.empty.\n");
    cmd.SaveInterpResult();
    cmd.SetCommandArg(0,stage);
    cmd.Eval();
    cmd.GetResultList(params);
    cmd.RestoreInterpResult();
  } else {
    fprintf(stderr,"YY_StageIsotropicMEL::ChangeDisplacementInitializer(): !filelist.empty.\n");
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

  //fprintf(stderr,"YY_StageIsotropicMEL::ChangeDisplacementInitializer(): OXS_GET_EXT_OBJECT, displacement_init.\n");
  OXS_GET_EXT_OBJECT(params,Oxs_VectorField,displacement_init);
//  CalculateMELField(state);
  working_stage = stage;
  stage_valid = 1;
}

void
YY_StageIsotropicMEL::FillDisplacementCache(const Oxs_SimState& state) const
{
  const Oxs_Mesh* mesh = state.mesh;
  // Displacement
  fprintf(stderr, "YY_StageIsotropicMEL::FillDisplacementCache(): FillMeshValue.\n");
  displacement_init->FillMeshValue(mesh,displacement);
  strain.SetDisplacement(state, displacement);

}

void YY_StageIsotropicMEL::UpdateCache(const Oxs_SimState& state) const
{
  // Update cache as necessary
  if(!stage_valid) {
    // The first go.
    fprintf(stderr, "YY_StageIsotropicMEL:UpdateCache(): !stage_valide.\n");
    mesh_id = 0;
    MELCoef_init->FillMeshValue(state.mesh,MELCoef);
    ChangeDisplacementInitializer(state.stage_number,state.mesh);
    FillDisplacementCache(state);
  } else if(working_stage != state.stage_number) {
    mesh_id = 0;
    fprintf(stderr, "YY_StageIsotropicMEL:UpdateCache(): ChangeDisplacementInitializer().\n");
    ChangeDisplacementInitializer(state.stage_number,state.mesh);
    fprintf(stderr, "YY_StageIsotropicMEL:UpdateCache(): FillDisplacementCache().\n");
    FillDisplacementCache(state);
    mesh_id = state.mesh->Id();
  } else if(mesh_id != state.mesh->Id()) {
    mesh_id = 0;
    fprintf(stderr, "YY_StageIsotropicMEL:UpdateCache(): in else, FillDisplacementCache().\n");
    FillDisplacementCache(state);
    mesh_id = state.mesh->Id();
  }
}


void YY_StageIsotropicMEL::GetEnergy
(const Oxs_SimState& state, Oxs_EnergyData& oed) const
{
  fprintf(stderr, "YY_StageIsotropicMEL::GetEnergy(): Beginning the method.\n");
  OC_INDEX size = state.mesh->Size();
  if(size<1) return;

  fprintf(stderr, "YY_StageIsotropicMEL::GetEnergy(): UpdateCache().\n");
  UpdateCache(state);

  fprintf(stderr, "YY_StageIsotropicMEL::GetEnergy(): CheckMesh().\n");
  if(!displacement.CheckMesh(state.mesh)) {
    throw Oxs_ExtError(this,"Programming error; displacement size"
        " incompatible with mesh.");
  }

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);

  stagefield.AdjustSize(state.mesh);

  fprintf(stderr, "YY_StageIsotropicMEL::GetEnergy(): Starting H_MEL calculation.\n");
  CalculateMELField(state);

  // Use supplied buffer space, and reflect that use in oed.
  oed.energy = oed.energy_buffer;
  oed.field = oed.field_buffer;
  Oxs_MeshValue<OC_REAL8m>& energy = *oed.energy_buffer;
  Oxs_MeshValue<ThreeVector>& field = *oed.field_buffer;
  energy.AdjustSize(state.mesh);
  field.AdjustSize(state.mesh);

  field = stagefield;

  // Calculate energy density
  for(OC_INDEX i=0; i<size; i++) {
    // TODO: replace this equation with a right one.
    energy[i] = -0.5 * MU0 * Ms[i] * (spin[i]*field[i]);
  }

}

void
YY_StageIsotropicMEL::Fill__B_MEL_output(const Oxs_SimState& state)
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

void YY_StageIsotropicMEL::CalculateMELField(
  const Oxs_SimState& state) const
{
  fprintf(stderr,"YY_StageIsotropicMEL::CalculateMELField()\n");
  // TODO: Check mesh validity
  const OC_INDEX size = state.mesh->Size();
  if(size<1) return;

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  const Oxs_MeshValue<OC_REAL8m>& Msi = *(state.Ms_inverse);
  const Oxs_RectangularMesh* mesh =
    dynamic_cast<const Oxs_RectangularMesh*>(state.mesh);
  const OC_INDEX xdim = mesh->DimX();
  const OC_INDEX ydim = mesh->DimX();
  const OC_INDEX zdim = mesh->DimX();
  const OC_INDEX xydim = xdim*ydim;

  // Compute MEL field
  fprintf(stderr,"Compute MEL field\n");
  for(OC_INDEX i=0; i<size; i++) {
    // stagefield[i]*strain.diag[i] returns a dot product. Don't use it.
    stagefield[i].x  = spin[i].x*strain.diag[i].x;
    stagefield[i].y  = spin[i].y*strain.diag[i].y;
    stagefield[i].z  = spin[i].z*strain.diag[i].z;
    stagefield[i].x += spin[i].y*strain.offdiag[i].z+spin[i].z*strain.offdiag[i].y;
    stagefield[i].y += spin[i].x*strain.offdiag[i].z+spin[i].z*strain.offdiag[i].x;
    stagefield[i].z += spin[i].x*strain.offdiag[i].y+spin[i].y*strain.offdiag[i].x;
    stagefield[i] *= -1/MU0*2*Msi[i]*MELCoef[i];
  }
  if(hmult != 1.0) stagefield *= hmult;

  // H-field
  fprintf(stderr, "Setting max field.\n");
  max_field.Set(0.,0.,0.);
  if(size>0) {
    fprintf(stderr, "size = %ld\n", size);
    OC_INDEX max_i = 0;
    OC_REAL8m max_magsq = stagefield[OC_INDEX(0)].MagSq();
    for(OC_INDEX i=1; i<size; i++) {
      OC_REAL8m magsq = stagefield[i].MagSq();
      if(magsq>max_magsq) {
        max_magsq = magsq;
        max_i = i;
      }
    }
    fprintf(stderr,"max_i = %ld\n",max_i);
    fprintf(stderr, "max_magsq = %e\n", max_magsq);
    max_field = stagefield[max_i];
  }

  //DisplayValues(state,6,8,14,16,1,3);

  // UpdateCache(state); // Do we need this?

}

void YY_StageIsotropicMEL::DisplayValues(
    const Oxs_SimState& state,
    OC_INDEX xmin, OC_INDEX xmax,
    OC_INDEX ymin, OC_INDEX ymax,
    OC_INDEX zmin, OC_INDEX zmax) const
{
  const Oxs_RectangularMesh* mesh =
    static_cast<const Oxs_RectangularMesh*>(state.mesh);
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);

  // MELCoef
  fprintf(stderr,"MELCoef:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",MELCoef[i]);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // strain.diag.x
  fprintf(stderr,"strain.diag.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",strain.diag[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // strain.diag.y
  fprintf(stderr,"strain.diag.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",strain.diag[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // strain.diag.z
  fprintf(stderr,"strain.diag.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",strain.diag[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // strain.offdiag.x
  fprintf(stderr,"strain.offdiag.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",strain.offdiag[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // strain.offdiag.y
  fprintf(stderr,"strain.offdiag.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",strain.offdiag[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // strain.offdiag.z
  fprintf(stderr,"strain.offdiag.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",strain.offdiag[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // stagefield.x
  fprintf(stderr,"stagefield.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",stagefield[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // stagefield.y
  fprintf(stderr,"stagefield.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",stagefield[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // stagefield.z
  fprintf(stderr,"stagefield.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",stagefield[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

}
