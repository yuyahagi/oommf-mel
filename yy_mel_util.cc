/* FILE: yy_MEL_util.cc                 -*-Mode: c++-*-
 *
 * Misc
 * 
 */

#include "nb.h"

#include "yy_MEL_util.h"

OC_USE_STRING;

/* End includes */

YY_Strain::YY_Strain():
  displacement_valid(0), strain_valid(0),
  diag(NULL), offdiag(NULL)
{
  // Set up initializers for the strain vectors
  // Use cmd to generate field initializer
  //vector<String> params;
  //cmd.SaveInterpResult();
  //cmd.SetCommandArg(0,stage);
  //cmd.Eval();
  //cmd.GetResultList(params);
  //cmd.RestoreInterpResult();

//  vector<String> params;
//  params.push_back("Oxs_UniformVectorField");
//  params.push_back("vector {1 0 0}");
//  OXS_GET_EXT_OBJECT(params,Oxs_VectorField,diag_init);

  //OXS_GET_EXT_OBJECT(params,Oxs_VectorField,diag_init);

}

void YY_Strain::Init(const Oxs_SimState& state)
{
  diag.Release();
  offdiag.Release();
}

void YY_Strain::CalculateStrain(const Oxs_SimState& state)
{
  // TODO: Check mesh validity
  fprintf(stderr, "YY_Strain::CalculateStrain(): Check mesh validity.\n");
  const OC_INDEX size = state.mesh->Size();
  if(size<1) return;
  // UpdateCache(state); // Do we need this?

  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_RectangularMesh* mesh =
    dynamic_cast<const Oxs_RectangularMesh*>(state.mesh);
  const OC_INDEX xdim = mesh->DimX();
  const OC_INDEX ydim = mesh->DimY();
  const OC_INDEX zdim = mesh->DimZ();
  const OC_INDEX xydim = xdim*ydim;
  const OC_REAL8m idelx = 1./mesh->EdgeLengthX();
  const OC_REAL8m idely = 1./mesh->EdgeLengthX();
  const OC_REAL8m idelz = 1./mesh->EdgeLengthX();

  diag.AdjustSize(mesh);
  offdiag.AdjustSize(mesh);

  fprintf(stderr,"ThreeVector::Set() test.\n");
  {
    ThreeVector temp;
    temp.Set(1.,2.,3.);
    fprintf(stderr,"temp = (%e,%e,%e)\n",temp.x,temp.y,temp.z);
    temp = u[7];
    fprintf(stderr,"temp = (%e,%e,%e)\n",temp.x,temp.y,temp.z);
  }

  // Compute du/dx
  fprintf(stderr, "YY_Strain::CalculateStrain(): Compute du/dx.\n");
  for(OC_INDEX z=0; z<zdim; z++) {
    for(OC_INDEX y=0; y<ydim; y++) {
      for(OC_INDEX x=0; x<xdim; x++) {
        OC_INDEX i = mesh->Index(x,y,z);  // Get base index
        ThreeVector du_dx;
        ThreeVector du_dy;
        ThreeVector du_dz;
        if(Ms[i]==0.0) {
  //fprintf(stderr, "Ms==0, diag[i].Set\n");
          diag[i].Set(0.0, 0.0, 0.0);
          offdiag[i].Set(0.0, 0.0, 0.0);
          continue;
        }

  //fprintf(stderr, "du/dx\n");
        if(x<xdim-1 && Ms[i+1]!=0.0) du_dx  = u[i+1];
        else                         du_dx  = u[i];
        if(x>0 && Ms[i-1]!=0.0)      du_dx -= u[i-1];
        else                         du_dx -= u[i];
        if(x<xdim-1 && Ms[i+1]!=0.0 && x>0 && Ms[i-1]!=0.0)
          du_dx *= 0.5*idelx;
        else
          du_dx *= idelx;

        if(y<ydim-1 && Ms[i+xdim]!=0.0) du_dy  = u[i+xdim];
        else                            du_dy  = u[i];
        if(y>0 && Ms[i-xdim]!=0.0)      du_dy -= u[i-xdim];
        else                            du_dy -= u[i];
        if(y<ydim-1 && Ms[i+xdim]!=0.0 && y>0 && Ms[i-xdim]!=0.0)
          du_dy *= 0.5*idely;
        else
          du_dy *= idely;

        if(z<zdim-1 && Ms[i+xydim]!=0.0) du_dz  = u[i+xydim];
        else                             du_dz  = u[i];
        if(z>0 && Ms[i-xydim]!=0.0)      du_dz -= u[i-xydim];
        else                             du_dz -= u[i];
        if(z<zdim-1 && Ms[i+xydim]!=0.0 && z>0 && Ms[i-xydim]!=0.0)
          du_dz *= 0.5*idelz;
        else
          du_dz *= idelz;

  //fprintf(stderr, "diag[i].Set() ");
        diag[i].Set(du_dx.x,du_dy.y,du_dz.z);
  //fprintf(stderr, "offdiag[i].Set()\n");
        offdiag[i].Set(
          0.5*(du_dz.y+du_dy.z),
          0.5*(du_dx.z+du_dz.x),
          0.5*(du_dy.x+du_dx.y)
        );
      }
    }
    strain_valid = 1;
  }

  // For debug, display values
  // arguments: state, xmin, xmax, ymin, ymax, zmin, zmax
  DisplayValues(state,4,6,0,2,0,2);

}

void YY_Strain::SetDisplacement(const Oxs_SimState& state,
    const Oxs_MeshValue<ThreeVector>& u_in)
{
  u = u_in;
  displacement_valid = 1;
  strain_valid = 0;
  fprintf(stderr,"u_in[i], u[i]:\n");
  for(OC_INDEX i=0; i<15; i++) {
    fprintf(stderr,"%e ",u_in[i].x);
    fprintf(stderr,"%e ",u[i].x);
  }

  fprintf(stderr, "YY_Strain::SetDisplacement(): CalculateStrain().\n");
  CalculateStrain(state);
}

void YY_Strain::DisplayValues(
    const Oxs_SimState& state,
    OC_INDEX xmin, OC_INDEX xmax,
    OC_INDEX ymin, OC_INDEX ymax,
    OC_INDEX zmin, OC_INDEX zmax) const
{
  const Oxs_RectangularMesh* mesh =
    static_cast<const Oxs_RectangularMesh*>(state.mesh);
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);

  // Indices
  fprintf(stderr,"Indices:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%ld %ld %ld %ld ",x,y,z,i);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // Ms
  fprintf(stderr,"Ms:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",Ms[i]);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // u.x
  fprintf(stderr,"u.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",u[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // u.y
  fprintf(stderr,"u.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",u[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // u.z
  fprintf(stderr,"u.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",u[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // diag.x
  fprintf(stderr,"diag.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",diag[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // diag.y
  fprintf(stderr,"diag.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",diag[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // diag.z
  fprintf(stderr,"diag.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",diag[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // offdiag.x
  fprintf(stderr,"offdiag.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",offdiag[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // offdiag.y
  fprintf(stderr,"offdiag.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",offdiag[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // offdiag.z
  fprintf(stderr,"offdiag.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",offdiag[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

}
