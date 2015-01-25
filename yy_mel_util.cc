/* FILE: yy_MEL_util.cc                 -*-Mode: c++-*-
 *
 * OOMMF magnetoelastic coupling extension module.
 * yy_mel_util.* contain YY_MELField class , which is to be used by other
 * YY_*MEL classes
 * 
 */

#include "yy_MEL_util.h"

#include "nb.h"
#include "energy.h"
#include "mesh.h"
#include "rectangularmesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "vectorfield.h"
#include "util.h"

/* End includes */

YY_MELField::YY_MELField():
  displacement_valid(0), strain_valid(0), 
  MELCoef1_valid(0), MELCoef2_valid(0)
{
}

void YY_MELField::Release()
{
  u.Release();
  e_diag.Release();
  e_offdiag.Release();
  MELCoef1.Release();
  MELCoef2.Release();
}

void YY_MELField::SetMELCoef(const Oxs_SimState& state,
    const Oxs_OwnedPointer<Oxs_ScalarField>& MELCoef1_init,
    const Oxs_OwnedPointer<Oxs_ScalarField>& MELCoef2_init)
{
  MELCoef1_init->FillMeshValue(state.mesh,MELCoef1);
  MELCoef2_init->FillMeshValue(state.mesh,MELCoef2);
  MELCoef1_valid = 1;
  MELCoef2_valid = 1;
}

void YY_MELField::SetDisplacement(const Oxs_SimState& state,
    const Oxs_OwnedPointer<Oxs_VectorField>& u_init)
{
  // Check mesh size
  const OC_INDEX size = state.mesh->Size();
  if(size<1) return;

  u_init->FillMeshValue(state.mesh,u);
  displacement_valid = 1;

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

  e_diag.AdjustSize(mesh);
  e_offdiag.AdjustSize(mesh);

  // Compute du/dx
  for(OC_INDEX z=0; z<zdim; z++) {
    for(OC_INDEX y=0; y<ydim; y++) {
      for(OC_INDEX x=0; x<xdim; x++) {
        OC_INDEX i = mesh->Index(x,y,z);  // Get base index
        ThreeVector du_dx;
        ThreeVector du_dy;
        ThreeVector du_dz;
        if(Ms[i]==0.0) {
          e_diag[i].Set(0.0, 0.0, 0.0);
          e_offdiag[i].Set(0.0, 0.0, 0.0);
          continue;
        }

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

        e_diag[i].Set(du_dx.x,du_dy.y,du_dz.z);
        e_offdiag[i].Set(
          0.5*(du_dz.y+du_dy.z),
          0.5*(du_dx.z+du_dz.x),
          0.5*(du_dy.x+du_dx.y)
        );
      }
    }
  }

  strain_valid = 1;

  // For debug, display values
  // arguments: state, xmin, xmax, ymin, ymax, zmin, zmax
#ifdef YY_DEBUG
  DisplayValues(state,4,6,0,2,0,2);
#endif

}

void YY_MELField::SetStrain(const Oxs_SimState& state,
    const Oxs_OwnedPointer<Oxs_VectorField>& e_diag_init,
    const Oxs_OwnedPointer<Oxs_VectorField>& e_offdiag_init)
{
  if(state.mesh->Size()<1) return;

  e_diag_init->FillMeshValue(state.mesh,e_diag);
  e_offdiag_init->FillMeshValue(state.mesh,e_offdiag);
  strain_valid = 1;

#ifdef YY_DEBUG
  DisplayValues(state,4,6,0,2,0,2);
#endif
}

void YY_MELField::CalculateMELField(
  const Oxs_SimState& state,
  OC_REAL8m hmult,
  Oxs_MeshValue<ThreeVector>& field_buf,
  Oxs_MeshValue<OC_REAL8m>& energy_buf) const
{
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

  // Compute MEL field and energy
  Oxs_MeshValue<ThreeVector> temp_field;
  temp_field.AdjustSize(mesh);
  for(OC_INDEX i=0; i<size; i++) {
    // field_buf[i]*e_diag[i] returns a dot product. Don't use it.
    field_buf[i].x  = spin[i].x*e_diag[i].x;
    field_buf[i].y  = spin[i].y*e_diag[i].y;
    field_buf[i].z  = spin[i].z*e_diag[i].z;
    field_buf[i] *= -1/MU0*2*Msi[i]*MELCoef1[i];
    temp_field[i].x = spin[i].y*e_offdiag[i].z+spin[i].z*e_offdiag[i].y;
    temp_field[i].y = spin[i].x*e_offdiag[i].z+spin[i].z*e_offdiag[i].x;
    temp_field[i].z = spin[i].x*e_offdiag[i].y+spin[i].y*e_offdiag[i].x;
    temp_field[i] *= -1/MU0*2*Msi[i]*MELCoef2[i];
    field_buf[i] += temp_field[i];
  }
  if(hmult != 1.0) field_buf *= hmult;

  // H-field
  max_field.Set(0.,0.,0.);
  if(size>0) {
    OC_INDEX max_i = 0;
    OC_REAL8m max_magsq = field_buf[OC_INDEX(0)].MagSq();
    for(OC_INDEX i=1; i<size; i++) {
      OC_REAL8m magsq = field_buf[i].MagSq();
      if(magsq>max_magsq) {
        max_magsq = magsq;
        max_i = i;
      }
    }
    max_field = field_buf[max_i];
  }

  // Energy density
  for(OC_INDEX i=0; i<size; i++) {
    energy_buf[i] = (-0.5*MU0*Ms[i])*(field_buf[i]*spin[i]);
  }
}

#ifdef YY_DEBUG
void YY_MELField::DisplayValues(
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

  if(displacement_valid) {
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
  }

  // e_diag.x
  fprintf(stderr,"e_diag.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",e_diag[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // e_diag.y
  fprintf(stderr,"e_diag.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",e_diag[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // e_diag.z
  fprintf(stderr,"e_diag.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",e_diag[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // e_offdiag.x
  fprintf(stderr,"e_offdiag.x:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",e_offdiag[i].x);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // e_offdiag.y
  fprintf(stderr,"e_offdiag.y:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",e_offdiag[i].y);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

  // e_offdiag.z
  fprintf(stderr,"e_offdiag.z:\n");
  for(OC_INDEX y=ymin; y<ymax+1; y++) {
    for(OC_INDEX z=zmin; z<zmax+1; z++) {
      for(OC_INDEX x=xmin; x<xmax+1; x++) {
        OC_INDEX i = mesh->Index(x,y,z);
        fprintf(stderr,"%e ",e_offdiag[i].z);
      }
      fprintf(stderr,"| ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");

}
#endif  // YY_DEBUG
