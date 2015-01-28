YY_MEL: OOMMF Magnetoelastic Extension Module
Release ver. 1.0.0
This file was last modified on 2015-01-28.


What's this:

YY_MEL module is an OOMMF extension to include the magnetoelastic copuling in
the micromagnetic simulation. The code is based on the work developed for
Yahagi et al, Phys. Rev. B, 90, 140405(R) (2014).


Basic principle:

The magnetoelastic (MEL) energy can be expressed in terms of the 
magnetoelastic coefficients B1 and B2, as a function of the magnetization and
elastic strain. By taking the variational derivative of the MEL energy
density, the MEL contribution to the effective field is obtained. This module 
takes either displacement or elastic strain, calculates MEL field and include 
it as an additional energy term in the micromagnetic simulation.


Installation and usage:

Copy .cc and .h files into the directory /app/oxs/local/ and run pimake as
>> tclsh oommf.tcl pimake

Include the following Oxs_Ext classes in your mif file and run it.


Oxs_Ext classes:

This extension contains three Oxs_Ext classes to include magnetoelastic (MEL)
energy terms in OOMMF simulations. The classes YY_FixedMEL and YY_StageMEL are
designed in the semblance of Oxs_FixedZeeman and Oxs_StageZeeman to account
for spatially nonuniform, non-time-varying strain. YY_TransformStageMEL
produces time-varying MEL field by applying a spatially uniform linear
transform, much in the way same as Oxs_TransformZeeman does. Unlike
Oxs_TransformZeeman which reads the input vector field only once at the
beginning of the simulation, YY_TransformStageMEL updates the untransformed 
vector field (either displacement or strain) at each stage and provides more
flexibility.

Common features:

In all three YY_*MEL classes, some initialization strings are always required.
B1 and B2 take a Oxs_ScalarField object to specify the MEL coefficients in the
unit of J/m^3. In most cases, Oxs_AtlasScalarField will suffice. In 
elastically isotropic media, B1 = B2 but you need to specify both of them 
anyway.

The elastic strain needs to be specified, either directly with a strain matrix
or indirectly with a displacement vector field. If the displacement is
specified, the strain is automatically calculated and is used for the
subsequent computation. Only one of them can be specified in each Specify
block.

The displacement u and strain e are handled as Oxs_VectorField objects. The
strain, however, is a symmetric real matrix and requires six elements to be
stored. In YY_*MEL modules, the diagonal and off-diagonal elements of the
strain are expressed with two Oxs_VectorField objects, together representing
a strain matrix. The indices of these vectors correspond to the matrix index
as depicted below.

 Strain matrix   Diag vec   Off-diag vec
  / 11 12 13 \   / 0     \   /   2 1 \
  | .. 22 23 | = |   1   | + |     0 |
  \ .. .. 33 /   \     2 /   \       /

YY_FixedMEL:

  Specify YY_FixedMEL {
    B1              scalar_field_spec
    B2              scalar_field_spec
    u_field         vector_field_spec
    e_diag_field    vector_field_spec
    e_offdiag_field vector_field_spec
    multiplier      number
  }

YY_FixedMEL takes a fixed displacement or strain and calculates the MEL field
and energy. Note that even though the strain remains constant throughout the
simulation, MEL field and energy are repeatedly updated reflecting the
magnetization evolution. The calculated MEL field and energy are multiplied by
the value of the multiplier parameter.

Either one of the displacement (u_field) or the strain (BOTH of e_diag_field
AND e_offdiag_field) should be specified.

  Example 1.
  Specify Oxs_UniformScalarField:B {
    value 7.85e6
  }
  Specify YY_FixedMEL:MELField {
    B1 :B
    B2 :B
    u_field :fileufield
  }

YY_StageMEL:

  Specify YY_FixedMEL {
    B1               scalar_field_spec
    B2               scalar_field_spec
    u_script         Tcl_script
    u_files          { list_of_files }
    e_diag_script    Tcl_script
    e_diag_files     { list_of_files }
    e_offdiag_script Tcl_script
    e_offdiag_files  { list_of_files }
    multiplier       number
    stage_count      number_of_stages
  }

For B1, B2, multiplier, see YY_FixedMEL. The stage_count functions as in
Oxs_StageZeeman and notifies the Oxs_Driver how many stages the YY_StageMEL
object wants. The default is 0, indicating no preference.

In YY_StageMEL, the displacement is specified as a Tcl script, which is
evaluated at each field to return an Oxs_VectorField object. See the section
for Oxs_StageZeeman of the OOMMF documentation for an example. Alternatively,
a list of ovf files may be specified as well.

The strain can also be specified in place of the displacement. Use
e_*diag_script and e_*diag_files parameters. Both the diagonal and
off-diagonal elements need to be specified.

YY_TransformStageMEL:

  Specify YY_TransformStageMEL {
    B1               scalar_field_spec
    B2               scalar_field_spec
    u_script         Tcl_script
    u_files          { list_of_files }
    e_diag_script    Tcl_script
    e_diag_files     { list_of_files }
    e_offdiag_script Tcl_script
    e_offdiag_files  { list_of_files }
    multiplier       number
    stage_count      number_of_stages

    type             transform_type
    script_args      { args_request }
    script           Tcl_script
  }

YY_TransformStageMEL is an extension of YY_StageMEL and has all the parameters
the latter has. See YY_StageMEL for descriptions for these parameters. It
applies a linear transformation to the displacement or the strain, whichever
is specified by the user, and the MEL field and energy is calculated at each
simulation step. 

The rest of the parameters, type, script_args, and script behave as they are
defined in Oxs_TransformZeeman. See the OOMMF documentation for the detailed
description. Note, however, that the transform matrix is applied differently
to the displacement and the strain. When M denotes the transformation matrix,
the displacement is transformed as the applied field is transformed in
Oxs_TransformZeeman,
  u_new = M*u_old
while the strain is transformed as
  e_new = M*e_old*transpose(M)
The type of the elastic input determines the type of linaer transformation.


Contact:
<yuyahagi2@gmail.com>
