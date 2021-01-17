#!/bin/sh

prefix=AEM-4501-S21-matlab-code

fileList="README.TXT\
         DISTMESH_COPYRIGHT.TXT\
         DISTMESH_FILE_LIST.TXT\
         LICENSE.TXT\
         PD_torsion.m\
         PD_torsion_poly.m\
         PD_truss_modes.m\
         PD_truss_static.m\
         PlotTruss.m\
         PlotTrussMode.m\
         README.TXT\
         assemble_PD_torsion_fem.m\
         assemble_PD_truss_fem.m\
         assemble_bar_fem.m\
         assemble_gmesh_bar_fem.m\
         assemble_gmesh_beam_fem.m\
         bar_dynamic.m\
         bar_forced.m\
         bar_impact.m\
         bar_modes.m\
         bar_static.m\
         beam_free_vib.m\
         beam_interpolate_results.m\
         boundedges.m\
         butterfly_torsion_example.m\
         create_gmesh_for_uniform_mesh.m\
         distmesh2d.m\
         dpoly.m\
         dsegment.m\
         dsegment.mexa64\
         dsegment.mexmaci64\
         dsegment.mexw32\
         dsegment.mexw64\
         ellipse_torsion_example.m\
         fixmesh.m\
         gmesh_bar_dynamic.m\
         gmesh_bar_modes.m\
         gmesh_bar_static.m\
         gmesh_beam_dynamic.m\
         gmesh_beam_modes.m\
         gmesh_beam_static.m\
         huniform.m\
         rectangle_torsion_example.m\
         simpplot.m\
         simpvol.m\
         spring_mass_forced.m\
         triangle_torsion_example.m\
         truss_2d_example.m\
         truss_example.m\
         truss_local_matrices.m"
fileList=`printf "${fileList}" | sed -e 's/[[:space:]][[:space:]]*/ /g'`


mkdir "${prefix}"
cp ${fileList} "${prefix}"/

zip -r "${prefix}.zip" "${prefix}"

rm -rf "${prefix}"
