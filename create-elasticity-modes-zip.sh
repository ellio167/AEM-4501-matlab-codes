#!/bin/sh

prefix=elasticity-modes-matlab

fileList="DISTMESH_COPYRIGHT.TXT\
         DISTMESH_FILE_LIST.TXT\
         boundedges.m\
         distmesh2d.m\
         distmeshnd.m\
         dpoly.m\
         dsegment.m\
         dsegment.mexa64\
         dsegment.mexmaci64\
         dsegment.mexw32\
         dsegment.mexw64\
         fixmesh.m\
         huniform.m\
         simpplot.m\
         simpvol.m\
         surftri.m\
         PD_elasticity_modes.m\
         assemble_PD_elasticity_fem.m\
         PD_3d_elasticity_modes.m\
         assemble_PD_3d_elasticity_fem.m\
         PlotMode_elasticity.m\
         PlotMode_elasticity_movie.m\
         PlotMode_elasticity_trajectory_movie.m\
         butterfly_elasticity_example.m\
         churro_elasticity_example.m\
         rectangle_elasticity_example.m\
         sphere_elasticity_example.m"
fileList=`printf "${fileList}" | sed -e 's/[[:space:]][[:space:]]*/ /g'`


mkdir "${prefix}"
cp ${fileList} "${prefix}"/

zip -r "${prefix}.zip" "${prefix}"

rm -rf "${prefix}"
