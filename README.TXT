AEM 4501 matlab files

place this directory in a convenient folder on your computer.

Use the 'pathtool' command to add this directory to Matlab's search path, so
that it can find these files.



== 1D Static Elastic Bar =======================================================

bar_static.m

    This script uses a uniform FEM mesh to solve the problem of a linear
    elastic bar subjected to an end load, Fz, and a distributed
    force-per-unit-length, qz.  The input is as simple as possible.  The output
    includes displacements, reactions, internal forces, and strains.

    Dependencies: assemble_bar_fem.m

assemble_bar_fem.m

    This script assembles the Mass and Stiffness matrices and the Force vector
    for the linear elastic bar subjected to a distributed load, qz, for a
    uniform mesh.  It uses 5 point Guass-quadrature to integrate the force
    vector.

gmesh_bar_static.m

    This script uses a general mesh (gmesh) to solve the problem of a linear
    elastic bar (tension or torsion) fixed at the Left and subjected to an end
    load, Fz, and a distributed load, qz.  The input uses a set of data
    structures to describe a general mesh.  The output includes displacements,
    reaction forces, internal forces, and strains.

    Dependencies: assemble_gmesh_bar_fem.m

assemble_gmesh_bar_fem.m

    This script assembles the Mass and Stiffness matrices and the Force vector
    for the linear elastic bar subjected to a distributed load, qz, for a
    general mesh.  It uses 5 point Guass-quadrature to integrate the force
    vector.

create_gmesh_for_uniform_mesh.m

    This script takes input similar to "bar_static" and provides output that is
    appropriate for use as input to "gmesh_bar_static".

The only real difference between bar_static and gmesh_bar_static is that the
gmesh versions implement the more general data structures.





== 1D Dynamic Elastic Bar ======================================================

bar_dynamic.m

    This script uses a uniform FEM mesh to integrate the equations of motion
    for a linear elastic bar subjected to an end load, Fz, and a distributed
    force-per-unit-length, qz.  This input is as simple as possible.  The
    output is an augmented 'sol' structure from matlab's ode45 solver.

    Dependencies: assemble_bar_fem.m


gmesh_bar_dynamic.m

    This script uses a general mesh (gmesh) to solve the problem of a linear
    elastic bar (tension or torsion) fixed at the Left and subjected to an end
    load, Fz, and a distributed load, qz.  The input uses a set of data
    structures to describe a general mesh.  The output includes displacements,
    reaction forces, internal forces, and strains.

    Dependencies: assemble_gmesh_bar_fem.m


bar_modes.m

    This script uses a uniform mesh to compute the natural frequencies and mode
    shapes for the elastic bar problems.  The output includes the frequencies
    and mode shapes.

    Dependencies: assemble_bar_fem.m


gmesh_bar_modes.m

    This script computes the natural frequencies and mode shapes for the
    elastic bar problems using a general mesh.  The output includes the
    frequencies and mode shapes.

    Dependencies: assemble_gmesh_bar_fem.m



spring_mass_forced.m

    An example script that uses a 2 node mesh to create a problem equivalent to
    a simple spring-mass oscillator subjected to harmonic forcing.  It uses
    "bar_modes.m" to determine a good time step.  It then creates two animated
    gif files plotting the results in various ways.

    Dependencies: bar_dynamic.m; bar_modes.m

bar_forced.m

    An example script that is the same as spring_mass_forced.m except it uses a
    20 node mesh for a bar.

    Dependencies: bar_dynamic.m; bar_modes.m


bar_impact.m

    An example script that is the same as bar_forced.m except it is subjected
    to an "impact" loading corresponding to an initial non-zero Right end
    velocity.

    Dependencies: bar_dynamic.m; bar_modes.m





== 2D Static Beam ==============================================================

gmesh_beam_static.m

     This script solves the static planar beam bending problem with build-in
     Left end and subjected to and shear force, Vy, and end bending moment, Mx,
     as well as a distributed shear loading, qy, using a general mesh.  The
     input allows for variable element size and material and geometric
     properties.  The output includes the solution in terms of nodal
     displacements and rotations, element moments and curvatures.  In addition,
     the output includes a list of displacements interior to each elements in
     order to make accurate plotting of the deformed configuration more
     straight forward.

     Dependencies: assemble_gmesh_beam_fem.m; beam_interpolate_results.m


assemble_gmesh_beam_fem.m

     This script assembles the Mass and Stiffness matrices and the Force vector
     for the planar beam bending problem subjected to a distributed shear
     loading, qy, using a general mesh.  It uses 5 point Gauss-quadrature to
     integrate the force vector.


beam_interpolate_results.m

     This script uses the beam shape functions to interpolate the solution
     vector (obtained from gmesh_beam_static.m) and generate the element
     moments and curvatures as well as the list of displacements interior to
     each element.





== 2D Dynamic Beam =============================================================


gmesh_beam_dynamic.m

     This script uses a general mesh (gmesh) to solve the dynamic problem of a
     planar beam bending problem with built-in Left end and subjected to and
     shear force, Vy, and end bending moment, Mx, as well as a distributed
     shear loading, qy, using a general mesh.  Initial conditions are given in
     terms of initial nodal displacements and velocities.  The output is an
     augmented 'sol' structure from matlab's ode45 solver.

     Dependencies: assemble_gmesh_beam_fem.m; beam_interpolate_results.m


gmesh_beam_modes.m

     This script computes the natural frequencies and mode shapes for the
     planar beam bending problem with built-in Left end.  The output includes
     a list of the requested number of natural frequencies, the corresponding
     mode shapes, mass and stiffness matrices.

     Dependencies: assemble_gmesh_beam_fem.m; beam_interpolate_results.m


beam_free_vib.m

     An example script that computes the free vibration response of a beam
     subjected to an initially deformed configuration corresponding to the
     equilibrium configuration associated with end loading of Vy=-5000,
     Mx=-150000.  This part of the script uses gmesh_beam_dynamic.m and creates
     the file beam_free_vib_saved.mat containing the results.  The next part of
     the script decomposes the initial conditions into natural mode components
     using the gmesh_beam_modes.m output.  The script then computes the
     time-dependent natural mode amplitudes from the explicit solution of the
     free vibration problem in terms of the beam's natural modes.  The script
     then plots the displacement of the Right end as a function of time for
     both solutions.  It also plots the full deformed configuration at 5
     different times using the full dynamic solution and a modal solution
     containing only the first 5 modes.  We see that these compare very well.
     You can rerun the last section of this script (cut-and-paste) using
     'order=x' with different values of x in order to see how the solution
     converges with the number of modes used.  (This script can take some time
     to run.)





== 3D Static Trusses ===========================================================


PD_truss_static.m

     This script solves the static 3D truss deformation problem.  A truss is
     specified in terms of nodes, elements, material sets, boundary condition
     types, and values.  Each nodal DOF is assigned a BCType (disp. or force)
     and a value (the numerical value of the disp. or force).  Input is
     encapsulated inside a matlab structure variable containing all the Problem
     Data (PD).  The script returns a copy of the original PD input structure
     augmented with output data including the nodal displacements, nodal
     reaction forces, element internal force values, and element stresses.

     Dependencies: assemble_PD_truss_fem.m; truss_local_matrices.m


assemble_PD_truss_fem.m

     This script assembles the Mass and Stiffness matrices for the linear 3D
     truss problem.


truss_local_matrices.m

     This script computes the element mass, stiffness, and coordinate
     transformation matrices for a 3D truss element.  The input consists of the
     two nodal coordinates for the element and a matlab structure containing
     the element's material and geometric properties.


PlotTruss.m

     This script generates an interactive 3D plot of the truss structure
     problem.  The input allows for specification of the scaling of the
     deformations and the plotting of labels.  If the deformed structure is not
     plotted colors indicate the material set of each element.  Dash marks at
     a node indicate applied boundary conditions (color coded for disp. or
     force type).  The plot may be rotated in 3D using the standard matlab plot
     view tools.





== 3D Dynamic Trusses ==========================================================


PD_truss_modes.m

     This script computes the free vibration natural modes for the 3D truss
     problem.  The output includes the number of requested mode frequencies and
     their shapes.

     Dependencies: assemble_PD_truss_fem.m


PlotTrussMode.m

     This script takes as input the output of PD_truss_modes.m and plots the
     deformed configuration of the requested natural mode.


truss_2d_example.m

     An example 2D truss with 2 material sets, 6 nodes, and 9 elements.  This
     is a good example to demonstrate the concept of "Mechanisms" by removing
     one or more of elements 3, 6, and 7.


truss_example.m

     An example 3D truss with 2 material sets, 7 nodes, and 9 elements.





== Prandtl Torsion =============================================================


PD_torsion_poly.m

     This script solves the Prandtl torsion problem for a polygonal cross
     section.  The input includes a list of vertexes, a bounding box, and a
     numerical value for the right hand side of the torsion equations.  Output
     consists of a copy of the input PD data structure augmented with lots of
     additional information as listed in the function documentation.  This
     script uses the DISTMESH package routines to construct the triangular FEM
     mesh.

     Dependencies: PD_torsion.m


PD_torsion.m

     This script solves the Prandtl torsion problem for a general cross section
     specified by a function that is zero on the boundary of the section,
     negative inside the section, and positive outside the section.  The input
     includes a function that describes the boundary, a bounding box, and a
     numerical value for the right hand side of the torsion equations.  Output
     consists of a copy of the input PD data structure augmented with lots of
     additional information as listed in the function documentation.  This
     script uses the DISTMESH package routines to construct the triangular FEM
     mesh.

     Dependencies: assemble_PD_torsion_fem.m; DISTMESH files


assemble_PD_torsion_fem.m

     This script assembles the Stiffness matrix and Force vector for the
     Prandtl torsion problem using the mesh generated by the DISTMESH utilities
     and the simple formula for the "constant-strain-triangle" element
     stiffness matrix and force vector.

ellipse_torsion_example.m
triangle_torsion_example.m
rectangle_torsion_example.m
butterfly_torsion_example.m

      Examples of setup and use of both PD_torsion.m and PD_torsion_poly.m


DIST MESH files:


     DISTMESH_COPYRIGHT.TXT
     DISTMESH_FILE_LIST.TXT
     boundedges.m
     distmesh2d.m
     dpoly.m
     dsegment.m
     dsegment.mexa64
     dsegment.mexmaci64
     dsegment.mexw32
     dsegment.mexw64
     fixmesh.m
     huniform.m
     simpplot.m
     simpvol.m
