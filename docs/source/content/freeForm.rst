Free-Form Deformation
=====================

Required inputs
^^^^^^^^^^^^^^^

To run the code an **input file** must be provided with the location of the VTK file to morph and the definition of all the Bernstein Boxes used for the transformation. 

Running the code
^^^^^^^^^^^^^^^^

To run the code once feMorph is compile type ::
 
  feMorph -j -f commandfile.txt

where *commandfile.txt* is the required input file.

Input File Format
^^^^^^^^^^^^^^^^^

An example of the valid input file is ::

  mesh.vtk
  1
  0.0 0.0 0.0
  3 3 3
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 5.0
  2
  5 0.0 1.0 0.0
  6 1.0 0.0 0.0
  ...

The first line of the file contains the **location of the VTK file to be morphed**. This file is an ASCII vtk legacy file. ::

  [line 1] mesh.vtk

The second line defines the total number of Bernstein Boxes to use. In this case a single box is selected ::

  [line 2] 1

The origin node is specified next ::

  [line 3] 0.0 0.0 0.0

The number of **control points** along three orthogonal axis is specified next. Control points are used to assign free-form deformations ::

  [line 4] 3 3 3 

The specification of the three axis follows. Note that the lenght of the axis is important to define the size of the Bernstein Box ::

  [line 5] 5.0 0.0 0.0
  [line 6] 0.0 5.0 0.0
  [line 7] 0.0 0.0 5.0

The next line specifies how many nodes in the box have displacements applied. ::

  [line 8] 2

In this case two nodes have displacements assigned and these are ::

  [line 9]  5 0.0 1.0 0.0
  [line 10] 6 1.0 0.0 0.0

Note that the displacements assigned to nodes 5 and 6 are referred to the axis system of the current Bernstein Box. 

In case more than one boxes are needed, the block from line 3 to 10 must be repeated.






