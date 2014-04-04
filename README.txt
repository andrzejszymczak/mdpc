
/* -------------------------------------------------------------------------- */

REQUIREMENTS

This implementation has been tested only on linux systems, but it 
should not be hard to make it run on other platforms. Dependencies:
- freeglut, http://freeglut.sourceforge.net/
- GLM, http://glm.g-truc.net/
- GLEW, http://glew.sourceforge.net/
- OpenGL 4 support (viewer only)
- graphviz, if MCG drawing is desired
- demos use eog (Eye Of Gnome) image viewer to view MCGs

/* -------------------------------------------------------------------------- */

HOW TO COMPILE

To compile: edit the makefile, adding system-specific access paths etc.
Make sure access paths to libraries and headers are provided if needed.
Make similar changes to the makefile in 'subd' subdirectory.

Then, run 'make'. In some cases, straightforward edits to #include 
statements may be needed.

/* -------------------------------------------------------------------------- */

DEMOS

Several demos are provided as bash scripts.
To run, compile the code, cd to demos and run any of the scripts.

/* -------------------------------------------------------------------------- */

EXECUTABLES

mdpc: computation of Morse decompositions, stable Morse decompositions
      and Morse connection graphs
msvis: visualization of the Morse decompositions and Morse connection
       graphs

In addition, vfsubd, which performs a PL-subdivision of the vector field,
is used in some of the provided deoms.

/* -------------------------------------------------------------------------- */

INPUT FILE FORMAT

Input files can be found in input directory. They are in a 
non-standard but very easy to read ASCII format. There are two 
types of input files the code currently accepts.

TYPE 1: FOR VECTOR FIELDS DEFINED ON A TRIANGLE MESH

These files contain only numbers (as opposed to files of type 2 described below 
that start with `n').

The first two numbers (both integers) are number of triangles (call it T) and 
vertices (V) in the mesh (the domain). Then, the file contains the triangle 
table, with each row specifying indices of vertices that form a triangle 
(note that vertices have indices in 0...V-1). There are T rows in the triangle 
table. Then, there comes the vertex table, coordinates of consecutive vertices 
(V rows, 3 floats per row). And finally, the vector data, i.e. vector values for
consecutive vertices (in the *V.t files; in this case, this part contains 
V tows with 3 floats in each row) OR triangles (for *T.t files; T rows, 3 floats
in each row). 

TYPE 2: FOR VECTOR FIELDS DEFINED ON POLYGONAL MESHES

The main goal of this format is to support 2D regular grids. mdpc will not work
correctly if there is a polygon with >8 edges.

The first character in the file is `n'.
The next two numbers (integers) are the number of 2D faces (F) and number of 
vertices (V).
Then, the file contains a F rows, each containing a tuple of vertex indices
(separated just by white spaces, no commas or parantheses!) bounding a 2D face.
The last entry of each row is -1, which plays the role of a separator.
Next, the file contains vertex coordinates (V rows, 3 floats per row).
Finally, the vector values are listed, one per vertex in *v.t files and one
per face in *f.t files.

/* -------------------------------------------------------------------------- */

DATA SOURCES

If you use the data sets in your research, please acknowledge their sources:

1) For the cooling jacket dataset, coolingJacket[TV].t:

- acknowledge: Robert S. Laramee at Swansea University

- cite: Robert S. Laramee, Christoph Garth, Helmut Doleisch, Juergen Schneider, 
Helwig Hauser, and Hans Hagen, 
"Visual Analysis and Exploration of Fluid Flow in a Cooling Jacket", 
in Proceedings of IEEE Visualization (IEEE Vis 2005), pages 623-630, 
October 23-28, 2005, Minneapolis, Minnesota.

2) For the gas and diesel engine dataset, gasEngine[TV].t and dieselEngine[TV].t:

- acknowledge: Robert S. Laramee at Swansea University

- cite: Robert S. Laramee, Daniel Weiskopf, Juergen Schneider, and Helwig Hauser, 
"Investigating Swirl and Tumble Flow with a Comparison of Visualization Techniques", 
in Proceedings of IEEE Visualization (IEEE Vis 2004), pages 51-58, October 15-19, 2004, 
Austin, Texas.

3) For most other files, bunny[12][TV].t, torus[12][TV].t, multicycles[TV].t, please cite

Guoning Chen, Konstantin Mischaikow, Robert S. Laramee, Pawel Pilarczyk, and Eugene Zhang,
"Vector Field Editing and Periodic Orbit Extraction Using Morse Decomposition", 
IEEE Transactions on Visualization and Computer Graphics, 2007, Vol 13(1), pp 769-785.

/* -------------------------------------------------------------------------- */

