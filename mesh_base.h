
/*
 * MDPC (Morse Decompositions for Piecewise Constant vector fields)
 * Copyright (c) 2012  Andrzej Szymczak
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"), 
 * to deal in the Software without restriction, including without limitation 
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 * and/or sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 * ANDRZEJ SZYMCZAK BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
 * SOFTWARE.
 *
*/

#ifndef __MESH_BASE_H
#define __MESH_BASE_H

#include <global.h>
#include <vector>

/* ------------------------------------------------------ */

class mesh_element {

 public:

  int ID;
  int dimension;

  int faces;
  mesh_element** face;   // ALL faces, in CCW order for 2D face 
  // for a 2D face 0D (vertices) are odd, 1D (edges) are even

  int cofaces;
  mesh_element** coface; // ALL cofaces; for a vertex, in order around it

  mesh_element ( int id, int dim );
  ~mesh_element();

  void set_faces ( int n, mesh_element **f );
  void set_cofaces ( int n, mesh_element **cf );
  bool isboundary();
  int find_face_index ( mesh_element *e );  // which face is this?

  bool isisolatedvertex();

  void print_out();
};

/* ------------------------------------------------------ */

// assumes manifold mesh, with planar 2D cells

class mesh_base {

  std::vector<std::vector<int>*> *workspace;

 protected:

  int vtcs, es, fcs;
  std::vector<mesh_element*> mel;  // mesh elements in decreasing dimension order (2,1,0)

  void add_2Dmel ( std::vector<int> *verts );   // adds a 2D mesh element
  void finalize();    // finalizes the datastructure; 
                      // in particular, computes faces, cofaces etc


 public:

  mesh_base ( );      
  ~mesh_base();

  mesh_element * get ( int i ); // get mesh element i
  mesh_element * getedge ( int i );
  mesh_element * getvertex ( int i );
  mesh_element * getface ( int i ); 

  void print_out();

  int vertices();
  int edges();
  int faces();
  int mesh_elements();
};

/* ------------------------------------------------------ */

#endif
