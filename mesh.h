

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

#ifndef __MESH_H
#define __MESH_H

#include <global.h>
#include <fstream>
#include <mesh_base.h>
#include <vec3d.h>

/* ------------------------------------------------------ */

class mesh : public mesh_base {

  void read_from_file ( std::ifstream &ifs, char format );
  // format: 't' - .t file
  // format: 'n' - polygonal mesh file
  // like .t but faces terminated with -1
  // file in n format starts with 'n'


  void compute_normals();

  vec3dd evec ( const mesh_element *f, int i, int j );

 protected:

  std::ifstream ifs;

  vec3dd *v;  // vertex coordinates
  vec3dd *n;  // unit normals for the faces

 public:

  mesh ( const char *name ); 
  ~mesh();

  void print_out();

  // returns point on edge eid 
  vec3dd edgepoint ( int eid, double t );  

  vec3dd vertex ( int i ); // vertex coordinate
  vec3dd normal ( int i ); // face normal
};

/* ------------------------------------------------------ */

#endif
