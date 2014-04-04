

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

#ifndef __TRACKBALL_H
#define __TRACKBALL_H

#include <global.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

/* ---------------------------------------------- */

class trackball {

  glm::dmat4 R;     // finished rotations
  glm::dmat4 R0;    // current rotation
  glm::dmat4 R0R;
  glm::dmat3 NM;

  int winsize[2];  // window size

  bool ismousedown;  // true if dragging
  glm::dvec3 last;  // last mouse down event
  int lasti, lastj;

  glm::dvec3 _ij2xyz ( int i, int j );

 public:

  trackball ( int sx, int sy );
  void resize ( int sx, int sy );
  
  void mousedown ( int i, int j );
  void mousemove ( int i, int j );
  void mouseup ( int i, int j );

  glm::dmat4 mv();  // modelview matrix
  glm::dmat3 nm();  // normal matrix

  double *mvpointer();   // returns the modelview matrix
  double *nmpointer();   // returns the normal matrix

  bool isactive();   // is mouse button down?
};

/* ---------------------------------------------- */

#endif
