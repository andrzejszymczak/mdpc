

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

#ifndef __PCVFDISPLAY_H
#define __PCVFDISPLAY_H

#include <global.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <pcvf.h>

/* ----------------------------------------------------------- */

class pcvfdisplay : public pcvf {

  int ts;
  float *vc;
  float *nc;
  float *vec;

  double center[3];
  double maxbboxspan;

  double dist; // used for projection
  double aspect;

  double fv;   // field of view

  glm::dmat4 ntrans;   // normalizing transformation
  glm::dmat4 ftrans;   // forward translation amount
  glm::dmat4 proj;     // projection matrix

 public:

  // type='v' or 'f'
  pcvfdisplay ( const char *name, double fov, char type );

  virtual ~pcvfdisplay();

  void render();  // this call does not include any transformations
  // attribute indices: 0: location, 1: normal, 2: vector

  glm::dmat4 pm();       // projection matrix
  double *pmpointer(); 

  glm::dmat4 normt();    // normalizing transform
  double *nmpointer();

  glm::dmat4 ftr();      // forward translation
  double *ftrpointer();

  void setfov ( double a );
  double getfov();
  void resize ( int sx, int sy );

  double boxsize();
};

/* ----------------------------------------------------------- */

#endif
