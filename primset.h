

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

#ifndef __PRIMSET_H
#define __PRIMSET_H

#include <global.h>
#include <primitive.h>
#include <vector>

/* ----------------------------------------------------------- */

class primset {

  std::vector<primitive*> l;
  std::vector<vec3dd> center;
  std::vector<double> diam;
  std::vector<bool> valid;

  int maxid;

  int ts;
  float *e2d;  // triangle array
  float *ce2d; // colors

  int es;
  float *e1d;  // edge array
  float *ce1d; // colors

  int vs;
  float *e0d;  // points
  float *ce0d; // colors

  float * color ( primitive *p );
  int code ( primitive *p );

 public:

  // bboxsize: bounding box size, sizethr: threshold for Morse set diameter estimate
  // for which points are drawn
  primset ( const char *name, double bboxsize, double sizethr );

  ~primset();

  // attribute 0: location, attribute 1: color
  void render ( int maxchunk = 65536 );

};

/* ----------------------------------------------------------- */

#endif
