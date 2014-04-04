

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

#ifndef __VFIELD_BASE_H
#define __VFIELD_BASE_H

#include <global.h>
#include <mesh.h>
#include <vec3d.h>
#include <fstream>

/* ------------------------------------------------------ */

class vfield_base : public mesh {
  
 public:

  // is vertex with ID i stationary?  
  virtual bool isstationary ( int i ) = 0;

  // fixed point index
  virtual int index ( int i ) = 0;

  // fixed point index assuming domain is repelling
  virtual int index2 ( int i ) = 0;

  // spiral flow (generically, true for spiral sinks/sources)
  virtual bool isspiral ( int i ) = 0;

  // does an edge i have flow up, i.e. vertex0->vertex1? down means the opposite
  virtual bool hasflowup ( int i ) = 0;
  virtual bool hasflowdown ( int i ) = 0;

  // is triangle stationary?
  virtual bool iststationary( int i ) = 0;

  // does edge or vertex attract/repel the flow?
  virtual bool attracts_flow ( int fce, int eix ) = 0;
  virtual bool repels_flow ( int fce, int eix ) = 0;

  // are two edges connected by flow through the face?
  virtual bool connects ( int fce, int eix1, int eix2, 
			  double s1 = 0, double e1 = 1, 
			  double s2 = 0, double e2 = 1 ) = 0;
  

  vfield_base ( const char *name );
  virtual ~vfield_base();
};

/* ------------------------------------------------------ */

/* ------------------------------------------------------ */

#endif
