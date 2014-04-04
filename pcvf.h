

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

#ifndef __PCVF_H
#define __PCVF_H

#include <global.h>
#include <vfield_base.h>

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

class pcvf : public vfield_base {

 protected:

  vec3dd *f;
  vec3dd *pvf;   // per vertex vector values; NULL if file contains only per-face values

  // precomputed stuff ...

  bool **aflow;           // attracts flow [edge rel face]?
  bool **rflow;           // attracts flow [edge rel face]?

  unsigned char *eflow;   // bit 0: up, bit 1: down

  // projections of vertices along face's vector
  double **proj0;
  double **proj1;

  // vertex data
  bool *isstat;  // stationary or not
  int *indx;     // fixed point index; for boundary vertices assumes flow converging toward the domain outside it
  int *indx2;    //  for boundary vertices assumes flow escaping from the domain
  bool *spiral;  // spiral or not

 public:

  // type='f' if per-face vector value,
  // type='v' if per vertex vector value
  // bd=true: don't allow flow across boundary edges
  //           'closed' dynamical system
  // bd=false: allow flow across boundary edges
  pcvf ( const char *name, char type = 'f', bool BD = true );

  virtual ~pcvf();

  vec3dd getvector ( int i );

  // see vfield_base for comments on fuctions below
  virtual bool isstationary ( int i );
  virtual int index ( int i );
  virtual int index2 ( int i );
  virtual bool hasflowup ( int i );
  virtual bool hasflowdown ( int i );
  virtual bool isspiral ( int i );
  virtual bool attracts_flow ( int fce, int eix );
  virtual bool repels_flow ( int fce, int eix );
  virtual bool iststationary ( int i );

  // are two edges connected by flow through the face?
  virtual bool connects ( int fce, int eix1, int eix2, 
			  double s1 = 0, double e1 = 1, 
			  double s2 = 0, double e2 = 1 );


  void print_out();
};

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

#endif
