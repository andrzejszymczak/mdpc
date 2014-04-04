
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

#ifndef __PCSTABLE_H
#define __PCSTABLE_H

#include <global.h>
#include <pcvf.h>

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

class pcstable : public pcvf
{
  double R;
  bool *fstat;  // stationary faces...

  // test vectors for fast connects() implementation
  vec3dd *testvec1;
  vec3dd *testvec2;

 public:

  virtual bool iststationary ( int i );
  
  pcstable ( double stability, const char *name, char type = 'f', bool BD=true );
  virtual ~pcstable();

  // are two edges connected by flow through the face?
  virtual bool connects ( int fce, int eix1, int eix2, 
			  double s1 = 0, double e1 = 1, 
			  double s2 = 0, double e2 = 1 );
  
};

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

#endif
