
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

#ifndef __MSTYPE_H
#define __MSTYPE_H

#include <global.h>

/* ------------------------------------------------------ */

class mstype {
  unsigned char stability;    // 1 = attracting, 2 = repelling, 3 = neither
  int index;
  int index2;
  bool bdry;
  int size;

 public:

  mstype ( unsigned char s = 0, int i1 = 0, int i2 = 0, int sz = 0 );

  bool isnontrivial();
  bool istrivial();
  bool issaddle();
  bool isapo();
  bool isrpo();
  bool issink();
  bool issource();
  bool isattracting();
  bool isrepelling();

  int getindex();
  int getindex2();
  unsigned char getstability();
  bool getbdry();
  int getsize();

  void addtoindex ( int i );
  void addtoindex2 ( int i );
  void orstability ( unsigned char msk );
  void setbdry();
  void incsize();
};

/* ------------------------------------------------------ */

#endif
