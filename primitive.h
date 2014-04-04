

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

#ifndef __IOMS_H
#define __IOMS_H

#include <global.h>
#include <fstream>
#include <vec3d.h>

/* ------------------------------------------------------ */

class primitive {
 public:
  vec3dd buf[4];
  unsigned char verts,stab;
  int id,index,index2;
  bool bdry;

  primitive ( int id, int index, int index2, 
	      unsigned char stab, bool bd,
	      vec3dd a );
  primitive ( int id, int index, int index2, 
	      unsigned char stab, bool bd,
	      vec3dd a, vec3dd b );
  primitive ( int id, int index, int index2,
	      unsigned char stab, bool bd,
	      vec3dd a, vec3dd b, vec3dd c );
  primitive ( int id, int index, int index2, 
	      unsigned char stab, bool bd,
	      vec3dd a, vec3dd b, vec3dd c, vec3dd d );
  primitive ( std::ifstream &ifs );
  void orient ( vec3dd n );   // orient consistently with n

  void save ( std::ofstream &ofs );

  operator bool();
};

/* ------------------------------------------------------ */

#endif
