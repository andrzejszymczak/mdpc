
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

#include <pcenv.h>
#include <iostream>

using namespace std;

/* ------------------------------------------------------ */

pcenv::~pcenv() {}

/* ------------------------------------------------------ */

pcenv::pcenv ( double weight, const char *name, bool BD ) :
  pchull(weight,name,'v',BD)
{
  F0 = new vector<vec3dd>[faces()];

  for ( int i=0; i<faces(); i++ )
    {
      mesh_element *ff = getface(i);
      for ( int j=1; j<ff->faces; j+=2 )
	{
	  vec3dd vv = pvf[ff->face[j]->ID];
	  F0[i].push_back(vv-(vv*n[i])*n[i]);
	}
    }
  initialize(BD);
}

/* ------------------------------------------------------ */
