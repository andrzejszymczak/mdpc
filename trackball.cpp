

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

#include <trackball.h>
#include <cmath>

using namespace glm;

/* ---------------------------------------------- */

trackball::trackball ( int sx, int sy ) : R(), R0(), R0R(), ismousedown(false)
{
  winsize[0] = sx;
  winsize[1] = sy;
  
}

/* ---------------------------------------------- */

void trackball::resize ( int sx, int sy )
{
  winsize[0] = sx;
  winsize[1] = sy;
}

/* ---------------------------------------------- */

dvec3 trackball::_ij2xyz ( int i, int j )
{
  float x = 2*((float)i)/winsize[0]-1;
  float y = 1-2*((float)j)/winsize[1];
  float x2y2 = 1-x*x-y*y;

  if (x2y2<0)
    {
      double l = sqrt(x*x+y*y);
      x = x/l;
      y = y/l;
      x2y2 = 0;
      
    }

  float z = sqrt(x2y2);

  return dvec3(x,y,z);      
}

/* ---------------------------------------------- */

void trackball::mousedown ( int i, int j )
{
  ismousedown = true;
  lasti = i;
  lastj = j;
  last = _ij2xyz(i,j);      
}

/* ---------------------------------------------- */

void trackball::mousemove ( int i, int j )
{
  if (i==lasti && j==lastj)
    {
      R0 = glm::mat4();
      return;
    }

  dvec3 cur = _ij2xyz(i,j);
  dvec3 axis = cross(last, cur);
  R0 = rotate(dmat4(),acos(dot(last,cur))/M_PI*180.0,axis);
  R0R = R0*R;
}

/* ---------------------------------------------- */

void trackball::mouseup ( int i, int j )
{
  if (!ismousedown)
    return;
  ismousedown = false;
  if (i==lasti && j==lastj)
    {
      R0 = mat4();
      return;
    }

  dvec3 cur = _ij2xyz(i,j);
  dvec3 axis = cross(last, cur);
  R0 = rotate(dmat4(),acos(dot(last,cur))/M_PI*180.0,axis);
  R = R0*R;
  R0 = dmat4();
  R0R = R;
}

/* ---------------------------------------------- */

double *trackball::mvpointer()
{
  return &R0R[0][0];
}

/* ---------------------------------------------- */

double *trackball::nmpointer()
{
  NM = dmat3(R0R);
  return &NM[0][0];
}

/* ---------------------------------------------- */

dmat4 trackball::mv()
{
  return R0R;
}

/* ---------------------------------------------- */

dmat3 trackball::nm()
{
  NM = dmat3(R0R);
  return NM;
}

/* ---------------------------------------------- */

bool trackball::isactive()
{
  return ismousedown;
}

/* ---------------------------------------------- */

