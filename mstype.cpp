

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

#include <mstype.h>


/* ------------------------------------------------------ */

mstype::mstype ( unsigned char s, int i1, int i2, int sz ) : stability(s), index(i1), index2(i2), size(sz), bdry(false)
{}

/* ------------------------------------------------------ */

bool mstype::isnontrivial()
{
  return !istrivial();
}

/* ------------------------------------------------------ */

bool mstype::istrivial()
{
  return stability==3 && index==0 && index2==0;
}

/* ------------------------------------------------------ */

bool mstype::issaddle()
{
  return stability==3 && index==-1 && index2==-1;
}

/* ------------------------------------------------------ */

bool mstype::isapo()
{
  return !index && !index2 && stability==1;
}

/* ------------------------------------------------------ */

bool mstype::isrpo()
{
  return !index && !index2 && stability==2;
}

/* ------------------------------------------------------ */

bool mstype::issource()
{
  return index==1 && index2==1 && stability==2;
}

/* ------------------------------------------------------ */

bool mstype::issink()
{
  return index==1 && index2==1 && stability==1;
}

/* ------------------------------------------------------ */

bool mstype::isrepelling()
{
  return stability==2;
}

/* ------------------------------------------------------ */

bool mstype::isattracting()
{
  return stability==1;
}

/* ------------------------------------------------------ */

int mstype::getindex()
{
  return index;
}

/* ------------------------------------------------------ */

int mstype::getindex2()
{
  return index2;
}

/* ------------------------------------------------------ */

void mstype::addtoindex ( int i )
{
  index += i;
}

/* ------------------------------------------------------ */

void mstype::addtoindex2 ( int i )
{
  index2 += i;
}

/* ------------------------------------------------------ */

void mstype::orstability ( unsigned char msk )
{
  stability |= msk;
}

/* ------------------------------------------------------ */

void mstype::setbdry()
{
  bdry = true;
}

/* ------------------------------------------------------ */

void mstype::incsize()
{
  size++;
}

/* ------------------------------------------------------ */

unsigned char mstype::getstability()
{
  return stability;
}

/* ------------------------------------------------------ */

bool mstype::getbdry()
{
  return bdry;
}

/* ------------------------------------------------------ */

int mstype::getsize()
{
  return size;
}

/* ------------------------------------------------------ */
