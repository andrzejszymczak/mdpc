

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

#include <primitive.h>

using namespace std;

/* ------------------------------------------------------ */

primitive::primitive ( int i, int ix, int ix2, unsigned char s, bool bb, vec3dd a ) : 
  verts(1), stab(s), index(ix), index2(ix2), id(i), bdry(bb)
{
  buf[0] = a;
}

primitive::primitive ( int i, int ix, int ix2, unsigned char s, bool bb, vec3dd a, vec3dd b ) : 
  verts(2), stab(s), index(ix), index2(ix2), id(i), bdry(bb)
{
  buf[0] = a;
  buf[1] = b;
}

primitive::primitive ( int i, int ix, int ix2, unsigned char s, bool bb, vec3dd a, vec3dd b, vec3dd c ) : 
  verts(3), stab(s), index(ix), index2(ix2), id(i), bdry(bb)
{
  buf[0] = a;
  buf[1] = b;
  buf[2] = c;
}

primitive::primitive ( int i, int ix, int ix2, unsigned char s, bool bb, vec3dd a, vec3dd b, vec3dd c, vec3dd d ) : 
  verts(4), stab(s), index(ix), index2(ix2), id(i), bdry(bb)
{
  buf[0] = a;
  buf[1] = b;
  buf[2] = c;
  buf[3] = d;
  // remove repetitions
  if (d==c || d==b || d==a)
    {
      verts--;
      return;
    }
  if (c==b || c==a)
    {
      buf[2] = buf[3];
      verts--;
      return;
    }
  if (a==b)
    {
      buf[1] = buf[3];
      verts--;
      return;
    }
}

/* ------------------------------------------------------ */

void primitive::orient ( vec3dd n )
{
  if (verts<=2) return;
  if (verts==3)
    if (((buf[1]-buf[0])^(buf[2]-buf[0]))*n<0)
      {
	vec3dd tmp = buf[1];
	buf[1] = buf[2];
	buf[2] = tmp;
      }
  /* -- old version
  if (verts==4)
    {
      if (((buf[1]-buf[0])^(buf[2]-buf[0]))*n<0)
	{
	  vec3dd tmp = buf[0];
	  buf[0] = buf[1];
	  buf[1] = tmp;
	}      
      if (((buf[2]-buf[1])^(buf[3]-buf[1]))*n<0)
	{
	  vec3dd tmp = buf[2];
	  buf[2] = buf[3];
	  buf[3] = tmp;
	}
    }
  */

  if (verts==4)
    {
      if ( ( ((buf[2]-buf[0])^(buf[3]-buf[0])) + ((buf[2]-buf[1])^(buf[3]-buf[1])) )*n < 0 )
	{
	  vec3dd tmp = buf[2];
	  buf[2] = buf[3];
	  buf[3] = tmp;
	}      
      if ( ( ((buf[0]-buf[2])^(buf[1]-buf[2])) + ((buf[0]-buf[3])^(buf[1]-buf[3])) )*n < 0)
	{
	  vec3dd tmp = buf[0];
	  buf[0] = buf[1];
	  buf[1] = tmp;
	}
     
    }
}

/* ------------------------------------------------------ */

void primitive::save ( ofstream &ofs )
{
  ofs.write((char*)&id,sizeof(int));
  ofs.write((char*)&index,sizeof(int));
  ofs.write((char*)&index2,sizeof(int));
  ofs.write((char*)&stab,sizeof(unsigned char));
  ofs.write((char*)&verts,sizeof(unsigned char));
  ofs.write((char*)&bdry,sizeof(bool));
  for ( int i=0; i<verts; i++ )
    ofs.write((char*)buf[i].pointer(),3*sizeof(double));
}

/* ------------------------------------------------------ */

primitive::primitive ( ifstream &ifs )
{
  ifs.read((char*)&id,sizeof(int));
  if (ifs.eof())
    {
      id = -1;
      return;
    }
  ifs.read((char*)&index,sizeof(int));
  ifs.read((char*)&index2,sizeof(int));
  ifs.read((char*)&stab,sizeof(unsigned char));
  ifs.read((char*)&verts,sizeof(unsigned char));
  ifs.read((char*)&bdry,sizeof(bool));
  assert(verts<=4);
  for ( int i=0; i<verts; i++ )
    ifs.read((char*)buf[i].pointer(),3*sizeof(double));
  
}

/* ------------------------------------------------------ */

primitive::operator bool()
{
  return id>=0;
}

/* ------------------------------------------------------ */
