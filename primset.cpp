

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

#include <GL/glew.h>
#include <GL/gl.h>
#include <primset.h>

using namespace std;

/* ----------------------------------------------------------- */

primset::primset ( const char *name, double bboxsize, double sizethr ) : 
  l(), e2d(NULL), ce2d(NULL), e1d(NULL), ce1d(NULL), e0d(NULL), ce0d(NULL),
  ts(0), es(0), vs(0), center(), valid(), diam(), maxid(-1)
{
  int i;
  ifstream ifs(name,ios::binary);

  double bound = (bboxsize>0) ? bboxsize*sizethr/100.0 : bboxsize;

  center.resize(1);
  valid.resize(1,false);
  diam.resize(1);

  if (ifs)
    {
      maxid = 0;
      primitive *p;
      do {
	p = new primitive(ifs);
	if (!*p)
	  break;
	l.push_back(p);
	if ((int)p->id>maxid)
	  {
	    maxid = p->id;
	    center.resize(maxid+1,vec3dd(0,0,0));
	    valid.resize(maxid+1,false);
	  }
	if (p->verts==4) 
	  ts+=2;
	else
	  if (p->verts==3)
	    ts++;
	  else
	    if (p->verts==2)
	      es++;
	    else
	      if (p->verts==1)
		{
		  center[p->id] = p->buf[0];
		  valid[p->id] = true;
		  vs++;
		}
      }
      while(1);
    }

  diam.resize(maxid+1,0.0);

  // now, build the arrays...

  e2d = new float[9*ts];
  ce2d = new float[9*ts];
  e1d = new float[6*es];
  ce1d = new float[6*es];
  e0d = new float[3*vs];
  ce0d = new float[3*vs];

  int ix2 = 0;
  int ix1 = 0;
  int ix0 = 0;

  for ( i=0; i<l.size(); ++i )
    if (valid[l[i]->id])
      {
	for ( int j=0; j<l[i]->verts; ++j )
	  {
	    double dst = (center[l[i]->id]-l[i]->buf[j]).norm();
	    if (dst>diam[l[i]->id])
	      diam[l[i]->id] = dst;
	  }
      }

  for ( i=0; i<l.size(); ++i )
    {
      float *c = color(l[i]);

      if (l[i]->verts==4)
	{
	  for ( int j=0; j<6; ++j )
	    {
	      ce2d[ix2+3*j+0] = c[0];
	      ce2d[ix2+3*j+1] = c[1];
	      ce2d[ix2+3*j+2] = c[2];
	    }

	  e2d[ix2+0] = l[i]->buf[0][0];
	  e2d[ix2+1] = l[i]->buf[0][1];
	  e2d[ix2+2] = l[i]->buf[0][2];
	  e2d[ix2+3] = l[i]->buf[1][0];
	  e2d[ix2+4] = l[i]->buf[1][1];
	  e2d[ix2+5] = l[i]->buf[1][2];
	  e2d[ix2+6] = l[i]->buf[2][0];
	  e2d[ix2+7] = l[i]->buf[2][1];
	  e2d[ix2+8] = l[i]->buf[2][2];
	  ix2 += 9;
	  e2d[ix2+0] = l[i]->buf[0][0];
	  e2d[ix2+1] = l[i]->buf[0][1];
	  e2d[ix2+2] = l[i]->buf[0][2];
	  e2d[ix2+3] = l[i]->buf[2][0];
	  e2d[ix2+4] = l[i]->buf[2][1];
	  e2d[ix2+5] = l[i]->buf[2][2];
	  e2d[ix2+6] = l[i]->buf[3][0];
	  e2d[ix2+7] = l[i]->buf[3][1];
	  e2d[ix2+8] = l[i]->buf[3][2];
	  ix2 += 9;
	}
      else
	if (l[i]->verts==3)
	  {
	    for ( int j=0; j<3; j++ )
	      {
		ce2d[ix2+3*j+0] = c[0];
		ce2d[ix2+3*j+1] = c[1];
		ce2d[ix2+3*j+2] = c[2];
	      }

	    e2d[ix2+0] = l[i]->buf[0][0];
	    e2d[ix2+1] = l[i]->buf[0][1];
	    e2d[ix2+2] = l[i]->buf[0][2];
	    e2d[ix2+3] = l[i]->buf[1][0];
	    e2d[ix2+4] = l[i]->buf[1][1];
	    e2d[ix2+5] = l[i]->buf[1][2];
	    e2d[ix2+6] = l[i]->buf[2][0];
	    e2d[ix2+7] = l[i]->buf[2][1];
	    e2d[ix2+8] = l[i]->buf[2][2];
	    ix2 += 9;
	  }
	else
	  if (l[i]->verts==2)
	    {
	      for ( int j=0; j<2; j++ )
		{
		  ce1d[ix1+3*j+0] = c[0];
		  ce1d[ix1+3*j+1] = c[1];
		  ce1d[ix1+3*j+2] = c[2];
		}
	      e1d[ix1+0] = l[i]->buf[0][0];
	      e1d[ix1+1] = l[i]->buf[0][1];
	      e1d[ix1+2] = l[i]->buf[0][2];
	      e1d[ix1+3] = l[i]->buf[1][0];
	      e1d[ix1+4] = l[i]->buf[1][1];
	      e1d[ix1+5] = l[i]->buf[1][2];
	      ix1 += 6;
	    }
	  else
	    if (l[i]->verts==1 && valid[l[i]->id] && diam[l[i]->id]<bound)
	      {
		
		ce0d[ix0+0] = c[0];
		ce0d[ix0+1] = c[1];
		ce0d[ix0+2] = c[2];
		e0d[ix0+0] = l[i]->buf[0][0];
		e0d[ix0+1] = l[i]->buf[0][1];
		e0d[ix0+2] = l[i]->buf[0][2];
		ix0 += 3;
	      }
    }

  vs = ix0/3;
}

/* ----------------------------------------------------------- */

float * primset::color ( primitive *p )
{
  float alpha = 1.0;

  static float c[10][4] =
    {
      {0.6,0.0,0.6,alpha},
      {0.0,0.7,0.0,alpha},
      {.3,.3,.3,alpha},  // originally A=alpha/3 magenta (0,0)
      {0,.8,0,alpha},        // green (1 -)
      {.8,0,0,alpha}, // red(1 +)
      {0,0,.8,alpha}, // blue (-1 0)
      {0,.5,0,alpha}, // green (0 -)
      {.5,0,0,alpha}, // red (0 +)
      {0,0,0,alpha},  // black other
      {.5*.542,.5*.27,.5*.0742, alpha}  // originally A=0.5
    };

  return c[code(p)];
}

/* ----------------------------------------------------------- */

int primset::code ( primitive *p )
{
  if (p->stab==255)
    return 9;
  if (p->stab==3 && p->index==0)
    return 2;
  if (p->stab==3 && p->index==-1) 
    return 5;
  if (p->stab==1)
    {
      if (p->index==0)
	return 6;
      if (p->index==1)
	return 3;
    }
  if (p->stab==2)
    {
      if (p->index==0)
	return 7;
      if (p->index==1)
	return 4;
    }
  return 8;  
}

/* ----------------------------------------------------------- */

void primset::render ( int maxchunk )
{
  if (ts)
    {
      glPolygonMode(GL_FRONT,GL_LINE);
      glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,e2d);
      glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,ce2d);
      glEnableVertexAttribArray(0);
      glEnableVertexAttribArray(1);

      int i = 0;

      while(1)
	{
	  glDrawArrays(GL_TRIANGLES,3*i*maxchunk,min(3*maxchunk,3*ts-3*i*maxchunk));
	  i++;
	  if (i*maxchunk>=ts)
	    break;
	}

      //      glDrawArrays(GL_TRIANGLES,0,3*ts);
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);

      glPolygonMode(GL_FRONT,GL_FILL);
      glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,e2d);
      glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,ce2d);
      glEnableVertexAttribArray(0);
      glEnableVertexAttribArray(1);

      i = 0;

      while(1)
	{
	  glDrawArrays(GL_TRIANGLES,3*i*maxchunk,min(3*maxchunk,3*ts-3*i*maxchunk));
	  i++;
	  if (i*maxchunk>=ts)
	    break;
	}

      //      glDrawArrays(GL_TRIANGLES,0,3*ts);
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);
    }

  if (es)
    {
      glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,e1d);
      glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,ce1d);
      glEnableVertexAttribArray(0);
      glEnableVertexAttribArray(1);

      int i = 0;

      while(1)
	{
	  glDrawArrays(GL_LINES,2*i*maxchunk,min(2*maxchunk,2*es-2*i*maxchunk));
	  i++;
	  if (i*maxchunk>=es)
	    break;
	}

      //      glDrawArrays(GL_LINES,0,2*es);
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);
    }

  if (vs)
    {
      glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,e0d);
      glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,ce0d);
      glEnableVertexAttribArray(0);
      glEnableVertexAttribArray(1);

      int i = 0;

      while(1)
	{
	  glDrawArrays(GL_POINTS,i*maxchunk,min(maxchunk,vs-i*maxchunk));
	  i++;
	  if (i*maxchunk>=vs)
	    break;
	}
      glDrawArrays(GL_POINTS,0,vs);
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);
    }
}

/* ----------------------------------------------------------- */

primset::~primset()
{
  if (e1d) delete[] e1d;
  if (ce1d) delete[] ce1d;
  if (e2d) delete[] e1d;
  if (ce2d) delete[] ce1d;
  if (e0d) delete[] e0d;
  if (ce0d) delete[] ce0d;
  for ( int i=0; i<l.size(); i++ )
    {
      delete l[i];
      l[i] = NULL;
    }
}

/* ----------------------------------------------------------- */
