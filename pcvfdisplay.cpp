
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
#include <pcvfdisplay.h>

using namespace glm;

/* ----------------------------------------------------------- */

void pcvfdisplay::render()
{
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,vc);  // vertex coord is attribute 0
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,nc);  // normals are 1
  glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,0,vec);  // vector is 2
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);
  glDrawArrays(GL_TRIANGLES,0,3*ts);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);  
}

/* ----------------------------------------------------------- */

pcvfdisplay::pcvfdisplay ( const char *name, double fov, char type ) : pcvf(name,type), fv(fov), aspect(1.0)
{
  // read mesh and build mesh data...
  int i;
  int trs = 0;

  double mx[3], mn[3], sz[3];
  for ( i=0; i<vertices(); i++ )
    {
      vec3dd v = vertex(i);
      if (!i)
	{
	  mx[0] = mn[0] = v[0];
	  mx[1] = mn[1] = v[1];
	  mx[2] = mn[2] = v[2];
	}
      else
	{
	  if (mx[0]<v[0]) mx[0] = v[0];
	  if (mx[1]<v[1]) mx[1] = v[1];
	  if (mx[2]<v[2]) mx[2] = v[2];
	  if (mn[0]>v[0]) mn[0] = v[0];
	  if (mn[1]>v[1]) mn[1] = v[1];
	  if (mn[2]>v[2]) mn[2] = v[2];	  
	}
    }

  center[0] = (mn[0]+mx[0])/2;
  center[1] = (mn[1]+mx[1])/2;
  center[2] = (mn[2]+mx[2])/2;

  sz[0] = mx[0]-center[0];
  sz[1] = mx[1]-center[1];
  sz[2] = mx[2]-center[2];

  if (sz[0]>=sz[1] && sz[0]>=sz[2])
    {
      maxbboxspan = sz[0];
    }
  if (sz[1]>=sz[0] && sz[1]>=sz[2])
    {
      maxbboxspan = sz[1];
    }
  if (sz[2]>=sz[1] && sz[2]>=sz[0])
    {
      maxbboxspan = sz[2];
    }
  
  // compute the transformations...

  dmat4 T = translate(glm::dmat4(),dvec3(-center[0],-center[1],-center[2]));
  dmat4 S = scale(glm::dmat4(),dvec3(1.0/maxbboxspan,1.0/maxbboxspan,1.0/maxbboxspan));
  ntrans = S*T;
  dist = 1/tan(fov/2.0/180.0*M_PI);
  proj = perspective(fov,aspect,dist-2,dist+2);
  ftrans = translate(glm::dmat4(),dvec3(0,0,-dist));

  vec3dd *pvn = new vec3dd[vertices()];  // per vertex normal
  bool *ma = new bool[vertices()];

  for ( i=0; i<vertices(); i++ )
    ma[i] = false;

  for ( i=0; i<faces(); i++ )
    {
      trs += getface(i)->faces/2-2;
      for ( int j=0; j<getface(i)->faces; j++ )
	if (getface(i)->face[j]->dimension==0)
	  pvn[getface(i)->face[j]->ID] += normal(getface(i)->ID);
    }

  for ( i=0; i<vertices(); i++ )
    pvn[i].normalize();

  for ( i=0; i<faces(); i++ )
    for ( int j=0; j<getface(i)->faces; j++ )
      if (getface(i)->face[j]->dimension==0)
	if (pvn[getface(i)->face[j]->ID]*normal(getface(i)->ID)<cos(M_PI/3))
	  ma[getface(i)->face[j]->ID] = true;
	
  vc = new float[9*trs];
  nc = new float[9*trs];
  vec = new float[9*trs];

  int cix = 0;
  for ( i=0; i<faces(); i++ )
    {
      for ( int j=5; j<getface(i)->faces; j+=2 )
	{
	  assert(cix<9*trs);
	  
	  vec3dd v1 = vertex(getface(i)->face[1]->ID);
	  vec3dd v2 = vertex(getface(i)->face[j-2]->ID);
	  vec3dd v3 = vertex(getface(i)->face[j]->ID);
	  vc[cix+0] = v2[0];
	  vc[cix+1] = v2[1];
	  vc[cix+2] = v2[2];
	  vc[cix+3] = v1[0];
	  vc[cix+4] = v1[1];
	  vc[cix+5] = v1[2];
	  vc[cix+6] = v3[0];
	  vc[cix+7] = v3[1];
	  vc[cix+8] = v3[2];

	  vec3dd n1 = ma[getface(i)->face[1]->ID] ? normal(i) : pvn[getface(i)->face[1]->ID];
	  vec3dd n2 = ma[getface(i)->face[j-2]->ID] ? normal(i) : pvn[getface(i)->face[j-2]->ID];
	  vec3dd n3 = ma[getface(i)->face[j]->ID] ? normal(i) : pvn[getface(i)->face[j]->ID]; 	
	  nc[cix+0] = n2[0];
	  nc[cix+1] = n2[1];
	  nc[cix+2] = n2[2];
	  nc[cix+3] = n1[0];
	  nc[cix+4] = n1[1];
	  nc[cix+5] = n1[2];
	  nc[cix+6] = n3[0];
	  nc[cix+7] = n3[1];
	  nc[cix+8] = n3[2];
	    
	  if (pvf)
	    {
	      vec3dd v1 = pvf[getface(i)->face[1]->ID];
	      vec3dd v2 = pvf[getface(i)->face[j-2]->ID];
	      vec3dd v3 = pvf[getface(i)->face[j]->ID];
	      vec3dd n1 = pvn[getface(i)->face[1]->ID];
	      vec3dd n2 = pvn[getface(i)->face[j-2]->ID];
	      vec3dd n3 = pvn[getface(i)->face[j]->ID];

	      v1 -= (v1*n1)*n1;
	      v2 -= (v2*n2)*n2;
	      v3 -= (v3*n3)*n3;

	      vec[cix+0] = v2[0];
	      vec[cix+1] = v2[1];
	      vec[cix+2] = v2[2];
	      vec[cix+3] = v1[0];
	      vec[cix+4] = v1[1];
	      vec[cix+5] = v1[2];
	      vec[cix+6] = v3[0];
	      vec[cix+7] = v3[1];
	      vec[cix+8] = v3[2];
	    }
	  else
	    {
	      vec3dd f = getvector(i);
	      vec[cix+0] = f[0];
	      vec[cix+1] = f[1];
	      vec[cix+2] = f[2];
	      vec[cix+3] = f[0];
	      vec[cix+4] = f[1];
	      vec[cix+5] = f[2];
	      vec[cix+6] = f[0];
	      vec[cix+7] = f[1];
	      vec[cix+8] = f[2];
	    }
	  cix += 9;
	}
    }

  ts = trs;
  assert(cix==9*trs);
}

/* ----------------------------------------------------------- */

dmat4 pcvfdisplay::pm()
{
  return proj;
} 

/* ----------------------------------------------------------- */

dmat4 pcvfdisplay::normt()
{
  return ntrans;
}

/* ----------------------------------------------------------- */

dmat4 pcvfdisplay::ftr()
{
  return ftrans;
}

/* ----------------------------------------------------------- */

double *pcvfdisplay::pmpointer()
{
  return &proj[0][0];
}

/* ----------------------------------------------------------- */

double *pcvfdisplay::nmpointer()
{
  return &ntrans[0][0];
}

/* ----------------------------------------------------------- */

double *pcvfdisplay::ftrpointer()
{
  return &ftrans[0][0];
}

/* ----------------------------------------------------------- */

void pcvfdisplay::resize ( int sx, int sy )
{
  aspect = float(sy)/sx;
  proj = perspective(fv,aspect,dist-2,dist+2);  
}

/* ----------------------------------------------------------- */

void pcvfdisplay::setfov ( double a )
{
  fv = a;
  proj = perspective(fv,aspect,dist-2,dist+2);    
}

/* ----------------------------------------------------------- */

double pcvfdisplay::getfov()
{
  return fv;
}

/* ----------------------------------------------------------- */

double pcvfdisplay::boxsize()
{
  return maxbboxspan;
}

/* ----------------------------------------------------------- */

pcvfdisplay::~pcvfdisplay()
{
  if (vc) delete[] vc;
  if (nc) delete[] nc;
  if (vec) delete[] vec;
}

/* ----------------------------------------------------------- */
