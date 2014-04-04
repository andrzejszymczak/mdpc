

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

#include <pcstable.h>

#define EPS 1e-8

/* ------------------------------------------------------ */

bool pcstable::iststationary ( int i )
{
  return fstat[i];
}

/* ------------------------------------------------------ */

pcstable::pcstable ( double stability, const char *name, char type, bool BD ) :
  pcvf(name,type,BD), R(stability)
{
  int i,j;

  fstat = new bool[faces()];
  for ( i=0; i<faces(); i++ )
    fstat[i] = (f[i].norm()<=R);

  // we need to update aflow, rflow, eflow

  // aflow and rflow first
  for ( i=0; i<faces(); i++ )
    {
      mesh_element *ff = getface(i);
      for ( j=0; j<ff->faces; j+=2 )
	{
	  vec3dd ev = v[ff->face[j+1]->ID] - v[ff->face[(j+ff->faces-1)%(ff->faces)]->ID];
	  vec3dd ei = n[i]^ev;
	  double tst = ei*(v[ff->face[(j+3)%ff->faces]->ID]-v[ff->face[j+1]->ID]);
	  assert(tst!=0);
	  if (tst<0)
	    ei = -ei;
	  ei.normalize();
	  aflow[i][j] |= fstat[i] || (ei*f[i]<=R);
	  rflow[i][j] |= fstat[i] || (ei*f[i]>=-R);
	}
      aflow[i][ff->faces-1] = aflow[i][ff->faces-2] && aflow[i][0];
      rflow[i][ff->faces-1] = rflow[i][ff->faces-2] && rflow[i][0];
      for ( j=1; j<ff->faces-2; j+=2 )
	{
	  aflow[i][j] = aflow[i][j-1] && aflow[i][j+1];
	  rflow[i][j] = rflow[i][j-1] && rflow[i][j+1];
	}      
    }

  // time for eflow
  for ( i=0; i<edges(); i++ )
    {
      mesh_element *ee = mesh::getedge(i);
      vec3dd evc = v[ee->face[1]->ID]-v[ee->face[0]->ID];
      evc.normalize();
      switch(ee->cofaces)
	{
	case 1:
	  {
	    // boundary edge 
	    if (BD)
	      {
		if (evc*f[ee->coface[0]->ID]>=-R || fstat[ee->coface[0]->ID])
		  eflow[i] |= 1;
		if (evc*f[ee->coface[0]->ID]<= R || fstat[ee->coface[0]->ID])
		  eflow[i] |= 2;
	      }
	  }
	  break;

	case 2:
	  {
	    if (fstat[ee->coface[0]->ID] || fstat[ee->coface[1]->ID])
	      {
		eflow[i] |= 3;
		continue;
	      }
	    vec3dd an = n[ee->coface[1]->ID]+n[ee->coface[0]->ID];
	    double shrink = an.norm()/2;
	    an.normalize();
	    vec3dd f1 = f[ee->coface[0]->ID]-(f[ee->coface[0]->ID]*an)*an;
	    vec3dd f2 = f[ee->coface[1]->ID]-(f[ee->coface[1]->ID]*an)*an;
	    vec3dd ev = v[ee->face[1]->ID] - v[ee->face[0]->ID];
	    ev.normalize();
	    vec3dd pp = ev^an;
	    double x2 = f2*ev;
	    double x1 = f1*ev;
	    double y2 = (f2*pp)/shrink;
	    double y1 = (f1*pp)/shrink;
	    
	    int ct = 0;

	    if (fabs(y1)<=R)
	      {
		double delta = sqrt(R*R-y1*y1);
		if (x1-delta<=0) 
		  eflow[i] |= 2;
		if (x1+delta>=0)
		  eflow[i] |= 1;
		ct++;
	      }

	    if (fabs(y2)<=R)
	      {
		double delta = sqrt(R*R-y2*y2);
		if (x2-delta<=0) 
		  eflow[i] |= 2;
		if (x2+delta>=0)
		  eflow[i] |= 1;
		ct++;
	      }

	    if (ct==2) continue;

	    if (ct==0 && y1*y2>0)
	      continue;

	    double x12 = x2-x1;
	    double y12 = y2-y1;

	    double px = -y12;
	    double py = x12;
	    double nm = sqrt(y12*y12+x12*x12);
	    if (nm==0) continue;

	    px *= R/nm;
	    py *= R/nm;

	    double a1 = x1+px;
	    double b1 = y1+py;
	    double c1 = x2+px;
	    double d1 = y2+py;

	    double a2 = x1-px;
	    double b2 = y1-py;
	    double c2 = x2-px;
	    double d2 = y2-py;

	    // check for intersections of intervals a1b1--c1d1 and a2b2--c2d2
	    // with the x-axis

	    if (b1*d1<0)
	      {
		double xis = fabs(b1)*c1+fabs(d1)*a1;
		if (xis<=0) 
		  eflow[i] |= 2;
		if (xis>=0)
		  eflow[i] |= 1;
	      }

	    if (b2*d2<0)
	      {
		double xis = fabs(b2)*c2+fabs(d2)*a2;
		if (xis<=0) 
		  eflow[i] |= 2;
		if (xis>=0)
		  eflow[i] |= 1;
	      }
	  }
	  break;

	default:
	  assert(0);
	  break;
	}
    }
  
  testvec1 = new vec3dd[faces()];
  testvec2 = new vec3dd[faces()];
  
  for ( i=0; i<faces(); i++ )
    {
      if (fstat[i])
	continue;
      double sinalpha = R/f[i].norm();
      if (sinalpha>1) sinalpha = 1;
      double cosalpha = sqrt(1-sinalpha*sinalpha);
      if (cosalpha<EPS) cosalpha = EPS;
      double tanalpha = sinalpha/cosalpha;
      vec3dd tv1 = f[i]+(/*sinalpha**/tanalpha)*(n[i]^f[i]);
      vec3dd tv2 = f[i]-(/*sinalpha**/tanalpha)*(n[i]^f[i]);
      vec3dd tst1 = tv1^n[i];
      vec3dd tst2 = tv2^n[i];
      if (tst1*tv2>0) 
	tst1 = -tst1;
      if (tst2*tv1>0)
	tst2 = -tst2;
      testvec1[i] = tst1;
      testvec2[i] = tst2;

      //      cout << f[i] << " " << n[i] << " / " << tv1 << " " << tv2 << " / " << testvec1[i] << " " << testvec2[i] << 
      //	sinalpha << " " << cosalpha << " " << tanalpha << endl;
    }
}

/* ------------------------------------------------------ */

bool pcstable::connects ( int fce, int eix1, int eix2, 
				  double s1, double e1, 
				  double s2, double e2 )
{
  bool res;

  // if fce stationary, return false
  if (fstat[fce])
    return false;

  if ( (eix1&1) && (eix2&1) )
    {
      vec3dd w = v[getface(fce)->face[eix2]->ID]-v[getface(fce)->face[eix1]->ID];
      res = (w*testvec1[fce]<=0) && (w*testvec2[fce]<=0);
    }

  if ( (eix1&1) && !(eix2&1) )
    {
      // vertex and edge piece

      vec3dd e2s = edgepoint(getface(fce)->face[eix2]->ID,s2);
      vec3dd e2e = edgepoint(getface(fce)->face[eix2]->ID,e2);
      vec3dd p1 = v[getface(fce)->face[eix1]->ID];
      vec3dd v1 = e2s-p1;
      vec3dd v2 = e2e-p1;
      
      res =  !(
	       (v1*testvec1[fce]>=0 && v2*testvec1[fce]>=0) 
	       ||
	       (v1*testvec2[fce]>=0 && v2*testvec2[fce]>=0)
	       );
    }

  if ( !(eix1&1) && (eix2&1) )
    {
      // vertex and edge piece

      vec3dd e1s = edgepoint(getface(fce)->face[eix1]->ID,s1);
      vec3dd e1e = edgepoint(getface(fce)->face[eix1]->ID,e1);
      vec3dd p2 = v[getface(fce)->face[eix2]->ID];
      vec3dd v1 = p2-e1s;
      vec3dd v2 = p2-e1e;
      
      res = !(
	      (v1*testvec1[fce]>=0 && v2*testvec1[fce]>=0) 
	      ||
	      (v1*testvec2[fce]>=0 && v2*testvec2[fce]>=0)
	      );
    }

  if ( !(eix1&1) && !(eix2&1) )
    {
      // both edges
      vec3dd e2s = edgepoint(getface(fce)->face[eix2]->ID,s2);
      vec3dd e2e = edgepoint(getface(fce)->face[eix2]->ID,e2);
      vec3dd e1s = edgepoint(getface(fce)->face[eix1]->ID,s1);
      vec3dd e1e = edgepoint(getface(fce)->face[eix1]->ID,e1);
      vec3dd v1 = e2s-e1s;
      vec3dd v2 = e2s-e1e;
      vec3dd v3 = e2e-e1s;
      vec3dd v4 = e2e-e1e;

      res = !(
	      (v1*testvec1[fce]>=0 && v2*testvec1[fce]>=0 && v3*testvec1[fce]>=0 && v4*testvec1[fce]>=0)
	      ||
	      (v1*testvec2[fce]>=0 && v2*testvec2[fce]>=0 && v3*testvec2[fce]>=0 && v4*testvec2[fce]>=0)
	      );
    }

  //  if (res) cout << "Y" << endl;

  return res;
}

/* ------------------------------------------------------ */

pcstable::~pcstable()
{
  if (fstat) delete[] fstat;
  if (testvec1) delete[] testvec1;
  if (testvec2) delete[] testvec2;
  testvec1 = testvec2 = NULL;
  fstat = NULL;  
}

/* ------------------------------------------------------ */

/* ------------------------------------------------------ */

/* ------------------------------------------------------ */

/* ------------------------------------------------------ */
