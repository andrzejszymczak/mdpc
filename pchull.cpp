
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

#include <pchull.h>

using namespace std;

/* ------------------------------------------------------ */

static bool _dotpositive ( const vec3dd &v, vector<vec3dd> *s )
{
  for ( int i=0; i<s->size(); i++ )
    if (v*(*s)[i]>=0)
      return true;
  return false;
}

static bool _dotnegative ( const vec3dd &v, vector<vec3dd> *s )
{
  for ( int i=0; i<s->size(); i++ )
    if (v*(*s)[i]<=0)
      return true;
  return false;
}

/* ------------------------------------------------------ */

static bool _samesign ( double a, double b, double c )
{
  return (a>=0 && b>=0 && c>=0) || (a<=0 && b<=0 && c<=0);
}

static bool _zeroinslab ( const vec3dd &n, const vec3dd &a, const vec3dd &b, const vec3dd &c )
{
  vec3dd Nab = (b-a)^n;
  vec3dd Nbc = (c-b)^n;
  vec3dd Nca = (a-c)^n;
  double zNab = (-a)*Nab;
  double zNbc = (-b)*Nbc;
  double zNca = (-c)*Nca;
  return _samesign(zNab,zNbc,zNca);
}

static bool _zeroinhull ( const vec3dd &n, vector<vec3dd> *s )
{
  int i,j,k;
  for ( int i=0; i<s->size(); i++ )
    for ( int j=i+1; j<s->size(); j++ )
      for ( int k=j+1; k<s->size(); k++ )
	if (_zeroinslab (n,(*s)[i],(*s)[j],(*s)[k]))
	  return true;
  return false;
}

/* ------------------------------------------------------ */

static unsigned char _testdir ( const vec3dd &tv, const vec3dd &p, vector<vec3dd> *s1, vector<vec3dd> *s2 )
{
  unsigned char res = 0;

  int i,j;

  vec3dd o = tv^p;

  // project vectors in s1 and s2...
  vector<vec3dd> p1,p2;
  for ( i=0; i<s1->size(); i++ )
    p1.push_back((*s1)[i]-((*s1)[i]*p)*p);
  for ( i=0; i<s2->size(); i++ )
    p2.push_back((*s2)[i]-((*s2)[i]*p)*p);

  for ( i=0; i<p1.size(); i++ )
    {
      double x1 = tv*p1[i];
      double y1 = o*p1[i];
      for ( j=0; j<p2.size(); j++ )
	{
	  double x2 = tv*p2[i];
	  double y2 = o*p2[i];
	  if (y1*y2<=0)
	    {
	      double xis = fabs(x1)*y2+fabs(x2)*y1;
	      if (xis>=0)
		res |= 1;
	      if (xis<=0)
		res |= 2;
	    }
	}
    }
}

/* ------------------------------------------------------ */

static void _findextreme ( const vec3dd &n, vector<vec3dd> *s, vec3dd *v, vec3dd *w )
{
  int i,j;
  int *ctp = new int[s->size()];
  for ( i=0; i<s->size(); i++ )
    ctp[i] = 0;

  for ( i=0; i<s->size(); i++ )
    {
      vec3dd tv = n^(*s)[i];
      for ( j=0; j<s->size(); j++ )
	if (i!=j)
	  if ((tv*(*s)[j])>=0)
	    ++(ctp[i]);
    }

  int min,max;
  min = max = 0;
  int minval,maxval;
  minval = maxval = ctp[0];

  for ( i=1; i<s->size(); i++ )
    {
      if (minval>ctp[i])
	{
	  minval = ctp[i];
	  min = i;
	}
      if (maxval<ctp[i])
	{
	  maxval = ctp[i];
	  max = i;
	}
    }

  *v = (*s)[max];
  *w = (*s)[min];

  delete[] ctp;
}

/* ------------------------------------------------------ */

bool pchull::iststationary ( int i )
{
  return fstat[i];
}

/* ------------------------------------------------------ */

pchull::pchull ( double wt, const char *name, char type, bool BD ) :
  pcvf(name,type,BD), weight(wt)
{
  int i,j;

  F0 = new vector<vec3dd>[faces()];

  for ( i=0; i<faces(); ++i )
    {
      mesh_element *ff = getface(i);
      for ( j=0; j<ff->faces; j+=2 )
	if (ff->face[j]->cofaces==2)
	  {
	    int af;  // adjacent face
	    vec3dd ev = v[ff->face[j+1]->ID] - v[ff->face[(j+ff->faces-1)%(ff->faces)]->ID];
	    vec3dd f0 = f[ af = ff->face[j]->coface[
						    (ff->face[j]->coface[0]==ff) ? 1 : 0
						    ]->ID];
	    F0[i].push_back(f[af]-(f[af]*n[i])*n[i]);
	  }
      F0[i].push_back(f[i]);

    }

  initialize(BD);
}

/* ------------------------------------------------------ */

void pchull::initialize ( bool BD )
{
  int i,j;

  // now compute F
  F = new vector<vec3dd>[faces()];

  for ( i=0; i<faces(); i++ )
    {
      for ( j=0; j<F0[i].size(); j++ )
	F[i].push_back(f[i]+weight*(F0[i][j]-f[i]));

      for ( j=0; j<F[i].size(); j++ )
	F[i][j] = F[i][j]-(F[i][j]*normal(i))*normal(i);
    }

  // fstat...
  fstat = new bool[faces()];
  for ( i=0; i<faces(); i++ )
    fstat[i] = _zeroinhull(normal(i),&F[i]);

  // update aflow and rflow 
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
	  aflow[i][j] |= fstat[i] || _dotnegative(ei,&F[i]);
	  rflow[i][j] |= fstat[i] || _dotpositive(ei,&F[i]);
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
		if (_dotpositive(evc,&F[ee->coface[0]->ID]) || fstat[ee->coface[0]->ID])
		  eflow[i] |= 1;
		if (_dotnegative(evc,&F[ee->coface[0]->ID]) || fstat[ee->coface[0]->ID])
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
	    vec3dd ev = v[ee->face[1]->ID] - v[ee->face[0]->ID];
	    an.normalize();
	    eflow[i] |=  _testdir(ev,an,&F[ee->coface[0]->ID],&F[ee->coface[1]->ID]);
	  }
	  break;

	default:
	  assert(0);
	  break;
	}
    }
    
  // finally, the test vectors...
  testvec1 = new vec3dd[faces()];
  testvec2 = new vec3dd[faces()];
  
  for ( i=0; i<faces(); i++ )
    {
      if (fstat[i])
	continue;

      vec3dd tv1;
      vec3dd tv2;
      _findextreme(n[i],&F[i],&tv1,&tv2);
      vec3dd tst1 = tv1^n[i];
      vec3dd tst2 = tv2^n[i];
      if (tst1*tv2>0) 
	tst1 = -tst1;
      if (tst2*tv1>0)
	tst2 = -tst2;
      testvec1[i] = tst1;
      testvec2[i] = tst2;
    }
}

/* ------------------------------------------------------ */

bool pchull::connects ( int fce, int eix1, int eix2, 
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

pchull::~pchull()
{
  if (fstat) delete[] fstat;
  if (testvec1) delete[] testvec1;
  if (testvec2) delete[] testvec2;
  if (F0) delete[] F0;
  if (F) delete[] F;
  F = NULL;
  F0 = NULL;
  fstat = NULL;
  testvec1 = testvec2 = NULL;
}

/* ------------------------------------------------------ */
