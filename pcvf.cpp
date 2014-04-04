
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

#include <cstdlib>
#include <pcvf.h>
#include <iostream>

using namespace std;

/* ------------------------------------------------------ */

pcvf::pcvf ( const char *name, char type, bool BD ) :
  vfield_base(name), pvf(NULL)
{
  int i;

  f = new vec3dd[faces()];

  if (type=='f')
    {
      for ( i=0; i<faces(); i++ )
	{
	  ifs >> f[i][0] >> f[i][1] >> f[i][2];   

	  if (ifs.eof())
	    {
	      cout << "Premature end of file: " << name << endl;
	      exit(1);
	    }
	}
    }
  else
    if (type=='v')
      {
	pvf = new vec3dd[vertices()];
	for ( i=0; i<vertices(); ++i )
	  {
	    ifs >> pvf[i][0] >> pvf[i][1] >> pvf[i][2];
	    if (ifs.eof())
	      {
		cout << "Premature end of file: " << name << endl;
		exit(1);
	      }
	  }

	for ( i=0; i<faces(); ++i )
	  {
	    for ( int j=1; j<getface(i)->faces; j+=2 )
	      f[i] += pvf[getface(i)->face[j]->ID];
	    f[i] *= (2.0/getface(i)->faces);
	  }
      }
    else
      {
	cout << "Unknown type in pcvf::pcvf" << endl;
	exit(1);
      }

  // project to face
  for ( i=0; i<faces(); i++ )
    f[i] -= (f[i]*n[i])*n[i];
  
  // fill aflow...
  aflow = new bool*[faces()];
  rflow = new bool*[faces()];
  for ( i=0; i<faces(); i++ )
    {
      int j;

      mesh_element *ff = getface(i);
      aflow[i] = new bool[ff->faces];
      rflow[i] = new bool[ff->faces];
      int crossings = 0;  // this is for consistency check...

      while(1)
	{
	  for ( j=0; j<ff->faces; j+=2 )
	    {
	      vec3dd ev = v[ff->face[j+1]->ID] - v[ff->face[(j+ff->faces-1)%(ff->faces)]->ID];
	      vec3dd ei = n[i]^ev;
	      double tst = ei*(v[ff->face[(j+3)%ff->faces]->ID]-v[ff->face[j+1]->ID]);
	      assert(tst!=0);
	      if (tst<0)
		ei = -ei;
	      aflow[i][j] = (ei*f[i]<0);
	      rflow[i][j] = !aflow[i][j];
	      if (j && (aflow[i][j]^aflow[i][j-2])) crossings++;
	    }
	  if (aflow[i][0]^aflow[i][ff->faces-2]) crossings++;

	  if (crossings!=2)
	    {
	      cout << "Face with more than 2 attract/repel switches! Trying a random perturbation..." << endl;
	      cout << "Index: " << i << " " << f[i];
	      f[i] = 1e-6*randomv3<double>();
	      cout << "--->" << f[i] << endl;
	    }
	  else
	    break;
	}

      aflow[i][ff->faces-1] = aflow[i][ff->faces-2] && aflow[i][0];
      rflow[i][ff->faces-1] = rflow[i][ff->faces-2] && rflow[i][0];
      for ( j=1; j<ff->faces-2; j+=2 )
	{
	  aflow[i][j] = aflow[i][j-1] && aflow[i][j+1];
	  rflow[i][j] = rflow[i][j-1] && rflow[i][j+1];
	}
    }

  // fill eflow...
  eflow = new unsigned char[edges()];
  for ( i=0; i<edges(); i++ )
    {
      mesh_element *ee = getedge(i);
      vec3dd evc = v[ee->face[1]->ID]-v[ee->face[0]->ID];
      switch(ee->cofaces)
	{
	case 1:
	  {
	    // boundary edge - has flow for sure!
	    eflow[i] = 0;
	    if (BD)
	      {
		if (evc*f[ee->coface[0]->ID]>0)
		  eflow[i] = 1;
		else
		  eflow[i] = 2;
	      }
	  }
	  break;

	case 2:
	  {
	    vec3dd an = n[ee->coface[1]->ID]+n[ee->coface[0]->ID];
	    an.normalize();
	    vec3dd f1 = f[ee->coface[0]->ID]-(f[ee->coface[0]->ID]*an)*an;
	    vec3dd f2 = f[ee->coface[1]->ID]-(f[ee->coface[1]->ID]*an)*an;
	    vec3dd ev = v[ee->face[1]->ID] - v[ee->face[0]->ID];
	    vec3dd pp = ev^an;
	    double w1 = f2*pp;
	    double w2 = f1*pp;

	    if (attracts_flow(ee->coface[0]->ID,ee->coface[0]->find_face_index(ee)) ^ 
		attracts_flow(ee->coface[1]->ID,ee->coface[1]->find_face_index(ee)))
	      eflow[i] = 0;
	    else
	      {
		assert(w1*w2<0);
		if ((fabs(w1)*f1+fabs(w2)*f2)*ev>0)
		  eflow[i] = 1;
		else
		  eflow[i] = 2;
	      }
	  }
	  break;

	default:
	  assert(0);
	  break;
	}
    }

  // fill proj
  proj0 = new double*[faces()];
  proj1 = new double*[faces()];
  for ( i=0; i<faces(); i++ )
    {
      int j;
      mesh_element *ff = getface(i);
      proj0[i] = new double[ff->faces>>1];
      proj1[i] = new double[ff->faces];
      vec3dd cmass(0,0,0);
      for ( j=1; j<ff->faces; j+=2 )
	cmass += v[ff->face[j]->ID];
      cmass *= (1.0/(ff->faces>>1));

      int check[2][2] = { {0,0}, {0,0} };
      for ( j=1; j<ff->faces; j+=2 )
	{
	  proj0[i][j>>1] = (f[i]^n[i])*(v[ff->face[j]->ID]-cmass);
	  if (j>1)
	    check[aflow[i][j-1] ? 0 : 1][proj0[i][j>>1]>=proj0[i][(j-2)>>1] ? 0 : 1] = 1;
	}
      check[aflow[i][0] ? 0 : 1][proj0[i][0]>=proj0[i][(j-2)>>1] ? 0 : 1] = 1;
      if (check[0][0]+check[1][0]!=1 || check[1][0]+check[1][1]!=1 || 
	  check[0][0]+check[1][0]!=1 || check[0][1]+check[1][1]!=1)
	{
	  cout << "Projection along vector field are inconsistent, ans so may be the output...." << endl;
	  cout << "Try applying a small random perturbation to the vector field." << endl;
	  cout << "If this does not work, check mesh for consistency." << endl;
	  assert(0);
	}

      for ( j=0; j<ff->faces; j+=2 )
	{
	  proj1[i][j] = (f[i]^n[i])*(v[ff->face[j]->face[0]->ID]-cmass);
	  proj1[i][j+1] = (f[i]^n[i])*(v[ff->face[j]->face[1]->ID]-cmass);
	}
    }

  // fill the vertex data now...
  
  indx = new int[vertices()];
  indx2 = new int[vertices()];
  isstat = new bool[vertices()];
  spiral = new bool[vertices()];

  for ( i=0; i<vertices(); i++ )
    {
      int j;
      mesh_element *cv = getvertex(i);

      bool boundary = (cv->cofaces&1);
      spiral[i] = false;

      // sector counts
      int sp = 0;
      int up = 0;
      int hy = 0;
      int el = 0;

      unsigned char *st = new unsigned char[cv->cofaces];

      for ( j=0; j<cv->cofaces; j++ )
	{
	  if (j&1)
	    st[j] = aflow[cv->coface[j]->ID][cv->coface[j]->find_face_index(cv->coface[j-1])] |
	      (aflow[cv->coface[j]->ID][cv->coface[j]->find_face_index(cv->coface[(j+1)%cv->cofaces])] << 1);
	  else
	    {
	      st[j] = (cv->coface[j]->face[0]==cv) ? eflow[cv->coface[j]->ID] : (eflow[cv->coface[j]->ID]^3);
	      if (st[j]==3) st[j]=0;
	    }
	}

      // st[i]: for even i, direction of flow along edge; 1: away, 2: toward the vertex
      // for odd i: whether flow is attracted in the respective face

      unsigned char lastdir = 0;   // 1: unstable, 2: stable
      unsigned char firstdir = 0;  // 1: unstable, 2; stable

      for ( j=0; j<cv->cofaces; j++ )
	{
	  if (!(j&1))
	    {
	      // incident edge

	      if (!st[j])
		continue;  // skip over if no flow
	      
	      if (!firstdir)
		{
		  lastdir = firstdir = st[j];
		}
	      else
		{
		  assert(j>0);
		  if (st[j]==lastdir)
		    continue;
		  if (lastdir==1)
		    {
		      up++;
		      if (st[j-1]>>1) el++;
		      else hy++;
		    }
		  if (lastdir==2)
		    {
		      sp++;
		      if (st[j-1]>>1) hy++;
		      else el++;
		    }
		  lastdir = st[j];
		}
	    }
	  else
	    {
	      // incident face

	      if (st[j]==1 || st[j]==2)
		continue;   // skip over if no (un)stable direction inside

	      unsigned char code = st[j] ? 2 : 1;

	      if (!firstdir)
		{
		  lastdir = firstdir = code;
		}
	      else
		{
		  if (code==lastdir)
		    continue;
		  if (lastdir==1)
		    {
		      up++;
		      hy++;
		    }
		  if (lastdir==2)
		    {
		      sp++;
		      hy++;
		    }
		  lastdir = code;
		}
	    }
	}

      if (!firstdir)
	{
	  // no (un)stable directions
	  // spiral sink/source
	  if (BD || !boundary)
	    {
	      indx[i] = indx2[i] = 1;
	      isstat[i] = true;
	      spiral[i] = true;
	      delete[] st;
	      continue;
	    }
	}
    
      if (boundary)
	{
	  if (BD)
	    {
	      if (lastdir==1) up++;
	      else sp++;

	      // add flow toward the boundary to get rid of it
	      //  i.e. reduce to internal vertex case
	      // adjust the sector counts 

	      int hy2 = hy;
	      int sp2 = sp;
	      int up2 = up;
	      int el2 = el;

	      bool force_stat = false;

	      switch(st[0] | (st[cv->cofaces-1]<<2))
		{
		case 5:
		  // both unstable...
		  hy+=2; sp++;
		  up2--;
		  if (!up2) up2 = 1;
		  force_stat = true;
		  break;
		case 6:
		case 9:
		  // one unstable, one stable
		  hy += 1;
		  hy2+= 1;
		  break;
		case 10:
		  // both stable
		  sp--;
		  if (!sp) sp = 1;
		  hy2+=2; up2++;
		  force_stat = true;
		  break;
		default:
		  assert(0);
		}

	      delete[] st;
	      indx[i] = 1+(el-hy)/2;
	      indx2[i] = 1+(el2-hy2)/2;
	      isstat[i] = (el || !(sp==1 && up==1 && hy==2) || el2 || !(sp2==1 && up2==1 && hy2==2) || force_stat);
	      continue;
	    }
	  else
	    {
	      isstat[i] = false;
	      indx[i] = indx2[i] = 0;
	      spiral[i] = false;
	      delete[] st;
	      continue;
	    }
	}
      else
	{
	  if (lastdir!=firstdir)
	    {
	      // need to add a hyperbolic or elliptic sector
	      if (firstdir==1)
		{
		  assert(lastdir==2);
		  sp++;
		  if (st[cv->cofaces-1]>>1)
		    hy++;
		  else
		    el++;
		}
	      else
		{
		  assert(lastdir==1);
		  assert(firstdir==2);
		  up++;
		  if (st[cv->cofaces-1]>>1)
		    el++;
		  else
		    hy++;
		}
	    }
	  else
	    {
	      if (sp==0 && up==0)
		if (lastdir==2)
		  sp++;
		else
		  up++;
	    }

	  indx2[i] = indx[i] = 1+(el-hy)/2;
	  isstat[i] = (el || !(sp==1 && up==1 && hy==2));
	  delete[] st;
	  continue;
	}
    }
}

/* ------------------------------------------------------ */

pcvf::~pcvf()
{
  if (f) delete[] f;
  f = NULL;
  if (pvf) delete[] pvf;
  pvf = NULL;
  if (isstat) delete[] isstat;
  isstat = NULL;
  if (indx) delete[] indx;
  indx = NULL;
  if (indx2) delete[] indx2;
  indx2 = NULL;
  if (spiral) delete[] spiral;
  spiral = NULL;
  if (eflow) delete[] eflow;
  eflow = NULL;
  for ( int i=0; i<faces(); i++ )
    {
      if (aflow[i]) delete[] aflow[i];
      aflow[i] = NULL;
      if (rflow[i]) delete[] rflow[i];
      rflow[i] = NULL;
      if (proj0[i]) delete[] proj0[i];
      proj0[i] = NULL;
      if (proj1[i]) delete[] proj1[i];
      proj1[i] = NULL;
    }
  if (aflow) delete[] aflow;
  aflow = NULL;
  if (rflow) delete[] rflow;
  rflow = NULL;
  if (proj0) delete[] proj0;
  proj0 = NULL;
  if (proj1) delete[] proj1;
  proj1 = NULL;
}

/* ------------------------------------------------------ */

bool pcvf::hasflowup ( int i )
{
  return eflow[i] & 1;
}

/* ------------------------------------------------------ */

bool pcvf::hasflowdown ( int i )
{
  return eflow[i] & 2;
}

/* ------------------------------------------------------ */

bool pcvf::attracts_flow ( int fce, int eix )
{
  return aflow[fce][eix];
}

/* ------------------------------------------------------ */

bool pcvf::repels_flow ( int fce, int eix )
{
  return rflow[fce][eix];
}

/* ------------------------------------------------------ */

bool pcvf::isstationary ( int i )
{
  return isstat[i];
}

/* ------------------------------------------------------ */

bool pcvf::isspiral ( int i )
{
  return spiral[i];
}

/* ------------------------------------------------------ */

int pcvf::index ( int i )
{
  return indx[i];
}

/* ------------------------------------------------------ */

vec3dd pcvf::getvector ( int i )
{
  return f[i];
}

/* ------------------------------------------------------ */

int pcvf::index2 ( int i )
{
  return indx2[i];
}

/* ------------------------------------------------------ */

bool pcvf::iststationary ( int i )
{
  return false; // technically, no triangle is in the PC case
}

/* ------------------------------------------------------ */

static void swap ( double &a, double &b )
{
  double tmp = a;
  a = b;
  b = tmp;
}

/* ------------------------------------------------------ */

bool pcvf::connects ( int fce, int eix1, int eix2, 
		      double s1, double e1, double s2, double e2 )
{
  if ( (eix1&1) && (eix2&1) )
    {
      // both are vertices
      return proj0[fce][eix1>>1]==proj0[fce][eix2>>1];
    }

  if ( (eix1&1) && !(eix2&1) )
    {
      // vertex and edge piece

      double a2,b2;
      a2 = proj1[fce][eix2];
      b2 = proj1[fce][eix2+1];
      double s = (1-s2)*a2+s2*b2;
      double e = (1-e2)*a2+e2*b2;
      if (s>e) swap(s,e);
      double p = proj0[fce][eix1>>1];
      return s<=p && p<=e;
    }

  if ( !(eix1&1) && (eix2&1) )
    {
      // vertex and edge piece

      double a1,b1;
      a1 = proj1[fce][eix1];
      b1 = proj1[fce][eix1+1];
      double s = (1-s1)*a1+s1*b1;
      double e = (1-e1)*a1+e1*b1;
      if (s>e) swap(s,e);
      double p = proj0[fce][eix2>>1];
      return s<=p && p<=e;
    }

  if ( !(eix1&1) && !(eix2&1) )
    {
      // both edges

      double a2,b2;
      a2 = proj1[fce][eix2];
      b2 = proj1[fce][eix2+1];
      double ss2 = (1-s2)*a2+s2*b2;
      double ee2 = (1-e2)*a2+e2*b2;
      if (ss2>ee2) swap(ss2,ee2);

      double a1,b1;
      a1 = proj1[fce][eix1];
      b1 = proj1[fce][eix1+1];
      double ss1 = (1-s1)*a1+s1*b1;
      double ee1 = (1-e1)*a1+e1*b1;
      if (ss1>ee1) swap(ss1,ee1);
      return ! (ee1<ss2 || ee2<ss1);
    }
}


/* ------------------------------------------------------ */

void pcvf::print_out()
{
  int i,j;
  mesh::print_out();
  cout << " ---------------------- " << endl;
  cout << "Vector field" << endl;
  for ( i=0; i<faces(); i++ )
    cout << i << " " << f[i] << endl;
  cout << " ---------------------- " << endl;
  cout << "Attracts/Repels flow" << endl;
  for ( i=0; i<faces(); i++ )
    {
      cout << i << " ";
      for ( j=0; j<getface(i)->faces; j++ )
	cout << (aflow[i][j] ? '1' : '0') << (rflow[i][j] ? '1' : '0') << " ";
      cout << endl;
    }
  cout << " ---------------------- " << endl;
  cout << "Edge flow" << endl;
  for ( i=0; i<edges(); i++ )
    cout << i << " " << (int)eflow[i] << endl;
  cout << " ---------------------- " << endl;
  cout << "Vertex data" << endl;
  for ( i=0; i<vertices(); i++ )
    cout << i << " " << indx[i] << " " << (isstat[i] ? '1' : '0') << " " << (spiral[i] ? 1 : 0) << endl;
}

/* ------------------------------------------------------ */
