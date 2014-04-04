
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

#include <mesh_base.h>
#include <cstdlib>
#include <cassert>
#include <iostream>

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

class i4 {
 public:
  int a,b,c,d;
  i4 ( int aa, int bb, int cc, int dd )
    {
      a = aa<bb ? aa : bb;
      b = aa<bb ? bb : aa;
      c = cc;
      d = dd;
    }
  i4() {}
};

class i3 {
 public:
  int a,b,c;
  i3 ( int aa, int bb, int cc )
    {
      a = aa;
      b = bb;
      c = cc;
    }
  i3() {}
};

int i4cmp ( const void *a0, const void *b0 )
{
  i4 *a = (i4*)a0;
  i4 *b = (i4*)b0;

  if (a->a<b->a)
    return -1;
  if (a->a>b->a)
    return 1;

  if (a->b<b->b)
    return -1;
  if (a->b>b->b)
    return 1;

  if (a->c<b->c)
    return -1;
  if (a->c>b->c)
    return 1;

  if (a->d<b->d)
    return -1;
  if (a->d>b->d)
    return 1;

  return 0;
}

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

mesh_element::mesh_element ( int id, int dim ) :
  ID(id), dimension(dim), faces(0), cofaces(0), face(NULL), coface(NULL)
{}

void mesh_element::set_faces ( int n, mesh_element **f )
{
  faces = n;
  face = f;
}

void mesh_element::set_cofaces ( int n, mesh_element **cf )
{
  cofaces = n;
  coface = cf;
}

bool mesh_element::isboundary()
{
  if (dimension==1)
    return cofaces==1;
  else
    if (dimension==0)
      return cofaces&1;
    else
      assert(0);
}

bool mesh_element::isisolatedvertex()
{
  return (dimension==0) && (cofaces==0);
}

mesh_element::~mesh_element()
{
  if (face)
    delete[] face;
  face = NULL;
  if (coface)
    delete[] coface;
  coface = NULL;
}

/* ------------------------------------------------------ */

mesh_base::mesh_base() : mel(), workspace(new std::vector<std::vector<int>*>)
{}

/* ------------------------------------------------------ */

mesh_element * mesh_base::get ( int i )
{
  assert(i>=0 && i<mel.size());
  return mel[i];
}

/* ------------------------------------------------------ */

mesh_element * mesh_base::getedge ( int i )
{
  assert(i>=0 && i<es);
  return mel[faces()+i];
}

/* ------------------------------------------------------ */

mesh_element * mesh_base::getvertex ( int i )
{
  assert(i>=0 && i<vtcs);
  return mel[faces()+es+i];
}

/* ------------------------------------------------------ */

mesh_element * mesh_base::getface ( int i )
{
  assert(i>=0 && i<fcs);
  return mel[i];
}

/* ------------------------------------------------------ */

void mesh_base::add_2Dmel ( std::vector<int> *verts )
{
  assert(workspace);
  workspace->push_back(verts);
}

/* ------------------------------------------------------ */

int mesh_element::find_face_index ( mesh_element *e )
{
  int i = -1;
  for ( i=0; i<faces; i++ )
    if (face[i]==e)
      break;
  assert(i>=0);
  return i;
}

/* ------------------------------------------------------ */

int mesh_base::edges()
{
  return es;
}

/* ------------------------------------------------------ */

int mesh_base::faces()
{
  return fcs;
}

/* ------------------------------------------------------ */

int mesh_base::vertices()
{
  return vtcs;
}

/* ------------------------------------------------------ */

int mesh_base::mesh_elements()
{
  return mel.size();
}

/* ------------------------------------------------------ */

void mesh_base::finalize()
{
  int i,j,k;
  int hedges = 0;
  int maxvid = -1; // max vertex id; number of vertices later on

  es = 0;
  fcs = workspace->size();

  // count half-edges

  for ( i=0; i<fcs; i++ )
    hedges += (*workspace)[i]->size();

  // need to construct half-edges and figure out their adjacency
  // first build an array of half-edges and sort it...
  // also, determine the max vertex id along the way
  
  i4 *e = new i4[hedges];
  k = 0;
  for ( i=0; i<fcs; i++ )
    {
      std::vector<int> *p = (*workspace)[i];
      int n = p->size();
      for ( j=0; j<p->size(); j++ )
	{
	  assert(k<hedges);
	  e[k++] = i4((*p)[j],(*p)[(j+1)%n],i,j);
	  if ((*p)[j]>maxvid)
	    maxvid = (*p)[j];
	}
    }
  assert(k==hedges);
  maxvid++;  // this is the number of vertices

  // sort the array...
  qsort(e,hedges,sizeof(i4),&i4cmp);

  // count edges
  for ( i=0; i<hedges; i++ )
    if (!i || e[i].a!=e[i-1].a || e[i].b!=e[i-1].b)
      es++;

  // allocate space for cells of different dimensions
  mesh_element **c2 = new mesh_element*[fcs];
  mesh_element **c1 = new mesh_element*[es];
  mesh_element **c0 = new mesh_element*[maxvid];
  for ( i=0; i<maxvid; i++ )
    {
      c0[i] = new mesh_element(i,0);
      c0[i]->set_faces(0,NULL);
      c0[i]->set_cofaces(0,NULL);
    }
  for ( i=0; i<es; i++ )
    c1[i] = new mesh_element(i,1);
  for ( i=0; i<fcs; i++ )
    {
      c2[i] = new mesh_element(i,2);
      c2[i]->set_cofaces(0,NULL);
      mesh_element **f = new mesh_element*[2*(*workspace)[i]->size()];
      for ( j=0; j<2*(*workspace)[i]->size(); j++ )
	f[j] = NULL;
      c2[i]->set_faces(2*(*workspace)[i]->size(),f);
    }

  // scan the sorted array and build edges/1D mesh elements
  for ( i=j=0; j<hedges; )
    {
      assert(i<es);
      if (j==0 || e[j].a!=e[j-1].a || e[j].b!=e[j-1].b)
	{
	  // new halfedge
	  mesh_element **f = new mesh_element*[2];
	  
	  f[0] = c0[e[j].a];
	  f[1] = c0[e[j].b];
	  c1[i]->set_faces(2,f);

	  assert(j+2>=hedges || !(e[j].a==e[j+2].a && e[j].b==e[j+2].b));
	  int cs = (e[j].a==e[j+1].a && e[j].b==e[j+1].b) ? 2 : 1;
	  mesh_element **c = new mesh_element*[cs];
	  c[0] = c2[e[j].c];
	  if (cs==2)
	    c[1] = c2[e[j+1].c];
	  c1[i]->set_cofaces(cs,c);

	  assert(!c2[e[j].c]->face[2*e[j].d]);
	  c2[e[j].c]->face[2*e[j].d] = c1[i];
	  if (cs==2)
	    {
	      assert(!c2[e[j+1].c]->face[2*e[j+1].d]);
	      c2[e[j+1].c]->face[2*e[j+1].d] = c1[i];
	    }
	  j += cs;
	  i++;
	}
    }

  assert(j==hedges);
  assert(i==es);

  // fill the 2D mesh element records
  for ( i=0; i<fcs; i++ )
    {
      int n = (*workspace)[i]->size();
      for ( j=0; j<n; j++ )
	{
	  mesh_element *m1 = c2[i]->face[2*j];
	  mesh_element *m2 = c2[i]->face[2*((j+1)%n)];
	  assert(m1); assert(m2);
	  assert(m1->dimension==1); assert(m2->dimension==1);
	  assert(m1->faces==2); assert(m2->faces==2);
	  int a1 = m1->face[0]->ID;
	  int b1 = m1->face[1]->ID;
	  int a2 = m2->face[0]->ID;
	  int b2 = m2->face[1]->ID;
	  if (a1==a2 || a1==b2)
	    c2[i]->face[2*j+1] = c0[a1];
	  else
	    if (b1==a2 || b1==b2)
	      c2[i]->face[2*j+1] = c0[b1];
	    else
	      assert(0);
	}
    }

  // all that remains is to fill the coface info for vertices
  // first, determine one incident edge for each vertex;
  // give preference to boundary edges
  std::vector<i3> *sec = new std::vector<i3>[maxvid];
  for ( i=0; i<fcs; i++ )
    {
      mesh_element *cc = c2[i];
      int n = cc->faces/2;
      for ( j=0; j<n; j++ )
	{
	  // std::cout << "L " << cc->face[2*j+1]->ID << " : " 
	  // << cc->face[2*j]->ID << " " << i << " " << cc->face[2*((j+2)%n)]->ID << std::endl;
	  sec[cc->face[2*j+1]->ID].push_back(i3(cc->face[2*j]->ID,i,cc->face[2*((j+1)%n)]->ID));
	}
    }

  int *pred = new int[es];
  int *succ = new int[es];
  int *fid = new int[es];

  for ( i=0; i<es; i++ )
    fid[i] = pred[i] = succ[i] = -1;

  for ( i=0; i<maxvid; i++ )
    {
      int one = -1;
      int firstone;

      //      std::cout << "Vertex " << i << std::endl;

      for ( j=0; j<sec[i].size(); j++ )
	{
	  //	  std::cout << sec[i][j].a << "--> " << sec[i][j].c << std::endl;
	  assert(pred[sec[i][j].c]==-1);
	  assert(succ[sec[i][j].a]==-1);
	  pred[sec[i][j].c] = sec[i][j].a;
	  one = succ[sec[i][j].a] = sec[i][j].c;
	  fid[sec[i][j].a] = sec[i][j].b;
	}
      
      assert(one>=0);

      firstone = one;

      // follow succ

      while(pred[one]>=0)
	{ 
	  one = pred[one];
	  if (one==firstone)
	    break;
	}

      firstone = one;

      if (pred[one]>=0)
	{
	  // internal vertex
	  mesh_element **c = new mesh_element*[2*sec[i].size()];

	  for ( j=0; j<sec[i].size(); j++ )
	    {
	      c[2*j] = c1[one];
	      c[2*j+1] = c2[fid[one]];
	      one = succ[one];
	    }
	  assert(one==firstone);

	  c0[i]->set_cofaces(2*sec[i].size(),c);
	}
      else
	{
	  // boundary vertex
	  mesh_element **c = new mesh_element*[2*sec[i].size()+1];

	  for ( j=0; j<sec[i].size(); j++ )
	    {
	      c[2*j] = c1[one];
	      c[2*j+1] = c2[fid[one]];
	      one = succ[one];
	    }
	  c[2*j] = c1[one];

	  c0[i]->set_cofaces(2*sec[i].size()+1,c);	  
	}

      for ( j=0; j<sec[i].size(); j++ )
	{
	  pred[sec[i][j].c] = -1;
	  succ[sec[i][j].a] = -1;
	}
    }

  // build list of all faces...
  for ( i=0; i<fcs; i++ )
    mel.push_back(c2[i]);
  for ( i=0; i<es; i++ )
    mel.push_back(c1[i]);
  for ( i=0; i<maxvid; i++ )
    mel.push_back(c0[i]);

  vtcs = maxvid;

  // clear workspace etc
  delete[] e;
  delete[] sec;
  delete[] fid;
  delete[] succ;
  delete[] pred;
  delete[] c0;
  delete[] c1;
  delete[] c2;
  for ( i=0; i<fcs; i++ )
    delete((*workspace)[i]);
  delete workspace;
  workspace = NULL;
}

/* ------------------------------------------------------ */

void mesh_element::print_out()
{
  int i;
  std::cout << "[[" << ID << " " << dimension << "]]  " << std::endl;
  std::cout << "     F(";
  for ( i=0; i<faces; i++ )
    std::cout << "[" << face[i]->ID << " " << face[i]->dimension << "]";
  std::cout << ")" << std::endl;
  std::cout << "     CF(";
  for ( i=0; i<cofaces; i++ )
    std::cout << "[" << coface[i]->ID << " " << coface[i]->dimension << "]";
  std::cout << ")" << std::endl;
  std::cout << " -------------------------------------------- " << std::endl;
}

void mesh_base::print_out()
{
  for ( int i=0; i<mel.size(); i++ )
    mel[i]->print_out();
}

/* ------------------------------------------------------ */

mesh_base::~mesh_base()
{
  if (workspace) delete workspace;
  for ( int i=0; i<mel.size(); i++ )
    {
      if (mel[i])
	{
	  delete mel[i];
	  mel[i] = NULL;
	}
    }
}

/* ------------------------------------------------------ */
