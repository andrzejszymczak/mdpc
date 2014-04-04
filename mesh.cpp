

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

#include <mesh.h>
#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

/* ------------------------------------------------------ */

mesh::mesh ( const char *name ) : mesh_base(), ifs(name)
{
  if (!ifs) 
    {
      cout << "Can't open " << name << endl;
      assert(ifs);
    }

  // detect format
  char c = ifs.get();
  char format;

  if (c=='n')
    format = 'n';
  else
    {
      ifs.putback(c);
      format = 't';
    }

  read_from_file(ifs,format);

  compute_normals();
}

/* ------------------------------------------------------ */

void mesh::read_from_file ( ifstream &ifs, char format )
{
  switch(format)
    {
    case 't':
      {
	int vs,ts,i;
	ifs >> ts >> vs;
	for ( i=0; i<ts; i++ )
	  {
	    int a,b,c;
	    ifs >> a >> b >> c;
	    vector<int> *t = new vector<int>;
	    t->push_back(a);
	    t->push_back(b);
	    t->push_back(c);
	    add_2Dmel(t);
	  }
	finalize();

	if (vs!=vtcs)
	  {
	    cout << "Warning: nused vertices detected" << endl;
	  }

	v = new vec3dd[vtcs];
	for ( i=0; i<vtcs; i++ )
	  ifs >> v[i][0] >> v[i][1] >> v[i][2];
      }
      break;

    case 'n':
      {
	int vs,fs,i;
	ifs >> fs >> vs;
	for ( i=0; i<fs; i++ )
	  {
	    int a;
	    vector<int> *t = new vector<int>;
	    while(1)
	      {
		ifs >> a;
		if (a==-1)
		  break;
		t->push_back(a);
	      }
	    add_2Dmel(t);
	  }
	finalize();

	if (vs!=vtcs)
	  {
	    cout << "Warning: nused vertices detected" << endl;
	  }

	v = new vec3dd[vtcs];
	for ( i=0; i<vtcs; i++ )
	  ifs >> v[i][0] >> v[i][1] >> v[i][2];
      }
    }
}

/* ------------------------------------------------------ */

vec3dd mesh::evec ( const mesh_element *f, int i, int j )
{
  return v[f->face[j]->ID]-v[f->face[i]->ID];
}

/* ------------------------------------------------------ */

vec3dd mesh::edgepoint ( int eid, double t )
{
  mesh_element *e = getedge(eid);
  return (1-t)*v[e->face[0]->ID]+t*v[e->face[1]->ID];
}

/* ------------------------------------------------------ */

vec3dd mesh::vertex ( int i )
{
  return v[i];
}

/* ------------------------------------------------------ */

vec3dd mesh::normal ( int i )
{
  return n[i];
}

/* ------------------------------------------------------ */

void mesh::compute_normals()
{
  n = new vec3dd[fcs];
  for ( int i=0; i<fcs; i++ )
    {
      vec3dd n0(0,0,0);
      const mesh_element *f = get(i);
      for ( int j=1; j<f->faces-4; j+=2 )
	n0 += evec(f,j,j+2)^evec(f,j,j+4);
      n0.normalize();
      n[i] = n0;
    }
}

/* ------------------------------------------------------ */

void mesh::print_out()
{
  int i;
  mesh_base::print_out();
  cout << " ---------------------- " << endl;
  cout << "Vertex coordinates" << endl;
  for ( i=0; i<vtcs; i++ )
    cout << i << " " << v[i] << endl;
  cout << " ---------------------- " << endl;
  cout << "Normals" << endl;
  for ( i=0; i<fcs; i++ )
    cout << i << " " << n[i] << endl;
}

/* ------------------------------------------------------ */

mesh::~mesh()
{
  if (v) delete[] v;
  if (n) delete[] n;
  v = n = NULL;
}

/* ------------------------------------------------------ */
