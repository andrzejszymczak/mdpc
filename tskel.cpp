

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

#include <tskel.h>
#include <fstream>
#include <cassert>

using namespace std;

/* ------------------------------------------------------------------ */

tsknode::tsknode ( ) : t(), in(), out() {}

/* ------------------------------------------------------------------ */

tsknode::tsknode ( int id, mstype bb ) : t(bb), in(), out() {}

/* ------------------------------------------------------------------ */

void tsknode::removein ( int a )
{
  for ( int i=0; i<in.size(); i++ )
    if (in[i]==a)
      {
	in[i] = in[in.size()-1];
	in.pop_back();
      }
}

/* ------------------------------------------------------------------ */

void tsknode::removeout ( int a )
{
  for ( int i=0; i<out.size(); i++ )
    if (out[i]==a)
      {
	out[i] = out[out.size()-1];
	out.pop_back();
      }
}

/* ------------------------------------------------------------------ */

bool tsknode::hasedgeto ( int i )
{
  for ( int j=0; j<out.size(); j++ )
    if (out[j]==i) return true;
  return false;
}

/* ------------------------------------------------------------------ */

bool tsknode::hasedgefrom ( int i )
{
  for ( int j=0; j<in.size(); j++ )
    if (in[j]==i) return true;
  return false;
}

/* ------------------------------------------------------------------ */

bool tsknode::isattracting()
{
  return out.size()==0;
}

/* ------------------------------------------------------------------ */

bool tsknode::isrepelling()
{
  return in.size()==0;
}

/* ------------------------------------------------------------------ */

tskel::tskel ( int nn, mstype *B ) : ns(nn)
{
  n = new tsknode[ns];
  for ( int i=0; i<ns; i++ )
    n[i] = tsknode(i,B[i]);
}

/* ------------------------------------------------------------------ */

void tskel::add_edge ( int i, int j )
{
  if (hasedge(i,j) || i==j) return;
  n[i].out.push_back(j);
  n[j].in.push_back(i);
}

/* ------------------------------------------------------------------ */

void tskel::remove_edge ( int i, int j )
{
  n[i].removeout(j);
  n[j].removein(i);
}

/* ------------------------------------------------------------------ */

void tskel::save ( const char *name )
{
  int i,j;
  ofstream ofs(name);

  ofs << "digraph G {" << endl;
  ofs << "node [ fontsize = 40, color = black, shape = box]" << endl;
  for ( i=0; i<ns; i++ )
    {
      if (!n[i].in.size() && !n[i].out.size())
	continue;
      ofs << i << " [ color = ";
      mstype b = n[i].t;
      if (b.issaddle())
	ofs << "blue";
      else
	if (b.issink())
	  ofs << "green";
	else
	  if (b.isapo())
	    ofs << "green, style = rounded";
	  else
	    if (b.issource())
	      ofs << "red";
	    else
	      if (b.isrpo())
		ofs << "red, style = rounded";
	      else
		if (b.istrivial())
		  ofs << "magenta";
		else
		  ofs << "grey";

      ofs << " ]" << endl;
    }
  for ( i=0; i<ns; i++ )
    for ( j=0; j<n[i].out.size(); j++ )
      {
	ofs << i << "->" << n[i].out[j];
	if (!iscertain(i,n[i].out[j]))
	  ofs << " [style=dotted]";
	ofs << endl;
      }
  ofs << "}" << endl;
}

/* ------------------------------------------------------------------ */

void tskel::cleanup()
{
  vector<int> a,b;

  for ( int i=0; i<ns; i++ )
    for ( int j=0; j<n[i].out.size(); j++ )
      for ( int k=0; k<n[n[i].out[j]].out.size(); k++ )
	{
	  a.push_back(i);
	  b.push_back(n[n[i].out[j]].out[k]);
	}

  for ( int l=0; l<a.size(); l++ )
    remove_edge(a[l],b[l]);
}

/* ------------------------------------------------------------------ */

bool tskel::hasedge ( int i, int j )
{
  assert(n[i].hasedgeto(j)==n[j].hasedgefrom(i));
  return n[i].hasedgeto(j);
}

/* ------------------------------------------------------------------ */

bool tskel::iscertain ( int i, int j )
{
  return !n[i].t.istrivial() && !n[j].t.istrivial() && (n[j].isattracting() || n[i].isrepelling());
}

/* ------------------------------------------------------------------ */

tskel::~tskel()
{
  if (n) delete[] n;
  n = NULL;
}

/* ------------------------------------------------------------------ */
