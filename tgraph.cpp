

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

#include <tgraph.h>
#include <iostream>
#include <primitive.h>

#define SIX0 2
#define SIX1 6
#define SIXC 10

using namespace std;

/* ------------------------------------------------------ */

int arc::arcs = 0;
int node::nodes = 0;

/* ------------------------------------------------------ */

node::node ( int id, mesh_element *o ) : 
  owner(o), in(), out(), left(NULL), right(NULL), s(0), e(0), scc(-1), flags(0), ID(id)
{
  nodes++;
}

/* ------------------------------------------------------ */

node::node ( int id, mesh_element *o, node *l, node *r, double ss, double ee ) :
  owner(o), in(), out(), left(l), right(r), s(ss), e(ee), scc(-1), flags(0), ID(id)
{
  nodes++;
}

/* ------------------------------------------------------ */

node::~node()
{
  while(in.size())
    delete in[0];
  while(out.size())
    delete out[0];
  if (left)
    {
      assert(left->right==this);
      left->right = NULL;
    }
  if (right)
    {
      assert(right->left==this);
      right->left = NULL;
    }
  nodes--;
}

/* ------------------------------------------------------ */

bool node::isepiece()
{
  return owner->dimension==1;
}

/* ------------------------------------------------------ */

void node::lockL()
{
  flags |= 32;
  if (left)
    left->flags |= 16;
}

/* ------------------------------------------------------ */

// this is an untested fix: originally left instead of right

void node::lockR()
{
  flags |= 16;
  if (right)
    right->flags |= 32;
}

/* ------------------------------------------------------ */

bool node::islockedR()
{
  return flags & 16;
}

/* ------------------------------------------------------ */

bool node::islockedL()
{
  return flags & 32;
}

/* ------------------------------------------------------ */

void node::markf ( unsigned char msk, vector<node*> *lst )
{
  stack<node*> S;

  S.push(this);
  while (!S.empty())
    {
      node *c = S.top();
      S.pop();
      unsigned char ofl = c->flags;
      c->flags |= msk;
      if (ofl==c->flags)
	continue;
      if (lst) lst->push_back(c);
      for ( int i=0; i<c->out.size(); i++ )
	S.push(c->out[i]->to);
    }
}

/* ------------------------------------------------------ */

void node::markb ( unsigned char msk, vector<node*> *lst )
{
  stack<node*> S;

  S.push(this);
  while (!S.empty())
    {
      node *c = S.top();
      S.pop();
      unsigned char ofl = c->flags;
      c->flags |= msk;
      if (ofl==c->flags)
	continue;
      if (lst) lst->push_back(c);
      for ( int i=0; i<c->in.size(); i++ )
	S.push(c->in[i]->from);
    }
}

/* ------------------------------------------------------ */

double node::span()
{
  return e-s;
}

/* ------------------------------------------------------ */

arc::arc ( node *f, node *t, mesh_element *carrier, int ix_orig, int ix_dest ) :
  from(f), to(t)
{
  assert(carrier->dimension==2);
  assert(ix_orig<16);
  assert(ix_dest<16);

  data = 2 | (ix_orig << SIX0) | (ix_dest << SIX1 ) | (carrier->ID << SIXC);

  t->in.push_back(this);
  f->out.push_back(this);

  arcs++;
}

/* ------------------------------------------------------ */

arc::arc ( node *f, node *t, mesh_element *carrier ) :
  from(f), to(t)
{
  switch(carrier->dimension)
    {
    case 0:
      data = (carrier->ID << SIXC);
      break;
    case 1:
      data = 1 | (carrier->ID << SIXC);
      break;
    case 2:
      {
	int ix_orig = carrier->find_face_index(f->owner);
	int ix_dest = carrier->find_face_index(t->owner);
	data = 2 | (ix_orig << SIX0) | (ix_dest << SIX1 ) | (carrier->ID << SIXC);
      }
    default:
      assert(0);
    }
  t->in.push_back(this);
  f->out.push_back(this);

  arcs++;
}

/* ------------------------------------------------------ */

arc::arc ( node *f, node *t, arc *a ) :
  from(f), to(t), data(a->data)
{
  t->in.push_back(this);
  f->out.push_back(this);

  arcs++;
}

/* ------------------------------------------------------ */

arc::~arc()
{
  int i;

  // remove arc from starting node
  for ( i=0; i<from->out.size(); i++ )
    if (this==from->out[i])
      break;
  assert(i<from->out.size());
  from->out[i] = from->out[from->out.size()-1];
  from->out.pop_back();

  // remove arc from end node
  for ( i=0; i<to->in.size(); i++ )
    if (this==to->in[i])
      break;
  assert(i<to->in.size());
  to->in[i] = to->in[to->in.size()-1];
  to->in.pop_back();

  arcs--;
}

/* ------------------------------------------------------ */

int arc::get_dimension()
{
  return data & ((1<<SIX0)-1);
}

int arc::get_ix()
{
  return data>>SIXC;
}

int arc::get_ixf()
{
  return (data>>SIX0) & ((1<<(SIX1-SIX0))-1);
}

int arc::get_ixt()
{
  return (data>>SIX1) & ((1<<(SIXC-SIX1))-1);
}

/* ------------------------------------------------------ */

int tgraph::_index ( mesh_element *e )
{
  assert(e->dimension==0 || e->dimension==1);
  return e->dimension ? (e->ID) : (e->ID+msh->edges());
}

node *tgraph::_getnode ( mesh_element *e )
{
  return n[_index(e)];
}

node *tgraph::_getedgenode ( int i )
{
  return n[i];
}

node *tgraph::_getvertexnode ( int i )
{
  return n[i+msh->edges()];
}

/* ------------------------------------------------------ */

void tgraph::remove_node ( int i )
{
  if (n[i]) 
    {
      delete n[i];
      n[i] = NULL;
    }
  if (i!=n.size()-1)
    {
      n[i] = n[n.size()-1];
      n[i]->ID = i;
    }
  n.pop_back();
}

/* ------------------------------------------------------ */

void tgraph::remove_all_nonSCC()
{
  int i,j;

  for ( i=0; i<n.size(); )
    if (n[i] && n[i]->scc==-1)
      {
	node *nn = n[i];
	bool remove = true;
	for ( j=0; j<nn->in.size(); j++ )
	  if (nn->in[j]->from->scc!=-1)
	    {
	      remove = false;
	      break;
	    }
	if (!remove)
	  {
	    i++;
	    continue;
	  }
	for ( j=0; j<nn->out.size(); j++ )
	  if (nn->out[j]->to->scc!=-1)
	    {
	      remove = false;
	      break;
	    }
	if (!remove) 
	  {
	    i++;
	    continue;
	  }
	remove_node(i);
      }
    else
      i++;
}

/* ------------------------------------------------------ */

tgraph::tgraph ( vfield_base *m ) : msh(m), mstp(NULL), sccs(-1), n()
{
  int i;

  // add edge pieces first
  for ( i=0; i<msh->mesh_elements(); i++ )
    {
      mesh_element *mel = msh->get(i);
      switch(mel->dimension)
	{
	case 2:
	  break;
	case 1:
	  n.push_back(new node(n.size(),mel,NULL,NULL,0,1));
	  break;
	case 0:
	  n.push_back(new node(n.size(),mel));
	  break;
	default:
	  assert(0);
	}
    }

  // now, construct arcs...

  for ( i=0; i<msh->vertices(); i++ )
    {
      if (msh->isstationary(i))
	new arc(_getvertexnode(i),_getvertexnode(i),msh->getvertex(i));
    }

  for ( i=0; i<msh->edges(); i++ )
    {
      if (msh->hasflowup(i))
	{
	  new arc(_getedgenode(i),_getvertexnode(msh->getedge(i)->face[1]->ID),msh->getedge(i));
	  new arc(_getvertexnode(msh->getedge(i)->face[0]->ID),_getedgenode(i),msh->getedge(i));
	}
      if (msh->hasflowdown(i))
	{
	  new arc(_getedgenode(i),_getvertexnode(msh->getedge(i)->face[0]->ID),msh->getedge(i));
	  new arc(_getvertexnode(msh->getedge(i)->face[1]->ID),_getedgenode(i),msh->getedge(i));
	}
    }

  for ( i=0; i<msh->faces(); i++ )
    {
      mesh_element *cf = msh->getface(i);
      for ( int j=0; j<cf->faces; j++ )
	for ( int k=0; k<cf->faces; k++ )
	  {
	    // originally, only k==j was used (wrong!)
	    if ((k==j) || (k==(j+1)%cf->faces) || (k==(j+cf->faces-1)%cf->faces)) continue;
	    if (msh->attracts_flow(i,k) && msh->repels_flow(i,j))
	      if (msh->connects(i,j,k))
		new arc(_getnode(cf->face[j]),_getnode(cf->face[k]),cf,j,k);
	  }
    }

  // we also need to add arcs from any isolated vertex to its incident edge and back

  for ( i=0; i<msh->vertices(); i++ )
    {
      if (!msh->isspiral(i))
	continue;
      node *vn = _getvertexnode(i);
      node *ince = _getedgenode(vn->owner->coface[0]->ID);
      new arc(ince,vn,ince->owner);
      new arc(vn,ince,ince->owner);
    }
}

/* ------------------------------------------------------ */

node * tgraph::merge ( int i, char side, bool move )
{
  int j;
 // for the new node
  float s,e; 
  node *l,*r; 

  if (side=='l' || side=='L')
    {
      assert(!n[i]->islockedL());
      assert(n[i]->left);
      j = n[i]->left->ID;
      s = n[i]->left->s;
      e = n[i]->e;
      l = n[j]->left;
      r = n[i]->right;
    }
  else
    if (side=='r' || side=='R')
      {
	assert(!n[i]->islockedR());
	assert(n[i]->right);
	j = n[i]->right->ID;
	s = n[i]->s;
	e = n[i]->right->e;
	l = n[i]->left;
	r = n[j]->right;
      }
    else
      assert(0);

  // need to merge nodes i and j

  node *nn = new node(i,n[i]->owner,l,r,s,e);

  if (side=='l' || side=='L')
    {
      if (n[i]->islockedR())
	nn->lockR();
      if (n[j]->islockedL())
	nn->lockL();
    }

  if (side=='r' || side=='R')
    {
      if (n[i]->islockedL())
	nn->lockL();
      if (n[j]->islockedR())
	nn->lockR();      
    }

  n[i]->flags |= (64|128);
  n[j]->flags |= (64|128);

  int k;
  for ( k=0; k<n[i]->out.size(); k++ )
    {
      if (!(n[i]->out[k]->to->flags & 64))
	{
	  new arc(nn,n[i]->out[k]->to,n[i]->out[k]);
	  n[i]->out[k]->to->flags |= 64;
	}
    }
  for ( k=0; k<n[j]->out.size(); k++ )
    {
      if (!(n[j]->out[k]->to->flags & 64))
	{
	  new arc(nn,n[j]->out[k]->to,n[j]->out[k]);
	  n[j]->out[k]->to->flags |= 64;
	}
    }
  for ( k=0; k<n[i]->in.size(); k++ )
    {
      if (!(n[i]->in[k]->from->flags & 128))
	{
	  new arc(n[i]->in[k]->from,nn,n[i]->in[k]);
	  n[i]->in[k]->from->flags |= 128;
	}
    }
  for ( k=0; k<n[j]->in.size(); k++ )
    {
      if (!(n[j]->in[k]->from->flags & 128))
	{
	  new arc(n[j]->in[k]->from,nn,n[j]->in[k]);
	  n[j]->in[k]->from->flags |= 128;
	}
    }

  // clear flags now...
  n[i]->flags &= ~(64|128);
  n[j]->flags &= ~(64|128);
  for ( k=0; k<nn->out.size(); k++ )
    nn->out[k]->to->flags &= ~(64|128);
  for ( k=0; k<nn->in.size(); k++ )
    nn->in[k]->from->flags &= ~(64|128);

  // set attributes of the new node
  nn->scc = n[i]->scc;
  nn->flags |= n[i]->flags & (~(16|32|64|128));

  // delete old nodes
  n[i]->right = n[j]->right = n[i]->left = n[j]->left = NULL;
  delete n[i];
  delete n[j];
  n[i] = n[j] = NULL;

  // update n
  n[i] = nn;
  if (move)
    remove_node(j);

  // update left and right
  if (n[i]->left)
    n[i]->left->right = n[i];
  if (n[i]->right)
    n[i]->right->left = n[i];

  return nn;
}

/* ------------------------------------------------------ */

void tgraph::subdivide ( int i, double spt )
{
  if (!(n[i] && n[i]->isepiece()))
    return;
  double mid = n[i]->s+spt*(n[i]->e-n[i]->s);
  node *nl = new node(i,n[i]->owner,n[i]->left,NULL,n[i]->s,mid);
  node *nr = new node(n.size(),n[i]->owner,nl,n[i]->right,mid,n[i]->e);
  nl->right = nr;
  if (n[i]->islockedL()) 
    nl->lockL();
  if (n[i]->islockedR())
    nr->lockR();

  // construct arcs now...
  if (msh->hasflowup(n[i]->owner->ID))
    new arc(nl,nr,n[i]->owner);
  if (msh->hasflowdown(n[i]->owner->ID))
    new arc(nr,nl,n[i]->owner);

  int j;
  for ( j=0; j<n[i]->in.size(); j++ )
    {
      arc *a = n[i]->in[j];
      node *src = a->from;
      switch(a->get_dimension())
	{
	case 1:
	  {
	    // connect only if starting and end nodes intersect
	    if (src->owner->dimension==0)
	      {
		if (src->owner==n[i]->owner->face[0])
		  {
		    assert(nl->s==0);
		    new arc(src,nl,n[i]->owner);
		  }
		else
		  {
		    assert(src->owner==n[i]->owner->face[1]);
		    assert(nr->e==1);
		    new arc(src,nr,n[i]->owner);
		  }
	      }
	    else
	      {
		assert(src->owner->dimension==1);
		if (src->e==nl->s || src->s==nl->e)
		  new arc(src,nl,n[i]->owner);
		else
		  if (src->e==nr->s || src->s==nr->e)
		    new arc(src,nr,n[i]->owner);
		  else
		    assert(0);
	      }
	  }
	  break;
	case 2:
	  {
	    bool yes = false;
	    if (msh->connects(a->get_ix(),a->get_ixf(),a->get_ixt(),src->s,src->e,nl->s,nl->e))
	      {
		new arc(src,nl,msh->getface(a->get_ix()),a->get_ixf(),a->get_ixt());
		yes = true;
	      }
	    if (msh->connects(a->get_ix(),a->get_ixf(),a->get_ixt(),src->s,src->e,nr->s,nr->e))
	      {
		new arc(src,nr,msh->getface(a->get_ix()),a->get_ixf(),a->get_ixt());
		yes = true;
	      }
	    //	    assert(yes);
	  }
	  break;
	default:
	  assert(0);
	  break;
	}
    }

  for ( j=0; j<n[i]->out.size(); j++ )
    {
      arc *a = n[i]->out[j];
      node *dst = a->to;
      switch(a->get_dimension())
	{
	case 1:
	  {
	    // connect only if starting and end nodes intersect
	    if (dst->owner->dimension==0)
	      {
		if (dst->owner==n[i]->owner->face[0])
		  {
		    assert(nl->s==0);
		    new arc(nl,dst,n[i]->owner);
		  }
		else
		  {
		    assert(dst->owner==n[i]->owner->face[1]);
		    assert(nr->e==1);
		    new arc(nr,dst,n[i]->owner);
		  }
	      }
	    else
	      {
		assert(dst->owner->dimension==1);
		if (dst->e==nl->s || dst->s==nl->e)
		  new arc(nl,dst,n[i]->owner);
		else
		  if (dst->e==nr->s || dst->s==nr->e)
		    new arc(nr,dst,n[i]->owner);
		  else
		    assert(0);
	      }
	  }
	  break;
	case 2:
	  {
	    bool yes = false;
	    if (msh->connects(a->get_ix(),a->get_ixf(),a->get_ixt(),nl->s,nl->e,dst->s,dst->e))
	      {
		new arc(nl,dst,msh->getface(a->get_ix()),a->get_ixf(),a->get_ixt());
		yes = true;
	      }
	    if (msh->connects(a->get_ix(),a->get_ixf(),a->get_ixt(),nr->s,nr->e,dst->s,dst->e))
	      {
		new arc(nr,dst,msh->getface(a->get_ix()),a->get_ixf(),a->get_ixt());
		yes = true;
	      }
	    //	    assert(yes);
	  }
	  break;
	default:
	  assert(0);
	  break;
	}
    }

  delete n[i];
  n[i] = nl;
  n.push_back(nr);

  if (nl->left)
    nl->left->right = nl;
  if (nl->right)
    nl->right->left = nl;
  if (nr->left)
    nr->left->right = nr;
  if (nr->right)
    nr->right->left = nr;
}

/* ------------------------------------------------------ */

void tgraph::subdivide_all()
{
  int num = n.size();

  for ( int i=0; i<num; i++ )
    subdivide(i);
}

/* ------------------------------------------------------ */

void tgraph::subdivide_scc_nodes()
{
  int num = n.size();
  for ( int i=0; i<num; i++ )
    if (n[i] && n[i]->scc>=0)
      subdivide(i);
}

/* ------------------------------------------------------ */

int tgraph::nodes()
{
  return node::nodes;
}

/* ------------------------------------------------------ */

int tgraph::arcs()
{
  return arc::arcs;
}

/* ------------------------------------------------------ */

int tgraph::arcs2()
{
  int res = 0;
  for ( int i=0; i<n.size(); i++ )
    if (n[i])
      for ( int j=0; j<n[i]->out.size(); j++ )
	if (n[i]->out[j]->get_dimension()==2) res++;
  return res;
}

/* ------------------------------------------------------ */

void tgraph::computeSCCs()
{
  int i;
  for ( i=0; i<n.size(); i++ )
    if (n[i]) n[i]->scc = 0;

  __S = new stack<node*>;
  __index = 1; 
  int nds = nodes();
  __c = nds-1;

  for ( i=0; i<n.size(); i++ )
    if (n[i] && !n[i]->scc)
      __visit(n[i]);

  delete __S;
  __S = NULL;

  // now, cleanup! Need to remove size-1 sccs except if represent a atationary vertex
  vector<node*> *nv = new vector<node*>[nds];
  for ( i=0; i<n.size(); i++ )
    if (n[i])
      {
	assert(n[i]->scc>=0 && n[i]->scc<nds);
	nv[n[i]->scc].push_back(n[i]);
      }
  
  int curid = 0;
  for ( i=0; i<nds; i++ )
    if (nv[i].size()>1 || (nv[i].size()==1 && nv[i][0]->owner->dimension==0 && msh->isstationary(nv[i][0]->owner->ID)))
      {
	for ( int j=0; j<nv[i].size(); j++ )
	  nv[i][j]->scc = curid;
	curid++;
      }
    else
      for ( int j=0; j<nv[i].size(); j++ )
	nv[i][j]->scc = -1;

  delete[] nv;
  sccs = curid;
}

/* ------------------------------------------------------ */

void tgraph::__visit ( node *v )
{
  bool root = true;
  v->scc = __index;
  __index++;

  for ( int i=0; i<v->out.size(); i++ )
    {
      node *w = v->out[i]->to;
      if (!w->scc) 
	__visit(w);
      if (w->scc<v->scc)
	{
	  v->scc=w->scc;
	  root = false;
	}
    }
  if (root)
    {
      __index--;
      while (!__S->empty() && v->scc<=__S->top()->scc)
	{
	  node *w = __S->top();
	  __S->pop();
	  w->scc = __c;
	  __index--;
	}
      v->scc = __c;
      __c--;
    }
  else
    __S->push(v);
}

/* ------------------------------------------------------ */

int tgraph::SCCs()
{
  return sccs;
}

/* ------------------------------------------------------ */

void tgraph::computeMSTypes()
{
  int i,j;

  if (mstp) delete[] mstp;

  if (sccs==0)
    {
      mstp = NULL;
      return;
    }

  mstp = new mstype[sccs];
  for ( i=0; i<n.size(); i++ )
    if (n[i] && n[i]->scc>=0)
      {
	mstp[n[i]->scc].incsize();

	if (n[i]->owner->dimension==0)
	  {
	    mstp[n[i]->scc].addtoindex(msh->index(n[i]->owner->ID));
	    mstp[n[i]->scc].addtoindex2(msh->index2(n[i]->owner->ID));
	    if (n[i]->owner->isboundary())
	      mstp[n[i]->scc].setbdry();
	  }

	for ( j=0; j<n[i]->in.size(); j++ )
	  if (n[i]->in[j]->from->scc!=n[i]->scc)
	    mstp[n[i]->scc].orstability(1);
	for ( j=0; j<n[i]->out.size(); j++ )
	  if (n[i]->out[j]->to->scc!=n[i]->scc)
	    mstp[n[i]->scc].orstability(2);
      }
}

/* ------------------------------------------------------ */

mstype tgraph::MStype ( int i )
{
  return mstp[i];
}

/* ------------------------------------------------------ */

void tgraph::saveSeparatrices ( const char *name )
{
  int i,j;
  ofstream ofs(name,ios::binary);

  if (!ofs)
    {
      cout << "Can't open file " << name << " for writing the output..." << endl;
      return;
    }

  mark();

  // start with vertices
  for ( i=0; i<n.size(); i++ )
    {
      if (!n[i] || !(n[i]->flags & (1|2)))
	continue;

      if (n[i]->owner->dimension==0 && n[i]->scc==-1)
	;
	//	primitive(0,0,0,255,false,
	//  msh->vertex(n[i]->owner->ID)).save(ofs);
      else
	if (n[i]->owner->dimension==1 && n[i]->scc==-1)
	  primitive(0,0,0,255,false,
		    msh->edgepoint(n[i]->owner->ID,n[i]->s),
		    msh->edgepoint(n[i]->owner->ID,n[i]->e)).save(ofs);

      for ( j=0; j<n[i]->out.size(); j++ )
	{
	  arc *a = n[i]->out[j];
	  if (a->from->scc>=0 && a->from->scc==a->to->scc)
	    continue;
	  if (!(a->from->flags & a->to->flags & (1|2)))
	    continue;
	  if (a->get_dimension()<2)
	    continue;
	  if (a->from->owner->dimension==0 && a->to->owner->dimension==0)
	    primitive(0,0,0,255,false,
		      msh->vertex(a->from->owner->ID),
		      msh->vertex(a->to->owner->ID)).save(ofs);
	  else
	    if (a->from->owner->dimension==1 && a->to->owner->dimension==0)
	      {
		primitive p(0,0,0,255,false,
			    msh->edgepoint(a->from->owner->ID,a->from->s),
			    msh->edgepoint(a->from->owner->ID,a->from->e),
			    msh->vertex(a->to->owner->ID));
		p.orient(msh->normal(a->get_ix()));
		p.save(ofs);
	      }
	    else
	      if (a->from->owner->dimension==0 && a->to->owner->dimension==1)
		{
		  primitive p(0,0,0,255,false,
			      msh->edgepoint(a->to->owner->ID,a->to->s),
			      msh->edgepoint(a->to->owner->ID,a->to->e),
			      msh->vertex(a->from->owner->ID));
		  p.orient(msh->normal(a->get_ix()));
		  p.save(ofs);
		}
	      else
		if (a->from->owner->dimension==1 && a->to->owner->dimension==1)
		  {
		    primitive p(0,0,0,255,false,
				msh->edgepoint(a->to->owner->ID,a->to->s),
				msh->edgepoint(a->to->owner->ID,a->to->e),
				msh->edgepoint(a->from->owner->ID,a->from->s),
				msh->edgepoint(a->from->owner->ID,a->from->e));
		    p.orient(msh->normal(a->get_ix()));
		    p.save(ofs);		
		  }
		else
		  assert(0);
	}
    }  
}

/* ------------------------------------------------------ */

void tgraph::saveMorseSets ( const char *name )
{
  int i,j;
  ofstream ofs(name,ios::binary);

  if (!ofs)
    {
      cout << "Can't open file " << name << " for writing the output..." << endl;
      return;
    }

  bool *saved = new bool[msh->faces()];
  for ( i=0; i<msh->faces(); i++ )
    saved[i] = false;
      
  // start with vertices
  for ( i=0; i<n.size(); i++ )
    {
      if (!n[i] || n[i]->scc==-1)
	continue;

      
      if (n[i]->owner->dimension==0)
	primitive(n[i]->scc,
		  mstp[n[i]->scc].getindex(), mstp[n[i]->scc].getindex2(),
		  mstp[n[i]->scc].getstability(), mstp[n[i]->scc].getbdry(),
		  msh->vertex(n[i]->owner->ID)).save(ofs);
      else
	if (n[i]->owner->dimension==1)
	  {
	    for ( j=0; j<n[i]->owner->cofaces; j++ )
	      if (msh->iststationary(n[i]->owner->coface[j]->ID))
		{
		  mesh_element *ff = n[i]->owner->coface[j];
		  if (saved[ff->ID])
		    continue;
		  saved[ff->ID] = true;
		  for ( int k=5; k<ff->faces; k+=2 )
		    primitive(n[i]->scc,
			      mstp[n[i]->scc].getindex(), mstp[n[i]->scc].getindex2(),
			      mstp[n[i]->scc].getstability(), mstp[n[i]->scc].getbdry(),
			      msh->vertex(ff->face[1]->ID),
			      msh->vertex(ff->face[k]->ID),
			      msh->vertex(ff->face[k-2]->ID)).save(ofs);
		}
	    primitive(n[i]->scc,
		      mstp[n[i]->scc].getindex(), mstp[n[i]->scc].getindex2(),
		      mstp[n[i]->scc].getstability(), mstp[n[i]->scc].getbdry(),
		      msh->edgepoint(n[i]->owner->ID,n[i]->s),
		      msh->edgepoint(n[i]->owner->ID,n[i]->e)).save(ofs);
	  }
	else
	  assert(0);

      for ( j=0; j<n[i]->out.size(); j++ )
	{
	  arc *a = n[i]->out[j];
	  if (a->from->scc!=a->to->scc)
	    continue;
	  if (a->get_dimension()<2)
	    continue;
	  if (a->from->owner->dimension==0 && a->to->owner->dimension==0)
	    {
	      primitive(n[i]->scc,
			mstp[n[i]->scc].getindex(), mstp[n[i]->scc].getindex2(),
			mstp[n[i]->scc].getstability(), mstp[n[i]->scc].getbdry(),
			msh->vertex(a->from->owner->ID),
			msh->vertex(a->to->owner->ID)).save(ofs);
	    }
	  else
	    if (a->from->owner->dimension==1 && a->to->owner->dimension==0)
	      {
		primitive p(n[i]->scc,
			    mstp[n[i]->scc].getindex(), mstp[n[i]->scc].getindex2(),
			    mstp[n[i]->scc].getstability(), mstp[n[i]->scc].getbdry(),
			    msh->edgepoint(a->from->owner->ID,a->from->s),
			    msh->edgepoint(a->from->owner->ID,a->from->e),
			    msh->vertex(a->to->owner->ID));
		p.orient(msh->normal(a->get_ix()));
		p.save(ofs);
	      }
	    else
	      if (a->from->owner->dimension==0 && a->to->owner->dimension==1)
		{
		  primitive p(n[i]->scc,
			      mstp[n[i]->scc].getindex(), mstp[n[i]->scc].getindex2(),
			      mstp[n[i]->scc].getstability(), mstp[n[i]->scc].getbdry(),
			      msh->edgepoint(a->to->owner->ID,a->to->s),
			      msh->edgepoint(a->to->owner->ID,a->to->e),
			      msh->vertex(a->from->owner->ID));
		  p.orient(msh->normal(a->get_ix()));
		  p.save(ofs);
		}
	      else
		if (a->from->owner->dimension==1 && a->to->owner->dimension==1)
		  {
		    primitive p(n[i]->scc,
				mstp[n[i]->scc].getindex(), mstp[n[i]->scc].getindex2(),
				mstp[n[i]->scc].getstability(), mstp[n[i]->scc].getbdry(),
				msh->edgepoint(a->to->owner->ID,a->to->s),
				msh->edgepoint(a->to->owner->ID,a->to->e),
				msh->edgepoint(a->from->owner->ID,a->from->s),
				msh->edgepoint(a->from->owner->ID,a->from->e));
		    p.orient(msh->normal(a->get_ix()));
      		    p.save(ofs);		    
		  }
		else
		  assert(0);
	}
    }
      
  delete[] saved;
}

/* ------------------------------------------------------ */

void tgraph::clear_flags()
{
  for ( int i=0; i<n.size(); i++ )
    if (n[i])
      n[i]->flags &= (16|32);
}

/* ------------------------------------------------------ */

void tgraph::mark()
{
  int i;

  clear_flags();
  
  for ( i=0; i<n.size(); i++ )
    if (n[i] && n[i]->scc>=0 && 
	!mstp[n[i]->scc].istrivial() && 
	!mstp[n[i]->scc].isattracting() && !mstp[n[i]->scc].isrepelling())
      {
	n[i]->markf(1);
	n[i]->markb(2);
      }
}

/* ------------------------------------------------------ */

void tgraph::traverse_and_subdivide ( int start, int minlevel, int maxlevel, char dir )
{
  int i;
  node *startnode = n[start];
  double tspan = 1.01/(1<<minlevel);
  double ltspan = 1.01/(1<<maxlevel);
  double cspan = tspan;
  vector<node*> traversed;
  int phase = 1;
  bool found;

  // NOTE: can probably just do phase 1 for backward traversals

  while(1)
    {
      traversed.clear();
      found = false;

      stack<node*> S;
      S.push(startnode);
      n[startnode->ID]->flags |= 1;

      vector<node*> newnodes;

      while (!S.empty())
	{
	  node *cn = S.top();
	  S.pop();
	  if (cn->scc!=-1 && cn->scc!=startnode->scc && 
	      !mstp[cn->scc].isattracting() && !mstp[cn->scc].isrepelling() &&
	      !mstp[cn->scc].istrivial())
	    found = true;

	  traversed.push_back(cn);	
	  
	  vector<node*> toBsubdivided;
	  for ( i=0; i<((dir=='f') ? cn->out : cn->in).size(); i++ )
	    if (((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->isepiece() && 
		((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->scc==-1 && 
		((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->span()>cspan && 
		!(((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->flags&2))
	      {
		toBsubdivided.push_back((dir=='f') ? cn->out[i]->to : cn->in[i]->from);
		((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->flags |= 2;   // newly added line
	      }
	  for ( i=0; i<toBsubdivided.size(); i++ )
	    {
	      int id = toBsubdivided[i]->ID;
	      subdivide(id);
	      n[id]->flags |= 2;
	      n[n.size()-1]->flags |= 2;
	      newnodes.push_back(n[id]);
	      newnodes.push_back(n[n.size()-1]);
	    }

	  for ( i=0; i<((dir=='f') ? cn->out : cn->in).size(); i++ )
	    if (!(((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->flags&1))
	      {
		S.push((dir=='f') ? cn->out[i]->to : cn->in[i]->from);
		((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->flags |= 1;
	      }
	}

      // clear 2 flag
      for ( i=0; i<newnodes.size(); i++ )
	newnodes[i]->flags &= ~2;

      if (phase==2 && (!found || !newnodes.size()))
	break;

      if (phase==1) 
	{
	  if (!newnodes.size())
	    if (found)
	      {
		phase = 2;
		cspan = ltspan;
	      }
	    else
	      break;
	}
	
      // clear 1 flag
      for ( i=0; i<traversed.size(); i++ )
	traversed[i]->flags &= ~1;
    }

  // OK, now merge and protect
  for ( i=0; i<traversed.size(); i++ )
    {
      node *cn = traversed[i];
      if (cn->scc>=0)
	continue;
      cn->flags |= 2;
      assert(cn->flags&1);
      int cnid = cn->ID;
      if (cn->left && (cn->left->flags&2) && !cn->islockedL())
	{
	  node *nn = merge(cn->ID,'l');
	  nn->flags |= (1|2);
	  cn = nn;
	}
      cnid = cn->ID;
      if (cn->right && (cn->right->flags&2) && !cn->islockedR())
	{
	  node *nn = merge(cn->ID,'r');
	  nn->flags |= (1|2);
	  cn = nn;
	}
      if (cn->right && !(cn->right->flags&1))
	cn->lockR();
      if (cn->left && !(cn->left->flags&1))
	cn->lockL();
    }

  // finally, reset all flags; cannot rely on any of the vectors though
  {
    stack<node*> S;
    S.push(startnode);
    while (!S.empty())
      {
	node *cn = S.top();
	S.pop();
	cn->flags &= ~(1|2);
	for ( i=0; i<((dir=='f') ? cn->out : cn->in).size(); i++ )
	  if (((dir=='f') ? cn->out[i]->to : cn->in[i]->from)->flags & (1|2))
	    S.push((dir=='f') ? cn->out[i]->to : cn->in[i]->from);
      }
  }
}

/* ------------------------------------------------------ */

void tgraph::prepare4MCG ( int minlevel, int maxlevel )
{
  int i;
  int count = 0;
  int all = 0;

  bool *done = new bool[SCCs()];
  for ( i=0; i<SCCs(); i++ )
    {
      done[i] = false;
      if (mstp[i].isattracting() || mstp[i].isrepelling())
	continue;
      if (mstp[i].istrivial())
	continue;
      all++;
    }

  cout << "Preparing graph for MCG computation: " << endl;

  for ( i=0; i<n.size(); i++ )
    {
      if (n[i]->scc==-1)
	continue;
      if (done[n[i]->scc])
	continue;
      if (mstp[n[i]->scc].isattracting() || mstp[n[i]->scc].isrepelling())
	continue;
      if (mstp[n[i]->scc].istrivial())
	continue;
      done[n[i]->scc] = true;

      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << count++ << " / " << all << flush;

      traverse_and_subdivide(i,minlevel,maxlevel,'f');
      traverse_and_subdivide(i,minlevel,maxlevel,'b');
    }

  delete[] done;
}

/* ------------------------------------------------------ */

void tgraph::coarsenMorseSets()
{
  bool cont;
  
  do
    {
      cont = false;

      for ( int i=0; i<n.size(); i++ )
	{
	  if (n[i]->left && !n[i]->islockedL() && n[i]->scc>=0 && n[i]->left->scc==n[i]->scc)
	    {
	      merge(n[i]->ID,'l');
	      cont = true;
	    }

	  if (i==n.size())
	    continue;

	  if (n[i]->right && !n[i]->islockedR() && n[i]->scc>=0 && n[i]->right->scc==n[i]->scc)
	    {
	      merge(n[i]->ID,'r');
	      cont = true;
	    }
	}
    }
  while(cont);
}

/* ------------------------------------------------------ */

void tgraph::coarsenMorseSetsP()
{
  for ( int e=0; e<n.size(); )
    {
      node *ni = n[e];
	      
      if (ni && ni->left && !ni->islockedL() && ni->scc>=0 && ni->left->scc==ni->scc)
	merge(ni->ID,'l',false);
      else
	e++;
    }
	  
  // need to compress n now...
  for ( int k=0; k<n.size(); )
    if (!n[k])
      {
	n[k] = n[n.size()-1];
	n.pop_back();
      }
    else
      k++;
	  
  for ( int k=0; k<n.size(); k++ )
    {
      assert(n[k]);
      n[k]->ID = k;
    }
}

/* ------------------------------------------------------ */

tskel *tgraph::MCG (  bool include_trivial )
{
  int i,j;

  vector<node*> traversed;
  tskel * res = new tskel(SCCs(),mstp);
  bool *done = new bool[SCCs()];

  for ( i=0; i<SCCs(); i++ )
    done[i] = false;

  // construct arcs

  for ( i=0; i<n.size(); i++ )
    if (n[i] && n[i]->scc>=0 && (!mstp[n[i]->scc].istrivial() || include_trivial) && !done[n[i]->scc])
      {
	bool *added = new bool[SCCs()];
	
	for ( j=0; j<SCCs(); j++ )
	  added[j] = false;

	// traverse...
	traversed.clear();
	stack<node*> S;
	S.push(n[i]);
	done[n[i]->scc] = true;

	while(!S.empty())
	  {
	    node *cn = S.top();
	    S.pop();
	    
	    cn->flags |= 8;
	    traversed.push_back(cn);

	    if (cn->scc>=0 && cn->scc!=n[i]->scc && !added[cn->scc] && (include_trivial || !mstp[cn->scc].istrivial()))
	      {
		res->add_edge(n[i]->scc,cn->scc);
		added[cn->scc] = true;
	      }
	    
	    for ( j=0; j<cn->out.size(); j++ )
	      if (!(cn->out[j]->to->flags & 8))
		S.push(cn->out[j]->to);
	  }
	
	for ( j=0; j<traversed.size(); j++ )
	  traversed[j]->flags &= ~8;

	delete[] added;
      }

  delete[] done;

  res->cleanup();

  return res;
}

/* ------------------------------------------------------ */

void node::print_out (  )
{
  int i;
  cout << "Node " << "[" << owner->ID << " " << owner->dimension << " " << s << " " << e << "]" << endl;
  if (left)
    cout << "L: " <<  "[" << left->owner->ID << " " << left->owner->dimension << " " 
	 << left->s << " " << left->e << "]"  << endl;
  if (right)
    cout << "R: " <<  "[" << right->owner->ID << " " << right->owner->dimension << " " 
	 << right->s << " " << right->e << "]" << endl;
  cout << "SCC: " << scc << endl;
  cout << "Arcs out: ";
  for ( i=0; i<out.size(); i++ )
    cout << "[" << out[i]->to->owner->ID  << " " << out[i]->to->owner->dimension << " " 
	 << out[i]->to->s << " " << out[i]->to->e << "] ";
  cout << endl << "Arcs in: ";
  for ( i=0; i<in.size(); i++ )
    cout << "[" << in[i]->from->owner->ID  << " " << in[i]->from->owner->dimension << " " 
	 << in[i]->from->s << " " << in[i]->from->e << "] ";
  cout << endl;
}

void tgraph::print_out()
{
  for ( int i=0; i<n.size(); i++ )
    if (n[i])
      n[i]->print_out();
  cout << "Number of SCCs: " << SCCs() << endl;
}

/* ------------------------------------------------------ */

double tgraph::usedspace()
{
  int totalc = 0;
  int totals = 0;

  for ( int i=0; i<n.size(); i++ )
    {
      totalc += n[i]->in.capacity();
      totals += n[i]->in.size();
      totalc += n[i]->out.capacity();
      totals += n[i]->out.size();
    }

  return double(totals)/totalc;
}

/* ------------------------------------------------------ */

tgraph::~tgraph()
{
  /*
  while(nodes())
    {
      remove_node(0);
    }
  */

  delete msh;
  msh = NULL;
  if (mstp) delete[] mstp;
  mstp = NULL;
}

/* ------------------------------------------------------ */
