

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

#ifndef __TGRAPH_H
#define __TGRAPH_H

#include <global.h>
#include <vfield_base.h>
#include <mstype.h>
#include <tskel.h>
#include <vector>
#include <stack>

class arc;
class node;
class tgraph;

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

// current limitation: maximum ID of a face, vertex or edge < 2^22

class node {

  friend class tgraph;
  friend class arc;

  std::vector<arc*> in;
  std::vector<arc*> out;
  mesh_element *owner;

  // properties below only for edge pieces
  node *left,*right; // neighboring half-edges
  float s,e;        // start and end parameter
  int scc;           // strongly connected component ID

  // bits 0, 1, 2 and 3 used for MCG computation, bits 6 and 7: when merging nodes
  // bits 4 and 5 used to protect from mergers
  unsigned char flags;

  int ID;

 public:

  static int nodes;

  node ( int id, mesh_element *o );  // use this only for vertex nodes
  node ( int id, mesh_element *o, node *l, node *r, double ss, double ee ); // only for edge pieces  
  ~node();

  bool isepiece();

  void print_out ( );

  // related to MCG computation
  // in both cases, return nodes not in scc but in (un)stable set
  void markf ( unsigned char msk, std::vector<node*> *lst = NULL );  // forward DFS
  void markb ( unsigned char msk, std::vector<node*> *lst = NULL );  // backward DFS

  double span();

  // protect endpoints
  void lockL();
  void lockR();
  bool islockedL();
  bool islockedR();
};

/* ------------------------------------------------------ */

class arc {

  friend class node;
  friend class tgraph;

  node *from;
  node *to;
  unsigned int data;

 public:

  static int arcs;

  arc ( node *f, node *t, mesh_element *carrier, int ix_orig, int ix_dest );  // 2D carrier only
  arc ( node *f, node *t, mesh_element *carrier );
  arc ( node *f, node *t, arc *a );   // copies data to the new arc
  ~arc();

  // extract carrier data from data field
  int get_dimension();
  int get_ixf();  // only for 2D carriers
  int get_ixt();  // only for 2D carriers
  int get_ix();
};

/* ------------------------------------------------------ */

class tgraph {

  std::vector<node*> n;   // all graph nodes

  // related to strongly connected component computation...
  std::stack<node*> *__S;
  int __index;
  int __c;
  void __visit ( node *v );
  int sccs;
  mstype *mstp;

  // index in the coarse graph ONLY
  int _index ( mesh_element *e ); 
  node *_getnode ( mesh_element *e );
  node *_getvertexnode ( int i );
  node *_getedgenode ( int i );

  void clear_flags();

  // MCG related calls
  void mark();   // mark nodes on generalized separatrices; assumes complete Morse set data
  void traverse_and_subdivide ( int start, int minlevel, int maxlevel, char dir );

  // node removal
  void remove_node ( int i );

  // refinement related calls
  void subdivide ( int i, double spt = 0.5 );

  // side can be 'l', 'r', 'L' or 'R' (to merge with left or right neighbor)
  // merge just copies the attrtibutes (scc, flags) from node i to the new node
  node* merge ( int i, char side, bool move = true );


 protected:
  vfield_base *msh;

 public:

  tgraph ( vfield_base *m );   // build a coarse graph
  ~tgraph();

  void print_out();
  double usedspace();  // size-to-capacity ratio for edge lists (just a statistics)

  void subdivide_all();
  void subdivide_scc_nodes();

  // merge edge pieces in the same Morse set
  void coarsenMorseSets();  
  void coarsenMorseSetsP();  

  // Morse set computation related calls
  void computeSCCs();
  void computeMSTypes();
  int SCCs();     // returns number of SCCs
  mstype MStype ( int i );
  void saveMorseSets ( const char *name ); // assumes up to date SCC and MS type info

  // MCG computation; assumes Morse sets and their types are up to date
  void prepare4MCG ( int minlevel, int maxlevel );  // uses bits 0 and 1 of flags
  void saveSeparatrices ( const char *name );
  tskel *MCG ( bool include_trivial = false );

  // remove all nodes except adjacent to an scc
  void remove_all_nonSCC();

  int nodes();
  int arcs();
  int arcs2();
};

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */

#endif
