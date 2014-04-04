

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

#if(!defined(__TSKEL_H))
#define __TSKEL_H

#include <global.h>
#include <vector>
#include <mstype.h>

/* ------------------------------------------------------------------ */

class tsknode {
 public:
  mstype t;
  std::vector<int> in;
  std::vector<int> out;

  tsknode ();
  tsknode ( int id, mstype ms );

  void removein ( int a );
  void removeout ( int a );

  bool hasedgeto ( int i );
  bool hasedgefrom ( int i ); 

  bool isattracting();
  bool isrepelling();
};


/* ------------------------------------------------------------------ */

class tskel {
 public:
  int ns;
  tsknode *n;

  tskel ( int nn, mstype *t );
  ~tskel();

  void add_edge ( int i, int j );
  void remove_edge ( int i, int j );
  void save ( const char *name );
  void cleanup();
  bool hasedge ( int i, int j );
  bool iscertain ( int i, int j );
};

/* ------------------------------------------------------------------ */

#endif
