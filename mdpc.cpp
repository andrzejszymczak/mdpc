
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
#include <pcenv.h>
#include <tgraph.h>

using namespace std;

/* ------------------------------------------------------ */

// basic Morse set types for the interior of the domain

int ix[6][3] = 
  {
    {  0,  0,  1 },  // attracting PO
    {  0,  0,  2 },  // repelling PO
    { -1, -1,  3 },  // saddle
    {  1,  1,  1 },  // sink
    {  1,  1,  2 },  // source
    {  0,  0,  3 }   // trivial
  };

const char *it[7] = 
  {
    "Attracting periodic trajectories: ",
    "Repelling periodic trajectories: ",
    "Saddles: ",
    "Sinks: ",
    "Sources: ",
    "Trivial: ",
    "Other: "
  };

int ic[7] = { 0, 0, 0, 0, 0, 0, 0 };

/* ------------------------------------------------------ */

// basic Morse sets touching the boundary of the domain

int bx[7][3] = 
  {
    {  1,  0,  1 },   // simple boundary sink
    {  0,  1,  2 },   // simple boundary source
    {  0, -1,  3 },   // converging half-saddle
    { -1,  0,  3 },   // diverging half-saddle
    {  0,  0,  3 },   // trivial
    {  0,  0,  2 },    // repelling periodic orbit
    {  0,  0,  1 }   // attracting periodic orbit
  };

const char *bt[8] = 
  {
      "Simple boundary sinks: ",
      "Simple boundary sources: ",
      "Converging half-saddles: ",
      "Diverging half-saddles: ",
      "Trivial: ",
      "Repelling periodic orbits:",
      "Attracting periodic orbits:",
      "Other: "
    };

int bc[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

/* ------------------------------------------------------ */

void print_statistics ()
{
  cout << "---------- MORSE SET COUNTS -----------" << endl;
  for ( int i=0; i<7; i++ )
    cout << " " << it[i] << ic[i] << endl;
  if (bc[0]+bc[1]+bc[2]+bc[3]+bc[4]+bc[5]+bc[6]+bc[7])
    {
      cout << "Morse sets touching the boundary: " << endl;
      for ( int i=0; i<8; i++ )
	cout << "   " << bt[i] << bc[i] << endl;
    }
}

/* ------------------------------------------------------ */

void new_MS ( mstype t )
{
  int j;

  if (t.getbdry())
    {
      for ( j=0; j<7; j++ )
	if (bx[j][0]==t.getindex() && bx[j][1]==t.getindex2() && bx[j][2]==t.getstability())
	  {
	    bc[j]++;
	    break;
	  }
      if (j==7)
	bc[j]++;
    }
  else
    {
      for ( j=0; j<6; j++ )
	if (ix[j][0]==t.getindex() && ix[j][1]==t.getindex2() && ix[j][2]==t.getstability())
	  {
	    ic[j]++;
	    break;
	  }
      if (j==6)
	ic[j]++;
    }
}

/* ------------------------------------------------------ */

void print_usage()
{
  cout << "Usage: " << endl;
  cout << " mdpc [options] <IN> <N> <OUT-MD> [<OUT-SEP> <OUT-DOT>] " << endl;
  cout << "   IN: input file name" << endl;
  cout << "   N : refinement iterations (when computing Morse decomposition)" << endl;
  cout << " Output file names: " << endl;
  cout << "   OUT-MD : geometric model of Morse decomposition" << endl;
  cout << "  If -c option is used: " << endl;
  cout << "   OUT-SEP: geometric model of separatrices" << endl;
  cout << "   OUT-DOT: dot file containing the MCG (graphviz can be used to draw the MCG)" << endl;
  cout << "  options: " << endl;
  cout << "   -c <MIN> <MAX>: compute MCG; MIN and MAX are refinement thresholds " << endl;
  cout << "   -s <R>        : compute stable Morse decomposition with radius R" << endl;
  cout << "   -h <WEIGHT>   : compute hull-based stable Morse decomposition" << endl;
  cout << "   -v            : vertex based input" << endl;
  cout << "   -e <WEIGHT>   : envelope (only with -v)" << endl;
  cout << "   -t            : include trivial Morse sets in the MCG" << endl;
  cout << "   -o            : treat the vector field as open system (allow flow into/out of domain)" << endl;
}

/* ------------------------------------------------------ */

static bool hopt = false;
static bool copt = false;
static bool eopt = false;
static bool sopt = false;
static char type = 'f';
static bool inct = false;
static bool osys = false;

static double R,wt;
static int minr, maxr;

/* ------------------------------------------------------ */

int main ( int argc, char *argv[] )
{

  int i = 1;
  int j;

  if (argc==1)
    {
      print_usage();
      return 0;
    }

  {
    // first, parse arguments
    do {
      if (argc<=i)
	{
	  print_usage();
	  return 0;
	}
      if (argv[i][0]=='-')
	{
	  switch(argv[i][1])
	    {
	    case 'c':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      copt = true;
	      if (argc<i+3)
		{
		  cout << "-c has to be followed by two nonnegative integers" << endl;
		  return 0;
		}
	      minr = atoi(argv[i+1]);
	      maxr = atoi(argv[i+2]);
	      if (minr<0 || maxr<0)
		{
		  cout << "-c has to be followed by two nonnegative integers" << endl;
		  return 0;
		}
		
	      i += 3;
	      break;

	    case 'h':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      hopt = true;
	      if (argc<i+2)
		{
		  cout << "-h has to be followed by a floating point number" << endl;
		  return 0;
		}
	      wt = atof(argv[i+1]);
	      i += 2;
	      break;
	      
	    case 's':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      sopt = true;
	      if (argc<i+2)
		{
		  cout << "-s has to be followed by a floating point number" << endl;
		  return 0;
		}
	      R = atof(argv[i+1]);
	      i += 2;
	      break;

	    case 'v':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      type = 'v';
	      i++;
	      break;


	    case 't':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      inct = true;
	      i++;
	      break;

	    case 'o':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      osys = true;
	      i++;
	      break;

	    case 'e':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      eopt = true;
	      if (argc<i+2)
		{
		  cout << "-e has to be followed by a floating point number" << endl;
		  return 0;
		}
	      wt = atof(argv[i+1]);
	      i += 2;
	      break;
	      
	    default:
	      cout << "Unknown option:" << argv[i] << endl;
	      print_usage();
	      return 0;
	      break;
	    }
	}
      else
	break;
    }
    while(1);
  }

  pcvf *m = NULL;

  if (argc<i+2)
    {
      cout << "Required args missing: input file or number of refinement iterations." << endl;
      return 0;
    }

  if (sopt)
    m = new pcstable(R,argv[i],type,!osys);
  else
    if (hopt)
      m = new pchull(wt,argv[i],type,!osys);
    else
      if (eopt)
	{
	  if (type=='f')
	    {
	      cout << "Envelope option is available only with vertex-based vector fields (-v option)" << endl;
	      return 0;
	    }
	  m = new pcenv(wt,argv[i],!osys);
	}
      else
	m = new pcvf(argv[i],type,!osys);

  tgraph t(m);
  t.computeSCCs();

  int iters = atoi(argv[i+1]);
  assert(iters>=0);

  cout << "Subdivision: " << endl;
  for ( j=1; j<=iters; j++ )
    {
      cout << j << ": " << flush;
      cout << " --> " << t.nodes() << "/" << t.arcs() << flush;
      t.subdivide_scc_nodes();
      cout << " --> " << t.nodes() << "/" << t.arcs() << flush;
      t.computeSCCs();
      if (!copt)
	t.remove_all_nonSCC();
      cout << " --> " << t.nodes() << "/" << t.arcs() << endl;
    }

  cout << endl;

  t.coarsenMorseSets();

  t.computeSCCs();
  t.computeMSTypes();

  int totali = 0;
  int totali2 = 0;

  for ( j=0; j<t.SCCs(); j++ )
    {
      totali += t.MStype(j).getindex();
      totali2 += t.MStype(j).getindex2();
      new_MS(t.MStype(j));
    }

  cout << "Total index: " << totali << " " << totali2 << endl;
  print_statistics();
  if (argc>i+2) t.saveMorseSets(argv[i+2]);

  if (copt)
    {
      cout << "Computing MCG..." << endl;
      t.prepare4MCG(minr,maxr);
      tskel *MCG = t.MCG(inct);
      if (argc>i+3) t.saveSeparatrices(argv[i+3]);
      if (argc>i+4) MCG->save(argv[i+4]);
      delete MCG;
    }

  cout << endl;

  return 0;
}

/* ------------------------------------------------------ */
