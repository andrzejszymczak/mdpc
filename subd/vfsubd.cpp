

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

#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>

#include <mesh.h>
#include <shared.h>

/* ------------------------------------------------ */

double ncf ( int deg )
{
  if (deg==3) return 3.0/32.0;
  return 3.0/(16.0*deg);
}

double tcf ( int deg )
{
  if (deg==3) return 7.0/16.0;
  return 5.0/8.0;
}

/* ------------------------------------------------ */

int main ( int argc, char *argv[] )
{
  FILE *fp;
  int i,j,eid;
  M_mesh msh;
  M_int_triple *enmb;

  int nnt,nnv;
  M_int_triple *ninc;
  M_point *nv;
  M_point *f, *nf;

  struct timeb start,finish;

  /* check if commandline fulfils syntax */

  if (argc!=3)
    {
      printf("Usage: subdiv <input .t file> <output .t file>\n");
      return 0;
    } 

  ftime(&start);

  /* read the input file */

  fp = fopen(argv[1],"r");
  if (!fp)
    {
      printf("Can't open input file %s\n",argv[1]);
      return 0;
    }
  msh = M_read_mesh(fp);
  f = (M_point*)malloc(msh.vertices*sizeof(M_point));

  for ( i=0; i<msh.vertices; i++ )
    fscanf(fp,"%lf %lf %lf",&f[i][0],&f[i][1],&f[i][2]);

  fclose(fp);

  /* number edges */

  enmb = (M_int_triple*)malloc(msh.triangles*sizeof(M_int_triple));

  for ( eid=i=0; i<msh.triangles; i++ )
    for ( j=0; j<3; j++ )
      if (msh.adj[i][j]>i)
	enmb[i][j] = enmb[msh.adj[i][j]][which(i,msh.adj[msh.adj[i][j]])] = eid++;

  /* now compute the new vertex data */

  nnv = eid+msh.vertices;
  nnt = 4*msh.triangles;

  nv = (M_point*)malloc(nnv*sizeof(M_point));
  nf = (M_point*)malloc(nnv*sizeof(M_point));
  ninc = (M_int_triple*)malloc(nnt*sizeof(M_int_triple));

  for ( i=0; i<nnv; i++ )
    {
      nv[i][0] = nv[i][1] = nv[i][2] = 0.0;
      nf[i][0] = nf[i][1] = nf[i][2] = 0.0;
    }

  for ( i=0; i<msh.triangles; i++ )
    {
      ninc[4*i+0][0] = msh.triangle[i][0];
      ninc[4*i+0][1] = msh.vertices+enmb[i][2];
      ninc[4*i+0][2] = msh.vertices+enmb[i][1];
      ninc[4*i+1][0] = msh.triangle[i][1];
      ninc[4*i+1][1] = msh.vertices+enmb[i][0];
      ninc[4*i+1][2] = msh.vertices+enmb[i][2];
      ninc[4*i+2][0] = msh.triangle[i][2];
      ninc[4*i+2][1] = msh.vertices+enmb[i][1];
      ninc[4*i+2][2] = msh.vertices+enmb[i][0];
      ninc[4*i+3][2] = msh.vertices+enmb[i][2]; 
      ninc[4*i+3][1] = msh.vertices+enmb[i][1];
      ninc[4*i+3][0] = msh.vertices+enmb[i][0];
    }

  for ( i=0; i<msh.triangles; i++ )
    {
      for ( j=0; j<3; j++ )
	{
	  M_point_scale_and_add(nv[msh.vertices+enmb[i][j]],3.0/16.0,msh.vertex[msh.triangle[i][(j+1)%3]]);
	  M_point_scale_and_add(nv[msh.vertices+enmb[i][j]],3.0/16.0,msh.vertex[msh.triangle[i][(j+2)%3]]);
	  M_point_scale_and_add(nv[msh.vertices+enmb[i][j]],1.0/8.0,msh.vertex[msh.triangle[i][j]]);
	  M_point_scale_and_add(nv[msh.triangle[i][j]],ncf(msh.deg[msh.triangle[i][j]]),
				msh.vertex[msh.triangle[i][(j+1)%3]]);
	  M_point_scale_and_add(nv[msh.triangle[i][j]],ncf(msh.deg[msh.triangle[i][j]]),
				msh.vertex[msh.triangle[i][(j+2)%3]]);

	  M_point_scale_and_add(nf[msh.vertices+enmb[i][j]],.25,f[msh.triangle[i][(j+1)%3]]);
	  M_point_scale_and_add(nf[msh.vertices+enmb[i][j]],.25,f[msh.triangle[i][(j+2)%3]]);
	}
    }

  for ( i=0; i<msh.vertices; i++ )
    {
      M_point_scale_and_add(nv[i],tcf(msh.deg[i]),msh.vertex[i]);
      M_point_scale_and_add(nf[i],1,f[i]);
    }

  /* save the output file */

  fp = fopen(argv[2],"w");
  if (!fp)
    {
      printf("Can't open output file %s\n",argv[2]);
      return 0;
    }

  fprintf(fp,"%d %d\n\n",nnt,nnv);
  for ( i=0; i<nnt; i++ )
    fprintf(fp,"%d %d %d\n",ninc[i][0],ninc[i][1],ninc[i][2]);
  for ( i=0; i<nnv; i++ )
    fprintf(fp,PFORMAT,nv[i][0],nv[i][1],nv[i][2]);
  for ( i=0; i<nnv; i++ )
    fprintf(fp,PFORMAT,nf[i][0],nf[i][1],nf[i][2]);

  fclose(fp);

  ftime(&finish);

  printf("time: %lf s\n",(double)(1000*(finish.time-start.time)+finish.millitm-start.millitm)/1000);

  return 0;
}
