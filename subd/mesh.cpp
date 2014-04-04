
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

/* --- mesh.c --- */

#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <stdio.h>

#include <mesh.h>
#include <shared.h>

/* -------------------------------- */

typedef int qdt[4];

/* -------------------------------- */

int which ( int i, M_int_triple t )
{
  if (i==t[0])
    return 0;
  if (i==t[1])
    return 1;
  assert(i==t[2]);
  return 2;
}

/* -------------------------------- */

void M_free_mesh ( M_mesh *m )
{
  int i;

  free(m->vertex);
  free(m->deg);
  free(m->adj);
  free(m->triangle);
  free(m->normal);

  for ( i=0; i<m->vertices; i++ )
    {
      free(m->inc[i]);
    }

  free(m->inc);

  m->triangles = m->vertices = 0;
}

/* -------------------------------- */

void M_point_subtract ( M_point t, M_point a, M_point b )
{
  t[0] = a[0] - b[0];
  t[1] = a[1] - b[1];
  t[2] = a[2] - b[2];
}

/* -------------------------------- */

double M_point_dot ( M_point a, M_point b )
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

/* -------------------------------- */

void M_point_add ( M_point t, M_point a, M_point b )
{
  t[0] = a[0] + b[0];
  t[1] = a[1] + b[1];
  t[2] = a[2] + b[2];
}

/* -------------------------------- */

void M_point_cross ( M_point t, M_point a, M_point b )
{
  t[0] = a[1]*b[2]-a[2]*b[1];
  t[1] = a[2]*b[0]-a[0]*b[2];
  t[2] = a[0]*b[1]-a[1]*b[0];
}

/* -------------------------------- */

void M_point_scale ( M_point t, double c, M_point a )
{
  t[0] = c*a[0];
  t[1] = c*a[1];
  t[2] = c*a[2];
}

/* -------------------------------- */

void M_point_scale_and_add ( M_point t, double c, M_point a )
{
  t[0] += c*a[0];
  t[1] += c*a[1];
  t[2] += c*a[2];
}

/* -------------------------------- */

double M_point_norm ( M_point p )
{
  return sqrt(M_point_dot(p,p));
}

/* -------------------------------- */

double M_point_area ( M_point p1, M_point p2, M_point p3 )
{
  M_point p12,p13,c;
  M_point_subtract(p12,p2,p1);
  M_point_subtract(p13,p3,p1);
  M_point_cross(c,p12,p13);
  return 0.5*M_point_norm(c);
}

/* -------------------------------- */

void M_point_normalize ( M_point p )
{
  M_point_scale(p,1/sqrt(M_point_dot(p,p)),p);
}

/* -------------------------------- */

static void qd_normalize ( qdt *q )
{
  if ((*q)[0]>(*q)[1])
    {
      int tmp = (*q)[0];
      (*q)[0] = (*q)[1];
      (*q)[1] = tmp;
    }
}

/* -------------------------------- */

static int qd_cmp ( const void *v0, const void *v1 )
{
  qdt *q0 = (qdt*)v0;
  qdt *q1 = (qdt*)v1;

  if ((*q0)[0]<(*q1)[0])
    return -1;
  if ((*q0)[0]>(*q1)[0])
    return 1;
  if ((*q0)[1]<(*q1)[1])
    return -1;
  if ((*q0)[1]>(*q1)[1])
    return 1;
  return 0;
}

/* -------------------------------- */

M_mesh M_read_mesh ( FILE *fp )
{
  M_mesh res;
  int i,j,di;
  qdt *q;
  int *w;

  fscanf(fp,"%d %d",&res.triangles,&res.vertices);

  res.vertex = (M_point*)malloc(res.vertices*sizeof(M_point));
  res.triangle = (M_int_triple*)malloc(res.triangles*sizeof(M_int_triple));
  res.adj = (M_int_triple*)malloc(res.triangles*sizeof(M_int_triple));
  res.deg = (int*)malloc(res.vertices*sizeof(int));

  memset(res.deg,0,res.vertices*sizeof(int));

  for ( i=0; i<res.triangles; i++ )
    {
      fscanf(fp,"%d %d %d",&res.triangle[i][0],
	     &res.triangle[i][1],&res.triangle[i][2]);

      assert(res.triangle[i][0]>=0 && res.triangle[i][0]<res.vertices);
      assert(res.triangle[i][1]>=0 && res.triangle[i][1]<res.vertices);
      assert(res.triangle[i][2]>=0 && res.triangle[i][2]<res.vertices);
      
      res.deg[res.triangle[i][0]]++;
      res.deg[res.triangle[i][1]]++;
      res.deg[res.triangle[i][2]]++;
    }

  res.inc = (M_pair**)malloc(res.vertices*sizeof(M_pair*));

  for ( i=0; i<res.vertices; i++ )
    res.inc[i] = res.deg ? (M_pair*)malloc(res.deg[i]*sizeof(M_pair)) : NULL;

  w = (int*)malloc(res.vertices*sizeof(int));
  memset(w,0,res.vertices*sizeof(int));

  for ( i=0; i<res.triangles; i++ )
    {
      res.inc[res.triangle[i][0]][w[res.triangle[i][0]]].t = i;
      res.inc[res.triangle[i][0]][w[res.triangle[i][0]]++].ix = 0;
      res.inc[res.triangle[i][1]][w[res.triangle[i][1]]].t = i;
      res.inc[res.triangle[i][1]][w[res.triangle[i][1]]++].ix = 1;
      res.inc[res.triangle[i][2]][w[res.triangle[i][2]]].t = i;
      res.inc[res.triangle[i][2]][w[res.triangle[i][2]]++].ix = 2;
    }

  for ( i=0; i<res.vertices; i++ )
    assert(w[i]==res.deg[i]);

  free(w);

  for ( i=0; i<res.vertices; i++ )
    fscanf(fp,"%lf %lf %lf",&res.vertex[i][0],
	  &res.vertex[i][1],&res.vertex[i][2]);

  q = (qdt*)malloc(3*res.triangles*sizeof(qdt));

  for ( j=i=0; i<res.triangles; i++ )
    {
      q[j][0] = res.triangle[i][0];
      q[j][1] = res.triangle[i][1];
      q[j][2] = i;
      q[j][3] = 2;
      qd_normalize(&q[j++]);

      q[j][0] = res.triangle[i][1];
      q[j][1] = res.triangle[i][2];
      q[j][2] = i;
      q[j][3] = 0;
      qd_normalize(&q[j++]);

      q[j][0] = res.triangle[i][2];
      q[j][1] = res.triangle[i][0];
      q[j][2] = i;
      q[j][3] = 1;
      qd_normalize(&q[j++]);
    }

  assert(j==3*res.triangles);

  qsort(q,j,sizeof(qdt),&qd_cmp);

  res.internal_edges = 0;
  res.external_edges = 0;

  for ( i=0; i<3*res.triangles; i+=di )
    {
      if (i && q[i][0]==q[i-1][0] && q[i][1]==q[i-1][1])
	{
	  printf("An edge with 3 incident triangles found - exiting...\n");
	  printf("triangle(%d %d) and triangle(%d %d)\n",q[i][2],q[i][3],
		 q[i-1][2],q[i-1][3]);

	  printf("(%lf %lf %lf) (%lf %lf %lf)\n",
		 res.vertex[res.triangle[q[i][2]][(q[i][3]+1)%3]][0],
		 res.vertex[res.triangle[q[i][2]][(q[i][3]+1)%3]][1],
		 res.vertex[res.triangle[q[i][2]][(q[i][3]+1)%3]][2],
		 res.vertex[res.triangle[q[i][2]][(q[i][3]+2)%3]][0],
		 res.vertex[res.triangle[q[i][2]][(q[i][3]+2)%3]][1],
		 res.vertex[res.triangle[q[i][2]][(q[i][3]+2)%3]][2]);
	  exit(1);
	}
      if ((i<3*res.triangles-1 && (q[i][0]!=q[i+1][0] || q[i][1]!=q[i+1][1])) 
	  || i==3*res.triangles-1) /* NEW!!! */
	{
	  res.adj[q[i][2]][q[i][3]] = -1;
	  res.external_edges++;
	  di = 1;
	}
      else
	{
	  res.adj[q[i][2]][q[i][3]] = q[i+1][2];
	  res.adj[q[i+1][2]][q[i+1][3]] = q[i][2];
	  res.internal_edges++;
	  di = 2;
	}
    }

  free(q);

  res.normal = (M_point*)malloc(res.triangles*sizeof(M_point));

  for ( i=0; i<res.triangles; i++ )
    {
      M_point v1,v2;

      M_point_subtract(v1,res.vertex[res.triangle[i][0]],res.vertex[res.triangle[i][1]]);
      M_point_subtract(v2,res.vertex[res.triangle[i][0]],res.vertex[res.triangle[i][2]]);
      M_point_cross(res.normal[i],v1,v2);
      M_point_normalize(res.normal[i]);
    }

  return res;
}
/* -------------------------------- */

M_mesh M_read_mesh_X ( FILE *fp, unsigned int flags )
{
  M_mesh res;
  int i,j,di;
  qdt *q;
  int *w;

  fscanf(fp,"%d %d",&res.triangles,&res.vertices);

  res.vertex = (M_point*)malloc(res.vertices*sizeof(M_point));

  if (flags & MF_ADJ)
    res.adj = (M_int_triple*)malloc(res.triangles*sizeof(M_int_triple));

  if (flags & MF_DEG)
    {
      res.deg = (int*)malloc(res.vertices*sizeof(int));
      memset(res.deg,0,res.vertices*sizeof(int));
    }

  res.triangle = (M_int_triple*)malloc(res.triangles*sizeof(M_int_triple));
  for ( i=0; i<res.triangles; i++ )
    {
      fscanf(fp,"%d %d %d",&res.triangle[i][0],
	     &res.triangle[i][1],&res.triangle[i][2]);
      
      assert(res.triangle[i][0]>=0 && res.triangle[i][0]<res.vertices);
      assert(res.triangle[i][1]>=0 && res.triangle[i][1]<res.vertices);
      assert(res.triangle[i][2]>=0 && res.triangle[i][2]<res.vertices);

      if (flags & MF_DEG)
	{
	  res.deg[res.triangle[i][0]]++;
	  res.deg[res.triangle[i][1]]++;
	  res.deg[res.triangle[i][2]]++;
	}
    }

  if (flags & MF_INC)
    {
      res.inc = (M_pair**)malloc(res.vertices*sizeof(M_pair*));

      for ( i=0; i<res.vertices; i++ )
	res.inc[i] = res.deg ? (M_pair*)malloc(res.deg[i]*sizeof(M_pair)) : NULL;

      w = (int*)malloc(res.vertices*sizeof(int));
      memset(w,0,res.vertices*sizeof(int));

      for ( i=0; i<res.triangles; i++ )
	{
	  res.inc[res.triangle[i][0]][w[res.triangle[i][0]]].t = i;
	  res.inc[res.triangle[i][0]][w[res.triangle[i][0]]++].ix = 0;
	  res.inc[res.triangle[i][1]][w[res.triangle[i][1]]].t = i;
	  res.inc[res.triangle[i][1]][w[res.triangle[i][1]]++].ix = 1;
	  res.inc[res.triangle[i][2]][w[res.triangle[i][2]]].t = i;
	  res.inc[res.triangle[i][2]][w[res.triangle[i][2]]++].ix = 2;
	}

      for ( i=0; i<res.vertices; i++ )
	assert(w[i]==res.deg[i]);

      free(w);
    }

  for ( i=0; i<res.vertices; i++ )
    fscanf(fp,"%lf %lf %lf",&res.vertex[i][0],
	  &res.vertex[i][1],&res.vertex[i][2]);

  if (flags & (MF_ADJ | MF_EDGE_STAT))
    {
      q = (qdt*)malloc(3*res.triangles*sizeof(qdt));

      for ( j=i=0; i<res.triangles; i++ )
	{
	  q[j][0] = res.triangle[i][0];
	  q[j][1] = res.triangle[i][1];
	  q[j][2] = i;
	  q[j][3] = 2;
	  qd_normalize(&q[j++]);
	  
	  q[j][0] = res.triangle[i][1];
	  q[j][1] = res.triangle[i][2];
	  q[j][2] = i;
	  q[j][3] = 0;
	  qd_normalize(&q[j++]);

	  q[j][0] = res.triangle[i][2];
	  q[j][1] = res.triangle[i][0];
	  q[j][2] = i;
	  q[j][3] = 1;
	  qd_normalize(&q[j++]);
	}

      assert(j==3*res.triangles);

      qsort(q,j,sizeof(qdt),&qd_cmp);

      res.internal_edges = 0;
      res.external_edges = 0;

      for ( i=0; i<3*res.triangles; i+=di )
	{
	  if (i && q[i][0]==q[i-1][0] && q[i][1]==q[i-1][1])
	    {
	      printf("An edge with 3 incident triangles found - exiting...\n");
	      printf("triangle(%d %d) and triangle(%d %d)\n",q[i][2],q[i][3],
		     q[i-1][2],q[i-1][3]);
	      printf("(%lf %lf %lf) (%lf %lf %lf)\n",
		     res.vertex[res.triangle[q[i][2]][(q[i][3]+1)%3]][0],
		     res.vertex[res.triangle[q[i][2]][(q[i][3]+1)%3]][1],
		     res.vertex[res.triangle[q[i][2]][(q[i][3]+1)%3]][2],
		     res.vertex[res.triangle[q[i][2]][(q[i][3]+2)%3]][0],
		     res.vertex[res.triangle[q[i][2]][(q[i][3]+2)%3]][1],
		     res.vertex[res.triangle[q[i][2]][(q[i][3]+2)%3]][2]);
	      exit(1);
	    }
	  if ((i<3*res.triangles-1 && (q[i][0]!=q[i+1][0] || q[i][1]!=q[i+1][1])) 
	      || i==3*res.triangles-1) /* NEW!!! */
	    {
	      res.adj[q[i][2]][q[i][3]] = -1;
	      res.external_edges++;
	      di = 1;
	    }
	  else
	    {
	      res.adj[q[i][2]][q[i][3]] = q[i+1][2];
	      res.adj[q[i+1][2]][q[i+1][3]] = q[i][2];
	      res.internal_edges++;
	      di = 2;
	    }
	}
      
      free(q);
    }

  if (flags & MF_NORMAL)
    {
      res.normal = (M_point*)malloc(res.triangles*sizeof(M_point));

      for ( i=0; i<res.triangles; i++ )
	{
	  M_point v1,v2;

	  M_point_subtract(v1,res.vertex[res.triangle[i][0]],res.vertex[res.triangle[i][1]]);
	  M_point_subtract(v2,res.vertex[res.triangle[i][0]],res.vertex[res.triangle[i][2]]);
	  M_point_cross(res.normal[i],v1,v2);
	  M_point_normalize(res.normal[i]);
	}
    }

  return res;
}

/* -------------------------------- */

void M_write_mesh ( FILE *fp, M_mesh m )
{
  int i;

  fprintf(fp,"%d %d\n\n",m.triangles,m.vertices);
  
  for ( i=0; i<m.triangles; i++ )
    fprintf(fp,"%d %d %d\n",m.triangle[i][0],m.triangle[i][1],m.triangle[i][2]);

  fprintf(fp,"\n");

  for ( i=0; i<m.vertices; i++ )
    fprintf(fp,PFORMAT,m.vertex[i][0],m.vertex[i][1],m.vertex[i][2]);
}

/* -------------------------------- */

void M_quantize ( M_mesh *m, double qsize )
{
  int i;

  for ( i=0; i<m->vertices; i++ )
    {
      m->vertex[i][0] = floor(m->vertex[i][0]/qsize);
      m->vertex[i][1] = floor(m->vertex[i][1]/qsize);
      m->vertex[i][2] = floor(m->vertex[i][2]/qsize);
    }
}

/* ---------------------------------- */

double M_point_dist ( M_point p1, M_point p2 )
{
  M_point p;
  M_point_subtract(p,p1,p2);
  return sqrt(M_point_dot(p,p));
}
