
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

#if(!defined(__VEC3D_H))
#define __VEC3D_H

#include <global.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <assert.h>

/* ---------------------------------------------------------------------------- */

template<class T>
class vec3d {
  T x[3];

 public:
  
  vec3d ( T a, T b, T c );
  vec3d ( T a );   // diagonal vector [a a a]
  vec3d();
  vec3d ( const vec3d<T> &v );
  vec3d & operator= ( vec3d<T> v );
  ~vec3d();
  
  template<class U>
    vec3d ( vec3d<U> w )
    {
      x[0] = w[0];
      x[1] = w[1];
      x[2] = w[2];
    }

  vec3d & operator+= ( vec3d<T> v );
  vec3d & operator-= ( vec3d<T> v );
  vec3d & operator^= ( vec3d<T> v );   // cross-prouct
  vec3d & operator*= ( T s );          // multiply by scalar
  vec3d & operator|= ( vec3d<T> v );
  vec3d & operator&= ( vec3d<T> v );

  T norm();
  T norm2();
  T min();
  T max();

  void normalize();
  void writeto ( std::ostream &o );

  T& operator[] ( int i );

  T* pointer();
};

/* ---------------------------------------------------------------------------- */

typedef vec3d<float> vec3df;
typedef vec3d<double> vec3dd;
typedef vec3d<int> vec3di;
typedef vec3d<unsigned short int> vec3dus;

/* ---------------------------------------------------------------------------- */

template<class T>
extern vec3d<T> operator+ ( vec3d<T> p, vec3d<T> q );

template<class T>
extern vec3d<T> operator- ( vec3d<T> p, vec3d<T> q );

template<class T>
extern vec3d<T> operator| ( vec3d<T> p, vec3d<T> q );

template<class T>
extern vec3d<T> operator& ( vec3d<T> p, vec3d<T> q );

template<class T>
extern vec3d<T> operator- ( vec3d<T> p );

template<class T>
extern vec3d<T> operator^ ( vec3d<T> p, vec3d<T> q );

template<class T>
extern vec3d<T> operator* ( T p, vec3d<T> q );

template<class T>
extern T operator* ( vec3d<T> p, vec3d<T> q );

template<class T>
std::ostream & operator<< ( std::ostream &o, vec3d<T> v );

template<class T>
vec3d<T> randomv();

template<class T>
extern bool operator== ( vec3d<T> p, vec3d<T> q );

template<class T>
extern bool operator!= ( vec3d<T> p, vec3d<T> q );

template<class T>
extern bool operator< (  vec3d<T> p, vec3d<T> q );

template<class T>
extern bool operator> (  vec3d<T> p, vec3d<T> q );

template<class T>
extern bool operator<= (  vec3d<T> p, vec3d<T> q );

template<class T>
extern bool operator>= (  vec3d<T> p, vec3d<T> q );

template<class T>
extern vec3d<T> average ( vec3d<T> u, vec3d<T> w );

/* ---------------------------------------------------------------------------- */
/* ------------------- IMPLEMENTATION ----------------------------------------- */
/* ---------------------------------------------------------------------------- */

template<class T>
T* vec3d<T>::pointer()
{
  return &x[0];
}

/* ---------------------------------------------------------------------------- */

template<class T>
T vec3d<T>::min()
{
  if (x[0]>x[1])
    return x[1]<x[2] ? x[1] : x[2];
  else
    return x[0]<x[2] ? x[0] : x[2];
}

/* ---------------------------------------------------------------------------- */


template<class T>
T vec3d<T>::max()
{
  if (x[0]<x[1])
    return x[1]>x[2] ? x[1] : x[2];
  else
    return x[0]>x[2] ? x[0] : x[2];
}


/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T>::~vec3d()
{
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T>::vec3d ( T a, T b, T c )
{
  x[0] = a;
  x[1] = b;
  x[2] = c;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T>::vec3d ( T a )
{
  x[0] = a;
  x[1] = a;
  x[2] = a;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T>::vec3d ( )
{
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T>::vec3d ( const vec3d<T> &v )
{
  x[0] = v.x[0];
  x[1] = v.x[1];
  x[2] = v.x[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> & vec3d<T>::operator= ( vec3d<T> v )
{
  x[0] = v.x[0];
  x[1] = v.x[1];
  x[2] = v.x[2];
  return *this;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> & vec3d<T>::operator+= ( vec3d<T> v )
{
  x[0] += v.x[0];
  x[1] += v.x[1];
  x[2] += v.x[2];
  return *this;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> & vec3d<T>::operator-= ( vec3d<T> v )
{
  x[0] -= v.x[0];
  x[1] -= v.x[1];
  x[2] -= v.x[2];
  return *this;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> & vec3d<T>::operator|= ( vec3d<T> v )
{
  if (x[0]<v.x[0]) x[0]=v.x[0];
  if (x[1]<v.x[1]) x[1]=v.x[1];
  if (x[2]<v.x[2]) x[2]=v.x[2];
  return *this;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> & vec3d<T>::operator&= ( vec3d<T> v )
{
  if (x[0]>v.x[0]) x[0]=v.x[0];
  if (x[1]>v.x[1]) x[1]=v.x[1];
  if (x[2]>v.x[2]) x[2]=v.x[2];
  return *this;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> & vec3d<T>::operator*= ( T s )
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return *this;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> & vec3d<T>::operator^= ( vec3d<T> v )
{
  T a = x[1]*v.x[2]-x[2]*v.x[1];
  T b = x[2]*v.x[0]-x[0]*v.x[2];
  T c = x[0]*v.x[1]-x[1]*v.x[0];
  x[0] = a;
  x[1] = b;
  x[2] = c;
  return *this;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> operator+ ( vec3d<T> p, vec3d<T> q )
{
  p+=q;
  return p;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> operator- ( vec3d<T> p, vec3d<T> q )
{
  p-=q;
  return p;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> operator- ( vec3d<T> p )
{
  p[0] = -p[0];
  p[1] = -p[1];
  p[2] = -p[2];
  return p;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> operator^ ( vec3d<T> p, vec3d<T> q )
{
  p^=q;
  return p;
}

/* ---------------------------------------------------------------------------- */

template<class T>
bool operator== ( vec3d<T> p, vec3d<T> q )
{
  return p[0]==q[0] && p[1]==q[1] && p[2]==q[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
bool operator> ( vec3d<T> p, vec3d<T> q )
{
  return p[0]>q[0] && p[1]>q[1] && p[2]>q[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
bool operator< ( vec3d<T> p, vec3d<T> q )
{
  return p[0]<q[0] && p[1]<q[1] && p[2]<q[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
bool operator>= ( vec3d<T> p, vec3d<T> q )
{
  return p[0]>=q[0] && p[1]>=q[1] && p[2]>=q[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
bool operator<= ( vec3d<T> p, vec3d<T> q )
{
  return p[0]<=q[0] && p[1]<=q[1] && p[2]<=q[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
bool operator!= ( vec3d<T> p, vec3d<T> q )
{
  return p[0]!=q[0] || p[1]!=q[1] || p[2]!=q[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> operator* ( T p, vec3d<T> q )
{
  q*=p;
  return q;
}

/* ---------------------------------------------------------------------------- */

template<class T>
T operator* ( vec3d<T> p, vec3d<T> q )
{
  return p[0]*q[0]+p[1]*q[1]+p[2]*q[2];
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> operator& ( vec3d<T> p, vec3d<T> q )
{
  vec3d<T> res = p;
  res &=q;
  return res;
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> operator| ( vec3d<T> p, vec3d<T> q )
{
  vec3d<T> res = p;
  res |=q;
  return res;
}

/* ---------------------------------------------------------------------------- */

template<class T>
std::ostream & operator<< ( std::ostream &o, vec3d<T> v )
{
  o << "[ " << v[0] << " ; " << v[1] << " ; " << v[2] << " ]";
  return o;
}

/* ---------------------------------------------------------------------------- */

template<class T>
T & vec3d<T>::operator[] ( int i )
{
  assert(i>=0 && i<3);
  return x[i];
}

/* ---------------------------------------------------------------------------- */

template<class T>
void vec3d<T>::normalize()
{
  T d = 1/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  x[0] *= d;
  x[1] *= d;
  x[2] *= d;
}

/* ---------------------------------------------------------------------------- */

template<class T>
T vec3d<T>::norm()
{
  return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

/* ---------------------------------------------------------------------------- */

template<class T>
T vec3d<T>::norm2()
{
  return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> randomv3()
{
  vec3d<T> res;

  do {
    res = vec3d<T>(1-2*drand48(),1-2*drand48(),1-2*drand48());
  }
  while (res*res<0.1);

  res.normalize();

  return res;
}

/* ---------------------------------------------------------------------------- */

template<class T>
void vec3d<T>::writeto ( std::ostream &o )
{
  o.write((char*)x,3*sizeof(T));
}

/* ---------------------------------------------------------------------------- */

template<class T>
vec3d<T> average ( vec3d<T> u, vec3d<T> w )
{
  return vec3d<T>((u[0]+w[0])/2,(u[1]+w[1])/2,(u[2]+w[2])/2);
}

/* ---------------------------------------------------------------------------- */

#endif
