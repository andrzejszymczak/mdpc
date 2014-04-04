
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

#ifndef __PROGRAM_H
#define __PROGRAM_H

#include <global.h>
#include <GL/glew.h>
#include <GL/gl.h>

/* ----------------------------------------------------------- */

class GLSLprogram {

  static void _PrintInfoLog ( GLuint obj );
  static GLchar *read_file ( const char *name );
  GLint p;

 public:

  // type='f': read from file; otherwise, just assume shader is in string
  GLSLprogram ( const char *vs, const char *gs, const char *fs, unsigned char type = 'f' );
  ~GLSLprogram();

  void SetUniform3fv ( const char *uname, GLfloat *val );
  void SetUniform4fv ( const char *uname, GLfloat *val );
  void SetUniform3x3fv ( const char *uname, GLfloat *val );
  void SetUniform4x4fv ( const char *uname, GLfloat *val );
  void SetUniformf ( const char *uname, GLfloat val );

  void SetUniform3dv ( const char *uname, GLdouble *val );
  void SetUniform4dv ( const char *uname, GLdouble *val );
  void SetUniform3x3dv ( const char *uname, GLdouble *val );
  void SetUniform4x4dv ( const char *uname, GLdouble *val );
  void SetUniformd ( const char *uname, GLdouble val );
    
  void SetTexture ( const char *uname, GLint texid, int texix );

  void UseMe();
  void StopMe();
  void PrintInfoLog();
};

/* ----------------------------------------------------------- */

#endif
