

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

#include <iostream>
#include <fstream>
#include <cassert>
#include <program.h>

using namespace std;

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniformf ( const char *uname, GLfloat val )
{
  glUniform1f(glGetUniformLocation(p,uname),val);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform3fv ( const char *uname, GLfloat *val )
{
  glUniform3fv(glGetUniformLocation(p,uname),1,val);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform4fv ( const char *uname, GLfloat *val )
{
  glUniform4fv(glGetUniformLocation(p,uname),1,val);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform3x3fv ( const char *uname, GLfloat *val )
{
  glUniformMatrix3fv(glGetUniformLocation(p,uname),1,GL_FALSE,val);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform4x4fv ( const char *uname, GLfloat *val )
{
  glUniformMatrix4fv(glGetUniformLocation(p,uname),1,GL_FALSE,val);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniformd ( const char *uname, GLdouble val )
{
  GLfloat valf = val;
  glUniform1f(glGetUniformLocation(p,uname),valf);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform3dv ( const char *uname, GLdouble *val )
{
  GLfloat valf[3] = {val[0], val[1], val[2] };
  glUniform3fv(glGetUniformLocation(p,uname),1,valf);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform4dv ( const char *uname, GLdouble *val )
{
  GLfloat valf[4] = {val[0], val[1], val[2], val[3] };  
  glUniform4fv(glGetUniformLocation(p,uname),1,valf);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform3x3dv ( const char *uname, GLdouble *val )
{
  GLfloat valf[9] = { val[0], val[1], val[2],
		      val[3], val[4], val[5],
		      val[6], val[7], val[8] };
  glUniformMatrix3fv(glGetUniformLocation(p,uname),1,GL_FALSE,valf);
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetUniform4x4dv ( const char *uname, GLdouble *val )
{
  GLfloat valf[16] = {  val[0], val[1], val[2], val[3], 
			val[4], val[5], val[6], val[7], 
			val[8], val[9], val[10],val[11],
                        val[12],val[13],val[14],val[15] };
  glUniformMatrix4fv(glGetUniformLocation(p,uname),1,GL_FALSE,valf);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::UseMe()
{
  glUseProgram(p);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::StopMe()
{
  glUseProgram(0);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

GLchar *GLSLprogram::read_file ( const char *name )
{
  int size = 0;
  {  
    ifstream ifs(name);
    assert(ifs);
    do {
      size++;
      char c;
      c = ifs.get();
    }
    while(!ifs.eof());
  }
  char *res = new char[size];
  ifstream ifs(name);
  ifs.read(res,size-1);
  res[size-1] = 0;
  return res;
}

/* ----------------------------------------------------------- */

void GLSLprogram::_PrintInfoLog ( GLuint obj )
{
  int infologLength = 0;
  int charsWritten  = 0;
  char *infoLog;
  
  glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

  if (infologLength > 0)
    {
      infoLog = new char[infologLength];
      glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
      cout << infoLog << endl;
      delete[] infoLog;
    }
}

/* ----------------------------------------------------------- */

void GLSLprogram::PrintInfoLog( )
{
  int infologLength = 0;
  int charsWritten  = 0;
  char *infoLog;
  
  glGetProgramiv(p, GL_INFO_LOG_LENGTH,&infologLength);
  
  if (infologLength > 0)
    {
      infoLog = new char[infologLength];
      glGetProgramInfoLog(p, infologLength, &charsWritten, infoLog);
      cout << infoLog << endl;
      delete[] infoLog;
    }
}

/* ----------------------------------------------------------- */

GLSLprogram::GLSLprogram( const char *vs, const char *gs, const char *fs, unsigned char type )
{
  // create program... this one is for shadow polygons
  GLint vid = glCreateShader(GL_VERTEX_SHADER);
  GLint fid = glCreateShader(GL_FRAGMENT_SHADER);
  GLint gid;
  if (gs) gid = glCreateShader(GL_GEOMETRY_SHADER_EXT);
  
  const GLchar *gsh = NULL;
  const GLchar *vsh = (type=='f') ? read_file(vs) : vs;
  const GLchar *fsh = (type=='f') ? read_file(fs) : fs;
  if (gs) gsh = (type=='f') ? read_file(gs) : gs;
  glShaderSource(vid,1,&vsh,NULL);
  glShaderSource(fid,1,&fsh,NULL);
  if (gs) glShaderSource(gid,1,&gsh,NULL);
  glCompileShader(vid);
  glCompileShader(fid);
  if (gs) glCompileShader(gid);
  _PrintInfoLog(vid);
  _PrintInfoLog(fid);
  if (gs) _PrintInfoLog(gid);
  p = glCreateProgram();
  glAttachShader(p,vid);
  glAttachShader(p,fid);
  if (gs) glAttachShader(p,gid);
  glLinkProgram(p);
  if (gsh && type=='f') delete[] gsh;
  if (vsh && type=='f') delete[] vsh;
  if (fsh && type=='f') delete[] fsh;
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void GLSLprogram::SetTexture ( const char *uname, GLint texid, int texix )
{
  glUniform1i(glGetUniformLocation(p,uname),texix);
  glActiveTexture(GL_TEXTURE0+texix);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,texid);
  assert(!glGetError());
}

/* ----------------------------------------------------------- */
