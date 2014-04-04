
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

#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <pcvfdisplay.h>
#include <program.h>
#include <trackball.h>
#include <primset.h>
#include <pcvf.h>

#include <cmath>
#include <cstdio>

#include <allshaders.cpp>

using namespace std;
using namespace glm;

/* ----------------------------------------------------------- */

GLint wid;
GLint width = 1024;
GLint height = 1024;
GLint initial_width = width;
GLint initial_height = height;

trackball *tb;
pcvfdisplay *vf;
primset *ms;
primset *sep;

GLuint depth_texture = 0;
GLuint vf_texture = 0;
GLuint noise_texture = 0;
GLuint fbo = 0;

/* ----------------------------------------------------------- */

// programs

GLSLprogram *p;
GLSLprogram *q;
GLSLprogram *r;

/* ----------------------------------------------------------- */

void setup_fbo()
{
  if (vf_texture==0)
    glGenTextures(1,&vf_texture);
  glBindTexture(GL_TEXTURE_2D,vf_texture);
  glTexImage2D(GL_TEXTURE_2D,0,4,width,height,0,GL_RGB,GL_FLOAT,NULL);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
  
  if (depth_texture==0)
    glGenTextures(1,&depth_texture);
  glBindTexture(GL_TEXTURE_2D,depth_texture);
  glTexImage2D(GL_TEXTURE_2D,0,GL_DEPTH_COMPONENT16,width,height,0,GL_DEPTH_COMPONENT,GL_FLOAT,NULL);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);

  if (fbo==0)
    glGenFramebuffers(1,&fbo);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,vf_texture,0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_DEPTH_ATTACHMENT_EXT,GL_TEXTURE_2D,depth_texture,0);
  assert(!glGetError());

  GLenum error;
  glGetIntegerv( GL_FRAMEBUFFER_BINDING_EXT, (GLint*)&fbo );

  // check the error status of this framebuffer */
  error = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);

  const char *fbName = "";
  // if error != GL_FRAMEBUFFER_COMPLETE_EXT, there's an error of some sort 
  {
    switch(error)
      {
      case GL_FRAMEBUFFER_COMPLETE_EXT:
	break;
      case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
	printf("Error!  %s missing a required image/buffer attachment!\n", fbName);
	break;
      case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
	printf("Error!  %s has no images/buffers attached!\n", fbName);
	break;
      case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
	printf("Error!  %s has mismatched image/buffer dimensions!\n", fbName);
	break;
      case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
	printf("Error!  %s's colorbuffer attachments have different types!\n", fbName);
	break;
      case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
	printf("Error!  %s trying to draw to non-attached color buffer!\n", fbName);
	break;
      case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
	printf("Error!  %s trying to read from a non-attached color buffer!\n", fbName);
	break;
      case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
	printf("Error!  %s format is not supported by current graphics card/driver!\n", fbName);
	break;
      default:
	printf("*UNKNOWN ERROR* reported from glCheckFramebufferStatusEXT() for %s!\n", fbName);
	break;
      }
  }
}

/* ----------------------------------------------------------- */

void setup_noise ( int w, int h )
{
  if (noise_texture==0)
    glGenTextures(1,&noise_texture);
  glBindTexture(GL_TEXTURE_2D,noise_texture);
  {
    float *arr = new float[w*h];
    for ( int i=0; i<w*h; i++ )
      arr[i] = drand48();
    glTexImage2D(GL_TEXTURE_2D,0,1,w,h,0,GL_LUMINANCE,GL_FLOAT,arr);
    delete[] arr;
  }
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
  assert(!glGetError());
}

/* ----------------------------------------------------------- */

void render_quad()
{
  static GLfloat vq[] = { 1.0, 1.0,
			  1.0,  -1.0,
			  -1.0,  1.0,
			  -1.0,  -1.0 };

  glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,vq);
  glEnableVertexAttribArray(0);
  glDrawArrays(GL_TRIANGLE_STRIP,0,8);
  glDisableVertexAttribArray(0);
}

/* ----------------------------------------------------------- */

void display()
{
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_POINT_SMOOTH);
  glDisable(GL_CULL_FACE);

  glPointSize(4.0);
  glutSetWindow(wid);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo);
  glClearColor(1.0,1.0,1.0,1.0);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  p->UseMe();
  p->SetUniform4x4dv("P",vf->pmpointer());
  p->SetUniform3x3dv("NM",tb->nmpointer());
  dmat4 mview = vf->ftr()*tb->mv()*vf->normt();
  p->SetUniform4x4dv("MV",&mview[0][0]);
  glViewport(0,0,initial_width,initial_height);
  vf->render();
  p->StopMe();

  glViewport(0,0,width,height);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  q->UseMe();
  q->SetTexture("noise",noise_texture,0);
  q->SetTexture("depth",depth_texture,1);
  q->SetTexture("vf",vf_texture,2);
  render_quad();

  r->UseMe();
  r->SetUniformf("depth_offset",-0.002);
  r->SetUniform4x4dv("P",vf->pmpointer());
  r->SetUniform4x4dv("MV",&mview[0][0]);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  if (ms) ms->render();
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  if (ms) ms->render();
  r->SetUniformf("depth_offset",-0.001);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  if (sep) sep->render();
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  if (sep) sep->render();

  glutSwapBuffers();
}

/* ----------------------------------------------------------- */

static bool zooming = false;
static int lasty;
static double ifov;

/* ----------------------------------------------------------- */

// function called when a mouse button is pressed or released

GLvoid mouse_button(GLint btn, GLint state, GLint x, GLint y)
{
  switch (btn)
    {
    case GLUT_LEFT_BUTTON:
      switch (state)
	{
	case GLUT_DOWN:
	  tb->mousedown(x,y);
	  break;
	case GLUT_UP:
	  tb->mouseup(x,y);
	  glutPostRedisplay();
	  break;
	}
      break;

    case GLUT_MIDDLE_BUTTON:
      switch (state)
	{
	case GLUT_DOWN:
	  zooming = true;
	  lasty = y;
	  ifov = vf->getfov();
	  break;
	case GLUT_UP:
	  zooming = false;
	  break;
	}
      break;
    }
}

/* --------------------------------------------- */

// function called when mouse is moving with a button down

GLvoid button_motion(GLint x, GLint y)
{
  if (tb->isactive())
    {
      tb->mousemove(x,y);
    }
  else
    if (zooming)
      {
	int dy = y-lasty;
	double nfov = ifov*exp(log(1.01)*dy);
	if (nfov>40)
	  nfov = 40;
	vf->setfov(nfov);
      }
  glutPostRedisplay(); // add display event to queue 
}

/* --------------------------------------------- */

/* handle keyboard events; here, just exit if ESC is hit */

GLvoid keyboard(GLubyte key, GLint x, GLint y)
{
  switch(key)
    {
    case 27:  /* ESC */
      exit(0);
      break;
      
    default:  
      break;
    }
}

/* ----------------------------------------------------------- */

/* handle resizing the glut window */

GLvoid reshape(GLint vpw, GLint vph)
{
  glutSetWindow(wid);
  width = vpw < vph ? vpw : vph;
  height = width;
  glViewport(0, 0, width, height);
  glutReshapeWindow(width, height);
  tb->resize(width,height);
  vf->resize(width,height);
  glutPostRedisplay();   // add display event to queue
}

/* ----------------------------------------------------------- */

void print_usage()
{
  cout << "Usage: " << endl;
  cout << " msvis [options] <VF> <MD> [<CONN>]" << endl;
  cout << "  VF: vector field (.t file) " << endl;
  cout << "  MD: Morse decomposition (.prs file)" << endl;
  cout << "  CONN: connection regions (.prs file)" << endl;
  cout << " options: " << endl;
  cout << "  -v: vector field is vertex based" << endl;
  cout << "  -f: vector is face based [default]" << endl;
}

/* ----------------------------------------------------------- */

int main ( int argc, char *argv[] )
{
  int i = 1;

  char type = 'f';

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
	      
	    case 'f':
	      if (argv[i][2])
		{
		  cout << "Unknown option: " << argv[i] << endl;
		  print_usage();
		  return 0;
		}
	      type = 'f';
	      i++;
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
  

  glutInit(&argc,argv);
  glutInitWindowSize(width,height);
  glutInitWindowPosition(10,10);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  wid = glutCreateWindow("Morse set visualization");    
   glewInit();
  if (glewIsSupported("GL_VERSION_4_0"))
    ;
  else 
    {
      cout << "OpenGL 4.0 not supported" << endl;;
      return 1;
    }

  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse_button);           // button press/release
  glutMotionFunc(button_motion);         // mouse motion w/ button down
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);

  // initialize trackball
  tb = new trackball(width,height);
  
  if (i>=argc)
    {
      print_usage();
      return 0;
    }

  // load vector field
  vf = new pcvfdisplay(argv[i],10,type);

  // initialize programn
  p = new GLSLprogram(shader1vsh,NULL,shader1fsh,'t');
  q = new GLSLprogram(shader2vsh,NULL,shader2fsh,'t');
  r = new GLSLprogram(shader3vsh,NULL,shader3fsh,'t');

  // initialize textures and fbos
  setup_noise(1024,1024);
  setup_fbo();

  // read primitive info
  if (argc>i+1)
    {
      ms = new primset(argv[i+1],2*vf->boxsize(),0.25);
    }
  else ms = NULL;
  if (argc>i+2) 
    {
      sep = new primset(argv[i+2],2*vf->boxsize(),0.25);
    }
  else sep = NULL;

  glutMainLoop();

  return 0;
}

/* ----------------------------------------------------------- */
