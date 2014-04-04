
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

const char shader1fsh[] = 
"#version 400 \n\
\n\
in vec4 wcoord;\n\
in vec4 wcoord_tip;\n\
in float shade;\n\
\n\
out vec4 color;\n\
\n\
void main()\n\
{\n\
  vec2 crp = (1.0/wcoord.w)*wcoord.xy;\n\
  vec2 ctp = (1.0/wcoord_tip.w)*wcoord_tip.xy;\n\
  vec2 vf = ctp-crp;\n\
  vf = (1.0/(length(vf)+0.0001))*vf;\n\
\n\
  color.rg = vec2(0.5,0.5)+0.5*vf;\n\
  color.b = shade;\n\
}\n\
";

const char shader1vsh[] =
"#version 400\n\
\n\
layout (location=0) in vec3 coord;\n\
layout (location=1) in vec3 normal;\n\
layout (location=2) in vec3 vf;\n\
\n\
uniform mat4 MV;\n\
uniform mat3 NM;\n\
uniform mat4 P;\n\
\n\
const float factor = 0.01;\n\
\n\
out vec4 wcoord;\n\
out vec4 wcoord_tip;\n\
out float shade;\n\
\n\
void main()\n\
{\n\
  vec4 b;\n\
  vec3 wnormal;\n\
\n\
  b.xyz = coord + factor/(length(vf)+0.001)*vf; \n\
  b.w = 1.0;\n\
\n\
  wcoord = MV*vec4(coord,1.0);\n\
  wcoord_tip = P*MV*b;\n\
\n\
  wnormal = NM*normal;\n\
  shade = -1.0/(length(wnormal)*length(wcoord.xyz))*dot(wnormal,wcoord.xyz);\n\
\n\
  wcoord = P*wcoord;\n\
\n\
  gl_Position = wcoord;\n\
}\n\
";


const char shader2fsh[] = 
"#version 400\n\
\n\
uniform sampler2D vf;\n\
uniform sampler2D depth;\n\
uniform sampler2D noise;\n\
\n\
in vec2 tcoord;\n\
\n\
out vec4 color;\n\
\n\
const int N = 64;\n\
const float scale = 0.001;\n\
const float eps = 0.1;\n\
\n\
float f ( float s, float t )\n\
{\n\
  return 0.6*s+0.4*t;\n\
}\n\
\n\
void main()\n\
{\n\
  vec2 p,dp,pnew;\n\
  int i;\n\
  float total = 0;\n\
  int j = 0;\n\
  float shade;\n\
  if (texture2D(depth,tcoord)==1.0) discard;\n\
  p = tcoord;\n\
  for ( i=0; i<N; i++ )\n\
    {\n\
      dp = scale*(texture2D(vf,p).rg-vec2(0.5,0.5));\n\
      pnew = p+dp;\n\
      if (abs(texture2D(depth,pnew).r-texture2D(depth,p).r)>eps)\n\
	break;\n\
      p = pnew;\n\
      total += texture2D(noise,p).r;\n\
      j++;\n\
    }\n\
  p = tcoord;\n\
  for ( i=0; i<=N; i++ )\n\
    {\n\
      dp = scale*(texture2D(vf,p).rg-vec2(0.5,0.5));\n\
      total += texture2D(noise,p).r;\n\
      j++;\n\
      pnew = p+dp;\n\
      if (abs(texture2D(depth,pnew).r-texture2D(depth,p).r)>eps)\n\
	break;\n\
      p = pnew;\n\
    }\n\
  shade = texture2D(vf,tcoord).b;\n\
  total = total/float(j);\n\
  color.rgba = vec4(f(shade,total),f(shade,total),f(shade,total),1);\n\
  gl_FragDepth = texture2D(depth,tcoord).r;\n\
}\n\
";




const char shader2vsh[] = 
"#version 400\n\
\n\
layout (location=0) in vec2 coord;\n\
\n\
out vec2 tcoord;\n\
\n\
void main()\n\
{\n\
  tcoord = 0.5*(coord+vec2(1.0,1.0));\n\
  gl_Position = vec4(coord,0.0,1.0);\n\
}\n\
";


const char shader3fsh[] = 
"#version 400\n\
\n\
in vec4 wcoord;\n\
in vec3 clr;\n\
\n\
uniform float depth_offset;\n\
\n\
out vec4 fcolor;\n\
\n\
void main()\n\
{\n\
  fcolor.rgb = clr.rgb;\n\
  gl_FragDepth = gl_FragCoord.z+depth_offset;\n\
}\n\
";

const char shader3vsh[] =
"#version 400\n\
\n\
layout (location=0) in vec3 coord;\n\
layout (location=1) in vec3 color;\n\
\n\
uniform mat4 MV;\n\
uniform mat4 P;\n\
\n\
out vec3 clr;\n\
\n\
void main()\n\
{\n\
  clr = color;\n\
  gl_Position = P*MV*vec4(coord,1);\n\
}\n\
";
