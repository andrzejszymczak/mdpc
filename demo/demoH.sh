
#!/bin/bash

ulimit -s 102400

echo ---------------------------------
echo A part of the ocean data http://svs.gsfc.nasa.gov/vis/a000000/a003800/a003827/
echo note that this vector field is face based not vertex based, hence no -v option
echo also, it is defined on a quad mesh -- file format is slightly different
echo ---------------------------------
echo ------ refinement level: 5, Morse sets only
echo ---------------------------------
echo Running: ../mdpc ../input/oceanf.t 5 o.prs 
../mdpc ../input/oceanf.t 5 o.prs 
echo Running: ../msvis ../input/oceanf.t o.prs 
echo Press ESC to exit
../msvis ../input/oceanf.t o.prs 
echo ---------------------------------
echo ------ refinement level: 10, Morse sets only
echo ---------------------------------
echo Running: ../mdpc ../input/oceanf.t 10 o.prs 
../mdpc ../input/oceanf.t 10 o.prs 
echo Running: ../msvis ../input/oceanf.t o.prs 
echo Press ESC to exit
../msvis ../input/oceanf.t o.prs 

