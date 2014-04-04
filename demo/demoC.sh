
#!/bin/bash

clear

echo Increasing the stack size
echo Executing:  ulimit -s 102400
ulimit -s 102400

echo ---------------------------------
echo A vector field on a torus, defined on a coarse grid
echo ---------------------------------
echo refinement level: 2 for connecting regions, 2 for Morse sets
echo ---------------------------------
echo Executing:  "../mdpc -v -c 2 2 ../input/torus1V.t 2 o.prs o.conn o.dot"
echo ---------------------------------
../mdpc -v -c 2 2 ../input/torus1V.t 2 o.prs o.conn o.dot
echo ---------------------------------
echo Executing:  "../msvis -v ../input/torus1V.t o.prs o.conn"
echo ---------------------------------
echo press ESC to exit
../msvis -v ../input/torus1V.t o.prs o.conn
echo ---------------------------------
echo Executing:  "dot -Tpng o.dot > o.png"
echo ---------------------------------
dot -Tpng o.dot > o.png
echo ---------------------------------
echo Executing:  "dot o.png"
echo ---------------------------------
echo press ESC to exit
eog o.png


clear


echo ---------------------------------
echo A vector field on a torus, defined on a coarse grid
echo ---------------------------------
echo refinement level: 5 for connecting regions, 5 for Morse sets
echo ---------------------------------
echo Executing:  "../mdpc -v -c 5 5 ../input/torus1V.t 5 o.prs o.conn o.dot"
echo ---------------------------------
../mdpc -v -c 5 5 ../input/torus1V.t 5 o.prs o.conn o.dot
echo ---------------------------------
echo Executing:  "../msvis -v ../input/torus1V.t o.prs o.conn"
echo ---------------------------------
echo press ESC to exit
../msvis -v ../input/torus1V.t o.prs o.conn
echo ---------------------------------
echo Executing:  "dot -Tpng o.dot > o.png"
echo ---------------------------------
dot -Tpng o.dot > o.png
echo ---------------------------------
echo Executing:  "dot o.png"
echo ---------------------------------
echo press ESC to exit
eog o.png


clear

echo ---------------------------------
echo A vector field on a torus, defined on a coarse grid
echo ---------------------------------
echo refinement level: 10 for connecting regions, 10 for Morse sets
echo ---------------------------------
echo Executing:  "../mdpc -v -c 10 10 ../input/torus1V.t 10 o.prs o.conn o.dot"
echo ---------------------------------
../mdpc -v -c 10 10 ../input/torus1V.t 10 o.prs o.conn o.dot
echo ---------------------------------
echo Executing:  "../msvis -v ../input/torus1V.t o.prs o.conn"
echo ---------------------------------
echo press ESC to exit
../msvis -v ../input/torus1V.t o.prs o.conn
echo ---------------------------------
echo Executing:  "dot -Tpng o.dot > o.png"
echo ---------------------------------
dot -Tpng o.dot > o.png
echo ---------------------------------
echo Executing:  "dot o.png"
echo ---------------------------------
echo press ESC to exit
eog o.png


clear

