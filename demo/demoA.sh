
#!/bin/bash

clear

echo First, increase the stack size
echo Executing:  ulimit -s 102400
ulimit -s 102400

echo ---------------------------------
echo A small planar vector field...
echo Note that it is vertex based so -v option is used for mdpc and msvis
echo ---------------------------------
echo Refinement depth: 2 for Morse sets and 2 for connection regions
echo ---------------------------------
echo Executing:  "../mdpc -v -c 2 2 ../input/multicyclesV.t 2 o.prs o.conn o.dot"
echo Note that one could also run "../mdpc -c 2 2 ../input/multicyclesT.t 2 o.prs o.conn o.dot"
echo to get the same output
../mdpc -v -c 2 2 ../input/multicyclesV.t 2 o.prs o.conn o.dot
echo
echo ---------------------------------
echo Executing:  "../msvis -v ../input/multicyclesV.t o.prs o.conn"
echo press ESC to exit
echo ---------------------------------
../msvis -v ../input/multicyclesV.t o.prs o.conn
echo ---------------------------------
echo Using graphviz to draw the MCG
echo Executing:  "dot -Tpng o.dot > o.png"
echo ---------------------------------
dot -Tpng o.dot > o.png
echo Executing:  "eog o.png"
echo press ESC to exit
eog o.png

clear

echo ---------------------------------
echo Refinement depth: 5 for Morse sets and 5 for connection regions
echo ---------------------------------
echo Executing:  "../mdpc -v -c 5 5 ../input/multicyclesV.t 5 o.prs o.conn o.dot"
../mdpc -v -c 5 5 ../input/multicyclesV.t 5 o.prs o.conn o.dot
echo
echo ---------------------------------
echo Executing:  "../msvis -v ../input/multicyclesV.t o.prs o.conn"
echo press ESC to exit
echo ---------------------------------
../msvis -v ../input/multicyclesV.t o.prs o.conn
echo ---------------------------------
echo Using graphviz to draw the MCG
echo Executing:  "dot -Tpng o.dot > o.png"
echo ---------------------------------
dot -Tpng o.dot > o.png
echo Executing:  "eog o.png"
echo press ESC to exit
eog o.png

clear

echo ---------------------------------
echo Refinement depth: 10 for Morse sets and 10 for connection regions
echo ---------------------------------
echo Executing:  "../mdpc -v -c 10 10 ../input/multicyclesV.t 10 o.prs o.conn o.dot"
../mdpc -v -c 10 10 ../input/multicyclesV.t 10 o.prs o.conn o.dot
echo
echo ---------------------------------
echo Executing:  "../msvis -v ../input/multicyclesV.t o.prs o.conn"
echo press ESC to exit
echo ---------------------------------
../msvis -v ../input/multicyclesV.t o.prs o.conn
echo ---------------------------------
echo Using graphviz to draw the MCG
echo Executing:  "dot -Tpng o.dot > o.png"
echo ---------------------------------
dot -Tpng o.dot > o.png
echo Executing:  "eog o.png"
echo press ESC to exit
eog o.png



