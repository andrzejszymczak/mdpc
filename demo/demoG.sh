#!/bin/bash

ulimit -s 102400

echo Morse connection graph for the diesel engine
echo Max subdivision depth for Morse sets: 5
echo Min/max depth for connecting regions: 5 5
echo Executing: ../mdpc -v -c 5 5 ../input/dieselEngineV.t 5 o.prs o.conn o.dot
../mdpc -v -c 5 5 ../input/dieselEngineV.t 5 o.prs o.conn o.dot
echo Executing: ../msvis -v ../input/dieselEngineV.t o.prs o.conn
echo Press ESC to exit
../msvis -v ../input/dieselEngineV.t o.prs o.conn
echo Executing: dot -Tpng o.dot > o.png
dot -Tpng o.dot > o.png
echo Executing: eog o.png
echo Press ESC to exit
eog o.png

clear

echo Morse connection graph for the diesel engine
echo Max subdivision depth for Morse sets: 10
echo Min/max depth for connecting regions: 10 10
echo Executing: ../mdpc -v -c 10 10 ../input/dieselEngineV.t 10 o.prs o.conn o.dot
../mdpc -v -c 10 10 ../input/dieselEngineV.t 10 o.prs o.conn o.dot
echo Executing: ../msvis -v ../input/dieselEngineV.t o.prs o.conn
echo Press ESC to exit
../msvis -v ../input/dieselEngineV.t o.prs o.conn
echo Executing: dot -Tpng o.dot > o.png
dot -Tpng o.dot > o.png
echo Executing: eog o.png
echo Press ESC to exit
eog o.png

clear

echo Morse connection graph for the gas engine
echo Max subdivision depth for Morse sets: 5
echo Min/max depth for connecting regions: 5 5
echo Executing: ../mdpc -v -c 5 5 ../input/gasEngineV.t 5 o.prs o.conn o.dot
../mdpc -v -c 5 5 ../input/gasEngineV.t 5 o.prs o.conn o.dot
echo Executing: ../msvis ../input/gasEngineV.t o.prs o.conn
echo Press ESC to exit
../msvis -v ../input/gasEngineV.t o.prs o.conn
echo Executing: dot -Tpng o.dot > o.png
dot -Tpng o.dot > o.png
echo Executing: eog o.png
echo Press ESC to exit
eog o.png

clear

echo Morse connection graph for the gas engine
echo Max subdivision depth for Morse sets: 9
echo Min/max depth for connecting regions: 10 12
echo Executing: ../mdpc -v -c 10 12 ../input/gasEngineV.t 9 o.prs o.conn o.dot
../mdpc -v -c 10 12 ../input/gasEngineV.t 9 o.prs o.conn o.dot
echo Executing ../msvis -v ../input/gasEngineV.t o.prs o.conn
echo Press ESC to exit
../msvis -v ../input/gasEngineV.t o.prs o.conn
echo Executing: dot -Tpng o.dot > o.png
dot -Tpng o.dot > o.png
echo Executing: eog o.png
echo Press ESC to exit
eog o.png

clear

echo Morse connection graph for the cooling jacket
echo Max subdivision depth for Morse sets: 10
echo Min/max depth for connecting regions: 10 14
echo Executing: ../mdpc -v -c 10 14 ../input/coolingJacketV.t 10 o.prs o.conn o.dot
../mdpc -v -c 10 14 ../input/coolingJacketV.t 10 o.prs o.conn o.dot
echo Executing ../msvis -v ../input/coolingJacketV.t o.prs o.conn
echo Press ESC to exit
../msvis -v ../input/coolingJacketV.t o.prs o.conn
echo Executing: dot -Tpng o.dot > o.png
dot -Tpng o.dot > o.png
echo Executing: gimp o.png
echo eog may have problems because the width of the image is very large
gimp o.png


