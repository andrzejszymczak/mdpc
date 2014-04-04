#!/bin/bash

ulimit -s 102400

clear

echo Diesel engine: stable Morse decomposition
echo 5 refinement iterations
echo stability: 4.2
echo Running: ../mdpc -v -s 4.2 ../input/dieselEngineV.t 5 o.prs
../mdpc -v -s 4.2 ../input/dieselEngineV.t 5 o.prs
echo Running: ../msvis -v ../input/dieselEngineV.t o.prs
echo Press ESC to quit
../msvis -v ../input/dieselEngineV.t o.prs

clear

echo Cooling jacket: stable Morse decomposition
echo 5 refinement iterations
echo Max perturbation: 0.15
echo Running: ../mdpc -v -s 0.15 ../input/coolingJacketV.t 5 o.prs
../mdpc -v -s 0.15 ../input/coolingJacketV.t 5 o.prs
echo Running: ../msvis -v ../input/coolingJacketV.t o.prs
echo Press ESC to quit
../msvis -v ../input/coolingJacketV.t o.prs

clear

echo Cooling jacket: stable Morse decomposition
echo 5 refinement iterations
echo Max perturbation: 0.45
echo Running: ../mdpc -v -s 0.45 ../input/coolingJacketV.t 5 o.prs
../mdpc -v -s 0.45 ../input/coolingJacketV.t 5 o.prs
echo Running: ../msvis -v ../input/coolingJacketV.t o.prs
echo Press ESC to quit
../msvis -v ../input/coolingJacketV.t o.prs

clear

echo Cooling jacket: Morse decomposition
echo 5 refinement iterations
echo stability: 0.9
echo Running: ../mdpc -v -s 0.9 ../input/coolingJacketV.t 5 o.prs
../mdpc -v -s 0.9 ../input/coolingJacketV.t 5 o.prs
echo Running: ../msvis -v ../input/coolingJacketV.t o.prs
echo Press ESC to quit
../msvis -v ../input/coolingJacketV.t o.prs

clear

echo Diesel engine: envelope
echo 7 refinement iterations
echo
echo Running: ../mdpc -v -e 1.0 ../input/dieselEngineV.t 7 o.prs
echo 
echo Note that the number that follows the -e option is the weight 
echo " used for the hull."
echo "The vector field used is f+weight*(H-f), where H is the envelope."
echo Weight=1 means that the envelope itself is used.
../mdpc -v -e 1.0 ../input/dieselEngineV.t 7 o.prs
echo Running: ../msvis -v ../input/dieselEngineV.t o.prs
echo Press ESC to quit
../msvis -v ../input/dieselEngineV.t o.prs

clear

echo Gas engine: envelope of the original dataset
echo 7 refinement iterations
echo
echo Running: ../mdpc -v -e 1.0 ../input/gasEngineV.t 7 o.prs
../mdpc -v -e 1.0 ../input/gasEngineV.t 7 o.prs
echo Running: ../msvis -v -v ../input/dieselEngineV.t o.prs
echo Press ESC to quit
../msvis -v -v ../input/gasEngineV.t o.prs


clear


echo Gas engine: envelope of the original dataset
echo after one PL subdivision
echo 7 refinement iterations
echo
echo Subdividing.
echo Running: ../subd/vsubd ../input/gasEngineV.t gas2V.t
../subd/vfsubd ../input/gasEngineV.t gas2V.t
echo Running: ../mdpc -v -e 1.0 gas2V.t 7 o.prs
../mdpc -v -e 1.0 gas2V.t 7 o.prs
echo Running: ../msvis -v -v gas2V.t o.prs
echo Press ESC to quit
../msvis -v -v gas2V.t o.prs


clear


echo Gas engine: envelope of the original dataset
echo after two PL subdivisions
echo 7 refinement iterations
echo
echo Subdividing.
echo Running: ../subd/vsubd gas2V.t gas3V.t
../subd/vfsubd gas2V.t gas3V.t
echo Running: ../mdpc -v -e 1.0 gas3V.t 7 o.prs
../mdpc -v -e 1.0 gas3V.t 7 o.prs
echo Running: ../msvis -v -v gas3V.t o.prs
echo Press ESC to quit
../msvis -v -v gas3V.t o.prs



