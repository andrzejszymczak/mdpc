
#!/bin/bash

#!/bin/bash

clear

echo Increasing the stack size
echo Executing:  ulimit -s 102400
ulimit -s 102400

echo ---------------------------------
echo A vector field on the bunny
echo 2 iterations of PL refinement is used to increase the resolution
echo ---------------------------------
echo Subdividing the input...
echo Ecexuting ../subd/vfsubd ../input/bunny2V.t level1.t
echo Executing ../subd/vfsubd level1.t level2.t
echo ---------------------------------
../subd/vfsubd ../input/bunny2V.t level1.t
../subd/vfsubd level1.t level2.t


echo refinement level: 2 for connecting regions, 2 for Morse sets
echo ---------------------------------
echo Executing:  "../mdpc -v -c 2 2 level2.t 2 o.prs o.conn o.dot"
echo ---------------------------------
../mdpc -v -c 2 2 level2.t 2 o.prs o.conn o.dot
echo ---------------------------------
echo Executing:  "../msvis -v level2.t o.prs o.conn"
echo ---------------------------------
echo press ESC to exit
../msvis -v level2.t o.prs o.conn
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
echo refinement level: 5 for connecting regions, 5 for Morse sets
echo ---------------------------------
echo Executing:  "../mdpc -v -c 5 5 level2.t 5 o.prs o.conn o.dot"
echo ---------------------------------
../mdpc -v -c 5 5 level2.t 5 o.prs o.conn o.dot
echo ---------------------------------
echo Executing:  "../msvis -v level2.t o.prs o.conn"
echo ---------------------------------
echo press ESC to exit
../msvis -v level2.t o.prs o.conn
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
echo refinement level: 12 for connecting regions, 12 for Morse sets
echo ---------------------------------
echo Executing:  "../mdpc -v -c 12 12 level2.t 12 o.prs o.conn o.dot"
echo ---------------------------------
../mdpc -v -c 12 12 level2.t 12 o.prs o.conn o.dot
echo ---------------------------------
echo Executing:  "../msvis -v level2.t o.prs o.conn"
echo ---------------------------------
echo press ESC to exit
../msvis -v level2.t o.prs o.conn
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

