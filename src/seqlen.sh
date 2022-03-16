#/bin/bash 

A=`grep -v \> $1 | wc -c`
B=`grep -v \> $1 | wc -l`
C=$(($A-$B))
echo $C