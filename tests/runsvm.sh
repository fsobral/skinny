#!/bin/bash

N=9
TOTAL=90 # N * 10

ulimit -St 3600 # 1h of CPU time per problem

for i in `seq 1 $N $TOTAL`; do
    echo $i
    awk '(NR == '$i'){print "'$2'"};(NR >= '$i' && NR < '$i'+'$N'){print $0}' $3 | $1
    cat runeng.out
done;
