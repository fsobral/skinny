#!/bin/bash

PROBS="1 2 3 4 5 6 7"

if (test $# = 3 ) then
    
    EXECNAME=$1;
    
    ulimit -St 3600 # 1 hour of CPU time per problem
    
    rm analysis-all -f
    
    for i in $PROBS; do
        echo -e "$2 $3\n $i" | ./$EXECNAME
#        echo -e "$i\n" | ./$EXECNAME
        cat runeng.out >> analysis-all
        rm runeng.out -f
    done;

else

    echo -e "\n\nUSO: ./rundftests.sh <NOMEEXEC>\n\n";

fi
