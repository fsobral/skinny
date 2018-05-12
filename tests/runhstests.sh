#!/bin/bash

PROBS="6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26 27 28 \
29 30 31 32 33 34 35 36 37 39 40 41 42 43 44 46 47 48 49 50 51 52 53 \
54 55 56 57 58 59 60 61 62 63 64 65 66 68 69 70 71 72 73 74 75 76 77 \
78 79 80 81 83 84 86 87 88 89 90 91 92 93 95 96 97 98 99 100 101 102 \
103 104 105 106 107 108 109 111 112 113 114 116 117 118 119"

SUBPROBS="16 17 18 30 33 54 65 70 75 84 86 87 96 97 98 101 102 103 \
104 105 106 107 108 109 111 112 113 114 116 117 118 119"

if (test $# = 2 ) then
    
    EXECNAME=$1;
    
    ulimit -St 1200 # 20 min of CPU time per problem
    
    rm analysis-all -f
    
    for i in $PROBS; do
#        echo $2 $3 $i | ./$EXECNAME
        echo $2 $i | ./$EXECNAME
        cat runhs.out >> analysis-all

        # Only for ALGENCAN's tests
        #
        #awk '{ printf "%15d %9.2f\n",$20,$25 }' algencan-tabline.out >> \
        #    analysis-all

        rm runhs.out -f
    done;

else

    echo -e "\n\nUSO: ./runhstests.sh <NOMEEXEC>\n\n";

fi
