#!/bin/bash

mv ../running.mp .
mpost running.mp
for i in `seq 1 $1`;
do echo -e "\\\\begin{figure} \\\\centering\n \\\\includegraphics{running.$i}\n \\\\end{figure}\n \\\\clearpage\n";
done > running.tex

echo -e "\\\\begin{figure} \\\\centering\n \\\\includegraphics{running.1000}\n \\\\end{figure}\n \\\\clearpage\n" >> running.tex;

latex figtest.tex
dvipdfm figtest.dvi 
