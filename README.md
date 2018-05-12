Skinny - An algorithm for constrained derivative-free problems
==============================================================

`Skinny` is a derivative-free algorithm based on pseudo-random direct
searches for solving constrained problems. The main feature of this
algorithm is the ability to handle equality constraints. The objective
function must be defined outside the feasible region, as the algorithm
works by *enlarging* the feasible set and iteratively reducing it.

This version was used to produce the numerical results presented in

  - J. M. Martínez and F. N. C. Sobral, "Constrained Derivative-Free
    Optimization on Thin Domains", Journal of Global Optimization,
    56(3), pp. 1217-1232, 2013. [Accepted
    version](http://link.springer.com/10.1007/s10898-012-9944-x)
    [Submitted
    version](http://www.optimization-online.org/DB_FILE/2011/08/3139.pdf).

Instructions for building and extending `Skinny` can be found in the
[author's
web page](http://fsobral.github.io/skinny/index.html). Complete results
can be downloaded [here](doc/tables.pdf).

The present branch is a developing branch in order to adapt `Skinny` to work with [Algencan version
3.1](https://www.ime.usp.br/~egbirgin/tango/). We are also including all test problems created during the development of the algorithm.
