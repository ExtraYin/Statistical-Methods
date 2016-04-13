# GMC
The generalized measures of correlation (GMC) is given by the following formula:
\left\{ GMC(Y|X), GMC(X|Y) \roght\} = \left\{ 1-frac{E[\left\{Y-E(Y|X)\right\}^2], 1-frac{E[\left\{X-E(X|Y)\right\}^2]}{Var(X)} \right}

GMC.R calculates this GMC value of given (X,Y).  -- This file is not written by me!

GMC.cpp does the same as GMC.R but much faster(because it is written in c++ :) ). Please note that you cannot compile GMC.cpp as a usual c++ file. Instead, you should call this function in R with: library(Rcpp);sourceCpp("GMC.cpp");
