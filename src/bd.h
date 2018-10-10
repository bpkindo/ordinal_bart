#ifndef GUARD_bd_h
#define GUARD_bd_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "info.h"
#include "tree.h"

bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi);

#endif
