
#include "acb_theta.h"

int arb_mat_is_nonsymmetric(const arb_mat_t m)
{
  arb_mat_t mt;
  slong nrows = arb_mat_nrows(m);
  int res;
  
  if (nrows != arb_mat_ncols(m)) return 1;

  arb_mat_init(mt, nrows, nrows);
  arb_mat_transpose(mt, m);
  res = !arb_mat_overlaps(mt, m);
  arb_mat_clear(mt);

  return res;
}
    
