
#include "acb_theta.h"

void arb_mat_randtest_sym_pos(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits)
{
  slong g = arb_mat_nrows(r);
  arb_mat_t tp;

  arb_mat_init(tp, g, g);
  arb_mat_randtest_cho(r, state, prec, mag_bits);
  arb_mat_transpose(tp, r);
  arb_mat_mul(r, tp, r, prec);

  arb_mat_clear(tp);
}
