
#include "acb_theta.h"


void arb_mat_randtest_cho(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits)
{
  slong g = arb_mat_nrows(r);
  slong k, j;
  
  arb_mat_zero(r);
  for (k = 0; k < g; k++)
    {
      arb_randtest_pos(arb_mat_entry(r, k, k), state, prec, mag_bits);
      for (j = k+1; j < g; j++)
	{
	  arb_randtest_precise(arb_mat_entry(r, k, j), state, prec, mag_bits);
	}
    }
}
