
#include "acb_theta.h"

void
acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)
{
    slong g = arb_mat_nrows(tau);
    arb_mat_t re, im;
    slong k, j;
  
    arb_mat_init(re, g, g);
    arb_mat_init(im, g, g);

    for (k = 0; k < g; k++)
    {
	for (j = k; j < g; j++)
	{
	    arb_randtest_precise(arb_mat_entry(re, k, j), state, prec, mag_bits);
	    arb_set(arb_mat_entry(re, j, k), arb_mat_entry(re, k, j));
	}
    }

    arb_mat_randtest_sym_pos(im, state, prec, mag_bits);
    acb_mat_set_arb_arb(tau, re, im);

    arb_mat_clear(re);
    arb_mat_clear(im);
}
