
#include "acb_theta.h"

void
acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec,
	slong mag_bits)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t mat;
    arb_t test;

    fmpz_mat_init(mat, 2*g, 2*g);
    arb_init(test);

    acb_siegel_randtest(tau, state, prec, mag_bits);
    acb_siegel_reduce(tau, mat, tau, prec);
    
    fmpz_mat_clear(mat);
    arb_clear(test);
}
