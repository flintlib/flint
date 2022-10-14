
#include "acb_theta.h"

static void get_fmpz_mat(fmpz_mat_t N, const arb_mat_t M, slong prec)
{
    slong j, k;
    
    for (j = 0; j < arb_mat_nrows(M); j++)
    {
	for (k = 0; k < arb_mat_ncols(M); k++)
	{
	    arf_get_fmpz_fixed_si(fmpz_mat_entry(N, j, k),
		    arb_midref(arb_mat_entry(M, j, k)),
		    -prec);	    
	}
    }
}

static void fmpz_mat_reduce(fmpz_mat_t N, fmpz_mat_t U)
{
    fmpz_lll_t fl;

    /* Default Flint LLL values, except Gram */
    fmpz_lll_context_init(fl, 0.99, 0.51, GRAM, EXACT);
    fmpz_mat_one(U);

    fmpz_lll(N, U, fl);    
}

void arb_mat_reduce(fmpz_mat_t U, const arb_mat_t M, slong prec)
{
    fmpz_mat_t N;
    slong g = acb_mat_nrows(M);

    
    fmpz_mat_one(U);
    
    /* Only proceed when M has finite entries */
    if (!arb_mat_is_finite(M))
    {
	return;
    }

    fmpz_mat_init(N, g, g);
    
    get_fmpz_mat(N, M, prec);
    fmpz_mat_reduce(N, U);

    fmpz_mat_clear(N);
}
