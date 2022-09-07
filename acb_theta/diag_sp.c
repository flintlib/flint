
#include "acb_theta.h"

void fmpz_mat_diag_sp(fmpz_mat_t mat, const fmpz_mat_t U)
{
    slong g = fmpz_mat_nrows(mat)/2;
    fmpz_mat_t D, zero;
    fmpz_t den;

    fmpz_mat_init(D, g, g);
    fmpz_mat_init(zero, g, g);
    fmpz_init(den);

    fmpz_mat_inv(D, den, U);
    fmpz_mat_transpose(D, D);
    if (!fmpz_is_one(den))
    {
	fmpz_neg(den, den);
	fmpz_mat_neg(D, D);
    }
    if (!fmpz_is_one(den))
    {
	flint_fprintf(stderr, "fmpz_mat_diag_sp: Error (not invertible)\n");
	flint_abort();
    }

    fmpz_mat_set_abcd(mat, U, zero, zero, D);

    fmpz_mat_clear(D);
    fmpz_mat_clear(zero);
    fmpz_clear(den);
}
