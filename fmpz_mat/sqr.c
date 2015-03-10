#include "fmpz_mat.h"

void
fmpz_mat_sqr(fmpz_mat_t C, const fmpz_mat_t A)
{
    slong dim, n;

    n = A->r;

    dim = n;

    if (dim <= 12)
    {
        fmpz_mat_mul(C, A, A);    
    }
    else
    {
        slong ab;

        ab = fmpz_mat_max_bits(A);
        ab = FLINT_ABS(ab);

        if (5*(ab + ab) > dim * dim )
        {
            fmpz_mat_sqr_bodrato(C, A);
        }
        else
        {
            fmpz_mat_mul(C, A, A);
        }

    }
}
