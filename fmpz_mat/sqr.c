#include "fmpz_mat.h"

#define E fmpz_mat_entry

void
fmpz_mat_sqr(fmpz_mat_t B, const fmpz_mat_t A)
{
        slong n = A->r;
        slong dim = n;

        slong ab, bits;

        ab = fmpz_mat_max_bits(A);
        ab = FLINT_ABS(ab);
        bits = ab + FLINT_BIT_COUNT(n) + 1;

        if (bits < 100)
            fmpz_mat_sqr_classical(B, A);
        else if (bits >= 100 && bits < 800)
        {
            if(dim < 15 || dim > 125)
                fmpz_mat_sqr_classical(B, A);
            else
                fmpz_mat_sqr_bodrato(B, A);
        }
        else
        {
            if(dim < 125)
                fmpz_mat_sqr_bodrato(B, A);
            else
                fmpz_mat_sqr_classical(B, A);
        }
}
