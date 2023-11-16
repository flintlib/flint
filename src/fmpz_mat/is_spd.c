/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "arb_mat.h"

static int
fmpz_mat_is_spd_arb(const fmpz_mat_t A)
{
    slong d = fmpz_mat_nrows(A);
    slong prec;
    slong maxprec = 32;
    arb_mat_t M, L;
    slong j, k;
    int res = 0;

    arb_mat_init(M, d, d);
    arb_mat_init(L, d, d);

    for (j = 0; j < d; j++)
    {
        for (k = 0; k <= j; k++)
        {
            maxprec = FLINT_MAX(maxprec, 32 + fmpz_bits(fmpz_mat_entry(A, j, k)));
        }
    }
    arb_mat_set_fmpz_mat(M, A);

    for (prec = 32; (prec < 4 * maxprec) && !res; prec *= 2)
    {
        res = arb_mat_ldl(L, M, prec);
    }

    arb_mat_clear(M);
    arb_mat_clear(L);
    return res;
}

static int
fmpz_mat_is_spd_charpoly(const fmpz_mat_t A)
{
    slong d = fmpz_mat_nrows(A);
    fmpz_poly_t pol;
    fmpz_t c;
    slong k;
    int res = 1;

    fmpz_poly_init(pol);
    fmpz_init(c);

    fmpz_mat_charpoly(pol, A);

    /* Descartes' rule of signs: pol has only positive roots iff (-1)^k a_k are
       all positive */
    for (k = 1; (k <= d) && res; k++)
    {
        fmpz_poly_get_coeff_fmpz(c, pol, d - k);
        if (k % 2 == 1)
        {
            fmpz_neg(c, c);
        }
        if (fmpz_cmp_si(c, 0) <= 0)
        {
            res = 0;
        }
    }

    fmpz_poly_clear(pol);
    return res;
}

int
fmpz_mat_is_spd(const fmpz_mat_t A)
{
    slong d = fmpz_mat_nrows(A);
    slong k, j;

    if (fmpz_mat_ncols(A) != d)
    {
        return 0;
    }

    for (j = 0; j < d; j++)
    {
        for (k = 0; k < j; k++)
        {
            if (!fmpz_equal(fmpz_mat_entry(A, j, k), fmpz_mat_entry(A, k, j)))
            {
                return 0;
            }
        }
    }

    if (fmpz_mat_is_spd_arb(A))
    {
        return 1;
    }
    else
    {
        return fmpz_mat_is_spd_charpoly(A);
    }
}
