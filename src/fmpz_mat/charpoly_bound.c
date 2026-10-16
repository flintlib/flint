/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz.h"
#include "fmpz_mat.h"
#include "mag.h"

void fmpz_mat_charpoly_bound(fmpz_t bound, const fmpz_mat_t A)
{
    slong n = fmpz_mat_nrows(A);
    slong i, j;
    mag_t t, rprod, cprod, rbound, cbound;
    mag_struct *rnorms, *cnorms, *binoms;

    if (n == 0)
    {
        fmpz_one(bound);
        return;
    }

    rnorms = _mag_vec_init(n);
    cnorms = _mag_vec_init(n);
    binoms = _mag_vec_init(n / 2 + 2);

    mag_init(t);
    mag_init(rprod);
    mag_init(cprod);
    mag_init(rbound);
    mag_init(cbound);

    /* Squared row and column norms */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            mag_set_fmpz(t, fmpz_mat_entry(A, i, j));
            mag_fast_addmul(rnorms + i, t, t);
            mag_fast_addmul(cnorms + j, t, t);
        }
    }

    for (i = 0; i < n; i++)
    {
        mag_sqrt(rnorms + i, rnorms + i);
        mag_sqrt(cnorms + i, cnorms + i);
    }

    /* Sort norms ascending; we will read in reverse */
    qsort(rnorms, n, sizeof(mag_struct),
        (int (*)(const void *, const void *)) mag_cmp);
    qsort(cnorms, n, sizeof(mag_struct),
        (int (*)(const void *, const void *)) mag_cmp);

    mag_one(binoms + 0);
    mag_set_ui(binoms + 1, n);
    for (i = 2; i <= n / 2; i++)
    {
        mag_mul_ui(binoms + i, binoms + i - 1, n - i + 1);
        mag_div_ui(binoms + i, binoms + i, i);
    }

    mag_one(rprod);
    mag_one(cprod);
    mag_one(rbound);
    mag_one(cbound);

    for (j = 0; j < n; j++)
    {
        slong k = j + 1;
        slong binom_idx = (k <= n / 2) ? k : n - k;

        mag_mul(rprod, rprod, rnorms + (n - 1 - j));
        mag_mul(cprod, cprod, cnorms + (n - 1 - j));

        mag_mul(t, binoms + binom_idx, rprod);
        mag_max(rbound, rbound, t);
        mag_mul(t, binoms + binom_idx, cprod);
        mag_max(cbound, cbound, t);
    }

    mag_min(t, rbound, cbound);
    mag_get_fmpz(bound, t);

    _mag_vec_clear(rnorms, n);
    _mag_vec_clear(cnorms, n);
    _mag_vec_clear(binoms, n / 2 + 2);

    mag_clear(t);
    mag_clear(rprod);
    mag_clear(cprod);
    mag_clear(rbound);
    mag_clear(cbound);
}

