/*
    Copyright (C) 2014 Alex J. Best
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

/*
    Compute SNF via HNF preprocessing followed by Iliopoulos.
    The product of the HNF pivot entries is a valid modulus for Iliopoulos
    (it equals the product of the nonzero invariant factors).
    We extract only the nonzero rows from the HNF before passing to
    Iliopoulos, since it cannot handle zero invariant factors (it computes
    gcd(0, mod) = mod instead of 0).
*/
static void
_fmpz_mat_snf_via_hnf(fmpz_mat_t S, const fmpz_mat_t A)
{
    fmpz_mat_t H, H_nz, S_nz;
    fmpz_t mod;
    slong i, j, r, m = A->r, n = A->c;

    fmpz_mat_init(H, m, n);
    fmpz_mat_hnf(H, A);

    /* Count nonzero rows (= rank) and compute modulus */
    fmpz_init(mod);
    fmpz_one(mod);
    r = 0;
    for (i = 0; i < m; i++)
    {
        if (fmpz_mat_is_zero_row(H, i))
            break;
        r++;
        for (j = 0; j < n; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(H, i, j)))
            {
                fmpz_mul(mod, mod, fmpz_mat_entry(H, i, j));
                break;
            }
        }
    }

    fmpz_mat_zero(S);

    if (r > 0)
    {
        /* Extract nonzero rows and compute SNF of the r x n submatrix */
        fmpz_mat_init(H_nz, r, n);
        fmpz_mat_init(S_nz, r, n);

        for (i = 0; i < r; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(H_nz, i, j),
                        fmpz_mat_entry(H, i, j));

        fmpz_mat_snf_iliopoulos(S_nz, H_nz, mod);

        /* Embed result into the output matrix */
        for (i = 0; i < r; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(S, i, j),
                        fmpz_mat_entry(S_nz, i, j));

        fmpz_mat_clear(S_nz);
        fmpz_mat_clear(H_nz);
    }

    fmpz_clear(mod);
    fmpz_mat_clear(H);
}

void
fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A)
{
    fmpz_t det;
    slong m = A->r, n = A->c, b = fmpz_mat_max_bits(A), cutoff = 9;

    if (b <= 2)
        cutoff = 15;
    else if (b <= 4)
        cutoff = 13;
    else if (b <= 8)
        cutoff = 13;
    else if (b <= 16)
        cutoff = 11;
    else if (b <= 32)
        cutoff = 11;
    else if (b <= 64)
        cutoff = 10;

    if (FLINT_MAX(m, n) < cutoff)
        fmpz_mat_snf_kannan_bachem(S, A);
    else if (m != n)
        _fmpz_mat_snf_via_hnf(S, A);
    else
    {
        fmpz_init(det);
        fmpz_mat_det(det, A);
        if (!fmpz_is_zero(det))
        {
            fmpz_abs(det, det);
            fmpz_mat_snf_iliopoulos(S, A, det);
        }
        else
        {
            _fmpz_mat_snf_via_hnf(S, A);
        }
        fmpz_clear(det);
    }
}
