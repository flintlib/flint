/*
    Copyright (C) 2010,2011,2018 Fredrik Johansson
    Copyright (C) 2016 Aaditya Thakkar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void _fmpz_mat_mul_1(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, br, bc;
    slong i, j, k;

    ar = fmpz_mat_nrows(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    fmpz_mat_zero(C);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            slong s = 0;

            for (k = 0; k < br; k++)
                s += *fmpz_mat_entry(A, i, k) * (*fmpz_mat_entry(B, k, j));

            *fmpz_mat_entry(C, i, j) = s;
        }
    }
}

void _fmpz_mat_mul_2a(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, br, bc;
    slong i, j, k;

    ar = fmpz_mat_nrows(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mp_limb_t hi, lo, shi, slo;
            slong x, y;

            shi = slo = 0;

            for (k = 0; k < br; k++)
            {
                x = *fmpz_mat_entry(A, i, k);
                y = *fmpz_mat_entry(B, k, j);

                smul_ppmm(hi, lo, x, y);
                add_ssaaaa(shi, slo, shi, slo, hi, lo);
            }

            fmpz_set_signed_uiui(fmpz_mat_entry(C, i, j), shi, slo);
        }
    }
}

void _fmpz_mat_mul_2b(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, br, bc;
    slong i, j, k;

    ar = fmpz_mat_nrows(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mp_limb_t hi, lo, cy, shh, shi, slo;
            slong x, y;

            shh = shi = slo = 0;

            for (k = 0; k < br; k++)
            {
                x = *fmpz_mat_entry(A, i, k);
                y = *fmpz_mat_entry(B, k, j);

                smul_ppmm(hi, lo, x, y);
                add_sssaaaaaa(cy, shi, slo, 0, shi, slo, 0, hi, lo);
                shh += (0 <= (slong) hi) ? cy : cy - 1;
            }

            fmpz_set_signed_uiuiui(fmpz_mat_entry(C, i, j), shh, shi, slo);
        }
    }
}

void
fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, br, bc;
    slong abits, bbits, bits;
    slong i, j, dim;

    ar = fmpz_mat_nrows(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    if (ar == 0 || br == 0 || bc == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        fmpz_mat_t T;
        fmpz_mat_init(T, ar, bc);
        fmpz_mat_mul(T, A, B);
        fmpz_mat_swap(C, T);
        fmpz_mat_clear(T);
        return;
    }

    if (br == 1)
    {
        for (i = 0; i < ar; i++)
            for (j = 0; j < bc; j++)
                fmpz_mul(fmpz_mat_entry(C, i, j),
                    fmpz_mat_entry(A, i, 0), fmpz_mat_entry(B, 0, j));
        return;
    }

    if (br == 2)
    {
        for (i = 0; i < ar; i++)
            for (j = 0; j < bc; j++)
                fmpz_fmma(fmpz_mat_entry(C, i, j),
                    fmpz_mat_entry(A, i, 0), fmpz_mat_entry(B, 0, j),
                    fmpz_mat_entry(A, i, 1), fmpz_mat_entry(B, 1, j));
        return;
    }

    abits = fmpz_mat_max_bits(A);
    bbits = fmpz_mat_max_bits(B);
    abits = FLINT_ABS(abits);
    bbits = FLINT_ABS(bbits);

    if (abits == 0 || bbits == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    bits = abits + bbits + FLINT_BIT_COUNT(br) + 1;
    dim = FLINT_MIN(ar, bc);
    dim = FLINT_MIN(dim, br);

    if (bits <= FLINT_BITS - 2)
    {
        if ((dim > 160 && abits + bbits <= 20) || dim > 600) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (dim > 160) /* tuning param */
            fmpz_mat_mul_strassen(C, A, B);
        else
            _fmpz_mat_mul_1(C, A, B);
    }
    else if (abits <= FLINT_BITS - 2 && bbits <= FLINT_BITS - 2)
    {
        if (dim > 400) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (bits <= 2 * FLINT_BITS - 1)
            _fmpz_mat_mul_2a(C, A, B);
        else
            _fmpz_mat_mul_2b(C, A, B);
    }
    else if (abits < 2 * FLINT_BITS && bbits < 2 * FLINT_BITS)
    {
        if (dim > 40) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else
            _fmpz_mat_mul_4(C, A, B);
    }
    else
    {
        if (dim >= 3 * FLINT_BIT_COUNT(bits))  /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (abits >= 500 && bbits >= 500 && dim >= 8)  /* tuning param */
            fmpz_mat_mul_strassen(C, A, B);
        else
            fmpz_mat_mul_classical_inline(C, A, B);
    }
}

