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
#if FLINT_USES_BLAS
#include "cblas.h"
#endif

static int
fmpz_get_sgnbit_mpn2(mp_ptr r, const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        slong v = *x;
        r[0] = FLINT_ABS(v);
        r[1] = 0;
        return v < 0;
    }
    else
    {
        __mpz_struct * p = COEFF_TO_PTR(*x);
        slong sz = p->_mp_size;
        r[0] = p->_mp_d[0];
        if (sz == 2 || sz == -2)
            r[1] = p->_mp_d[1];
        else
            r[1] = 0;
        return sz < 0;
    }
}

#define nn_mul_2x1(r2, r1, r0, a1, a0, b0)                  \
    do {                                                    \
        mp_limb_t t1;                                       \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
    } while (0)

#define nn_mul_2x2(r3, r2, r1, r0, a1, a0, b1, b0)          \
    do {                                                    \
        mp_limb_t t1, t2, t3;                               \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
        umul_ppmm(t1, t2, a0, b1);                          \
        umul_ppmm(r3, t3, a1, b1);                          \
        add_ssaaaa(r3, t1, r3, t1, 0, t3);                  \
        add_ssaaaa(r2, r1, r2, r1, t1, t2);                 \
        r3 += r2 < t1;                                      \
} while (0)

FLINT_DLL void
fmpz_mat_mul_1(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
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

FLINT_DLL void
fmpz_mat_mul_2a(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
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

FLINT_DLL void
fmpz_mat_mul_2b(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
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

FLINT_DLL void
fmpz_mat_mul_4(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, ac, br, bc;
    slong i, j, k;
    mp_ptr AL, BL;
    char *AS, *BS;
    TMP_INIT;

    ar = fmpz_mat_nrows(A);
    ac = fmpz_mat_ncols(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    TMP_START;

    AL = TMP_ALLOC(2 * sizeof(mp_limb_t) * ar * ac);
    BL = TMP_ALLOC(2 * sizeof(mp_limb_t) * br * bc);
    AS = TMP_ALLOC(sizeof(char) * ar * ac);
    BS = TMP_ALLOC(sizeof(char) * br * bc);

    for (i = 0; i < ar; i++)
        for (j = 0; j < ac; j++)
            AS[i * ac + j] = fmpz_get_sgnbit_mpn2(AL + 2 * i * ac + 2 * j,
                fmpz_mat_entry(A, i, j));

    for (i = 0; i < br; i++)
        for (j = 0; j < bc; j++)
            BS[i * bc + j] = fmpz_get_sgnbit_mpn2(BL + 2 * i * bc + 2 * j,
                fmpz_mat_entry(B, i, j));

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mp_limb_t s[4];
            mp_limb_t t[4];
            mp_limb_t u[4];

            flint_mpn_zero(s, 4);
            flint_mpn_zero(t, 4);
            flint_mpn_zero(u, 4);

            for (k = 0; k < br; k++)
            {
                mp_limb_t ah, al, bh, bl;
                mp_srcptr aptr, bptr;

                aptr = AL + 2 * i * ac + 2 * k;
                bptr = BL + 2 * k * bc + 2 * j;

                al = aptr[0];
                ah = aptr[1];
                bl = bptr[0],
                bh = bptr[1];

                if (ah == 0 && bh == 0)
                {
                    if (al == 0 || bl == 0)
                        continue;

                    umul_ppmm(t[1], t[0], al, bl);

                    if (AS[i * ac + k] == BS[k * bc + j])
                        add_sssaaaaaa(u[2], u[1], u[0], u[2], u[1], u[0], 0, t[1], t[0]);
                    else
                        sub_dddmmmsss(u[2], u[1], u[0], u[2], u[1], u[0], 0, t[1], t[0]);
                }
                else
                {
                    if (ah == 0)
                    {
                        nn_mul_2x1(t[2], t[1], t[0], bh, bl, al);
                        t[3] = 0;
                    }
                    else if (bh == 0)
                    {
                        nn_mul_2x1(t[2], t[1], t[0], ah, al, bl);
                        t[3] = 0;
                    }
                    else
                    {
                        nn_mul_2x2(t[3], t[2], t[1], t[0], ah, al, bh, bl);
                    }

                    if (AS[i * ac + k] == BS[k * bc + j])
                        mpn_add_n(s, s, t, 4);
                    else
                        mpn_sub_n(s, s, t, 4);
                }
            }

            if (((mp_limb_signed_t) u[2]) >= 0)
            {
                s[3] += mpn_add_n(s, s, u, 3);
            }
            else
            {
                sub_dddmmmsss(u[2], u[1], u[0], 0, 0, 0, u[2], u[1], u[0]);
                s[3] -= mpn_sub_n(s, s, u, 3);
            }

            if (((mp_limb_signed_t) s[3]) >= 0)
            {
                fmpz_set_ui_array(fmpz_mat_entry(C, i, j), s, 4);
            }
            else
            {
                mpn_neg_n(s, s, 4);
                fmpz_set_ui_array(fmpz_mat_entry(C, i, j), s, 4);
                fmpz_neg(fmpz_mat_entry(C, i, j), fmpz_mat_entry(C, i, j));
            }
        }
    }

    TMP_END;
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

    /*make more use of strassen on 1 thread and mul_blas on more than one thread */
    if (
#if FLINT_USES_BLAS
        (openblas_get_num_threads() == 1 && abits + bbits >= 10) &&
#endif
        ((abits + bbits >= 20  && abits + bbits < 60 && dim > 150) || (dim > 500)))
            fmpz_mat_mul_strassen(C, A, B);
#if FLINT_USES_BLAS
    else if ((openblas_get_num_threads() > 1 && abits + bbits >= 15) &&
	((abits + bbits >= 400 && dim > 40) ||
	 (abits + bbits >= 120 && dim > 60) || dim > 500))
	    fmpz_mat_mul_blas(C, A, B);
#endif
    else if (bits <= FLINT_BITS - 2)
    {
        if ((abits + bbits < 18 && dim > 100) ||
#if FLINT_USES_BLAS
	    (abits + bbits < 50 && openblas_get_num_threads() > 1 && dim > 100) ||
#endif
	    ((dim > 160 && abits + bbits <= 20) || dim > 600)) /* tuning param */
                _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (dim > 160) /* tuning param */
            fmpz_mat_mul_strassen(C, A, B);
        else
            fmpz_mat_mul_1(C, A, B);
    }
    else if (abits <= FLINT_BITS - 2 && bbits <= FLINT_BITS - 2)
    {
        if (dim > 400) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (bits <= 2 * FLINT_BITS - 1)
            fmpz_mat_mul_2a(C, A, B);
        else
            fmpz_mat_mul_2b(C, A, B);
    }
    else if (abits <= 2 * FLINT_BITS && bbits <= 2 * FLINT_BITS
                                     && bits <= 4 * FLINT_BITS - 1)
    {
        if (dim > 40) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else
            fmpz_mat_mul_4(C, A, B);
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

