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
_fmpz_mat_mul_4(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
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

