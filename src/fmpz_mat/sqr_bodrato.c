/*
    Copyright (C) 2012, 2024 Fredrik Johansson
    Copyright (C) 2015 Anubhav Srivastava
    Copyright (C) 2024 Marco Bodrato

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

#define E fmpz_mat_entry

static void
local_fmma(fmpz_t f, mp_limb_t a, fmpz b,
           fmpz c, fmpz d)
{
    mp_limb_t sh, sl, th, tl;

    smul_ppmm(sh, sl, a, b);
    smul_ppmm(th, tl, c, d);
    add_ssaaaa(sh, sl, sh, sl, th, tl);

    fmpz_set_signed_uiui(f, sh, sl);
}

void
fmpz_mat_sqr_bodrato(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong n = A->r;

    if (n == 0)
    {
        return;
    }
    else if (n == 1)
    {
        fmpz_mul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));
    }
    else if (n == 2)
    {
        const fmpz a = *E(A, 0, 0);
        const fmpz b = *E(A, 0, 1);
        const fmpz c = *E(A, 1, 0);
        const fmpz d = *E(A, 1, 1);

        if (!COEFF_IS_MPZ(a) && !COEFF_IS_MPZ(b) &&
            !COEFF_IS_MPZ(c) && !COEFF_IS_MPZ(d))
        {
            mp_limb_t s, t, u, v;

            smul_ppmm(s, t, a, a);
            smul_ppmm(u, v, b, c);
            add_ssaaaa(s, t, s, t, u, v);
            fmpz_set_signed_uiui(E(B, 0, 0), s, t);
            smul_ppmm(s, t, d, d);
            add_ssaaaa(s, t, s, t, u, v);
            fmpz_set_signed_uiui(E(B, 1, 1), s, t);
            smul_ppmm(s, t, a + d, b);
            fmpz_set_signed_uiui(E(B, 0, 1), s, t);
            smul_ppmm(s, t, a + d, c);
            fmpz_set_signed_uiui(E(B, 1, 0), s, t);
        }
        else
        {
            fmpz_t t;

            fmpz_init(t);

            fmpz_mul(t, E(A, 0, 1), E(A, 1, 0));

            fmpz_mul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));
            fmpz_add(E(B, 0, 0), E(B, 0, 0), t);

            fmpz_mul(E(B, 1, 1), E(A, 1, 1), E(A, 1, 1));
            fmpz_add(E(B, 1, 1), E(B, 1, 1), t);

            fmpz_add(t, E(A, 0, 0), E(A, 1, 1));
            fmpz_mul(E(B, 0, 1), E(A, 0, 1), t);
            fmpz_mul(E(B, 1, 0), E(A, 1, 0), t);

            fmpz_clear(t);
        }
    }
    else if (n == 3)
    {
      if (!COEFF_IS_MPZ(*E(A, 0, 0)) && !COEFF_IS_MPZ(*E(A, 1, 0)) &&
          !COEFF_IS_MPZ(*E(A, 0, 1)) && !COEFF_IS_MPZ(*E(A, 1, 1)) &&
          !COEFF_IS_MPZ(*E(A, 0, 2)) && !COEFF_IS_MPZ(*E(A, 1, 2)) &&
          !COEFF_IS_MPZ(*E(A, 2, 0)) && !COEFF_IS_MPZ(*E(A, 2, 1)) &&
          !COEFF_IS_MPZ(*E(A, 2, 2)))
          {
            mp_limb_t s, t, u, v, j, k;

            smul_ppmm(s, t, *E(A, 0, 2), *E(A, 2, 0));
            smul_ppmm(u, v, *E(A, 0, 1), *E(A, 1, 0));

            smul_ppmm(j, k, *E(A, 0, 0), *E(A, 0, 0));
            add_ssaaaa(j, k, j, k, u, v);
            add_ssaaaa(j, k, j, k, s, t);
            fmpz_set_signed_uiui(E(B, 0, 0), j, k);

            smul_ppmm(j, k, *E(A, 1, 2), *E(A, 2, 1));
            add_ssaaaa(u, v, u, v, j, k);
            add_ssaaaa(s, t, s, t, j, k);

            smul_ppmm(j, k, *E(A, 1, 1), *E(A, 1, 1));
            add_ssaaaa(j, k, j, k, u, v);
            fmpz_set_signed_uiui(E(B, 1, 1), j, k);

            smul_ppmm(j, k, *E(A, 2, 2), *E(A, 2, 2));
            add_ssaaaa(j, k, j, k, s, t);
            fmpz_set_signed_uiui(E(B, 2, 2), j, k);

            s = *E(A, 0, 0) + *E(A, 1, 1);
            local_fmma(E(B, 0, 1), s, *E(A, 0, 1), *E(A, 0, 2), *E(A, 2, 1));
            local_fmma(E(B, 1, 0), s, *E(A, 1, 0), *E(A, 2, 0), *E(A, 1, 2));

            t = *E(A, 0, 0) + *E(A, 2, 2);
            local_fmma(E(B, 0, 2), t, *E(A, 0, 2), *E(A, 0, 1), *E(A, 1, 2));
            local_fmma(E(B, 2, 0), t, *E(A, 2, 0), *E(A, 2, 1), *E(A, 1, 0));

            u = *E(A, 1, 1) + *E(A, 2, 2);
            local_fmma(E(B, 1, 2), u, *E(A, 1, 2), *E(A, 1, 0), *E(A, 0, 2));
            local_fmma(E(B, 2, 1), u, *E(A, 2, 1), *E(A, 0, 1), *E(A, 2, 0));
          }
        else
      {
        slong bits;
        int sign;

        bits = fmpz_mat_max_bits(A);
        sign = (bits < 0);
        bits = FLINT_ABS(bits);

        if (bits + sign <= 2 * FLINT_BITS)
        {
            bits = 2 * bits + FLINT_BIT_COUNT(n);

            /* todo: specialize for squaring */
            _fmpz_mat_mul_double_word_internal(B, A, A, sign, bits);
        }
        else
        {
            fmpz_t temp;

            fmpz_mul(E(B, 2, 2), E(A, 2, 2), E(A, 2, 2));
            fmpz_mul(E(B, 1, 1), E(A, 1, 1), E(A, 1, 1));
            fmpz_mul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));

            fmpz_init(temp);
            fmpz_mul(temp, E(A, 0, 1), E(A, 1, 0));
            fmpz_add(E(B, 0, 0), E(B, 0, 0), temp);
            fmpz_add(E(B, 1, 1), E(B, 1, 1), temp);

            fmpz_mul(temp, E(A, 1, 2), E(A, 2, 1));
            fmpz_add(E(B, 1, 1), E(B, 1, 1), temp);
            fmpz_add(E(B, 2, 2), E(B, 2, 2), temp);

            fmpz_mul(temp, E(A, 0, 2), E(A, 2, 0));
            fmpz_add(E(B, 2, 2), E(B, 2, 2), temp);
            fmpz_add(E(B, 0, 0), E(B, 0, 0), temp);

            fmpz_add(temp, E(A, 0, 0), E(A, 1, 1));
            fmpz_fmma(E(B, 0, 1), temp, E(A, 0, 1), E(A, 0, 2), E(A, 2, 1));
            fmpz_fmma(E(B, 1, 0), temp, E(A, 1, 0), E(A, 2, 0), E(A, 1, 2));

            fmpz_add(temp, E(A, 0, 0), E(A, 2, 2));
            fmpz_fmma(E(B, 0, 2), temp, E(A, 0, 2), E(A, 0, 1), E(A, 1, 2));
            fmpz_fmma(E(B, 2, 0), temp, E(A, 2, 0), E(A, 2, 1), E(A, 1, 0));

            fmpz_add(temp, E(A, 1, 1), E(A, 2, 2));
            fmpz_fmma(E(B, 1, 2), temp, E(A, 1, 2), E(A, 1, 0), E(A, 0, 2));
            fmpz_fmma(E(B, 2, 1), temp, E(A, 2, 1), E(A, 0, 1), E(A, 2, 0));

            fmpz_clear(temp);
        }
      }
    }
    else
    {
        slong anr;

        fmpz_mat_t A11, A12, A21, A22;
        fmpz_mat_t B11, B12, B21, B22;
        fmpz_mat_t X1, X2;

        anr = n / 2;

        fmpz_mat_window_init(A11, A, 0, 0, anr, anr);
        fmpz_mat_window_init(A12, A, 0, anr, anr, 2*anr);
        fmpz_mat_window_init(A21, A, anr, 0, 2*anr, anr);
        fmpz_mat_window_init(A22, A, anr, anr, 2*anr, 2*anr);

        fmpz_mat_window_init(B11, B, 0, 0, anr, anr);
        fmpz_mat_window_init(B12, B, 0, anr, anr, 2*anr);
        fmpz_mat_window_init(B21, B, anr, 0, 2*anr, anr);
        fmpz_mat_window_init(B22, B, anr, anr, 2*anr, 2*anr);

        fmpz_mat_init(X2, anr, anr);

        fmpz_mat_add(X2, A22, A12);
        fmpz_mat_sqr(B21, X2);

        fmpz_mat_sub(X2, A22, A21);
        fmpz_mat_sqr(B22, X2);

        fmpz_mat_add(X2, X2, A12);
        fmpz_mat_sqr(B11, X2);

        fmpz_mat_sub(X2, X2, A11);
        fmpz_mat_mul(B12, X2, A12);

        fmpz_mat_init(X1, anr, anr);
        fmpz_mat_mul(X1, A12, A21);
        fmpz_mat_add(B11, B11, X1);
        fmpz_mat_sub(B12, B11, B12);
        fmpz_mat_sub(B11, B21, B11);
        fmpz_mat_mul(B21, A21, X2);

        fmpz_mat_clear(X2);

        fmpz_mat_sub(B21, B11, B21);
        fmpz_mat_sub(B12, B12, B22);
        fmpz_mat_add(B22, B22, B11);
        fmpz_mat_sqr(B11, A11);

        fmpz_mat_add(B11, X1, B11);

        fmpz_mat_clear(X1);

        fmpz_mat_window_clear(A11);
        fmpz_mat_window_clear(A12);
        fmpz_mat_window_clear(A21);
        fmpz_mat_window_clear(A22);

        fmpz_mat_window_clear(B11);
        fmpz_mat_window_clear(B12);
        fmpz_mat_window_clear(B21);
        fmpz_mat_window_clear(B22);

        if (n > 2*anr) {
          {
            fmpz_mat_t Ac, Bc;
            fmpz_mat_window_init(Ac, A, 0, 2*anr, n, n);
            fmpz_mat_window_init(Bc, B, 0, 2*anr, n, n);
            fmpz_mat_mul(Bc, A, Ac);
            fmpz_mat_window_clear(Ac);
            fmpz_mat_window_clear(Bc);
          }
          {
            fmpz_mat_t Ar, Br;
            fmpz_mat_t As;
            fmpz_mat_window_init(Ar, A, 2*anr, 0, n, n);
            fmpz_mat_window_init(Br, B, 2*anr, 0, n, 2*anr);
            fmpz_mat_window_init(As, A, 0, 0, n, 2*anr);
            fmpz_mat_mul(Br, Ar, As);
            fmpz_mat_window_clear(As);
            fmpz_mat_window_clear(Ar);
            fmpz_mat_window_clear(Br);
          }
#if 0
          {
            fmpz_mat_t Ac, Ar, Bb, tmp;

            fmpz_mat_window_init(Ac, A, 0, 2*anr, 2*anr, n);
            fmpz_mat_window_init(Ar, A, 2*anr, 0, n, 2*anr);
            fmpz_mat_window_init(Bb, B, 0, 0, 2*anr, 2*anr);

            fmpz_mat_init(tmp, 2*anr, 2*anr);
            fmpz_mat_mul(tmp, Ac, Ar);
            fmpz_mat_add(Bb, Bb, tmp);
            fmpz_mat_clear(tmp);
            fmpz_mat_window_clear(Ac);
            fmpz_mat_window_clear(Ar);
            fmpz_mat_window_clear(Bb);
          }
#else
          {
            slong i, j;

            for (i = 0; i < 2*anr; ++i)
              {
                const fmpz tmptr = *E(A, i, n - 1);
                for (j = 0; j < 2*anr; ++j)
                  {
                    fmpz_addmul(E(B, i, j), &tmptr, E(A, n - 1, j));
                  }
              }
          }
#endif
        }
    }
}