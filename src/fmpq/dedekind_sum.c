/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    fmpz_t i, j, q1, q2, r2;

    if (fmpz_is_zero(k))
    {
        fmpq_zero(s);
        return;
    }

    fmpz_init(i);
    fmpz_init(j);
    fmpz_init(q1);
    fmpz_init(q2);
    fmpz_init(r2);

    fmpz_zero(fmpq_numref(s));

    for (fmpz_one(i); fmpz_cmp(i, k) < 0; fmpz_add_ui(i, i, 1))
    {
        fmpz_fdiv_q(q1, i, k);

        fmpz_mul(j, h, i);
        fmpz_fdiv_qr(q2, r2, j, k);
        if (fmpz_is_zero(r2))
            continue;

        fmpz_mul(q1, q1, k);
        fmpz_sub(q1, i, q1);
        fmpz_mul_ui(q1, q1, 2);
        fmpz_sub(q1, q1, k);

        fmpz_mul(q2, q2, k);
        fmpz_sub(q2, j, q2);
        fmpz_mul_ui(q2, q2, 2);
        fmpz_sub(q2, q2, k);

        fmpz_addmul(fmpq_numref(s), q1, q2);
    }

    fmpz_mul(fmpq_denref(s), k, k);
    fmpz_mul_ui(fmpq_denref(s), fmpq_denref(s), 4);
    fmpq_canonicalise(s);

    fmpz_clear(i);
    fmpz_clear(j);
    fmpz_clear(q1);
    fmpz_clear(q2);
    fmpz_clear(r2);
}

#define _UI_MAT22_RMUL_ELEM(m11, m12, m21, m22, q) \
  do {                                             \
    mp_limb_t __t1 = m12 + q*m11;                  \
    mp_limb_t __t2 = m22 + q*m21;                  \
    m12 = m11;                                     \
    m22 = m21;                                     \
    m11 = __t1;                                    \
    m21 = __t2;                                    \
  } while (0)

void fmpq_dedekind_sum(fmpq_t sum, const fmpz_t h, const fmpz_t k)
{
    if (fmpz_cmp_ui(k, 2) <= 0 || fmpz_is_zero(h))
    {
        fmpq_zero(sum);
        return;
    }
    else if (fmpz_fits_si(k))
    {
        /*
            since the alternating sum of quotients
                s1 - s2 + - ... + sn - 3, n odd
                s1 - s2 + - ... - sn,     n even
            is < k in absolute value, we have a version for machine k
        */
        ulong a = fmpz_get_ui(k);
        ulong b = fmpz_fdiv_ui(h, a);
        ulong m11 = 1, m12 = 0, m21 = 0, m22 = 1;
        ulong sum_hi, sum_lo, t = 0, q, r;

        while (b != 0)
        {
            udiv_qrnnd(q, r, 0, a, b);
            a = b;
            b = r;
            t += q;
            _UI_MAT22_RMUL_ELEM(m11, m12, m21, m22, q);

            if (b == 0)
            {
                /* break with det -1 */
                t -= 3;
                smul_ppmm(sum_hi, sum_lo, t, m11);
                add_ssaaaa(sum_hi, sum_lo, sum_hi, sum_lo, 0, m21 + m12);
                goto set_sum;
            }

            udiv_qrnnd(q, r, 0, a, b);
            a = b;
            b = r;
            t -= q;
            _UI_MAT22_RMUL_ELEM(m11, m12, m21, m22, q);
        }

        /* break with det +1 */
        smul_ppmm(sum_hi, sum_lo, t, m11);
        t = m21 - m12;
        add_ssaaaa(sum_hi, sum_lo, sum_hi, sum_lo, FLINT_SIGN_EXT(t), t);

set_sum:

        fmpz_set_signed_uiui(fmpq_numref(sum), sum_hi, sum_lo);
        fmpz_set_ui(fmpq_denref(sum), m11);
    }
    else
    {
        _fmpq_cfrac_list_t s;
        _fmpz_mat22_t M;
        _fmpq_ball_t x;

        _fmpq_cfrac_list_init(s);
        s->length = -1;
        s->want_alt_sum = 1;

        _fmpz_mat22_init(M);
        _fmpz_mat22_one(M);

        _fmpq_ball_init(x);
        x->exact = 1;
        fmpz_set(x->left_num, k);
        fmpz_fdiv_r(x->left_den, h, k);

        if (!fmpz_is_zero(x->left_den))
        {
            _fmpq_ball_get_cfrac(s, M, 1, x);

            /* exactly one extra iteration needed if get_cfrac is working */
            FLINT_ASSERT(fmpz_divisible(x->left_num, x->left_den));

            do {
                fmpz_fdiv_qr(x->right_num, x->left_num, x->left_num, x->left_den);
                _fmpz_mat22_rmul_elem(M, x->right_num);
                _fmpq_cfrac_list_push_back(s, x->right_num);
                fmpz_swap(x->left_num, x->left_den);
            } while (!fmpz_is_zero(x->left_den));
        }

        FLINT_ASSERT(!fmpz_is_zero(M->_11));

        FLINT_ASSERT(s->want_alt_sum == M->det);
        if (M->det == 1)
        {
            fmpz_sub(fmpq_numref(sum), M->_21, M->_12);
        }
        else
        {
            FLINT_ASSERT(M->det == -1);
            fmpz_sub_ui(s->alt_sum, s->alt_sum, 3);
            fmpz_add(fmpq_numref(sum), M->_21, M->_12);
        }

        fmpz_swap(fmpq_denref(sum), M->_11);
        fmpz_addmul(fmpq_numref(sum), s->alt_sum, fmpq_denref(sum));

        _fmpq_ball_clear(x);
        _fmpq_cfrac_list_clear(s);
        _fmpz_mat22_clear(M);
    }

    fmpz_mul_ui(fmpq_denref(sum), fmpq_denref(sum), 12);
    fmpq_canonicalise(sum); /* extra gcd seems to be unavoidable */
}
