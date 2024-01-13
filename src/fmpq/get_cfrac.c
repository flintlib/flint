/*
    Copyright (C) 2011 Fredrik Johansson
    Copywrite (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"


slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t f, slong limit)
{
    slong i;
    int cmp, num_sign, den_sign;
    _fmpq_ball_t x;
    _fmpq_cfrac_list_t s;
    _fmpz_mat22_t M;

#if FLINT_WANT_ASSERT
    int input_is_canonical;
#endif

    num_sign = fmpz_sgn(fmpq_denref(f));
    den_sign = fmpz_sgn(fmpq_denref(f));

    if (limit <= 0 || den_sign == 0)
    {
        if (num_sign < 0)
        {
            fmpz_neg(fmpq_numref(rem), fmpq_numref(f));
            fmpz_neg(fmpq_denref(rem), fmpq_denref(f));
        }
        else
        {
            fmpz_set(fmpq_numref(rem), fmpq_numref(f));
            fmpz_set(fmpq_denref(rem), fmpq_denref(f));
        }
        fmpz_swap(fmpq_numref(rem), fmpq_denref(rem));
        return 0;
    }

#if FLINT_WANT_ASSERT
    input_is_canonical = fmpq_is_canonical(f);
#endif

    _fmpz_mat22_init(M);
    _fmpz_mat22_one(M);

    _fmpq_ball_init(x);
    if (den_sign > 0)
    {
        fmpz_set(x->left_num, fmpq_numref(f));
        fmpz_set(x->left_den, fmpq_denref(f));
    }
    else
    {
        fmpz_neg(x->left_num, fmpq_numref(f));
        fmpz_neg(x->left_den, fmpq_denref(f));
    }
    x->exact = 1;

    _fmpq_cfrac_list_init(s);
    s->limit = limit;

    cmp = fmpz_cmp(x->left_num, x->left_den);
    if (cmp > 0)
    {
        /* 1 < x */
        _fmpq_ball_get_cfrac(s, M, 0, x);
    }
    else
    {
        /* x <= 1 */
        _fmpq_cfrac_list_push_back_zero(s);
        if (cmp >= 0 || fmpz_sgn(x->left_num) < 0)
        {
            /* x == 1 or x < 0 */
            fmpz_fdiv_qr(s->array + 0, x->left_num, x->left_num, x->left_den);
        }

        fmpz_swap(x->left_num, x->left_den);

        if (!fmpz_is_zero(x->left_den))
            _fmpq_ball_get_cfrac(s, M, 0, x);
    }

    FLINT_ASSERT(fmpz_is_zero(x->left_den) || _fmpq_ball_gt_one(x));
    FLINT_ASSERT(x->exact);

    while (s->length < s->limit && !fmpz_is_zero(x->left_den))
    {
        _fmpq_cfrac_list_push_back_zero(s);
        fmpz_fdiv_qr(s->array + s->length - 1, x->left_num,
                                                     x->left_num, x->left_den);
        fmpz_swap(x->left_num, x->left_den);
    }

    /* write remainder */
    FLINT_ASSERT(!fmpz_is_zero(x->left_num));
    fmpz_swap(fmpq_numref(rem), x->left_den);
    fmpz_swap(fmpq_denref(rem), x->left_num);

    FLINT_ASSERT(!input_is_canonical || fmpq_is_canonical(rem));

    /* write terms */
    FLINT_ASSERT(s->length <= limit);
    for (i = 0; i < s->length; i++)
        fmpz_swap(c + i, s->array + i);

    _fmpz_mat22_clear(M);
    _fmpq_ball_clear(x);
    _fmpq_cfrac_list_clear(s);

    return i;
}


slong fmpq_get_cfrac_naive(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)
{
    fmpz_t p, q;
    slong i;

    fmpz_init(p);
    fmpz_init(q);

    fmpz_set(p, fmpq_numref(x));
    fmpz_set(q, fmpq_denref(x));

    for (i = 0; i < n && !fmpz_is_zero(q); i++)
    {
        fmpz_fdiv_qr(c + i, p, p, q);
        fmpz_swap(p, q);
    }

    fmpz_set(fmpq_numref(rem), q);
    fmpz_set(fmpq_denref(rem), p);
    fmpq_canonicalise(rem);

    fmpz_clear(p);
    fmpz_clear(q);

    return i;
}
