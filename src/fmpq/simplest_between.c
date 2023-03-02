/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"


void _fmpq_simplest_between(fmpz_t mid_num, fmpz_t mid_den,
                                        const fmpz_t l_num, const fmpz_t l_den,
                                        const fmpz_t r_num, const fmpz_t r_den)
{
    fmpz_t q, r;
    _fmpq_cfrac_list_t s;
    _fmpz_mat22_t M;
    _fmpq_ball_t x;

    FLINT_ASSERT(fmpz_sgn(l_den) > 0);
    FLINT_ASSERT(fmpz_sgn(r_den) > 0);
    FLINT_ASSERT(_fmpq_cmp(l_num, l_den, r_num, r_den) <= 0);

    fmpz_init(q);
    fmpz_init(r);

    _fmpq_cfrac_list_init(s);
    s->length = -1; /* no write */

    _fmpz_mat22_init(M);
    _fmpz_mat22_one(M);

    _fmpq_ball_init(x);
    fmpz_set(x->left_num, l_num);
    fmpz_set(x->left_den, l_den);
    fmpz_set(x->right_num, r_num);
    fmpz_set(x->right_den, r_den);
    x->exact = 0;

    if (fmpz_cmp(x->left_num, x->left_den) > 0)
    {
        /* 1 < x */
        _fmpq_ball_get_cfrac(s, M, 1, x);
    }
    else if (fmpz_sgn(x->left_num) > 0
                 && fmpz_cmp(x->right_num, x->right_den) < 0)
    {
        /* 0 < x < 1 */
        fmpz_swap(x->left_den, x->right_num);
        fmpz_swap(x->left_num, x->right_den);
        _fmpq_ball_get_cfrac(s, M, 1, x);
        fmpz_zero(q);
        _fmpz_mat22_lmul_elem(M, q);
    }
    else
    {
        _fmpq_ball_t y;

        _fmpq_ball_init(y);
        y->exact = 0;

        fmpz_fdiv_qr(q, r, x->left_num, x->left_den);
        fmpz_set(y->right_num, x->left_den);
        fmpz_swap(y->right_den, r);
        fmpz_set(y->left_num, x->right_den);
        fmpz_set(y->left_den, x->right_num);
        fmpz_submul(y->left_den, x->right_den, q);

        if (_fmpq_ball_gt_one(y))
        {
            _fmpq_ball_swap(x, y);
            _fmpq_ball_get_cfrac(s, M, 1, x);
            _fmpz_mat22_lmul_elem(M, q);
        }
        _fmpq_ball_clear(y);
    }

    fmpz_cdiv_q(q, x->left_num, x->left_den);

    FLINT_ASSERT(_fmpq_cmp_fmpz(x->right_num, x->right_den, q) >= 0);

    FLINT_ASSERT(M->det == 1 || M->det == -1);
    fmpz_swap(mid_num, M->_12);
    fmpz_addmul(mid_num, M->_11, q);
    fmpz_swap(mid_den, M->_22);
    fmpz_addmul(mid_den, M->_21, q);

    fmpz_clear(q);
    fmpz_clear(r);
    _fmpq_cfrac_list_clear(s);
    _fmpq_ball_clear(x);
    _fmpz_mat22_clear(M);

    FLINT_ASSERT(_fmpq_is_canonical(mid_num, mid_den));

    return;
}


void fmpq_simplest_between(fmpq_t mid, const fmpq_t l, const fmpq_t r)
{
    if (fmpq_cmp(l, r) <= 0)
    {
        _fmpq_simplest_between(fmpq_numref(mid), fmpq_denref(mid),
                               fmpq_numref(l), fmpq_denref(l),
                               fmpq_numref(r), fmpq_denref(r));
    }
    else
    {
        _fmpq_simplest_between(fmpq_numref(mid), fmpq_denref(mid),
                               fmpq_numref(r), fmpq_denref(r),
                               fmpq_numref(l), fmpq_denref(l));
    }
}

