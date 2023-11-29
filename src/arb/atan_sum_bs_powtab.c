/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "fmpz_vec.h"
#include "arb.h"

slong _arb_compute_bs_exponents(slong * tab, slong n);

slong _arb_get_exp_pos(const slong * tab, slong step);

static void
bsplit(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const slong * xexp,
    const fmpz * xpow, flint_bitcnt_t r, slong a, slong b)
{
    if (b - a == 1)
    {
        fmpz_set(T, xpow);

        if (a % 2 == 0)
            fmpz_neg_ui(Q, 2 * a + 3);
        else
            fmpz_set_ui(Q, 2 * a + 3);

        *Qexp = 2 * r;
    }
    else if (b - a == 2)
    {
        fmpz_mul_ui(T, xpow, 2 * a + 5);
        fmpz_mul_2exp(T, T, 2 * r);
        fmpz_submul_ui(T, xpow + 1, 2 * a + 3);

        if (a % 2 == 1)
            fmpz_neg(T, T);

        fmpz_neg_ui(Q, 2 * a + 3);
        fmpz_mul_ui(Q, Q, 2 * a + 5);
        *Qexp = 4 * r;
    }
    else
    {
        slong step, m, i;
        flint_bitcnt_t Q2exp[1];
        fmpz_t Q2, T2;

        step = (b - a) / 2;
        m = a + step;

        fmpz_init(Q2);
        fmpz_init(T2);

        bsplit(T,  Q,  Qexp,  xexp, xpow, r, a, m);
        bsplit(T2, Q2, Q2exp, xexp, xpow, r, m, b);

        fmpz_mul(T, T, Q2);
        fmpz_mul_2exp(T, T, *Q2exp);

        /* find x^step in table */
        i = _arb_get_exp_pos(xexp, step);
        fmpz_mul(T2, T2, Q);
        fmpz_addmul(T, xpow + i, T2);
        fmpz_clear(T2);

        fmpz_mul(Q, Q, Q2);
        *Qexp = *Qexp + *Q2exp;
        fmpz_clear(Q2);
    }
}

typedef struct
{
    fmpz_t T;
    fmpz_t Q;
    flint_bitcnt_t Qexp;
    slong a;
    slong b;
}
atan_bsplit_struct;

typedef atan_bsplit_struct atan_bsplit_t[1];

static void atan_bsplit_init(atan_bsplit_t x, void * args)
{
    fmpz_init(x->T);
    fmpz_init(x->Q);
}

static void atan_bsplit_clear(atan_bsplit_t x, void * args)
{
    fmpz_clear(x->T);
    fmpz_clear(x->Q);
}

typedef struct
{
    const slong * xexp;
    const fmpz * xpow;
    flint_bitcnt_t r;
}
atan_bsplit_args;

static void
atan_bsplit_merge(atan_bsplit_t res, atan_bsplit_t L, atan_bsplit_t R, atan_bsplit_args * args)
{
    slong i, step;

    slong a = L->a;
    slong b = R->b;

    step = (b - a) / 2;

    fmpz_mul(res->T, L->T, R->Q);
    fmpz_mul_2exp(res->T, res->T, R->Qexp);

    /* find x^step in table */
    i = _arb_get_exp_pos(args->xexp, step);

    fmpz_mul(R->T, R->T, L->Q);
    fmpz_addmul(res->T, args->xpow + i, R->T);
    fmpz_zero(R->T);

    fmpz_mul(res->Q, L->Q, R->Q);
    res->Qexp = L->Qexp + R->Qexp;

    res->a = L->a;  /* actually a no-op because of aliasing */
    res->b = R->b;
}

static void
atan_bsplit_basecase(atan_bsplit_t res, slong a, slong b, atan_bsplit_args * args)
{
    bsplit(res->T, res->Q, &(res->Qexp), args->xexp, args->xpow, args->r, a, b);
    res->a = a;
    res->b = b;
}

static void
bsplit2(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const slong * xexp,
    const fmpz * xpow, flint_bitcnt_t r, slong a, slong b)
{
    atan_bsplit_t s;
    atan_bsplit_args args;
    slong max_threads;
    slong prec_hint;

    args.xexp = xexp;
    args.xpow = xpow;
    args.r = r;

    *s->T = *T;
    *s->Q = *Q;

    max_threads = flint_get_num_threads();

    prec_hint = 2 * (b - a) * FLINT_MAX(r, 1);

    if (prec_hint < 30000)
        max_threads = 1;
    else if (prec_hint < 1000000)
        max_threads = FLINT_MIN(2, max_threads);
    else if (prec_hint < 5000000)
        max_threads = FLINT_MIN(4, max_threads);
    else
        max_threads = FLINT_MIN(8, max_threads);

    flint_parallel_binary_splitting(s,
        (bsplit_basecase_func_t) atan_bsplit_basecase,
        (bsplit_merge_func_t) atan_bsplit_merge,
        sizeof(atan_bsplit_struct),
        (bsplit_init_func_t) atan_bsplit_init,
        (bsplit_clear_func_t) atan_bsplit_clear,
        &args, a, b, 4, max_threads, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);

    *T = *s->T;
    *Q = *s->Q;
    *Qexp = s->Qexp;
}

void
_arb_atan_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N)
{
    slong * xexp;
    slong length, i;
    fmpz * xpow;

    /* compute the powers of x^2 that will appear (at least x^2) */
    xexp = flint_calloc(2 * FLINT_BITS, sizeof(slong));
    length = _arb_compute_bs_exponents(xexp, N);

    xpow = _fmpz_vec_init(length);
    fmpz_mul(xpow, x, x);

    /* build x^i table */
    for (i = 1; i < length; i++)
    {
        if (xexp[i] == 2 * xexp[i-1])
        {
            fmpz_mul(xpow + i, xpow + i - 1, xpow + i - 1);
        }
        else if (xexp[i] == 2 * xexp[i-2]) /* prefer squaring if possible */
        {
            fmpz_mul(xpow + i, xpow + i - 2, xpow + i - 2);
        }
        else if (xexp[i] == 2 * xexp[i-1] + 1)
        {
            fmpz_mul(xpow + i, xpow + i - 1, xpow + i - 1);
            fmpz_mul(xpow + i, xpow + i, xpow);
        }
        else if (xexp[i] == 2 * xexp[i-2] + 1)
        {
            fmpz_mul(xpow + i, xpow + i - 2, xpow + i - 2);
            fmpz_mul(xpow + i, xpow + i, xpow);
        }
        else
        {
            flint_throw(FLINT_ERROR, "power table has the wrong structure!\n");
        }
    }

    if (flint_get_num_available_threads() == 1)
        bsplit(T, Q, Qexp, xexp, xpow, r, 0, N);
    else
        bsplit2(T, Q, Qexp, xexp, xpow, r, 0, N);

    _fmpz_vec_clear(xpow, length);
    flint_free(xexp);
}

