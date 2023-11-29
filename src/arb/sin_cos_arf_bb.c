/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "fmpz_vec.h"
#include "arb.h"
#include "acb.h"

slong _arb_compute_bs_exponents(slong * tab, slong n);
slong _arb_get_exp_pos(const slong * tab, slong step);

static void
bsplit(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const slong * xexp,
    const fmpz * xpow, flint_bitcnt_t r, slong a, slong b)
{
    int cc;

    if (b - a == 1)
    {
        cc = flint_ctz((2 * a + 2));
        fmpz_neg_ui(Q, (2 * a + 2) >> cc);
        fmpz_mul_ui(Q, Q, 2 * a + 3);
        *Qexp = 2 * r + cc;

        fmpz_set(T, xpow);
    }
    else if (b - a == 2)
    {
        fmpz_mul2_uiui(T, xpow, (2 * a + 4), (2 * a + 5));
        fmpz_mul_2exp(T, T, 2 * r);
        fmpz_neg(T, T);
        fmpz_add(T, T, xpow + 1);

        cc = flint_ctz((2 * a + 4));
        fmpz_neg_ui(Q, (2 * a + 4) >> cc);
        fmpz_mul_ui(Q, Q, 2 * a + 5);
        *Qexp = 2 * r + cc;

        cc = flint_ctz((2 * a + 2));
        fmpz_mul2_uiui(Q, Q, (2 * a + 2) >> cc, (2 * a + 3));
        fmpz_neg(Q, Q);
        *Qexp += 2 * r + cc;
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
cos_bsplit_struct;

typedef cos_bsplit_struct cos_bsplit_t[1];

static void cos_bsplit_init(cos_bsplit_t x, void * args)
{
    fmpz_init(x->T);
    fmpz_init(x->Q);
}

static void cos_bsplit_clear(cos_bsplit_t x, void * args)
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
cos_bsplit_args;

static void
cos_bsplit_merge(cos_bsplit_t res, cos_bsplit_t L, cos_bsplit_t R, cos_bsplit_args * args)
{
    slong i, step;

    slong a = L->a;
    slong b = R->b;

    step = (b - a) / 2;

    fmpz_mul(res->T, L->T, R->Q);
    fmpz_mul_2exp(res->T, res->T, R->Qexp);

    /* find x^step in table */
    i = _arb_get_exp_pos(args->xexp, step);
    fmpz_addmul(res->T, args->xpow + i, R->T);
    fmpz_zero(R->T);

    fmpz_mul(res->Q, L->Q, R->Q);
    res->Qexp = L->Qexp + R->Qexp;

    res->a = L->a;  /* actually a no-op because of aliasing */
    res->b = R->b;
}

static void
cos_bsplit_basecase(cos_bsplit_t res, slong a, slong b, cos_bsplit_args * args)
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
    cos_bsplit_t s;
    cos_bsplit_args args;
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
        (bsplit_basecase_func_t) cos_bsplit_basecase,
        (bsplit_merge_func_t) cos_bsplit_merge,
        sizeof(cos_bsplit_struct),
        (bsplit_init_func_t) cos_bsplit_init,
        (bsplit_clear_func_t) cos_bsplit_clear,
        &args, a, b, 4, max_threads, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);

    *T = *s->T;
    *Q = *s->Q;
    *Qexp = s->Qexp;
}

/* todo: also allow computing cos, using the same table... */
void
_arb_sin_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
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

/*
Determine N such that the error is bounded by 2^-prec when summing the
Taylor series of sin(x) up to term x^(2N+1) inclusive. We choose an N with
many trailing zeros to improve efficiency of the binary splitting.
*/
static slong
bs_num_terms(slong mag, slong prec)
{
    slong N;

    N = _arb_exp_taylor_bound(mag, prec);
    N = N / 2 - 1;
    N = FLINT_MAX(N, 1);

    if (N > 10000)
        while (N % 128 != 0)
            N++;

    if (N > 1000)
        while (N % 16 != 0)
            N++;

    if (N > 100)
        while (N % 2 != 0)
            N++;

    return N;
}

void
arb_sin_cos_fmpz_div_2exp_bsplit(arb_t wsin, arb_t wcos, const fmpz_t x, flint_bitcnt_t r, slong prec)
{
    fmpz_t T, Q;
    slong N, xmag;
    flint_bitcnt_t Qexp[1];

    /* slightly reduce memory usage at very high precision */
    arb_zero(wsin);
    arb_zero(wcos);

    fmpz_init(T);
    fmpz_init(Q);

    if (r > prec)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    /* Binary splitting (+1 fixed-point ulp truncation error). */
    xmag = fmpz_bits(x) - r;
    N = bs_num_terms(xmag, prec);
   _arb_sin_sum_bs_powtab(T, Q, Qexp, x, r, N);

    /* we still need to multiply and add x/2^r to get sine */
    fmpz_mul(T, T, x);
    Qexp[0] += r;

    /* T = T / Q  (+1 fixed-point ulp error). */
    if (Qexp[0] >= prec)
        fmpz_tdiv_q_2exp(T, T, Qexp[0] - prec);
    else
        fmpz_mul_2exp(T, T, prec - Qexp[0]);

    arb_fmpz_divapprox(T, T, Q);

    fmpz_mul_2exp(Q, x, prec - r);
    fmpz_add(T, T, Q);

    /* T = sin(u) with at most 2 fixed-point ulp error. */
    arf_set_fmpz(arb_midref(wsin), T);
    arf_mul_2exp_si(arb_midref(wsin), arb_midref(wsin), -prec);
    mag_set_ui_2exp_si(arb_radref(wsin), 2, -prec);

    /* compute cos from sin */
    arb_mul(wcos, wsin, wsin, prec);
    arb_sub_ui(wcos, wcos, 1, prec);
    arb_neg(wcos, wcos);
    arb_sqrt(wcos, wcos, prec);

    fmpz_clear(T);
    fmpz_clear(Q);
}

typedef struct
{
    acb_ptr vs;
    fmpz * u;
    slong * r;
    slong wp;
}
work_t;

static void
worker(slong iter, work_t * work)
{
    arb_sin_cos_fmpz_div_2exp_bsplit(acb_imagref(work->vs + iter), acb_realref(work->vs + iter), work->u + iter, work->r[iter], work->wp);
}

/* parallel product of complex numbers; destructive (overwrites input) */

typedef struct
{
    acb_ptr vec;
    slong prec;
}
pwork_t;

static void
pbasecase(acb_t res, slong a, slong b, pwork_t * work)
{
    if (b - a == 0)
    {
        acb_one(res);
    }
    else if (b - a == 1)
    {
        acb_swap(res, work->vec + a);
    }
    else
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }
}

static void
pmerge(acb_t res, acb_t a, acb_t b, pwork_t * work)
{
    arb_t tmp1;
    arb_ptr zsin, zcos, wsin, wcos;
    slong wp = work->prec;

    zcos = acb_realref(res);
    zsin = acb_imagref(res);
    wcos = acb_realref(b);
    wsin = acb_imagref(b);

    arb_init(tmp1);

    arb_add(tmp1, zsin, zcos, wp);
    arb_mul(zcos, zcos, wcos, wp);
    arb_add(wcos, wcos, wsin, wp);
    arb_mul(wsin, wsin, zsin, wp);
    arb_mul(wcos, tmp1, wcos, wp);
    arb_zero(tmp1);
    arb_sub(zsin, wcos, wsin, wp);
    arb_zero(wcos);
    arb_sub(zsin, zsin, zcos, wp);
    arb_sub(zcos, zcos, wsin, wp);
    arb_zero(wsin);

    arb_clear(tmp1);
}

void
_acb_vec_prod_bsplit_threaded(acb_t res, acb_ptr vec, slong len, slong prec)
{
    pwork_t work;

    work.vec = vec;
    work.prec = prec;

    flint_parallel_binary_splitting(res,
        (bsplit_basecase_func_t) pbasecase,
        (bsplit_merge_func_t) pmerge,
        sizeof(acb_struct),
        (bsplit_init_func_t) acb_init,
        (bsplit_clear_func_t) acb_clear,
        &work, 0, len, 1, -1, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);
}

void
arb_sin_cos_arf_bb(arb_t zsin, arb_t zcos, const arf_t x, slong prec)
{
    slong k, bits, r, xmag, q, wp;
    slong argred_bits, start_bits;
    int inexact, negative;
    fmpz_t t, u;
    arb_t wcos, wsin, tmp1;

    if (zsin == NULL)
    {
        arb_init(tmp1);
        arb_sin_cos_arf_bb(tmp1, zcos, x, prec);
        arb_clear(tmp1);
        return;
    }

    if (zcos == NULL)
    {
        arb_init(tmp1);
        arb_sin_cos_arf_bb(zsin, tmp1, x, prec);
        arb_clear(tmp1);
        return;
    }

    if (arf_is_zero(x))
    {
        arb_zero(zsin);
        arb_one(zcos);
        return;
    }

    xmag = arf_abs_bound_lt_2exp_si(x);
    negative = arf_sgn(x) < 0;

    /* We assume that this function only gets called with something
       reasonable as input (huge/tiny input will be handled by
       the main sin/cos wrapper). */
    if (arf_is_special(x) || arf_cmpabs_d(x, 3.15) > 0 || xmag < -2 * prec - 100)
    {
        flint_throw(FLINT_ERROR, "arb_sin_cos_arf_bb: unexpectedly large/small input\n");
    }

    argred_bits = 24;
    start_bits = argred_bits * 3;

    q = FLINT_MAX(0, xmag + argred_bits);
    if (q <= 2)
        q = 0;

    wp = prec + 10 + 2 * (q - xmag) + 2 * FLINT_BIT_COUNT(prec);

    fmpz_init(t);
    fmpz_init(u);
    arb_init(wcos);
    arb_init(wsin);
    arb_init(tmp1);

    /* Convert x/2^q to a fixed-point number. */
    inexact = arf_get_fmpz_fixed_si(t, x, -wp + q);
    fmpz_abs(t, t);

    /* Aliasing of z and x is safe now that only use t. */
    /* Start with z = 1. */
    arb_one(zcos);
    arb_zero(zsin);

    /* We have two ways to parallelize the BB algorithm: run
       the main loop serially and rely on parallel binary splitting,
       or compute all the sines/cosines in parallel. The latter is
       more efficient (empirically about 1.7x) but uses more memory,
       so we fall back on a serial main loop at high enough precision. */
    if (flint_get_num_available_threads() == 1 || prec >= 4e8)
    {
        /* Bit-burst loop. */
        for (bits = start_bits; !fmpz_is_zero(t); bits *= 3)
        {
            /* Extract bits. */
            r = FLINT_MIN(bits, wp);
            fmpz_tdiv_q_2exp(u, t, wp - r);

            arb_sin_cos_fmpz_div_2exp_bsplit(wsin, wcos, u, r, wp);

            /* Remove used bits. */
            fmpz_mul_2exp(u, u, wp - r);
            fmpz_sub(t, t, u);

            /* zsin, zcos = zsin wcos + zcos wsin, zcos wcos - zsin wsin */
            /* using karatsuba */
            arb_add(tmp1, zsin, zcos, wp);
            arb_mul(zcos, zcos, wcos, wp);
            arb_add(wcos, wcos, wsin, wp);
            arb_mul(wsin, wsin, zsin, wp);
            arb_mul(tmp1, tmp1, wcos, wp);
            arb_sub(zsin, tmp1, wsin, wp);
            arb_sub(zsin, zsin, zcos, wp);
            arb_sub(zcos, zcos, wsin, wp);
            arb_zero(tmp1);  /* slightly reduce memory usage */
        }
    }
    else
    {
        acb_ptr vs;
        fmpz * us;
        slong * rs;
        slong num = 0;

        vs = _acb_vec_init(FLINT_BITS);
        us = _fmpz_vec_init(FLINT_BITS);
        rs = flint_malloc(sizeof(slong) * FLINT_BITS);

        /* Bit-burst loop. */
        for (bits = start_bits; !fmpz_is_zero(t); bits *= 3)
        {
            /* Extract bits. */
            r = FLINT_MIN(bits, wp);
            fmpz_tdiv_q_2exp(u, t, wp - r);

            if (!fmpz_is_zero(u))
            {
                fmpz_set(us + num, u);
                rs[num] = r;
                num++;
            }

            /* Remove used bits. */
            fmpz_mul_2exp(u, u, wp - r);
            fmpz_sub(t, t, u);
        }

        /* todo: only allocate as many temporaries as threads,
           reducing memory */
        {
            work_t work;

            work.vs = vs;
            work.u = us;
            work.r = rs;
            work.wp = wp;

            flint_parallel_do((do_func_t) worker, &work, num, -1, FLINT_PARALLEL_STRIDED);
        }

        {
            acb_t z;
            *acb_realref(z) = *zcos;
            *acb_imagref(z) = *zsin;

            _acb_vec_prod_bsplit_threaded(z, vs, num, wp);

            *zcos = *acb_realref(z);
            *zsin = *acb_imagref(z);
        }

        _acb_vec_clear(vs, FLINT_BITS);
        _fmpz_vec_clear(us, FLINT_BITS);
        flint_free(rs);
    }

    /* Initial fixed-point truncation error. */
    if (inexact)
    {
        arb_add_error_2exp_si(zcos, -wp);
        arb_add_error_2exp_si(zsin, -wp);
    }

    if (q != 0)
    {
        /* cos(x) = 2 cos(x/2)^2 - 1 */
        for (k = 0; k < q; k++)
        {
            arb_mul(zcos, zcos, zcos, wp);
            arb_mul_2exp_si(zcos, zcos, 1);
            arb_sub_ui(zcos, zcos, 1, wp);
        }

        arb_mul(tmp1, zcos, zcos, wp);
        arb_sub_ui(tmp1, tmp1, 1, wp);
        arb_neg(tmp1, tmp1);
        arb_sqrt(zsin, tmp1, wp);
    }

    if (negative)
        arb_neg(zsin, zsin);

    arb_set_round(zsin, zsin, prec);
    arb_set_round(zcos, zcos, prec);

    fmpz_clear(t);
    fmpz_clear(u);
    arb_clear(wcos);
    arb_clear(wsin);
    arb_clear(tmp1);
}

