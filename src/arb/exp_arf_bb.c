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

/*
Determine N such that the error is bounded by 2^-prec when summing the
Taylor series of exp(x) up to term x^N inclusive. We choose an N with
many trailing zeros to improve efficiency of the binary splitting.
*/
static slong
bs_num_terms(slong mag, slong prec)
{
    slong N;

    N = _arb_exp_taylor_bound(mag, prec);
    /* Convert from N exclusive to N inclusive. */
    N--;

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

typedef struct
{
    arb_ptr w;
    fmpz * u;
    slong * r;
    slong wp;
}
work_t;

static void
worker(slong iter, work_t * work)
{
    slong mag, wp, N;
    flint_bitcnt_t Qexp[1];
    fmpz_t T, Q;

    fmpz_init(T);
    fmpz_init(Q);

    wp = work->wp;

    /* Binary splitting (+1 fixed-point ulp truncation error). */
    mag = fmpz_bits(work->u + iter) - work->r[iter];
    N = bs_num_terms(mag, wp);

   _arb_exp_sum_bs_powtab(T, Q, Qexp, work->u + iter, work->r[iter], N);

    /* T = T / Q  (+1 fixed-point ulp error). */
    if (*Qexp >= wp)
        fmpz_tdiv_q_2exp(T, T, *Qexp - wp);
    else
        fmpz_mul_2exp(T, T, wp - *Qexp);

    arb_fmpz_divapprox(T, T, Q);

    /* T = 1 + T */
    fmpz_one(Q);
    fmpz_mul_2exp(Q, Q, wp);
    fmpz_add(T, T, Q);

    /* Now T = exp(u) with at most 2 fixed-point ulp error. */
    /* Set z = z * T. */
    arf_set_fmpz(arb_midref(work->w + iter), T);
    arf_mul_2exp_si(arb_midref(work->w + iter), arb_midref(work->w + iter), -wp);
    mag_set_ui_2exp_si(arb_radref(work->w + iter), 2, -wp);

    fmpz_clear(T);
    fmpz_clear(Q);
}

typedef struct
{
    arb_srcptr vec;
    slong prec;
}
pwork_t;

static void
pbasecase(arb_t res, slong a, slong b, pwork_t * work)
{
    if (b - a == 0)
    {
        arb_one(res);
    }
    else if (b - a == 1)
    {
        arb_set(res, work->vec + a);
    }
    else if (b - a == 2)
    {
        arb_mul(res, work->vec + a, work->vec + a + 1, work->prec);
    }
    else if (b - a == 3)
    {
        arb_mul(res, work->vec + a, work->vec + a + 1, work->prec);
        arb_mul(res, res, work->vec + a + 2, work->prec);
    }
    else
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }
}

static void
pmerge(arb_t res, const arb_t a, const arb_t b, pwork_t * work)
{
    arb_mul(res, a, b, work->prec);
}

void
_arb_vec_prod_bsplit_threaded(arb_t res, arb_srcptr vec, slong len, slong prec)
{
    pwork_t work;

    work.vec = vec;
    work.prec = prec;

    flint_parallel_binary_splitting(res,
        (bsplit_basecase_func_t) pbasecase,
        (bsplit_merge_func_t) pmerge,
        sizeof(arb_struct),
        (bsplit_init_func_t) arb_init,
        (bsplit_clear_func_t) arb_clear,
        &work, 0, len, 3, -1, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);
}

void
arb_exp_arf_bb(arb_t z, const arf_t x, slong prec, int minus_one)
{
    slong k, bits, r, mag, q, wp, N;
    slong argred_bits, start_bits;
    slong num_threads;
    flint_bitcnt_t Qexp[1];
    int inexact;
    fmpz_t t, u, T, Q;
    arb_t w;

    if (arf_is_zero(x))
    {
        if (minus_one)
            arb_zero(z);
        else
            arb_one(z);
        return;
    }

    if (arf_is_special(x))
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }

    mag = arf_abs_bound_lt_2exp_si(x);

    /* We assume that this function only gets called with something
       reasonable as input (huge/tiny input will be handled by
       the main exp wrapper). */
    if (mag > 200 || mag < -2 * prec - 100)
    {
        flint_throw(FLINT_ERROR, "arb_exp_arf_bb: unexpectedly large/small input\n");
    }

    if (prec < 100000000)
    {
        argred_bits = 16;
        start_bits = 32;
    }
    else
    {
        argred_bits = 32;
        start_bits = 64;
    }

    /* Argument reduction: exp(x) -> exp(x/2^q). This improves efficiency
       of the first iteration in the bit-burst algorithm. */
    q = FLINT_MAX(0, mag + argred_bits);

    /* Determine working precision. */
    wp = prec + 10 + 2 * q + 2 * FLINT_BIT_COUNT(prec);
    if (minus_one && mag < 0)
        wp += (-mag);

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(Q);
    fmpz_init(T);
    arb_init(w);

    /* Convert x/2^q to a fixed-point number. */
    inexact = arf_get_fmpz_fixed_si(t, x, -wp + q);

    /* Aliasing of z and x is safe now that only use t. */
    /* Start with z = 1. */
    arb_one(z);

    num_threads = flint_get_num_available_threads();

    /* We have two ways to parallelize the BB algorithm: run
       the main loop serially and rely on parallel binary splitting,
       or compute all the exponentials in parallel. The latter is
       more efficient (empirically about 1.3x) but uses more memory,
       so we fall back on a serial main loop at high enough precision. */
    if (num_threads == 1 || prec >= 1e9)
    {
        /* Bit-burst loop. */
        for (bits = start_bits; !fmpz_is_zero(t); bits *= 2)
        {
            /* Extract bits. */
            r = FLINT_MIN(bits, wp);
            fmpz_tdiv_q_2exp(u, t, wp - r);

            /* Binary splitting (+1 fixed-point ulp truncation error). */
            mag = fmpz_bits(u) - r;
            N = bs_num_terms(mag, wp);

           _arb_exp_sum_bs_powtab(T, Q, Qexp, u, r, N);

            /* T = T / Q  (+1 fixed-point ulp error). */
            if (*Qexp >= wp)
            {
                fmpz_tdiv_q_2exp(T, T, *Qexp - wp);
                arb_fmpz_divapprox(T, T, Q);
            }
            else
            {
                fmpz_mul_2exp(T, T, wp - *Qexp);
                arb_fmpz_divapprox(T, T, Q);
            }

            /* T = 1 + T */
            fmpz_one(Q);
            fmpz_mul_2exp(Q, Q, wp);
            fmpz_add(T, T, Q);

            /* Now T = exp(u) with at most 2 fixed-point ulp error. */
            /* Set z = z * T. */
            arf_set_fmpz(arb_midref(w), T);
            arf_mul_2exp_si(arb_midref(w), arb_midref(w), -wp);
            mag_set_ui_2exp_si(arb_radref(w), 2, -wp);
            arb_mul(z, z, w, wp);

            /* Remove used bits. */
            fmpz_mul_2exp(u, u, wp - r);
            fmpz_sub(t, t, u);
        }
    }
    else
    {
        arb_ptr ws;
        fmpz * us;
        slong * rs;
        slong num = 0;

        ws = _arb_vec_init(FLINT_BITS);
        us = _fmpz_vec_init(FLINT_BITS);
        rs = flint_malloc(sizeof(slong) * FLINT_BITS);

        /* Bit-burst loop. */
        for (bits = start_bits; !fmpz_is_zero(t); bits *= 2)
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

        num_threads = FLINT_MIN(num_threads, num);
        num_threads = FLINT_MAX(num_threads, 1);

        /* todo: only allocate as many temporaries as threads,
           reducing memory */
        {
            work_t work;

            work.w = ws;
            work.u = us;
            work.r = rs;
            work.wp = wp;

            flint_parallel_do((do_func_t) worker, &work, num, -1, FLINT_PARALLEL_STRIDED);
        }

        /* Parallel product. */
        _arb_vec_prod_bsplit_threaded(z, ws, num, wp);

        _arb_vec_clear(ws, FLINT_BITS);
        _fmpz_vec_clear(us, FLINT_BITS);
        flint_free(rs);
    }

    /* We have exp(x + eps) - exp(x) < 2*eps (by assumption that the argument
       reduction is large enough). */
    if (inexact)
        arb_add_error_2exp_si(z, -wp + 1);

    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(Q);
    fmpz_clear(T);
    arb_clear(w);

    /* exp(x) = exp(x/2^q)^(2^q) */
    for (k = 0; k < q; k++)
        arb_mul(z, z, z, wp);

    if (minus_one)
        arb_sub_ui(z, z, 1, wp);

    arb_set_round(z, z, prec);
}

