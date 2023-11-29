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
Determine N such that the error is bounded by 2^-prec.
We choose an N with many trailing zeros to improve efficiency
of the binary splitting.

With N = 0, 1, 2, ... the highest included term is x, x^3, x^5, ...
so the error is bounded by x^3, x^5, x^7, ... = x^(2N+3)
*/
static slong
bs_num_terms(slong mag, slong prec)
{
    slong N;

    if (mag >= 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    N = 0;

    while (mag * (2 * N + 3) > -prec)
        N++;

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

/* Argument reduction: apply atan(x) = 2 atan(x/(1+sqrt(1+x^2))) a total
   of r times, and convert result to a fixed-point number res
   together with absolute error err. With an initial inversion
   if xmag > 0.
*/
void
arb_atan_bb_reduce(fmpz_t res, mag_t err, const arf_t x, slong xmag, slong r, slong prec)
{
    int inexact;

    if (r == 0)
    {
        if (xmag <= 0)
        {
            inexact = arf_get_fmpz_fixed_si(res, x, -prec);
            mag_set_ui_2exp_si(err, inexact, -prec);
        }
        else
        {
            slong wp;
            arb_t t;

            wp = FLINT_MAX(8, prec - xmag);

            arb_init(t);
            arb_set_arf(t, x);
            arb_set_round(t, t, wp);
            arb_ui_div(t, 1, t, wp);

            mag_set(err, arb_radref(t));
            inexact = arf_get_fmpz_fixed_si(res, arb_midref(t), -prec);
            mag_add_ui_2exp_si(err, err, inexact, -prec);

            arb_clear(t);
        }
    }
    else
    {
        slong k;
        arb_t p, p2, q, q2;

        arb_init(p);
        arb_init(p2);
        arb_init(q);
        arb_init(q2);

        if (xmag <= 0)
        {
            arb_set_arf(p, x);
            arb_set_round(p, p, prec);
            arb_mul(p2, p, p, prec);
            arb_add_ui(q, p2, 1, prec);
            arb_sqrt(q, q, prec);
            arb_add_ui(q, q, 1, prec);

            for (k = 1; k < r; k++)
            {
                if (k == 1)
                {
                    arb_mul_2exp_si(q2, q, 1);
                    arb_add(q2, q2, p2, prec);
                }
                else
                {
                    arb_mul(q2, q, q, prec);
                }

                arb_add(q2, p2, q2, prec);
                arb_sqrt(q2, q2, prec);
                arb_add(q, q, q2, prec);
            }
        }
        else
        {
            arb_one(p);
            arb_one(p2);
            arb_set_arf(q, x);
            arb_set_round(q, q, prec);

            for (k = 0; k < r; k++)
            {
                arb_mul(q2, q, q, prec);
                arb_add(q2, p2, q2, prec);
                arb_sqrt(q2, q2, prec);
                arb_add(q, q, q2, prec);
            }
        }

        arb_div(p, p, q, prec);

        mag_set(err, arb_radref(p));
        inexact = arf_get_fmpz_fixed_si(res, arb_midref(p), -prec);
        mag_add_ui_2exp_si(err, err, inexact, -prec);

        arb_clear(p);
        arb_clear(p2);
        arb_clear(q);
        arb_clear(q2);
    }
}

typedef struct
{
    fmpz * s;
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
    fmpz * s;
    fmpz * u;
    slong r;
    fmpz_t P, Q;

    fmpz_init(P);
    fmpz_init(Q);

    s = work->s + iter;
    u = work->u + iter;
    r = work->r[iter];
    wp = work->wp;

    /* Binary splitting (+1 fixed-point ulp truncation error). */
    mag = fmpz_bits(u) - r;

    N = bs_num_terms(mag, wp);

    if (N != 0)
    {
        _arb_atan_sum_bs_powtab(P, Q, Qexp, u, r, N);

        /* multiply by u/2^r */
        fmpz_mul(P, P, u);
        *Qexp += r;

        /* T = T / Q  (+1 fixed-point ulp error). */
        if (*Qexp >= wp)
        {
            fmpz_tdiv_q_2exp(P, P, *Qexp - wp);
            fmpz_tdiv_q(P, P, Q);
        }
        else
        {
            fmpz_mul_2exp(P, P, wp - *Qexp);
            fmpz_tdiv_q(P, P, Q);
        }

        fmpz_add(s, s, P);
    }

    /* add u/2^r */
    fmpz_mul_2exp(Q, u, wp - r);
    fmpz_add(s, s, Q);

    fmpz_clear(P);
    fmpz_clear(Q);
}

void
arb_atan_arf_bb(arb_t z, const arf_t x, slong prec)
{
    slong iter, bits, r, mag, q, wp, N;
    slong argred_bits, start_bits;
    flint_bitcnt_t Qexp[1];
    int inverse;
    mag_t inp_err;
    fmpz_t s, t, u, P, Q, err;

    if (arf_is_zero(x))
    {
        arb_zero(z);
        return;
    }

    if (arf_is_special(x))
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }

    if (ARF_SGNBIT(x))
    {
        arf_t y;
        arf_init_neg_shallow(y, x);
        arb_atan_arf_bb(z, y, prec);
        arb_neg(z, z);
        return;
    }

    mag = arf_abs_bound_lt_2exp_si(x);

    /* We assume that this function only gets called with something
       reasonable as input (huge/tiny input will be handled by
       the main atan wrapper). */
    if (FLINT_ABS(mag) > 2 * prec + 100)
    {
        flint_throw(FLINT_ERROR, "arb_atan_arf_bb: unexpectedly large/small input\n");
    }

    /* approximate by x - x^3 / 3 or pi/2 - 1/x + (1/3)/x^3 */
    if (mag < -prec / 4 - 2 || (mag-1) > prec / 5 + 3)
    {
        arb_t t;
        arb_init(t);
        arb_set_arf(t, x);

        if (mag < 0)
        {
            arb_mul(t, t, t, prec);
            arb_mul_arf(t, t, x, prec);
            arb_div_ui(t, t, 3, prec);
            arb_sub_arf(t, t, x, prec);
            arb_neg(z, t);
            /* error is bounded by x^5 */
            mag_add_ui_2exp_si(arb_radref(z), arb_radref(z), 1, 5 * mag);
        }
        else
        {
            arb_ui_div(t, 1, t, prec);
            arb_mul(z, t, t, prec);
            arb_mul(z, z, t, prec);
            arb_div_ui(z, z, 3, prec);
            arb_sub(z, t, z, prec);

            arb_const_pi(t, prec + 2);
            arb_mul_2exp_si(t, t, -1);

            arb_sub(z, t, z, prec);
            /* error is bounded by 1/x^5, and 1/x <= 2^(1-mag) */
            mag_add_ui_2exp_si(arb_radref(z), arb_radref(z), 1, 5 * (1-mag));
        }

        arb_clear(t);
        return;
    }

    argred_bits = 8;
    start_bits = 16;

    /* Argument reduction q times. */
    q = FLINT_MAX(0, argred_bits - FLINT_ABS(mag));

    /* Determine working precision. */
    wp = prec + 10 + 2 * q + 2 * FLINT_BIT_COUNT(prec);
    if (mag < 0)
        wp += (-mag);

    fmpz_init(s);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(Q);
    fmpz_init(P);
    fmpz_init(err);  /* in fixed-point ulp */
    mag_init(inp_err); /* absolute error */

    arb_atan_bb_reduce(t, inp_err, x, mag, q, wp);
    inverse = mag > 0; /* todo: compute in function, or pass to it */

    /* s = 0, t = x */

    if (flint_get_num_available_threads() == 1 || prec >= 1000000000)
    {
        for (iter = 0, bits = start_bits; !fmpz_is_zero(t);
            iter++, bits *= 2)
        {
            /* Extract bits. */
            r = FLINT_MIN(bits, wp);
            fmpz_tdiv_q_2exp(u, t, wp - r);

            if (!fmpz_is_zero(u))
            {
                /* Binary splitting (+1 fixed-point ulp truncation error). */
                mag = fmpz_bits(u) - r;

                N = bs_num_terms(mag, wp);

                if (N != 0)
                {
                    _arb_atan_sum_bs_powtab(P, Q, Qexp, u, r, N);

                    /* multiply by u/2^r */
                    fmpz_mul(P, P, u);
                    *Qexp += r;

                    /* T = T / Q  (+1 fixed-point ulp error). */
                    if (*Qexp >= wp)
                    {
                        fmpz_tdiv_q_2exp(P, P, *Qexp - wp);
                        fmpz_tdiv_q(P, P, Q);
                    }
                    else
                    {
                        fmpz_mul_2exp(P, P, wp - *Qexp);
                        fmpz_tdiv_q(P, P, Q);
                    }

                    fmpz_add(s, s, P);
                }

                /* add u/2^r */
                fmpz_mul_2exp(Q, u, wp - r);
                fmpz_add(s, s, Q);

                /* 1 ulp from the division,
                   1 ulp from truncating the Taylor series */
                fmpz_add_ui(err, err, 2);
            }

            /* atan(t) = atan(u/2^r) + atan((t 2^r - u)/(2^r + u t)) */
            fmpz_mul_2exp(P, t, r);
            fmpz_mul_2exp(Q, u, wp);
            fmpz_sub(P, P, Q);

            fmpz_one(Q);
            fmpz_mul_2exp(Q, Q, r + wp);
            fmpz_addmul(Q, t, u);

            fmpz_mul_2exp(P, P, wp);
            fmpz_tdiv_q(t, P, Q);

            /* 1 ulp error from the division */
            fmpz_add_ui(err, err, 1);
        }
    }
    else
    {
        fmpz * ws;
        fmpz * us;
        slong * rs;
        slong num = 0;

        ws = _fmpz_vec_init(FLINT_BITS);
        us = _fmpz_vec_init(FLINT_BITS);
        rs = flint_malloc(sizeof(slong) * FLINT_BITS);

        /* Bit-burst loop. */
        for (iter = 0, bits = start_bits; !fmpz_is_zero(t);
            iter++, bits *= 2)
        {
            /* Extract bits. */
            r = FLINT_MIN(bits, wp);
            fmpz_tdiv_q_2exp(u, t, wp - r);

            if (!fmpz_is_zero(u))
            {
                fmpz_set(us + num, u);
                rs[num] = r;
                num++;

                /* 1 ulp from the division,
                   1 ulp from truncating the Taylor series */
                fmpz_add_ui(err, err, 2);
            }

            /* atan(t) = atan(u/2^r) + atan((t 2^r - u)/(2^r + u t)) */
            fmpz_mul_2exp(P, t, r);
            fmpz_mul_2exp(Q, u, wp);
            fmpz_sub(P, P, Q);

            fmpz_one(Q);
            fmpz_mul_2exp(Q, Q, r + wp);
            fmpz_addmul(Q, t, u);

            fmpz_mul_2exp(P, P, wp);
            fmpz_tdiv_q(t, P, Q);

            /* 1 ulp error from the division */
            fmpz_add_ui(err, err, 1);
        }

        /* todo: only allocate as many temporaries as threads,
           reducing memory */
        {
            work_t work;

            work.s = ws;
            work.u = us;
            work.r = rs;
            work.wp = wp;

            flint_parallel_do((do_func_t) worker, &work, num, -1, FLINT_PARALLEL_STRIDED);
        }

        for (iter = 0; iter < num; iter++)
            fmpz_add(s, s, ws + iter);

        _fmpz_vec_clear(ws, FLINT_BITS);
        _fmpz_vec_clear(us, FLINT_BITS);
        flint_free(rs);
    }

    /* add both err and inp_err */
    arf_set_fmpz(arb_midref(z), s);
    mag_set_fmpz(arb_radref(z), err);
    arb_mul_2exp_si(z, z, -wp);
    mag_add(arb_radref(z), arb_radref(z), inp_err);

    /* argument reduction: atan(x) = 2^q atan(x') */
    arb_mul_2exp_si(z, z, q);

    /* outmost argument reduction: atan(x) = +/-pi/2 - atan(1/x) */
    if (inverse)
    {
        arb_t pi2;
        arb_init(pi2);
        arb_const_pi(pi2, wp);
        arb_mul_2exp_si(pi2, pi2, -1);
        arb_sub(z, pi2, z, wp);
        arb_clear(pi2);
    }

    arb_set_round(z, z, prec);

    fmpz_clear(s);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(Q);
    fmpz_clear(P);
    fmpz_clear(err);
    mag_clear(inp_err);
}

