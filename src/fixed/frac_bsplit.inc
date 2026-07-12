/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Binary splitting for atan(p/q) and atanh(p/q)
   (arb_atan1_frac_bsplit, valid for p = 1: the leaves build the
   terms from q alone) and for log(1 + p/q) (arb_log1_frac_bsplit),
   carried from the newfastexp2 prototype; the alternating
   (non-hyperbolic) leaf signs were added here, everything else is
   verbatim.  Compared with arb_atan_frac_bsplit these fold pairs of
   terms at the leaves and strip two-powers from the denominators as
   they go.  TODO: reimplement natively with mpn arithmetic. */

#ifndef FIXED_FRAC_BSPLIT_INC
#define FIXED_FRAC_BSPLIT_INC

#include <math.h>
#include "mag.h"

#ifndef FIXED_LOG2_CONST
#define FIXED_LOG2_CONST 0.69314718055994530942
#endif

static void
atan1_bsplit(fmpz_t p1, fmpz_t q1, fmpz_t r1,
        const fmpz_t p, const fmpz_t q,
        const fmpz_t ppow2, const fmpz_t qpow2, const fmpz_t qpow4,
        int alternate, slong a, slong b, int cont)
{
    if (b - a == 1)
    {
        /* alternating: the terms carry (-1)^a absolutely, so the
           merges below need no sign logic */
        if (alternate && (a & 1))
            fmpz_set_si(p1, -1);
        else
            fmpz_one(p1);
        fmpz_mul_ui(q1, a == 0 ? q : qpow2, 2*a+1);
        fmpz_set_ui(r1, 2*a+1);
    }
    else if (b - a == 2 && a < (WORD(1) << (FLINT_BITS / 2 - 2)))
    {
        fmpz_mul_ui(p1, qpow2, 2 * a + 3);
        if (alternate)
        {
            fmpz_sub_ui(p1, p1, 2 * a + 1);
            if (a & 1)
                fmpz_neg(p1, p1);
        }
        else
            fmpz_add_ui(p1, p1, 2 * a + 1);

        if (a == 0)
        {
            fmpz_mul(q1, q, qpow2);
            fmpz_mul_ui(q1, q1, (2 * a + 1) * (2 * a + 3));
        }
        else
        {
            fmpz_mul_ui(q1, qpow4, (2 * a + 1) * (2 * a + 3));
        }

        fmpz_set_ui(r1, (2 * a + 1) * (2 * a + 3));
    }
    else
    {
        fmpz_t p2, q2, r2;
        slong m;

        m = (a + b) / 2;

        fmpz_init(p2);
        fmpz_init(q2);
        fmpz_init(r2);

        atan1_bsplit(p1, q1, r1, p, q, ppow2, qpow2, qpow4, alternate, a, m, 1);
        atan1_bsplit(p2, q2, r2, p, q, ppow2, qpow2, qpow4, alternate, m, b, cont);

        fmpz_mul(p1, p1, q2);
        fmpz_addmul(p1, r1, p2);
        fmpz_mul(q1, q1, q2);

        if (cont)
            fmpz_mul(r1, r1, r2);

        fmpz_clear(p2);
        fmpz_clear(q2);
        fmpz_clear(r2);
    }
}

static void
arb_atan1_frac_bsplit(arb_t s, const fmpz_t p, const fmpz_t q, int hyperbolic, slong prec)
{
    fmpz_t P, Q, R, p2, q2, q4;
    double logqp;
    slong N;
    mag_t err;

    if (fmpz_is_zero(p))
    {
        arb_zero(s);
        return;
    }

    if (fmpz_is_zero(q))
    {
        arb_indeterminate(s);
        return;
    }

    if (fmpz_sgn(p) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, p);
        arb_atan1_frac_bsplit(s, t, q, hyperbolic, prec);
        arb_neg(s, s);
        fmpz_clear(t);
        return;
    }

    fmpz_init(P);
    fmpz_init(Q);
    fmpz_init(R);
    fmpz_init(p2);
    fmpz_init(q2);
    fmpz_init(q4);
    mag_init(err);

    /* If N >= 1 and p/q <= 1/2, the error is bounded by (p/q)^(2N+1).
    For error <= 2^-prec, it is sufficient to pick
    N >= (1/2) * (prec * log(2) / log(q/p) - 1). */

    {
        slong qexp, pexp;
        double logp, logq;

        logp = fmpz_get_d_2exp(&pexp, p);
        logq = fmpz_get_d_2exp(&qexp, q);

        logqp = fabs(logq) / logp;
        logqp = mag_d_log_lower_bound(logqp) * (1.0 - 1e-12);
        logqp = logqp + (qexp - pexp) * FIXED_LOG2_CONST * (1.0 - 1e-12);
    }

    N = ceil((prec * (0.5 * FIXED_LOG2_CONST) / logqp) * (1.0 + 1e-12));
    N = FLINT_MAX(N, 1);
    N = FLINT_MIN(N, 4 * prec);
    if (N % 2)
        N++;

    fmpz_mul(p2, p, p);
    fmpz_mul(q2, q, q);
    fmpz_mul(q4, q2, q2);
    atan1_bsplit(P, Q, R, p, q, p2, q2, q4, !hyperbolic, 0, N, 0);

    mag_set_fmpz(err, p);
    mag_div_fmpz(err, err, q);
    mag_geom_series(err, err, 2 * N + 1);
    mag_div_ui(err, err, 2 * N + 1);

    arb_fmpz_div_fmpz(s, P, Q, prec);
    arb_add_error_mag(s, err);

    fmpz_clear(p2);
    fmpz_clear(q2);
    fmpz_clear(q4);
    fmpz_clear(P);
    fmpz_clear(Q);
    fmpz_clear(R);
    mag_clear(err);
}

#ifndef FIXED_FRAC_BSPLIT_NO_LOG1

static void
log1_bsplit(fmpz_t p1, fmpz_t q1, slong * q1exp, fmpz_t r1,
        const fmpz_t p, const fmpz_t q,
        int alternate, slong a, slong b, int cont)
{
    if (b - a == 1)
    {
        fmpz_set(p1, p);
        fmpz_set(q1, q);
        fmpz_mul_ui(q1, q1, a+1);
        fmpz_mul_ui(r1, p1, a+1);

        q1exp[0] = fmpz_val2(q1);
        fmpz_tdiv_q_2exp(q1, q1, q1exp[0]);

        if (alternate)
        {
            fmpz_neg(p1, p1);
            fmpz_neg(q1, q1);
        }
    }
    else
    {
        fmpz_t p2, q2, r2;
        slong q2exp[1];
        slong m;

        m = (a + b) / 2;

        fmpz_init(p2);
        fmpz_init(q2);
        fmpz_init(r2);

        log1_bsplit(p1, q1, q1exp, r1, p, q, alternate, a, m, 1);
        log1_bsplit(p2, q2, q2exp, r2, p, q, alternate, m, b, cont);

        fmpz_mul(p1, p1, q2);
        fmpz_mul_2exp(p1, p1, q2exp[0]);
        fmpz_addmul(p1, r1, p2);
        fmpz_mul(q1, q1, q2);
        q1exp[0] += q2exp[0];
        if (cont)
            fmpz_mul(r1, r1, r2);

        fmpz_clear(p2);
        fmpz_clear(q2);
        fmpz_clear(r2);
    }
}

static void
arb_log1_frac_bsplit(arb_t s, const fmpz_t p, const fmpz_t q, slong prec)
{
    fmpz_t P, Q, R;
    double logqp;
    slong N, Qexp;
    mag_t err;

    fmpz_init(P);
    fmpz_init(Q);
    fmpz_init(R);
    mag_init(err);

    {
        slong qexp, pexp;
        double logp, logq;

        logp = fmpz_get_d_2exp(&pexp, p);
        logq = fmpz_get_d_2exp(&qexp, q);

        logqp = fabs(logq) / logp;
        logqp = mag_d_log_lower_bound(logqp) * (1.0 - 1e-12);
        logqp = logqp + (qexp - pexp) * FIXED_LOG2_CONST * (1.0 - 1e-12);
    }

    N = ceil((prec * (FIXED_LOG2_CONST) / logqp) * (1.0 + 1e-12));
    N = FLINT_MAX(N, 1);
    N = FLINT_MIN(N, 4 * prec);

    log1_bsplit(P, Q, &Qexp, R, p, q, 1, 0, N, 0);

    mag_set_fmpz(err, p);
    mag_div_fmpz(err, err, q);
    mag_geom_series(err, err, N + 1);
    mag_div_ui(err, err, N + 1);

    arb_fmpz_div_fmpz(s, P, Q, prec);
    arb_mul_2exp_si(s, s, -Qexp);
    arb_add_error_mag(s, err);

    fmpz_clear(P);
    fmpz_clear(Q);
    fmpz_clear(R);
    mag_clear(err);
}

#endif /* FIXED_FRAC_BSPLIT_NO_LOG1 */

#endif /* FIXED_FRAC_BSPLIT_INC */
