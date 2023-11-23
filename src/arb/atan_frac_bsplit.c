/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

#ifdef __GNUC__
# define ceil __builtin_ceil
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

/* todo: arb arithmetic when sizes exceed prec */

static void
bsplit(fmpz_t p1, fmpz_t q1, fmpz_t r1,
        const fmpz_t p, const fmpz_t q,
        const fmpz_t ppow2, const fmpz_t qpow2,
        int alternate, slong a, slong b, int cont)
{
    if (b - a == 1)
    {
        if (a == 0)
        {
            fmpz_set(p1, p);
            fmpz_mul_ui(q1, q, 2*a+1);
            fmpz_mul_ui(r1, p, 2*a+1);
        }
        else
        {
            fmpz_set(p1, ppow2);
            fmpz_set(q1, qpow2);
            fmpz_mul_ui(q1, q1, 2*a+1);
            fmpz_mul_ui(r1, p1, 2*a+1);
        }

        if (alternate)
        {
            fmpz_neg(p1, p1);
            fmpz_neg(q1, q1);
        }
    }
    else
    {
        fmpz_t p2, q2, r2;
        slong m;

        m = (a + b) / 2;

        fmpz_init(p2);
        fmpz_init(q2);
        fmpz_init(r2);

        bsplit(p1, q1, r1, p, q, ppow2, qpow2, alternate, a, m, 1);
        bsplit(p2, q2, r2, p, q, ppow2, qpow2, alternate, m, b, cont);

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

#define LOG2 0.69314718055994530942

void
arb_atan_frac_bsplit(arb_t s, const fmpz_t p, const fmpz_t q, int hyperbolic, slong prec)
{
    fmpz_t P, Q, R, p2, q2;
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
        arb_atan_frac_bsplit(s, t, q, hyperbolic, prec);
        arb_neg(s, s);
        fmpz_clear(t);
        return;
    }

    fmpz_init(P);
    fmpz_init(Q);
    fmpz_init(R);
    fmpz_init(p2);
    fmpz_init(q2);
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
        logqp = logqp + (qexp - pexp) * LOG2 * (1.0 - 1e-12);
    }

    N = ceil((prec * (0.5 * LOG2) / logqp) * (1.0 + 1e-12));
    N = FLINT_MAX(N, 1);
    N = FLINT_MIN(N, 4 * prec);

    fmpz_mul(p2, p, p);
    fmpz_mul(q2, q, q);
    bsplit(P, Q, R, p, q, p2, q2, !hyperbolic, 0, N, 0);

    mag_set_fmpz(err, p);
    mag_div_fmpz(err, err, q);
    mag_geom_series(err, err, 2 * N + 1);
    mag_div_ui(err, err, 2 * N + 1);

    arb_fmpz_div_fmpz(s, P, Q, prec);
    arb_add_error_mag(s, err);

    fmpz_clear(p2);
    fmpz_clear(q2);
    fmpz_clear(P);
    fmpz_clear(Q);
    fmpz_clear(R);
    mag_clear(err);
}
