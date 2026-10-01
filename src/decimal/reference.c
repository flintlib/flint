/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "decimal.h"
#include "fmpq.h"
#include "gr.h"

/* Not performance-critical: optimize for size. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/*
    Reference implementation of correctly rounded conversion of a rational
    number to prec digits, using only fmpz arithmetic. Used as an oracle by
    the test code. Exponent limits of the context are applied through
    _decfloat_finalize.
*/
int
decfloat_set_round_fmpq_reference(decfloat_t res, const fmpq_t q, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz_t a, b, N, R, U, t, E;
    slong da, db, Es, k;
    int negative, up, cmp, sticky, status;

    if (fmpq_is_zero(q))
        return decfloat_zero(res, ctx);

    if (prec == DECIMAL_PREC_EXACT)
    {
        /* exact only if the denominator is 2^i 5^j */
        fmpz_t d, two, five;
        slong v2, v5;
        int ok;

        fmpz_init(d);
        fmpz_init_set_ui(two, 2);
        fmpz_init_set_ui(five, 5);
        fmpz_set(d, fmpq_denref(q));
        v2 = fmpz_remove(d, d, two);
        v5 = fmpz_remove(d, d, five);
        ok = fmpz_is_one(d);

        if (ok)
        {
            fmpz_t m, e;
            slong kk = FLINT_MAX(v2, v5);
            fmpz_init(m);
            fmpz_init(e);
            fmpz_ui_pow_ui(m, 10, kk);
            fmpz_mul(m, m, fmpq_numref(q));
            fmpz_divexact(m, m, fmpq_denref(q));
            fmpz_set_si(e, -kk);
            status = decfloat_set_round_fmpz_10exp_fmpz(res, m, e, DECIMAL_PREC_EXACT, rnd, ctx);
            fmpz_clear(m);
            fmpz_clear(e);
        }
        else
            status = GR_UNABLE;

        fmpz_clear(d);
        fmpz_clear(two);
        fmpz_clear(five);
        return status;
    }

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(N);
    fmpz_init(R);
    fmpz_init(U);
    fmpz_init(t);
    fmpz_init(E);

    fmpz_abs(a, fmpq_numref(q));
    fmpz_set(b, fmpq_denref(q));
    negative = fmpz_sgn(fmpq_numref(q)) < 0;

    /* E = floor(log10(a/b)) */
    da = fmpz_sizeinbase(a, 10);
    db = fmpz_sizeinbase(b, 10);
    Es = da - db;
    /* adjust so that 10^Es <= a/b < 10^(Es+1) (sizeinbase may be off by one) */
    for (;;)
    {
        int c;
        fmpz_ui_pow_ui(t, 10, FLINT_ABS(Es));
        if (Es >= 0)
        {
            fmpz_mul(t, t, b);
            c = fmpz_cmp(a, t);        /* a >= b 10^Es ? */
        }
        else
        {
            fmpz_mul(t, t, a);
            c = fmpz_cmp(t, b);        /* a 10^-Es >= b ? */
        }
        if (c < 0) { Es--; continue; }

        fmpz_ui_pow_ui(t, 10, FLINT_ABS(Es + 1));
        if (Es + 1 >= 0)
        {
            fmpz_mul(t, t, b);
            c = fmpz_cmp(a, t);        /* a < b 10^(Es+1) ? */
        }
        else
        {
            fmpz_mul(t, t, a);
            c = fmpz_cmp(t, b);
        }
        if (c >= 0) { Es++; continue; }
        break;
    }
    /* sanity: now 10^Es <= a/b < 10^(Es+1) */

    /* N = floor(a / (b 10^(Es - prec + 1))) */
    k = Es - prec + 1;

    if (k >= 0)
    {
        fmpz_ui_pow_ui(t, 10, k);
        fmpz_mul(t, t, b);
        fmpz_fdiv_qr(N, R, a, t);
    }
    else
    {
        fmpz_ui_pow_ui(t, 10, -k);
        fmpz_mul(t, t, a);
        fmpz_fdiv_qr(N, R, t, b);
        fmpz_set(t, b);
    }

    /* now t is the divisor, R the remainder: compare 2R with t */
    sticky = !fmpz_is_zero(R);
    fmpz_mul_2exp(R, R, 1);
    cmp = fmpz_cmp(R, t);

    up = 0;
    if (sticky)
    {
        switch (rnd)
        {
            case DECIMAL_RND_DOWN: up = 0; break;
            case DECIMAL_RND_UP: up = 1; break;
            case DECIMAL_RND_FLOOR: up = negative; break;
            case DECIMAL_RND_CEIL: up = !negative; break;
            case DECIMAL_RND_NEAR: up = (cmp > 0) || (cmp == 0 && fmpz_is_odd(N)); break;
            case DECIMAL_RND_NEAR_AWAY: up = (cmp >= 0); break;
            case DECIMAL_RND_NEAR_ZERO: up = (cmp > 0); break;
            default: flint_abort();
        }
    }

    if (up)
        fmpz_add_ui(N, N, 1);
    if (negative)
        fmpz_neg(N, N);

    fmpz_set_si(E, k);
    status = decfloat_set_round_fmpz_10exp_fmpz(res, N, E, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(N);
    fmpz_clear(R);
    fmpz_clear(U);
    fmpz_clear(t);
    fmpz_clear(E);

    return status;
}

POP_OPTIONS
