/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static void
bsplit(arb_t P, arb_t Q, const fmpz_t n, const fmpz_t a, const fmpz_t b, slong prec)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_sub(t, b, a);

    if (fmpz_sgn(t) <= 0)
    {
        arb_zero(P);
        arb_one(Q);
    }
    else if (fmpz_cmp_ui(t, 20) < 0)
    {
        slong steps, k;
        arb_t u;
        arb_init(u);

        arb_zero(P);
        arb_one(Q);

        steps = fmpz_get_si(t);

        for (k = steps - 1; k >= 0; k--)
        {
            fmpz_add_ui(t, a, k);

            arb_set_round_fmpz(u, t, prec);
            arb_pow_fmpz(u, u, n, prec);
            arb_addmul(P, Q, u, prec);

            if (!fmpz_is_zero(t))
                arb_mul_fmpz(Q, Q, t, prec);
        }

        arb_clear(u);
    }
    else
    {
        fmpz_t m;
        arb_t P1, Q2;

        fmpz_init(m);
        arb_init(P1);
        arb_init(Q2);

        fmpz_add(m, a, b);
        fmpz_tdiv_q_2exp(m, m, 1);

        bsplit(P1, Q, n, a, m, prec);
        bsplit(P, Q2, n, m, b, prec);

        arb_mul(Q, Q, Q2, prec);
        arb_addmul(P, P1, Q2, prec);

        fmpz_clear(m);
        arb_clear(P1);
        arb_clear(Q2);
    }

    fmpz_clear(t);
}

void
arb_bell_sum_bsplit(arb_t res, const fmpz_t n,
    const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec)
{
    if (fmpz_cmp(a, b) >= 0)
    {
        arb_zero(res);
    }
    else
    {
        slong wp;
        arb_t P, Q;

        wp = _fmpz_sub_small(b, a);
        wp = FLINT_BIT_COUNT(FLINT_ABS(wp));
        wp = prec + fmpz_bits(n) + fmpz_bits(a) + wp;

        arb_init(P);
        arb_init(Q);

        bsplit(P, Q, n, a, b, wp);
        arb_div(res, P, Q, wp);

        if (!fmpz_is_zero(a))
        {
            arb_gamma_fmpz(P, a, wp);
            arb_div(res, res, P, wp);
        }

        arb_set_round(res, res, prec);

        arb_clear(P);
        arb_clear(Q);
    }
}

