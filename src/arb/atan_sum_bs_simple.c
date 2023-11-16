/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static void
bsplit(fmpz_t P, fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong a, slong b)
{
    if (b - a == 1)
    {
        fmpz_mul(P, x, x);
        fmpz_set(T, P);

        if (a % 2 == 0)
            fmpz_neg_ui(Q, 2 * a + 3);
        else
            fmpz_set_ui(Q, 2 * a + 3);

        *Qexp = 2 * r;
    }
    else
    {
        slong step, m;
        flint_bitcnt_t Q2exp[1];
        fmpz_t P2, Q2, T2;

        step = (b - a) / 2;
        m = a + step;

        fmpz_init(P2);
        fmpz_init(Q2);
        fmpz_init(T2);

        bsplit(P, T,  Q,  Qexp,  x, r, a, m);
        bsplit(P2, T2, Q2, Q2exp, x, r, m, b);

        fmpz_mul(T, T, Q2);
        fmpz_mul_2exp(T, T, *Q2exp);
        fmpz_mul(T2, T2, Q);
        fmpz_addmul(T, P, T2);
        fmpz_mul(P, P, P2);
        fmpz_mul(Q, Q, Q2);
        *Qexp = *Qexp + *Q2exp;

        fmpz_clear(P2);
        fmpz_clear(Q2);
        fmpz_clear(T2);
    }
}

void
_arb_atan_sum_bs_simple(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N)
{
    fmpz_t P;
    fmpz_init(P);
    bsplit(P, T, Q, Qexp, x, r, 0, N);
    fmpz_clear(P);
}

