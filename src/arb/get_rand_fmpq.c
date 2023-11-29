/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb.h"

void
_arb_get_rand_fmpq(fmpz_t num, fmpz_t den, flint_rand_t state,
    const fmpz_t den_mult, const arb_t x)
{
    fmpz_t a, b, exp;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(exp);

    arb_get_interval_fmpz_2exp(a, b, exp, x);

    if (COEFF_IS_MPZ(*exp))
    {
        flint_throw(FLINT_ERROR, "exception: arb_get_rand_fmpq: too large exponent\n");
    }

    if (*exp >= 0)
    {
        fmpz_mul_2exp(a, a, *exp);
        fmpz_mul_2exp(b, b, *exp);
    }

    /* generate random integer in [a*den, b*den] */
    fmpz_mul(a, a, den_mult);
    fmpz_mul(b, b, den_mult);
    fmpz_add_ui(b, b, UWORD(1));
    fmpz_sub(b, b, a);

    /* return one endpoint with high probability (used for stress
       testing rounding) */
    if (n_randint(state, 6) == 0)
    {
        if (n_randint(state, 2))
            fmpz_zero(num);
        else
            fmpz_sub_ui(num, b, UWORD(1));
    }
    else
    {
        fmpz_randtest_mod(num, state, b);
    }

    fmpz_add(num, num, a);

    fmpz_set(den, den_mult);

    if (*exp < 0)
        fmpz_mul_2exp(den, den, -(*exp));

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(exp);
}

void
arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, slong bits)
{
    /* there is only one rational */
    if (arb_is_exact(x))
    {
        arf_get_fmpq(q, arb_midref(x));
        return;
    }

    /* pick a denominator */
    fmpz_randbits(fmpq_denref(q), state, n_randint(state, bits + 1));
    fmpz_abs(fmpq_denref(q), fmpq_denref(q));
    if (fmpz_is_zero(fmpq_denref(q)))
        fmpz_one(fmpq_denref(q));

    _arb_get_rand_fmpq(fmpq_numref(q), fmpq_denref(q), state, fmpq_denref(q), x);
    fmpq_canonicalise(q);
}

