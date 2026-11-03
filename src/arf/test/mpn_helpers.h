/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARF_MPN_TEST_HELPERS_H
#define ARF_MPN_TEST_HELPERS_H

#include "test_helpers.h"
#include "fmpz.h"
#include "arf.h"

/* a random mantissa of n limbs with the top bit set, sometimes with
   structure (long runs of ones or zeros) */
static void
_arf_mpn_test_randtest_mantissa(nn_ptr d, slong n, flint_rand_t state)
{
    slong i;

    for (i = 0; i < n; i++)
        d[i] = n_randtest(state);
    d[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
}

/* a random sign-free rounding mode */
static arf_rnd_t
_arf_mpn_test_rnd_signfree(flint_rand_t state)
{
    switch (n_randint(state, 3))
    {
        case 0: return ARF_RND_DOWN;
        case 1: return ARF_RND_UP;
        default: return ARF_RND_NEAR;
    }
}

/* the reference rounding of the positive integer X to prec bits: Q, the
   exponent of the result (bit length) and the inexact flag */
static int
_arf_mpn_test_round_fmpz(fmpz_t q, slong * bits, const fmpz_t X, slong prec, arf_rnd_t rnd)
{
    slong L = fmpz_bits(X), tp = L - prec;
    fmpz_t r, half;
    int inexact, increment;

    if (tp <= 0)
    {
        fmpz_mul_2exp(q, X, -tp);
        *bits = prec;
        return 0;
    }

    fmpz_init(r);
    fmpz_init(half);
    fmpz_tdiv_q_2exp(q, X, tp);
    fmpz_tdiv_r_2exp(r, X, tp);
    inexact = !fmpz_is_zero(r);

    if (rnd == ARF_RND_DOWN)
        increment = 0;
    else if (rnd == ARF_RND_UP)
        increment = inexact;
    else
    {
        FLINT_ASSERT(rnd == ARF_RND_NEAR);
        int c;
        fmpz_one(half);
        fmpz_mul_2exp(half, half, tp - 1);
        c = fmpz_cmp(r, half);
        increment = (c > 0) || (c == 0 && fmpz_is_odd(q));
    }

    fmpz_add_ui(q, q, increment);
    *bits = prec;
    if (fmpz_bits(q) > prec)
    {
        /* a carry to the next power of two */
        fmpz_tdiv_q_2exp(q, q, 1);
        *bits = prec + 1;
    }

    fmpz_clear(r);
    fmpz_clear(half);
    return inexact;
}

/* whether a 2^ea == b 2^eb */
static int
_arf_mpn_test_equal_2exp(const fmpz_t a, slong ea, const fmpz_t b, slong eb)
{
    fmpz_t t;
    int r;
    fmpz_init(t);
    if (ea >= eb)
    {
        fmpz_mul_2exp(t, a, ea - eb);
        r = fmpz_equal(t, b);
    }
    else
    {
        fmpz_mul_2exp(t, b, eb - ea);
        r = fmpz_equal(t, a);
    }
    fmpz_clear(t);
    return r;
}

#endif
