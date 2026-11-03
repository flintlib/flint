/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

/* Division, inversion and (reciprocal) square roots of exact balls at
   precisions where the Newton iterations of mp_real are used inside the
   mpn kernels (above the former Newton cutoffs of arf): the radius must
   be zero exactly when the result is exact, and the ball must contain the
   exact result. Uses its own deterministically seeded generator, since
   the inexact cases that were reported as exact are rare (about 1 in
   300 operands) and occur at specific operands. The precisions are 8 bits
   below a multiple of 64 bits, where the 8 guard bits of the Newton
   iterations end at a limb boundary (the case where the bug showed). */

#define PREC_ABOVE(cutoff) (64 * ((cutoff) / 64 + 2) - 8)

/* the former cutoffs of arf, in bits */
#define ARF_INV_NEWTON_CUTOFF 5000
#define ARF_DIV_NEWTON_CUTOFF 12000
#define ARF_SQRT_NEWTON_CUTOFF 100000
#define ARF_RSQRT_NEWTON_CUTOFF 4000

/* y random in (1, 2) with a full prec-bit mantissa */
static void
_rand_full(arf_t y, flint_rand_t state, slong prec)
{
    do {
        arf_urandom(y, state, prec, ARF_RND_DOWN);
        arf_add_ui(y, y, 1, prec, ARF_RND_DOWN);
    } while (arf_bits(y) < prec - 20);
}

TEST_FUNCTION_START(arb_div_newton, state)
{
    flint_rand_t st;
    arb_t y, z, t, u;
    arf_t a;
    slong i, prec;

    flint_rand_init(st);
    arb_init(y); arb_init(z); arb_init(t); arb_init(u);
    arf_init(a);

    /* inversion and reciprocal square root: 1/y and 1/sqrt(y) are not
       dyadic */
    prec = PREC_ABOVE(FLINT_MAX(ARF_INV_NEWTON_CUTOFF, ARF_RSQRT_NEWTON_CUTOFF));
    for (i = 0; i < 600; i++)
    {
        _rand_full(arb_midref(y), st, prec);
        mag_zero(arb_radref(y));

        arb_inv(z, y, prec);
        arb_mul(t, z, y, 3 * prec);
        if (arb_is_exact(z) || !arb_contains_si(t, 1))
            TEST_FUNCTION_FAIL("inv: i = %wd\ny = %{arb}\nz = %{arb}\n", i, y, z);

        arb_rsqrt(z, y, prec);
        arb_mul(t, z, z, 3 * prec);
        arb_mul(t, t, y, 3 * prec);
        if (arb_is_exact(z) || !arb_contains_si(t, 1))
            TEST_FUNCTION_FAIL("rsqrt: i = %wd\ny = %{arb}\nz = %{arb}\n", i, y, z);
    }

    /* division by a full-length divisor, exact and inexact */
    prec = PREC_ABOVE(ARF_DIV_NEWTON_CUTOFF);
    for (i = 0; i < 300; i++)
    {
        _rand_full(arb_midref(y), st, prec);
        mag_zero(arb_radref(y));

        arb_set_ui(u, 3);
        arb_div(z, u, y, prec);
        arb_mul(t, z, y, 3 * prec);
        if (arb_is_exact(z) || !arb_contains(t, u))
            TEST_FUNCTION_FAIL("div: i = %wd\ny = %{arb}\nz = %{arb}\n", i, y, z);

        /* u = y * a with a of at most prec bits: y / u is exact */
        arf_urandom(a, st, 1 + n_randint(st, prec), ARF_RND_DOWN);
        if (arf_is_zero(a))
            arf_one(a);
        arf_mul(arb_midref(u), arb_midref(y), a, ARF_PREC_EXACT, ARF_RND_DOWN);
        mag_zero(arb_radref(u));
        arb_div(z, u, y, prec);
        if (!arb_is_exact(z) || !arf_equal(arb_midref(z), a))
            TEST_FUNCTION_FAIL("div (exact): i = %wd\ny = %{arb}\nz = %{arb}\n", i, y, z);
    }

    /* square roots, exact and inexact */
    prec = PREC_ABOVE(ARF_SQRT_NEWTON_CUTOFF);
    for (i = 0; i < 10; i++)
    {
        _rand_full(arb_midref(y), st, prec);
        mag_zero(arb_radref(y));

        arb_sqrt(z, y, prec);
        arb_mul(t, z, z, 3 * prec);
        if (arb_is_exact(z) || !arb_contains(t, y))
            TEST_FUNCTION_FAIL("sqrt: i = %wd\n", i);

        arf_mul(arb_midref(u), arb_midref(y), arb_midref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
        mag_zero(arb_radref(u));
        arb_sqrt(z, u, prec);
        if (!arb_is_exact(z) || !arb_equal(z, y))
            TEST_FUNCTION_FAIL("sqrt (exact): i = %wd\n", i);
    }

    arb_clear(y); arb_clear(z); arb_clear(t); arb_clear(u);
    arf_clear(a);
    flint_rand_clear(st);

    TEST_FUNCTION_END(state);
}
