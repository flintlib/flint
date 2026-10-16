/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fft_small.h"
#include "machine_vectors.h"

/* fail printer */
static void
mulmid_fail(const char* msg, ulong an, ulong bn, ulong lo, ulong hi)
{
    flint_printf("\nFAILED: %s\n", msg);
    flint_printf("an = %wu, bn = %wu, lo = %wu, hi = %wu\n", an, bn, lo, hi);
    fflush(stdout);
    flint_abort();
}

/*
    Check the limb window [lo, hi) of a*b returned by _mpn_ctx_mpn_mul_range.

    'full' must already hold _mpn_ctx_mpn_mul_range(R, full, lo, an+bn, a, an, b, bn),
    i.e. the same window but with the maximal top, so that the low limbs of the
    requested window can be matched against it without any high-truncation
    ambiguity.
*/
static void
check_window(mpn_ctx_t R, const ulong* a, ulong an, const ulong* b, ulong bn,
             const ulong* full, ulong lo, ulong hi, ulong* scratch)
{
    ulong zn = an + bn;

    _mpn_ctx_mpn_mul_range(R, scratch, lo, hi, a, an, b, bn);

    /* limbs [lo, min(hi, zn)) must agree with the maximal-top run */
    for (ulong k = lo; k < hi; k++)
    {
        ulong got = scratch[k - lo];
        ulong exp = (k < zn) ? full[k - lo] : UWORD(0);
        if (got != exp)
            mulmid_fail("window disagrees with full-top run / zero pad", an, bn, lo, hi);
    }
}

void test_mulmid(mpn_ctx_t R, ulong maxsize, ulong nreps, flint_rand_t state)
{
    ulong *a, *b, *c, *full, *scratch;
    fmpz_t fa, fb, fv, fvlo, fd, fdef, fbound;

    a       = FLINT_ARRAY_ALLOC(maxsize, ulong);
    b       = FLINT_ARRAY_ALLOC(maxsize, ulong);
    c       = FLINT_ARRAY_ALLOC(maxsize, ulong);
    full    = FLINT_ARRAY_ALLOC(maxsize + 2, ulong);
    scratch = FLINT_ARRAY_ALLOC(maxsize + 2, ulong);

    fmpz_init(fa);
    fmpz_init(fb);
    fmpz_init(fv);
    fmpz_init(fvlo);
    fmpz_init(fd);
    fmpz_init(fdef);
    fmpz_init(fbound);

    for (ulong rep = 0; rep < nreps; rep++)
    {
        ulong an = 2 + n_randint(state, maxsize - 4);
        ulong bn = 1 + n_randint(state, n_min(an, maxsize - an));
        ulong zn = an + bn;
        ulong lo, hi;

        for (ulong i = 0; i < maxsize; i++)
        {
            a[i] = n_randlimb(state);
            b[i] = n_randlimb(state);
        }

        /* a random low cut, biased to also hit 0 and the extremes */
        switch (n_randint(state, 4))
        {
            case 0:  lo = 0; break;
            case 1:  lo = n_randint(state, zn); break;
            case 2:  lo = (zn > 4) ? zn - 1 - n_randint(state, 4) : 0; break;
            default: lo = n_randint(state, n_max(1, zn/2)); break;
        }

        /* exact product and exact integers */
        mpn_mul(c, a, an, b, bn);
        fmpz_set_ui_array(fa, a, an);
        fmpz_set_ui_array(fb, b, bn);
        fmpz_mul(fv, fa, fb);

        /* maximal-top window [lo, zn): a lower approximation of floor(V/2^(64 lo)) */
        _mpn_ctx_mpn_mul_range(R, full, lo, zn, a, an, b, bn);

        /* --- lower-bound + deficit-bound, as full integers (no truncation) --- */
        fmpz_fdiv_q_2exp(fvlo, fv, FLINT_BITS*lo);          /* floor(V / 2^(64 lo)) */
        fmpz_set_ui_array(fd, full, zn - lo);               /* our approximation    */
        fmpz_sub(fdef, fvlo, fd);                           /* deficit = exact - ours */

        if (fmpz_sgn(fdef) < 0)
            mulmid_fail("not a lower approximation", an, bn, lo, zn);

        /* deficit must be <= min(an, bn, lo) * 2^64 */
        fmpz_set_ui(fbound, n_min(n_min(an, bn), lo));
        fmpz_mul_2exp(fbound, fbound, FLINT_BITS);
        if (fmpz_cmp(fdef, fbound) > 0)
            mulmid_fail("deficit exceeds min(an,bn,lo)*2^64", an, bn, lo, zn);

        /* with lo == 0 there is no dropped low part: the window is exact */
        if (lo == 0 && !fmpz_is_zero(fdef))
            mulmid_fail("low product must be exact (lo == 0)", an, bn, lo, zn);

        /* --- high-truncation consistency against the maximal-top run --- */
        for (int t = 0; t < 3; t++)
        {
            switch (n_randint(state, 4))
            {
                case 0:  hi = lo + 1; break;                       /* tiny window */
                case 1:  hi = lo + 1 + n_randint(state, zn - lo);  /* random      */
                         break;
                case 2:  hi = zn + 1 + n_randint(state, 3); break; /* past the top */
                default: hi = zn; break;
            }

            check_window(R, a, an, b, bn, full, lo, hi, scratch);
        }
    }

    fmpz_clear(fa);
    fmpz_clear(fb);
    fmpz_clear(fv);
    fmpz_clear(fvlo);
    fmpz_clear(fd);
    fmpz_clear(fdef);
    fmpz_clear(fbound);

    flint_free(a);
    flint_free(b);
    flint_free(c);
    flint_free(full);
    flint_free(scratch);
}

TEST_FUNCTION_START(_mpn_ctx_mpn_mul_range, state)
{
    {
        mpn_ctx_t R;
        mpn_ctx_init(R, UWORD(0x0003f00000000001));
        test_mulmid(R, 1000, 200 * flint_test_multiplier(), state);
        test_mulmid(R, 10000, 40 * flint_test_multiplier(), state);
        test_mulmid(R, 50000, 4 * flint_test_multiplier(), state);

        /* exercise the threading paths */
        flint_set_num_threads(1 + n_randint(state, 8));
        test_mulmid(R, 20000, 10 * flint_test_multiplier(), state);

        mpn_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}
