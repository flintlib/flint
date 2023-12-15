/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mpn_extras.h"

/* See verify_taylor.py for code to generate tables and
   proof of correctness */

#define TMP_ALLOC_LIMBS(size) TMP_ALLOC((size) * sizeof(mp_limb_t))

#define FACTORIAL_TAB_SIZE 288

FLINT_DLL extern const mp_limb_t factorial_tab_numer[FACTORIAL_TAB_SIZE];
FLINT_DLL extern const mp_limb_t factorial_tab_denom[FACTORIAL_TAB_SIZE];

void _arb_sin_cos_taylor_rs(mp_ptr ysin, mp_ptr ycos,
    mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N,
    int sinonly, int alternating)
{
    mp_ptr s, t, xpow;
    mp_limb_t new_denom, old_denom, c;
    slong power, k, m;
    int cosorsin;

    TMP_INIT;
    TMP_START;

    if (2 * N >= FACTORIAL_TAB_SIZE - 1)
    {
        flint_throw(FLINT_ERROR, "_arb_sin_cos_taylor_rs: N too large!\n");
    }

    if (N <= 1)
    {
        if (N == 0)
        {
            flint_mpn_zero(ysin, xn);
            if (!sinonly) flint_mpn_zero(ycos, xn);
            error[0] = 0;
        }
        else if (N == 1)
        {
            flint_mpn_copyi(ysin, x, xn);
            if (!sinonly) flint_mpn_store(ycos, xn, LIMB_ONES);
            error[0] = 1;
        }
    }
    else
    {
        /* Choose m ~= sqrt(num_terms) (m must be even, >= 2) */
        m = 2;
        while (m * m < N)
            m += 2;

        /* todo: merge allocations */
        xpow = TMP_ALLOC_LIMBS((m + 1) * xn);
        s = TMP_ALLOC_LIMBS(xn + 2);
        t = TMP_ALLOC_LIMBS(2 * xn + 2);     /* todo: 1 limb too much? */

        /* higher index ---> */
        /*        | ---xn--- | */
        /* xpow = |  <temp>  | x^m | x^(m-1) | ... | x^2 | x | */

#define XPOW_WRITE(__k) (xpow + (m - (__k)) * xn)
#define XPOW_READ(__k) (xpow + (m - (__k) + 1) * xn)

        flint_mpn_sqr(XPOW_WRITE(1), x, xn);
        flint_mpn_sqr(XPOW_WRITE(2), XPOW_READ(1), xn);

        for (k = 4; k <= m; k += 2)
        {
            flint_mpn_mul_n(XPOW_WRITE(k - 1), XPOW_READ(k / 2), XPOW_READ(k / 2 - 1), xn);
            flint_mpn_sqr(XPOW_WRITE(k), XPOW_READ(k / 2), xn);
        }

        for (cosorsin = sinonly; cosorsin < 2; cosorsin++)
        {
            flint_mpn_zero(s, xn + 1);

            /* todo: skip one nonscalar multiplication (use x^m)
               when starting on x^0 */
            power = (N - 1) % m;

            for (k = N - 1; k >= 0; k--)
            {
                c = factorial_tab_numer[2 * k + cosorsin];
                new_denom = factorial_tab_denom[2 * k + cosorsin];
                old_denom = factorial_tab_denom[2 * k + cosorsin + 2];

                /* change denominators */
                if (new_denom != old_denom && k < N - 1)
                {
                    if (alternating && (k % 2 == 0))
                        s[xn] += old_denom;

                    mpn_divrem_1(s, 0, s, xn + 1, old_denom);

                    if (alternating && (k % 2 == 0))
                        s[xn] -= 1;
                }

                if (power == 0)
                {
                    /* add c * x^0 -- only top limb is affected */
                    if (alternating & k)
                        s[xn] -= c;
                    else
                        s[xn] += c;

                    /* Outer polynomial evaluation: multiply by x^m */
                    if (k != 0)
                    {
                        flint_mpn_mul(t, s, xn + 1, XPOW_READ(m), xn);
                        flint_mpn_copyi(s, t + xn, xn + 1);
                    }

                    power = m - 1;
                }
                else
                {
                    if (alternating & k)
                        s[xn] -= mpn_submul_1(s, XPOW_READ(power), xn, c);
                    else
                        s[xn] += mpn_addmul_1(s, XPOW_READ(power), xn, c);

                    power--;
                }
            }

            /* finally divide by denominator */
            if (cosorsin == 0)
            {
                mpn_divrem_1(t, 0, s, xn + 1, factorial_tab_denom[0]);

                /* perturb down to a number < 1 if necessary. note that this
                   does not invalidate the error bound: 1 - ulp is either
                   1 ulp too small or must be closer to the exact value */
                if (t[xn] == 0)
                    flint_mpn_copyi(ycos, t, xn);
                else
                    flint_mpn_store(ycos, xn, LIMB_ONES);
            }
            else
            {
                mpn_divrem_1(s, 0, s, xn + 1, factorial_tab_denom[0]);
                flint_mpn_mul(t, s, xn + 1, x, xn);
                flint_mpn_copyi(ysin, t + xn, xn);
            }
        }

        /* error bound (ulp) */
        error[0] = 2;
    }

    TMP_END;
}

