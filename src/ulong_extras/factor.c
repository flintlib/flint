/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "factor_tables.h"

static int n_is_oddprime_small2(ulong n)
{
    ulong q = n / 2;
    ulong x = (q & 63);
    return (flint_odd_prime_lookup[q / 64] & (((uint64_t) 1) << x)) >> x;
}

#define APPEND_FACTOR(prime, exponent) \
    do \
    { \
        fac->p[fac->num] = (prime); \
        fac->exp[fac->num] = (exponent); \
        fac->num++; \
    } while (0)

#define REMOVE_ODD_TRIAL_PRIME(i) \
    do { \
        ulong inv1 = flint_trial_primes_inv1[i]; \
        ulong inv2 = flint_trial_primes_inv2[i]; \
        if (n_divisible_odd_gm(x, inv1, inv2)) \
        { \
            exp = 1; \
            x *= inv1; \
            while ((n_divisible_odd_gm(x, inv1, inv2))) \
            { \
                exp++; \
                x *= inv1; \
            } \
            APPEND_FACTOR(flint_trial_primes[i], exp); \
            if (x == 1) \
                return; \
            if (x < SMALL_ODDPRIME_LIMIT && n_is_oddprime_small2(x)) \
            { \
                APPEND_FACTOR(x, 1); \
                return; \
            } \
        } \
    } while (0)

/* Currently P+1 does not seem to save anything. */
#define USE_PP1 0

static void n_factor_after_trial(n_factor_t * factors, ulong n, ulong cutoff)
{
    ulong factor_arr[FLINT_MAX_FACTORS_IN_LIMB];
    ulong exp_arr[FLINT_MAX_FACTORS_IN_LIMB];
    ulong factors_left;
    ulong exp;
    ulong cofactor, factor;

    if (n == 1)
        return;

    if (n_is_prime_odd_no_trial(n))
    {
        n_factor_insert(factors, n, 1);
        return;
    }

    cofactor = n;

    factor_arr[0] = cofactor;
    factors_left = 1;
    exp_arr[0] = 1;

    while (factors_left > 0)
    {
        factor = factor_arr[factors_left - 1];

        if (factor >= cutoff)
        {
            if ((cofactor = n_factor_power235(&exp, factor)))
            {
                exp_arr[factors_left - 1] *= exp;
                factor_arr[factors_left - 1] = factor = cofactor;
            }

            if ((factor >= cutoff) && !n_is_prime_odd_no_trial(factor))
            {
                if ((
#if FLINT64
                    (factor < FLINT_FACTOR_ONE_LINE_MAX) &&
#endif
                     (cofactor = n_factor_one_line(factor, FLINT_FACTOR_ONE_LINE_ITERS)))
#if USE_PP1
                  || (cofactor = n_factor_pp1_wrapper2(factor))
#endif
                  || (cofactor = n_factor_SQUFOF(factor, FLINT_FACTOR_SQUFOF_ITERS)))
                {
                   exp_arr[factors_left] = exp_arr[factors_left - 1];
                       factor_arr[factors_left] = cofactor;
                       factor_arr[factors_left - 1] /= cofactor;
                       factors_left++;
                }
                else
                {
                    flint_throw(FLINT_ERROR, "Exception (n_factor). Failed to factor %wd.\n", n);
                }
            }
            else
            {
                n_factor_insert(factors, factor, exp_arr[factors_left - 1]);
                factors_left--;
            }
        }
        else
        {
            n_factor_insert(factors, factor, exp_arr[factors_left - 1]);
            factors_left--;
        }
    }
}

void
n_factor(n_factor_t * fac, ulong x, int FLINT_UNUSED(proved))
{
    slong i, exp;

    fac->num = 0;

    if (x <= 1)
        return;

    if (x % 2 == 0)
    {
        exp = flint_ctz(x);
        x >>= exp;
        APPEND_FACTOR(2, exp);
        if (x == 1)
            return;
    }

    if (x < SMALL_ODDPRIME_LIMIT && n_is_oddprime_small2(x))
    {
        APPEND_FACTOR(x, 1);
        return;
    }

    /* Simplify unrolling by 2 */
    FLINT_ASSERT(FLINT_FACTOR_TRIAL_PRIMES_BEFORE_PRIMALITY_TEST % 2 == 0);
    FLINT_ASSERT(FLINT_FACTOR_TRIAL_PRIMES % 2 == 0);

    i = 1;
    if (x % 3 == 0)
        REMOVE_ODD_TRIAL_PRIME(1);

    for (i = 2; i < FLINT_FACTOR_TRIAL_PRIMES_BEFORE_PRIMALITY_TEST; i += 2)
    {
        /* Manually unroll by two for performance */
        REMOVE_ODD_TRIAL_PRIME(i);
        REMOVE_ODD_TRIAL_PRIME(i + 1);

        if (x <= (ulong) flint_trial_primes[i+1] * (ulong) flint_trial_primes[i+1])
        {
            APPEND_FACTOR(x, 1);
            return;
        }
    }

    if (n_is_prime_odd_no_trial(x))
    {
        APPEND_FACTOR(x, 1);
        return;
    }

    for ( ; i < FLINT_FACTOR_TRIAL_PRIMES; i += 2)
    {
        REMOVE_ODD_TRIAL_PRIME(i);
        REMOVE_ODD_TRIAL_PRIME(i + 1);

        if (x <= (ulong) flint_trial_primes[i+1] * (ulong) flint_trial_primes[i+1])
        {
            APPEND_FACTOR(x, 1);
            return;
        }
    }

    n_factor_after_trial(fac, x, FLINT_FACTOR_TRIAL_CUTOFF);
}

