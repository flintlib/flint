/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Contents and LCMs of factorial type, as in manipulating exponential-type
   power series over Q: the previous chained algorithms against
   _fmpz_vec_content_chained, _fmpz_vec_lcm and _fmpq_vec_get_fmpz_vec_fmpz */

#include <string.h>
#include <time.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "fmpq_vec.h"

/* previous _fmpz_vec_content_chained */
static void
content_chained_old(fmpz_t res, const fmpz * vec, slong len, const fmpz_t in)
{
    while (len > 0 && fmpz_is_zero(vec + 0))
    {
        len--;
        vec++;
    }

    while (len > 1 && fmpz_is_zero(vec + len - 1))
        len--;

    if (len == 0)
    {
        fmpz_abs(res, in);
        return;
    }

    if (len == 1)
    {
        fmpz_gcd(res, vec + 0, in);
        return;
    }

    if (fmpz_is_pm1(in) || fmpz_is_pm1(vec + 0) || fmpz_is_pm1(vec + len - 1))
    {
        fmpz_one(res);
        return;
    }

    fmpz_gcd3(res, vec + 0, vec + len - 1, in);
    vec++;
    len -= 2;

    while (len >= 2 && !fmpz_is_one(res))
    {
        fmpz_gcd3(res, vec + 0, vec + len - 1, res);
        vec++;
        len -= 2;
    }

    if (len != 0 && !fmpz_is_one(res))
        fmpz_gcd(res, res, vec + 0);
}

/* previous _fmpz_vec_lcm and _fmpq_vec_get_fmpz_vec_fmpz */
static void
lcm_old(fmpz_t res, const fmpz * vec, slong len)
{
    slong i;
    fmpz_one(res);
    for (i = 0; i < len && !fmpz_is_zero(res); i++)
        fmpz_lcm(res, res, vec + i);
    fmpz_abs(res, res);
}

static void
fmpq_vec_get_old(fmpz * num, fmpz_t den, const fmpq * a, slong len)
{
    slong i;

    fmpz_lcm(den, fmpq_denref(a + 0), fmpq_denref(a + 1));
    for (i = 2; i < len; i++)
        fmpz_lcm(den, den, fmpq_denref(a + i));

    for (i = 0; i < len; i++)
    {
        fmpz_divexact(num + i, den, fmpq_denref(a + i));
        fmpz_mul(num + i, num + i, fmpq_numref(a + i));
    }
}

static double
now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

#define TIMEIT(expr, out)                                       \
    do {                                                        \
        double __t0, __t, __best = 1e300;                       \
        slong __reps = 1, __i, __k;                             \
        for (__k = 0; __k < 3; __k++)                           \
        {                                                       \
            while (1)                                           \
            {                                                   \
                __t0 = now();                                   \
                for (__i = 0; __i < __reps; __i++)              \
                {                                               \
                    expr;                                       \
                }                                               \
                __t = now() - __t0;                             \
                if (__t > 0.02)                                 \
                    break;                                      \
                __reps *= 2;                                    \
            }                                                   \
            __best = FLINT_MIN(__best, __t / __reps);           \
        }                                                       \
        (out) = __best;                                         \
    } while (0)

/* the series sum m_k x^k / k!, k < n, with integer numerator vector num and
   denominator den = (n-1)! (not canonical: as obtained before
   canonicalisation) */
static void
exp_type(fmpz * num, fmpz_t den, slong n, flint_bitcnt_t mbits, flint_rand_t state)
{
    slong k;
    fmpz_t f;

    fmpz_init(f);
    fmpz_fac_ui(den, n - 1);

    /* num[k] = (n-1)! / k! * m_k, from the top */
    fmpz_one(f);
    for (k = n - 1; k >= 0; k--)
    {
        if (mbits == 0)
            fmpz_set(num + k, f);
        else
        {
            fmpz_randbits(num + k, state, mbits);
            if (fmpz_is_zero(num + k))
                fmpz_one(num + k);
            fmpz_mul(num + k, num + k, f);
        }
        if (k > 0)
            fmpz_mul_ui(f, f, k);
    }

    fmpz_clear(f);
}

static void
report(const char * what, slong n, double t1, double t2, const fmpz_t r)
{
    flint_printf("%-34s %5wd %11.4g %11.4g %7.2f   (result %wu bits)\n",
        what, n, t1, t2, t1 / t2, fmpz_bits(r));
    fflush(stdout);
}

int main(void)
{
    static const slong ns[] = { 10, 30, 100, 300, 1000, 0 };
    flint_rand_t state;
    slong ni, n, k;
    double t1, t2;
    fmpz_t r1, r2, den, den2, t;

    flint_rand_init(state);
    fmpz_init(r1);
    fmpz_init(r2);
    fmpz_init(den);
    fmpz_init(den2);
    fmpz_init(t);

    flint_printf("%-34s %5s %11s %11s %7s\n", "case", "n", "old", "new", "speedup");

    for (ni = 0; ns[ni] != 0; ni++)
    {
        fmpz * a, * b, * c;
        fmpq * q;
        n = ns[ni];

        a = _fmpz_vec_init(n);
        b = _fmpz_vec_init(n);
        c = _fmpz_vec_init(n);
        q = _fmpq_vec_init(n);

        /* contents ********************************************************/

        /* exp(x) * exp(x) before canonicalisation: numerators of the
           product with denominator (n-1)!^2 */
        exp_type(a, den, n, 0, state);
        _fmpz_poly_mullow(c, a, n, a, n, n);
        fmpz_mul(den2, den, den);
        TIMEIT(content_chained_old(r1, c, n, den2), t1);
        TIMEIT(_fmpz_vec_content_chained(r2, c, n, den2), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("content: exp*exp, chained with den", n, t1, t2, r2);

        /* f g with f, g of the form sum m_k x^k / k! (10-bit m_k) */
        exp_type(a, den, n, 10, state);
        exp_type(b, den, n, 10, state);
        _fmpz_poly_mullow(c, a, n, b, n, n);
        TIMEIT(content_chained_old(r1, c, n, den2), t1);
        TIMEIT(_fmpz_vec_content_chained(r2, c, n, den2), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("content: f*g, chained with den", n, t1, t2, r2);

        /* f + g with different exp-type series, numerators over (n-1)! */
        _fmpz_vec_add(c, a, b, n);
        TIMEIT(content_chained_old(r1, c, n, den), t1);
        TIMEIT(_fmpz_vec_content_chained(r2, c, n, den), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("content: f+g, chained with den", n, t1, t2, r2);

        /* the derivative of exp-type f: numerators k m_k (n-1)!/k!, k >= 1 */
        for (k = 1; k < n; k++)
            fmpz_mul_ui(c + k - 1, a + k, k);
        TIMEIT(content_chained_old(r1, c, n - 1, den), t1);
        TIMEIT(_fmpz_vec_content_chained(r2, c, n - 1, den), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("content: f', chained with den", n, t1, t2, r2);

        /* (n-1)! times random 20-bit entries: factorial content */
        fmpz_fac_ui(t, n - 1);
        for (k = 0; k < n; k++)
        {
            fmpz_randbits(c + k, state, 20);
            fmpz_mul(c + k, c + k, t);
        }
        TIMEIT(content_chained_old(r1, c, n, t), t1);
        TIMEIT(_fmpz_vec_content_chained(r2, c, n, t), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("content: (n-1)! r_k, chained with (n-1)!", n, t1, t2, r2);

        /* k! (n-1-k)! r_k: binomial type */
        for (k = 0; k < n; k++)
        {
            fmpz_fac_ui(t, k);
            fmpz_fac_ui(c + k, n - 1 - k);
            fmpz_mul(c + k, c + k, t);
            fmpz_mul_ui(c + k, c + k, 1 + n_randint(state, 1000));
        }
        TIMEIT(content_chained_old(r1, c, n, den), t1);
        TIMEIT(_fmpz_vec_content_chained(r2, c, n, den), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("content: k!(n-1-k)! r_k, with (n-1)!", n, t1, t2, r2);

        /* lcms ************************************************************/

        /* denominators k!, ascending and descending */
        for (k = 0; k < n; k++)
            fmpz_fac_ui(c + k, k);
        TIMEIT(lcm_old(r1, c, n), t1);
        TIMEIT(_fmpz_vec_lcm(r2, c, n), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("lcm: k!, ascending", n, t1, t2, r2);

        for (k = 0; k < n; k++)
            fmpz_fac_ui(c + k, n - 1 - k);
        TIMEIT(lcm_old(r1, c, n), t1);
        TIMEIT(_fmpz_vec_lcm(r2, c, n), t2);
        if (!fmpz_equal(r1, r2)) flint_abort();
        report("lcm: k!, descending", n, t1, t2, r2);

        /* sum x^k / k! as an fmpq vector, and sum m_k x^k / k! (canonical,
           so the denominators are divisors of k!) */
        for (k = 0; k < n; k++)
        {
            fmpz_one(fmpq_numref(q + k));
            fmpz_fac_ui(fmpq_denref(q + k), k);
        }
        TIMEIT(fmpq_vec_get_old(a, r1, q, n), t1);
        TIMEIT(_fmpq_vec_get_fmpz_vec_fmpz(b, r2, q, n), t2);
        if (!fmpz_equal(r1, r2) || !_fmpz_vec_equal(a, b, n)) flint_abort();
        report("fmpq -> fmpz: 1/k!", n, t1, t2, r2);

        for (k = 0; k < n; k++)
        {
            fmpz_randbits(fmpq_numref(q + k), state, 20);
            fmpz_fac_ui(fmpq_denref(q + k), k);
            fmpq_canonicalise(q + k);
        }
        TIMEIT(fmpq_vec_get_old(a, r1, q, n), t1);
        TIMEIT(_fmpq_vec_get_fmpz_vec_fmpz(b, r2, q, n), t2);
        if (!fmpz_equal(r1, r2) || !_fmpz_vec_equal(a, b, n)) flint_abort();
        report("fmpq -> fmpz: m_k/k! (20-bit m_k)", n, t1, t2, r2);

        /* 1/(k+1)! - the integral of exp */
        for (k = 0; k < n; k++)
        {
            fmpz_one(fmpq_numref(q + k));
            fmpz_fac_ui(fmpq_denref(q + k), k + 1);
        }
        TIMEIT(fmpq_vec_get_old(a, r1, q, n), t1);
        TIMEIT(_fmpq_vec_get_fmpz_vec_fmpz(b, r2, q, n), t2);
        if (!fmpz_equal(r1, r2) || !_fmpz_vec_equal(a, b, n)) flint_abort();
        report("fmpq -> fmpz: 1/(k+1)!", n, t1, t2, r2);

        _fmpz_vec_clear(a, n);
        _fmpz_vec_clear(b, n);
        _fmpz_vec_clear(c, n);
        _fmpq_vec_clear(q, n);
        flint_printf("\n");
    }

    fmpz_clear(r1);
    fmpz_clear(r2);
    fmpz_clear(den);
    fmpz_clear(den2);
    fmpz_clear(t);
    flint_rand_clear(state);
    return 0;
}
