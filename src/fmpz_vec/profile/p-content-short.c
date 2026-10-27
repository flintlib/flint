/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Contents of short vectors (lengths 1 to 10), as in number field and
   short fmpq_poly arithmetic: the previous algorithms against
   _fmpz_vec_content and _fmpz_vec_content_chained (with a denominator, as
   in _fmpq_poly_canonicalise) */

#include <string.h>
#include <time.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

/* The previous _fmpz_vec_content (chained gcds from both ends), kept for
   comparison. */
static void
content_old(fmpz_t res, const fmpz * vec, slong len)
{
    while (len > 0 && fmpz_is_zero(vec + 0))
    {
        len--;
        vec++;
    }

    while (len > 1 && fmpz_is_zero(vec + len - 1))
        len--;

    if (len <= 2)
    {
        if (len == 0)
            fmpz_zero(res);
        else if (len == 1)
            fmpz_abs(res, vec + 0);
        else
            fmpz_gcd(res, vec + 0, vec + 1);
        return;
    }

    if (fmpz_is_pm1(vec + 0) || fmpz_is_pm1(vec + len - 1))
    {
        fmpz_one(res);
        return;
    }

    fmpz_gcd3(res, vec + 0, vec + 1, vec + len - 1);
    vec += 2;
    len -= 3;

    while (len >= 2 && !fmpz_is_one(res))
    {
        fmpz_gcd3(res, vec + 0, vec + len - 1, res);
        vec++;
        len -= 2;
    }

    if (len != 0 && !fmpz_is_one(res))
        fmpz_gcd(res, res, vec + 0);
}

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
                if (__t > 0.002)                                 \
                    break;                                      \
                __reps *= 2;                                    \
            }                                                   \
            __best = FLINT_MIN(__best, __t / __reps);           \
        }                                                       \
        (out) = __best;                                         \
    } while (0)


static const char * shapes[] = {
    "random",   /* random numerators and denominator */
    "common",   /* numerators and denominator share a factor of half the bits */
    "small",    /* ... a 12-bit factor */
    NULL
};

/* numerators f[0..len) and denominator d of the given bits */
static void
gen(fmpz * f, fmpz_t d, slong len, flint_bitcnt_t bits, const char * shape, flint_rand_t state)
{
    fmpz_t c;
    flint_bitcnt_t cb;
    slong i;

    fmpz_init(c);
    cb = !strcmp(shape, "common") ? bits / 2 : !strcmp(shape, "small") ? FLINT_MIN(12, bits - 1) : 0;
    if (cb == 0)
        fmpz_one(c);
    else
    {
        fmpz_randbits(c, state, cb);
        fmpz_abs(c, c);
        if (fmpz_is_zero(c))
            fmpz_one(c);
    }

    for (i = 0; i < len; i++)
    {
        fmpz_randbits(f + i, state, bits - cb);
        if (fmpz_is_zero(f + i))
            fmpz_one(f + i);
        fmpz_mul(f + i, f + i, c);
    }
    fmpz_randbits(d, state, bits - cb);
    fmpz_abs(d, d);
    if (fmpz_is_zero(d))
        fmpz_one(d);
    fmpz_mul(d, d, c);

    fmpz_clear(c);
}

int main(void)
{
    static const flint_bitcnt_t bitss[] = { 30, 62, 128, 320, 1000, 0 };
    flint_rand_t state;
    slong len, si, bi, trial;
    flint_bitcnt_t bits;
    fmpz_t r1, r2, d;
    fmpz * f;
    double t1, t2, s1, s2;
    int chained;

    flint_rand_init(state);
    fmpz_init(r1);
    fmpz_init(r2);
    fmpz_init(d);
    f = _fmpz_vec_init(10);

    for (chained = 1; chained >= 0; chained--)
    {
        flint_printf("\n%s: speedup (old time / new time), averaged over 20 vectors\n",
            chained ? "_fmpz_vec_content_chained with the denominator" : "_fmpz_vec_content");
        flint_printf("%7s %5s ", "shape", "bits");
        for (len = 1; len <= 10; len++)
            flint_printf("%6wd", len);
        flint_printf("   (new ns, len 1..10)\n");

        for (si = 0; shapes[si] != NULL; si++)
        {
            for (bi = 0; bitss[bi] != 0; bi++)
            {
                double ns[11];
                bits = bitss[bi];
                flint_printf("%7s %5wu ", shapes[si], bits);

                for (len = 1; len <= 10; len++)
                {
                    s1 = s2 = 0;
                    for (trial = 0; trial < 20; trial++)
                    {
                        gen(f, d, len, bits, shapes[si], state);
                        if (chained)
                        {
                            TIMEIT(content_chained_old(r1, f, len, d), t1);
                            TIMEIT(_fmpz_vec_content_chained(r2, f, len, d), t2);
                        }
                        else
                        {
                            TIMEIT(content_old(r1, f, len), t1);
                            TIMEIT(_fmpz_vec_content(r2, f, len), t2);
                        }
                        if (!fmpz_equal(r1, r2))
                            flint_abort();
                        s1 += t1;
                        s2 += t2;
                    }
                    ns[len] = s2 / 20 * 1e9;
                    flint_printf("%6.2f", s1 / s2);
                    fflush(stdout);
                }
                flint_printf("   (%.0f .. %.0f)\n", ns[1], ns[10]);
            }
        }
    }

    _fmpz_vec_clear(f, 10);
    fmpz_clear(r1);
    fmpz_clear(r2);
    fmpz_clear(d);
    flint_rand_clear(state);
    return 0;
}
