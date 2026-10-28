/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "thread_support.h"
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

static double
now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

typedef void (*content_fn)(fmpz_t, const fmpz *, slong);

static double
timeit(content_fn f, fmpz_t res, const fmpz * vec, slong len)
{
    double t0, t, best = 1e300;
    slong reps = 1, i, k;

    for (k = 0; k < 3; k++)
    {
        while (1)
        {
            t0 = now();
            for (i = 0; i < reps; i++)
                f(res, vec, len);
            t = now() - t0;
            if (t > 0.01)
                break;
            reps *= 2;
        }
        best = FLINT_MIN(best, t / reps);
    }

    return best;
}

/* random signed integer of exactly bits bits (bits >= 1) */
static void
randbits(fmpz_t x, flint_rand_t state, flint_bitcnt_t bits)
{
    fmpz_randbits(x, state, bits);
    if (fmpz_is_zero(x))
        fmpz_one(x);
}

static const char * kinds[] = {
    "triv",     /* random entries */
    "small",    /* 16-bit content */
    "bal",      /* content of half the bits */
    "large",    /* 6-bit cofactors */
    "pow2",     /* content 2^k times a small odd number */
    "end",      /* half-size content except for the last entry */
    "mid",      /* half-size content except for the middle entry */
    "grow",     /* 24-bit content, sizes growing along the vector */
    "shrink",   /* 24-bit content, sizes shrinking along the vector */
    "sparse",   /* 90% zeros, 16-bit content */
    NULL
};

static void
gen(fmpz * vec, slong len, flint_bitcnt_t bits, const char * kind, flint_rand_t state)
{
    fmpz_t c, t;
    slong i;
    flint_bitcnt_t cb = 0, b;

    fmpz_init(c);
    fmpz_init(t);

    if (!strcmp(kind, "triv"))
        cb = 0;
    else if (!strcmp(kind, "small") || !strcmp(kind, "sparse"))
        cb = FLINT_MIN(16, bits - 1);
    else if (!strcmp(kind, "bal") || !strcmp(kind, "end") || !strcmp(kind, "mid"))
        cb = bits / 2;
    else if (!strcmp(kind, "large"))
        cb = bits > 7 ? bits - 6 : 1;
    else if (!strcmp(kind, "grow") || !strcmp(kind, "shrink"))
        cb = FLINT_MIN(24, bits - 1);

    if (!strcmp(kind, "pow2"))
    {
        fmpz_set_ui(c, n_randint(state, 1000) | 1);
        fmpz_mul_2exp(c, c, FLINT_MIN(bits / 2, 70));
        cb = fmpz_bits(c);
    }
    else if (cb == 0)
        fmpz_one(c);
    else
    {
        randbits(c, state, cb);
        fmpz_abs(c, c);
    }

    for (i = 0; i < len; i++)
    {
        b = bits;

        if (!strcmp(kind, "grow"))
            b = FLINT_MAX(cb + 1, bits * (i + 1) / len);
        else if (!strcmp(kind, "shrink"))
            b = FLINT_MAX(cb + 1, bits * (len - i) / len);

        b = (b > cb) ? b - cb : 1;
        randbits(t, state, b);
        fmpz_mul(vec + i, t, c);

        if (!strcmp(kind, "sparse") && n_randint(state, 10) != 0 && i != 0)
            fmpz_zero(vec + i);
    }

    if (!strcmp(kind, "end"))
        randbits(vec + len - 1, state, bits);
    else if (!strcmp(kind, "mid"))
        randbits(vec + len / 2, state, bits);

    fmpz_clear(c);
    fmpz_clear(t);
}

int main(int argc, char ** argv)
{
    static const slong lens[] = { 4, 8, 16, 32, 100, 1000, 10000, 0 };
    static const flint_bitcnt_t bitss[] = { 20, 62, 128, 640, 6400, 64000, 0 };
    flint_rand_t state;
    slong li, bi, ki, len;
    flint_bitcnt_t bits;
    const char * only = (argc > 1) ? argv[1] : NULL;
    double maxsize = (argc > 2) ? atof(argv[2]) : 3e7;
    slong threads = (argc > 3) ? atol(argv[3]) : 1;   /* for the new version */
    double t1, t2;
    fmpz_t r1, r2;

    flint_rand_init(state);
    fmpz_init(r1);
    fmpz_init(r2);

    /* the old version does not use threads */
    flint_set_num_threads(threads);

    flint_printf("%7s %6s %6s %12s %12s %7s\n", "kind", "len", "bits", "old", "new", "speedup");

    for (ki = 0; kinds[ki] != NULL; ki++)
    {
        if (only != NULL && strcmp(only, "all") && strcmp(only, kinds[ki]))
            continue;

        for (bi = 0; bitss[bi] != 0; bi++)
        {
            bits = bitss[bi];

            for (li = 0; lens[li] != 0; li++)
            {
                fmpz * vec;
                len = lens[li];

                if ((double) len * bits > maxsize)
                    continue;

                vec = _fmpz_vec_init(len);
                gen(vec, len, bits, kinds[ki], state);

                t1 = timeit(content_old, r1, vec, len);
                t2 = timeit(_fmpz_vec_content, r2, vec, len);

                if (!fmpz_equal(r1, r2))
                {
                    flint_printf("FAIL: %s %wd %wu\n", kinds[ki], len, bits);
                    flint_abort();
                }

                flint_printf("%7s %6wd %6wu %12.4g %12.4g %7.2f\n", kinds[ki], len, bits, t1, t2, t1 / t2);
                fflush(stdout);

                _fmpz_vec_clear(vec, len);
            }
        }
    }

    fmpz_clear(r1);
    fmpz_clear(r2);
    flint_rand_clear(state);
    return 0;
}
