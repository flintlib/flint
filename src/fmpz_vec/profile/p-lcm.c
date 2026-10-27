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
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"

/* The previous implementations, kept for comparison */
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

    if (len < 1)
    {
        fmpz_one(den);
        return;
    }
    else if (len == 1)
    {
        fmpz_set(num + 0, fmpq_numref(a + 0));
        fmpz_set(den, fmpq_denref(a + 0));
        return;
    }

    fmpz_lcm(den, fmpq_denref(a + 0), fmpq_denref(a + 1));
    for (i = 2; i < len; i++)
        fmpz_lcm(den, den, fmpq_denref(a + i));

    if (fmpz_is_one(den))
    {
        for (i = 0; i < len; i++)
            fmpz_set(num + i, fmpq_numref(a + i));
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            fmpz_divexact(num + i, den, fmpq_denref(a + i));
            fmpz_mul(num + i, num + i, fmpq_numref(a + i));
        }
    }
}

static void
colwise_old(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)
{
    fmpz_t t, lcm;
    slong i, j;

    fmpz_init(t);
    fmpz_init(lcm);

    for (j = 0; j < mat->c; j++)
    {
        fmpz_set(lcm, fmpq_mat_entry_den(mat, 0, j));
        for (i = 1; i < mat->r; i++)
            fmpz_lcm(lcm, lcm, fmpq_mat_entry_den(mat, i, j));
        fmpz_set(den + j, lcm);
        if (fmpz_is_one(lcm))
        {
            for (i = 0; i < mat->r; i++)
                fmpz_set(fmpz_mat_entry(num, i, j), fmpq_mat_entry_num(mat, i, j));
        }
        else
        {
            for (i = 0; i < mat->r; i++)
            {
                fmpz_divexact(t, lcm, fmpq_mat_entry_den(mat, i, j));
                fmpz_mul(fmpz_mat_entry(num, i, j), fmpq_mat_entry_num(mat, i, j), t);
            }
        }
    }

    fmpz_clear(t);
    fmpz_clear(lcm);
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
                if (__t > 0.01)                                 \
                    break;                                      \
                __reps *= 2;                                    \
            }                                                   \
            __best = FLINT_MIN(__best, __t / __reps);           \
        }                                                       \
        (out) = __best;                                         \
    } while (0)

static const char * kinds[] = {
    "one",      /* all denominators 1 */
    "same",     /* all equal to a random D of the given bits */
    "mostly1",  /* 90% ones, else D */
    "div",      /* random divisors of D, a product of 12-bit primes */
    "small",    /* random 1 to 8 bit denominators */
    "harm",     /* d_i = i + 1 */
    "fact",     /* d_i = (i + 1)! (bits ignored) */
    "rand",     /* random denominators of the given bits */
    NULL
};

static void
gen_dens(fmpz * d, slong len, flint_bitcnt_t bits, const char * kind, flint_rand_t state)
{
    fmpz_t D, t;
    slong i, j, np;
    fmpz primes[64];

    fmpz_init(D);
    fmpz_init(t);

    fmpz_randbits(D, state, bits);
    fmpz_abs(D, D);
    if (fmpz_is_zero(D))
        fmpz_one(D);

    /* D as a product of 12-bit primes for "div" */
    np = FLINT_MAX(1, FLINT_MIN(64, bits / 12));
    for (j = 0; j < np; j++)
    {
        fmpz_init(primes + j);
        fmpz_set_ui(primes + j, n_nextprime(2048 + n_randint(state, 2048), 1));
    }

    for (i = 0; i < len; i++)
    {
        if (!strcmp(kind, "one"))
            fmpz_one(d + i);
        else if (!strcmp(kind, "same"))
            fmpz_set(d + i, D);
        else if (!strcmp(kind, "mostly1"))
        {
            if (n_randint(state, 10) == 0)
                fmpz_set(d + i, D);
            else
                fmpz_one(d + i);
        }
        else if (!strcmp(kind, "div"))
        {
            fmpz_one(d + i);
            for (j = 0; j < np; j++)
                if (n_randint(state, 3) != 0)
                    fmpz_mul(d + i, d + i, primes + j);
        }
        else if (!strcmp(kind, "small"))
            fmpz_set_ui(d + i, 1 + n_randint(state, 1 + n_randint(state, 256)));
        else if (!strcmp(kind, "harm"))
            fmpz_set_ui(d + i, i + 1);
        else if (!strcmp(kind, "fact"))
        {
            if (i == 0)
                fmpz_one(d + i);
            else
                fmpz_mul_ui(d + i, d + i - 1, i + 1);
        }
        else
        {
            fmpz_randbits(d + i, state, bits);
            fmpz_abs(d + i, d + i);
            if (fmpz_is_zero(d + i))
                fmpz_one(d + i);
        }
    }

    for (j = 0; j < np; j++)
        fmpz_clear(primes + j);
    fmpz_clear(D);
    fmpz_clear(t);
}

int main(int argc, char ** argv)
{
    static const slong lens[] = { 4, 16, 100, 1000, 10000, 0 };
    static const flint_bitcnt_t bitss[] = { 16, 62, 128, 640, 0 };
    const char * only = (argc > 1) ? argv[1] : NULL;
    int mode = (argc > 2) ? atoi(argv[2]) : 0;   /* 0: _fmpz_vec_lcm, 1: _fmpq_vec_get_fmpz_vec_fmpz, 2: fmpq_mat_get_fmpz_mat_colwise (len x 4) */
    flint_rand_t state;
    slong li, bi, ki, len, i;
    flint_bitcnt_t bits;
    double t1, t2;
    fmpz_t r1, r2;

    flint_rand_init(state);
    fmpz_init(r1);
    fmpz_init(r2);

    flint_printf("%8s %6s %5s %12s %12s %7s   (%s)\n", "kind", "len", "bits", "old", "new", "speedup",
        mode == 0 ? "_fmpz_vec_lcm" : mode == 1 ? "_fmpq_vec_get_fmpz_vec_fmpz" : "fmpq_mat_get_fmpz_mat_colwise, len x 4");

    for (ki = 0; kinds[ki] != NULL; ki++)
    {
        if (only != NULL && strcmp(only, "all") && strcmp(only, kinds[ki]))
            continue;

        for (bi = 0; bitss[bi] != 0; bi++)
        {
            bits = bitss[bi];

            if ((!strcmp(kinds[ki], "harm") || !strcmp(kinds[ki], "fact") || !strcmp(kinds[ki], "one")
                    || !strcmp(kinds[ki], "small")) && bi != 0)
                continue;

            for (li = 0; lens[li] != 0; li++)
            {
                fmpz * d;
                len = lens[li];

                /* keep the lcm (or the factorials) of moderate size */
                if (!strcmp(kinds[ki], "rand") && (double) len * bits > 3e5)
                    continue;
                if (!strcmp(kinds[ki], "fact") && len > 1000)
                    continue;
                if (!strcmp(kinds[ki], "harm") && len > 10000)
                    continue;

                d = _fmpz_vec_init(len);
                gen_dens(d, len, bits, kinds[ki], state);

                if (mode == 0)
                {
                    TIMEIT(lcm_old(r1, d, len), t1);
                    TIMEIT(_fmpz_vec_lcm(r2, d, len), t2);
                }
                else if (mode == 2)
                {
                    fmpq_mat_t A;
                    fmpz_mat_t N1, N2;
                    fmpz * D1 = _fmpz_vec_init(4), * D2 = _fmpz_vec_init(4);
                    slong j;
                    fmpq_mat_init(A, len, 4);
                    fmpz_mat_init(N1, len, 4);
                    fmpz_mat_init(N2, len, 4);
                    for (j = 0; j < 4; j++)
                    {
                        gen_dens(d, len, bits, kinds[ki], state);
                        for (i = 0; i < len; i++)
                        {
                            fmpz_randbits(fmpq_mat_entry_num(A, i, j), state, bits);
                            fmpz_set(fmpq_mat_entry_den(A, i, j), d + i);
                            fmpq_canonicalise(fmpq_mat_entry(A, i, j));
                        }
                    }
                    TIMEIT(colwise_old(N1, D1, A), t1);
                    TIMEIT(fmpq_mat_get_fmpz_mat_colwise(N2, D2, A), t2);
                    if (!fmpz_mat_equal(N1, N2) || !_fmpz_vec_equal(D1, D2, 4))
                    {
                        flint_printf("FAIL (colwise): %s %wd %wu\n", kinds[ki], len, bits);
                        flint_abort();
                    }
                    fmpz_set(r1, D1);
                    fmpz_set(r2, D2);
                    fmpq_mat_clear(A);
                    fmpz_mat_clear(N1);
                    fmpz_mat_clear(N2);
                    _fmpz_vec_clear(D1, 4);
                    _fmpz_vec_clear(D2, 4);
                }
                else
                {
                    fmpq * a = _fmpq_vec_init(len);
                    fmpz * n1 = _fmpz_vec_init(len);
                    fmpz * n2 = _fmpz_vec_init(len);
                    for (i = 0; i < len; i++)
                    {
                        fmpz_randbits(fmpq_numref(a + i), state, bits);
                        fmpz_set(fmpq_denref(a + i), d + i);
                        fmpq_canonicalise(a + i);
                    }
                    TIMEIT(fmpq_vec_get_old(n1, r1, a, len), t1);
                    TIMEIT(_fmpq_vec_get_fmpz_vec_fmpz(n2, r2, a, len), t2);
                    if (!_fmpz_vec_equal(n1, n2, len))
                    {
                        flint_printf("FAIL (numerators): %s %wd %wu\n", kinds[ki], len, bits);
                        flint_abort();
                    }
                    _fmpq_vec_clear(a, len);
                    _fmpz_vec_clear(n1, len);
                    _fmpz_vec_clear(n2, len);
                }

                if (!fmpz_equal(r1, r2))
                {
                    flint_printf("FAIL: %s %wd %wu\n", kinds[ki], len, bits);
                    flint_abort();
                }

                flint_printf("%8s %6wd %5wu %12.4g %12.4g %7.2f\n", kinds[ki], len, bits, t1, t2, t1 / t2);
                fflush(stdout);

                _fmpz_vec_clear(d, len);
            }
        }
    }

    fmpz_clear(r1);
    fmpz_clear(r2);
    flint_rand_clear(state);
    return 0;
}
