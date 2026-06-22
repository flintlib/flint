/*
   Copyright (C) 2025 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod.h"
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_vec.h"

static inline ulong
_get_c_from_ctype(int ctype, flint_rand_t state, ulong modn)
{
    if (ctype == 0)
        return modn - 1;
    else if (ctype == 1)
        return 1;
    else
        return n_randint(state, modn);
}

typedef struct
{
    flint_bitcnt_t bits;
    slong length;
    ulong n;
    int ctype;  /* 0: c == -1; 1: c == 1; else: general*/
} info_t;

void sample_rem(void * arg, ulong count)
{
    ulong modn;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    ulong n = info->n;
    ulong c;

    slong i, j;

    nmod_poly_t poly, div, rem;

    FLINT_TEST_INIT(state);

    for (i = 0; i < (slong)count; i++)
    {
        modn = n_randprime(state, bits, 1);

        nmod_poly_init2(poly, modn, length);
        _nmod_poly_set_length(poly, length);
        for (j = 0; j < length; j++)
            poly->coeffs[j] = n_randint(state, modn);
        _nmod_poly_normalise(poly);

        c = _get_c_from_ctype(info->ctype, state, modn);

        nmod_poly_init_mod(rem, poly->mod);
        nmod_poly_init_mod(div, poly->mod);

        /* divisor xnmc:  b == x**n - c */
        nmod_poly_set_coeff_ui(div, n, 1);
        nmod_poly_set_coeff_ui(div, 0, n_negmod(c, modn));

        prof_start();
        nmod_poly_rem(rem, poly, div);
        prof_stop();

        nmod_poly_clear(rem);
        nmod_poly_clear(poly);
        nmod_poly_clear(div);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_interface(void * arg, ulong count)
{
    ulong modn;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    ulong n = info->n;
    ulong c;

    slong i, j;

    nmod_poly_t poly, rem;

    FLINT_TEST_INIT(state);

    for (i = 0; i < (slong)count; i++)
    {
        modn = n_randprime(state, bits, 1);

        nmod_poly_init2(poly, modn, length);
        _nmod_poly_set_length(poly, length);
        for (j = 0; j < length; j++)
            poly->coeffs[j] = n_randint(state, modn);
        _nmod_poly_normalise(poly);

        c = _get_c_from_ctype(info->ctype, state, modn);

        nmod_poly_init_mod(rem, poly->mod);

        prof_start();
        nmod_poly_rem_xnmc(rem, poly, n, c);
        prof_stop();

        nmod_poly_clear(rem);
        nmod_poly_clear(poly);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_generic(void * arg, ulong count)
{
    ulong modn;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    ulong n = info->n;
    ulong c;

    slong i;

    nn_ptr r, a;

    FLINT_TEST_INIT(state);

    for (i = 0; i < (slong)count; i++)
    {
        modn = n_randprime(state, bits, 1);
        nmod_init(&mod, modn);

        r = _nmod_vec_init(n);
        a = _nmod_vec_init(length);
        _nmod_vec_rand(a, state, length, mod);

        c = _get_c_from_ctype(info->ctype, state, modn);

        prof_start();
        _nmod_poly_rem_xnmc(r, a, length, n, c, mod);
        prof_stop();

        _nmod_vec_clear(r);
        _nmod_vec_clear(a);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_precomp(void * arg, ulong count)
{
    ulong modn;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    ulong n = info->n;
    ulong c;

    slong i;

    nn_ptr r, a;

    FLINT_TEST_INIT(state);

    for (i = 0; i < (slong)count; i++)
    {
        modn = n_randprime(state, bits, 1);
        nmod_init(&mod, modn);

        r = _nmod_vec_init(n);
        a = _nmod_vec_init(length);
        _nmod_vec_rand(a, state, length, mod);

        c = _get_c_from_ctype(info->ctype, state, modn);

        prof_start();
        const ulong c_precomp = n_mulmod_precomp_shoup(c, modn);
        _nmod_poly_rem_xnmc_precomp(r, a, length, n, c, c_precomp, modn);
        prof_stop();

        _nmod_vec_clear(a);
        _nmod_vec_clear(r);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_precomp_lazy(void * arg, ulong count)
{
    ulong modn;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    ulong n = info->n;
    ulong c;

    slong i;

    nn_ptr r, a;

    FLINT_TEST_INIT(state);

    for (i = 0; i < (slong)count; i++)
    {
        modn = n_randprime(state, bits, 1);
        nmod_init(&mod, modn);

        r = _nmod_vec_init(n);
        a = _nmod_vec_init(length);
        _nmod_vec_rand(a, state, length, mod);

        c = _get_c_from_ctype(info->ctype, state, modn);

        prof_start();
        const ulong c_precomp = n_mulmod_precomp_shoup(c, modn);
        _nmod_poly_rem_xnmc_precomp_lazy(r, a, length, n, c, c_precomp, modn);
        prof_stop();

        _nmod_vec_clear(r);
        _nmod_vec_clear(a);
    }

    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    const int nblengths = 15;
    slong lengths[15] = {5, 40, 80, 160, 320,
                        640, 1280, 2560, 5120, 10240,
                        20480, 40960, 81920, 163840, 327680};

    double min, max;
    double mins_rem;
    double mins_interface;
    double mins_generic;
    double mins_precomp;
    double mins_precomp_lazy;
    info_t info;
    flint_bitcnt_t i;

    flint_printf("profiled:\n");
    flint_printf("\t1. rem\n\t2. poly interface\n\t3. vec general\n");
    flint_printf("\t4. vec precomp\n\t5. precomp_lazy (no excess correction)\n");
    flint_printf("column c: 0 -> c == -1, 1 -> c == 1, 2 -> c general\n");
    flint_printf("note: variant 2. expected to be faster than variants 3.4.5. when c == 1 or c == -1\n\n");
    flint_printf("bit c len    n      1.      2.      3.      4.      5.\n");

    for (i = 62; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        for (int ct = 0; ct <= 2; ct++)
        {
            info.ctype = ct;
            for (int len = 0; len < nblengths; ++len)
            {
                info.length = lengths[len];

                info.n = 1;
                while (info.n < (ulong)lengths[len])
                {
                    prof_repeat(&min, &max, sample_rem, (void *) &info);
                    mins_rem = min/(double)FLINT_CLOCK_SCALE_FACTOR;
                    prof_repeat(&min, &max, sample_interface, (void *) &info);
                    mins_interface = min/(double)FLINT_CLOCK_SCALE_FACTOR;
                    prof_repeat(&min, &max, sample_generic, (void *) &info);
                    mins_generic = min/(double)FLINT_CLOCK_SCALE_FACTOR;
                    prof_repeat(&min, &max, sample_precomp, (void *) &info);
                    mins_precomp = min/(double)FLINT_CLOCK_SCALE_FACTOR;
                    prof_repeat(&min, &max, sample_precomp_lazy, (void *) &info);
                    mins_precomp_lazy = min/(double)FLINT_CLOCK_SCALE_FACTOR;

                    if (i < FLINT_BITS-1)
                    {
                        flint_printf("%-4d%-2d%-7d%-7d%.1e %.1e %.1e %.1e %.1e",
                                    info.bits, info.ctype, info.length, info.n,
                                    mins_rem, mins_interface, mins_generic, mins_precomp, mins_precomp_lazy);
                        flint_printf("\n");
                    }
                    else if (i < FLINT_BITS)
                    {
                        flint_printf("%-4d%-2d%-7d%-7d%.1e %.1e %.1e %.1e  na",
                                    info.bits, info.ctype, info.length, info.n,
                                    mins_rem, mins_interface, mins_generic, mins_precomp);
                        flint_printf("\n");
                    }
                    else  /* i == FLINT_BITS */
                    {
                        flint_printf("%-4d%-2d%-7d%-7d%.1e %.1e %.1e  na     na",
                                    info.bits, info.ctype, info.length, info.n,
                                    mins_rem, mins_interface, mins_generic);
                        flint_printf("\n");
                    }

                    info.n = 5 * info.n;
                }
            }
        }
    }

    return 0;
}
