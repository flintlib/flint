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

void sample_divrem(void * arg, ulong count)
{
    ulong modn;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    ulong n = info->n;
    ulong c;

    slong i, j;

    nmod_poly_t poly, div, quo, rem;

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

        nmod_poly_init_mod(quo, poly->mod);
        nmod_poly_init_mod(rem, poly->mod);
        nmod_poly_init_mod(div, poly->mod);

        /* divisor xnmc:  b == x**n - c */
        nmod_poly_set_coeff_ui(div, n, 1);
        nmod_poly_set_coeff_ui(div, 0, n_negmod(c, modn));

        prof_start();
        nmod_poly_divrem(quo, rem, poly, div);
        prof_stop();

        nmod_poly_clear(quo);
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

    nmod_poly_t poly, quo, rem;

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

        nmod_poly_init_mod(quo, poly->mod);
        nmod_poly_init_mod(rem, poly->mod);

        prof_start();
        nmod_poly_divrem_xnmc(quo, rem, poly, n, c);
        prof_stop();

        nmod_poly_clear(quo);
        nmod_poly_clear(rem);
        nmod_poly_clear(poly);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_generic(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < 100; j++)
            _nmod_poly_evaluate_nmod(poly, length, pt, mod);
        prof_stop();
    }

    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

void sample_precomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < 100; j++)
        {
            const ulong pt_precomp = n_mulmod_precomp_shoup(pt, n);
            _nmod_poly_evaluate_nmod_precomp(poly, length, pt, pt_precomp, n);
        }
        prof_stop();
    }

    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

void sample_precomp_lazy(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < 100; j++)
        {
            const ulong pt_precomp = n_mulmod_precomp_shoup(pt, n);
            ulong val = _nmod_poly_evaluate_nmod_precomp_lazy(poly, length, pt, pt_precomp, n);
            if (val >= 2*n)
                val -= 2*n;
            else if (val >= n)
                val -= n;
        }
        prof_stop();
    }

    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    const int nblengths = 17;
    slong lengths[17] = {5, 10, 20, 40, 80,
                        160, 320, 640, 1280, 2560,
                        5120, 10240, 20480, 40960, 81920,
                        163840, 327680};

    double min, max;
    double mins_divrem[nblengths];
    double mins_interface[nblengths];
    double mins_generic[nblengths];
    double mins_precomp[nblengths];
    double mins_precomp_lazy[nblengths];
    info_t info;
    flint_bitcnt_t i;

    flint_printf("profiled: divrem | interface | generic | precomp | precomp_lazy\n");

    for (i = 62; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        for (int ct = 0; ct <= 2; ct++)
        {
            info.ctype = ct;
            printf("nbits = %ld, c type = %d\n", i, ct);
            for (int len = 0; len < nblengths; ++len)
            {
                info.length = lengths[len];
                info.n = FLINT_MAX(1, lengths[len] / 5);

                prof_repeat(&min, &max, sample_divrem, (void *) &info);
                mins_divrem[len-1] = min;
                prof_repeat(&min, &max, sample_interface, (void *) &info);
                mins_interface[len-1] = min;
                /* prof_repeat(&min, &max, sample_generic, (void *) &info); */
                mins_generic[len-1] = min;
                /* prof_repeat(&min, &max, sample_precomp, (void *) &info); */
                mins_precomp[len-1] = min;
                /* prof_repeat(&min, &max, sample_precomp_lazy, (void *) &info); */
                mins_precomp_lazy[len-1] = min;
            }

            if (i < FLINT_BITS-1)
            {
                for (int len = 0; len < nblengths; ++len)
                {
                    flint_printf("   len %ld\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e",
                                 lengths[len],
                                 mins_divrem[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_interface[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_precomp[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_precomp_lazy[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR);
                    flint_printf("\n");
                }
            }
            else if (i < FLINT_BITS)
            {
                for (int len = 0; len < nblengths; ++len)
                {
                    flint_printf("   len %ld\t%.1e\t%.1e\t%.1e\t%.1e\t na",
                                 lengths[len],
                                 mins_divrem[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_interface[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_precomp[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR);
                    flint_printf("\n");
                }
            }
            else  // i == FLINT_BITS
            {
                for (int len = 0; len < nblengths; ++len)
                {
                    flint_printf("   len %ld\t%.1e\t%.1e\t%.1e\t na \t na",
                                 lengths[len],
                                 mins_divrem[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_interface[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR,
                                 mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR);
                    flint_printf("\n");
                }
            }
        }

        flint_printf("\n");
    }

    return 0;
}
