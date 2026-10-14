/*
   Copyright (C) 2026 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
   Geometric evaluation, interpolation and extrapolation with a precomputation
   G of length len, built without (function = 7) and with (function = 15)
   stored transforms: times of one call of each, and of the two inits.

   usage: p-geometric_progression_fft [mod ...]

   where each mod is a bit length (a random prime of that size is used),
   "fft" (a prime of the default fft_small context, for which fft_small
   needs a single prime), or "ntt" (998244353 = 119*2^23 + 1, which
   fft_small uses directly, without chinese remaindering).
   Shapes: evaluation of a polynomial of length len at len points,
   interpolation from len points, forward extrapolation from len/2 values to
   the next len/2 points.
*/

#include <stdlib.h>
#include <string.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/* minimum over 3 runs of at least 30 ms of the time per call, in us */
#define TIME_MIN(t, stmt)                                                 \
    do {                                                                  \
        timeit_t __T;                                                     \
        slong __reps = 1;                                                 \
        (t) = 1e300;                                                      \
        for (int __r = 0; __r < 3; )                                      \
        {                                                                 \
            timeit_start(__T);                                            \
            for (slong __k = 0; __k < __reps; __k++)                      \
            {                                                             \
                stmt;                                                     \
            }                                                             \
            timeit_stop(__T);                                             \
            if (__T->wall < 30)                                           \
                __reps *= 2;                                              \
            else                                                          \
            {                                                             \
                (t) = FLINT_MIN((t), __T->wall * 1000.0 / __reps);        \
                __r++;                                                    \
            }                                                             \
        }                                                                 \
    } while (0)

static const slong lens[] = {16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512,
    700, 1024, 1400, 2048, 3000, 4096, 6000, 8192, 12000, 16384};

int main(int argc, char * argv[])
{
    const char * defaults[] = {"20", "30", "fft", "ntt", "60", "64"};
    const char ** mods = (argc > 1) ? (const char **) argv + 1 : defaults;
    slong nmods = (argc > 1) ? argc - 1 : 6;
    flint_rand_t state;

    flint_rand_init(state);

    for (slong m = 0; m < nmods; m++)
    {
        ulong p;
        nmod_t mod;

        if (!strcmp(mods[m], "fft"))
            p = UWORD(1108307720798209);  /* 63*2^44 + 1 */
        else if (!strcmp(mods[m], "ntt"))
            p = UWORD(998244353);
        else
            p = n_randprime(state, atoi(mods[m]), 1);

        nmod_init(&mod, p);
        flint_printf("\nmodulus %s (%wu), times in us, ratio = without / with transforms\n", mods[m], p);
        flint_printf("%6s | %9s %9s %5s | %9s %9s %5s | %9s %9s %5s | %9s %9s\n", "len",
                     "eval", "eval+tf", "ratio", "interp", "interp+tf", "ratio",
                     "extrap", "extrap+tf", "ratio", "init", "init+tf");

        for (ulong l = 0; l < sizeof(lens) / sizeof(lens[0]); l++)
        {
            slong len = lens[l], h = len / 2;
            ulong r = n_primitive_root_prime(p);
            nmod_geometric_progression_t G0, G1;
            nn_ptr a = _nmod_vec_init(len), b = _nmod_vec_init(len), c = _nmod_vec_init(len);
            double te[2], ti[2], tx[2], tinit[2];

            if (p - 1 < (ulong) 2 * len)
                continue;

            _nmod_vec_rand(a, state, len, mod);

            TIME_MIN(tinit[0], nmod_geometric_progression_t T;
                     _nmod_geometric_progression_init_function(T, r, len, mod, 7);
                     nmod_geometric_progression_clear(T));
            TIME_MIN(tinit[1], nmod_geometric_progression_t T;
                     _nmod_geometric_progression_init_function(T, r, len, mod, 15);
                     nmod_geometric_progression_clear(T));

            _nmod_geometric_progression_init_function(G0, r, len, mod, 7);
            _nmod_geometric_progression_init_function(G1, r, len, mod, 15);

            for (int k = 0; k < 2; k++)
            {
                nmod_geometric_progression_struct * G = (k == 0) ? G0 : G1;
                TIME_MIN(te[k], _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(b, a, len, G, len, mod));
                TIME_MIN(ti[k], _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(b, a, G, len, mod));
                TIME_MIN(tx[k], nmod_poly_extrapolate_geometric_precomp(b, h, a, h, h, G));
            }

            /* sanity check */
            nmod_poly_extrapolate_geometric_precomp(b, h, a, h, h, G0);
            nmod_poly_extrapolate_geometric_precomp(c, h, a, h, h, G1);
            if (!_nmod_vec_equal(b, c, h))
                flint_abort();

            flint_printf("%6wd | %9.2f %9.2f %5.2f | %9.2f %9.2f %5.2f | %9.2f %9.2f %5.2f | %9.2f %9.2f\n",
                         len, te[0], te[1], te[0] / te[1], ti[0], ti[1], ti[0] / ti[1],
                         tx[0], tx[1], tx[0] / tx[1], tinit[0], tinit[1]);
            fflush(stdout);

            nmod_geometric_progression_clear(G0);
            nmod_geometric_progression_clear(G1);
            _nmod_vec_clear(a);
            _nmod_vec_clear(b);
            _nmod_vec_clear(c);
        }
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
