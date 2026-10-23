/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpn_extras.h"

/*
    Usage: p-multi_crt [pbits] [chunk_bits] [mod_base_bits] [preinv_cutoff]

    Times flint_mpn_multi_mod and flint_mpn_multi_crt for primes of pbits
    bits and increasing total sizes, optionally with explicit tuning
    parameters (0 = default).
*/

int main(int argc, char ** argv)
{
    slong bits_list[] = {20, 32, 50, 62, 64};
    slong totbits_list[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072};
    slong nb, ntot;
    slong only_bits = argc > 1 ? atol(argv[1]) : 0;
    slong chunk_bits = argc > 2 ? atol(argv[2]) : 0;
    slong mod_base_bits = argc > 3 ? atol(argv[3]) : 0;
    slong preinv = argc > 4 ? atol(argv[4]) : 0;
    flint_rand_t state;

    flint_rand_init(state);

    flint_printf("%6s %6s %8s | %10s %10s\n", "pbits", "n", "totbits", "mod (ns)", "crt (ns)");

    for (nb = 0; nb < 5; nb++)
    {
        slong bits = bits_list[nb];

        if (only_bits && bits != only_bits)
            continue;

        for (ntot = 0; ntot < 12; ntot++)
        {
            slong n = FLINT_MAX(1, totbits_list[ntot] / bits), i;
            ulong * primes = flint_malloc(n * sizeof(ulong));
            ulong * res = flint_malloc(n * sizeof(ulong));
            ulong * out, * x, * tmp;
            ulong p;
            flint_mpn_crt_t C;
            fmpz_t X;
            slong xn;
            double t_mod, t_crt;

            p = (bits == 64) ? n_nextprime(UWORD(1) << 63, 1) : n_nextprime(UWORD(1) << (bits - 1), 1);
            for (i = 0; i < n; i++)
            {
                primes[i] = p;
                p = n_nextprime(p, 1);
            }

            if (chunk_bits || mod_base_bits || preinv)
                flint_mpn_crt_init_tuned(C, primes, n, FLINT_MPN_CRT_MOD | FLINT_MPN_CRT_CRT,
                    chunk_bits ? chunk_bits : 2048,
                    mod_base_bits ? mod_base_bits : 2048,
                    preinv ? preinv : 16);
            else
                flint_mpn_crt_init(C, primes, n);

            fmpz_init(X);
            fmpz_set_ui_array(X, C->prod, C->prod_len);
            fmpz_randm(X, state, X);
            xn = fmpz_size(X);
            x = flint_malloc((xn + 1) * sizeof(ulong));
            fmpz_get_ui_array(x, xn + 1, X);
            out = flint_malloc(C->prod_len * sizeof(ulong));
            tmp = flint_malloc(C->tmp_limbs * sizeof(ulong));

            TIMEIT_START;
            flint_mpn_multi_mod(res, x, xn, C, tmp);
            TIMEIT_STOP_VALUES(t_mod, t_mod);

            TIMEIT_START;
            flint_mpn_multi_crt(out, res, C, 1, tmp);
            TIMEIT_STOP_VALUES(t_crt, t_crt);

            flint_printf("%6wd %6wd %8wd | %10.0f %10.0f\n", bits, n, n * bits, t_mod * 1e9, t_crt * 1e9);
            fflush(stdout);

            flint_free(primes);
            flint_free(res);
            flint_free(out);
            flint_free(x);
            flint_free(tmp);
            fmpz_clear(X);
            flint_mpn_crt_clear(C);
        }
    }

    flint_rand_clear(state);
    return 0;
}
