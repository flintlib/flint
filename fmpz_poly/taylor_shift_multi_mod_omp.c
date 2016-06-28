/*
    Copyright (C) 2014 Fredrik Johansson
    Authored 2016 Daniel S. Roche

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_taylor_shift_multi_mod_omp(fmpz * poly, const fmpz_t c, slong len)
{
    slong xbits, ybits, num_primes, i;
    mp_ptr primes;
    mp_ptr * residues;

    if (len <= 1 || fmpz_is_zero(c))
        return;

    xbits = _fmpz_vec_max_bits(poly, len);

    if (xbits == 0)
        return;

    /* If poly has degree D and coefficients at most |C|, the
       output has coefficient at most D * |C| * 2^D * c^D */
    xbits = FLINT_ABS(xbits) + 1;
    ybits = xbits + len + FLINT_BIT_COUNT(len);

    if (!fmpz_is_pm1(c))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_pow_ui(t, c, len);
        ybits += fmpz_bits(t);
        fmpz_clear(t);
    }

    /* Use primes greater than 2^(FLINT_BITS-1) */
    num_primes = (ybits + (FLINT_BITS - 1) - 1) / (FLINT_BITS - 1);
    primes = flint_malloc(sizeof(mp_limb_t) * num_primes);
    primes[0] = n_nextprime(UWORD(1) << (FLINT_BITS - 1), 1);
    for (i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i-1], 1);

    /* Space for poly reduced modulo the primes */
    residues = flint_malloc(sizeof(mp_ptr) * num_primes);
    for (i = 0; i < num_primes; i++)
        residues[i] = flint_malloc(sizeof(mp_limb_t) * len);

#pragma omp parallel private(i)
    {
        mp_ptr tmp;
        slong j;

        fmpz_comb_t comb;
        fmpz_comb_temp_t comb_temp;

        tmp = flint_malloc(sizeof(mp_limb_t) * num_primes);
        fmpz_comb_init(comb, primes, num_primes);
        fmpz_comb_temp_init(comb_temp, comb);

        /* _fmpz_vec_multi_mod_ui_threaded(residues, poly, len, primes, num_primes, 0); */
#pragma omp for schedule(static)
        for (i = 0; i < len; i++)
        {
            fmpz_multi_mod_ui(tmp, poly + i, comb, comb_temp);
            for (j = 0; j < num_primes; j++)
                residues[j][i] = tmp[j];
        }

        /* _fmpz_poly_multi_taylor_shift_threaded(residues, len, c, primes, num_primes); */
#pragma omp for
        for (i = 0; i < num_primes; i++)
        {
            nmod_t mod;
            mp_limb_t p, cm;

            p = primes[i];
            nmod_init(&mod, p);
            cm = fmpz_fdiv_ui(c, p);
            _nmod_poly_taylor_shift(residues[i], cm, len, mod);
        }

        /* _fmpz_vec_multi_mod_ui_threaded(residues, poly, len, primes, num_primes, 1); */
#pragma omp for schedule(static) nowait
        for (i = 0; i < len; i++)
        {
            for (j = 0; j < num_primes; j++)
                tmp[j] = residues[j][i];
            fmpz_multi_CRT_ui(poly + i, tmp, comb, comb_temp, 1);
        }

        flint_free(tmp);
        fmpz_comb_clear(comb);
        fmpz_comb_temp_clear(comb_temp);

        flint_parallel_cleanup();
    }

    for (i = 0; i < num_primes; i++)
        flint_free(residues[i]);
    flint_free(residues);
    flint_free(primes);
}
