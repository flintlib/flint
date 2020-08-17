/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "arith.h"

#define CRT_MAX_RESOLUTION 16

void
arith_bell_number_vec_multi_mod(fmpz * res, slong n)
{
    fmpz_comb_t comb[CRT_MAX_RESOLUTION];
    fmpz_comb_temp_t temp[CRT_MAX_RESOLUTION];
    mp_ptr primes, residues;
    mp_ptr * polys;
    nmod_t mod;
    slong i, j, k, num_primes, num_primes_k, resolution;
    flint_bitcnt_t size, prime_bits;

    if (n < 1)
        return;

    resolution = FLINT_MAX(1, FLINT_MIN(CRT_MAX_RESOLUTION, n / 16));

    size = arith_bell_number_size(n);
    prime_bits = FLINT_BITS - 1;
    num_primes = (size + prime_bits - 1) / prime_bits;

    primes = flint_malloc(num_primes * sizeof(mp_limb_t));
    residues = flint_malloc(num_primes * sizeof(mp_limb_t));
    polys = flint_malloc(num_primes * sizeof(mp_ptr));

    /* Compute Bell numbers mod p */
    primes[0] = n_nextprime(UWORD(1)<<prime_bits, 0);
    for (k = 1; k < num_primes; k++)
        primes[k] = n_nextprime(primes[k-1], 0);

    for (k = 0; k < num_primes; k++)
    {
        /* flint_printf("prime %wd of %wd\n", k, num_primes); */
        polys[k] = _nmod_vec_init(n);
        nmod_init(&mod, primes[k]);
        arith_bell_number_nmod_vec(polys[k], n, mod);
    }

    /* Init CRT comb */
    for (i = 0; i < resolution; i++)
    {
        fmpz_comb_init(comb[i], primes, num_primes * (i + 1) / resolution);
        fmpz_comb_temp_init(temp[i], comb[i]);
    }

    /* Reconstruction */
    for (k = 0; k < n; k++)
    {
        size = arith_bell_number_size(k);
        /* Use only as large a comb as needed */
        num_primes_k = (size + prime_bits - 1) / prime_bits;
        for (i = 0; i < resolution; i++)
        {
            if (comb[i]->num_primes >= num_primes_k)
                break;
        }
        num_primes_k = comb[i]->num_primes;
        for (j = 0; j < num_primes_k; j++)
            residues[j] = polys[j][k];
        fmpz_multi_CRT_ui(res + k, residues, comb[i], temp[i], 0);
    }

    /* Cleanup */
    for (k = 0; k < num_primes; k++)
        _nmod_vec_clear(polys[k]);

    for (i = 0; i < resolution; i++)
    {
        fmpz_comb_temp_clear(temp[i]);
        fmpz_comb_clear(comb[i]);
    }

    flint_free(primes);
    flint_free(residues);
    flint_free(polys);
}
