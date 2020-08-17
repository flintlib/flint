/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void
arith_bell_number_multi_mod(fmpz_t res, ulong n)
{
    fmpz_comb_temp_t temp;
    fmpz_comb_t comb;
    nmod_t mod;
    mp_ptr primes, residues;
    slong k, num_primes;
    flint_bitcnt_t size, prime_bits;

    size = arith_bell_number_size(n);
    prime_bits = FLINT_BITS - 1;
    num_primes = (size + prime_bits - 1) / prime_bits;

    primes = flint_malloc(num_primes * sizeof(mp_limb_t));
    residues = flint_malloc(num_primes * sizeof(mp_limb_t));

    primes[0] = n_nextprime(UWORD(1) << prime_bits, 0);
    for (k = 1; k < num_primes; k++)
        primes[k] = n_nextprime(primes[k-1], 0);

    for (k = 0; k < num_primes; k++)
    {
        nmod_init(&mod, primes[k]);
        residues[k] = arith_bell_number_nmod(n, mod);
    }

    fmpz_comb_init(comb, primes, num_primes);
    fmpz_comb_temp_init(temp, comb);

    fmpz_multi_CRT_ui(res, residues, comb, temp, 0);

    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);

    flint_free(primes);
    flint_free(residues);
}
