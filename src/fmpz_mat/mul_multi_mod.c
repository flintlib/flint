/*
    Copyright (C) 2010, 2018, 2026 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

void _fmpz_mat_mul_multi_mod(
    fmpz_mat_t C,
    const fmpz_mat_t A,
    const fmpz_mat_t B,
    int sign,
    flint_bitcnt_t bits)
{
    slong i;
    slong m, k, n;
    flint_bitcnt_t primes_bits;
    fmpz_comb_t comb;
    fmpz_comb_struct * comb_ptr;
    ulong first_prime;
    int squaring = (A == B);
    nmod_mat_t * mod_A, * mod_B, * mod_C;
    nn_ptr primes;
    slong num_primes;

    m = A->r;
    k = A->c;
    n = B->c;

    if (m < 1 || n < 1 || k < 1)
    {
        fmpz_mat_zero(C);
        return;
    }

    FLINT_ASSERT(sign == 0 || sign == 1);
    bits += sign;

    primes_bits = NMOD_MAT_OPTIMAL_MODULUS_BITS;

    if (bits < primes_bits || bits <= FLINT_BITS - 1)
    {
        num_primes = 1;
        first_prime = UWORD(1) << bits;
    }
    else
    {
        num_primes = 1 + (bits - (FLINT_BITS - 1) + primes_bits - 1) / primes_bits;
        first_prime = UWORD(1) << (FLINT_BITS - 1);
    }

    primes = FLINT_ARRAY_ALLOC(num_primes, ulong);
    primes[0] = first_prime;
    if (num_primes > 1)
    {
        primes[1] = n_nextprime(UWORD(1) << primes_bits, 0);
        for (i = 2; i < num_primes; i++)
            primes[i] = n_nextprime(primes[i - 1], 0);
    }

    mod_A = FLINT_ARRAY_ALLOC(num_primes, nmod_mat_t);
    mod_B = squaring ? NULL : FLINT_ARRAY_ALLOC(num_primes, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(num_primes, nmod_mat_t);

    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_init(mod_A[i], A->r, A->c, primes[i]);
        if (!squaring)
            nmod_mat_init(mod_B[i], B->r, B->c, primes[i]);
        nmod_mat_init(mod_C[i], C->r, C->c, primes[i]);
    }

    if (num_primes > FMPZ_MAT_MOD_PRIMES_COMB_CUTOFF ||
        num_primes > FMPZ_MAT_CRT_PRIMES_COMB_CUTOFF)
    {
        fmpz_comb_init(comb, primes, num_primes);
        comb_ptr = comb;
    }
    else
    {
        comb_ptr = NULL;
    }

    fmpz_mat_multi_mod_2_ui_precomp(mod_A, squaring ? NULL : mod_B, num_primes,
            A, squaring ? NULL : B, num_primes > FMPZ_MAT_MOD_PRIMES_COMB_CUTOFF ? comb_ptr : NULL);

    for (i = 0; i < num_primes; i++)
        nmod_mat_mul(mod_C[i], mod_A[i], squaring ? mod_A[i] : mod_B[i]);

    fmpz_mat_multi_CRT_ui_precomp(C, mod_C, num_primes,
            num_primes > FMPZ_MAT_CRT_PRIMES_COMB_CUTOFF ? comb_ptr : NULL, sign);

    /* Cleanup */
    if (comb_ptr != NULL)
        fmpz_comb_clear(comb);

    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_clear(mod_A[i]);
        if (!squaring)
            nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }

    flint_free(mod_A);
    if (!squaring)
        flint_free(mod_B);
    flint_free(mod_C);
    flint_free(primes);
}


void
fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong Abits, Bbits;
    int sign = 0;
    flint_bitcnt_t Cbits;

    Abits = fmpz_mat_max_bits(A);
    Bbits = fmpz_mat_max_bits(B);

    if (Abits < 0)
    {
        sign = 1;
        Abits = -Abits;
    }

    if (Bbits < 0)
    {
        sign = 1;
        Bbits = -Bbits;
    }

    Cbits = Abits + Bbits + FLINT_BIT_COUNT(A->c);

    _fmpz_mat_mul_multi_mod(C, A, B, sign, Cbits);
}
