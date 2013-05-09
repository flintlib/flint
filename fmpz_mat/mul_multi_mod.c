/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


void
_fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B,
    len_t bits)
{
    len_t i, j;

    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;

    len_t num_primes;
    len_t primes_bits;
    mp_limb_t * primes;
    mp_limb_t * residues;

    nmod_mat_t * mod_C;
    nmod_mat_t * mod_A;
    nmod_mat_t * mod_B;

    primes_bits = NMOD_MAT_OPTIMAL_MODULUS_BITS;

    if (bits < primes_bits)
    {
        primes_bits = bits;
        num_primes = 1;
    }
    else
    {
        /* Round up in the division */
        num_primes = (bits + primes_bits - 1) / primes_bits;
    }

    /* Initialize */
    primes = flint_malloc(sizeof(mp_limb_t) * num_primes);
    primes[0] = n_nextprime(1UL << primes_bits, 0);
    for (i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i-1], 0);

    residues = flint_malloc(sizeof(mp_limb_t) * num_primes);

    mod_A = flint_malloc(sizeof(nmod_mat_t) * num_primes);
    mod_B = flint_malloc(sizeof(nmod_mat_t) * num_primes);
    mod_C = flint_malloc(sizeof(nmod_mat_t) * num_primes);
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_init(mod_A[i], A->r, A->c, primes[i]);
        nmod_mat_init(mod_B[i], B->r, B->c, primes[i]);
        nmod_mat_init(mod_C[i], C->r, C->c, primes[i]);
    }

    fmpz_comb_init(comb, primes, num_primes);
    fmpz_comb_temp_init(comb_temp, comb);

    /* Calculate residues of A */
    for (i = 0; i < A->r * A->c; i++)
    {
        fmpz_multi_mod_ui(residues, &A->entries[i], comb, comb_temp);
        for (j = 0; j < num_primes; j++)
            mod_A[j]->entries[i] = residues[j];
    }

    /* Calculate residues of B */
    for (i = 0; i < B->r * B->c; i++)
    {
        fmpz_multi_mod_ui(residues, &B->entries[i], comb, comb_temp);
        for (j = 0; j < num_primes; j++)
            mod_B[j]->entries[i] = residues[j];
    }

    /* Multiply */
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_mul(mod_C[i], mod_A[i], mod_B[i]);
    }

    /* Chinese remaindering */
    for (i = 0; i < C->r * C->c; i++)
    {
        for (j = 0; j < num_primes; j++)
            residues[j] = mod_C[j]->entries[i];
        fmpz_multi_CRT_ui(&C->entries[i], residues, comb, comb_temp, 1);
    }

    /* Cleanup */
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_clear(mod_A[i]);
        nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }

    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);

    fmpz_comb_temp_clear(comb_temp);
    fmpz_comb_clear(comb);

    flint_free(residues);
    flint_free(primes);
}

void
fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    len_t A_bits;
    len_t B_bits;

    A_bits = fmpz_mat_max_bits(A);
    B_bits = fmpz_mat_max_bits(B);

    _fmpz_mat_mul_multi_mod(C, A, B, FLINT_ABS(A_bits) + FLINT_ABS(B_bits) \
        + FLINT_BIT_COUNT(A->c) + 1);
}
