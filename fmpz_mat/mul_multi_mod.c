/*
    Copyright (C) 2010, 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
_fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B,
    mp_bitcnt_t bits)
{
    slong i, j, k;
    slong num_primes;
    mp_bitcnt_t primes_bits;
    mp_ptr primes;

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
    primes[0] = n_nextprime(UWORD(1) << primes_bits, 0);
    for (i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i-1], 0);

    mod_A = flint_malloc(sizeof(nmod_mat_t) * num_primes);
    mod_B = flint_malloc(sizeof(nmod_mat_t) * num_primes);
    mod_C = flint_malloc(sizeof(nmod_mat_t) * num_primes);
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_init(mod_A[i], A->r, A->c, primes[i]);
        nmod_mat_init(mod_B[i], B->r, B->c, primes[i]);
        nmod_mat_init(mod_C[i], C->r, C->c, primes[i]);
    }

    /* Basecase reduction & CRT */
    if (num_primes < 500)
    {
        /* Calculate residues of A */
        for (i = 0; i < A->r; i++)
        {
            for (j = 0; j < A->c; j++)
            {
                for (k = 0; k < num_primes; k++)
                    nmod_mat_entry(mod_A[k], i, j) =
                        fmpz_fdiv_ui(fmpz_mat_entry(A, i, j), primes[k]);
            }
        }

        /* Calculate residues of B */
        for (i = 0; i < B->r; i++)
        {
            for (j = 0; j < B->c; j++)
            {
                for (k = 0; k < num_primes; k++)
                    nmod_mat_entry(mod_B[k], i, j) =
                        fmpz_fdiv_ui(fmpz_mat_entry(B, i, j), primes[k]);
            }
        }

        /* Multiply */
        for (i = 0; i < num_primes; i++)
        {
            nmod_mat_mul(mod_C[i], mod_A[i], mod_B[i]);
        }

        /* Chinese remaindering */
        if (num_primes == 1)
        {
            mp_limb_t r, t, p = primes[0];

            for (i = 0; i < C->r; i++)
            {
                for (j = 0; j < C->c; j++)
                {
                    r = nmod_mat_entry(mod_C[0], i, j);
                    t = p - r;
                    if (t < r)
                        fmpz_neg_ui(fmpz_mat_entry(C, i, j), t);
                    else
                        fmpz_set_ui(fmpz_mat_entry(C, i, j), r);
                }
            }
        }
        else if (num_primes == 2)
        {
            mp_limb_t c1, c2, r1, r2, m1, m2, c1m2[2], c2m1[2], t[3], u[3], M[2];
            m1 = primes[0];
            m2 = primes[1];
            c1 = n_invmod(m2 % m1, m1);
            c2 = n_invmod(m1, m2);  /* Assumes m1 < m2 */
            umul_ppmm(M[1], M[0], m1, m2);
            umul_ppmm(c1m2[1], c1m2[0], c1, m2);
            umul_ppmm(c2m1[1], c2m1[0], c2, m1);

            for (i = 0; i < C->r; i++)
            {
                for (j = 0; j < C->c; j++)
                {
                    r1 = nmod_mat_entry(mod_C[0], i, j);
                    r2 = nmod_mat_entry(mod_C[1], i, j);

                    /* Assumes no overflow (fine with 60-bit moduli) */
                    t[2] = mpn_mul_1(t, c1m2, 2, r1);
                    t[2] += mpn_addmul_1(t, c2m1, 2, r2);

                    /* Assumes M[1] != 0 (fine with 60-bit moduli) */
                    /* todo: write a preinv 3by2 division function */
                    mpn_tdiv_qr(u, t, 0, t, 3, M, 2);

                    sub_ddmmss(u[1], u[0], M[1], M[0], t[1], t[0]);
                    if (u[1] < t[1] || (u[1] == t[1] && u[0] < t[0]))
                        fmpz_neg_uiui(fmpz_mat_entry(C, i, j), u[1], u[0]);
                    else
                        fmpz_set_uiui(fmpz_mat_entry(C, i, j), t[1], t[0]);
                }
            }
        }
        else
        {
            mp_ptr M, Ns, T, U;
            mp_size_t Msize, Nsize;
            mp_limb_t cy, ri;

            M = flint_malloc(sizeof(mp_limb_t) * (num_primes + 1));

            M[0] = primes[0];
            Msize = 1;
            for (i = 1; i < num_primes; i++)
            {
                M[Msize] = cy = mpn_mul_1(M, M, Msize, primes[i]);
                Msize += (cy != 0);
            }

            /* We add terms with Msize + 1 limbs, with one extra limb for the
               carry accumulation. todo: reduce Nsize by 1 when the carries
               do not require an extra limb. */
            Nsize = Msize + 2;

            Ns = flint_malloc(sizeof(mp_limb_t) * Nsize * num_primes);
            T = flint_malloc(sizeof(mp_limb_t) * Nsize);
            U = flint_malloc(sizeof(mp_limb_t) * Nsize);

            for (i = 0; i < num_primes; i++)
            {
                Ns[i * Nsize + (Nsize - 1)] = 0;
                Ns[i * Nsize + (Nsize - 2)] = 0;
                mpn_divrem_1(Ns + i * Nsize, 0, M, Msize, primes[i]);
                ri = mpn_mod_1(Ns + i * Nsize, Msize, primes[i]);
                ri = n_invmod(ri, primes[i]);
                Ns[i * Nsize + Msize] = mpn_mul_1(Ns + i * Nsize, Ns + i * Nsize, Msize, ri);
            }

            for (i = 0; i < C->r; i++)
            {
                for (j = 0; j < C->c; j++)
                {
                    ri = nmod_mat_entry(mod_C[0], i, j);
                    T[Nsize - 1] = mpn_mul_1(T, Ns, Nsize - 1, ri);

                    for (k = 1; k < num_primes; k++)
                    {
                        ri = nmod_mat_entry(mod_C[k], i, j);
                        T[Nsize - 1] += mpn_addmul_1(T, Ns + k * Nsize, Nsize - 1, ri);
                    }

                    mpn_tdiv_qr(U, T, 0, T, Nsize, M, Msize);
                    mpn_sub_n(U, M, T, Msize);

                    if (mpn_cmp(U, T, Msize) < 0)
                    {
                        fmpz_set_ui_array(fmpz_mat_entry(C, i, j), U, Msize);
                        fmpz_neg(fmpz_mat_entry(C, i, j), fmpz_mat_entry(C, i, j));
                    }
                    else
                    {
                        fmpz_set_ui_array(fmpz_mat_entry(C, i, j), T, Msize);
                    }
                }
            }

            flint_free(M);
            flint_free(Ns);
            flint_free(T);
            flint_free(U);
        }
    }
    else   /* Use comb */
    {
        fmpz_comb_t comb;
        fmpz_comb_temp_t comb_temp;
        mp_ptr residues;

        fmpz_comb_init(comb, primes, num_primes);
        fmpz_comb_temp_init(comb_temp, comb);
        residues = flint_malloc(sizeof(mp_limb_t) * num_primes);

        /* Calculate residues of A */
        for (i = 0; i < A->r; i++)
        {
            for (j = 0; j < A->c; j++)
            {
                fmpz_multi_mod_ui(residues, &(A->rows[i][j]), comb, comb_temp);
                for (k = 0; k < num_primes; k++)
                    mod_A[k]->rows[i][j] = residues[k];
            }
        }

        /* Calculate residues of B */
        for (i = 0; i < B->r; i++)
        {
            for (j = 0; j < B->c; j++)
            {
                fmpz_multi_mod_ui(residues, &B->rows[i][j], comb, comb_temp);
                for (k = 0; k < num_primes; k++)
                    mod_B[k]->rows[i][j] = residues[k];
            }
        }

        /* Multiply */
        for (i = 0; i < num_primes; i++)
        {
            nmod_mat_mul(mod_C[i], mod_A[i], mod_B[i]);
        }

        /* Chinese remaindering */
        for (i = 0; i < C->r; i++)
        {
            for (j = 0; j < C->c; j++)
            {
                for (k = 0; k < num_primes; k++)
                    residues[k] = mod_C[k]->rows[i][j];
                fmpz_multi_CRT_ui(&C->rows[i][j], residues, comb, comb_temp, 1);
            }
        }

        fmpz_comb_temp_clear(comb_temp);
        fmpz_comb_clear(comb);
        flint_free(residues);
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
    flint_free(primes);
}

void
fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong A_bits;
    slong B_bits;

    A_bits = fmpz_mat_max_bits(A);
    B_bits = fmpz_mat_max_bits(B);

    _fmpz_mat_mul_multi_mod(C, A, B, FLINT_ABS(A_bits) + FLINT_ABS(B_bits) \
        + FLINT_BIT_COUNT(A->c) + 1);
}
