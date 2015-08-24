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
    Copyright (C) 2015 Nitin Kumar
******************************************************************************/

#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

#include <inttypes.h>
#define _STDC_FORMAT_MACROS
#include <time.h>

/*
   Returns a factor of n.
   Assumes n is not prime and not a perfect power.
*/

mp_limb_t qsieve_factor(fmpz_t n, fmpz_factor_t factors)
{
    clock_t start = clock(), diff;
    qs_t qs_inf;
    mp_limb_t small_factor, factor = 0, t, delta;
    ulong exp = 0;
    unsigned char * sieve;
    slong ncols, nrows, i, j, k = 0, count, relation = 0, num_primes;
    uint64_t * nullrows;
    uint64_t mask;
    flint_rand_t state;
    fmpz_t temp, X, Y;
    fmpz_init(temp);
    fmpz_init(X);
    fmpz_init(Y);

    /**************************************************************************
        INITIALISATION:
        Initialise the qs_t structure.
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nstart\n");
#endif

    qsieve_init(qs_inf, n);

#if QS_DEBUG
    flint_printf("factoring ");
    fmpz_print(qs_inf->n);
    flint_printf(" of %wu bits\n", qs_inf->bits);
#endif

    /**************************************************************************
        KNUTH SCHROEPPEL:
        Try to compute a multiplier k such that there are a lot of small primes
        which are quadratic residues modulo kn. If a small factor of n is found
        during this process it is returned.
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nKnuth_Schroeppel\n");
#endif

    small_factor = qsieve_knuth_schroeppel(qs_inf);

    if (small_factor)
    {

#if QS_DEBUG
        flint_printf("found small factor %wu in Knuth-Schroeppel\n", small_factor);
#endif

        while (fmpz_fdiv_ui(qs_inf->n, small_factor) == 0)
        {
	        fmpz_divexact_ui(temp, qs_inf->n, small_factor);
	        fmpz_init_set(qs_inf->n, temp);
	        exp++;
        }

        _fmpz_factor_append_ui(factors, small_factor, exp);
        exp = 0;

        return 0;
    }

    /* compute kn */
    fmpz_mul_ui(qs_inf->kn, qs_inf->n, qs_inf->k);

    /* refine qs_inf->bits */
    qs_inf->bits = fmpz_bits(qs_inf->kn);

    /**************************************************************************
        COMPUTE FACTOR BASE:
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nCompute factor-base\n");
#endif

    /* compute factor base primes and associated data*/
    small_factor = qsieve_primes_init(qs_inf);

    if (small_factor)
    {

#if QS_DEBUG
        flint_printf("found small factor %wu while generating factor base\n", small_factor);
#endif

        while (fmpz_fdiv_ui(qs_inf->n, small_factor) == 0)
        {
            fmpz_divexact_ui(temp, qs_inf->n, small_factor);
            fmpz_init_set(qs_inf->n, temp);
            exp++;
        }

        _fmpz_factor_append_ui(factors, small_factor, exp);
        exp = 0;

        return 0;
    }

    /**************************************************************************
        INITIALISE RELATION/LINAXYXYXYXY DATA:
        Create space for all the relations and matrix information
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nInitializing Relations and Linear Algbera\n");
#endif

    qsieve_linalg_init(qs_inf);

    /**************************************************************************
        POLYNOMIAL INITIALIZATION AND SIEVEING:
        Sieve for relations
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nPolynomial Initializaton and Sieveing\n");
#endif

    sieve = flint_malloc(qs_inf->sieve_size + sizeof(ulong));

    qs_inf->sieve_bits = 64;
    qs_inf->sieve_fill = 8;
    qs_inf->q_idx = qs_inf->num_primes;
    qs_inf->siqs = fopen("siqs.dat", "w");

    for (j = qs_inf->small_primes; j < qs_inf->num_primes; j++)
    {
        if (qs_inf->factor_base[j].p > BLOCK_SIZE)
            break;
    }

    qs_inf->second_prime = j;

    while(1)
    {
        if (qs_inf->s)
            qsieve_re_init_A0(qs_inf);
        else
            qsieve_init_A0(qs_inf);

        do
        {
            qsieve_compute_pre_data(qs_inf);

            for (j = qs_inf->num_primes; j < qs_inf->num_primes + qs_inf->ks_primes; j++)
            {

                qs_inf->q_idx  = j;
                relation += qsieve_collect_relations(qs_inf, sieve);

                qs_inf->num_cycles = qs_inf->edges + qs_inf->components - qs_inf->vertices;

                if (qs_inf->full_relation  + qs_inf->num_cycles    >= (qs_inf->num_primes + qs_inf->ks_primes + 1.5 * qs_inf->extra_rels))
                {
                    fclose(qs_inf->siqs);

                    qsieve_process_relation(qs_inf);

    /**************************************************************************
        REDUCE MATRIX:
        Perform some light filtering on the matrix
    **************************************************************************/

                    num_primes = qs_inf->num_primes;
                    qs_inf->num_primes += qs_inf->ks_primes;

                    ncols = qs_inf->num_primes + qs_inf->extra_rels;
                    nrows = qs_inf->num_primes;

                    reduce_matrix(qs_inf, &nrows, &ncols, qs_inf->matrix);

   /**************************************************************************
        BLOCK LANCZOS:
        Find extra_rels nullspace vectors (if they exist)
    **************************************************************************/
#if QS_DEBUG
                    flint_printf("\nBlock Lanczos\n");
#endif
                    flint_randinit(state); /* initialise the random generator */

                    do /* repeat block lanczos until it succeeds */
                    {
                        nullrows = block_lanczos(state, nrows, 0, ncols, qs_inf->matrix);
                    } while (nullrows == NULL);

                    for (i = 0, mask = 0; i < ncols; i++) /* create mask of nullspace vectors */
                        mask |= nullrows[i];

                    for (i = count = 0; i < 64; i++) /* count nullspace vectors found */
                    {
                        if (mask & ((uint64_t)(1) << i))
                            count++;
                    }

                    flint_randclear(state); /* clean up random state */

    /**************************************************************************
        SQUARE ROOT:
        Compute the square root and take the GCD of X-Y with N
    **************************************************************************/

#if QS_DEBUG
                    flint_printf("\nSquare Root\n");
#endif

                    fmpz_init(X);
                    fmpz_init(Y);

                    for (count = 0; count < 64; count++)
                    {
                        if (mask & ((uint64_t)(1) << count))
                        {

                            qsieve_square_root(X, Y, qs_inf, nullrows, ncols, count, qs_inf->kn);

                            fmpz_sub(X, X, Y);
                            fmpz_gcd(X, X, qs_inf->n);

                            if (fmpz_cmp(X, qs_inf->n) != 0 && fmpz_cmp_ui(X, 1) != 0) /* have a factor */
                            {
                                if (fmpz_size(X)!= 1)
                                    fmpz_fdiv_q(X, qs_inf->n, X); /* take smaller of two factors */

                                diff = clock() - start;
                                flint_printf("factor found = ");
                                fmpz_print(X);
                                flint_printf(" time taken = %f \n", diff / 1000.0);
                                goto cleanup;
                            }
                        }
                    }

                    qs_inf->siqs = fopen("siqs.dat", "w");
                    qsieve_linalg_re_init(qs_inf);
                    qs_inf->num_primes = num_primes;
                    relation = 0;

                }
            }

        } while (qsieve_next_A0(qs_inf));

        delta = qs_inf->num_primes / 10;
        num_primes = qs_inf->num_primes + delta;
        qs_inf->num_primes += qs_inf->ks_primes;

#if QS_DEBUG
        flint_printf("\nfactor base increment\n");
#endif
        small_factor = qsieve_primes_increment(qs_inf, delta);
        qs_inf->num_primes = num_primes;

        if (small_factor)
        {

#if QS_DEBUG
            flint_printf("found small factor %wu while incrementing factor base\n, small_factor");
#endif
            goto cleanup;
        }

        qsieve_linalg_re_alloc(qs_inf);
        relation = 0;
    }

#if QS_DEBUG
    flint_printf("\nCleanup\n");
#endif

cleanup:
    flint_free(sieve);
    fmpz_clear(X);
    fmpz_clear(Y);
    flint_free(nullrows);
    fmpz_clear(temp);
    qsieve_clear(qs_inf);

    return 1;
}
