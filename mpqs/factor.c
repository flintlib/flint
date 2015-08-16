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

    Built upon existing FLINT siqs
    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpqs.h"

/********************* qsieve factor **************************/

mp_limb_t
mpqs_factor(fmpz_t n, fmpz_factor_t factors)
{
    mpqs_t mpqs_inf;
    mp_limb_t small_factor, factor = 0;
    ulong exp = 0;
    unsigned char * sieve;
    slong ncols, nrows, i, count, relations = 0;
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

    mpqs_init(mpqs_inf, n);

#if QS_DEBUG
    flint_printf("factoring ");
    fmpz_print(mpqs_inf->n);
    flint_printf(" of %wu bits\n", mpqs_inf->bits);
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

    small_factor = mpqs_knuth_schroeppel(mpqs_inf);

    flint_printf("\n knuth is : %wu\n", mpqs_inf->k);

    if (small_factor)
    {

#if QS_DEBUG
        flint_printf("found small factor %wu in Knuth-Schroeppel\n", small_factor);
#endif
        while (fmpz_fdiv_ui(mpqs_inf->n, small_factor) == 0)
        {
            fmpz_divexact_ui(temp, mpqs_inf->n, small_factor);
            fmpz_init_set(mpqs_inf->n, temp);
            exp++;
        }

        _fmpz_factor_append_ui(factors, small_factor, exp);
        exp = 0;

        return 0;
    }

        /* compute kn */
        fmpz_mul_ui(mpqs_inf->kn, mpqs_inf->n, mpqs_inf->k);

        /* refine qs_inf->bits */
        mpqs_inf->bits = fmpz_bits(mpqs_inf->kn);

    /**************************************************************************
        COMPUTE FACTOR BASE:
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nCompute factor-base\n");
#endif

    /* compute factor base primes and associated data*/
    small_factor = mpqs_primes_init(mpqs_inf);

    if (small_factor)
    {

#if QS_DEBUG
        flint_printf("found small factor %wu while generating factor base\n", small_factor);
#endif

        while (fmpz_fdiv_ui(mpqs_inf->n, small_factor) == 0)
        {
            fmpz_divexact_ui(temp, mpqs_inf->n, small_factor);
            fmpz_init_set(mpqs_inf->n, temp);
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

    mpqs_linalg_init(mpqs_inf);

    /**************************************************************************
        POLYNOMIAL INITIALIZATION AND SIEVEING:
        Sieve for relations
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nPolynomial Initializaton and Sieveing\n");
#endif

    sieve = flint_malloc(mpqs_inf->sieve_size + sizeof(ulong));
    mpqs_inf->sieve_bits = 64;

    mpqs_poly_init(mpqs_inf); /* computes first poly */

    while (mpqs_inf->columns <= mpqs_inf->num_primes + mpqs_inf->extra_rels)
    {
        relations += mpqs_collect_relations(mpqs_inf, sieve);

#if QS_DEBUG
        flint_printf("%wd/%wd relations.\n", relations, mpqs_inf->num_primes + mpqs_inf->extra_rels);
#endif

        mpqs_compute_next_poly(mpqs_inf);
    }

    flint_free(sieve);


    /************************************************************************
        REDUCE MATRIX:
        
        Perform some light filtering on the matrix
    ************************************************************************/

    ncols = mpqs_inf->num_primes + mpqs_inf->extra_rels;
    nrows = mpqs_inf->num_primes;

#if QS_DEBUG
    flint_printf("Reduce matrix:\n");
#endif

    mpqs_reduce_matrix(mpqs_inf, &nrows, &ncols, mpqs_inf->matrix); 
 
    /************************************************************************
        BLOCK LANCZOS:
        
        Find extra_rels nullspace vectors (if they exist)
    ************************************************************************/

#if QS_DEBUG
    flint_printf("Block lanczos:\n");
#endif

    flint_randinit(state); /* initialise the random generator */
   
    do /* repeat block lanczos until it succeeds */
    {
        nullrows = mpqs_block_lanczos(state, nrows, 0, ncols, mpqs_inf->matrix);
    } while (nullrows == NULL); 
        
    for (i = 0, mask = 0; i < ncols; i++) /* create mask of nullspace vectors */
        mask |= nullrows[i];

    for (i = count = 0; i < 64; i++) /* count nullspace vectors found */
    {
      if (mask & ((uint64_t)(1) << i))
        count++;
    }

    flint_randclear(state); /* clean up random state */

#if QS_DEBUG
    flint_printf("%wd nullspace vectors found\n", count);
#endif

    /************************************************************************
        SQUARE ROOT:
        
        Compute the square root and take the GCD of X-Y with N
    ************************************************************************/

#if QS_DEBUG
    flint_printf("Square root:\n");
#endif
  
    fmpz_fdiv_q_ui(mpqs_inf->kn, mpqs_inf->kn, mpqs_inf->k); /* divide kn by multiplier */
   
    fmpz_init(X);
    fmpz_init(Y);

    for (count = 0; count < 64; count++)
    {
        if (mask & ((uint64_t)(1) << count))
        {
            mpqs_square_root(X, Y, mpqs_inf, nullrows, ncols, count, mpqs_inf->kn); 
            fmpz_sub(X, X, Y);
            fmpz_gcd(X, X, mpqs_inf->kn);
         
            if (fmpz_cmp(X, mpqs_inf->kn) != 0 && fmpz_cmp_ui(X, 1) != 0) /* have a factor */
            {
                if (fmpz_size(X)!= 1)
                    fmpz_fdiv_q(X, mpqs_inf->kn, X); /* take smaller of two factors */
                factor = fmpz_get_ui(X); 
                break;
            }
        }
    }

#if QS_DEBUG
    flint_printf("\nCleanup\n");
#endif

    mpqs_linalg_clear(mpqs_inf);
    flint_free(sieve);
    fmpz_clear(X);
    fmpz_clear(Y);
    fmpz_clear(temp);
    mpqs_clear(mpqs_inf);
    flint_free(nullrows);

    return 1;
}
