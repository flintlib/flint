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

    Copyright (C) 2006, 2011 William Hart

******************************************************************************/

#undef ulong /* avoid clash with stdlib */
#include <stdio.h>
#define ulong unsigned long 

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

/* 
   Factor n = (hi, lo). Returns a factor of n. 
   Assumes n is not prime and not a perfect power.
   Returns 0 if n is too large to be factored by this function.
*/
mp_limb_t qsieve_ll_factor(mp_limb_t hi, mp_limb_t lo)
{
    qs_t qs_inf;
    mp_limb_t factor = 0, t;
    len_t rels_found = 0;
    char * sieve;
    len_t ncols, nrows, i, count;
    uint64_t * nullrows;
    uint64_t mask;
    flint_rand_t state;
    fmpz_t X, Y;

    /************************************************************************
        INITIALISATION:
          
        Initialise the qs_t structure. 
    ************************************************************************/
#if QS_DEBUG
    printf("\nStart:\n");
#endif

    qsieve_ll_init(qs_inf, hi, lo);

#if QS_DEBUG /* print some diagnostic information */
    {
       mpz_t _n;
       mpz_init(_n);
       mpz_set_ui(_n, hi);
       mpz_mul_2exp(_n, _n, FLINT_BITS);
       mpz_add_ui(_n, _n, lo);
       gmp_printf("Factoring %Zd of %ld bits\n", _n, qs_inf->bits);
       mpz_clear(_n);
    }
#endif

    /************************************************************************
        KNUTH SCHROEPPEL:
        
        Try to compute a multiplier k such that there are a lot of small primes
        which are quadratic residues modulo kn. If a small factor of n is found
        during this process it is returned.
    ************************************************************************/
#if QS_DEBUG
    printf("\nKnuth-Schroeppel:\n");
#endif

    factor = qsieve_ll_knuth_schroeppel(qs_inf); 
    if (factor) 
	{
#if QS_DEBUG
        printf("Found small factor %ld in Knuth-Schroeppel\n", factor);
#endif
		qsieve_ll_clear(qs_inf);
        return factor;
	}

    /* compute kn */
    fmpz_set_ui(qs_inf->kn, hi);
    fmpz_mul_2exp(qs_inf->kn, qs_inf->kn, FLINT_BITS);
    fmpz_add_ui(qs_inf->kn, qs_inf->kn, lo);
    fmpz_mul_ui(qs_inf->kn, qs_inf->kn, qs_inf->k);

    /* refine qs_inf->bits */
    qs_inf->bits = fmpz_bits(qs_inf->kn);
    if (qs_inf->bits > 2*FLINT_BITS)
    {
		qsieve_ll_clear(qs_inf);
        return 0; /* kn is too large to factor */
    }

    /************************************************************************
        COMPUTE FACTOR BASE:
        
        Try to compute a multiplier k such that there are a lot of small primes
        which are quadratic residues modulo kn. If a small factor of n is found
        during this process it is returned.
    ************************************************************************/
#if QS_DEBUG
    printf("\nCompute factor base:\n");
#endif

    /* compute factor base primes and associated data*/
    factor = qsieve_ll_primes_init(qs_inf);
	if (factor) 
    {
#if QS_DEBUG
        printf("Found small factor %ld whilst generating factor base\n", factor);
#endif
		qsieve_ll_clear(qs_inf);
        return factor;
	}
    
    /************************************************************************
        INITIALISE POLYNOMIAL DATA:
        
        Create space for all the polynomial information
    ************************************************************************/
#if QS_DEBUG
    printf("\nInitialise poly:\n");
#endif

    /* set (hi, lo) to kn */
    umul_ppmm(t, lo, lo, qs_inf->k);
    hi = hi*qs_inf->k + t;
    qs_inf->hi = hi;
    qs_inf->lo = lo;

    qsieve_ll_poly_init(qs_inf);

    /************************************************************************
        INITIALISE RELATION/LINALG DATA:
        
        Create space for all the relations and matrix information
    ************************************************************************/
#if QS_DEBUG
    printf("\nInitialise relations and linear algebra:\n");
#endif

    qsieve_ll_linalg_init(qs_inf);

    /************************************************************************
        SIEVE:
        
        Sieve for relations
    ************************************************************************/
#if QS_DEBUG
    printf("\nSieve:\n");
#endif

    sieve = flint_malloc(qs_inf->sieve_size + sizeof(ulong));

    while (rels_found < qs_inf->num_primes + qs_inf->extra_rels)
    {
        rels_found += qsieve_ll_collect_relations(qs_inf, sieve);

#if (QS_DEBUG & 128)
        printf("%ld/%ld relations.\n", rels_found, qs_inf->num_primes + qs_inf->extra_rels);
#endif
    }

    flint_free(sieve);

    /************************************************************************
        REDUCE MATRIX:
        
        Perform some light filtering on the matrix
    ************************************************************************/

    ncols = qs_inf->num_primes + qs_inf->extra_rels;
    nrows = qs_inf->num_primes;

#if QS_DEBUG
    printf("Reduce matrix:\n");
#endif

	reduce_matrix(qs_inf, &nrows, &ncols, qs_inf->matrix); 
 
    /************************************************************************
        BLOCK LANCZOS:
        
        Find extra_rels nullspace vectors (if they exist)
    ************************************************************************/

#if QS_DEBUG
    printf("Block lanczos:\n");
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

#if QS_DEBUG
    printf("%ld nullspace vectors found\n", count);
#endif

    /************************************************************************
        SQUARE ROOT:
        
        Compute the square root and take the GCD of X-Y with N
    ************************************************************************/

#if QS_DEBUG
    printf("Square root:\n");
#endif
	
    fmpz_fdiv_q_ui(qs_inf->kn, qs_inf->kn, qs_inf->k); /* divide kn by multiplier */
   
    fmpz_init(X);
    fmpz_init(Y);

    for (count = 0; count < 64; count++)
    {
        if (mask & ((uint64_t)(1) << count))
        {
            qsieve_ll_square_root(X, Y, qs_inf, nullrows, ncols, count, qs_inf->kn); 
            fmpz_sub(X, X, Y);
            fmpz_gcd(X, X, qs_inf->kn);
         
            if (fmpz_cmp(X, qs_inf->kn) != 0 && fmpz_cmp_ui(X, 1) != 0) /* have a factor */
            {
                if (fmpz_size(X)!= 1)
                    fmpz_fdiv_q(X, qs_inf->kn, X); /* take smaller of two factors */
                factor = fmpz_get_ui(X); 
                break;
            }
        }
    }

    fmpz_clear(X);
    fmpz_clear(Y);
    flint_free(nullrows);

    /************************************************************************
        CLEAN UP:
        
        Free all used memory
    ************************************************************************/

#if QS_DEBUG
    printf("\nClean up:\n");
#endif

    qsieve_ll_clear(qs_inf);

#if QS_DEBUG
    printf("\nDone.\n");
#endif

    return factor;
}
