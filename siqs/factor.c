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

/*
   Returns a factor of n.
   Assumes n is not prime and not a perfect power.
*/

mp_limb_t qsieve_factor(fmpz_t n, fmpz_factor_t factors)
{
    qs_t qs_inf;
    mp_limb_t small_factor, factor = 0, t;
    ulong exp = 0;
    char * sieve;
    slong ncols, nrows, i, count;
    uint64_t * nullrows;
    uint64_t mask;
    flint_rand_t state;
    fmpz_t temp, X, Y;
    fmpz_init(temp);
    fmpz_init(X);
    fmpz_init(Y);

    /************************************************************************
        INITIALISATION:

        Initialise the qs_t structure.
    ************************************************************************/

    qsieve_init(qs_inf, n);

    /************************************************************************
        KNUTH SCHROEPPEL:

        Try to compute a multiplier k such that there are a lot of small primes
        which are quadratic residues modulo kn. If a small factor of n is found
        during this process it is returned.
    ************************************************************************/

    small_factor = qsieve_knuth_schroeppel(qs_inf);

    if (small_factor)
    {
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

    /************************************************************************
        COMPUTE FACTOR BASE:

        Try to compute a multiplier k such that there are a lot of small primes
        which are quadratic residues modulo kn. If a small factor of n is found
        during this process it is returned.
    ************************************************************************/

    /* compute factor base primes and associated data*/
    small_factor = qsieve_primes_init(qs_inf);

    if (small_factor)
    {
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

  /************************************************************************
        INITIALISE POLYNOMIAL DATA:

        Create space for all the polynomial information
   ************************************************************************/

    qsieve_poly_init(qs_inf);

    flint_printf("poly initialized\n");

  /************************************************************************
        INITIALISE RELATION/LINAXYXYXYXY DATA:

        Create space for all the relations and matrix information
    ************************************************************************/

    qsieve_linalg_init(qs_inf);

    flint_printf("linear algebra initialized\n");

    /************************************************************************
        SIEVE:

        Sieve for relations
    ************************************************************************/

    sieve = (char *) flint_malloc(qs_inf->sieve_size + sizeof(ulong));

    qs_inf->sieve_bits = 30;

    qsieve_collect_relations(qs_inf, sieve);

    flint_printf("relation collected\n");

    flint_free(sieve);

    /************************************************************************
        REDUCE MATRIX:

        Perform some light filtering on the matrix
    ************************************************************************/

    ncols = qs_inf->num_primes + qs_inf->extra_rels;
    nrows = qs_inf->num_primes;

	reduce_matrix(qs_inf, &nrows, &ncols, qs_inf->matrix);

    /************************************************************************
        BLOCK LANCZOS:

        Find extra_rels nullspace vectors (if they exist)
    ************************************************************************/

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

    /************************************************************************
        SQUARE ROOT:

        Compute the square root and take the GCD of X-Y with N
    ************************************************************************/

    fmpz_fdiv_q_ui(qs_inf->kn, qs_inf->kn, qs_inf->k); /* divide kn by multiplier */

    fmpz_init(X);
    fmpz_init(Y);

    for (count = 0; count < 64; count++)
    {
        if (mask & ((uint64_t)(1) << count))
        {
            qsieve_square_root(X, Y, qs_inf, nullrows, ncols, count, qs_inf->kn);
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
    fmpz_clear(temp);
    qsieve_clear(qs_inf);

    return 1;
}
