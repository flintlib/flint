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

#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

/**************************IGNORE*********************************************/

#include "C:\Users\measure\Documents\GitHub\flint2\siqs\qsieve.h"

/*****************************************************************************/

/*
   Returns a factor of n.
   Assumes n is not prime and not a perfect power.
*/

mp_limb_t qsieve_factor(fmpz_t n)
{
    qs_t qs_inf;
    mp_limb_t factor = 0, t;
    slong rels_found = 0;
    char * sieve;
    slong ncols, nrows, i, count;
    uint64_t * nullrows;
    uint64_t mask;
    flint_rand_t state;
    fmpz_t X, Y;

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

    factor = qsieve_knuth_schroeppel(qs_inf);

    if (factor)
	{
	    qsieve_clear(qs_inf);
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

    /* compute factor base primes and associated data*/
    factor = qsieve_ll_primes_init(qs_inf);
	if (factor)
    {
        qsieve_ll_clear(qs_inf);
        return factor;
	}

    return factor;
}
