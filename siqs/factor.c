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
#include "qsieve.h"

/*
   Returns a factor of n.
   Assumes n is not prime and not a perfect power.
*/

mp_limb_t qsieve_factor(fmpz_t n, fmpz_factor_t factors)
{
    qs_t qs_inf;
    mp_limb_t small_factor;
    ulong exp = 0;
    fmpz_t x;
    fmpz_init(x);

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
	        fmpz_divexact_ui(x, qs_inf->n, small_factor);
	        fmpz_init_set(qs_inf->n, x);
	        exp++;
        }

        _fmpz_factor_append_ui(factors, small_factor, exp);
	    exp = 0;
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
            fmpz_divexact_ui(x, qs_inf->n, small_factor);
            fmpz_init_set(qs_inf->n, x);
            exp++;
        }

        _fmpz_factor_append_ui(factors, small_factor, exp);
	    exp = 0;
    }

    if (factors->num) return 0;

    while (small_factor == 0)
        small_factor = qsieve_primes_increment(qs_inf, qs_inf->num_primes / 10);

    if (small_factor)
    {
        while (fmpz_fdiv_ui(qs_inf->n, small_factor) == 0)
        {
            fmpz_divexact_ui(x, qs_inf->n, small_factor);
            fmpz_init_set(qs_inf->n, x);
            exp++;
        }

        _fmpz_factor_append_ui(factors, small_factor, exp);
	    exp = 0;
    }

    fmpz_clear(x);
    qsieve_clear(qs_inf);

    return 0;
}
