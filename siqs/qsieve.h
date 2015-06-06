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

#ifndef QSIEVE_H
#define QSIEVE_H

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

#ifdef __cplusplus
 extern "C" {
#endif

#if FLINT_BITS==64
   #ifndef uint64_t
   #define uint64_t ulong
   #endif
#else
   #include <stdint.h>
#endif

typedef struct prime_t
{
   mp_limb_t pinv; /* precomputed inverse */
   int p; /* prime */
   char size;
} prime_t;

typedef struct qs_s
{
   fmpz_t n; /* Number to factor */

   mp_bitcnt_t bits; /* Number of bits of n */

   ulong ks_primes; /* number of Knuth-Schroeppel primes */

   mp_limb_t k; /* Multiplier */
   fmpz_t kn; /* kn as a multiprecision integer */

   slong num_primes; /* number of factor base primes including k and 2 */

   prime_t * factor_base; /* data about factor base primes */

   int * sqrts; /* square roots of kn mod factor base primes */

   slong small_primes; /* number of primes to not sieve with */
   slong sieve_size; /* size of sieve to use */

   /***************************************************************************
                       POLYNOMIAL DATA
    **************************************************************************/

   mp_limb_t A;
   mp_limb_t A0;  /* value of coeff A excluding the non-factor-base prime */
   mp_limb_t q0;   /* primes immediately following largest factor-base prime */
   mp_limb_t B;
   fmpz_t C;

   mp_limb_t * A_ind;   /* indices of factor base primes dividing A */
   mp_limb_t * A_modp;  /* (A / p) mod p for each prime dividing A */
   mp_limb_t * A0_divp;  /* (A0 / p) for each prime dividing A0 */
   mp_limb_t * B_terms;

   mp_limb_t * B0_terms;  /* similar to B_terms for A0 */

   mp_limb_t * A_inv;      /* A^(-1) mod p */
   mp_limb_t ** A_inv2B;   /* A_inv2B[j][i] = 2 * B_terms[j] * A^(-1)  mod p */

   mp_limb_t * soln1;     /* first root of poly */
   mp_limb_t * soln2;     /* second root of poly */

   mp_limb_t target_A;  /* approximate target value for A coeff of poly */

   slong s;


} qs_s;

typedef qs_s qs_t[1];

/*
   Tuning parameters { bits, ks_primes, fb_primes, small_primes }
   for qsieve_ll_factor where:
     * bits is the number of bits of n
     * ks_primes is the max number of primes to try in Knuth-Schroeppel algo
     * fb_primes is the number of factor base primes to use (including k and 2)
     * small_primes is the number of small primes to not factor with (including k and 2)
     * sieve_size is the size of the sieve to use
*/
static const mp_limb_t qsieve_tune[][5] =
{
    {0, 50, 80, 2, 14000 },
    {30, 50, 80, 2, 16000 },
    {40, 50, 100, 3, 18000 },
    {50, 50, 120, 3, 20000 },
    {60, 50, 140, 4, 22000 },
    {70, 50, 160, 4, 24000 },
    {80, 100, 180, 5, 26000 },
    {90, 100, 200, 5, 28000 },
    {100, 100, 250, 6, 30000 },
    {110, 100, 300, 6, 34000 },
    {120, 100, 500, 7, 40000 },
    {130, 100, 550, 7, 60000 }
};

/* number of entries in the tuning table */
#define QS_TUNE_SIZE (sizeof(qsieve_tune)/(5*sizeof(mp_limb_t)))

#define P_GOODNESS 5

#define P_GOODNESS2 200

/* need to replace with function prototype */

/******************************************************************************/

#include "test_helpers.c"
#include "clear.c"
#include "init.c"
#include "knuth_schroeppel.c"
#include "primes_init.c"
#include "factor.c"
#include "compute_poly_data.c"

/******************************************************************************/


#ifdef __cplusplus
}
#endif

#endif
