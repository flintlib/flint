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

#ifndef MPQS_H
#define MPQS_H

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


typedef struct fac_t /* struct for factors of relations */
{
   slong ind;
   slong exp;
} fac_t;

typedef struct la_col_t  /* matrix column */
{
    slong * data;       /* The list of occupied rows in this column */
    slong weight;       /* Number of nonzero entries in this column */
    slong orig;         /* Original relation number */
} la_col_t;

typedef struct mpqs_s
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

   unsigned char sieve_bits; /* sieve threshold */

   /***************************************************************************
                       POLYNOMIAL DATA
    **************************************************************************/

    fmpz_t A;                /* current value of coefficient A */
    fmpz_t B;                /* B values corresponding to current value of A */
    fmpz_t C;                /* value of coefficient 'C' for current 'A' & 'B' */
    
    mp_limb_t * A_ind;       /* indices of factor base primes dividing A0 */
    fmpz_t * A_divp;        /* (A0 / p) for each prime dividing A0 */
    fmpz_t * inv_p2;     /* newton inverse of p^2 for each factor p of A */

    fmpz_t * B_terms;        /* B_terms[i] = A_divp[i] * (B0_terms[i] * q0^(-1))%p,
                               where 'p' is a prime factor of 'A0' */

    mp_limb_t * A_inv;      /* A0^(-1) mod p, for factor base prime p */
    mp_limb_t ** A_inv2B;    /* A_inv2B[j][i] = 2 * B_terms[j] * A^(-1)  mod p */
    mp_limb_t * soln1;       /* first root of poly */
    mp_limb_t * soln2;       /* second root of poly */

    fmpz_t target_A;         /* approximate target value for A coeff of poly */

    slong s; /* number of prime factors of A coeff */
    slong min; /* minimum FB prime that can appear as factor of A */
    slong span; /* size of set of possible prime factors of A */
    slong fact; /* middle of set of possible prime factors of A */
    slong mid; /* start of range for middle factor */
    slong high; /* end of range for middle factor */

   /***************************************************************************
                       RELATION DATA
   ***************************************************************************/

   slong qsort_rels; /* number of relations to accumulate before sorting */
   slong extra_rels; /* number of extra relations beyond num_primes */
   slong max_factors; /* maximum number of factors a relation can have */

   slong * small; /* exponents of small prime factors in relations */
   fac_t * factor; /* factors for a relation */
   fmpz * Y_arr; /* array of Y's corresponding to relations */
   slong * curr_rel; /* current relation in array of relations */
   slong * relation; /* relation array */

   slong buffer_size; /* size of buffer of relations */
   slong num_relations; /* number of relations so far */

   slong num_factors; /* number of factors found in a relation */

   /***************************************************************************
                       LINEAR ALGEBRA DATA
   ***************************************************************************/

   la_col_t * matrix; /* the main matrix over GF(2) in sparse format */
   la_col_t * unmerged; /* unmerged matrix columns */
   la_col_t ** qsort_arr; /* array of columns ready to be sorted */

   slong num_unmerged; /* number of columns unmerged */
   slong columns; /* number of columns in matrix so far */

   /***************************************************************************
                       SQUARE ROOT DATA
   ***************************************************************************/

   slong * prime_count; /* counts of the exponents of primes appearing in the square */

} mpqs_s;

typedef mpqs_s mpqs_t[1];

/*
   Tuning parameters { bits, ks_primes, fb_primes, small_primes, sieve_size}
   for qsieve_factor where:
     * bits is the number of bits of n
     * ks_primes is the max number of primes to try in Knuth-Schroeppel function
     * fb_primes is the number of factor base primes to use (including k and 2)
     * small_primes is the number of small primes to not factor with (including k and 2)
     * sieve_size is the size of the sieve to use
*/

static const mp_limb_t mpqs_tune[][5] =
{
   {40,   50,    60,  5,   2 *   3000}, /* 12 digit, change factor base to 60 from 50 */
   {50,   50,    80,  5,   2 *   3500}, /* 15 digit, change factor base to 100 from 80 */
   {60,   50,   100,  5,   2 *   4000}, /* 18 digit */
   {70,   50,   300,  6,   2 *   6000}, /* 21 digit */
   {80,   50,   400,  6,   2 *   8000}, /* 24 digit */
   {90,   50,   500,  7,   2 *  10000}, 
   {100, 100,   650,  7,   2 *  13000}, 
   {110, 100,   800,  7,   2 *  15000}, /* 31 digits */
   {120, 100,  1000,  7,   2 *  20000}, 
   {130, 100,  3600,  9,   2 *  32000}, /* 41 digits, changed from 800 */
   {140, 100,  1200,  8,   2 *  28000}, 
   {150, 100,  1800,  8,   2 *  32000}, /* 45 digit */
   {160, 150,  2000,  8,   2 *  40000}, 
   {170, 150,  2200,  9,   2 *  64000}, /* 50 digits */
   {180, 150,  2400,  9,   2 *  64000}, 
   {190, 150,  5400, 10,   2 *  64000}, /* factor base size changed from 2700 to x */
   {200, 150,  3600, 10,   2 *  64000}, /* 60 digits */
   {210, 150,  6000, 12,   2 *  64000}, 
   {220, 200,  7500, 15,   2 *  64000}, 
   {230, 200,  8500, 17,   2 *  64000}, /* 70 digits */
   {240, 200, 18000, 19,   2 *  64000}, 
   {250, 200, 24000, 19,   2 *  64000}, /* 75 digits */
   {260, 200, 55000, 25,   2 * 128000}, /* 80 digits */
   {270, 200, 64000, 27,   2 * 128000}
};

/* number of entries in the tuning table */

#define MPQS_TUNE_SIZE (sizeof(mpqs_tune) / (5 * sizeof(mp_limb_t)))

void mpqs_init(mpqs_t mpqs_inf, fmpz_t m);

void mpqs_clear(mpqs_t mpqs_inf);

mp_limb_t mpqs_knuth_schroeppel(mpqs_t mpqs_inf);

mp_limb_t mpqs_primes_init(mpqs_t mpqs_inf);

prime_t * mpqs_compute_factor_base(mp_limb_t * small_factor, mpqs_t mpqs_inf,
                                   slong num_primes);

void mpqs_linalg_init(mpqs_t mpqs_inf);

void mpqs_linalg_re_init(mpqs_t mpqs_inf);

void mpqs_linalg_re_alloc(mpqs_t mpqs_inf);

#ifdef __cplusplus
}
#endif

#endif
