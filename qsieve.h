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

/* 
    Debug verbosity (bits are as follows):
     >0   - print headings for each phase of the sieve and some 
            simple data about the number of FB primes used, etc. 
      2   - print poly A coeffs
      4   - print polynomials used
      8   - print relations
      16  - print statistics for sieve contents
      32  - print X values for each sieve entry that passes initial sieve test
      64  - print information about the number of duplicate relations
      128 - print lanczos information and errors
*/
#define QS_DEBUG 0

#define CACHE_SIZE 65536 /* size of L1 cache */

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

typedef struct la_col_t /* matrix column */
{
  slong * data;   /* The list of occupied rows in this column */
  slong weight;   /* Number of nonzero entries in this column */
  slong orig;         /* Original relation number */
} la_col_t;

typedef struct qs_s
{
   mp_limb_t hi; /* Number to factor */
   mp_limb_t lo;

   mp_bitcnt_t bits; /* Number of bits of n */
   
   ulong ks_primes; /* number of Knuth-Schroeppel primes */

   mp_limb_t k; /* Multiplier */
   fmpz_t kn; /* kn as a multiprecision integer */

   slong num_primes; /* number of factor base primes including k and 2 */
   slong small_primes; /* number of primes to not sieve with */
   slong sieve_size; /* size of sieve to use */

   prime_t * factor_base; /* data about factor base primes */

   int * sqrts; /* square roots of kn mod factor base primes */

   char sieve_bits; /* number of bits to exceed in sieve */

   /******************
     Polynomial data
   ******************/

   mp_limb_t A; /* coefficient A */
   mp_limb_t B; /* coefficient B */
   fmpz_t C; /* coefficient C */

   mp_limb_t * A_ind; /* indices of factor base primes dividing A */
   mp_limb_t * A_modp; /* (A/p) mod p for each prime p dividing A */
   mp_limb_t * inv_p2; /* newton inverse of p^2 for each factor p of A */

   mp_limb_t * B_terms; /* 
                           Let A_i = (A/p) mod p where p is the i-th prime
                           which is a factor of A, then B_terms[i] is
                           {p^(1/2) / A_i} mod p * (A/p) where we take 
                           the smaller square root of p
                        */
 
   mp_limb_t * A_inv; /* A^(-1) mod p */

   mp_limb_t ** A_inv2B; /* A_inv[j][i] = 2*B_terms[j]*A^(-1) mod p */

   mp_limb_t * soln1; /* first root of poly */
   mp_limb_t * soln2; /* second root of poly */

   mp_limb_t target_A; /* approximate target value for A coeff of poly */

   slong s; /* number of prime factors of A coeff */
   slong min; /* minimum FB prime that can appear as factor of A */
   slong span; /* size of set of possible prime factors of A */
   slong fact; /* middle of set of possible prime factors of A */
   slong mid; /* start of range for middle factor */
   slong high; /* end of range for middle factor */


   /*********************
     Relations data
   **********************/

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

   /*********************
     Linear algebra data
   **********************/

   la_col_t * matrix; /* the main matrix over GF(2) in sparse format */
   la_col_t * unmerged; /* unmerged matrix columns */
   la_col_t ** qsort_arr; /* array of columns ready to be sorted */

   slong num_unmerged; /* number of columns unmerged */
   slong columns; /* number of columns in matrix so far */

   /*********************
     Square root data
   **********************/

   slong * prime_count; /* counts of the exponents of primes appearing in the square */

   /*********************
     Statistics
   **********************/

#if (QS_DEBUG & 16)
   slong * sieve_tally;
#endif

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
static const mp_limb_t qsieve_ll_tune[][5] =
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
#define QS_LL_TUNE_SIZE (sizeof(qsieve_ll_tune)/(5*sizeof(mp_limb_t)))

#define P_GOODNESS 100 /* within what factor of target_A must A be */
#define P_GOODNESS2 200 /* within what factor of target_A must A be when s = 2 */

#define BITS_ADJUST 10 /* no. bits less than f(X) to qualify for trial division */

FLINT_DLL void qsieve_ll_init(qs_t qs_inf, mp_limb_t hi, mp_limb_t lo);

FLINT_DLL void qsieve_ll_clear(qs_t qs_inf);

FLINT_DLL mp_limb_t qsieve_ll_knuth_schroeppel(qs_t qs_inf);

FLINT_DLL mp_limb_t qsieve_ll_primes_init(qs_t qs_inf);

FLINT_DLL mp_limb_t qsieve_ll_poly_init(qs_t qs_inf);

FLINT_DLL void qsieve_ll_linalg_init(qs_t qs_inf);

FLINT_DLL void qsieve_ll_compute_poly_data(qs_t qs_inf);

FLINT_DLL void qsieve_ll_compute_A_factor_offsets(qs_t qs_inf);

FLINT_DLL void qsieve_ll_compute_C(qs_t qs_inf);

FLINT_DLL slong qsieve_ll_collect_relations(qs_t qs_inf, char * sieve);

FLINT_DLL slong qsieve_ll_merge_sort(qs_t qs_inf);
      
FLINT_DLL slong qsieve_ll_merge_relations(qs_t qs_inf);

FLINT_DLL slong qsieve_ll_insert_relation(qs_t qs_inf, fmpz_t Y);

FLINT_DLL mp_limb_t qsieve_ll_factor(mp_limb_t hi, mp_limb_t lo);

static __inline__ void insert_col_entry(la_col_t * col, slong entry)
{
   if (((col->weight >> 4) << 4) == col->weight) /* need more space */
   {
       if (col->weight != 0) col->data = 
           (slong *) flint_realloc(col->data, (col->weight + 16)*sizeof(slong));
       else col->data = (slong *) flint_malloc(16*sizeof(slong));
   }
   
   col->data[col->weight] = entry;
   col->weight++;
}

static __inline__ void copy_col(la_col_t * col2, la_col_t * col1)
{
   col2->weight = col1->weight;
   col2->data = col1->data;
   col2->orig = col1->orig;
}

static __inline__ void swap_cols(la_col_t * col2, la_col_t * col1)
{
   la_col_t temp;
   
   temp.weight = col1->weight;
   temp.data = col1->data;
   temp.orig = col1->orig;
   
   col1->weight = col2->weight;
   col1->data = col2->data;
   col1->orig = col2->orig;
   
   col2->weight = temp.weight;
   col2->data = temp.data;
   col2->orig = temp.orig;
}

static __inline__ void clear_col(la_col_t * col)
{
   col->weight = 0;
}

static __inline__ void free_col(la_col_t * col)
{
   if (col->weight) flint_free(col->data);
}

FLINT_DLL uint64_t get_null_entry(uint64_t * nullrows, slong i, slong l);

FLINT_DLL void reduce_matrix(qs_t qs_inf, slong * nrows, slong * ncols, la_col_t * cols);

uint64_t * block_lanczos(flint_rand_t state, slong nrows, slong dense_rows, 
                                                       slong ncols, la_col_t *B);

FLINT_DLL void qsieve_ll_square_root(fmpz_t X, fmpz_t Y, qs_t qs_inf,
                             uint64_t * nullrows, slong ncols, slong l, fmpz_t N);

#ifdef __cplusplus
}
#endif

#endif
