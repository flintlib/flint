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


#define QS_DEBUG 128

typedef struct prime_t
{
   mp_limb_t pinv;     /* precomputed inverse */
   int p;              /* prime */
   char size;
} prime_t;


typedef struct fac_t    /* struct for factors of relations */
{
   slong ind;
   slong exp;
} fac_t;

typedef struct la_col_t  /* matrix column */
{
	slong * data;		/* The list of occupied rows in this column */
	slong weight;		/* Number of nonzero entries in this column */
	slong orig;         /* Original relation number */
} la_col_t;


typedef struct hash_t   /* entry in hash table */
{
    mp_limb_t prime;    /* value of prime */
    mp_limb_t next;     /* next prime which have same hash value as 'prime' */
    mp_limb_t count;    /* number of occurrence of 'prime' */
} hash_t;

typedef struct relation_t  /* format for relation */
{
    mp_limb_t lp;          /* large prime, is 1, if relation is full */
    slong num_factors;     /* number of factors, excluding small factor */
    slong * small;         /* exponent of small factors */
    fac_t * factor;        /* factor of relation */
    fmpz_t Y;              /* square root of sieve value for relation */
} relation_t;

typedef struct qs_s
{
   fmpz_t n;               /* Number to factor */

   mp_bitcnt_t bits;       /* Number of bits of n */

   ulong ks_primes;        /* number of Knuth-Schroeppel primes */

   mp_limb_t k;            /* Multiplier */
   fmpz_t kn;              /* kn as a multiprecision integer */

   slong num_primes;       /* number of factor base primes including k and 2 */

   prime_t * factor_base;  /* data about factor base primes */

   int * sqrts;            /* square roots of kn mod factor base primes */

   slong small_primes;     /* number of primes to not sieve with */
   slong second_prime;
   slong sieve_size;       /* size of sieve to use */

   unsigned char sieve_bits;  /* sieve threshold */
   unsigned char sieve_fill;  /* for biasing sieve values */

   /***************************************************************************
                       POLYNOMIAL DATA
    **************************************************************************/

   fmpz_t A;                /* current value of coefficient A */
   fmpz_t A0;               /* possible candidate for A0 i.e. value of
                               coefficient A excluding the non-factor-base
                               prime  */
   slong q_idx;             /* offset of q0 in factor base */

   fmpz_t B;                /* B values corresponding to current value of A */
   fmpz_t C;                /* value of coefficient 'C' for current 'A' & 'B' */
   mp_limb_t * A_ind;       /* indices of factor base primes dividing A0 */
   fmpz_t * A0_divp;        /* (A0 / p) for each prime dividing A0 */
   fmpz_t * B_terms;        /* B_terms[i] = A_divp[i] * (B0_terms[i] * q0^(-1))%p,
                               where 'p' is a prime factor of 'A0' */

   mp_limb_t * B0_terms;    /* B0_terms[i] = (sqrt(kn) * (A0_divp[i])^(-1)) modulo p,
                               where 'p' is a prime factor of 'A0' */

   mp_limb_t * A0_inv;      /* A0^(-1) mod p, for factor base prime p */
   mp_limb_t ** A_inv2B;    /* A_inv2B[j][i] = 2 * B_terms[j] * A^(-1)  mod p */
   mp_limb_t * soln1;       /* first root of poly */
   mp_limb_t * soln2;       /* second root of poly */

   mp_limb_t * xr1;
   mp_limb_t * xr2;

   fmpz_t target_A;         /* approximate target value for A coeff of poly */

   fmpz_t upp_bound;
   fmpz_t low_bound;

   slong curr_poly;         /* rank of polynomial according to gray-code
                               formula */
   slong s;                 /* number of prime factor of A0 */
   slong low;               /* minimum offset in factor base,
                               for possible factor of 'A0'*/

   slong high;              /* maximum offset in factor base,
                               for possible factor of 'A0' */
   slong span;              /* total number of possible factor for 'A0' */

   /* parameters for calculating next subset of possible factor of 'A0' */

   slong h;
   slong m;
   mp_limb_t * curr_subset;

   /***************************************************************************
                       RELATION DATA
   ***************************************************************************/

   FILE * siqs;          /* pointer to file for storing relations */

   slong full_relation;  /* number of full relation */
   slong num_cycles;     /* number of possible full relation from partials */

   slong vertices;       /* number of different prime in partials */
   slong components;     /* equals to 1 */
   slong edges;          /* total number of partials */

   slong table_size;     /* size of table */
   hash_t * table;       /* store 'prime' occurring in partial */
   mp_limb_t * hash_table;  /* to keep track of location of primes in 'table' */

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

} qs_s;

typedef qs_s qs_t[1];

/*
   Tuning parameters { bits, ks_primes, fb_primes, small_primes, sieve_size}
   for qsieve_factor where:
     * bits is the number of bits of n
     * ks_primes is the max number of primes to try in Knuth-Schroeppel function
     * fb_primes is the number of factor base primes to use (including k and 2)
     * small_primes is the number of small primes to not factor with (including k and 2)
     * sieve_size is the size of the sieve to use
*/
static const mp_limb_t qsieve_tune[][5] =
{
   {40,   50,    60,  5,   2 *   3000}, // 12 digit,
   {50,   50,    80,  5,   2 *   3500}, // 15 digit,
   {60,   50,   100,  5,   2 *   4000}, // 18 digit
   {70,   50,   300,  6,   2 *   6000}, // 21 digit
   {80,   50,   400,  6,   2 *   8000}, // 24 digit
   {90,   50,   500,  7,   2 *  10000}, //
   {100, 100,   650,  7,   2 *  13000}, //
   {110, 100,   800,  7,   2 *  15000}, // 31 digits
   {120, 100,  1000,  7,   2 *  20000}, //
   {130, 100,   800,  9,   2 *  32000}, // 41 digits,
   {140, 100,  1200,  8,   2 *  28000}, //
   {150, 100,  1800,  8,   2 *  32000}, // 45 digit
   {160, 150,  2000,  8,   2 *  40000}, //
   {170, 150,  2200,  9,   2 *  64000}, // 50 digits
   {180, 150,  2400,  9,   2 *  64000}, //
   {190, 150,  4000, 10,   2 *  65536}, //
   {200, 150,  4000, 10,   2 *  65536}, // 60 digits
   {210, 150,  6000, 12,   2 *  64000}, //
   {220, 200,  7500, 15,   2 *  64000}, //
   {230, 200,  8500, 17,   2 *  64000}, // 70 digits
   {240, 200, 18000, 19,   2 *  64000}, //
   {250, 200, 24000, 19,   2 *  64000}, // 75 digits
   {260, 200, 55000, 25,   2 * 128000}, // 80 digits
   {270, 200, 64000, 27,   2 * 128000}
};

/* number of entries in the tuning table */
#define QS_TUNE_SIZE (sizeof(qsieve_tune)/(5*sizeof(mp_limb_t)))

#define BITS_ADJUST 10

#define BLOCK_SIZE 65536

void qsieve_init(qs_t qs_inf, fmpz_t n);

mp_limb_t qsieve_knuth_schroeppel(qs_t qs_inf);

void qsieve_clear(qs_t qs_inf);

prime_t *
compute_factor_base(mp_limb_t * small_factor, qs_t qs_inf, slong num_primes);

mp_limb_t qsieve_primes_init(qs_t qs_inf);

mp_limb_t qsieve_primes_increment(qs_t qs_inf, mp_limb_t delta);

mp_limb_t qsieve_poly_init(qs_t qs_inf);

mp_limb_t qsieve_next_A0(qs_t qs_inf);

void qsieve_re_init_A0(qs_t qs_inf);

void qsieve_init_A0(qs_t qs_inf);

void qsieve_compute_pre_data(qs_t qs_inf);

void qsieve_init_poly_first(qs_t qs_inf);

void qsieve_init_poly_next(qs_t qs_inf);

void qsieve_compute_C(qs_t qs_inf);

void qsieve_poly_clear(qs_t qs_inf);

void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve);

void qsieve_do_sieving2(qs_t qs_inf, unsigned char * sieve);

slong qsieve_evaluate_candidate(qs_t qs_inf, slong i, unsigned char * sieve);

slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve);

slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve);

void qsieve_linalg_init(qs_t qs_inf);

void qsieve_linalg_re_init(qs_t qs_inf);

void qsieve_linalg_re_alloc(qs_t qs_inf);

void qsieve_linalg_clear(qs_t qs_inf);

int qsieve_relations_cmp(const void * a, const void * b);

int qsieve_relations_cmp2(const void * a, const void * b);

slong qsieve_merge_sort(qs_t qs_inf);

slong qsieve_merge_relations(qs_t qs_inf);

slong qsieve_insert_relation(qs_t qs_inf, fmpz_t Y);

void qsieve_write_to_file(qs_t qs_inf, mp_limb_t prime, fmpz_t Y);

hash_t * qsieve_get_table_entry(qs_t qs_inf, mp_limb_t prime);

void qsieve_add_to_hashtable(qs_t qs_inf, mp_limb_t prime);

relation_t qsieve_parse_relation(qs_t qs_inf, char * str);

relation_t qsieve_merge_relation(qs_t qs_inf, relation_t  a, relation_t  b);

int qsieve_compare_relation(const void * a, const void * b);

int qsieve_remove_duplicates(relation_t * rel_list, slong num_relations);

void qsieve_insert_relation2(qs_t qs_inf, relation_t * rel_list, slong num_relations);

void qsieve_process_relation(qs_t qs_inf);


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

uint64_t get_null_entry(uint64_t * nullrows, slong i, slong l);

void reduce_matrix(qs_t qs_inf, slong *nrows, slong *ncols, la_col_t *cols);

uint64_t * block_lanczos(flint_rand_t state, slong nrows,
			slong dense_rows, slong ncols, la_col_t *B);

void qsieve_square_root(fmpz_t X, fmpz_t Y, qs_t qs_inf,
   uint64_t * nullrows, slong ncols, slong l, fmpz_t N);


#ifdef __cplusplus
}
#endif

#endif
