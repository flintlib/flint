/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef QSIEVE_H
#define QSIEVE_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"

#if HAVE_OPENMP
#include <omp.h> /* must include flint.h first */
#endif

#ifdef __cplusplus
 extern "C" {
#endif

#define QS_DEBUG 0

#define BITS_ADJUST 25 /* added to sieve entries to compensate for approximations */

#define BLOCK_SIZE 65536 /* size of sieving cache block */

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
   slong small_primes;   /* number of small factors */
   slong * small;         /* exponent of small factors */
   fac_t * factor;        /* factor of relation */
   fmpz_t Y;              /* square root of sieve value for relation */
} relation_t;

typedef struct qs_poly_s
{
   fmpz_t B;          /* current B coeff of poly */
   int * soln1;       /* first start position in sieve per prime */
   int * soln2;       /* second start position in sieve per prime */
   int * posn1;       /* temp space for sieving */
   int * posn2;       /* temp space for sieving */
   slong * small;     /* exponents of small prime factors in relations */
   fac_t * factor;    /* factors for a relation */
   slong num_factors; /* number of factors found in a relation */
} qs_poly_s;

typedef qs_poly_s qs_poly_t[1];

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
   slong second_prime;     /* index of first prime bigger than block size */
   slong sieve_size;       /* size of sieve to use */

   unsigned char sieve_bits;  /* sieve threshold */
   unsigned char sieve_fill;  /* for biasing sieve values */

   /***************************************************************************
                       POLYNOMIAL DATA
    **************************************************************************/

   fmpz_t A;                /* current value of coefficient A */
   fmpz_t A0;               /* coefficient A excluding the non-factor-base
                               prime  */
   slong q_idx;             /* offset of q0 in factor base */

   fmpz_t B;                /* B values corresponding to current value of A */
   mp_limb_t * A_ind;       /* indices of factor base primes dividing A0 */
   fmpz_t * A0_divp;        /* (A0 / p) for each prime dividing A0 */
   fmpz_t * B_terms;        /* B_terms[i] = A_divp[i] * (B0_terms[i] * q0^(-1)) % p,
                               where 'p' is a prime factor of 'A0' */

   mp_limb_t * B0_terms;    /* B0_terms[i] = (sqrt(kn) * (A0_divp[i])^(-1)) modulo p,
                               where 'p' is a prime factor of 'A0' */

   mp_limb_t * A0_inv;      /* A0^(-1) mod p, for factor base primes p */
   mp_limb_t ** A_inv2B;    /* A_inv2B[j][i] = 2 * B_terms[j] * A^(-1)  mod p */
   int * soln1;       /* first root of poly */
   int * soln2;       /* second root of poly */

   fmpz_t target_A;         /* approximate target value for A coeff of poly */

   fmpz_t upp_bound;
   fmpz_t low_bound;

   slong s;                 /* number of prime factors of A0 */
   slong low;               /* minimum offset in factor base,
                               for possible factors of 'A0' */

   slong high;              /* maximum offset in factor base,
                               for possible factors of 'A0' */
   slong span;              /* total number of possible factors for 'A0' */

   /* parameters for calculating next subset of possible factor of 'A0' */

   slong h;
   slong m;
   mp_limb_t * curr_subset;

#if QS_DEBUG
   slong poly_count;         /* keep track of the number of polynomials used */
#endif

   qs_poly_s * poly;         /* poly data per thread */

   /***************************************************************************
                       RELATION DATA
   ***************************************************************************/

   FILE * siqs;          /* pointer to file for storing relations */

   slong full_relation;  /* number of full relations */
   slong num_cycles;     /* number of possible full relations from partials */

   slong vertices;       /* number of different primes in partials */
   slong components;     /* equal to 1 */
   slong edges;          /* total number of partials */

   slong table_size;     /* size of table */
   hash_t * table;       /* store 'prime' occurring in partial */
   mp_limb_t * hash_table;  /* to keep track of location of primes in 'table' */

   slong extra_rels;     /* number of extra relations beyond num_primes */
   slong max_factors;    /* maximum number of factors a relation can have */

   fmpz * Y_arr;         /* array of Y's corresponding to relations */
   slong * curr_rel;     /* current relation in array of relations */
   slong * relation;     /* relation array */

   slong buffer_size;    /* size of buffer of relations */
   slong num_relations;  /* number of relations so far */

   ulong small_factor;   /* small factor found when merging relations */

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
     * sieve_bits - sieve_fill
*/

#if HAVE_OPENMP

static const mp_limb_t qsieve_tune[][6] =
{
   {10,   50,   100,  5,   2 *  2000,  30}, /* */
   {20,   50,   120,  6,   2 *  2500,  30}, /* */
   {30,   50,   150,  6,   2 *  2000,  31}, /* */
   {40,   50,   150,  8,   2 *  3000,  32}, /* 12 digits */
   {50,   50,   150,  8,   2 *  3000,  34}, /* 15 digits */
   {60,   50,   150,  9,   2 *  3500,  36}, /* 18 digits */
   {70,  100,   200,  9,   2 *  4000,  42}, /* 21 digits */
   {80,  100,   200,  9,   2 *  6000,  44}, /* 24 digits */
   {90,  100,   200,  9,   2 *  6000,  50}, /* */
   {100, 100,   300,  9,   2 *  7000,  54}, /* */
   {110, 100,   500,  9,   2 *  25000, 62}, /* 31 digits */
   {120, 100,   800,  9,   2 *  30000, 64}, /* */
   {130, 100,  1000,  9,   2 *  30000, 64}, /* 41 digits */
   {140, 100,  1200,  9,   2 *  30000, 66}, /* */
   {150, 100,  1500, 10,   2 *  32000, 68}, /* 45 digit */
   {160, 150,  1800, 11,   2 *  32000, 70}, /* */
   {170, 150,  2000, 12,   2 *  32000, 72}, /* 50 digits */
   {180, 150,  2500, 12,   2 *  32000, 73}, /* */
   {190, 150,  2800, 12,   2 *  32000, 76}, /* */
   {200, 200,  4000, 12,   2 *  32000, 80}, /* 60 digits */
   {210, 100,  3600, 12,   2 *  32000, 83}, /* */
   {220, 300,  6000, 15,   2 *  65536, 87}, /* */
   {230, 350,  8500, 17,   3 *  65536, 90}, /* 70 digits */
   {240, 400, 10000, 19,   4 *  65536, 93}, /* */
   {250, 500, 15000, 19,   4 *  65536, 97}, /* 75 digits */
   {260, 600, 25000, 25,   4 *  65536, 100}, /* 80 digits */
   {270, 800, 35000, 27,   5 *  65536, 104}  /* */
};

#else /* currently tuned for four threads */

static const mp_limb_t qsieve_tune[][6] =
{
   {10,   50,   100,  5,   2 *  2000,  30}, /* */
   {20,   50,   120,  6,   2 *  2500,  30}, /* */
   {30,   50,   150,  6,   2 *  2000,  31}, /* */
   {40,   50,   150,  8,   2 *  3000,  32}, /* 12 digits */
   {50,   50,   150,  8,   2 *  4000,  34}, /* 15 digits */
   {60,   50,   150,  9,   2 *  5000,  36}, /* 18 digits */
   {70,  100,   200,  9,   2 *  6000,  42}, /* 21 digits */
   {80,  100,   200,  9,   2 *  8000,  44}, /* 24 digits */
   {90,  100,   200,  9,   2 *  9000,  50}, /* */
   {100, 100,   300,  9,   2 *  10000,  54}, /* */
   {110, 100,   500,  9,   2 *  30000, 62}, /* 31 digits */
   {120, 100,   800,  9,   2 *  40000, 64}, /* */
   {130, 100,  1000,  9,   2 *  50000, 64}, /* 41 digits */
   {140, 100,  1200,  9,   2 *  65536, 66}, /* */
   {150, 100,  1500, 10,   2 *  65536, 68}, /* 45 digit */
   {160, 150,  1800, 11,   3 *  65536, 70}, /* */
   {170, 150,  2000, 12,   4 *  65536, 72}, /* 50 digits */
   {180, 150,  2500, 12,   5 *  65536, 73}, /* */
   {190, 150,  2800, 12,   6 *  65536, 76}, /* */
   {200, 200,  4000, 12,   6 *  65536, 80}, /* 60 digits */
   {210, 100,  3600, 12,   7 *  65536, 83}, /* */
   {220, 300,  6000, 15,   9 *  65536, 87}, /* */
   {230, 350,  8500, 17,   10 *  65536, 90}, /* 70 digits */
   {240, 400, 10000, 19,   12 *  65536, 93}, /* */
   {250, 500, 15000, 19,   14 *  65536, 97}, /* 75 digits */
   {260, 600, 25000, 25,   15 *  65536, 100}, /* 80 digits */
   {270, 800, 35000, 27,   16 *  65536, 104}  /* */
};

#endif

/* number of entries in the tuning table */
#define QS_TUNE_SIZE (sizeof(qsieve_tune)/(6*sizeof(mp_limb_t)))

void qsieve_init(qs_t qs_inf, const fmpz_t n);

mp_limb_t qsieve_knuth_schroeppel(qs_t qs_inf);

void qsieve_clear(qs_t qs_inf);

void qsieve_factor(fmpz_factor_t factors, const fmpz_t n);

prime_t * compute_factor_base(mp_limb_t * small_factor, qs_t qs_inf,
                                                             slong num_primes);

mp_limb_t qsieve_primes_init(qs_t qs_inf);

mp_limb_t qsieve_primes_increment(qs_t qs_inf, mp_limb_t delta);

mp_limb_t qsieve_poly_init(qs_t qs_inf);

mp_limb_t qsieve_next_A0(qs_t qs_inf);

void qsieve_re_init_A0(qs_t qs_inf);

int qsieve_init_A0(qs_t qs_inf);

void qsieve_compute_pre_data(qs_t qs_inf);

void qsieve_init_poly_first(qs_t qs_inf);

void qsieve_init_poly_next(qs_t qs_inf, slong i);

void qsieve_compute_C(fmpz_t C, qs_t qs_inf, qs_poly_t poly);

void qsieve_poly_copy(qs_poly_t poly, qs_t qs_inf);

void qsieve_poly_clear(qs_t qs_inf);

void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly);

void qsieve_do_sieving2(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly);

slong qsieve_evaluate_candidate(qs_t qs_inf, ulong i, unsigned char * sieve, qs_poly_t poly);

slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly);

slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve);

void qsieve_linalg_init(qs_t qs_inf);

void qsieve_linalg_re_init(qs_t qs_inf);

void qsieve_linalg_re_alloc(qs_t qs_inf);

void qsieve_linalg_clear(qs_t qs_inf);

int qsieve_relations_cmp(const void * a, const void * b);

slong qsieve_merge_relations(qs_t qs_inf);

slong qsieve_insert_relation(qs_t qs_inf, fmpz_t Y);

void qsieve_write_to_file(qs_t qs_inf, mp_limb_t prime, fmpz_t Y, qs_poly_t poly);

hash_t * qsieve_get_table_entry(qs_t qs_inf, mp_limb_t prime);

void qsieve_add_to_hashtable(qs_t qs_inf, mp_limb_t prime);

relation_t qsieve_parse_relation(qs_t qs_inf, char * str);

relation_t qsieve_merge_relation(qs_t qs_inf, relation_t  a, relation_t  b);

int qsieve_compare_relation(const void * a, const void * b);

int qsieve_remove_duplicates(relation_t * rel_list, slong num_relations);

void qsieve_insert_relation2(qs_t qs_inf, relation_t * rel_list,
                                                          slong num_relations);

int qsieve_process_relation(qs_t qs_inf);

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
