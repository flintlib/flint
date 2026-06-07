/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef QSIEVE_H
#define QSIEVE_H

#include <stdint.h>
#include "fmpz_types.h"

#if FLINT_USES_PTHREAD
# include <pthread.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Windows systems may define `small` macro, which leads to conflicts */
#undef small

#define QS_DEBUG 0 /* level of debug information printed, 0 = none */

#define BITS_ADJUST 25 /* add to sieve entries to compensate approximations */

#define BLOCK_SIZE (4*65536) /* size of sieving cache block */

typedef struct
{
   ulong pinv;     /* precomputed inverse */
   ulong pinv2;    /* alternative precomputed inverse */
   int p;              /* prime */
   char size;
} prime_t;

typedef struct          /* struct for factors of relations */
{
   slong ind;
   slong exp;
} fac_t;

typedef struct           /* matrix column */
{
   slong * data;		/* The list of occupied rows in this column */
   slong weight;		/* Number of nonzero entries in this column */
   slong orig;         /* Original relation number */
} la_col_t;


typedef struct          /* entry in hash table */
{
   ulong prime;    /* value of prime */
   ulong next;     /* next prime which have same hash value as 'prime' */
   ulong count;    /* number of occurrence of 'prime' */
} hash_t;

typedef struct             /* format for relation */
{
   ulong lp;          /* large prime, is 1, if relation is full */
   slong num_factors;     /* number of factors, excluding small factor */
   slong small_primes;   /* number of small factors */
   slong * small;         /* exponent of small factors */
   fac_t * factor;        /* factor of relation */
   fmpz_t Y;              /* square root of sieve value for relation */
} relation_t;

typedef struct
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

typedef struct
{
   volatile slong index_j;
#if FLINT_USES_PTHREAD
   pthread_mutex_t mutex;
#endif
   thread_pool_handle * handles;
   slong num_handles;

   fmpz_t n;               /* Number to factor */

   flint_bitcnt_t bits;    /* Number of bits of n */

   ulong ks_primes;        /* number of Knuth-Schroeppel primes */

   slong fb_primes;        /* number of factor base primes to use (incl. k and 2) */

   ulong k;            /* Multiplier */
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

   fmpz_t A;                /* current value of coeff A of poly Ax^2 + Bx + C */
   fmpz_t B;                /* B value of poly */

   ulong * A_ind;       /* indices of factor base primes dividing A */

   fmpz_t * A_divp;         /* A_divp[i] = (A/p_i),
                               where the p_i are the prime factors of A */
   ulong * B0_terms;    /* B0_terms[i] = min(gamma_i, p - gamma_i) where
                               gamma_i = (sqrt(kn)*(A_divp[i])^(-1)) mod p_i,
                               where the p_i are the prime factors of A */

   fmpz_t * B_terms;        /* B_terms[i] = A_divp[i]*B0_terms[i] (multprec) */

   ulong * A_inv;       /* A_inv[k] = A^(-1) mod p_k, for FB prime p_k */
   ulong ** A_inv2B;    /* A_inv2B[i][k] = 2 * B_terms[i] * A^(-1) mod p_k
                               for FB prime p_k */

   int * soln1;             /* soln1[k] = first poly root mod FB prime p_k */
   int * soln2;             /* soln2[k] = second poly root mod FB prime p_k */

   fmpz_t target_A;         /* approximate target value for A coeff of poly */

   fmpz_t upp_bound;        /* upper bound on desired A = 2*target_A */
   fmpz_t low_bound;        /* lower bound on desired A = target_A/2 */

   slong s;                 /* number of prime factors of A */
   slong low;               /* minimum offset in factor base,
                               for possible factors of A */
   slong high;              /* maximum offset in factor base,
                               for possible factors of A */
   slong span;              /* total number of possible factors for A */

   /*
      parameters for calculating next subset of possible factors of A, giving
      a lexicographic ordering of all such tuples
   */
   slong h; /* tuple entry we just set, numbered from 1 at end of tuple */
   slong m; /* last value we just set a tuple entry to */
   slong A_ind_diff; /* diff. between indices of (s-1) and (s-2)-th A-factor */
   ulong * curr_subset; /* current tuple */
   ulong * first_subset; /* first tuple, in case of restart */
   ulong j; /* index of s-th factor of first A, if s > 3 */

#if QS_DEBUG
   slong poly_count;         /* keep track of the number of polynomials used */
#endif

   qs_poly_s * poly;         /* poly data per thread */

   /***************************************************************************
                       RELATION DATA
   ***************************************************************************/

   FLINT_FILE * siqs;           /* pointer to file for storing relations */
   char * fname;          /* name of file used for relations */

   slong full_relation;   /* number of full relations */
   slong num_cycles;      /* number of possible full relations from partials */

   slong vertices;        /* number of different primes in partials */
   slong components;      /* equal to 1 */
   slong edges;           /* total number of partials */

   slong table_size;      /* size of table */
   hash_t * table;        /* store 'prime' occurring in partial */
   ulong * hash_table;  /* to keep track of location of primes in 'table' */

   slong extra_rels;      /* number of extra relations beyond num_primes */
   slong max_factors;     /* maximum number of factors a relation can have */

   fmpz * Y_arr;          /* array of Y's corresponding to relations */
   slong * curr_rel;      /* current relation in array of relations */
   slong * relation;      /* relation array */

   slong buffer_size;     /* size of buffer of relations */
   slong num_relations;   /* number of relations so far */

   ulong small_factor;    /* small factor found when merging relations */

   /***************************************************************************
                       LINEAR ALGEBRA DATA
   ***************************************************************************/

   la_col_t * matrix; /* the main matrix over GF(2) in sparse format */
   la_col_t ** qsort_arr; /* array of columns ready to be sorted */

   slong columns; /* number of columns in matrix so far */

   /***************************************************************************
                       SQUARE ROOT DATA
   ***************************************************************************/

   slong * prime_count; /* counts of the exponents of primes appearing in the square */

} qs_s;

typedef qs_s qs_t[1];

/*
   Tuning parameters { bits, ks_primes, fb_primes, small_primes, sieve_size}
   for qsieve_factor_threaded where:
     * bits is the number of bits of n
     * ks_primes is the max number of primes to try in Knuth-Schroeppel function
     * fb_primes is the number of factor base primes to use (including k and 2)
     * small_primes is the number of small primes to not factor with (including k and 2)
     * sieve_size is the size of the sieve to use
     * sieve_bits - sieve_fill
*/

static const ulong qsieve_tune[][6] =
{
    /* Generated with qsieve-tune */
    {10,  49,     84,  4,    4300,  25},  /* 3 digits, 0.00 ms/semiprime */
    {20,  10,    120,  5,    4800,  26},  /* 6 digits, 0.01 ms/semiprime */
    {30,  14,     25,  3,    4500,  19},  /* 9 digits, 0.72 ms/semiprime */
    {40,  28,     73,  9,    5700,  24},  /* 12 digits, 1.05 ms/semiprime */
    {50,  25,     85,  7,    5400,  28},  /* 15 digits, 1.25 ms/semiprime */
    {60,  27,    100, 14,   11000,  28},  /* 18 digits, 1.45 ms/semiprime */
    {65,  26,    110, 19,   13000,  28},  /* 19 digits, 1.86 ms/semiprime */
    {70,  26,    180,  3,   12000,  43},  /* 21 digits, 2.26 ms/semiprime */
    {75,  26,    190,  4,   12000,  45},  /* 22 digits, 2.61 ms/semiprime */
    {80,  27,    190,  5,   12000,  45},  /* 24 digits, 3.16 ms/semiprime */
    {85,  28,    150,  5,   10000,  45},  /* 25 digits, 2.34 ms/semiprime */
    {90,  34,    150,  5,   10000,  45},  /* 27 digits, 2.68 ms/semiprime */
    {100, 40,    210,  6,   14000,  48},  /* 30 digits, 4.20 ms/semiprime */
    {110, 46,    220,  7,   13000,  48},  /* 33 digits, 6.50 ms/semiprime */
    {120, 50,    250,  7,   13000,  50},  /* 36 digits, 11.58 ms/semiprime */
    {130, 47,    310,  9,   16000,  50},  /* 39 digits, 19.37 ms/semiprime */
    {135, 39,    450, 10,   20000,  53},  /* 40 digits, 25.39 ms/semiprime */
    {140, 47,    580,  8,   17000,  56},  /* 41 digits, 33.47 ms/semiprime */
    {145, 47,    610,  9,   16000,  56},  /* 43 digits, 46.47 ms/semiprime */
    {150, 45,    610,  8,   17000,  59},  /* 45 digits, 59.05 ms/semiprime */
    {155, 43,    860, 10,   21000,  59},  /* 46 digits, 92.20 ms/semiprime */
    /* Legacy tuning values that seem reasonable */
    {160, 150,  2000, 11,   4 *  65536, 73}, /* 49 digit */
    {170, 150,  2000, 12,   4 *  65536, 75}, /* 52 digits */
    {180, 150,  3000, 12,   4 *  65536, 76}, /* 55 digits */
    {190, 150,  3000, 13,   4 *  65536, 78}, /* 58 digit */
    {200, 200,  4500, 14,   4 *  65536, 81}, /* 61 digits */
    {210, 100,  8000, 14,   12 *  65536, 84}, /* 64 digits */
    {220, 300, 10000, 15,   12 *  65536, 88}, /* 67 digits */
    {230, 400, 20000, 17,   20 *  65536, 90}, /* 70 digits */
    {240, 450, 20000, 19,   20 *  65536, 93}, /* 73 digis */
    {250, 500, 22000, 22,   24 *  65536, 97}, /* 76 digits */
    {260, 600, 25000, 25,   24 *  65536, 100}, /* 79 digits */
    {270, 800, 35000, 27,   28 *  65536, 102}, /* 82 digits */
    {280, 900, 40000, 29,   28 *  65536, 104}, /* 85 digits */
    {290, 1000, 60000, 29,  32 *  65536, 106}, /* 88 digits */
    {300, 1100, 140000, 30,  32 * 65536, 108} /* 91 digits */
};

/* number of entries in the tuning table */
#define QS_TUNE_SIZE (sizeof(qsieve_tune)/(6*sizeof(ulong)))

void qsieve_init(qs_t qs_inf, const fmpz_t n);

/* As qsieve_init, but with the tuning parameters supplied explicitly rather
   than read from the qsieve_tune table.  qsieve_init calls this with the
   table defaults for the bit-size of n. */
void qsieve_init_with_tune(qs_t qs_inf, const fmpz_t n, ulong ks_primes,
        slong fb_primes, slong small_primes, slong sieve_size, ulong sieve_bits);

ulong qsieve_knuth_schroeppel(qs_t qs_inf);

void qsieve_clear(qs_t qs_inf);

void qsieve_factor(fmpz_factor_t factors, const fmpz_t n);

/* As qsieve_factor, but with the tuning parameters supplied explicitly. */
void qsieve_factor_with_tune(fmpz_factor_t factors, const fmpz_t n,
        ulong ks_primes, slong fb_primes, slong small_primes,
        slong sieve_size, ulong sieve_bits);

prime_t * compute_factor_base(ulong * small_factor, qs_t qs_inf,
                                                             slong num_primes);

ulong qsieve_primes_init(qs_t qs_inf);

ulong qsieve_primes_increment(qs_t qs_inf, ulong delta);

ulong qsieve_poly_init(qs_t qs_inf);

int qsieve_init_A(qs_t qs_inf);

void qsieve_reinit_A(qs_t qs_inf);

int qsieve_next_A(qs_t qs_inf);

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

void qsieve_linalg_realloc(qs_t qs_inf);

void qsieve_linalg_clear(qs_t qs_inf);

int qsieve_relations_cmp(const void * a, const void * b);

slong qsieve_merge_relations(qs_t qs_inf);

void qsieve_write_to_file(qs_t qs_inf, ulong prime,
                                                     const fmpz_t Y, const qs_poly_t poly);

hash_t * qsieve_get_table_entry(qs_t qs_inf, ulong prime);

void qsieve_add_to_hashtable(qs_t qs_inf, ulong prime);

relation_t qsieve_parse_relation(qs_t qs_inf);

relation_t qsieve_merge_relation(qs_t qs_inf, relation_t  a, relation_t  b);

int qsieve_compare_relation(const void * a, const void * b);

int qsieve_remove_duplicates(relation_t * rel_list, slong num_relations);

void qsieve_insert_relation(qs_t qs_inf, relation_t * rel_list,
                                                          slong num_relations);

int qsieve_process_relation(qs_t qs_inf);

static inline void insert_col_entry(la_col_t * col, slong entry)
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

static inline void swap_cols(la_col_t * col2, la_col_t * col1)
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

static inline void clear_col(la_col_t * col)
{
   col->weight = 0;
}

static inline void free_col(la_col_t * col)
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
