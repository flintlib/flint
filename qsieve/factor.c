/*
    Copyright (C) 2006, 2011, 2016, 2020 William Hart
    Copyright (C) 2015 Nitin Kumar
    Copyright (C) 2020 Dan Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qsieve.h"
#include "fmpz_factor.h"
#include "thread_support.h"

#include <inttypes.h>
#define _STDC_FORMAT_MACROS
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)
#include <unistd.h>
#endif
#if defined (__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#endif

int compare_facs(const void * a, const void * b)
{
   fmpz * x = (fmpz *) a;
   fmpz * y = (fmpz *) b;

   return fmpz_cmp(x, y);
}

/*
   Finds at least one nontrivial factor of n using the self initialising
   multiple polynomial quadratic sieve with single large prime variation.
   Assumes n is not prime and not a perfect power.
*/
void qsieve_factor(fmpz_factor_t factors, const fmpz_t n)
{
    qs_t qs_inf;
    mp_limb_t small_factor, delta;
    ulong expt = 0;
    unsigned char * sieve;
    slong ncols, nrows, i, j = 0, count, num_primes;
    uint64_t * nullrows = NULL;
    uint64_t mask;
    flint_rand_t state;
    fmpz_t temp, temp2, X, Y;
    slong num_facs;
    fmpz * facs;
    int nchars;

    if (fmpz_sgn(n) < 0)
    {
       fmpz_t n2;

       fmpz_init(n2);
       fmpz_abs(n2, n);

       factors->sign *= -1;
       
       qsieve_factor(factors, n2);

       fmpz_clear(n2);
       
       return;
    }

    /**************************************************************************
        INITIALISATION:
        Initialise the qs_t structure.
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nstart\n");
#endif

    qsieve_init(qs_inf, n);

#if QS_DEBUG
    flint_printf("factoring ");
    fmpz_print(qs_inf->n);
    flint_printf(" of %wu bits\n", qs_inf->bits);
#endif

    /**************************************************************************
        KNUTH SCHROEPPEL:
        Try to compute a multiplier k such that there are a lot of small primes
        which are quadratic residues modulo kn. If a small factor of n is found
        during this process it is returned.
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nKnuth_Schroeppel\n");
#endif

    small_factor = qsieve_knuth_schroeppel(qs_inf);

    if (small_factor)
    {

#if QS_DEBUG
        flint_printf("found small factor %wu in Knuth-Schroeppel\n", small_factor);
#endif
        fmpz_init_set_ui(temp, small_factor);

        expt += fmpz_remove(temp, qs_inf->n, temp);

        _fmpz_factor_append_ui(factors, small_factor, expt);
        
        qsieve_clear(qs_inf);

        fmpz_factor_no_trial(factors, temp);

        fmpz_clear(temp);
        
        return;
    }

    /* compute kn */
    fmpz_mul_ui(qs_inf->kn, qs_inf->n, qs_inf->k);

    /* refine qs_inf->bits */
    qs_inf->bits = fmpz_bits(qs_inf->kn);

#if QS_DEBUG
    flint_printf("kn bits = %wd\n", qs_inf->bits);
#endif

    /**************************************************************************
        COMPUTE FACTOR BASE:
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nCompute factor-base\n");
#endif

    /* compute factor base primes and associated data */
    small_factor = qsieve_primes_init(qs_inf);

    if (small_factor)
    {

#if QS_DEBUG
        flint_printf("found small factor %wu while generating factor base\n", small_factor);
#endif

        fmpz_init_set_ui(temp, small_factor);

        expt += fmpz_remove(temp, qs_inf->n, temp);

        _fmpz_factor_append_ui(factors, small_factor, expt);
        
        qsieve_clear(qs_inf);

        fmpz_factor_no_trial(factors, temp);

        fmpz_clear(temp);

        return;
    }

    fmpz_init(temp);
    fmpz_init(temp2);
    fmpz_init(X);
    fmpz_init(Y);

    /**************************************************************************
        INITIALISE RELATION/LINEAR ALGEBRA DATA:
        Create space for all the relations and matrix information
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nInitializing Relations and Linear Algebra\n");
#endif

    qsieve_linalg_init(qs_inf);

    /**************************************************************************
        POLYNOMIAL INITIALIZATION AND SIEVING:
        Sieve for relations
    **************************************************************************/

#if QS_DEBUG
    flint_printf("\nPolynomial Initialisation and Sieving\n");
#endif

    qs_inf->num_handles = flint_request_threads(&qs_inf->handles, flint_get_num_threads());

    /* ensure cache lines don't overlap if num_handles > 0 */
    sieve = flint_malloc((qs_inf->sieve_size + sizeof(ulong)
               + (qs_inf->num_handles > 0 ? 64 : 0))*(qs_inf->num_handles + 1));

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&qs_inf->mutex, NULL);
#endif
    
#if defined (__WIN32) && !defined(__CYGWIN__)
    srand((int) GetCurrentProcessId());
#else
    srand((int) getpid());
#endif
    nchars = sprintf(qs_inf->fname, "%d", (int) rand());
    strcat(qs_inf->fname + nchars, "siqs.dat");

    qs_inf->siqs = fopen(qs_inf->fname, "w");

    for (j = qs_inf->small_primes; j < qs_inf->num_primes; j++)
    {
        if (qs_inf->factor_base[j].p > BLOCK_SIZE)
            break;
    }

    qs_inf->second_prime = j;

#if QS_DEBUG
    flint_printf("second prime index = %wd\n", qs_inf->second_prime);
#endif

    while (1)
    {
        if (qs_inf->s) /* we have already tried factoring, so restart */
            qsieve_reinit_A(qs_inf);
        else
        {
            if (!qsieve_init_A(qs_inf))
                goto more_primes; /* initialisation failed, increase FB */
        }

        do
        {           
            qsieve_collect_relations(qs_inf, sieve);
                
            qs_inf->num_cycles = qs_inf->edges + qs_inf->components - qs_inf->vertices;

#if QS_DEBUG
            flint_printf("full relations = %wd, num cycles = %wd, ks_primes = %wd, "
                         "extra rels = %wd, poly_count = %wd, num_primes = %wd\n", qs_inf->full_relation,
                          qs_inf->num_cycles, qs_inf->ks_primes,
                          qs_inf->extra_rels, qs_inf->poly_count, qs_inf->num_primes);
#endif
 
            if (qs_inf->full_relation + qs_inf->num_cycles >= 
                ((slong) (1.10*qs_inf->num_primes) + qs_inf->ks_primes + qs_inf->extra_rels))
            {
                int ok;

                fclose(qs_inf->siqs);

                ok = qsieve_process_relation(qs_inf);

                if (ok == -1)
                {
                    small_factor = qs_inf->small_factor;
#if QS_DEBUG
                    flint_printf("found small factor %wu while incrementing factor base\n, small_factor");
#endif
                    goto found_small_factor;
                }

                if (ok)
                {

    /**************************************************************************
        REDUCE MATRIX:
        Perform some light filtering on the matrix
    **************************************************************************/

                    num_primes = qs_inf->num_primes;
                    qs_inf->num_primes += qs_inf->ks_primes;

                    ncols = qs_inf->num_primes + qs_inf->extra_rels;
                    nrows = qs_inf->num_primes;

                    reduce_matrix(qs_inf, &nrows, &ncols, qs_inf->matrix);


   /**************************************************************************
        BLOCK LANCZOS:
        Find extra_rels nullspace vectors (if they exist)
    **************************************************************************/

#if QS_DEBUG
                    flint_printf("\nBlock Lanczos\n");
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

    /**************************************************************************
        SQUARE ROOT:
        Compute the square root and take the GCD of X-Y with N
    **************************************************************************/

#if QS_DEBUG
                    flint_printf("\nSquare Root\n");
                    flint_printf("Found %ld kernel vectors\n", count);
#endif

                    facs = _fmpz_vec_init(100);
                    num_facs = 0;

                    for (i = 0; i < 64; i++)
                    {
                        if (mask & ((uint64_t)(1) << i))
                        {
                            qsieve_square_root(X, Y, qs_inf, nullrows, ncols, i, qs_inf->kn);

                            fmpz_sub(X, X, Y);
                            fmpz_gcd(X, X, qs_inf->n);

                            if (fmpz_cmp(X, qs_inf->n) != 0 && fmpz_cmp_ui(X, 1) != 0) /* have a factor */
                                fmpz_set(facs + num_facs++, X);
                        }
                    }

                    flint_free(nullrows);

                    if (num_facs > 0)
                    {
                        _fmpz_factor_append(factors, qs_inf->n, 1);

                        qsort((void *) facs, num_facs, sizeof(fmpz), compare_facs);

                        for (i = 0; i < num_facs; i++)
                        {
                            fmpz_gcd(temp, factors->p + factors->num - 1, facs + i);
                            if (!fmpz_is_one(temp))
                            {
                                factors->exp[factors->num - 1] = fmpz_remove(temp2, factors->p + factors->num - 1, temp);
                                fmpz_set(factors->p + factors->num - 1, temp);
                                
                                if (fmpz_is_one(temp2))
                                   break;
                                else
                                   _fmpz_factor_append(factors, temp2, 1);
                             }  
                        }

                        _fmpz_vec_clear(facs, 100);

                        goto cleanup;
                    }

                    _fmpz_vec_clear(facs, 100);

                    qs_inf->siqs = fopen(qs_inf->fname, "w");
                    qs_inf->num_primes = num_primes; /* linear algebra adjusts this */
                    goto more_primes; /* factoring failed, may need more primes */
                }
            }
        } while (qsieve_next_A(qs_inf));

more_primes: /* ran out of A's in init/sieving of linalg failed, increase FB */

#if QS_DEBUG
        printf("Increasing factor base.\n");
#endif

        delta = qs_inf->num_primes / 10;
        delta = FLINT_MAX(delta, 100); /* add at least 100 more primes */
        
#if QS_DEBUG
        flint_printf("\nfactor base increment\n");
#endif
        qsieve_poly_clear(qs_inf);

        small_factor = qsieve_primes_increment(qs_inf, delta);

        for (j = qs_inf->small_primes; j < qs_inf->num_primes; j++)
        {
            if (qs_inf->factor_base[j].p > BLOCK_SIZE)
               break;
        }

        qs_inf->second_prime = j;

        qs_inf->s = 0; /* indicate polynomials need setting up again */

#if QS_DEBUG
        printf("Now %ld primes\n", qs_inf->num_primes);
#endif

        if (small_factor)
        {

#if QS_DEBUG
            flint_printf("found small factor %wu while incrementing factor base\n, small_factor");
#endif

found_small_factor:

            fmpz_set_ui(temp, small_factor);

            expt += fmpz_remove(temp, qs_inf->n, temp);

            _fmpz_factor_append_ui(factors, small_factor, expt);
        
            fmpz_factor_no_trial(factors, temp);
            
            goto cleanup;
        }

        qsieve_linalg_realloc(qs_inf);
    }

    /**************************************************************************
        CLEANUP:
        Clean up allocated memory
    **************************************************************************/

cleanup:

#if QS_DEBUG
    flint_printf("\nCleanup\n");
#endif

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&qs_inf->mutex);
#endif

    flint_give_back_threads(qs_inf->handles, qs_inf->num_handles);

    flint_free(sieve);
    remove(qs_inf->fname);
    qsieve_clear(qs_inf);
    qsieve_linalg_clear(qs_inf);
    qsieve_poly_clear(qs_inf);
    fmpz_clear(X);
    fmpz_clear(Y);
    fmpz_clear(temp);
    fmpz_clear(temp2);
}
