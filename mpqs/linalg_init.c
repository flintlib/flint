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
    Built upon existing FLINT qsieve
    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpqs.h"

void
mpqs_linalg_init(mpqs_t mpqs_inf)
{
    slong i, num_primes;

    mpqs_inf->extra_rels = 64; /* number of opportunities to factor n */
    mpqs_inf->max_factors = 30; /* maximum number of factors a relation can have */

    /* allow as many dups as relations */
    num_primes = mpqs_inf->num_primes;
    mpqs_inf->num_primes += mpqs_inf->ks_primes;
    mpqs_inf->buffer_size = 2 * (mpqs_inf->num_primes + mpqs_inf->extra_rels +
                                 mpqs_inf->qsort_rels);
    mpqs_inf->small     = flint_malloc(mpqs_inf->small_primes * sizeof(mp_limb_t));
    mpqs_inf->factor    = flint_malloc(mpqs_inf->max_factors * sizeof(fac_t));
    mpqs_inf->matrix    = flint_malloc((mpqs_inf->buffer_size + mpqs_inf->qsort_rels)
                                        * sizeof(la_col_t));
    mpqs_inf->unmerged  = mpqs_inf->matrix + mpqs_inf->buffer_size;
    mpqs_inf->Y_arr     = flint_malloc(mpqs_inf->buffer_size * sizeof(fmpz));
    mpqs_inf->curr_rel  = mpqs_inf->relation
                        = flint_malloc(2 * mpqs_inf->buffer_size *
                                       mpqs_inf->max_factors * sizeof(slong));
    mpqs_inf->qsort_arr = flint_malloc(mpqs_inf->qsort_rels * sizeof(la_col_t *));

    for (i = 0; i < mpqs_inf->buffer_size; i++)
    {
        fmpz_init(mpqs_inf->Y_arr + i);
        mpqs_inf->matrix[i].weight = 0;
        mpqs_inf->matrix[i].data = NULL;
    }

    for (i = 0; i < mpqs_inf->qsort_rels; i++)
    {
        mpqs_inf->unmerged[i].weight = 0;
        mpqs_inf->unmerged[i].data = NULL;
    }

    mpqs_inf->prime_count = flint_malloc(mpqs_inf->num_primes*sizeof(slong));

    mpqs_inf->num_primes = num_primes;
    mpqs_inf->num_unmerged = 0;
    mpqs_inf->columns = 0;
    mpqs_inf->num_relations = 0;
}

/* re-initialize all the linear algebra parameter */

void
mpqs_linalg_re_init(mpqs_t mpqs_inf)
{
    slong i;
    mpqs_inf->curr_rel = mpqs_inf->relation;

    for (i = 0; i < mpqs_inf->buffer_size; i++)
    {
        fmpz_init(mpqs_inf->Y_arr + i);
        mpqs_inf->matrix[i].weight = 0;
        mpqs_inf->matrix[i].data = NULL;
    }

    for (i = 0; i < mpqs_inf->qsort_rels; i++)
    {
        mpqs_inf->unmerged[i].weight = 0;
        mpqs_inf->unmerged[i].data = NULL;
    }

    mpqs_inf->num_unmerged = 0;
    mpqs_inf->columns = 0;
    mpqs_inf->num_relations = 0;
}

/* increase size of different array after factor base increment*/

void
mpqs_linalg_re_alloc(mpqs_t mpqs_inf)
{
    slong num_primes;

    num_primes = mpqs_inf->num_primes;
    mpqs_inf->num_primes += mpqs_inf->ks_primes;
    mpqs_inf->buffer_size = 2 * (mpqs_inf->num_primes + mpqs_inf->extra_rels +
                                 mpqs_inf->qsort_rels);
    mpqs_inf->matrix = flint_realloc(mpqs_inf->matrix, (mpqs_inf->buffer_size + 
                                     mpqs_inf->qsort_rels) * sizeof(la_col_t));
    mpqs_inf->unmerged = mpqs_inf->matrix + mpqs_inf->buffer_size;
    mpqs_inf->Y_arr = flint_realloc(mpqs_inf->Y_arr, mpqs_inf->buffer_size
                                    * sizeof(fmpz));
    mpqs_inf->curr_rel = mpqs_inf->relation
                     = flint_realloc(mpqs_inf->relation, 2 * mpqs_inf->buffer_size
                                     * mpqs_inf->max_factors * sizeof(slong));
    mpqs_inf->prime_count = flint_realloc(mpqs_inf->prime_count,
                                          mpqs_inf->num_primes * sizeof(slong));
    mpqs_inf->num_primes = num_primes;
    mpqs_linalg_init(mpqs_inf);
}
