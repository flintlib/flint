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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

void qsieve_ll_linalg_init(qs_t qs_inf)
{
    long i;
    
    qs_inf->extra_rels = 64; /* number of opportunities to factor n */
    qs_inf->max_factors = 30; /* maximum number of factors a relation can have */

    /* allow as many dups as relations */
    qs_inf->buffer_size = 2*(qs_inf->num_primes + qs_inf->extra_rels + qs_inf->qsort_rels);

    qs_inf->small = flint_malloc(qs_inf->small_primes*sizeof(mp_limb_t));
    qs_inf->factor = flint_malloc(qs_inf->max_factors*sizeof(fac_t));
    qs_inf->matrix = flint_malloc((qs_inf->buffer_size + qs_inf->qsort_rels)*sizeof(la_col_t));
    qs_inf->unmerged = qs_inf->matrix + qs_inf->buffer_size;
    qs_inf->Y_arr = flint_malloc(qs_inf->buffer_size*sizeof(fmpz));
    qs_inf->curr_rel = qs_inf->relation
                     = flint_malloc(2*qs_inf->buffer_size*qs_inf->max_factors*sizeof(long));
    qs_inf->qsort_arr = flint_malloc(qs_inf->qsort_rels*sizeof(la_col_t *));

    for (i = 0; i < qs_inf->buffer_size; i++)
    {
        fmpz_init(qs_inf->Y_arr + i);
        qs_inf->matrix[i].weight = 0;
        qs_inf->matrix[i].data = NULL;
    }

    for (i = 0; i < qs_inf->qsort_rels; i++)
    {
        qs_inf->unmerged[i].weight = 0;
        qs_inf->unmerged[i].data = NULL;
    }
    
    qs_inf->prime_count = flint_malloc(qs_inf->num_primes*sizeof(long));

    qs_inf->num_unmerged = 0;
    qs_inf->columns = 0;
    qs_inf->num_relations = 0;
}
