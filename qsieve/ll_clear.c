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

void qsieve_ll_clear(qs_t qs_inf)
{
    long i;
    
    fmpz_clear(qs_inf->kn);
    fmpz_clear(qs_inf->C);
   
    flint_free(qs_inf->factor_base);
    flint_free(qs_inf->sqrts);
    flint_free(qs_inf->B_terms);
    flint_free(qs_inf->A_inv);
    
    if (qs_inf->A_inv2B != NULL)
        flint_free(qs_inf->A_inv2B[0]);

    flint_free(qs_inf->A_inv2B);
    
    qs_inf->factor_base = NULL;
    qs_inf->sqrts       = NULL;
    qs_inf->B_terms     = NULL;
    qs_inf->A_inv       = NULL;
    qs_inf->A_inv2B     = NULL;

    flint_free(qs_inf->small);
    flint_free(qs_inf->factor);
    flint_free(qs_inf->relation);
    flint_free(qs_inf->qsort_arr);
    
    if (qs_inf->matrix != NULL)
    {
        for (i = 0; i < qs_inf->buffer_size + qs_inf->num_unmerged; i++)
        {
            la_col_t * col = qs_inf->matrix + i;
            if (col->weight)
                flint_free(col->data);
        }

        flint_free(qs_inf->matrix);
    }
    
    if (qs_inf->Y_arr != NULL)
    {
        for (i = 0; i < qs_inf->buffer_size; i++)
            fmpz_clear(qs_inf->Y_arr + i);    
        flint_free(qs_inf->Y_arr);
    }

    flint_free(qs_inf->prime_count);
     
    qs_inf->small       = NULL;
    qs_inf->factor      = NULL;
    qs_inf->matrix      = NULL;
    qs_inf->Y_arr       = NULL;
    qs_inf->relation    = NULL;
    qs_inf->qsort_arr   = NULL;
    
    qs_inf->prime_count = NULL;

    qs_inf->A = 0;

#if (QS_DEBUG & 16)
    flint_free(qs_inf->sieve_tally);
    qs_inf->sieve_tally = NULL;
#endif
}
