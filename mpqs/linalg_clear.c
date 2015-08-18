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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpqs.h"

void
mpqs_linalg_clear(mpqs_t mpqs_inf)
{
    slong i;

<<<<<<< HEAD:mpqs/linalg_clear.c
    flint_free(mpqs_inf->small);
    flint_free(mpqs_inf->factor);
    flint_free(mpqs_inf->relation);
    flint_free(mpqs_inf->qsort_arr);
=======
    flint_free(qs_inf->small);
    flint_free(qs_inf->factor);
    flint_free(qs_inf->relation);
    flint_free(qs_inf->qsort_arr);
    flint_free(qs_inf->hash_table);
    flint_free(qs_inf->table);
>>>>>>> 5903f39ba98692427f9f0b120c0ac156fb229a3c:siqs/linalg_clear.c

    if (mpqs_inf->matrix != NULL)
    {
        for (i = 0; i < mpqs_inf->buffer_size + mpqs_inf->num_unmerged; i++)
        {
            la_col_t * col = mpqs_inf->matrix + i;

            if (col->weight)
                flint_free(col->data);
        }

        flint_free(mpqs_inf->matrix);
    }

    if (mpqs_inf->Y_arr != NULL)
    {
        for (i = 0; i < mpqs_inf->buffer_size; i++)
            fmpz_clear(mpqs_inf->Y_arr + i);
        flint_free(mpqs_inf->Y_arr);
    }

    flint_free(mpqs_inf->prime_count);

<<<<<<< HEAD:mpqs/linalg_clear.c
    mpqs_inf->small = NULL;
    mpqs_inf->factor = NULL;
    mpqs_inf->relation = NULL;
    mpqs_inf->qsort_arr = NULL;
    mpqs_inf->matrix = NULL;
    mpqs_inf->Y_arr = NULL;
    mpqs_inf->prime_count = NULL;
=======
    qs_inf->small = NULL;
    qs_inf->factor = NULL;
    qs_inf->relation = NULL;
    qs_inf->qsort_arr = NULL;
    qs_inf->matrix = NULL;
    qs_inf->Y_arr = NULL;
    qs_inf->prime_count = NULL;
    qs_inf->hash_table = NULL;
    qs_inf->table = NULL;
>>>>>>> 5903f39ba98692427f9f0b120c0ac156fb229a3c:siqs/linalg_clear.c
}
