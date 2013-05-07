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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("init/clear....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        mpfr_mat_t a;
        long j, k;
        long rows = n_randint(state, 100);
        long cols = n_randint(state, 100);
        mp_prec_t prec = n_randint(state, 200) + MPFR_PREC_MIN;

        mpfr_mat_init(a, rows, cols, prec);

        for (j = 0; j < rows; j++)
            for (k = 0; k < cols; k++)
                mpfr_set_ui(a->rows[j] + k, 0, GMP_RNDN);

        mpfr_mat_clear(a);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
