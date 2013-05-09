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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "fmpz.h"
#include "ulong_extras.h"


int main(void)
{
    len_t k;
    flint_rand_t state;

    printf("zeta_ui_bsplit....");
    fflush(stdout);

    flint_randinit(state);

    for (k = 0; k < 100; k++)
    {
        len_t prec, n;
        mpfr_t x, y;

        n = 2 + n_randint(state, 20);
        prec = 2 + n_randint(state, 10000);

        mpfr_init2(x, prec);
        mpfr_init2(y, prec);

        mpfr_zeta_ui(x, n, MPFR_RNDN);
        mpfr_zeta_ui_bsplit(y, n, MPFR_RNDN);

        if (!mpfr_equal_p(x, y))
        {
            printf("FAIL:\n");
            printf("Wrong value at n = %ld, prec = %ld\n", n, prec);
            printf("x:\n");
            mpfr_out_str(stdout, 10, 0, x, MPFR_RNDN);
            printf("\n");
            printf("x:\n");
            mpfr_out_str(stdout, 10, 0, y, MPFR_RNDN);
            printf("\n");
            abort();
        }

        mpfr_clear(x);
        mpfr_clear(y);
    }

    flint_randclear(state);
    mpfr_free_cache();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
