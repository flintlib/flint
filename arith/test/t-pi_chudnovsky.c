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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"


int main(void)
{
    len_t k;

    printf("pi_chudnovsky....");
    fflush(stdout);

    for (k = 2; k < 20; k++)
    {
        mpfr_t x, y;

        mpfr_init2(x, 1L << k);
        mpfr_init2(y, 1L << k);

        mpfr_const_pi(x, MPFR_RNDN);
        mpfr_pi_chudnovsky(y, MPFR_RNDN);

        if (!mpfr_equal_p(x, y))
        {
            printf("FAIL:\n");
            printf("Wrong value at prec = %ld\n", 1L << k);
            abort();
        }

        mpfr_clear(x);
        mpfr_clear(y);
    }

    mpfr_free_cache();
    printf("PASS\n");
    return 0;
}
