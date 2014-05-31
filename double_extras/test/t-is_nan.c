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

   Copyright (C) 2012 Fredrik Johansson
   Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "ulong_extras.h"
#include "double_extras.h"

int
main(void)
{
    double x;
    slong iter;

    FLINT_TEST_INIT(state);

    flint_printf("is_nan....");
    fflush(stdout);

    /* check non-zero value returned if x == NaN */
    x = D_NAN;
    if (!d_is_nan(x))
    {
        flint_printf("FAIL\n");
        flint_printf("0 returned for %g\n", x);
        abort();
    }

    /* check 0 returned if x != NaN */
    x = D_INF;
    if (d_is_nan(x))
    {
        flint_printf("FAIL\n");
        flint_printf("Non-zero returned for %g\n", x);
        abort();
    }
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest_signed(state, 0, 0);
        if (d_is_nan(x))
        {
            flint_printf("FAIL\n");
            flint_printf("Non-zero returned for %g\n", x);
            abort();
        }
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
