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

    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/


#ifdef T

#include "templates.h"

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("concat_horizontal....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
    	TEMPLATE(T, ctx_t) ctx;
    	TEMPLATE(T, mat_t) A, B, C;
    	TEMPLATE(T, mat_t) window1, window2;
        slong c1, c2, r1;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        c1 = n_randint(state, 50);
        c2 = n_randint(state, 50);
        r1 = n_randint(state, 50);

        TEMPLATE(T, mat_init) (A, r1, c1, ctx);
        TEMPLATE(T, mat_init) (B, r1, c2, ctx);
        TEMPLATE(T, mat_init) (C, r1, (c1 + c2), ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);
        
        TEMPLATE(T, mat_randtest) (C, state, ctx);

        TEMPLATE(T, mat_concat_horizontal) (C, A, B, ctx);

        TEMPLATE(T, mat_window_init) (window1, C, 0, 0, r1, c1, ctx);
        TEMPLATE(T, mat_window_init) (window2, C, 0, c1, r1, (c1 + c2), ctx);


        if (!(TEMPLATE(T, mat_equal) (window1, A, ctx) && TEMPLATE(T, mat_equal) (window2, B, ctx)))
        {
            flint_printf("FAIL: results not equal\n");
            abort();
        }

        
        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);

        TEMPLATE(T, mat_window_clear) (window1, ctx);
        TEMPLATE(T, mat_window_clear) (window2, ctx);

    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


#endif
