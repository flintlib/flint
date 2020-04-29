/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include "ulong_extras.h"

int main(void)
{
    slong rep, r1, r2, c, nreps = 100;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, B, C;
    TEMPLATE(T, sparse_mat_t) window1, window2;
    FLINT_TEST_INIT(state);


    flint_printf("concat_vertical....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        r1 = n_randint(state, 100);
        r2 = n_randint(state, 100);
        c = n_randint(state, 100);
        TEMPLATE(T, sparse_mat_init) (A, r1, c, ctx);
        TEMPLATE(T, sparse_mat_init) (B, r2, c, ctx);
        TEMPLATE(T, sparse_mat_init) (C, r1+r2, c, ctx);

        TEMPLATE(T, sparse_mat_randtest) (A, state, 0, c, ctx);
        TEMPLATE(T, sparse_mat_randtest) (B, state, 0, c, ctx);
        TEMPLATE(T, sparse_mat_randtest) (C, state, 0, c, ctx);

        TEMPLATE(T, sparse_mat_concat_vertical) (C, A, B, ctx);
        
        TEMPLATE(T, sparse_mat_window_init) (window1, C, 0, 0, r1, c, ctx);
        TEMPLATE(T, sparse_mat_window_init) (window2, C, r1, 0, r1+r2, c, ctx);

        if (!(TEMPLATE(T, sparse_mat_equal) (window1, A, ctx) && TEMPLATE(T, sparse_mat_equal) (window2, B, ctx)))
        {
            flint_printf("A = \n");
            TEMPLATE(T, sparse_mat_print_pretty) (A, ctx);
            flint_printf("B = \n");
            TEMPLATE(T, sparse_mat_print_pretty) (B, ctx);
            flint_printf("A concat_vertical B = \n");
            TEMPLATE(T, sparse_mat_print_pretty) (C, ctx);
            flint_printf("FAIL: results not equal\n");
            abort();
        }

        TEMPLATE(T, sparse_mat_window_clear) (window1, ctx);
        TEMPLATE(T, sparse_mat_window_clear) (window2, ctx);

        TEMPLATE(T, sparse_mat_init) (window1, r1, c, ctx);
        TEMPLATE(T, sparse_mat_init) (window2, r2, c, ctx);
        TEMPLATE(T, sparse_mat_split_vertical) (window1, window2, C, r1, ctx);
        
        if (!TEMPLATE(T, sparse_mat_equal) (window1, A, ctx) || !TEMPLATE(T, sparse_mat_equal) (window2, B, ctx))
        {
            flint_printf("A = \n");
            TEMPLATE(T, sparse_mat_print_pretty) (A, ctx);
            flint_printf("B = \n");
            TEMPLATE(T, sparse_mat_print_pretty) (B, ctx);
            flint_printf("A concat_vertical B = \n");
            TEMPLATE(T, sparse_mat_print_pretty) (C, ctx);
            flint_printf("FAIL: results not equal\n");
            abort();
        }

        TEMPLATE(T, sparse_mat_window_clear) (window1, ctx);
        TEMPLATE(T, sparse_mat_window_clear) (window2, ctx);
        
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, sparse_mat_clear) (B, ctx);
        TEMPLATE(T, sparse_mat_clear) (C, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
