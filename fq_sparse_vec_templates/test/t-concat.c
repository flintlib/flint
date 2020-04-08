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
#include <gmp.h>
#include "ulong_extras.h"

int main(void)
{
    slong rep, len, nnz;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_vec_t) u, v, w;
    TEMPLATE(T, sparse_vec_t) window1, window2;
    FLINT_TEST_INIT(state);

    flint_printf("concat....");
    fflush(stdout);

    for (rep = 0; rep < 100; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);

        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);

        TEMPLATE(T, sparse_vec_init) (u, ctx);
        TEMPLATE(T, sparse_vec_init) (v, ctx);
        TEMPLATE(T, sparse_vec_init) (w, ctx);

        TEMPLATE(T, sparse_vec_randtest) (u, state, nnz, len, ctx);
        TEMPLATE(T, sparse_vec_randtest) (v, state, nnz, len, ctx);
        TEMPLATE(T, sparse_vec_randtest) (w, state, nnz, len, ctx);

        TEMPLATE(T, sparse_vec_concat) (w, u, v, len, ctx);        
        
        TEMPLATE(T, sparse_vec_window_init) (window1, w, 0, len, ctx);
        TEMPLATE(T, sparse_vec_window_init)(window2, w, len, 2*len, ctx);

        if (!(TEMPLATE(T, sparse_vec_equal) (window1, u, 0, ctx) && TEMPLATE(T, sparse_vec_equal) (window2, v, len, ctx)))
        {
            flint_printf("u = ");
            TEMPLATE(T, sparse_vec_print_pretty) (u, 0, len, ctx);
            flint_printf("v = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (v, 0, len, ctx);
            flint_printf("u | v = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (w, 0, len, ctx);
            flint_printf("window1 = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (window1, 0, len, ctx);
            flint_printf("window2 = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (window2, len, len, ctx);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
        TEMPLATE(T, sparse_vec_window_clear) (window1, ctx);
        TEMPLATE(T, sparse_vec_window_clear) (window2, ctx);

        TEMPLATE(T, sparse_vec_init) (window1, ctx);
        TEMPLATE(T, sparse_vec_init) (window2, ctx);
        TEMPLATE(T, sparse_vec_split) (window1, window2, w, len, ctx);
        if (!(TEMPLATE(T, sparse_vec_equal) (window1, u, 0, ctx) && TEMPLATE(T, sparse_vec_equal) (window2, v, 0, ctx)))
        {
            flint_printf("u = ");
            TEMPLATE(T, sparse_vec_print_pretty) (u, 0, len, ctx);
            flint_printf("v = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (v, 0, len, ctx);
            flint_printf("u | v = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (w, 0, len, ctx);
            flint_printf("window1 = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (window1, 0, len, ctx);
            flint_printf("window2 = \n");
            TEMPLATE(T, sparse_vec_print_pretty) (window2, 0, len, ctx);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
        TEMPLATE(T, sparse_vec_window_clear) (window1, ctx);
        TEMPLATE(T, sparse_vec_window_clear) (window2, ctx);
        TEMPLATE(T, sparse_vec_clear) (u, ctx);
        TEMPLATE(T, sparse_vec_clear) (v, ctx);
        TEMPLATE(T, sparse_vec_clear) (w, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif