/*
    Copyright (C) 2010 Fredrik Johansson

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

int
main(void)
{
    slong rep, r, c, k;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A;
    TEMPLATE(T, struct) *x, *y, *y2;
    TEMPLATE(T, mat_t) B, X, Y, Y2;
    FLINT_TEST_INIT(state);
    

    flint_printf("multiplication by vec/mat....");
    fflush(stdout);

    for (rep = 0; rep < 10; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        r = n_randint(state, 200);
        c = n_randint(state, 200);
        k = n_randint(state, 200);
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, mat_init) (B, r, c, ctx);
        TEMPLATE(T, mat_init) (X, c, k, ctx);
        TEMPLATE(T, mat_init) (Y, r, k, ctx);
        TEMPLATE(T, mat_init) (Y2, r, k, ctx);
        x = _TEMPLATE(T, vec_init) (c, ctx);
        y = _TEMPLATE(T, vec_init) (r, ctx);
        y2 = _TEMPLATE(T, vec_init) (r, ctx);

        TEMPLATE(T, sparse_mat_randtest) (A, state, 0, c, ctx);
        TEMPLATE(T, sparse_mat_to_dense) (B, A, ctx);
        _TEMPLATE(T, vec_randtest) (x, state, c, ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (y, A, x, ctx);
        TEMPLATE(T, mat_mul_vec) (y2, B, x, ctx);

        if (!_TEMPLATE(T, vec_equal) (y, y2, A->r, ctx))
        {
            flint_printf("FAIL: y != y2\n");
            abort();
        }

        TEMPLATE(T, mat_randtest) (X, state, ctx);
        TEMPLATE(T, sparse_mat_mul_mat) (Y, A, X, ctx);
        TEMPLATE(T, mat_mul) (Y2, B, X, ctx);

        if (!TEMPLATE(T, mat_equal) (Y, Y2, ctx))
        {
            flint_printf("Fail: Y != Y2\n");
            abort();
        }

        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        _TEMPLATE(T, vec_clear) (x, c, ctx);
        _TEMPLATE(T, vec_clear) (y, r, ctx);
        _TEMPLATE(T, vec_clear) (y2, r, ctx);
        TEMPLATE(T, mat_clear) (X, ctx);
        TEMPLATE(T, mat_clear) (Y, ctx);
        TEMPLATE(T, mat_clear) (Y2, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
     }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
