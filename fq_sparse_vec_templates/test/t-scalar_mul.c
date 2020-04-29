/*
    wopyright (w) 2011 Fredrik Johansson

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
    slong rep, len, nnz;
    TEMPLATE(T, t) c;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_vec_t) u, v, w, x;
    FLINT_TEST_INIT(state);
    
    flint_printf("scalar_mul and muladd....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, init) (c, ctx);
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        TEMPLATE(T, sparse_vec_init) (u, ctx);
        TEMPLATE(T, sparse_vec_init) (v, ctx);
        TEMPLATE(T, sparse_vec_init) (w, ctx);
        TEMPLATE(T, sparse_vec_init) (x, ctx);

        TEMPLATE(T, sparse_vec_randtest) (u, state, nnz, len, ctx);
        TEMPLATE(T, sparse_vec_randtest) (v, state, nnz, len, ctx);
        TEMPLATE(T, randtest) (c, state, ctx);

        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (w, u, v, c, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (x, v, c, ctx);
        TEMPLATE(T, sparse_vec_add) (x, u, x, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (w, x, 0, ctx))
        {
            flint_printf("FAIL: u + c*v != u + (c*v)\n");
            abort();
        }

        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T)) (w, u, v, c, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (x, v, c, ctx);
        TEMPLATE(T, sparse_vec_sub) (x, u, x, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (w, x, 0, ctx))
        {
            flint_printf("FAIL: u + c*v != u + (c*v)\n");
            abort();
        }

        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (w, u, v, c, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (u, u, v, c, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (u, w, 0, ctx))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (w, u, v, c, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (v, u, v, c, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (v, w, 0, ctx))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T)) (w, u, v, c, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T)) (u, u, v, c, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (u, w, 0, ctx))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T)) (w, u, v, c, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T)) (v, u, v, c, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (v, w, 0, ctx))
        {
            flint_printf("FAIL: u + c*v != (u += c*v)\n");
            abort();
        }

        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (x, v, c, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (v, v, c, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (v, x, 0, ctx))
        {
            flint_printf("FAIL: c*v != (c *= v)\n");
            abort();
        }

        TEMPLATE(T, sparse_vec_clear) (u, ctx);
        TEMPLATE(T, sparse_vec_clear) (v, ctx);
        TEMPLATE(T, sparse_vec_clear) (w, ctx);
        TEMPLATE(T, sparse_vec_clear) (x, ctx);
        TEMPLATE(T, clear) (c, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif