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
#include <sys/time.h>

int
main(void)
{
    slong rep, r, c;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, Ai;
    TEMPLATE(T, mat_t) dA, dAiA;
    FLINT_TEST_INIT(state);
    
    flint_printf("inverting A....");
    fflush(stdout);
    
    for (rep = 0; rep < 200; rep++)
    {
        if (rep % 5==0) {flint_printf("."); fflush(stdout);}
        TEMPLATE(T, ctx_randtest) (ctx, state);
        do r = n_randint(state, 100), c = n_randint(state, 100);
        while (r==UWORD(0) || c==UWORD(0));
        
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (Ai, r, r, ctx);
        TEMPLATE(T, sparse_mat_randtest) (A, state, 0, c, ctx);
        TEMPLATE(T, mat_init) (dA, r, c, ctx);
        TEMPLATE(T, mat_init) (dAiA, r, c, ctx);
        TEMPLATE(T, sparse_mat_to_dense) (dA, A, ctx);

        TEMPLATE(T, sparse_mat_inv) (Ai, A, ctx);
        TEMPLATE(T, sparse_mat_mul_mat) (dAiA, Ai, dA, ctx);
        TEMPLATE(T, mat_rref) (dA, ctx);
        if (!TEMPLATE(T, mat_equal) (dAiA, dA, ctx)) 
        {
            flint_printf("FAIL!\n");
            flint_printf("A^-1 x A = ");
            TEMPLATE(T, mat_print_pretty) (dAiA, ctx);
            flint_printf("rref(A) = ");
            TEMPLATE(T, mat_print_pretty) (dA, ctx);
            abort();
        }

        TEMPLATE(T, sparse_mat_clear) (Ai, ctx);
        TEMPLATE(T, mat_clear) (dA, ctx);
        TEMPLATE(T, mat_clear) (dAiA, ctx);
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
