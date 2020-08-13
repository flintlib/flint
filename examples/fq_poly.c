/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program to demonstrate some use of the fq_poly module.
*/

#include <time.h>

#include "fq_poly.h"

int main(void)
{
    fmpz_t p;
    slong d, i;
    fq_ctx_t ctx;
    clock_t c0, c1;
    double c;
    fq_poly_t f, g, h;

    FLINT_TEST_INIT(state);

    fq_poly_init(f, ctx);
    fq_poly_init(g, ctx);
    fq_poly_init(h, ctx);

    printf("Polynomial multiplication over GF(q)\n");
    printf("------------------------------------\n");

    {
        printf("1)  Two length-10,000 polynomials over GF(3^2)\n");

        fmpz_init_set_ui(p, 3);
        d = 2;
        fq_ctx_init_conway(ctx, p, d, "X");

        fq_poly_randtest(g, state, 10000, ctx);
        fq_poly_randtest(h, state, 10000, ctx);

        c0 = clock();
        fq_poly_mul_classical(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("Classical: %fs\n", c);

        c0 = clock();
        for (i = 0; i < 100; i++)
            fq_poly_mul_reorder(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("Reorder: %fms\n", 10 * c);

        c0 = clock();
        for (i = 0; i < 100; i++)
            fq_poly_mul_KS(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("KS: %fms\n", 10 * c);

        fq_ctx_clear(ctx);
        fmpz_clear(p);
    }
    {
        printf("2)  Two length-500 polynomials over GF(3^263)\n");

        fmpz_init_set_ui(p, 3);
        d = 263;
        fq_ctx_init_conway(ctx, p, d, "X");

        fq_poly_randtest(g, state, 500, ctx);
        fq_poly_randtest(h, state, 500, ctx);

        c0 = clock();
        fq_poly_mul_classical(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("Classical: %fs\n", c);

        c0 = clock();
        fq_poly_mul_reorder(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("Reorder: %fs\n", c);

        c0 = clock();
        for (i = 0; i < 100; i++)
            fq_poly_mul_KS(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("KS: %fms\n", 10 * c);

        fq_ctx_clear(ctx);
        fmpz_clear(p);
    }
    {
        printf("3)  Two length-5 polynomials over GF(109987^4)\n");

        fmpz_init_set_ui(p, 109987);
        d = 4;
        fq_ctx_init_conway(ctx, p, d, "X");

        fq_poly_randtest(g, state, 4, ctx);
        fq_poly_randtest(h, state, 4, ctx);

        c0 = clock();
        for (i = 0; i < 1000 * 100; i++)
            fq_poly_mul_classical(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("Classical: %f\xb5s\n", 10 * c);

        c0 = clock();
        for (i = 0; i < 1000 * 100; i++)
            fq_poly_mul_reorder(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("Reorder: %f\xb5s\n", 10 * c);

        c0 = clock();
        for (i = 0; i < 1000 * 100; i++)
            fq_poly_mul_KS(f, g, h, ctx);
        c1 = clock();
        c  = (double) (c1 - c0) / CLOCKS_PER_SEC;
        printf("KS: %f\xb5s\n", 10 * c);

        fq_ctx_clear(ctx);
        fmpz_clear(p);
    }

    fq_poly_clear(f, ctx);
    fq_poly_clear(g, ctx);
    fq_poly_clear(h, ctx);
    FLINT_TEST_CLEANUP(state);

    return EXIT_SUCCESS;
}

