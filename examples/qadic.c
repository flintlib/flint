/*
    Copyright (C) 2011, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program to demonstrate some use of the qadic module.
*/

#include <stdlib.h>
#include <stdio.h>

#include "qadic.h"

int main(void)
{
    fmpz_t p;
    slong d, N;
    qadic_ctx_t ctx;

    qadic_t a, b, c;

    int ans;

    /*************************************************************************/

    fmpz_t e = {WORD(4)};
    fmpz_t nine = {WORD(9)};

    fmpz_init_set_ui(p, 3);
    d = 2;
    N = 5;
    qadic_ctx_init_conway(ctx, p, d, 0, N, "a", PADIC_SERIES);

    qadic_init2(a, N);
    qadic_init2(b, N);
    qadic_init2(c, N);

    flint_printf("Compute a power and a sum\n");
    padic_poly_fit_length(a, 2);
    fmpz_one(a->coeffs + 0);
    fmpz_set_ui(a->coeffs + 1, 2);
    a->val = 0;
    _padic_poly_set_length(a, 2);

    qadic_print_pretty(a, ctx); flint_printf("\n");

    qadic_pow(a, a, e, ctx);
    padic_poly_set_ui(b, 3249, &ctx->pctx);
    qadic_add(c, a, b, ctx);

    qadic_print_pretty(a, ctx); flint_printf("\n");
    qadic_print_pretty(b, ctx); flint_printf("\n");
    qadic_print_pretty(c, ctx); flint_printf("\n");
    flint_printf("\n");

    flint_printf("Compute a Teichmuller lift\n");
    padic_poly_fit_length(a, 2);
    fmpz_one(a->coeffs + 0);
    fmpz_set_ui(a->coeffs + 1, 2);
    a->val = 0;
    _padic_poly_set_length(a, 2);

    qadic_teichmuller(b, a, ctx);
    qadic_pow(c, b, nine, ctx);

    qadic_print_pretty(a, ctx); flint_printf("\n");
    qadic_print_pretty(b, ctx); flint_printf("\n");
    qadic_print_pretty(c, ctx); flint_printf("\n");
    flint_printf("\n");

    flint_printf("Compute an inverse\n");
    qadic_set(a, b, ctx);
    qadic_inv(b, a, ctx);
    qadic_mul(c, a, b, ctx);

    qadic_print_pretty(a, ctx); flint_printf("\n");
    qadic_print_pretty(b, ctx); flint_printf("\n");
    qadic_print_pretty(c, ctx); flint_printf("\n");
    flint_printf("\n");

    qadic_clear(a);
    qadic_clear(b);
    qadic_clear(c);

    qadic_ctx_clear(ctx);
    fmpz_clear(p);

    /*************************************************************************/

    flint_printf("Compute a Frobenius image\n");

    fmpz_init_set_ui(p, 3);
    d = 2;
    N = 5;
    qadic_ctx_init_conway(ctx, p, d, 0, N, "X", PADIC_TERSE);

    qadic_init2(a, N);
    qadic_init2(b, N);

    padic_poly_fit_length(a, 2);
    a->coeffs[0] = WORD(78);
    a->coeffs[1] = WORD(45);
    a->val = 0;
    _padic_poly_set_length(a, 2);

    qadic_frobenius(b, a, 1, ctx);
    flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
    flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
    flint_printf("Context:\n"), qadic_ctx_print(ctx);
    flint_printf("\n");

    qadic_clear(a);
    qadic_clear(b);

    qadic_ctx_clear(ctx);
    fmpz_clear(p);

    /*************************************************************************/

    flint_printf("Compute a square root\n");

    fmpz_init_set_ui(p, 2);
    d = 3;
    N = 2;
    qadic_ctx_init_conway(ctx, p, d, 0, N, "X", PADIC_SERIES);

    qadic_init2(a, N);
    qadic_init2(b, N);

    padic_poly_fit_length(a, d);
    a->coeffs[0] = WORD(1);
    a->coeffs[1] = WORD(3);
    a->coeffs[2] = WORD(1);
    a->val = 0;
    _padic_poly_set_length(a, d);

    ans = qadic_sqrt(b, a, ctx);
    flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
    flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
    flint_printf("ans = %d\n", ans);
    flint_printf("Context:\n"), qadic_ctx_print(ctx);
    flint_printf("\n");

    qadic_clear(a);
    qadic_clear(b);

    qadic_ctx_clear(ctx);
    fmpz_clear(p);


    return EXIT_SUCCESS;
}

