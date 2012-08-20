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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

/*
    Demo FLINT program to demonstrate some use of the qadic module.
*/

#include <stdlib.h>
#include <stdio.h>

#include "qadic.h"

int main(void)
{
    fmpz_t p;
    qadic_ctx_t ctx;
    qadic_t a, b, c;
    fmpz_t e = {4L};
    fmpz nine[1] = {9L};

    fmpz_init_set_ui(p, 3);
    qadic_ctx_init_conway(ctx, p, 2, 5, "a", PADIC_SERIES);

    qadic_init(a);
    qadic_init(b);
    qadic_init(c);

    printf("Compute a power and a sum\n");
    padic_poly_fit_length(a, 2);
    fmpz_one(a->coeffs + 0);
    fmpz_set_ui(a->coeffs + 1, 2);
    a->val = 0;
    _padic_poly_set_length(a, 2);

    qadic_print_pretty(a, ctx); printf("\n");

    qadic_pow(a, a, e, ctx);
    padic_poly_set_ui(b, 3249, &ctx->pctx);
    qadic_add(c, a, b, ctx);

    qadic_print_pretty(a, ctx); printf("\n");
    qadic_print_pretty(b, ctx); printf("\n");
    qadic_print_pretty(c, ctx); printf("\n");

    printf("Compute a Teichmuller lift\n");
    padic_poly_fit_length(a, 2);
    fmpz_one(a->coeffs + 0);
    fmpz_set_ui(a->coeffs + 1, 2);
    a->val = 0;
    _padic_poly_set_length(a, 2);

    qadic_teichmuller(b, a, ctx);
    qadic_pow(c, b, nine, ctx);

    qadic_print_pretty(a, ctx); printf("\n");
    qadic_print_pretty(b, ctx); printf("\n");
    qadic_print_pretty(c, ctx); printf("\n");

    printf("Compute an inverse\n");
    qadic_set(a, b);
    qadic_inv(b, a, ctx);
    qadic_mul(c, a, b, ctx);

    qadic_print_pretty(a, ctx); printf("\n");
    qadic_print_pretty(b, ctx); printf("\n");
    qadic_print_pretty(c, ctx); printf("\n");

    qadic_clear(a);
    qadic_clear(b);
    qadic_clear(c);

    qadic_ctx_clear(ctx);
    fmpz_clear(p);

    /*************************************************************************/

    printf("Compute a Frobenius image\n");

    fmpz_init_set_ui(p, 3);
    qadic_ctx_init_conway(ctx, p, 2, 5, "X", PADIC_TERSE);

    qadic_init(a);
    qadic_init(b);

    padic_poly_fit_length(a, 2);
    a->coeffs[0] = 78L;
    a->coeffs[1] = 45L;
    a->val = 0;
    _padic_poly_set_length(a, 2);

    qadic_frobenius(b, a, 1, ctx);

    printf("Context:\n"), qadic_ctx_print(ctx), printf("\n");

    printf("a = "), qadic_print_pretty(a, ctx), printf("\n");
    printf("b = "), qadic_print_pretty(b, ctx), printf("\n");

    qadic_clear(a);
    qadic_clear(b);

    qadic_ctx_clear(ctx);
    fmpz_clear(p);

    return EXIT_SUCCESS;
}

