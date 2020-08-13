/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program to demonstrate some use of the padic module.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "flint.h"
#include "padic.h"

int main(void)
{
    fmpz_t p;
    padic_ctx_t ctx;

    char *str;
    padic_t x, y;

    flint_printf("Output:\n\n");

    /* Case 1 */
    flint_printf("Positive integer:  x = 127 mod 7^10\n");

    fmpz_init(p);
    fmpz_set_ui(p, 7);
    padic_ctx_init(ctx, p, 8, 12, PADIC_TERSE);

    padic_init2(x, 10);
    padic_set_ui(x, 127, ctx);

    ctx->mode = PADIC_TERSE;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_SERIES;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_VAL_UNIT;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    padic_clear(x);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    /* Case 2 */
    flint_printf("Positive integer larger than p^N:  x = 1057 mod 2^10\n");

    fmpz_init(p);
    fmpz_set_ui(p, 2);
    padic_ctx_init(ctx, p, 10, 12, PADIC_TERSE);

    padic_init2(x, 10);
    padic_set_ui(x, 1057, ctx);

    ctx->mode = PADIC_TERSE;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_SERIES;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_VAL_UNIT;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    padic_clear(x);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    /* Case 3 */
    flint_printf("Negative integer:  x = -127 mod 3^10\n");

    fmpz_init(p);
    fmpz_set_ui(p, 3);
    padic_ctx_init(ctx, p, 10, 12, PADIC_TERSE);

    padic_init2(x, 10);
    padic_set_si(x, -127, ctx);

    ctx->mode = PADIC_TERSE;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_VAL_UNIT;
    str = padic_get_str(NULL, x, ctx);
    flint_printf("print:   "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("get_str: %s\n", str);
    flint_free(str);

    padic_clear(x);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    /* Log */
    flint_printf("Log of 7380996 mod 5^20\n");

    fmpz_init(p);
    fmpz_set_ui(p, 5);
    padic_ctx_init(ctx, p, 10, 25, PADIC_SERIES);

    padic_init(x);
    padic_init(y);
    padic_set_ui(x, 7380996, ctx);

    padic_log(y, x, ctx);

    flint_printf("x = "), padic_print(x, ctx), flint_printf("\n");
    flint_printf("y = "), padic_print(y, ctx), flint_printf("\n");

    padic_clear(x);
    padic_clear(y);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    return EXIT_SUCCESS;
}

