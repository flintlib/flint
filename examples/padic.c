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
    Demo FLINT program to demonstrate some use of the padic module.
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>

#include "flint.h"
#include "padic.h"

int main(void)
{
    fmpz_t p;
    padic_ctx_t ctx;

    char *str;
    padic_t x, y;

    printf("Output:\n\n");

    /* Case 1 */
    printf("Positive integer:  x = 127 mod 7^10\n");

    fmpz_init(p);
    fmpz_set_ui(p, 7);
    padic_ctx_init(ctx, p, 10, PADIC_TERSE);

    _padic_init(x);
    _padic_set_ui(x, 127, ctx);

    ctx->mode = PADIC_TERSE;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_SERIES;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_VAL_UNIT;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    _padic_clear(x);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    /* Case 2 */
    printf("Positive integer larger than p^N:  x = 1057 mod 2^10\n");

    fmpz_init(p);
    fmpz_set_ui(p, 2);
    padic_ctx_init(ctx, p, 10, PADIC_TERSE);

    _padic_init(x);
    _padic_set_ui(x, 1057, ctx);

    ctx->mode = PADIC_TERSE;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_SERIES;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_VAL_UNIT;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    _padic_clear(x);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    /* Case 3 */
    printf("Negative integer:  x = -127 mod 3^10\n");

    fmpz_init(p);
    fmpz_set_ui(p, 3);
    padic_ctx_init(ctx, p, 10, PADIC_TERSE);

    _padic_init(x);
    _padic_set_si(x, -127, ctx);

    ctx->mode = PADIC_TERSE;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    ctx->mode = PADIC_VAL_UNIT;
    str = _padic_get_str(NULL, x, ctx);
    printf("print:   "), padic_print(x, ctx), printf("\n");
    printf("get_str: %s\n", str);
    flint_free(str);

    _padic_clear(x);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    /* Log */
    printf("Log of 7380996 mod 5^10\n");

    fmpz_init(p);
    fmpz_set_ui(p, 5);
    padic_ctx_init(ctx, p, 10, PADIC_SERIES);

    _padic_init(x);
    _padic_init(y);
    _padic_set_ui(x, 7380996, ctx);

    padic_log(y, x, ctx);

    printf("x = "), padic_print(x, ctx), printf("\n");
    printf("y = "), padic_print(y, ctx), printf("\n");

    _padic_clear(x);
    _padic_clear(y);

    padic_ctx_clear(ctx);
    fmpz_clear(p);

    return EXIT_SUCCESS;
}

