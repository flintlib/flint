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

    padic_t x;

    fmpz_init(p);
    fmpz_set_ui(p, 7);
    padic_ctx_init(ctx, p, 10, PADIC_TERSE);

    _padic_init(x);
    _padic_set_ui(x, 127, ctx);

    ctx->mode = PADIC_TERSE;
    _padic_print(x, ctx), printf("\n");

    ctx->mode = PADIC_SERIES;
    _padic_print(x, ctx), printf("\n");

    ctx->mode = PADIC_VAL_UNIT;
    _padic_print(x, ctx), printf("\n");

    _padic_clear(x);
    padic_ctx_clear(ctx);
    fmpz_clear(p);

    return EXIT_SUCCESS;
}

