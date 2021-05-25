/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Example program for the fmpz_mod_poly module.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz_mod_poly.h"

int main(int argc, char* argv[])
{
    fmpz_t n;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t x, y;

    fmpz_init_set_ui(n, 7);
    fmpz_mod_ctx_init(ctx, n);
    fmpz_mod_poly_init(x, ctx);
    fmpz_mod_poly_init(y, ctx);
    fmpz_mod_poly_set_coeff_ui(x, 3, 5, ctx);
    fmpz_mod_poly_set_coeff_ui(x, 0, 6, ctx);
    fmpz_mod_poly_sqr(y, x, ctx);
    fmpz_mod_poly_print(x, ctx); flint_printf("\n");
    fmpz_mod_poly_print(y, ctx); flint_printf("\n");
    fmpz_mod_poly_clear(x, ctx);
    fmpz_mod_poly_clear(y, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(n);

    return EXIT_SUCCESS;
}

