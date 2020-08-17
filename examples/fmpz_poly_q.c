/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Simple example demonstrating the use of the fmpz_poly_q module.
 */

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_q.h"

int main(int argc, char* argv[])
{
    char *str, *strf, *strg;
    fmpz_poly_q_t f, g;

    fmpz_poly_q_init(f);
    fmpz_poly_q_init(g);
    fmpz_poly_q_set_str(f, "2  1 3/1  2");
    fmpz_poly_q_set_str(g, "1  3/2  2 7");
    strf = fmpz_poly_q_get_str_pretty(f, "t");
    strg = fmpz_poly_q_get_str_pretty(g, "t");
    fmpz_poly_q_mul(f, f, g);
    str  = fmpz_poly_q_get_str_pretty(f, "t");
    flint_printf("%s * %s = %s\n", strf, strg, str);
    flint_free(str);
    flint_free(strf);
    flint_free(strg);
    fmpz_poly_q_clear(f);
    fmpz_poly_q_clear(g);

    return EXIT_SUCCESS;
}

