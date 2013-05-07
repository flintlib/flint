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

    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int result;
    char *str;
    fmpz_poly_t a;

    printf("get_str_pretty....");
    fflush(stdout);

    fmpz_poly_init(a);

    str = fmpz_poly_get_str_pretty(a, "t");
    result = strcmp(str, "0") == 0;
    if (!result)
    {
        printf("FAIL:\n");
        printf("a = "), fmpz_poly_print(a), printf("\n");
        printf("str(a) = {%s}\n", str);
        abort();
    }
    flint_free(str);

    fmpz_poly_set_si(a, -2);
    str = fmpz_poly_get_str_pretty(a, "t");
    result = strcmp(str, "-2") == 0;
    if (!result)
    {
        printf("FAIL:\n");
        printf("a = "), fmpz_poly_print(a), printf("\n");
        printf("str(a) = {%s}\n", str);
        abort();
    }
    flint_free(str);

    fmpz_poly_set_coeff_si(a, 3, 1);
    str = fmpz_poly_get_str_pretty(a, "t");
    result = strcmp(str, "t^3-2") == 0;
    if (!result)
    {
        printf("FAIL:\n");
        printf("a = "), fmpz_poly_print(a), printf("\n");
        printf("str(a) = {%s}\n", str);
        abort();
    }
    flint_free(str);
    fmpz_poly_clear(a);

    printf("PASS\n");
    return 0;
}
