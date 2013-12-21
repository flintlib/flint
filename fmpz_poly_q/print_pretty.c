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

    Copyright (C) 2010, 2011 Sebastian Pancratz
   
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "fmpz_poly_q.h"

int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x)
{
    char *str;

    str = fmpz_poly_q_get_str_pretty(op, x);
    flint_printf("%s", str);
    flint_free(str);

    return 1;
}
