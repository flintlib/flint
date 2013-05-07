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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"

void
fmpz_factor_print(const fmpz_factor_t factor)
{
    long i;

    if (factor->sign == 0)
    {
        printf("0");
        return;
    }

    if (factor->sign == -1)
    {
        if (factor->num)
            printf("-1 * ");
        else
            printf("-1");
    }

    for (i = 0; i < factor->num; i++)
    {
        fmpz_print(factor->p + i);

        if (factor->exp[i] != 1UL)
            printf("^%lu", factor->exp[i]);

        if (i != factor->num - 1)
            printf(" * ");
    }
}
