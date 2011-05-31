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

#include <limits.h>

#include "padic.h"

/*
    Returns the number of digits in the base $b$ representation 
    of the integer $n$.

    Assumes that $b \geq 2$.
 */
static long sizeinbase_si(long n, long b)
{
    long c;

    if (n > 0)
        c = 0;
    else if (n > LONG_MIN)
    {
        n = -n;
        c = 1;
    }
    else  /* n == LONG_MIN */
    {
        if (n % b)
        {
            n = - (n + 1);
            c = 1;
        }
        else
        {
            n = - (n / b);
            c = 2;
        }
    }

    for ( ; n > 0; n /= b, c++) ;

    return c;
}

char * padic_get_str(const padic_t op, const padic_ctx_t ctx)
{
    if (1)
    {
        printf("ERROR (padic_get_str).  Not implemented yet.\n");
        abort();
    }

    return NULL;
}

