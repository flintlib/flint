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

#include <stdlib.h>
#include <limits.h>
#include "flint.h"

size_t z_sizeinbase(slong n, int b)
{
    slong c = 0;

    if (n == 0)
    {
        return 1;
    }

    if (n <= 0)
    {
        if (n > WORD_MIN)
        {
            n = -n;
        }
        else  /* n == WORD_MIN */
        {
            if (n % b)
            {
                n = - (n + 1);
            }
            else
            {
                n = - (n / b);
                c = 1;
            }
        }
    }

    for ( ; n > 0; n /= b, c++) ;

    return c;
}

