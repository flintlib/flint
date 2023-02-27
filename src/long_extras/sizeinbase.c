/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <limits.h>
#include "flint.h"
#include "long_extras.h"

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

