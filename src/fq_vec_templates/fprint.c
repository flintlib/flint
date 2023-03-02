/*
    Copyright (C) 2008, 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdlib.h>
#include <stdio.h>

/*
    Recall the return value conventions for fputc (of type int) 

    ``If there are no errors, the same character that has been written is 
    returned.  If an error occurs, EOF is returned and the error indicator 
    is set''

    where the EOF macro expands to a negative int, and fprintf (of type int)

    ``On success, the total number of characters written is returned.
    On failure, a negative number is returned.''
 */

int _TEMPLATE(T, vec_fprint) (FILE * file,
                              const TEMPLATE(T, struct) * vec,
                              slong len, const TEMPLATE(T, ctx_t) ctx)
{
    int r;
    slong i;

    r = fprintf(file, "%li", len);
    if ((len > 0) && (r > 0))
    {
        r = fputc(' ', file);
        for (i = 0; (i < len) && (r > 0); i++)
        {
            r = fputc(' ', file);
            if (r > 0)
                r = TEMPLATE(T, fprint) (file, vec + i, ctx);
        }
    }

    return r;
}


#endif
