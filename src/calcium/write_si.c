/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "calcium.h"

void
calcium_write_si(calcium_stream_t out, slong x)
{
    if (out->fp != NULL)
    {
        flint_fprintf((FILE *) out->fp, "%wd", x);
    }
    else
    {
        char tmp[22];
        sprintf(tmp, WORD_FMT "d", x);
        calcium_write(out, tmp);
    }
}
