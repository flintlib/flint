/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

char *
fq_zech_get_str(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    char *s = flint_malloc(n_clog(op->value, 10) * sizeof(char));
    flint_sprintf(s, "%wd", op->value);
    return s;
}
