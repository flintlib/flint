/*
    Copyright (C) 2019 Julian Rüth

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>

#include "arf.h"

int
arf_dump_file(FILE* stream, const arf_t x)
{
    int nwrite;
    char* data = arf_dump_str(x);
    nwrite = fputs(data, stream);
    if (nwrite == EOF)
        return nwrite;

    flint_free(data);
    return 0;
}
