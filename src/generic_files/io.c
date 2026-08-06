/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdarg.h>
#include <stdio.h>
#include "flint.h"

int flint_printf(const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vfprintf(stdout, str, vlist);
   va_end(vlist);

   return ret;
}

int flint_vprintf(const char * str, va_list vlist)
{
    return flint_vfprintf(stdout, str, vlist);
}
