/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define PADIC_INLINES_C

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

#include "padic.h"

#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

void padic_get_unit(fmpz_t f, padic_t p)
{
   fmpz_set(f, padic_unit(p));
}
