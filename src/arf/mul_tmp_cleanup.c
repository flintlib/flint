/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

FLINT_TLS_PREFIX mp_ptr __arf_mul_tmp = NULL;
FLINT_TLS_PREFIX slong __arf_mul_alloc = 0;

void _arf_mul_tmp_cleanup(void)
{
    flint_free(__arf_mul_tmp);
    __arf_mul_tmp = NULL;
    __arf_mul_alloc = 0;
}
