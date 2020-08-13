/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_nth_prime(ulong n)
{
    if (n == 0)
    {
        flint_printf("Exception (n_nth_prime). n_nth_prime(0) is undefined.\n");
        flint_abort();
    }

    return n_primes_arr_readonly(n)[n-1];
}

