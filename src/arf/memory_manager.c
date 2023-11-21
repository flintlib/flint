/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

#define ARF_USE_CACHE 0

#define ARF_MAX_CACHE_LIMBS 64

FLINT_TLS_PREFIX mp_ptr * arf_free_arr = NULL;
FLINT_TLS_PREFIX ulong arf_free_num = 0;
FLINT_TLS_PREFIX ulong arf_free_alloc = 0;
FLINT_TLS_PREFIX int arf_have_registered_cleanup = 0;

void _arf_cleanup(void)
{
    slong i;
    for (i = 0; i < arf_free_num; i++)
        flint_free(arf_free_arr[i]);

    flint_free(arf_free_arr);

    arf_free_arr = NULL;
    arf_free_num = 0;
    arf_free_alloc = 0;
}

void
_arf_promote(arf_t x, mp_size_t n)
{
    if (ARF_USE_CACHE && n <= ARF_MAX_CACHE_LIMBS && arf_free_num != 0)
    {
        mp_ptr ptr;
        mp_size_t alloc;

        ptr = arf_free_arr[--arf_free_num];
        alloc = ptr[0];

        if (alloc >= n)
        {
            ARF_PTR_ALLOC(x) = ptr[0];
            ARF_PTR_D(x) = ptr;
        }
        else
        {
            ptr = flint_realloc(ptr, n * sizeof(mp_limb_t));
            ARF_PTR_ALLOC(x) = n;
            ARF_PTR_D(x) = ptr;
        }
    }
    else
    {
        ARF_PTR_ALLOC(x) = n;
        ARF_PTR_D(x) = flint_malloc(n * sizeof(mp_limb_t));
    }
}

void
_arf_demote(arf_t x)
{
    mp_ptr ptr;
    mp_size_t alloc;

    alloc = ARF_PTR_ALLOC(x);
    ptr = ARF_PTR_D(x);

    if (ARF_USE_CACHE && alloc <= ARF_MAX_CACHE_LIMBS)
    {
        if (arf_free_num == arf_free_alloc)
        {
            if (!arf_have_registered_cleanup)
            {
                flint_register_cleanup_function(_arf_cleanup);
                arf_have_registered_cleanup = 1;
            }

            arf_free_alloc = FLINT_MAX(64, arf_free_alloc * 2);
            arf_free_arr = flint_realloc(arf_free_arr,
                arf_free_alloc * sizeof(mp_ptr));
        }

        ptr[0] = alloc;
        arf_free_arr[arf_free_num++] = ptr;
    }
    else
    {
        flint_free(ptr);
    }
}

