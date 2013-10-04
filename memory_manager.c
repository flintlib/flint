/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "flint.h"

static void flint_memory_error()
{
    flint_printf("Exception (FLINT memory_manager). Unable to allocate memory.\n");
    abort();
}

void * flint_malloc(size_t size)
{
    void * ptr = malloc(size);

    if (ptr == NULL)
        flint_memory_error();

    return ptr;
}

void * flint_realloc(void * ptr, size_t size)
{
    void * ptr2 = realloc(ptr, size);

    if (ptr2 == NULL)
        flint_memory_error();

    return ptr2;
}

void * flint_calloc(size_t num, size_t size)
{
    void * ptr = calloc(num, size);

    if (ptr == NULL)
        flint_memory_error();

    return ptr;
}

void flint_free(void * ptr)
{
    free(ptr);
}


FLINT_TLS_PREFIX size_t flint_num_cleanup_functions = 0;

FLINT_TLS_PREFIX flint_cleanup_function_t * flint_cleanup_functions = NULL;

void flint_register_cleanup_function(flint_cleanup_function_t cleanup_function)
{
    flint_cleanup_functions = flint_realloc(flint_cleanup_functions,
        (flint_num_cleanup_functions + 1) * sizeof(flint_cleanup_function_t));

    flint_cleanup_functions[flint_num_cleanup_functions] = cleanup_function;

    flint_num_cleanup_functions++;
}

void _fmpz_cleanup();

void flint_cleanup()
{
    size_t i;

    for (i = 0; i < flint_num_cleanup_functions; i++)
        flint_cleanup_functions[i]();

    flint_free(flint_cleanup_functions);
    flint_cleanup_functions = NULL;
    flint_num_cleanup_functions = 0;

    mpfr_free_cache();
    _fmpz_cleanup();
}


