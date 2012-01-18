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
    printf("FLINT error: unable to allocate memory\n");
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

