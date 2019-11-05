/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"


void _fmpz_vector_init(_fmpz_vector_t v)
{
    v->array = NULL;
    v->length = 0;
    v->alloc = 0;
    v->limit = WORD_MAX;
}


void _fmpz_vector_init_nowrite(_fmpz_vector_t v)
{
    v->array = NULL;
    v->length = -1;
    v->alloc = 0;
    v->limit = WORD_MAX;
}


void _fmpz_vector_clear(_fmpz_vector_t v)
{
    slong i;

    for (i = 0; i < v->alloc; i++)
        fmpz_clear(v->array + i);

    if (v->array)
        flint_free(v->array);
}


void _fmpz_vector_fit_length(_fmpz_vector_t v, slong len)
{
    if (len <= v->alloc)
        return;

    if (v->alloc > 0)
    {
        len = FLINT_MAX(len, v->alloc + v->alloc/2);

        v->array = (fmpz *) flint_realloc(v->array, len * sizeof(fmpz));
        FLINT_ASSERT(len > v->alloc);
        flint_mpn_zero((mp_ptr) (v->array + v->alloc), len - v->alloc);
    }
    else
    {
        v->array = (fmpz *) flint_calloc(len, sizeof(fmpz));
    }

    v->alloc = len;
}


void _fmpz_vector_push_back(_fmpz_vector_t v, const fmpz_t a)
{
    if (v->length < 0)
        return;

    _fmpz_vector_fit_length(v, v->length + 1);
    fmpz_set(v->array + v->length, a);
    v->length++;
    FLINT_ASSERT(v->length <= v->limit);
}


void _fmpz_vector_push_back_zero(_fmpz_vector_t v)
{
    if (v->length < 0)
        return;

    _fmpz_vector_fit_length(v, v->length + 1);
    fmpz_zero(v->array + v->length);
    v->length++;

    FLINT_ASSERT(v->length <= v->limit);
}


void _fmpz_vector_append_ui(_fmpz_vector_t v, const ulong * a, slong n)
{
    slong i;

    if (v->length < 0)
        return;

    _fmpz_vector_fit_length(v, v->length + n);
    for (i = 0; i < n; i++)
        fmpz_set_ui(v->array + v->length + i, a[i]);
    v->length += n;

    FLINT_ASSERT(v->length <= v->limit);
}

