/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include <string.h>

void
fmpz_vec_init(fmpz_vec_t vec, slong len)
{
    vec->length = vec->alloc = len;

    if (len == 0)
    {
        vec->entries = NULL;
    }
    else
    {
        vec->entries = _fmpz_vec_init(len);
    }
}

void
fmpz_vec_clear(fmpz_vec_t vec)
{
    _fmpz_vec_clear(vec->entries, vec->alloc);
}

void
fmpz_vec_fit_length(fmpz_vec_t vec, slong len)
{
    slong alloc = vec->alloc;

    if (len > alloc)
    {
        if (len < 2 * alloc)
            len = 2 * alloc;

        vec->entries = flint_realloc(vec->entries, len * sizeof(fmpz));
        memset(vec->entries + alloc, 0, (len - alloc) * sizeof(fmpz));
        vec->alloc = len;
    }
}

void
fmpz_vec_set_length(fmpz_vec_t vec, slong len)
{
    if (vec->length > len)
    {
        _fmpz_vec_zero(vec->entries + len, vec->length - len);
    }
    else if (vec->length < len)
    {
        fmpz_vec_fit_length(vec, len);
        _fmpz_vec_zero(vec->entries + vec->length, len - vec->length);
    }

    vec->length = len;
}

void
fmpz_vec_set(fmpz_vec_t res, const fmpz_vec_t src)
{
    if (res == src)
        return;

    fmpz_vec_set_length(res, src->length);
    _fmpz_vec_set(res->entries, src->entries, src->length);
}

void
fmpz_vec_append(fmpz_vec_t vec, const fmpz_t f)
{
    fmpz_vec_fit_length(vec, vec->length + 1);
    vec->length++;
    fmpz_set(vec->entries + vec->length - 1, f);
}

void
fmpz_vec_append_ui(fmpz_vec_t vec, ulong f)
{
    fmpz_vec_fit_length(vec, vec->length + 1);
    vec->length++;
    fmpz_set_ui(vec->entries + vec->length - 1, f);
}
