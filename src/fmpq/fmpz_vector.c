/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

/*
    after initialization one may:

        (1) set length to < 0: terms should not be stored, or
        (2) set want_alt_sum to +-1: terms should be add/sub to alt_sum, or
        (3) set limit to anything >= 0: stop generating terms after limit
        (4) set length to < 0 and set want_alt_sum to +-1: combination of (1) and (2)

    different combinations of these settings may or may not work
*/
void _fmpq_cfrac_list_init(_fmpq_cfrac_list_t v)
{
    v->array = NULL;
    v->length = 0;
    v->alloc = 0;
    v->limit = WORD_MAX;
    v->want_alt_sum = 0;
    fmpz_init(v->alt_sum);
}


void _fmpq_cfrac_list_clear(_fmpq_cfrac_list_t v)
{
    slong i;

    for (i = 0; i < v->alloc; i++)
        fmpz_clear(v->array + i);

    if (v->array)
        flint_free(v->array);

    fmpz_clear(v->alt_sum);
}


void _fmpq_cfrac_list_fit_length(_fmpq_cfrac_list_t v, slong len)
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


void _fmpq_cfrac_list_push_back(_fmpq_cfrac_list_t v, const fmpz_t a)
{
    if (v->want_alt_sum)
    {
        v->want_alt_sum *= -1;
        if (v->want_alt_sum > 0)
            fmpz_sub(v->alt_sum, v->alt_sum, a);
        else
            fmpz_add(v->alt_sum, v->alt_sum, a);
    }

    if (v->length < 0)
        return;

    _fmpq_cfrac_list_fit_length(v, v->length + 1);
    fmpz_set(v->array + v->length, a);
    v->length++;
    FLINT_ASSERT(v->length <= v->limit);
}


void _fmpq_cfrac_list_push_back_zero(_fmpq_cfrac_list_t v)
{
    v->want_alt_sum *= -1;

    if (v->length < 0)
        return;

    _fmpq_cfrac_list_fit_length(v, v->length + 1);
    fmpz_zero(v->array + v->length);
    v->length++;
    FLINT_ASSERT(v->length <= v->limit);
}


void _fmpq_cfrac_list_append_ui(_fmpq_cfrac_list_t v, const ulong * a, slong n)
{
    slong i;

    if (v->want_alt_sum)
    {
        ulong hi = 0, lo = 0;
        for (i = 0; i + 2 <= n; i += 2)
        {
            add_ssaaaa(hi,lo, hi,lo, 0,a[i]);
            sub_ddmmss(hi,lo, hi,lo, 0,a[i + 1]);
        }

        if (i < n)
        {
            add_ssaaaa(hi,lo, hi,lo, 0,a[i]);
        }

        if (v->want_alt_sum < 0)
        {
            hi = -hi - (lo != 0);
            lo = -lo;
        }

        if (i < n)
        {
            v->want_alt_sum *= -1;
        }

        if (hi == 0)
        {
            fmpz_add_ui(v->alt_sum, v->alt_sum, lo);
        }
        else if (lo != 0 && hi == -UWORD(1))
        {
            fmpz_sub_ui(v->alt_sum, v->alt_sum, -lo);
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_set_signed_uiui(t, hi, lo);
            fmpz_add(v->alt_sum, v->alt_sum, t);
            fmpz_clear(t);
        }
    }

    if (v->length < 0)
        return;

    _fmpq_cfrac_list_fit_length(v, v->length + n);
    for (i = 0; i < n; i++)
        fmpz_set_ui(v->array + v->length + i, a[i]);
    v->length += n;

    FLINT_ASSERT(v->length <= v->limit);
}

