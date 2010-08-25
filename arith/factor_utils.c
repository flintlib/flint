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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "mpn_extras.h"
#include "ulong_extras.h"



void fmpz_factor_init(fmpz_factor_t factor)
{
    factor->sign = 0;
    factor->p = NULL;
    factor->exp = NULL;
    factor->length = 0;
    factor->alloc = 0;
}

void fmpz_factor_clear(fmpz_factor_t factor)
{
    _fmpz_vec_clear(factor->p, factor->alloc);
    _fmpz_vec_clear(factor->exp, factor->alloc);
}

void _fmpz_factor_fit_length(fmpz_factor_t factor, long len)
{
    if (len > factor->alloc)
    {
        if (len < 2 * factor->alloc)
            len = 2 * factor->alloc;
        factor->p = (fmpz *) realloc(factor->p, len * sizeof(fmpz));
        factor->exp = (fmpz *) realloc(factor->exp, len * sizeof(fmpz));
        if (len > factor->alloc)
        {
            mpn_zero((mp_ptr)(factor->p + factor->alloc), len-factor->alloc);
            mpn_zero((mp_ptr)(factor->exp + factor->alloc), len-factor->alloc);
        }
        factor->alloc = len;
    }
}

void _fmpz_factor_append_ui(fmpz_factor_t factor, ulong p, ulong exp)
{
    _fmpz_factor_fit_length(factor, factor->length + 1);
    fmpz_set_ui(factor->p + factor->length, p);
    fmpz_set_ui(factor->exp + factor->length, exp);
    factor->length++;
}

void _fmpz_factor_set_length(fmpz_factor_t factor, long newlen)
{
    if (factor->length > newlen)
    {
        long i;
        for (i = newlen; i < factor->length; i++)
        {
            _fmpz_demote(factor->p + i); 
            _fmpz_demote(factor->exp + i); 
        }
    }
    factor->length = newlen;
}

void _fmpz_factor_si(fmpz_factor_t factor, long n)
{
    _fmpz_factor_set_length(factor, 0);
    if (n < 0)
    {
        _fmpz_factor_extend_factor_n(factor, -n);
        factor->sign = -1;
        return;
    }
    factor->sign = 1;
    _fmpz_factor_extend_factor_n(factor, n);
}

void _fmpz_factor_extend_factor_n(fmpz_factor_t factor, ulong n)
{
    int i, len;
    n_factor_t nfac;

    if (n == 0)
    {
        _fmpz_factor_set_length(factor, 0);
        factor->sign = 0;
        return;
    }

    n_factor_init(&nfac);
    n_factor(&nfac, n, 0);

    len = factor->length;

    _fmpz_factor_fit_length(factor, len + nfac.num);
    _fmpz_factor_set_length(factor, len + nfac.num);

    for (i = 0; i < nfac.num; i++)
    {
        fmpz_set_ui(factor->p + len + i, nfac.p[i]);
        fmpz_set_ui(factor->exp + len + i, nfac.exp[i]);
    }
}

void fmpz_factor_print(fmpz_factor_t factor)
{
    int i;

    if (factor->sign == 0)
    {
        printf("0");
        return;
    }

    if (factor->sign == -1)
    {
        if (factor->length)
            printf("-1 * ");
        else
            printf("-1");
    }

    for (i = 0; i < factor->length; i++)
    {
        fmpz_print(factor->p + i);
        if (factor->exp[i] != 1UL)
        {
            printf("^");
            fmpz_print(factor->exp + i);
        }
        if (i != factor->length - 1)
            printf(" * ");
    }
}
