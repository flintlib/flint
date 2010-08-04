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

    Copyright (C) 2010, William Hart

    We acknowledge the generous assistance of Roman Pearce in helping us 
    understand his heap multiplication algorithm, and in particular for
    providing us with a reference implementation to look at.

******************************************************************************/

#ifndef FMPZ_MPOLY_H
#define FMPZ_MPOLY_H

#include <mpir.h>
#include "fmpz.h"

typedef struct
{
    fmpz * coeffs;
    fmpz * exps;
    long alloc;
    long length;
    long vars;  /* number of variables */
    mp_bitcnt_t ebits;  /* width of exponent bitfield (per variable) */
} fmpz_mpoly_struct;

typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

typedef struct fmpz_mpoly_entry_t fmpz_mpoly_entry_t;

struct fmpz_mpoly_entry_t
{
    fmpz_mpoly_entry_t * next;  /* linked list of entries */
    unsigned int i1;  /* index of coefficient in first poly */
    unsigned int i2;  /* index of coefficient in second poly */
};

typedef struct
{
    fmpz exp;
    fmpz_mpoly_entry_t * entry; 
} fmpz_mpoly_heap_t;

void fmpz_mpoly_init(fmpz_mpoly_t poly, long vars, mp_bitcnt_t ebits);

void fmpz_mpoly_init2(fmpz_mpoly_t poly, long alloc, 
                                        long vars, mp_bitcnt_t ebits);

void fmpz_mpoly_realloc(fmpz_mpoly_t poly, long alloc);

void fmpz_mpoly_fit_length(fmpz_mpoly_t poly, long len);

void fmpz_mpoly_clear(fmpz_mpoly_t poly);

static __inline__
void _fmpz_mpoly_truncate(fmpz_mpoly_t poly, long len)
{
    if (poly->length > len)  /* only truncate if necessary */
    {
        long i;
        for (i = len; i < poly->length; i++)
        {
            _fmpz_demote(poly->coeffs + i);
            _fmpz_demote(poly->exps + i);
        }
        poly->length = len;
    }
}

void fmpz_mpoly_reheapify(fmpz_mpoly_heap_t * heap, ulong * n);

void fmpz_mpoly_heap_insert(fmpz_mpoly_heap_t * heap, ulong * n, 
                                 fmpz_mpoly_entry_t * entry, fmpz exp);

void fmpz_mpoly_mul_heap(fmpz_mpoly_t res, fmpz_mpoly_t poly1, 
                                                   fmpz_mpoly_t poly2);

#endif

