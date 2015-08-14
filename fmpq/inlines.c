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

    Copyright (C) 2009 William Hart

******************************************************************************/

#define FMPQ_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpq.h"

void fmpq_numerator(fmpz_t n, const fmpq_t q)
{
   fmpz_set(n, fmpq_numref(q));
}

void fmpq_denominator(fmpz_t n, const fmpq_t q)
{
   fmpz_set(n, fmpq_denref(q));
}

fmpz * fmpq_numerator_ptr(fmpq_t q)
{
   return fmpq_numref(q);
}

fmpz * fmpq_denominator_ptr(fmpq_t q)
{
   return fmpq_denref(q);
}

int fmpq_equal_si(fmpq_t q, slong n)
{
   return q->num == n && q->den == WORD(1);
}

int fmpq_equal_fmpz(fmpq_t q, fmpz_t n)
{
   return fmpz_equal(fmpq_numref(q), n) && q->den == WORD(1);
}

