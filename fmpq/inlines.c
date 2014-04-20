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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#include <mpfr.h>
#define ulong mp_limb_t
#include "flint.h"
#include "fmpq.h"

void __fmpq_init(fmpq_t x)
{
    x->num = WORD(0);
    x->den = WORD(1);
}

void __fmpq_clear(fmpq_t x)
{
    fmpz_clear(fmpq_numref(x));
    fmpz_clear(fmpq_denref(x));
}

fmpq * __fmpq_vec_init(slong n)
{
    fmpq * v = (fmpq *) flint_malloc(sizeof(fmpq) * n);
    slong i;

    for (i = 0; i < n; i++)
        fmpq_init(v + i);

    return v;
}

void __fmpq_vec_clear(fmpq * vec, slong n)
{
    _fmpz_vec_clear((fmpz *) vec, 2 * n);
}

void __fmpq_zero(fmpq_t res)
{
    fmpz_zero(fmpq_numref(res));
    fmpz_one(fmpq_denref(res));
}

void __fmpq_one(fmpq_t res)
{
    fmpz_one(fmpq_numref(res));
    fmpz_one(fmpq_denref(res));
}

int __fmpq_equal(const fmpq_t x, const fmpq_t y)
{
    return fmpz_equal(fmpq_numref(x), fmpq_numref(y)) &&
           fmpz_equal(fmpq_denref(x), fmpq_denref(y));
}

int __fmpq_sgn(const fmpq_t x)
{
    return fmpz_sgn(fmpq_numref(x));
}

int __fmpq_is_zero(const fmpq_t x)
{
    return fmpz_is_zero(fmpq_numref(x));
}

int __fmpq_is_one(const fmpq_t x)
{
    return fmpz_is_one(fmpq_numref(x)) && fmpz_is_one(fmpq_denref(x));
}

void __fmpq_set(fmpq_t dest, const fmpq_t src)
{
    fmpz_set(fmpq_numref(dest), fmpq_numref(src));
    fmpz_set(fmpq_denref(dest), fmpq_denref(src));
}

void __fmpq_swap(fmpq_t op1, fmpq_t op2)
{
    fmpz_swap(fmpq_numref(op1), fmpq_numref(op2));
    fmpz_swap(fmpq_denref(op1), fmpq_denref(op2));
}

void __fmpq_neg(fmpq_t dest, const fmpq_t src)
{
    fmpz_neg(fmpq_numref(dest), fmpq_numref(src));
    fmpz_set(fmpq_denref(dest), fmpq_denref(src));
}

void __fmpq_abs(fmpq_t dest, const fmpq_t src)
{
    fmpz_abs(fmpq_numref(dest), fmpq_numref(src));
    fmpz_set(fmpq_denref(dest), fmpq_denref(src));
}

void __fmpq_set_mpq(fmpq_t dest, const mpq_t src)
{
    fmpz_set_mpz(fmpq_numref(dest), mpq_numref(src));
    fmpz_set_mpz(fmpq_denref(dest), mpq_denref(src));
}

void __fmpq_get_mpq(mpq_t dest, const fmpq_t src)
{
    fmpz_get_mpz(mpq_numref(dest), fmpq_numref(src));
    fmpz_get_mpz(mpq_denref(dest), fmpq_denref(src));
}
