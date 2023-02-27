/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FMPZ_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
 
fmpz * __new_fmpz()
{
    return flint_calloc(sizeof(fmpz), 1);
}

void __free_fmpz(fmpz * f)
{
   _fmpz_demote(f);
   flint_free(f);
}   

void __fmpz_set_si(fmpz_t f, slong val)
{
    if (val < COEFF_MIN || val > COEFF_MAX) /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        flint_mpz_set_si(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

void __fmpz_set_ui(fmpz_t f, ulong val)
{
    if (val > COEFF_MAX)        /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        flint_mpz_set_ui(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

void __fmpz_init(fmpz_t f)
{
	(*f) = WORD(0);
}

void __fmpz_init_set_ui(fmpz_t f, ulong g)
{
    if (g <= COEFF_MAX)
    {
        *f = g;
    }
    else
    {
        __mpz_struct *ptr;

        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        flint_mpz_set_ui(ptr, g);
    }
}

void __fmpz_clear(fmpz_t f)
{
	_fmpz_demote(f);
}

int __fmpz_lt(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) < 0;
}

int __fmpz_gt(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) > 0;
}

int __fmpz_lte(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) <= 0;
}

int __fmpz_gte(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) >= 0;
}

int __fmpz_eq(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) == 0;
}

int __fmpz_neq(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) != 0;
}

void __fmpz_init_set(fmpz_t f, const fmpz_t g)
{
    if (!COEFF_IS_MPZ(*g))
    {
        *f = *g;
    }
    else
    {
        __mpz_struct *ptr;

        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        mpz_set(ptr, COEFF_TO_PTR(*g));
    }
}

void __fmpz_neg(fmpz_t f1, const fmpz_t f2)
{
    if (!COEFF_IS_MPZ(*f2))     /* coeff is small */
    {
        fmpz t = -*f2;          /* Need to save value in case of aliasing */
        _fmpz_demote(f1);
        *f1 = t;
    }
    else                        /* coeff is large */
    {
        /* No need to retain value in promotion, as if aliased, both already large */
        __mpz_struct * mf1 = _fmpz_promote(f1);
        mpz_neg(mf1, COEFF_TO_PTR(*f2));
    }
}
