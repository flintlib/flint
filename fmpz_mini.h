/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MINI_H
#define FMPZ_MINI_H

#ifdef FMPZ_MINI_INLINES_C
#define FMPZ_MINI_INLINE FLINT_DLL
#else
#define FMPZ_MINI_INLINE static __inline__
#endif

#include "flint.h"
#include "fmpz-conversions.h"

#ifdef __cplusplus
extern "C" {
#endif

FMPZ_MINI_INLINE
void fmpz_init(fmpz_t f)
{
    (*f) = WORD(0);
}

FLINT_DLL void _fmpz_clear_mpz(fmpz f);

FMPZ_MINI_INLINE
void fmpz_clear(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f))
        _fmpz_clear_mpz(*f);
}

FMPZ_MINI_INLINE
void _fmpz_demote(fmpz_t f)
{
    /* 
       warning, if fmpz_demote changes, fmpz_zero must
       also be changed to match
    */
    if (COEFF_IS_MPZ(*f)) 
    {
        _fmpz_clear_mpz(*f);
        (*f) = WORD(0);
    }
}

FMPZ_MINI_INLINE
void fmpz_zero(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f))
        _fmpz_clear_mpz(*f);
    *f = WORD(0);
}

FMPZ_MINI_INLINE
void fmpz_one(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f)) 
        _fmpz_clear_mpz(*f);
    *f = WORD(1);
}

FMPZ_MINI_INLINE
int fmpz_is_zero(const fmpz_t f)
{
    return (*f == 0);
}

FMPZ_MINI_INLINE
int fmpz_is_one(const fmpz_t f)
{
    return (*f == 1);
}

FLINT_DLL int fmpz_equal(const fmpz_t f, const fmpz_t g);

#ifdef __cplusplus
}
#endif

#endif
