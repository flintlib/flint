/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* This is an interface intended for basic fmpz functionality as it reduces the
 * header size by a substantially. */

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

/* memory management **********************************************************/

FMPZ_MINI_INLINE
void fmpz_init(fmpz_t f)
{
    (*f) = WORD(0);
}

FLINT_DLL void _fmpz_clear_mpz(fmpz f);

FLINT_DLL mpz_mock_ptr _fmpz_new_mpz(void);

FMPZ_MINI_INLINE
void fmpz_clear(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f))
        _fmpz_clear_mpz(*f);
}

FMPZ_MINI_INLINE
void _fmpz_demote(fmpz_t f)
{
    if (COEFF_IS_MPZ(*f)) 
    {
        _fmpz_clear_mpz(*f);
        (*f) = WORD(0);
    }
}

FLINT_DLL mpz_mock_ptr _fmpz_promote_val(fmpz_t f);

FLINT_DLL void _fmpz_demote_val(fmpz_t f);

FLINT_DLL mpz_mock_ptr _fmpz_promote(fmpz_t f);

FMPZ_MINI_INLINE
void fmpz_init_set_ui(fmpz_t f, ulong g)
{
    if (g <= COEFF_MAX)
    {
        *f = g;
    }
    else
    {
        mpz_mock_ptr ptr;
        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        ptr->_mp_size = 1;
        ptr->_mp_d[0] = g;
    }
}

FMPZ_MINI_INLINE
void fmpz_init_set_si(fmpz_t f, slong g)
{
    if (COEFF_MIN <= g && g <= COEFF_MAX)
    {
        *f = g;
    }
    else
    {
        mpz_mock_ptr ptr;
        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        if (g < 0)
        {
            ptr->_mp_size = -1;
            ptr->_mp_d[0] = -g;
        }
        else
        {
            ptr->_mp_size = 1;
            ptr->_mp_d[0] = g;
        }
    }
}

/* basic assignment and manipulation ******************************************/

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

FMPZ_MINI_INLINE
int fmpz_is_pm1(const fmpz_t f)
{
   return (*f == 1 || *f == -1);
}

FLINT_DLL int fmpz_equal(const fmpz_t f, const fmpz_t g);
FLINT_DLL int fmpz_equal_si(const fmpz_t f, slong g);
FLINT_DLL int fmpz_equal_ui(const fmpz_t f, ulong g);

FLINT_DLL void fmpz_set(fmpz_t f, const fmpz_t g);

FMPZ_MINI_INLINE void
fmpz_set_si(fmpz_t f, slong g)
{
    if (g < COEFF_MIN || g > COEFF_MAX)
    {
        mpz_mock_ptr mf = _fmpz_promote(f);
        if (g < 0)
        {
            mf->_mp_size = -1;
            mf->_mp_d[0] = -g;
        }
        else
        {
            mf->_mp_size = 1;
            mf->_mp_d[0] = g;
        }
    }
    else
    {
        _fmpz_demote(f);
        *f = g;
    }
}

FMPZ_MINI_INLINE void
fmpz_set_ui(fmpz_t f, ulong g)
{
    if (g > COEFF_MAX)
    {
        mpz_mock_ptr mf = _fmpz_promote(f);
        mf->_mp_size = 1;
        mf->_mp_d[0] = g;
    }
    else
    {
        _fmpz_demote(f);
        *f = g;
    }
}

FMPZ_MINI_INLINE void
fmpz_neg_ui(fmpz_t f, ulong g)
{
    if (g > COEFF_MAX)
    {
        mpz_mock_ptr mf = _fmpz_promote(f);
        mf->_mp_size = -1;
        mf->_mp_d[0] = g;
    }
    else
    {
        _fmpz_demote(f);
        *f = -(slong) g;
    }
}

FMPZ_MINI_INLINE
void fmpz_swap(fmpz_t f, fmpz_t g)
{
    if (f != g)
    {
        fmpz t = *f;
        *f = *g;
        *g = t;
    }
}

FLINT_DLL slong fmpz_get_si(const fmpz_t f);
FLINT_DLL ulong fmpz_get_ui(const fmpz_t f);

FLINT_DLL void fmpz_abs(fmpz_t f1, const fmpz_t f2);

FLINT_DLL void fmpz_neg(fmpz_t f1, const fmpz_t f2);

FLINT_DLL void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);
FLINT_DLL void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x);

FLINT_DLL void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h);
FLINT_DLL void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x);

FMPZ_MINI_INLINE
void fmpz_add_si(fmpz_t f, const fmpz_t g, slong x)
{
    if (x >= 0)
        fmpz_add_ui(f, g, (ulong) x);
    else
        fmpz_sub_ui(f, g, (ulong) -x);
}

FMPZ_MINI_INLINE
void fmpz_sub_si(fmpz_t f, const fmpz_t g, slong x)
{
    if (x >= 0)
        fmpz_sub_ui(f, g, (ulong) x);
    else
        fmpz_add_ui(f, g, (ulong) -x);
}

FLINT_DLL void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);
FLINT_DLL void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x);
FLINT_DLL void fmpz_mul_si(fmpz_t f, const fmpz_t g, slong x);

FLINT_DLL void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp);

FLINT_DLL void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h);
FLINT_DLL void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slong h);
FLINT_DLL void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h);

/* comparison *****************************************************************/

FLINT_DLL int fmpz_cmp(const fmpz_t f, const fmpz_t g);
FLINT_DLL int fmpz_cmp_ui(const fmpz_t f, ulong g);
FLINT_DLL int fmpz_cmp_si(const fmpz_t f, slong g);

FMPZ_MINI_INLINE
int fmpz_sgn(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
        return FLINT_SGN(*f);
    else
    {
        mpz_mock_ptr mf = COEFF_TO_PTR(*f);
        return (mf->_mp_size > 0) ? 1 : -1;
    }
}

FLINT_DLL int fmpz_fits_si(const fmpz_t f);
FLINT_DLL int fmpz_abs_fits_ui(const fmpz_t f);

FLINT_DLL mp_mock_size_t fmpz_size(const fmpz_t f);
FLINT_DLL flint_bitcnt_t fmpz_bits(const fmpz_t f);

/* miscellaneous **************************************************************/

FLINT_DLL int fmpz_divisible(const fmpz_t f, const fmpz_t g);

FLINT_DLL void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h);

FLINT_DLL void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h);

/* NOTE: Only for padic.h and padic_poly.h. Hopefully there is a way to remove
 * this from here. */
FLINT_DLL slong _fmpz_remove(fmpz_t x, const fmpz_t f, double finv);
FLINT_DLL slong fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f);

#ifdef __cplusplus
}
#endif

#endif
