/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef MPFR_VEC_H
#define MPFR_VEC_H

#ifdef MPFR_VEC_INLINES_C
#define MPFR_VEC_INLINE FLINT_DLL
#else
#define MPFR_VEC_INLINE static __inline__
#endif

#include <gmp.h>
#include <mpfr.h> 

#ifdef __cplusplus
 extern "C" {
#endif

FLINT_DLL mpfr * _mpfr_vec_init(slong length, mp_bitcnt_t prec);

FLINT_DLL void _mpfr_vec_clear(mpfr * vec, slong length);

FLINT_DLL void _mpfr_vec_randtest(mpfr * f, flint_rand_t state, slong len);

FLINT_DLL void _mpfr_vec_zero(mpfr * vec, slong length);

FLINT_DLL void _mpfr_vec_set(mpfr * vec1, const mpfr * vec2, slong length);

FLINT_DLL int _mpfr_vec_equal(const mpfr * vec1, const mpfr * vec2, slong len);

FLINT_DLL void _mpfr_vec_add(mpfr * res, const mpfr * vec1, const mpfr * vec2, slong length);

FLINT_DLL void _mpfr_vec_scalar_mul_2exp(mpfr * res, const mpfr * vec, slong length, mp_bitcnt_t exp);

FLINT_DLL void _mpfr_vec_scalar_mul_mpfr(mpfr * res, const mpfr * vec, slong length, mpfr_t c);

FLINT_DLL void _mpfr_vec_scalar_product(mpfr_t res, const mpfr * vec1, const mpfr * vec2, slong length);

#ifdef __cplusplus
}
#endif

#endif






