/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPFR_VEC_H
#define MPFR_VEC_H

#ifdef MPFR_VEC_INLINES_C
#define MPFR_VEC_INLINE
#else
#define MPFR_VEC_INLINE static inline
#endif

#include <mpfr.h>
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Soon to be deprecated */
#define flint_mpfr __mpfr_struct

mpfr_ptr _mpfr_vec_init(slong length, flint_bitcnt_t prec);

void _mpfr_vec_clear(mpfr_ptr vec, slong length);

void _mpfr_vec_randtest(mpfr_ptr f, flint_rand_t state, slong len);

void _mpfr_vec_zero(mpfr_ptr vec, slong length);

void _mpfr_vec_set(mpfr_ptr vec1, mpfr_srcptr vec2, slong length);

int _mpfr_vec_equal(mpfr_srcptr vec1, mpfr_srcptr vec2, slong len);

void _mpfr_vec_add(mpfr_ptr res, mpfr_srcptr vec1, mpfr_srcptr vec2, slong length);

void _mpfr_vec_scalar_mul_2exp(mpfr_ptr res, mpfr_srcptr vec, slong length, flint_bitcnt_t exp);

void _mpfr_vec_scalar_mul_mpfr(mpfr_ptr res, mpfr_srcptr vec, slong length, mpfr_t c);

void _mpfr_vec_scalar_product(mpfr_t res, mpfr_srcptr vec1, mpfr_srcptr vec2, slong length);

#ifdef __cplusplus
}
#endif

#endif
