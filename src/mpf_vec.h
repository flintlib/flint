/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPF_VEC_H
#define MPF_VEC_H

#ifdef MPF_VEC_INLINES_C
#define MPF_VEC_INLINE FLINT_DLL
#else
#define MPF_VEC_INLINE static __inline__
#endif

#include "flint.h"

typedef __mpf_struct mpf;

#ifdef __cplusplus
 extern "C" {
#endif

/*  Memory management  *******************************************************/

FLINT_DLL mpf * _mpf_vec_init(slong len, flint_bitcnt_t prec);

FLINT_DLL void _mpf_vec_clear(mpf * vec, slong len);

/*  Randomisation  ***********************************************************/

FLINT_DLL void _mpf_vec_randtest(mpf * f, flint_rand_t state, 
                        slong len, flint_bitcnt_t bits);

/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void _mpf_vec_zero(mpf * vec, slong len);

FLINT_DLL void _mpf_vec_set(mpf * vec1, const mpf * vec2, slong len2);

/*  Comparison  **************************************************************/

FLINT_DLL int _mpf_vec_equal(const mpf * vec1, const mpf * vec2, slong len);

FLINT_DLL int _mpf_vec_approx_equal(const mpf * vec1, const mpf * vec2, slong len, flint_bitcnt_t bits);

FLINT_DLL int _mpf_vec_is_zero(const mpf * vec, slong len);

/*  Addition  ****************************************************************/

FLINT_DLL void _mpf_vec_add(mpf * res, const mpf * vec1, const mpf * vec2, slong len2);

FLINT_DLL void _mpf_vec_sub(mpf * res, const mpf * vec1, const mpf * vec2, slong len2);

/*  Scalar multiplication  **************************************/

FLINT_DLL void _mpf_vec_scalar_mul_2exp(mpf * res, const mpf * vec, slong len, flint_bitcnt_t exp);

FLINT_DLL void _mpf_vec_scalar_mul_mpf(mpf * res, const mpf * vec, slong len, mpf_t c);

/*  Dot product and norm  **************************************/

FLINT_DLL void _mpf_vec_dot(mpf_t res, const mpf * vec1, const mpf * vec2, slong len2);

FLINT_DLL void _mpf_vec_norm(mpf_t res, const mpf * vec, slong len);

FLINT_DLL int _mpf_vec_dot2(mpf_t res, const mpf * vec1, const mpf * vec2, slong len2, flint_bitcnt_t prec);

FLINT_DLL void _mpf_vec_norm2(mpf_t res, const mpf * vec, slong len, flint_bitcnt_t prec);

#ifdef __cplusplus
}
#endif

#endif

