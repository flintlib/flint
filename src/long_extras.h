/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONG_EXTRAS_H
#define LONG_EXTRAS_H

#ifdef LONG_EXTRAS_INLINES_C
#define LONG_EXTRAS_INLINE
#else
#define LONG_EXTRAS_INLINE static inline
#endif

#include "longlong.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Properties ****************************************************************/

size_t z_sizeinbase(slong n, int b);

/* Checked arithmetic ********************************************************/

LONG_EXTRAS_INLINE int z_mul_checked(slong * a, slong b, slong c)
{
#if defined(__GNUC__)
    return __builtin_mul_overflow(b, c, a);
#else
	ulong ahi, alo;
	smul_ppmm(ahi, alo, b, c);
	*a = alo;
	return FLINT_SIGN_EXT(alo) != ahi;
#endif
}

LONG_EXTRAS_INLINE int z_add_checked(slong * a, slong b, slong c)
{
#if defined(__GNUC__)
    return __builtin_add_overflow(b, c, a);
#else
    int of = (b > 0 && c > WORD_MAX - b) || (b < 0 && c < WORD_MIN - b);
    *a = b + c;
    return of;
#endif
}

LONG_EXTRAS_INLINE
int z_mat22_det_is_negative(slong m11, slong m12, slong m21, slong m22)
{
    ulong t1, t2, t3, t4;
    smul_ppmm(t1, t2, m11, m22);
    smul_ppmm(t3, t4, m12, m21);
    sub_ddmmss(t1, t2, t1, t2, t3, t4);
    return FLINT_SIGN_EXT(t1) != 0;
}

/* Randomisation  ************************************************************/

slong z_randtest(flint_rand_t state);

slong z_randtest_not_zero(flint_rand_t state);

slong z_randint(flint_rand_t state, ulong limit);

/*****************************************************************************/

int z_kronecker(slong a, slong n);

#ifdef __cplusplus
}
#endif

#endif
