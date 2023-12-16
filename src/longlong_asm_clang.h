/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONGLONG_ASM_H
#define LONGLONG_ASM_H

#if FLINT64
# define _STOR_ULONG __uint128_t
# define _STOR_SLONG __int128_t
#else
# define _STOR_ULONG unsigned long long int
# define _STOR_SLONG long long int
#endif

#if defined(_LONG_LONG_LIMB)
# define _FLINT_ADC __builtin_addcll
# define _FLINT_SBB __builtin_subcll
#else
# define _FLINT_ADC __builtin_addcl
# define _FLINT_SBB __builtin_subcl
#endif

#define add_ssaaaa(s1, s0, a1, a0, b1, b0) \
do \
{ \
    ulong _carry; \
    (s0) = _FLINT_ADC(a0, b0, 0, &_carry); \
    (s1) = _FLINT_ADC(a1, b1, _carry, &_carry); \
} while (0)

#define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
do \
{ \
    ulong _carry; \
    (s0) = _FLINT_ADC(a0, b0, 0, &_carry); \
    (s1) = _FLINT_ADC(a1, b1, _carry, &_carry); \
    (s2) = _FLINT_ADC(a2, b2, _carry, &_carry); \
} while (0)

#define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
do \
{ \
    ulong _carry; \
    (s0) = _FLINT_ADC(a0, b0, 0, &_carry); \
    (s1) = _FLINT_ADC(a1, b1, _carry, &_carry); \
    (s2) = _FLINT_ADC(a2, b2, _carry, &_carry); \
    (s3) = _FLINT_ADC(a3, b3, _carry, &_carry); \
} while (0)

#define sub_ddmmss(s1, s0, a1, a0, b1, b0) \
do \
{ \
    ulong _carry; \
    (s0) = _FLINT_SBB(a0, b0, 0, &_carry); \
    (s1) = _FLINT_SBB(a1, b1, _carry, &_carry); \
} while (0)

#define sub_dddmmmsss(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
do \
{ \
    ulong _carry; \
    (s0) = _FLINT_SBB(a0, b0, 0, &_carry); \
    (s1) = _FLINT_SBB(a1, b1, _carry, &_carry); \
    (s2) = _FLINT_SBB(a2, b2, _carry, &_carry); \
} while (0)

#define _mul_ppmm(big_type, small_type, r1, r0, u, v) \
do \
{ \
    big_type __r = (big_type) ((small_type) u) * (big_type) ((small_type) v); \
    (r0) = (small_type) __r; \
    (r1) = (small_type) (__r >> FLINT_BITS); \
} while (0)

#define umul_ppmm(r1, r0, u, v) _mul_ppmm(_STOR_ULONG, ulong, r1, r0, u, v)
#define smul_ppmm(r1, r0, u, v) _mul_ppmm(_STOR_SLONG, slong, r1, r0, u, v)

#endif
