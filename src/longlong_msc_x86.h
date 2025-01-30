/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONGLONG_MSC_H
#define LONGLONG_MSC_H

#include <stdlib.h>
#include <intrin.h>
#include <immintrin.h>

#if defined(_M_IX86) || FLINT_BITS == 32

/* Trailing and leading zeros */
# define flint_clz _lzcnt_u32
# define flint_ctz _tzcnt_u32

/* Addition and subtraction */
# define _FLINT_ADC _addcarry_u32
# define _FLINT_SBB _subborrow_u32

/* Multiplication */
# define umul_ppmm(r1, r0, u, v) \
do \
{ \
    unsigned __int64 _tmp = __emulu(u, v); \
    (r0) = _tmp; \
    (r1) = _tmp >> 32; \
} while (0)

# define smul_ppmm(r1, r0, u, v) \
do \
{ \
    __int64 _tmp = __emul(u, v); \
    (r0) = _tmp; \
    (r1) = _tmp >> 32; \
} while (0)

/* Division */
# define _FLINT_DIV _udiv64
# define _FLINT_IDIV _div64

#else

/* Trailing and leading zeros */
# define flint_clz _lzcnt_u64
# define flint_ctz _tzcnt_u64

/* Addition and subtraction */
# define _FLINT_ADC _addcarry_u64
# define _FLINT_SBB _subborrow_u64

/* Multiplication */
#define umul_ppmm(r1, r0, u, v) do { (r0) = _umul128(u, v, (unsigned __int64 *) &(r1)); } while (0)
#define smul_ppmm(r1, r0, u, v) do { (r0) = _mul128(u, v, (__int64 *) &(r1)); } while (0)

/* Division */
# define _FLINT_DIV _udiv128
# define _FLINT_IDIV _div128

#endif

/* Addition and subtraction */
#define add_ssaaaa(s1, s0, a1, a0, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_ADC(0, a0, b0, &(s0)); \
    _FLINT_ADC(_carry, a1, b1, &(s1)); \
} while (0)

#define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_ADC(0, a0, b0, &(s0)); \
    _carry = _FLINT_ADC(_carry, a1, b1, &(s1)); \
    _FLINT_ADC(_carry, a2, b2, &(s2)); \
} while (0)

#define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_ADC(0, a0, b0, &(s0)); \
    _carry = _FLINT_ADC(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_ADC(_carry, a2, b2, &(s2)); \
    _FLINT_ADC(_carry, a3, b3, &(s3)); \
} while (0)

#define add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_ADC(0, a0, b0, &(s0)); \
    _carry = _FLINT_ADC(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_ADC(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_ADC(_carry, a3, b3, &(s3)); \
    _FLINT_ADC(_carry, a4, b4, &(s4)); \
} while (0)

#define add_ssssssaaaaaaaaaaaa(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_ADC(0, a0, b0, &(s0)); \
    _carry = _FLINT_ADC(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_ADC(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_ADC(_carry, a3, b3, &(s3)); \
    _carry = _FLINT_ADC(_carry, a4, b4, &(s4)); \
    _FLINT_ADC(_carry, a5, b5, &(s5)); \
} while (0)

#define add_sssssssaaaaaaaaaaaaaa(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_ADC(0, a0, b0, &(s0)); \
    _carry = _FLINT_ADC(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_ADC(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_ADC(_carry, a3, b3, &(s3)); \
    _carry = _FLINT_ADC(_carry, a4, b4, &(s4)); \
    _carry = _FLINT_ADC(_carry, a5, b5, &(s5)); \
    _FLINT_ADC(_carry, a6, b6, &(s6)); \
} while (0)

#define add_ssssssssaaaaaaaaaaaaaaaa(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_ADC(0, a0, b0, &(s0)); \
    _carry = _FLINT_ADC(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_ADC(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_ADC(_carry, a3, b3, &(s3)); \
    _carry = _FLINT_ADC(_carry, a4, b4, &(s4)); \
    _carry = _FLINT_ADC(_carry, a5, b5, &(s5)); \
    _carry = _FLINT_ADC(_carry, a6, b6, &(s6)); \
    _FLINT_ADC(_carry, a7, b7, &(s7)); \
} while (0)


#define sub_ddmmss(s1, s0, a1, a0, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_SBB(0, a0, b0, &(s0)); \
    _FLINT_SBB(_carry, a1, b1, &(s1)); \
} while (0)

#define sub_dddmmmsss(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_SBB(0, a0, b0, &(s0)); \
    _carry = _FLINT_SBB(_carry, a1, b1, &(s1)); \
    _FLINT_SBB(_carry, a2, b2, &(s2)); \
} while (0)

#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_SBB(0, a0, b0, &(s0)); \
    _carry = _FLINT_SBB(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_SBB(_carry, a2, b2, &(s2)); \
    _FLINT_SBB(_carry, a3, b3, &(s3)); \
} while (0)

#define sub_dddddmmmmmsssss(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_SBB(0, a0, b0, &(s0)); \
    _carry = _FLINT_SBB(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_SBB(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_SBB(_carry, a3, b3, &(s3)); \
    _FLINT_SBB(_carry, a4, b4, &(s4)); \
} while (0)

#define sub_ddddddmmmmmmssssss(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_SBB(0, a0, b0, &(s0)); \
    _carry = _FLINT_SBB(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_SBB(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_SBB(_carry, a3, b3, &(s3)); \
    _carry = _FLINT_SBB(_carry, a4, b4, &(s4)); \
    _FLINT_SBB(_carry, a5, b5, &(s5)); \
} while (0)

#define sub_dddddddmmmmmmmsssssss(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_SBB(0, a0, b0, &(s0)); \
    _carry = _FLINT_SBB(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_SBB(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_SBB(_carry, a3, b3, &(s3)); \
    _carry = _FLINT_SBB(_carry, a4, b4, &(s4)); \
    _carry = _FLINT_SBB(_carry, a5, b5, &(s5)); \
    _FLINT_SBB(_carry, a6, b6, &(s6)); \
} while (0)


#define sub_ddddddddmmmmmmmmssssssss(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) \
do \
{ \
    unsigned char _carry; \
    _carry = _FLINT_SBB(0, a0, b0, &(s0)); \
    _carry = _FLINT_SBB(_carry, a1, b1, &(s1)); \
    _carry = _FLINT_SBB(_carry, a2, b2, &(s2)); \
    _carry = _FLINT_SBB(_carry, a3, b3, &(s3)); \
    _carry = _FLINT_SBB(_carry, a4, b4, &(s4)); \
    _carry = _FLINT_SBB(_carry, a5, b5, &(s5)); \
    _carry = _FLINT_SBB(_carry, a6, b6, &(s6)); \
    _FLINT_SBB(_carry, a7, b7, &(s7)); \
} while (0)

/* Division */
#define udiv_qrnnd(q, r, n1, n0, dx) do { (q) = _FLINT_DIV(n1, n0, dx, &(r)); } while (0)
#define sdiv_qrnnd(q, r, n1, n0, dx) do { (q) = _FLINT_IDIV(n1, n0, dx, &(r)); } while (0)

#endif
