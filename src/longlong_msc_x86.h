/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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

/* Byte swap */
# define byte_swap(x) do { (x) = _byteswap_ulong(x); } while (0)

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

/* Byte swap */
# define byte_swap(x) do { (x) = _byteswap_uint64(x); } while (0)

/* Addition and subtraction */
# define _FLINT_ADC _addcarry_u64
# define _FLINT_SBB _subborrow_u64

/* Multiplication */
#define umul_ppmm(r1, r0, u, v) do { (r0) = _umul128(u, v, &(r1)); } while (0)
#define smul_ppmm(r1, r0, u, v) do { (r0) = _mul128(u, v, &(r1)); } while (0)

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

/* Division */
#define udiv_qrnnd(q, r, n1, n0, dx) do { (q) = _FLINT_DIV(n1, n0, dx, &(r)); } while (0)
#define sdiv_qrnnd(q, r, n1, n0, dx) do { (q) = _FLINT_IDIV(n1, n0, dx, &(r)); } while (0)

#endif
