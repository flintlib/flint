/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MP_REAL_TEST_MAX_ERR_H
#define MP_REAL_TEST_MAX_ERR_H

/* The documented maxima of the err outputs (for the bitwise functions,
   r = 0 stands for the tuned default, at most 768). */
#define DOC_R(r) ((slong) ((r) ? (r) : 768))
#define DOC_EXP_RS_MAX_ERR(n) 10
#define DOC_SIN_RS_MAX_ERR(n) 15
#define DOC_COS_RS_MAX_ERR(n) 15
#define DOC_SIN_COS_RS_MAX_ERR(n) 15
#define DOC_SINH_RS_MAX_ERR(n) 15
#define DOC_COSH_RS_MAX_ERR(n) 15
#define DOC_SINH_COSH_RS_MAX_ERR(n) 15
#define DOC_ATAN_RS_MAX_ERR(n) 15
#define DOC_ATANH_RS_MAX_ERR(n) 15
#define DOC_EXP_BITWISE_RS_MAX_ERR(n, r) (9 * DOC_R(r) + 100)
#define DOC_LOG1P_BITWISE_RS_MAX_ERR(n, r) (3 * DOC_R(r) + 64)
#define DOC_SIN_COS_BITWISE_RS_MAX_ERR(n, r) (6 * DOC_R(r) + 128)
#define DOC_ATAN_BITWISE_RS_MAX_ERR(n, r) (4 * DOC_R(r) + 64)
#define DOC_TAN_BITWISE_RS_MAX_ERR(n, r) (8 * DOC_R(r) + 256)
#define DOC_EXP_REDUCED_MAX_ERR 96
#define DOC_SIN_COS_REDUCED_MAX_ERR 96
#define DOC_EXP_NOTAB_MAX_ERR 128
#define DOC_SIN_COS_NOTAB_MAX_ERR 128
#define DOC_EXP_DIOPHANTINE_MAX_ERR 3
#define DOC_SIN_COS_DIOPHANTINE_MAX_ERR 4
#define DOC_NEGLOG_NEWTON_MAX_ERR 2
#define DOC_ATAN_NEWTON_MAX_ERR 2

/* the err output of the last call under test; a comparison against
   TEST_ERR(documented maximum) fails if the returned bound exceeds the
   maximum or the measured error exceeds the returned bound (the
   measurements go through mag_t upper bounds with 30-bit mantissas,
   hence the relative allowance of 2^-20 for their own rounding) */
static ulong test_err_out;
#define TEST_ERR(m) (test_err_out <= (ulong) (m) \
    ? (double) test_err_out * (1.0 + 0x1p-20) : -1.0)

#endif
