/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef CALCIUM_H
#define CALCIUM_H

#ifdef CALCIUM_INLINES_C
#define CALCIUM_INLINE
#else
#define CALCIUM_INLINE static __inline__
#endif

#include <stdio.h>
#include "flint/flint.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Global library definitions */

const char * calcium_version(void);

double calcium_test_multiplier(void);

/* Triple-valued logic */

typedef enum
{
    T_TRUE,
    T_FALSE,
    T_UNKNOWN
} truth_t;

CALCIUM_INLINE void truth_print(truth_t t)
{
    if (t == T_TRUE) flint_printf("T_TRUE");
    if (t == T_FALSE) flint_printf("T_FALSE");
    if (t == T_UNKNOWN) flint_printf("T_UNKNOWN");
}

/* IDs for builtin mathematical functions and constants */
typedef enum
{
    /* Special case for representing qqbar instances */
    CA_QQBar,
    /* Arithmetic */
    CA_Neg,
    CA_Add,
    CA_Sub,
    CA_Mul,
    CA_Div,
    /* Roots */
    CA_Sqrt,
    CA_Cbrt,
    CA_Root,
    /* Complex parts */
    CA_Floor,
    CA_Ceil,
    CA_Abs,
    CA_Sign,
    CA_Re,
    CA_Im,
    CA_Arg,
    CA_Conjugate,
    /* Elementary constants */
    CA_Pi,
    /* Elementary functions */
    CA_Exp,
    CA_Log,
    CA_Pow,
    CA_Cos,
    CA_Sin,
    CA_Tan,
    CA_Cosh,
    CA_Sinh,
    CA_Tanh,
    CA_Atan,
    CA_Acos,
    CA_Asin,
    CA_Atanh,
    CA_Acosh,
    CA_Asinh,
    /* Euler's constant */
    CA_Euler,
    /* Gamma and related functions */
    CA_Gamma,
    CA_LogGamma,
    CA_Psi,
    CA_RiemannZeta,
    CA_HurwitzZeta,
    CA_FUNC_CODE_LENGTH
} calcium_func_code;

const char * calcium_func_name(calcium_func_code func);

#ifdef __cplusplus
}
#endif

#endif

