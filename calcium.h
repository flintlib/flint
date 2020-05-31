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

/* todo: IDs for builtin mathematical functions and constants */

#define CA_Neg   10
#define CA_Add   11
#define CA_Sub   12
#define CA_Mul   13
#define CA_Div   14
#define CA_Sqrt  15
#define CA_Root  16

#define CA_Pi    50

#define CA_Exp   100
#define CA_Log   101
#define CA_Pow   102
#define CA_Cos   103
#define CA_Sin   104
#define CA_Tan   105
#define CA_Atan  106
#define CA_Acos  107
#define CA_Asin  108

#ifdef __cplusplus
}
#endif

#endif

