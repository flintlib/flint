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
#include "flint/fmpz.h"
#include "acb.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Global library definitions */

const char * calcium_version(void);

#define __CALCIUM_VERSION 0
#define __CALCIUM_VERSION_MINOR 2
#define __CALCIUM_VERSION_PATCHLEVEL 0
#define CALCIUM_VERSION "0.2.0"
#define __CALCIUM_RELEASE (__CALCIUM_VERSION * 10000 + \
                         __CALCIUM_VERSION_MINOR * 100 + \
                         __CALCIUM_VERSION_PATCHLEVEL)

double calcium_test_multiplier(void);

/* Input and output */

typedef struct
{
    FILE * fp;
    char * s;
    slong len;
    slong alloc;
}
calcium_stream_struct;

typedef calcium_stream_struct calcium_stream_t[1];

CALCIUM_INLINE
void calcium_stream_init_file(calcium_stream_t out, FILE * fp)
{
    out->fp = fp;
    out->s = NULL;
}

CALCIUM_INLINE
void calcium_stream_init_str(calcium_stream_t out)
{
    out->fp = NULL;
    out->s = NULL;
    out->len = 0;
    out->alloc = 0;
}

void calcium_write(calcium_stream_t out, const char * s);
void calcium_write_si(calcium_stream_t out, slong x);
void calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags);

CALCIUM_INLINE
void calcium_write_free(calcium_stream_t out, char * s)
{
    calcium_write(out, s);
    flint_free(s);
}

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
    CA_Erf,
    CA_Erfc,
    CA_Erfi,
    CA_RiemannZeta,
    CA_HurwitzZeta,
    CA_FUNC_CODE_LENGTH
} calcium_func_code;

const char * calcium_func_name(calcium_func_code func);

/* Flint extras */

/* slower alternative: fmpz_fdiv_ui(x 1000000007) */
CALCIUM_INLINE ulong calcium_fmpz_hash(const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
        return *x;
    else
    {
        __mpz_struct * z = COEFF_TO_PTR(*x);
        return (z->_mp_size > 0) ? z->_mp_d[0] : -z->_mp_d[0];
    }
}

#ifdef __cplusplus
}
#endif

#endif

