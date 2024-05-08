/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_TYPES_H
#define CA_TYPES_H

#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* streams *******************************************************************/

/* TODO: Deprecate these and simply replace with gr_stream_struct */
#define calcium_stream_struct gr_stream_struct
#define calcium_stream_t gr_stream_t

/* builtin functions and constants *******************************************/

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
    CA_Sin,
    CA_Cos,
    CA_Exp,
    CA_Log,
    CA_Pow,
    CA_Tan,
    CA_Cot,
    CA_Cosh,
    CA_Sinh,
    CA_Tanh,
    CA_Coth,
    CA_Atan,
    CA_Acos,
    CA_Asin,
    CA_Acot,
    CA_Atanh,
    CA_Acosh,
    CA_Asinh,
    CA_Acoth,
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
}
calcium_func_code;

#ifdef __cplusplus
}
#endif

#endif /* CA_TYPES_H */
