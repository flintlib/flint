/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "calcium.h"

const char * calcium_func_name(calcium_func_code func)
{
    switch (func)
    {
        /* Arithmetic */
        case CA_Neg: return "Neg";
        case CA_Add: return "Add";
        case CA_Sub: return "Sub";
        case CA_Mul: return "Mul";
        case CA_Div: return "Div";
        /* Roots */
        case CA_Sqrt: return "Sqrt";
        case CA_Cbrt: return "Cbrt";
        case CA_Root: return "Root";
        /* Complex parts */
        case CA_Floor: return "Floor";
        case CA_Ceil: return "Ceil";
        case CA_Abs: return "Abs";
        case CA_Sign: return "Sign";
        case CA_Re: return "Re";
        case CA_Im: return "Im";
        case CA_Arg: return "Arg";
        case CA_Conjugate: return "Conjugate";
        /* Elementary constants */
        case CA_Pi: return "Pi";
        /* Elementary functions */
        case CA_Exp: return "Exp";
        case CA_Log: return "Log";
        case CA_Pow: return "Pow";
        case CA_Cos: return "Cos";
        case CA_Sin: return "Sin";
        case CA_Tan: return "Tan";
        case CA_Cosh: return "Cosh";
        case CA_Sinh: return "Sinh";
        case CA_Tanh: return "Tanh";
        case CA_Atan: return "Atan";
        case CA_Acos: return "Acos";
        case CA_Asin: return "Asin";
        case CA_Atanh: return "Atanh";
        case CA_Acosh: return "Acosh";
        case CA_Asinh: return "Asinh";
        /* Euler's constant */
        case CA_Euler: return "Euler";
        /* Gamma and related functions */
        case CA_Gamma: return "Gamma";
        case CA_LogGamma: return "LogGamma";
        case CA_Psi: return "Psi";
        case CA_Erf: return "Erf";
        case CA_Erfc: return "Erfc";
        case CA_Erfi: return "Erfi";
        case CA_RiemannZeta: return "RiemannZeta";
        case CA_HurwitzZeta: return "HurwitzZeta";
        default: return "<unknown function>";
    }
}

