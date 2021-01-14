/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

const fexpr_symbol_info fexpr_builtins[FEXPR_BUILTIN_LENGTH] = {
    { FEXPR_Abs, "Abs" },
    { FEXPR_Acos, "Acos" },
    { FEXPR_Acosh, "Acosh" },
    { FEXPR_Add, "Add" },
    { FEXPR_AiryAi, "AiryAi" },
    { FEXPR_AiryBi, "AiryBi" },
    { FEXPR_Arg, "Arg" },
    { FEXPR_Asin, "Asin" },
    { FEXPR_Asinh, "Asinh" },
    { FEXPR_Atan, "Atan" },
    { FEXPR_Atanh, "Atanh" },
    { FEXPR_BesselI, "BesselI" },
    { FEXPR_BesselJ, "BesselJ" },
    { FEXPR_BesselK, "BesselK" },
    { FEXPR_BesselY, "BesselY" },
    { FEXPR_Ceil, "Ceil" },
    { FEXPR_Conjugate, "Conjugate" },
    { FEXPR_Cos, "Cos" },
    { FEXPR_Cosh, "Cosh" },
    { FEXPR_Div, "Div" },
    { FEXPR_Erf, "Erf" },
    { FEXPR_Erfc, "Erfc" },
    { FEXPR_Erfi, "Erfi" },
    { FEXPR_Euler, "Euler" },
    { FEXPR_Exp, "Exp" },
    { FEXPR_Floor, "Floor" },
    { FEXPR_Gamma, "Gamma" },
    { FEXPR_HurwitzZeta, "HurwitzZeta" },
    { FEXPR_I, "I" },
    { FEXPR_Im, "Im" },
    { FEXPR_JacobiTheta, "JacobiTheta" },
    { FEXPR_LambertW, "LambertW" },
    { FEXPR_Log, "Log" },
    { FEXPR_LogGamma, "LogGamma" },
    { FEXPR_Mul, "Mul" },
    { FEXPR_Neg, "Neg" },
    { FEXPR_Pi, "Pi" },
    { FEXPR_Pos, "Pos" },
    { FEXPR_Pow, "Pow" },
    { FEXPR_Psi, "Psi" },
    { FEXPR_Re, "Re" },
    { FEXPR_RiemannZeta, "RiemannZeta" },
    { FEXPR_Root, "Root" },
    { FEXPR_RootOfUnity, "RootOfUnity" },
    { FEXPR_Sign, "Sign" },
    { FEXPR_Sin, "Sin" },
    { FEXPR_Sinh, "Sinh" },
    { FEXPR_Sqrt, "Sqrt" },
    { FEXPR_Sub, "Sub" },
    { FEXPR_Tan, "Tan" },
    { FEXPR_Tanh, "Tanh" },
};

slong
fexpr_get_builtin_str(const char * s)
{
    slong a, mid, b;
    int cmp;

    a = 0;
    b = FEXPR_BUILTIN_LENGTH - 1;

    while (a <= b)
    {
        mid = (a + b) / 2;
        cmp = strcmp(fexpr_builtins[mid].string, s);

        if (cmp == 0)
            return mid;
        else if (cmp > 0)
            b = mid - 1;
        else
            a = mid + 1;
    }

    return -1;
}
