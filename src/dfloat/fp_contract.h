/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Included by every dfloat source file. The error-free transformations
    break if the compiler contracts a separate multiply and add into an
    fma. The build passes -ffp-contract=off to compilers that take it
    (GCC ignores the standard pragma below); MSVC before Visual Studio
    2022 contracts under /fp:precise with /arch:AVX2, and clang-cl may
    not be given the flag, so both get a pragma as well.
*/
#if defined(_MSC_VER) && !defined(__clang__)
# pragma fp_contract(off)
#elif defined(__clang__)
# pragma STDC FP_CONTRACT OFF
#endif
