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
    fma. The build passes -ffp-contract=off -DDFLOAT_FP_CONTRACT_OFF to
    compilers that take the flag (GCC ignores the standard pragma
    below); MSVC (which contracts under /fp:precise with /arch:AVX2
    before Visual Studio 2022) and clang (also as clang-cl, which may
    not be given the flag) get a pragma. With any other compiler that
    does not take the flag, contraction cannot be ruled out, and
    dfloat_is_supported() returns 0.
*/
#if defined(_MSC_VER) && !defined(__clang__)
# pragma fp_contract(off)
# define DFLOAT_NO_CONTRACTION 1
#elif defined(__clang__)
# pragma STDC FP_CONTRACT OFF
# define DFLOAT_NO_CONTRACTION 1
#elif defined(DFLOAT_FP_CONTRACT_OFF)
# define DFLOAT_NO_CONTRACTION 1
#else
# define DFLOAT_NO_CONTRACTION 0
#endif
