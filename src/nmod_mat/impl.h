/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MAT_IMPL_H
#define NMOD_MAT_IMPL_H

#include "flint.h"

/*
    Whether nmod_mat_mul_u52 has a kernel: it needs AVX512-IFMA (with F and
    DQ) at compile time and 64-bit words; otherwise it declines every
    multiplication. Both mul_u52.c and the dispatch in mul.c test this.
*/
#if FLINT_BITS == 64 && defined(__AVX512F__) && defined(__AVX512DQ__) \
        && defined(__AVX512IFMA__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_MAT_HAVE_MUL_U52 1
#else
# define NMOD_MAT_HAVE_MUL_U52 0
#endif

#endif
