/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#  if HAVE_ALLOCA_H
#   include <alloca.h>
#  else
#   if _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    ifdef __DECC
#     define alloca(x) __ALLOCA(x)
#    else
#     ifdef BSD
#      include <stdlib.h>
#     else
#      error Could not find alloca
#     endif
#    endif
#   endif
#  endif
# endif
#endif

#include "flint-impl.h"
#include "nmod_mat.h"

void nmod_mat_mul_nmod_vec_ptr(
    mp_limb_t * const * c,
    const nmod_mat_t A,
    const mp_limb_t * const * b, slong blen)
{
    slong i;
    slong len = FLINT_MIN(A->c, blen);
    slong nrows = A->r;
    mp_limb_t * bb, * cc;
    TMP_INIT;

    TMP_START;

    bb = TMP_ARRAY_ALLOC(len, mp_limb_t);
    cc = TMP_ARRAY_ALLOC(nrows, mp_limb_t);

    for (i = 0; i < len; i++)
        bb[i] = b[i][0];

    nmod_mat_mul_nmod_vec(cc, A, bb, len);

    for (i = 0; i < nrows; i++)
        c[i][0] = cc[i];

    TMP_END;
}

