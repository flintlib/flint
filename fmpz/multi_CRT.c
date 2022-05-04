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
#include "fmpz.h"

int fmpz_multi_CRT(
    fmpz_t output,
    const fmpz * moduli,
    const fmpz * values,
    slong len,
    int sign)
{
    int success;
    slong i;
    fmpz_multi_crt_t P;
    fmpz * out;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    fmpz_multi_CRT_init(P);
    success = fmpz_multi_CRT_precompute(P, moduli, len);

    out = TMP_ARRAY_ALLOC(P->localsize, fmpz);
    for (i = 0; i < P->localsize; i++)
        fmpz_init(out + i);

    fmpz_swap(out + 0, output);
    _fmpz_multi_CRT_precomp(out, P, values, sign);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
        fmpz_clear(out + i);

    fmpz_multi_CRT_clear(P);

    TMP_END;

    return success;
}

