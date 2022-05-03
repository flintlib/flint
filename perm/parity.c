/*
    Copyright (C) 2011 Fredrik Johansson

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

#include "flint.h"
#include "flint-impl.h"

int
_perm_parity(const slong *vec, slong n)
{
    slong i, k;
    int * encountered;
    int parity;
    TMP_INIT;

    if (n <= 1)
        return 0;

    TMP_START;

    parity = 0;
    encountered = (int *) TMP_ALLOC(n*sizeof(int));

    for (i = 0; i < n; i++)
       encountered[i] = 0;

    for (i = 0; i < n; i++)
    {
        if (encountered[i] != 0)
        {
            parity ^= 1;
        }
        else
        {
            k = i;
            do
            {
                k = vec[k];
                encountered[k] = 1;
            }
            while (k != i);
        }
    }

    TMP_END;

    return parity;
}
