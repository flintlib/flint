/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ULONG_EXTRAS_MINI_H
#define ULONG_EXTRAS_MINI_H

#ifdef ULONG_EXTRAS_MINI_INLINES_C
#define ULONG_EXTRAS_MINI_INLINE FLINT_DLL
#else
#define ULONG_EXTRAS_MINI_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL ulong n_sqrt(ulong a);

FLINT_DLL ulong n_gcd(ulong x, ulong y);
FLINT_DLL ulong n_gcdinv(ulong_ptr a, ulong x, ulong y);

ULONG_EXTRAS_MINI_INLINE
ulong n_invmod(ulong x, ulong y)
{
   ulong r, g;

   g = n_gcdinv(&r, x, y);
   if (g != 1)
      flint_throw(FLINT_IMPINV, "Cannot invert modulo in n_invmod\n");

   return r;
}

ULONG_EXTRAS_MINI_INLINE
ulong n_submod(ulong x, ulong y, ulong n)
{
    FLINT_ASSERT(x < n);
    FLINT_ASSERT(y < n);
    FLINT_ASSERT(n != 0);

    return (y > x ? x - y + n : x - y);
}

ULONG_EXTRAS_MINI_INLINE
ulong n_preinvert_limb(ulong n)
{
   ulong norm, ninv;

   count_leading_zeros(norm, n);
   invert_limb(ninv, n << norm);

   return ninv;
}

FLINT_DLL ulong n_ll_mod_preinv(ulong a_hi, ulong a_lo, ulong n, ulong ninv);

FLINT_DLL ulong n_mod2_preinv(ulong a, ulong n, ulong ninv);

ULONG_EXTRAS_MINI_INLINE
ulong n_mulmod2_preinv(ulong a, ulong b, ulong n, ulong ninv)
{
    ulong p1, p2;

    FLINT_ASSERT(n != 0);

    umul_ppmm(p1, p2, a, b);
    return n_ll_mod_preinv(p1, p2, n, ninv);
}

FLINT_DLL ulong n_powmod2_ui_preinv(ulong a, ulong exp, ulong n, ulong ninv);

#ifdef __cplusplus
}
#endif

#endif
