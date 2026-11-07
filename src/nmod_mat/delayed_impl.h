/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Reductions of the unreduced two- and three-limb entries of Gaussian
    elimination with delayed reduction (lu_classical_delayed.c,
    ple_classical_delayed.c).
*/

#ifndef NMOD_MAT_DELAYED_IMPL_H
#define NMOD_MAT_DELAYED_IMPL_H

#include "longlong.h"
#include "mpn_extras.h"
#include "nmod.h"

FLINT_FORCE_INLINE
void mullo_2x1(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b0)
{
    ulong t0, t1;
    umul_ppmm(t1, t0, a0, b0);
    t1 += a1 * b0;
    *r0 = t0;
    *r1 = t1;
}

/* two words, n < 2^(FLINT_BITS - 1): n_ll_barrett_precomp_preinv and
   n_ll_mod_barrett of ulong_extras */

FLINT_FORCE_INLINE ulong
n_lll_rem_l_fullword_limited(ulong y2, ulong y1, ulong y0, nmod_t mod, ulong alpha2, ulong alpha1)
{
    ulong c1, c0, t1, t0;
    ulong xhi, xlo;

    FLINT_ASSERT(mod.n >= (UWORD(1) << (FLINT_BITS - 1)));
    FLINT_ASSERT(mod.n < (UWORD(1) << (FLINT_BITS - 1)) + (UWORD(1) << (FLINT_BITS / 2 - 2)));

    umul_ppmm(t1, t0, y2, alpha2);
    umul_ppmm(c1, c0, y1, alpha1);
    add_ssaaaa(xhi, xlo, t1, t0, c1, c0);
    add_ssaaaa(xhi, xlo, xhi, xlo, 0, y0);

    NMOD_RED2_FULLWORD(xlo, xhi, xlo, mod);

    return xlo;
}

FLINT_FORCE_INLINE ulong
n_lll_rem_l(ulong y2, ulong y1, ulong y0, nmod_t mod, ulong alpha2, ulong alpha1)
{
    ulong c1, c0, t1, t0;
    ulong xhi, xlo;

    umul_ppmm(t1, t0, y2, alpha2);
    umul_ppmm(c1, c0, y1, alpha1);
    add_ssaaaa(xhi, xlo, t1, t0, c1, c0);
    add_ssaaaa(xhi, xlo, xhi, xlo, 0, y0);

    if (xhi >= mod.n) xhi -= mod.n;
    NMOD_RED2(xlo, xhi, xlo, mod);

    return xlo;
}

#endif
