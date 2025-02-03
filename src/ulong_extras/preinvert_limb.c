/*
    Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
    2004, 2005 Free Software Foundation, Inc.

    Copyright 2009, 2015 William Hart
    Copyright 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong n_preinvert_limb(ulong n)
{
    ulong ninv, r;
    unsigned int norm;

    norm = flint_clz(n);
    n <<= norm;

    udiv_qrnnd(ninv, r, ~n, ~UWORD(0), n);
    return ninv;
}

ulong n_preinvert_limb_prenorm(ulong n)
{
    ulong ninv, r;

    udiv_qrnnd(ninv, r, ~n, ~UWORD(0), n);
    return ninv;
}
