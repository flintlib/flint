/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_is_probabprime_fermat(mp_limb_t n, mp_limb_t i)
{
    if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
        return (n_powmod(i, n - 1, n) == UWORD(1));
    else
        return n_powmod2_ui_preinv(i, n - 1, n, n_preinvert_limb(n)) == UWORD(1);
}
