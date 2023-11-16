/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "dlog.h"

void
dlog_rho_init(dlog_rho_t t, ulong a, ulong mod, ulong n)
{
    t->a = a;
    nmod_init(&t->n, n);
    nmod_init(&t->mod, mod);
    t->nisprime = n_is_prime(n);
}
