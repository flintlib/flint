/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef BERNOULLI_IMPL_H
#define BERNOULLI_IMPL_H

#include "fmpz_types.h"
#include "bernoulli.h"

ulong _bernoulli_n_muldivrem_precomp(ulong * q, ulong a, ulong b, ulong n, double bnpre);
ulong _bernoulli_mod_p_harvey_pow2(ulong p, ulong pinv, ulong k);
ulong _bernoulli_mod_p_harvey_powg(ulong p, ulong pinv, ulong k);

/* bernoulli_rev_init with test flags: the generic table kernel even when
   a vector kernel is available, and j^2 and n (n - 1) applied as two
   factors (as on 32-bit machines for large n) */
#define BERNOULLI_REV_TEST_GENERIC 1
#define BERNOULLI_REV_TEST_SPLIT 2
/* the sieve of the iterator (by default by {3, 5} for large n) */
#define BERNOULLI_REV_SIEVE_3 4
#define BERNOULLI_REV_SIEVE_15 8
#define BERNOULLI_REV_NO_SIEVE 16

void _bernoulli_rev_init_flags(bernoulli_rev_t iter, ulong nmax, int flags);

#endif
