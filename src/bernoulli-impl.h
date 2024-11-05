/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef BERNOULLI_IMPL_H
#define BERNOULLI_IMPL_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

void
_arb_tree_crt(fmpz_t r, fmpz_t m, nn_srcptr residues, nn_srcptr primes, slong len);

ulong _bernoulli_n_muldivrem_precomp(ulong * q, ulong a, ulong b, ulong n, double bnpre);

ulong _bernoulli_mod_p_harvey_powg(ulong p, ulong pinv, ulong k);
ulong _bernoulli_mod_p_harvey_pow2(ulong p, ulong pinv, ulong k);

#ifdef __cplusplus
}
#endif

#endif /* BERNOULLI_IMPL_H */
