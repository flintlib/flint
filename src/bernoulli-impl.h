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

#ifdef __cplusplus
}
#endif

#endif /* BERNOULLI_IMPL_H */
