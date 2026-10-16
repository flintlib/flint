/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PARTITIONS_IMPL_H
#define PARTITIONS_IMPL_H

#include "flint.h"

slong partitions_hrr_needed_terms(double n);

/* p(n) for n < 128 */
extern const unsigned int partitions_lookup[128];

/* v[n] = p(n) mod 2^FLINT_BITS for 0 <= n < len, by the pentagonal
   recurrence (exact while p(n) < 2^FLINT_BITS) */
void _partitions_vec_ui(nn_ptr v, slong len);

#endif
