/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_MODULAR_IMPL_H
#define ACB_MODULAR_IMPL_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

slong
acb_modular_rs_optimal_m(const int * best_ms, const int * num_residues, slong N);

#ifdef __cplusplus
}
#endif

#endif /* ACB_MODULAR_IMPL_H */
