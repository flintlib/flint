/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_DFT_IMPL_H
#define ACB_DFT_IMPL_H

#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void
acb_dft_rad2_reorder(acb_ptr v, slong n);

#ifdef __cplusplus
}
#endif

#endif /* ACB_DFT_IMPL_H */
