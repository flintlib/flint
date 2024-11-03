/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_MPOLY_IMPL_H
#define FQ_ZECH_MPOLY_IMPL_H

#include "fq_zech_mpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

fq_zech_mpoly_struct * _fq_zech_mpolyu_get_coeff(fq_zech_mpolyu_t A,
                                     ulong pow, const fq_zech_mpoly_ctx_t uctx);

#ifdef __cplusplus
}
#endif

#endif /* FQ_ZECH_MPOLY_IMPL_H */
