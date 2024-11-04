/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MPOLY_IMPL_H
#define FQ_NMOD_MPOLY_IMPL_H

#include "fq_nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void fq_nmod_poly_product_roots(fq_nmod_poly_t P, fq_nmod_struct * r,
                                            slong n, const fq_nmod_ctx_t fqctx);

#ifdef __cplusplus
}
#endif

#endif /* FQ_NMOD_MPOLY_IMPL_H */
