/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_SERIES_H
#define GR_SERIES_H

#include "flint.h"
#include "gr.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    gr_ctx_struct * base_ring;
    slong n;
    char * var;
}
series_mod_ctx_t;

typedef struct
{
    gr_ctx_struct * base_ring;
    slong prec;     /* default approximate truncation */
    char * var;
}
series_ctx_t;

#define GR_SERIES_CTX(ring_ctx) ((series_ctx_t *)((ring_ctx)))
#define GR_SERIES_ELEM_CTX(ring_ctx) (GR_SERIES_CTX(ring_ctx)->base_ring)
#define GR_SERIES_PREC(ring_ctx) (GR_SERIES_CTX(ring_ctx)->prec)

#define GR_SERIES_MOD_CTX(ring_ctx) ((series_mod_ctx_t *)((ring_ctx)))
#define GR_SERIES_MOD_ELEM_CTX(ring_ctx) (GR_SERIES_MOD_CTX(ring_ctx)->base_ring)
#define GR_SERIES_MOD_N(ring_ctx) (GR_SERIES_MOD_CTX(ring_ctx)->n)


#define GR_SERIES_ERR_EXACT WORD_MAX
#define GR_SERIES_ERR_MAX WORD_MAX / 4

typedef struct
{
    gr_poly_struct poly;
    slong error;
}
gr_series_struct;

typedef gr_series_struct gr_series_t[1];

#define GR_SERIES_POLY(x) (&((x)->poly))
#define GR_SERIES_ERROR(x) ((x)->error)



#ifdef __cplusplus
}
#endif

#endif /* GR_SERIES_H */
