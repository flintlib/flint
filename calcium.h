/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef CALCIUM_H
#define CALCIUM_H

#ifdef CALCIUM_INLINES_C
#define CALCIUM_INLINE
#else
#define CALCIUM_INLINE static __inline__
#endif

#include <stdio.h>
#include "flint/flint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    slong dummy;
}
ca_ctx_struct;

typedef ca_ctx_struct ca_ctx_t[1];

void ca_ctx_init(ca_ctx_t ctx);

void ca_ctx_clear(ca_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

