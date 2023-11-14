/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_VEC_H
#define FQ_ZECH_VEC_H

#ifdef FQ_ZECH_VEC_INLINES_C
#define FQ_VEC_TEMPLATES_INLINE
#define FQ_ZECH_VEC_INLINE
#else
#define FQ_VEC_TEMPLATES_INLINE static inline
#define FQ_ZECH_VEC_INLINE static inline
#endif

#include "fq_zech_types.h"

#define FQ_ZECH_VEC_NORM(vec, i, ctx)                    \
do {                                                \
    while ((i) && fq_zech_is_zero((vec) + (i) - 1, ctx)) \
        (i)--;                                      \
} while (0)


#define FQ_ZECH_VEC_SWAP(vec1, len1, vec2, len2)   \
do {                                          \
    fq_zech_struct *__t;                           \
    slong __tn;                               \
    __t    = (vec1);                          \
    (vec1) = (vec2);                          \
    (vec2) = __t;                             \
    __tn   = (len1);                          \
    (len1) = (len2);                          \
    (len2) = __tn;                            \
} while (0);


#define T fq_zech
#define CAP_T FQ_ZECH
#include "fq_vec_templates.h"
#undef CAP_T
#undef T

#endif
