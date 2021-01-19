/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FEXPR_EVAL_H
#define FEXPR_EVAL_H

#ifdef FEXPR_EVAL_INLINES_C
#define FEXPR_EVAL_INLINE
#else
#define FEXPR_EVAL_INLINE static __inline__
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "fexpr.h"

#define FEXPR_EVAL_SUCCESS                0
#define FEXPR_EVAL_ERR_NOT_IMPLEMENTED    1
#define FEXPR_EVAL_ERR_SYNTAX             2
#define FEXPR_EVAL_ERR_SYMBOL             3
#define FEXPR_EVAL_ERR_TYPE               4
#define FEXPR_EVAL_ERR_OVERFLOW           5
#define FEXPR_EVAL_ERR_PRECISION          6

typedef struct
{
    fexpr_vec_struct bound_variables;         /* variables bound to values */
    fexpr_vec_struct bound_values;            /* the corresponding values */
    fexpr_vec_struct assumptions;             /* assumptions on free variables */
    fexpr_vec_struct assumptions_variables;   /* free variables appearing in assumptions */
}
fexpr_eval_struct;

typedef fexpr_eval_struct fexpr_eval_t[1];

void fexpr_eval_init(fexpr_eval_t eval);
void fexpr_eval_clear(fexpr_eval_t eval);

void fexpr_eval_push_def(fexpr_eval_t eval, const fexpr_t symbol, const fexpr_t value);
void fexpr_eval_pop_def(fexpr_eval_t eval);

void fexpr_eval_push_assumptions(fexpr_eval_t eval, const fexpr_t assumptions);
void fexpr_eval_pop_assumptions(fexpr_eval_t eval);

#ifdef __cplusplus
}
#endif

#endif

