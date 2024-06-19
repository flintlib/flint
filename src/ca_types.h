/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_TYPES_H
#define CA_TYPES_H

#include "gr_types.h"
#include "fmpz_types.h"
#include "fmpz_mpoly.h" /* FIXME: Move fmpz_mpoly_vec_struct to fmpz_types? */
#include "qqbar.h" /* FIXME: Would be nice to remove */
#include "nf_elem.h" /* FIXME: Would be nice to remove */

#ifdef __cplusplus
extern "C" {
#endif

/* I/O ***********************************************************************/

#ifndef calcium_stream_struct
# define calcium_stream_struct gr_stream_struct
# define calcium_stream_t gr_stream_t
#endif

/* IDs for builtin mathematical functions and constants **********************/

typedef enum
{
    /* Special case for representing qqbar instances */
    CA_QQBar,
    /* Arithmetic */
    CA_Neg,
    CA_Add,
    CA_Sub,
    CA_Mul,
    CA_Div,
    /* Roots */
    CA_Sqrt,
    CA_Cbrt,
    CA_Root,
    /* Complex parts */
    CA_Floor,
    CA_Ceil,
    CA_Abs,
    CA_Sign,
    CA_Re,
    CA_Im,
    CA_Arg,
    CA_Conjugate,
    /* Elementary constants */
    CA_Pi,
    /* Elementary functions */
    CA_Sin,
    CA_Cos,
    CA_Exp,
    CA_Log,
    CA_Pow,
    CA_Tan,
    CA_Cot,
    CA_Cosh,
    CA_Sinh,
    CA_Tanh,
    CA_Coth,
    CA_Atan,
    CA_Acos,
    CA_Asin,
    CA_Acot,
    CA_Atanh,
    CA_Acosh,
    CA_Asinh,
    CA_Acoth,
    /* Euler's constant */
    CA_Euler,
    /* Gamma and related functions */
    CA_Gamma,
    CA_LogGamma,
    CA_Psi,
    CA_Erf,
    CA_Erfc,
    CA_Erfi,
    CA_RiemannZeta,
    CA_HurwitzZeta,
    CA_FUNC_CODE_LENGTH
} calcium_func_code;

/* element object ************************************************************/

typedef union
{
    fmpq q;                           /* rational number */
    nf_elem_struct nf;                /* algebraic number field element */
    fmpz_mpoly_q_struct * mpoly_q;    /* generic field element */
}
ca_elem_struct;

typedef struct
{
    ulong field;
    ca_elem_struct elem;
}
ca_struct;

typedef ca_struct ca_t[1];
typedef ca_struct * ca_ptr;
typedef const ca_struct * ca_srcptr;

/* extension object **********************************************************/

typedef struct
{
    qqbar_struct x;        /* qqbar_t element */
    nf_struct * nf;        /* antic number field for fast arithmetic */
}
ca_ext_qqbar;

typedef struct
{
    ca_struct * args;       /* Function arguments x1, ..., xn. */
    slong nargs;            /* Number of function arguments n. */
    acb_struct enclosure;   /* Cached numerical enclosure of f(x1,...,xn) */
    slong prec;             /* Working precision of cached enclosure */
    qqbar_struct * qqbar;   /* Cached qqbar */
}
ca_ext_func_data;

typedef struct
{
    calcium_func_code head;   /* f = F_Pi, F_Exp, ... */
    ulong hash;
    slong depth;
    union {
        ca_ext_qqbar qqbar;
        ca_ext_func_data func_data;
    } data;
} ca_ext_struct;

typedef ca_ext_struct ca_ext_t[1];
typedef ca_ext_struct * ca_ext_ptr;
typedef const ca_ext_struct * ca_ext_srcptr;

typedef struct
{
    ca_ext_struct ** items;
    slong length;
    slong alloc;

    slong hash_size;
    slong * hash_table;
}
ca_ext_cache_struct;

typedef ca_ext_cache_struct ca_ext_cache_t[1];

/* field object **************************************************************/

typedef struct
{
    slong length;                /* Number of generators              */
    ca_ext_struct ** ext;        /* Generators                        */
    fmpz_mpoly_vec_struct ideal; /* Algebraic relations for reduction */
    ulong hash;
}
ca_field_struct;

typedef ca_field_struct ca_field_t[1];
typedef ca_field_struct * ca_field_ptr;
typedef const ca_field_struct * ca_field_srcptr;

typedef struct
{
    ca_field_struct ** items;
    slong length;
    slong alloc;

    slong hash_size;
    slong * hash_table;
}
ca_field_cache_struct;

typedef ca_field_cache_struct ca_field_cache_t[1];

/* context object ************************************************************/

typedef struct
{
    ca_ext_cache_struct ext_cache;              /* Cached extension objects */
    ca_field_cache_struct field_cache;          /* Cached extension fields  */
    ca_field_struct * field_qq;                 /* Quick access to QQ      */
    ca_field_struct * field_qq_i;               /* Quick access to QQ(i)   */
    fmpz_mpoly_ctx_struct ** mctx;              /* Cached contexts for multivariate polys */
    slong mctx_len;
    slong * options;
}
ca_ctx_struct;

typedef ca_ctx_struct ca_ctx_t[1];

/* factor object *************************************************************/

typedef struct
{
    ca_ptr base;
    ca_ptr exp;
    slong length;
    slong alloc;
}
ca_factor_struct;

typedef ca_factor_struct ca_factor_t[1];

/* vector object *************************************************************/

typedef struct
{
    ca_ptr entries;
    slong alloc;
    slong length;
}
ca_vec_struct;

typedef ca_vec_struct ca_vec_t[1];

/* matrix object *************************************************************/

typedef struct
{
    ca_ptr entries;
    slong r;
    slong c;
    ca_ptr * rows;
}
ca_mat_struct;

typedef ca_mat_struct ca_mat_t[1];

/* polynomial object *********************************************************/

typedef struct
{
    ca_struct * coeffs;
    slong alloc;
    slong length;
}
ca_poly_struct;

typedef ca_poly_struct ca_poly_t[1];

#ifdef __cplusplus
}
#endif

#endif /* CA_TYPES_H */
