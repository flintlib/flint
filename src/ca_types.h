/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_TYPES_H
#define CA_TYPES_H

#include "fmpz_types.h"
#include "mpoly_types.h"
#include "acb_types.h"
#include "gr_types.h"
#include "nf_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* macros ********************************************************************/

/* Use the low two bits of the field pointer to encode special values. */
/* The field pointer with the mask removed is NULL for
   Unknown/Undefined/Uinf, and a normal field pointer for signed
   infinity (encoding the sign). */
#define CA_UNKNOWN        UWORD(1)
#define CA_UNDEFINED      UWORD(2)
#define CA_INF            UWORD(3)
#define CA_SPECIAL        (CA_UNKNOWN | CA_UNDEFINED | CA_INF)

/* streams *******************************************************************/

/* TODO: Deprecate these and simply replace with gr_stream_struct */
#define calcium_stream_struct gr_stream_struct
#define calcium_stream_t gr_stream_t

/* builtin functions and constants *******************************************/

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
}
calcium_func_code;

/* symbolic expressions ******************************************************/

typedef struct
{
    ulong * data;
    slong alloc;
}
fexpr_struct;

typedef fexpr_struct fexpr_t[1];
typedef fexpr_struct * fexpr_ptr;
typedef const fexpr_struct * fexpr_srcptr;

typedef struct
{
    fexpr_struct * entries;
    slong alloc;
    slong length;
}
fexpr_vec_struct;

typedef fexpr_vec_struct fexpr_vec_t[1];

#define fexpr_vec_entry(vec, i) ((vec)->entries + (i))

/* numbers objects ***********************************************************/

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

#define CA_FMPQ(x)         (&((x)->elem.q))
#define CA_MPOLY_Q(x)      (&(((x)->elem.mpoly_q)[0]))
#define CA_NF_ELEM(x)      (&((x)->elem.nf))
#define CA_FMPQ_NUMREF(x)  (fmpq_numref(CA_FMPQ(x)))
#define CA_FMPQ_DENREF(x)  (fmpq_denref(CA_FMPQ(x)))

#define CA_FIELD(x, ctx)     ((ca_field_ptr) ((x)->field))
#define CA_FIELD_ULONG(x)    ((x)->field)

#define CA_IS_SPECIAL(x)       (CA_FIELD_ULONG(x) & CA_SPECIAL)
#define CA_IS_UNKNOWN(x)       (CA_FIELD_ULONG(x) == CA_UNKNOWN)
#define CA_IS_UNDEFINED(x)     (CA_FIELD_ULONG(x) == CA_UNDEFINED)
#define CA_IS_INF(x)           ((CA_FIELD_ULONG(x) & CA_SPECIAL) == CA_INF)
#define CA_IS_UNSIGNED_INF(x)  (CA_FIELD_ULONG(x) == CA_INF)
#define CA_IS_SIGNED_INF(x)    (CA_IS_INF(x) && !CA_IS_UNSIGNED_INF(x))

/* We always allocate QQ and QQ(i) */
#define CA_IS_QQ(x, ctx) (CA_FIELD(x, ctx) == (ctx)->field_qq)
#define CA_IS_QQ_I(x, ctx) (CA_FIELD(x, ctx) == (ctx)->field_qq_i)

#define CA_FIELD_UNSPECIAL(x, ctx) ((ca_field_ptr) (CA_FIELD_ULONG(x) & ~CA_SPECIAL))

/* algebraic numbers via minimal polynomials *********************************/

typedef struct
{
    fmpz_poly_struct poly;
    acb_struct enclosure;
}
qqbar_struct;

typedef qqbar_struct qqbar_t[1];
typedef qqbar_struct * qqbar_ptr;
typedef const qqbar_struct * qqbar_srcptr;

/* extension objects *********************************************************/

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

#define CA_EXT_HEAD(x) ((x)->head)
#define CA_EXT_HASH(x) ((x)->hash)
#define CA_EXT_DEPTH(x) ((x)->depth)

#define CA_EXT_IS_QQBAR(x) ((x)->head == CA_QQBar)

#define CA_EXT_QQBAR(_x) (&((_x)->data.qqbar.x))
#define CA_EXT_QQBAR_NF(_x) ((_x)->data.qqbar.nf)

#define CA_EXT_FUNC_ARGS(x) ((x)->data.func_data.args)
#define CA_EXT_FUNC_NARGS(x) ((x)->data.func_data.nargs)
#define CA_EXT_FUNC_ENCLOSURE(x) (&((x)->data.func_data.enclosure))
#define CA_EXT_FUNC_PREC(x) ((x)->data.func_data.prec)

/* field objects *************************************************************/

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

#define CA_FIELD_LENGTH(K) ((K)->length)
#define CA_FIELD_EXT(K) ((K)->ext)
#define CA_FIELD_EXT_ELEM(K, i) ((K)->ext[i])
#define CA_FIELD_HASH(K) ((K)->hash)

#define CA_FIELD_IS_QQ(K) ((K)->length == 0)
#define CA_FIELD_IS_NF(K) ((K)->ideal.length == -1)
#define CA_FIELD_IS_GENERIC(K) (!CA_FIELD_IS_QQ(K) && !CA_FIELD_IS_NF(K))

#define CA_FIELD_NF(K) (((K)->ext[0]->data.qqbar.nf))
#define CA_FIELD_NF_QQBAR(K) (&((K)->ext[0]->data.qqbar.x))

#define CA_FIELD_IDEAL(K) (&((K)->ideal))
#define CA_FIELD_IDEAL_ELEM(K, i) fmpz_mpoly_vec_entry(CA_FIELD_IDEAL(K), i)
#define CA_FIELD_IDEAL_LENGTH(K) ((K)->ideal.length)
#define CA_FIELD_IDEAL_ALLOC(K) ((K)->ideal.alloc)
#define CA_FIELD_IDEAL_P(K) ((K)->ideal.p)

#define CA_MCTX_1(ctx) ((ctx)->mctx[0])
#define CA_FIELD_MCTX(K, ctx) ((ctx)->mctx[CA_FIELD_LENGTH(K) - 1])

/* context objects ***********************************************************/

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

#ifdef __cplusplus
}
#endif

#endif /* CA_TYPES_H */
