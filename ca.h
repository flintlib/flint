/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef CA_H
#define CA_H

#ifdef CA_INLINES_C
#define CA_INLINE
#else
#define CA_INLINE static __inline__
#endif

#include <stdio.h>
#include "flint/flint.h"
#include "nf.h"
#include "nf_elem.h"
#include "calcium.h"
#include "qqbar.h"
#include "fmpz_mpoly_q.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Number object *************************************************************/

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

#define CA_FMPQ(x)         (&((x)->elem.q))
#define CA_MPOLY_Q(x)      (((x)->elem->mpoly_q))
#define CA_NF_ELEM(x)      (&((x)->elem.nf))
#define CA_FMPQ_NUMREF(x)  (fmpq_numref(CA_FMPQ(x)))
#define CA_FMPQ_DENREF(x)  (fmpq_denref(CA_FMPQ(x)))

/* Extension object **********************************************************/

/* There are currently two kinds of extension elements: algebraic numbers,
   and symbolic functions. */
typedef enum
{
    CA_EXT_QQBAR,
    CA_EXT_FUNCTION
}
ca_extension_type;

typedef struct
{
    qqbar_struct x;     /* qqbar element */
    nf_struct nf;       /* antic number field for fast arithmetic */
}
ca_extension_data_qqbar;

typedef struct
{
    ulong func;             /* f = F_Pi, F_Exp, ... */
    slong num_args;         /* n */
    ca_struct * args;       /* x1, ..., xn */
    acb_struct enclosure;   /* Numerical enclosure of f(x1,...,xn) */
}
ca_extension_data_function;

typedef union
{
    ca_extension_data_qqbar qqbar;
    ca_extension_data_function function;
}
ca_extension_data_struct;

typedef struct
{
    ca_extension_type type;
    char * string;
    ca_extension_data_struct data;
}
ca_extension_struct;

typedef ca_extension_struct ca_extension_t[1];

/* Field object **************************************************************/

typedef enum
{
    CA_FIELD_QQ,       /* field elements are represented as fmpq_t */
    CA_FIELD_NF,       /* field elements are represented as nf_elem_t */
    CA_FIELD_MPOLY_Q   /* field elements are represented as fmpz_mpoly_q_t */
}
ca_field_type_t;

typedef struct
{
    fmpz_mpoly_ctx_struct mctx;  /* todo: should perhaps be a reference to a fixed table of precomputed contexts */
    ca_field_type_t type;
    ca_extension_struct ** ext;
    slong len;
    fmpz_mpoly_struct * ideal;
    slong ideal_len;
}
ca_field_struct;

typedef ca_field_struct ca_field_t[1];

#define CA_FIELD_MCTX(K) (&((K)->mctx))

/* Context object ************************************************************/

typedef struct
{
    slong dummy;
}
ca_ctx_struct;

typedef ca_ctx_struct ca_ctx_t[1];

/* Memory management */

void ca_ctx_init(ca_ctx_t ctx);

void ca_ctx_clear(ca_ctx_t ctx);

void ca_init(ca_t x, ca_ctx_t ctx);

void ca_clear(ca_t x, ca_ctx_t ctx);

/* Extension and field methods */

void ca_extension_init_qqbar(ca_extension_t ext, const qqbar_t x);

void ca_extension_init_const(ca_extension_t ext, ulong func);

void ca_extension_init_fx(ca_extension_t ext, ulong func, const ca_t x);

void ca_extension_clear(ca_extension_t ext);

void ca_extension_print(const ca_extension_t ext);

void ca_field_init_qq(ca_field_t K);

void ca_field_init_mpoly_q(ca_field_t K, slong len);

void ca_field_clear(ca_field_t K);

void ca_field_set_ext(ca_field_t K, slong i, ca_extension_struct * ext);

void ca_field_print(const ca_field_t K);

#ifdef __cplusplus
}
#endif

#endif

