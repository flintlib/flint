/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PADIC_NMOD_H
#define PADIC_NMOD_H

#ifdef PADIC_NMOD_INLINES_C
#define PADIC_NMOD_INLINE
#else
#define PADIC_NMOD_INLINE static inline
#endif

#include "gr.h"
#include "nmod.h"
#include "padic_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PADIC_EMIN (-(WORD_MAX / 4))
#define PADIC_EMAX ((WORD_MAX / 4))

#define PADIC_NMOD_CTX(ctx) \
((_padic_nmod_ctx_struct *)(GR_CTX_DATA_AS_PTR(ctx)))
#define PADIC_NMOD_CTX_P(ctx) (PADIC_NMOD_CTX(ctx)->p)
#define PADIC_NMOD_CTX_PINV(ctx) (PADIC_NMOD_CTX(ctx)->pinv)
#define PADIC_NMOD_CTX_N(ctx) (PADIC_NMOD_CTX(ctx)->n)
#define PADIC_NMOD_CTX_POW(ctx) (PADIC_NMOD_CTX(ctx)->pow)
#define PADIC_NMOD_CTX_P_MOD(ctx) (PADIC_NMOD_CTX(ctx)->p_mod)
#define PADIC_NMOD_CTX_PN_MOD(ctx) (PADIC_NMOD_CTX(ctx)->pn_mod)

/* Context *******************************************************************/

int padic_nmod_ctx_init(gr_ctx_t ctx, ulong p, slong n);
void padic_nmod_ctx_clear(gr_ctx_t ctx);

/* Memory management *********************************************************/

void padic_nmod_init(padic_nmod_t res, gr_ctx_t ctx);

PADIC_NMOD_INLINE
void _padic_nmod_canonicalise(padic_nmod_t x, gr_ctx_t ctx)
{
    if (x->u != 0)
    {
        if (x->u > 0)
        {
            x->v += n_remove2_precomp(&x->u, PADIC_NMOD_CTX_P(ctx),
                                        PADIC_NMOD_CTX_PINV(ctx));
        }

        else
        {
            ulong z = - x->u;
            slong e = n_remove2_precomp(&z, PADIC_NMOD_CTX_P(ctx),
                                        PADIC_NMOD_CTX_PINV(ctx));

            if (e > 0)
            {
                x->u = - (slong) z;
            }

            x->v += e;
        }
    }

    else
    {
        x->v = 0;
    }
}

/* Randomisation *************************************************************/

int padic_nmod_randtest(padic_nmod_t rop, flint_rand_t state, gr_ctx_t ctx);

int padic_nmod_randtest_not_zero(padic_nmod_t rop, flint_rand_t state,
                                 gr_ctx_t ctx);

/* Assignments and conversions ***********************************************/

int padic_nmod_set(padic_nmod_t res, const padic_nmod_t x, gr_ctx_t ctx);

int padic_nmod_set_ui(padic_nmod_t res, ulong x, gr_ctx_t ctx);

PADIC_NMOD_INLINE
int padic_nmod_zero(padic_nmod_t res, gr_ctx_t ctx)
{
    res->u = 0;
    res->v = PADIC_EMAX;

    return GR_SUCCESS;
}

PADIC_NMOD_INLINE 
int padic_nmod_one(padic_nmod_t res, gr_ctx_t ctx)
{
    res->u = 1;
    res->v = 0;

    return GR_SUCCESS;
}

/* Comparison ****************************************************************/

PADIC_NMOD_INLINE
truth_t padic_nmod_is_zero(const padic_nmod_t x, gr_ctx_t ctx)
{
    return (x->u == 0) ? T_TRUE : T_FALSE;
}

PADIC_NMOD_INLINE
truth_t padic_nmod_is_one(const padic_nmod_t x, gr_ctx_t ctx)
{
    return (x->u == 1 && x->v == 0) ? T_TRUE : T_FALSE;
}

/* Arithmetic operations *****************************************************/

int padic_nmod_add(padic_nmod_t res, const padic_nmod_t a,
                   const padic_nmod_t b, gr_ctx_t ctx);
int padic_nmod_div(padic_nmod_t res, const padic_nmod_t a,
                   const padic_nmod_t b, gr_ctx_t ctx);
int padic_nmod_inv(padic_nmod_t res, const padic_nmod_t a, gr_ctx_t ctx);
int padic_nmod_mul(padic_nmod_t res, const padic_nmod_t a,
                   const padic_nmod_t b, gr_ctx_t ctx);
int padic_nmod_neg(padic_nmod_t res, const padic_nmod_t a, gr_ctx_t ctx);
int padic_nmod_sub(padic_nmod_t res, const padic_nmod_t a,
                   const padic_nmod_t b, gr_ctx_t ctx);

/* Input and output **********************************************************/

void padic_nmod_println(const padic_nmod_t x, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
