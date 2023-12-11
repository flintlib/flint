/*
    Copyright (C) 2021 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_DEFAULT_POLY_FACTOR_H
#define FQ_DEFAULT_POLY_FACTOR_H

#ifdef FQ_DEFAULT_POLY_FACTOR_INLINES_C
#define FQ_DEFAULT_POLY_FACTOR_INLINE
#else
#define FQ_DEFAULT_POLY_FACTOR_INLINE static inline
#endif

#include "thread_pool.h"
#include "nmod_poly_factor.h"
#include "fmpz_mod_poly_factor.h"
#include "fq_poly_factor.h"
#include "fq_nmod_poly_factor.h"
#include "fq_zech_poly_factor.h"
#include "fq_default_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef union fq_default_poly_factor_struct
{
    fq_poly_factor_t fq;
    fq_nmod_poly_factor_t fq_nmod;
    fq_zech_poly_factor_t fq_zech;
    nmod_poly_factor_t nmod;
    fmpz_mod_poly_factor_t fmpz_mod;
} fq_default_poly_factor_struct;

typedef fq_default_poly_factor_struct fq_default_poly_factor_t[1];


FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_init(fq_default_poly_factor_t fac,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_init(fac->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_init(fac->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_init(fac->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        fmpz_mod_poly_factor_init(fac->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_init(fac->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_clear(fq_default_poly_factor_t fac,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_clear(fac->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_clear(fac->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_clear(fac->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_clear(fac->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_clear(fac->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_realloc(fq_default_poly_factor_t fac,
                                       slong alloc, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_realloc(fac->fq_zech, alloc, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_realloc(fac->fq_nmod, alloc, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_realloc(fac->nmod, alloc);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_realloc(fac->fmpz_mod, alloc, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_realloc(fac->fq, alloc, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_fit_length(fq_default_poly_factor_t fac,
                                         slong len, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_fit_length(fac->fq_zech, len, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_fit_length(fac->fq_nmod, len, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_fit_length(fac->nmod, len);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_fit_length(fac->fmpz_mod, len, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_fit_length(fac->fq, len, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
slong fq_default_poly_factor_length(fq_default_poly_factor_t fac,
                                                   const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fac->fq_zech->num;
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fac->fq_nmod->num;
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return fac->nmod->num;
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        return fac->fmpz_mod->num;
    }
    else
    {
        return fac->fq->num;
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
slong fq_default_poly_factor_exp(fq_default_poly_factor_t fac, slong i,
                                                   const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fac->fq_zech->exp[i];
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fac->fq_nmod->exp[i];
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return fac->nmod->exp[i];
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        return fac->fmpz_mod->exp[i];
    }
    else
    {
        return fac->fq->exp[i];
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_set(fq_default_poly_factor_t res,
                const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_set(res->fq_zech, fac->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_set(res->fq_nmod, fac->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_set(res->nmod, fac->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_set(res->fmpz_mod, fac->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_set(res->fq, fac->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_insert(fq_default_poly_factor_t fac,
           const fq_default_poly_t poly, slong exp, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_insert(fac->fq_zech,
                                         poly->fq_zech, exp, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_insert(fac->fq_nmod,
                                         poly->fq_nmod, exp, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_insert(fac->nmod, poly->nmod, exp);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_insert(fac->fmpz_mod,
                                   poly->fmpz_mod, exp, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_insert(fac->fq, poly->fq, exp, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_get_poly(fq_default_poly_t poly,
       const fq_default_poly_factor_t fac, slong i, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_get_poly(poly->fq_zech,
                                            fac->fq_zech, i, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_get_poly(poly->fq_nmod,
                                            fac->fq_nmod, i, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_get_poly(poly->nmod, fac->nmod, i);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_get_poly(poly->fmpz_mod,
                                      fac->fmpz_mod, i, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_get_poly(poly->fq, fac->fq, i, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_print(const fq_default_poly_factor_t fac,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_print(fac->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_print(fac->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_print(fac->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_print(fac->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_print(fac->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_print_pretty(const fq_default_poly_factor_t fac,
                                  const char * var, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_print_pretty(fac->fq_zech, var, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_print_pretty(fac->fq_nmod, var, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_print_pretty(fac->nmod, var);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_print_pretty(fac->fmpz_mod, var,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_print_pretty(fac->fq, var, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_concat(fq_default_poly_factor_t res,
                const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_concat(res->fq_zech, fac->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_concat(res->fq_nmod, fac->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_concat(res->nmod, fac->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_concat(res->fmpz_mod, fac->fmpz_mod,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_concat(res->fq, fac->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_pow(fq_default_poly_factor_t fac,
                                         slong exp, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_pow(fac->fq_zech, exp, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_pow(fac->fq_nmod, exp, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_pow(fac->nmod, exp);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_pow(fac->fmpz_mod, exp, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_pow(fac->fq, exp, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
int fq_default_poly_is_squarefree(const fq_default_poly_t f,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_is_squarefree(f->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_is_squarefree(f->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return nmod_poly_is_squarefree(f->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_is_squarefree(f->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        return fq_poly_is_squarefree(f->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_squarefree(fq_default_poly_factor_t res,
                         const fq_default_poly_t f, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_squarefree(res->fq_zech,
		                                 f->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_squarefree(res->fq_nmod,
		                                 f->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_squarefree(res->nmod, f->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_squarefree(res->fmpz_mod, f->fmpz_mod,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_squarefree(res->fq, f->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
int fq_default_poly_is_irreducible(const fq_default_poly_t f,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_is_irreducible(f->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_is_irreducible(f->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return nmod_poly_is_irreducible(f->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_is_irreducible(f->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        return fq_poly_is_irreducible(f->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_distinct_deg(fq_default_poly_factor_t res,
                        const fq_default_poly_t poly, slong * const * degs,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_distinct_deg(res->fq_zech,
                                        poly->fq_zech, degs, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_distinct_deg(res->fq_nmod,
                                        poly->fq_nmod, degs, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_distinct_deg(res->nmod, poly->nmod, degs);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_distinct_deg(res->fmpz_mod, poly->fmpz_mod, degs,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_distinct_deg(res->fq, poly->fq, degs, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_equal_deg(fq_default_poly_factor_t factors,
              const fq_default_poly_t pol, slong d, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_equal_deg(factors->fq_zech,
                                            pol->fq_zech, d, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_equal_deg(factors->fq_nmod,
                                            pol->fq_nmod, d, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_factor_equal_deg(factors->nmod, pol->nmod, d);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor_equal_deg(factors->fmpz_mod, pol->fmpz_mod, d,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_factor_equal_deg(factors->fq, pol->fq, d, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor(fq_default_poly_factor_t result,
            fq_default_t leading_coeff, const fq_default_poly_t input,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor(result->fq_zech,
                     leading_coeff->fq_zech, input->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor(result->fq_nmod,
                     leading_coeff->fq_nmod, input->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        leading_coeff->nmod = nmod_poly_factor(result->nmod, input->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_factor(result->fmpz_mod, input->fmpz_mod,
                                                        ctx->ctx.fmpz_mod.mod);
        if (input->fmpz_mod->length < 1)
            fmpz_zero(leading_coeff->fmpz_mod);
        else
            fmpz_set(leading_coeff->fmpz_mod, input->fmpz_mod->coeffs +
                                                  input->fmpz_mod->length - 1);
    }
    else
    {
        fq_poly_factor(result->fq, leading_coeff->fq, input->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_factor_split_single(fq_default_poly_t linfactor,
                     const fq_default_poly_t input, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_factor_split_single(linfactor->fq_zech,
                                             input->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_factor_split_single(linfactor->fq_nmod,
                                             input->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        flint_throw(FLINT_ERROR, "operation not implemented");
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        flint_throw(FLINT_ERROR, "operation not implemented");
    }
    else
    {
        fq_poly_factor_split_single(linfactor->fq, input->fq, ctx->ctx.fq);
    }
}

FQ_DEFAULT_POLY_FACTOR_INLINE
void fq_default_poly_roots(fq_default_poly_factor_t r,
                    const fq_default_poly_t f, int with_multiplicity,
                                                    const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_roots(r->fq_zech,
                              f->fq_zech, with_multiplicity, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_roots(r->fq_nmod,
                              f->fq_nmod, with_multiplicity, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_roots(r->nmod, f->nmod, with_multiplicity);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_roots(r->fmpz_mod, f->fmpz_mod, with_multiplicity,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_poly_roots(r->fq, f->fq, with_multiplicity, ctx->ctx.fq);
    }
}

#ifdef __cplusplus
}
#endif

#endif
