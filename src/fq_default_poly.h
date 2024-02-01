/*
    Copyright (C) 2021 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_DEFAULT_POLY_H
#define FQ_DEFAULT_POLY_H

#ifdef FQ_DEFAULT_POLY_INLINES_C
#define FQ_DEFAULT_POLY_INLINE
#else
#define FQ_DEFAULT_POLY_INLINE static inline
#endif

#include "fq_default.h"
#include "fq_vec.h"
#include "fq_poly.h"
#include "fq_nmod_poly.h"
#include "fq_zech_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef union fq_default_poly_struct
{
    fq_poly_t fq;
    fq_nmod_poly_t fq_nmod;
    fq_zech_poly_t fq_zech;
    nmod_poly_t nmod;
    fmpz_mod_poly_t fmpz_mod;
}
fq_default_poly_struct;

typedef fq_default_poly_struct fq_default_poly_t[1];

/*  Memory management ********************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_init(fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_init(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_init(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_init(poly->nmod, FQ_DEFAULT_CTX_NMOD(ctx).n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_init(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_init(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_init2(fq_default_poly_t poly,
                                       slong alloc, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_init2(poly->fq_zech, alloc, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_init2(poly->fq_nmod, alloc, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_init2(poly->nmod, FQ_DEFAULT_CTX_NMOD(ctx).n, alloc);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_init2(poly->fmpz_mod, alloc, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_init2(poly->fq, alloc, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_realloc(fq_default_poly_t poly,
                                       slong alloc, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_realloc(poly->fq_zech, alloc, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_realloc(poly->fq_nmod, alloc, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_realloc(poly->nmod, alloc);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_realloc(poly->fmpz_mod, alloc, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_realloc(poly->fq, alloc, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_truncate(fq_default_poly_t poly,
                                         slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_truncate(poly->fq_zech, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_truncate(poly->fq_nmod, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_truncate(poly->nmod, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_truncate(poly->fmpz_mod, len, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_truncate(poly->fq, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_set_trunc(fq_default_poly_t poly1,
                  fq_default_poly_t poly2, slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_set_trunc(poly1->fq_zech,
		                        poly2->fq_zech, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_set_trunc(poly1->fq_nmod,
		                        poly2->fq_nmod, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_set_trunc(poly1->nmod, poly2->nmod, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_set_trunc(poly1->fmpz_mod,
                                  poly2->fmpz_mod, len, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_set_trunc(poly1->fq, poly2->fq, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_fit_length(fq_default_poly_t poly,
                                         slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_fit_length(poly->fq_zech, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_fit_length(poly->fq_nmod, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_fit_length(poly->nmod, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_fit_length(poly->fmpz_mod, len, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_fit_length(poly->fq, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void _fq_default_poly_set_length(fq_default_poly_t poly,
                                         slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        _fq_zech_poly_set_length(poly->fq_zech, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        _fq_nmod_poly_set_length(poly->fq_nmod, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        _nmod_poly_set_length(poly->nmod, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        _fmpz_mod_poly_set_length(poly->fmpz_mod, len);
    }
    else
    {
        _fq_poly_set_length(poly->fq, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_clear(fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_clear(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_clear(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_clear(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_clear(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_clear(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Polynomial parameters  ***************************************************/

FQ_DEFAULT_POLY_INLINE slong
fq_default_poly_length(const fq_default_poly_t poly,
		                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_length(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_length(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_length(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_length(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_length(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE slong
fq_default_poly_degree(const fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_degree(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_degree(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_degree(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_degree(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_degree(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Randomisation  ***********************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_randtest(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_randtest(f->fq_zech, state, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_randtest(f->fq_nmod, state, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_randtest(f->nmod, state, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_randtest(f->fmpz_mod, state, len, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_randtest(f->fq, state, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_randtest_not_zero(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_randtest_not_zero(f->fq_zech, state, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_randtest_not_zero(f->fq_nmod, state, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_randtest_not_zero(f->nmod, state, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_randtest_not_zero(f->fmpz_mod, state, len,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_randtest_not_zero(f->fq, state, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_randtest_monic(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_randtest_monic(f->fq_zech, state, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_randtest_monic(f->fq_nmod, state, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_randtest_monic(f->nmod, state, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_randtest_monic(f->fmpz_mod, state, len,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_randtest_monic(f->fq, state, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_randtest_irreducible(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_randtest_irreducible(f->fq_zech,
                                                 state, len, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_randtest_irreducible(f->fq_nmod,
                                                 state, len, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_randtest_irreducible(f->nmod, state, len);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_randtest_irreducible(f->fmpz_mod, state, len,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_randtest_irreducible(f->fq, state, len, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Assignment and basic manipulation  ***************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_set(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_set(rop->fq_zech, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_set(rop->fq_nmod, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_set(rop->nmod, op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_set(rop->fmpz_mod, op->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_set(rop->fq, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_set_fq_default(fq_default_poly_t poly,
                              const fq_default_t c, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_set_fq_zech(poly->fq_zech, c->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_set_fq_nmod(poly->fq_nmod, c->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_zero(poly->nmod);
        nmod_poly_set_coeff_ui(poly->nmod, 0, c->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_set_fmpz(poly->fmpz_mod, c->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_set_fq(poly->fq, c->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_swap(fq_default_poly_t op1,
                               fq_default_poly_t op2, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_swap(op1->fq_zech, op2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_swap(op1->fq_nmod, op2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_swap(op1->nmod, op2->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_swap(op1->fmpz_mod, op2->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_swap(op1->fq, op2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_zero(fq_default_poly_t poly, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_zero(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_zero(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_zero(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_zero(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_zero(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_one(fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_one(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_one(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_one(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_one(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_one(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_gen(fq_default_poly_t f,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_gen(f->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_gen(f->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_one(f->nmod);
        nmod_poly_shift_left(f->nmod, f->nmod, 1);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_gen(f->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_gen(f->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_make_monic(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_make_monic(rop->fq_zech, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_make_monic(rop->fq_nmod, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_make_monic(rop->nmod, op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_make_monic(rop->fmpz_mod, op->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_make_monic(rop->fq, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_reverse(fq_default_poly_t res,
             const fq_default_poly_t poly, slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_reverse(res->fq_zech, poly->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_reverse(res->fq_nmod, poly->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_reverse(res->nmod, poly->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_reverse(res->fmpz_mod, poly->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_reverse(res->fq, poly->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
ulong fq_default_poly_deflation(const fq_default_poly_t input,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_deflation(input->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_deflation(input->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_deflation(input->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_deflation(input->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_deflation(input->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_deflate(fq_default_poly_t result,
    const fq_default_poly_t input, ulong deflation, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_deflate(result->fq_zech,
                                  input->fq_zech, deflation, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_deflate(result->fq_nmod,
                                  input->fq_nmod, deflation, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_deflate(result->nmod, input->nmod, deflation);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_deflate(result->fmpz_mod, input->fmpz_mod, deflation,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_deflate(result->fq, input->fq, deflation, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_inflate(fq_default_poly_t result,
    const fq_default_poly_t input, ulong inflation, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_inflate(result->fq_zech,
                                  input->fq_zech, inflation, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_inflate(result->fq_nmod,
                                  input->fq_nmod, inflation, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_inflate(result->nmod, input->nmod, inflation);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_inflate(result->fmpz_mod, input->fmpz_mod, inflation,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_inflate(result->fq, input->fq, inflation, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Getting and setting coefficients  ****************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_get_coeff(fq_default_t x,
             const fq_default_poly_t poly, slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_get_coeff(x->fq_zech, poly->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_get_coeff(x->fq_nmod, poly->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        x->nmod = nmod_poly_get_coeff_ui(poly->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_get_coeff_fmpz(x->fmpz_mod, poly->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_get_coeff(x->fq, poly->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_set_coeff(fq_default_poly_t poly,
                     slong n, const fq_default_t x, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_set_coeff(poly->fq_zech, n, x->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_set_coeff(poly->fq_nmod, n, x->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_set_coeff_ui(poly->nmod, n, x->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_set_coeff_fmpz(poly->fmpz_mod, n, x->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_set_coeff(poly->fq, n, x->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_set_coeff_fmpz(fq_default_poly_t poly,
                           slong n, const fmpz_t x, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_set_coeff_fmpz(poly->fq_zech, n, x, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_set_coeff_fmpz(poly->fq_nmod, n, x, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_set_coeff_ui(poly->nmod, n, fmpz_get_nmod(x, FQ_DEFAULT_CTX_NMOD(ctx)));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_set_coeff_fmpz(poly->fmpz_mod, n, x,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_set_coeff_fmpz(poly->fq, n, x, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_set_nmod_poly(fq_default_poly_t rop,
                          const     nmod_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_set_nmod_poly(rop->fq_zech, op, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_set_nmod_poly(rop->fq_nmod, op, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_set(rop->nmod, op);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_set_nmod_poly(rop->fmpz_mod, op);
    }
    else
    {
        fq_poly_set_nmod_poly(rop->fq, op, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_set_fmpz_mod_poly(fq_default_poly_t rop,
                          const fmpz_mod_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_set_fmpz_mod_poly(rop->fq_zech, op, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_set_fmpz_mod_poly(rop->fq_nmod, op, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        fmpz_mod_poly_get_nmod_poly(rop->nmod, op);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_set(rop->fmpz_mod, op, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_set_fmpz_mod_poly(rop->fq, op, FQ_DEFAULT_CTX_FQ(ctx));
    }
}


void fq_default_poly_set_fmpz_poly(fq_default_poly_t rop,
                             const fmpz_poly_t op, const fq_default_ctx_t ctx);

/*  Comparison  **************************************************************/

FQ_DEFAULT_POLY_INLINE
int fq_default_poly_equal(const fq_default_poly_t poly1,
                     const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_equal(poly1->fq_zech, poly2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_equal(poly1->fq_nmod, poly2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_equal(poly1->nmod, poly2->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_equal(poly1->fmpz_mod, poly2->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_equal(poly1->fq, poly2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
int fq_default_poly_equal_trunc(const fq_default_poly_t poly1,
            const fq_default_poly_t poly2, slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_equal_trunc(poly1->fq_zech,
                                          poly2->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_equal_trunc(poly1->fq_nmod,
                                          poly2->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_equal_trunc(poly1->nmod, poly2->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_equal_trunc(poly1->fmpz_mod, poly2->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_equal_trunc(poly1->fq, poly2->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_zero(const fq_default_poly_t poly,
		                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_is_zero(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_is_zero(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_is_zero(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_is_zero(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_is_zero(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_one(const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_is_one(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_is_one(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_is_one(op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_is_one(op->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
       return fq_poly_is_one(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_unit(const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_is_unit(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_is_unit(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_is_unit(op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_is_unit(op->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_is_unit(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_gen(const fq_default_poly_t poly,
		                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_is_gen(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_is_gen(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_is_gen(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_is_gen(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_is_gen(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_equal_fq_default(const fq_default_poly_t poly,
                              const fq_default_t c, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_equal_fq_zech(poly->fq_zech,
		                                 c->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_equal_fq_nmod(poly->fq_nmod,
		                                 c->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        if (c->nmod == 0)
            return poly->nmod->length == 0;
        else
            return poly->nmod->length == 1 && poly->nmod->coeffs[0] == c->nmod;
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        if (fmpz_is_zero(c->fmpz_mod))
            return poly->fmpz_mod->length == 0;
        else
            return poly->fmpz_mod->length == 1 &&
                   fmpz_equal(poly->fmpz_mod->coeffs + 0, c->fmpz_mod);
    }
    else
    {
        return fq_poly_equal_fq(poly->fq, c->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Addition and subtraction  ************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_add(fq_default_poly_t rop, const fq_default_poly_t op1,
                       const fq_default_poly_t op2, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_add(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_add(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_add(rop->nmod, op1->nmod, op2->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_add(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_add(rop->fq, op1->fq, op2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_add_si(fq_default_poly_t rop,
              const fq_default_poly_t op1, slong c, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_add_si(rop->fq_zech, op1->fq_zech, c, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_add_si(rop->fq_nmod, op1->fq_nmod, c, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_add_ui(rop->nmod, op1->nmod, nmod_set_si(c, FQ_DEFAULT_CTX_NMOD(ctx)));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_add_si(rop->fmpz_mod, op1->fmpz_mod, c,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_add_si(rop->fq, op1->fq, c, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_add_series(fq_default_poly_t rop,
		const fq_default_poly_t op1, const fq_default_poly_t op2,
                                           slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_add_series(rop->fq_zech,
                              op1->fq_zech, op2->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_add_series(rop->fq_nmod,
                              op1->fq_nmod, op2->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_add_series(rop->nmod, op1->nmod, op2->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_add_series(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_add_series(rop->fq, op1->fq, op2->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_sub(fq_default_poly_t rop,
                const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_sub(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_sub(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_sub(rop->nmod, op1->nmod, op2->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_sub(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_sub(rop->fq, op1->fq, op2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_sub_series(fq_default_poly_t rop,
		const fq_default_poly_t op1, const fq_default_poly_t op2,
                                           slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_sub_series(rop->fq_zech,
                              op1->fq_zech, op2->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_sub_series(rop->fq_nmod,
                              op1->fq_nmod, op2->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_sub_series(rop->nmod, op1->nmod, op2->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_sub_series(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_sub_series(rop->fq, op1->fq, op2->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_neg(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_neg(rop->fq_zech, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_neg(rop->fq_nmod, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_neg(rop->nmod, op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_neg(rop->fmpz_mod, op->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_neg(rop->fq, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Scalar multiplication and division  **************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_mul_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_scalar_mul_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_scalar_mul_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_scalar_mul_nmod(rop->nmod, op->nmod, x->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_scalar_mul_fmpz(rop->fmpz_mod, op->fmpz_mod, x->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_scalar_mul_fq(rop->fq, op->fq, x->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_div_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_scalar_div_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_scalar_div_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_scalar_mul_nmod(rop->nmod, op->nmod,
                                         nmod_inv(x->nmod, FQ_DEFAULT_CTX_NMOD(ctx)));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mod_inv(t, x->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
        fmpz_mod_poly_scalar_mul_fmpz(rop->fmpz_mod, op->fmpz_mod, t,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
        fmpz_clear(t);
    }
    else
    {
        fq_poly_scalar_div_fq(rop->fq, op->fq, x->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_addmul_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_scalar_addmul_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_scalar_addmul_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_scalar_addmul_nmod(rop->nmod, op->nmod, x->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_scalar_addmul_fmpz(rop->fmpz_mod, op->fmpz_mod,
                                           x->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_scalar_addmul_fq(rop->fq, op->fq, x->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_submul_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_scalar_submul_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_scalar_submul_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_scalar_addmul_nmod(rop->nmod, op->nmod,
                                         nmod_neg(x->nmod, FQ_DEFAULT_CTX_NMOD(ctx)));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mod_neg(t, x->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
        fmpz_mod_poly_scalar_addmul_fmpz(rop->fmpz_mod, op->fmpz_mod, t,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
        fmpz_clear(t);
    }
    else
    {
        fq_poly_scalar_submul_fq(rop->fq, op->fq, x->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Multiplication  **********************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mul(fq_default_poly_t rop,
                 const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_mul(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_mul(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_mul(rop->nmod, op1->nmod, op2->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_mul(rop->fmpz_mod, op1->fmpz_mod,
                                         op2->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_mul(rop->fq, op1->fq, op2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mullow(fq_default_poly_t rop,
                 const fq_default_poly_t op1, const fq_default_poly_t op2,
                                           slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_mullow(rop->fq_zech,
                              op1->fq_zech, op2->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_mullow(rop->fq_nmod,
                              op1->fq_nmod, op2->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_mullow(rop->nmod, op1->nmod, op2->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_mullow(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_mullow(rop->fq, op1->fq, op2->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mulhigh(fq_default_poly_t rop,
                 const fq_default_poly_t op1, const fq_default_poly_t op2,
                                       slong start, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_mulhigh(rop->fq_zech,
                          op1->fq_zech, op2->fq_zech, start, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_mulhigh(rop->fq_nmod,
                          op1->fq_nmod, op2->fq_nmod, start, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_mulhigh(rop->nmod, op1->nmod, op2->nmod, start);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_mulhigh(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod,
                                                 start, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_mulhigh(rop->fq, op1->fq, op2->fq, start, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mulmod(fq_default_poly_t res,
                 const fq_default_poly_t poly1, const fq_default_poly_t poly2,
                         const fq_default_poly_t f, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_mulmod(res->fq_zech, poly1->fq_zech, poly2->fq_zech,
                                                 f->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_mulmod(res->fq_nmod, poly1->fq_nmod, poly2->fq_nmod,
                                                 f->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_mulmod(res->nmod, poly1->nmod, poly2->nmod, f->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_mulmod(res->fmpz_mod, poly1->fmpz_mod, poly2->fmpz_mod,
                                           f->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_mulmod(res->fq, poly1->fq, poly2->fq, f->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Squaring ******************************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_sqr(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_sqr(rop->fq_zech, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_sqr(rop->fq_nmod, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_mul(rop->nmod, op->nmod, op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_sqr(rop->fmpz_mod, op->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_sqr(rop->fq, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Powering  ****************************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_pow(fq_default_poly_t rop,
                                      const fq_default_poly_t op, ulong e,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_pow(rop->fq_zech, op->fq_zech, e, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_pow(rop->fq_nmod, op->fq_nmod, e, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_pow(rop->nmod, op->nmod, e);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_pow(rop->fmpz_mod, op->fmpz_mod, e, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_pow(rop->fq, op->fq, e, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_pow_trunc(fq_default_poly_t res,
                                    const fq_default_poly_t poly, ulong e,
				       slong trunc, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_pow_trunc(res->fq_zech,
                                    poly->fq_zech, e, trunc, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_pow_trunc(res->fq_nmod,
                                    poly->fq_nmod, e, trunc, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_pow_trunc(res->nmod, poly->nmod, e, trunc);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_pow_trunc(res->fmpz_mod, poly->fmpz_mod, e, trunc,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_pow_trunc(res->fq, poly->fq, e, trunc, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_powmod_fmpz_binexp(fq_default_poly_t res,
                             const fq_default_poly_t poly, const fmpz_t e,
                         const fq_default_poly_t f, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_powmod_fmpz_binexp(res->fq_zech,
                               poly->fq_zech, e, f->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_powmod_fmpz_binexp(res->fq_nmod,
                               poly->fq_nmod, e, f->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_powmod_fmpz_binexp(res->nmod, poly->nmod, (fmpz*) e, f->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_powmod_fmpz_binexp(res->fmpz_mod, poly->fmpz_mod, e,
                                           f->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_powmod_fmpz_binexp(res->fq, poly->fq, e, f->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_powmod_ui_binexp(fq_default_poly_t res,
                         const fq_default_poly_t poly, ulong e,
                         const fq_default_poly_t f, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_powmod_ui_binexp(res->fq_zech,
                               poly->fq_zech, e, f->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_powmod_ui_binexp(res->fq_nmod,
                               poly->fq_nmod, e, f->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_powmod_ui_binexp(res->nmod, poly->nmod, e, f->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_powmod_ui_binexp(res->fmpz_mod, poly->fmpz_mod, e,
                                           f->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_powmod_ui_binexp(res->fq, poly->fq, e, f->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Shifting  ****************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_shift_left(fq_default_poly_t rop,
               const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_shift_left(rop->fq_zech, op->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_shift_left(rop->fq_nmod, op->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_shift_left(rop->nmod, op->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_shift_left(rop->fmpz_mod, op->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_shift_left(rop->fq, op->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_shift_right(fq_default_poly_t rop,
               const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_shift_right(rop->fq_zech, op->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_shift_right(rop->fq_nmod, op->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_shift_right(rop->nmod, op->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_shift_right(rop->fmpz_mod, op->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_shift_right(rop->fq, op->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Norms  *******************************************************************/

FQ_DEFAULT_POLY_INLINE
slong fq_default_poly_hamming_weight(const fq_default_poly_t op,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_hamming_weight(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_hamming_weight(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_hamming_weight(op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_hamming_weight(op->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_hamming_weight(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Greatest common divisor  *************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_gcd(fq_default_poly_t rop,
                  const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_gcd(rop->fq_zech,
		                 op1->fq_zech, op2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_gcd(rop->fq_nmod,
		                 op1->fq_nmod, op2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_gcd(rop->nmod, op1->nmod, op2->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_gcd(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_gcd(rop->fq, op1->fq, op2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_xgcd(fq_default_poly_t G,
                         fq_default_poly_t S, fq_default_poly_t T,
                  const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_xgcd(G->fq_zech,
             S->fq_zech, T->fq_zech, A->fq_zech, B->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_xgcd(G->fq_nmod,
             S->fq_nmod, T->fq_nmod, A->fq_nmod, B->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_xgcd(G->nmod, S->nmod, T->nmod, A->nmod, B->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_xgcd(G->fmpz_mod, S->fmpz_mod, T->fmpz_mod, A->fmpz_mod,
                                           B->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_xgcd(G->fq, S->fq, T->fq, A->fq, B->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Euclidean division  ******************************************************/

FQ_DEFAULT_POLY_INLINE ulong fq_default_poly_remove(fq_default_poly_t f,
                         const fq_default_poly_t g, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_remove(f->fq_zech, g->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_remove(f->fq_nmod, g->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_remove(f->nmod, g->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_remove(f->fmpz_mod, g->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_remove(f->fq, g->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_divrem(fq_default_poly_t Q, fq_default_poly_t R,
                   const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_divrem(Q->fq_zech,
                         R->fq_zech, A->fq_zech, B->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_divrem(Q->fq_nmod,
                         R->fq_nmod, A->fq_nmod, B->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_divrem(Q->nmod, R->nmod, A->nmod, B->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_divrem(Q->fmpz_mod, R->fmpz_mod, A->fmpz_mod,
                                           B->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_divrem(Q->fq, R->fq, A->fq, B->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_rem(fq_default_poly_t R,
             const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_rem(R->fq_zech, A->fq_zech, B->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_rem(R->fq_nmod, A->fq_nmod, B->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_rem(R->nmod, A->nmod, B->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_rem(R->fmpz_mod, A->fmpz_mod, B->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_rem(R->fq, A->fq, B->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_inv_series(fq_default_poly_t Qinv,
                               const fq_default_poly_t Q, slong n,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_inv_series(Qinv->fq_zech, Q->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_inv_series(Qinv->fq_nmod, Q->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_inv_series(Qinv->nmod, Q->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_inv_series(Qinv->fmpz_mod, Q->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_inv_series(Qinv->fq, Q->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_div_series(fq_default_poly_t Q,
               const fq_default_poly_t A, const fq_default_poly_t B,
                                           slong n, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_div_series(Q->fq_zech,
                                  A->fq_zech, B->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_div_series(Q->fq_nmod,
                                  A->fq_nmod, B->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_div_series(Q->nmod, A->nmod, B->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_div_series(Q->fmpz_mod, A->fmpz_mod, B->fmpz_mod, n,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_div_series(Q->fq, A->fq, B->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Divisibility testing  ****************************************************/

FQ_DEFAULT_POLY_INLINE int fq_default_poly_divides(fq_default_poly_t Q,
                      const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_divides(Q->fq_zech,
                                     A->fq_zech, B->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_divides(Q->fq_nmod,
                                     A->fq_nmod, B->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_divides(Q->nmod, A->nmod, B->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_divides(Q->fmpz_mod, A->fmpz_mod, B->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_divides(Q->fq, A->fq, B->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Derivative  **************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_derivative(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_derivative(rop->fq_zech, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_derivative(rop->fq_nmod, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_derivative(rop->nmod, op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_derivative(rop->fmpz_mod, op->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_derivative(rop->fq, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Square root ***************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_invsqrt_series(fq_default_poly_t rop,
                    const fq_default_poly_t op, slong n, fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_invsqrt_series(rop->fq_zech, op->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_invsqrt_series(rop->fq_nmod, op->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_invsqrt_series(rop->nmod, op->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_invsqrt_series(rop->fmpz_mod, op->fmpz_mod, n, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_invsqrt_series(rop->fq, op->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_sqrt_series(fq_default_poly_t rop,
                    const fq_default_poly_t op, slong n, fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_sqrt_series(rop->fq_zech, op->fq_zech, n, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_sqrt_series(rop->fq_nmod, op->fq_nmod, n, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_sqrt_series(rop->nmod, op->nmod, n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_sqrt_series(rop->fmpz_mod, op->fmpz_mod, n, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_sqrt_series(rop->fq, op->fq, n, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
int fq_default_poly_sqrt(fq_default_poly_t rop,
                              const fq_default_poly_t op, fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_sqrt(rop->fq_zech, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_sqrt(rop->fq_nmod, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_sqrt(rop->nmod, op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_sqrt(rop->fmpz_mod, op->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_sqrt(rop->fq, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Evaluation  **************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_evaluate_fq_default(fq_default_t res,
                           const fq_default_poly_t f, const fq_default_t a,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_evaluate_fq_zech(res->fq_zech,
		                     f->fq_zech, a->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_evaluate_fq_nmod(res->fq_nmod,
		                     f->fq_nmod, a->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        res->nmod = nmod_poly_evaluate_nmod(f->nmod, a->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        fmpz_mod_poly_evaluate_fmpz(res->fmpz_mod, f->fmpz_mod, a->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_evaluate_fq(res->fq, f->fq, a->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Composition  *************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_compose(fq_default_poly_t rop,
            const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_compose(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_compose(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_compose(rop->nmod, op1->nmod, op2->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_compose(rop->fmpz_mod, op1->fmpz_mod, op2->fmpz_mod,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_compose(rop->fq, op1->fq, op2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_compose_mod(fq_default_poly_t res,
            const fq_default_poly_t poly1, const fq_default_poly_t poly2,
                     const fq_default_poly_t poly3, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_poly_compose_mod(res->fq_zech,
             poly1->fq_zech, poly2->fq_zech, poly3->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_poly_compose_mod(res->fq_nmod,
             poly1->fq_nmod, poly2->fq_nmod, poly3->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_poly_compose_mod(res->nmod, poly1->nmod, poly2->nmod, poly3->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_compose_mod(res->fmpz_mod, poly1->fmpz_mod,
                      poly2->fmpz_mod, poly3->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_poly_compose_mod(res->fq,
                                 poly1->fq, poly2->fq, poly3->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/*  Input and output  ********************************************************/

#ifdef FLINT_HAVE_FILE
int fq_default_poly_fprint(FILE * file, const fq_default_poly_t poly, const fq_default_ctx_t ctx);
int fq_default_poly_fprint_pretty(FILE * file, const fq_default_poly_t poly, const char *x, const fq_default_ctx_t ctx);
#endif

int fq_default_poly_print(const fq_default_poly_t poly, const fq_default_ctx_t ctx);
int fq_default_poly_print_pretty(const fq_default_poly_t poly, const char *x, const fq_default_ctx_t ctx);

FQ_DEFAULT_POLY_INLINE
char * fq_default_poly_get_str_pretty(const fq_default_poly_t poly,
                                     const char *x, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_get_str_pretty(poly->fq_zech, x, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_get_str_pretty(poly->fq_nmod, x, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_get_str_pretty(poly->nmod, x);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_get_str_pretty(poly->fmpz_mod, x,
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_get_str_pretty(poly->fq, x, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_POLY_INLINE
char * fq_default_poly_get_str(const fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_get_str(poly->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_get_str(poly->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return nmod_poly_get_str(poly->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_get_str(poly->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        return fq_poly_get_str(poly->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

#include "fq_default_mat.h"

/* Characteristic polynomial *************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_mat_charpoly(fq_default_poly_t p,
                          const fq_default_mat_t M, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_mat_charpoly(p->fq_zech, M->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_mat_charpoly(p->fq_nmod, M->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_mat_charpoly(p->nmod, M->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_mat_charpoly(p->fmpz_mod, M->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_mat_charpoly(p->fq, M->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Minimal polynomial ********************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_mat_minpoly(fq_default_poly_t p,
                          const fq_default_mat_t X, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_mat_minpoly(p->fq_zech, X->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_mat_minpoly(p->fq_nmod, X->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        nmod_mat_minpoly(p->nmod, X->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_mat_minpoly(p->fmpz_mod, X->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_mat_minpoly(p->fq, X->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

#ifdef __cplusplus
}
#endif

#endif
