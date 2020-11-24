/*
    Copyright (C) 2019-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_H
#define FMPZ_MOD_MPOLY_H

#ifdef FMPZ_MOD_MPOLY_INLINES_C
#define FMPZ_MOD_MPOLY_INLINE FLINT_DLL
#else
#define FMPZ_MOD_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_mod.h"
#include "fmpz_mpoly.h"
#include "mpoly.h"
#include "n_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    fmpz_mod_ctx_t ffinfo;
    mpoly_ctx_t minfo;
} fmpz_mod_mpoly_ctx_struct;

typedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1];


/*
    fmpz_mod_mpolyn_t
    sparse multivariates with fmpz_mod_poly_t coefficients
        with LEX ordering
*/
typedef struct
{
   fmpz_mod_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fmpz_mod_mpolyn_struct;

typedef fmpz_mod_mpolyn_struct fmpz_mod_mpolyn_t[1];

typedef struct
{
    fmpz_mod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;
} fmpz_mod_mpolyun_struct;

typedef fmpz_mod_mpolyun_struct fmpz_mod_mpolyun_t[1];


/* Context object ************************************************************/

FLINT_DLL void fmpz_mod_mpoly_ctx_init(fmpz_mod_mpoly_ctx_t ctx, 
                      slong nvars, const ordering_t ord, const fmpz_t modulus);

FLINT_DLL void fmpz_mod_mpoly_ctx_clear(fmpz_mod_mpoly_ctx_t ctx);


/* mpolyn and mpolyun ********************************************************/

FLINT_DLL void fmpz_mod_mpolyn_init(fmpz_mod_mpolyn_t A, flint_bitcnt_t bits,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE void fmpz_mod_mpolyn_swap(fmpz_mod_mpolyn_t A,
                           fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx)
{
   fmpz_mod_mpolyn_struct t = *A;
   *A = *B;
   *B = t;    
}

FMPZ_MOD_MPOLY_INLINE void fmpz_mod_mpolyn_zero(fmpz_mod_mpolyn_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    A->length = 0;
}

FMPZ_MOD_MPOLY_INLINE int fmpz_mod_mpolyn_is_zero(fmpz_mod_mpolyn_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length == 0;
}

FMPZ_MOD_MPOLY_INLINE fmpz_mod_poly_struct * fmpz_mod_mpolyn_leadcoeff_poly(
                     const fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}


FLINT_DLL void fmpz_mod_mpolyn_fit_length(fmpz_mod_mpolyn_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_fit_length(fmpz_mod_mpolyun_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL fmpz * fmpz_mod_mpolyn_leadcoeff_last_ref(const fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL fmpz_mod_poly_struct * fmpz_mod_mpolyun_leadcoeff_ref(
                         fmpz_mod_mpolyun_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_content_poly(fmpz_mod_poly_t a,
                    const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_content_last(fmpz_mod_poly_t a,
                   const fmpz_mod_mpolyun_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_divexact_poly(fmpz_mod_mpolyn_t A,
                      const fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_divexact_last(fmpz_mod_mpolyun_t A,
                      const fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL ulong nmod_mpolyn_bidegree(const nmod_mpolyn_t A);

FLINT_DLL ulong fmpz_mod_mpolyn_bidegree(const fmpz_mod_mpolyn_t A);

FLINT_DLL slong fmpz_mod_mpolyn_lastdeg(const fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mod_mpolyun_lastdeg(const fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_init(fmpz_mod_mpolyun_t A, flint_bitcnt_t bits,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_clear(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_clear(fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL fmpz * fmpz_mod_mpolyun_leadcoeff_last_ref(const fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_one(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_one(fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_mul_poly(fmpz_mod_mpolyn_t A,
                            fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_mul_last(fmpz_mod_mpolyun_t A, fmpz_mod_poly_t b,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(fmpz_mod_mpolyn_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_scalar_mul_fmpz_mod(fmpz_mod_mpolyun_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_equal(const fmpz_mod_mpolyn_t A,
                   const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyun_equal(const fmpz_mod_mpolyun_t A,
                   const fmpz_mod_mpolyun_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_gcd_brown_bivar(
     fmpz_mod_mpolyn_t G, fmpz_mod_mpolyn_t Abar, fmpz_mod_mpolyn_t Bbar,
     fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_print_pretty(const fmpz_mod_mpolyn_t poly,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_print_pretty(const fmpz_mod_mpolyun_t poly,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

/* interp ********************************************************************/

FLINT_DLL void fmpz_mod_mpolyn_intp_reduce_sm_poly(fmpz_mod_poly_t E,
                            const fmpz_mod_mpolyn_t A, const fmpz_t alpha,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_intp_lift_sm_poly(fmpz_mod_mpolyn_t A,
                      const fmpz_mod_poly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_intp_crt_sm_poly(slong * lastdeg_,
             fmpz_mod_mpolyn_t F, fmpz_mod_mpolyn_t T, const fmpz_mod_poly_t A,
                            const fmpz_mod_poly_t modulus, const fmpz_t alpha,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* helpers *******************************************************************/

typedef struct {
    nmod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} nmod_bma_mpoly_struct;

typedef nmod_bma_mpoly_struct nmod_bma_mpoly_t[1];

typedef struct {
    fmpz_mod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} fmpz_mod_bma_mpoly_struct;

typedef fmpz_mod_bma_mpoly_struct fmpz_mod_bma_mpoly_t[1];

typedef struct
{
    slong * degbounds;
    ulong * subdegs;
    fmpz_mod_discrete_log_pohlig_hellman_t dlogenv;
    nmod_discrete_log_pohlig_hellman_t dlogenv_sp;
} mpoly_bma_interpolate_ctx_struct;
typedef mpoly_bma_interpolate_ctx_struct mpoly_bma_interpolate_ctx_t[1];


FLINT_DLL void nmod_bma_mpoly_init(nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_reset_prime(nmod_bma_mpoly_t A, nmod_t fpctx);

FLINT_DLL void nmod_bma_mpoly_clear(nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_print(const nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_fit_length(nmod_bma_mpoly_t A, slong length,
                                                                 nmod_t fpctx);

FLINT_DLL void nmod_bma_mpoly_zero(nmod_bma_mpoly_t L);

FLINT_DLL int nmod_bma_mpoly_reduce(nmod_bma_mpoly_t L);

FLINT_DLL void nmod_bma_mpoly_add_point(
    nmod_bma_mpoly_t L,
    const n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx_sp);

FLINT_DLL int nmod_bma_mpoly_get_fmpz_mpolyu(fmpz_mpolyu_t A,
      const fmpz_mpoly_ctx_t ctx, ulong alphashift, const nmod_bma_mpoly_t L,
                         const mpoly_bma_interpolate_ctx_t Ictx, nmod_t fpctx);



FLINT_DLL void fmpz_mod_bma_mpoly_init(fmpz_mod_bma_mpoly_t A);

FLINT_DLL void fmpz_mod_bma_mpoly_clear(fmpz_mod_bma_mpoly_t A,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_bma_mpoly_fit_length(fmpz_mod_bma_mpoly_t A,
                                     slong length, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_bma_mpoly_zero(fmpz_mod_bma_mpoly_t L);

FLINT_DLL void fmpz_mod_bma_mpoly_add_point(fmpz_mod_bma_mpoly_t L,
                  const fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx_mp);

FLINT_DLL int fmpz_mod_bma_mpoly_get_fmpz_mpolyu(fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx, const fmpz_t alphashift, const fmpz_mod_bma_mpoly_t L,
           const mpoly_bma_interpolate_ctx_t Ictx, const fmpz_mod_ctx_t fpctx);

FMPZ_MOD_MPOLY_INLINE
void mpoly_bma_interpolate_ctx_init(mpoly_bma_interpolate_ctx_t I, slong nvars)
{
    I->degbounds = (slong *) flint_malloc(nvars*sizeof(slong));
    I->subdegs   = (ulong *) flint_malloc(nvars*sizeof(ulong));
    fmpz_mod_discrete_log_pohlig_hellman_init(I->dlogenv);
    nmod_discrete_log_pohlig_hellman_init(I->dlogenv_sp);
}

FMPZ_MOD_MPOLY_INLINE
void mpoly_bma_interpolate_ctx_clear(mpoly_bma_interpolate_ctx_t I)
{
    flint_free(I->degbounds);
    flint_free(I->subdegs);
    fmpz_mod_discrete_log_pohlig_hellman_clear(I->dlogenv);
    nmod_discrete_log_pohlig_hellman_clear(I->dlogenv_sp);
}

FLINT_DLL int nmod_mpoly_bma_get_fmpz_mpoly(fmpz_mpoly_t A,
     const fmpz_mpoly_ctx_t ctx, ulong alphashift, nmod_berlekamp_massey_t I,
                         const mpoly_bma_interpolate_ctx_t Ictx, nmod_t fpctx);

FLINT_DLL int fmpz_mod_bma_get_fmpz_mpoly(fmpz_mpoly_t A,
     const fmpz_mpoly_ctx_t ctx, const fmpz_t alphashift, fmpz_mod_berlekamp_massey_t I,
           const mpoly_bma_interpolate_ctx_t Ictx, const fmpz_mod_ctx_t fpctx);


FLINT_DLL void nmod_mpoly_bma_interpolate_alpha_powers(mp_limb_t * out,
                            ulong w, const mpoly_bma_interpolate_ctx_t Ictx,
                                     const fmpz_mpoly_ctx_t ctx, nmod_t fpctx);


FLINT_DLL void fmpz_mod_mpoly_bma_interpolate_alpha_powers(fmpz * out,
                     const fmpz_t w, const mpoly_bma_interpolate_ctx_t Ictx,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

/* skel */

FLINT_DLL void fmpz_mod_mpoly_red_skel(fmpz_mpolyc_t Ared, const fmpz_mpoly_t A,
                                                      const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyu_red_skel(fmpz_mpolycu_t Ared, const fmpz_mpolyu_t A,
                                                      const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpoly_copy_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S);

FLINT_DLL void fmpz_mod_mpolyu_copy_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S);

FLINT_DLL void fmpz_mod_mpoly_pow_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S,
                                          ulong k, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyu_pow_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S,
                                          ulong k, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_skel(fmpz_mpolyc_t S,
                      const fmpz_mod_mpoly_ctx_t ctx_mp, const fmpz_mpoly_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyu_set_skel(fmpz_mpolycu_t S,
                     const fmpz_mod_mpoly_ctx_t ctx_mp, const fmpz_mpolyu_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_use_skel_mul(fmpz_t eval, fmpz_mpolyc_t Ared,
                               fmpz_mpolyc_t Avar, const fmpz_mpolyc_t Ainc,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyuu_use_skel_mul(fmpz_mod_mpolyn_t E,
        const fmpz_mpolyu_t A, const fmpz_mpolycu_t Ared, fmpz_mpolycu_t Acur,
                 const fmpz_mpolycu_t Ainc, const fmpz_mod_mpoly_ctx_t ctx_mp);


/* eval */

FLINT_DLL mp_limb_t fmpz_mpoly_eval_nmod(nmod_t fpctx,
    const fmpz_mpoly_t A, const mp_limb_t * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_eval_fmpz_mod(fmpz_t eval,
                        const fmpz_mod_ctx_t fpctx, const fmpz_mpoly_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyuu_eval_nmod(
    n_bpoly_t E,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyuu_eval_fmpz_mod(fmpz_mod_mpolyn_t E,
                 const fmpz_mod_mpoly_ctx_t ctx_mp, const fmpz_mpolyu_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _mpoly_monomial_evals_nmod(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    const mpoly_ctx_t mctx,
    nmod_t fctx);

FLINT_DLL void nmod_mpolyuu_eval_step2(
    n_bpoly_t E,
    n_bpoly_t Acur,
    const n_polyun_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp);


#ifdef __cplusplus
}
#endif

#endif

