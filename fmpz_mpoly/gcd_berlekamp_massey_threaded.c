/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mod_mpoly.h"
#include "thread_pool.h"
#include "n_poly.h"

#define LOW_HALF_MASK ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2))

typedef struct
{
    n_bpoly_t nGeval_sp, nAbareval_sp, nBbareval_sp;
    ulong GevaldegXY;
    int success;
}
_eval_sp_ret_struct;

typedef struct
{
    fmpz_mod_mpolyn_t nGeval_mp, nAbareval_mp, nBbareval_mp;
    ulong GevaldegXY;
    int success;
}
_eval_mp_ret_struct;

/*
    Most of the objects used by fmpz_mpolyuu_gcd_berlekamp_massey_threaded
    need to be made visible to other threads. They are all passed via a
    _base_t, of which there is only one instance on the stack of 
    fmpz_mpolyuu_gcd_berlekamp_massey_threaded, and this instance is
    succinctly referred to as "w".

    When calculations are done modulo a machine prime, the "_sp" suffix is
    used, otherwise "_mp" is used. The _base_t w has two array pointers
    evals_{sp|mp} of length evals_{sp|mp}_alloc for returning evaluations.
*/
typedef struct
{
    volatile slong index;
    volatile int zip_find_coeffs_no_match, zip_find_coeffs_non_invertible;
    volatile int Gchanged, Abarchanged, Bbarchanged, Hchanged, failed;
    int which_check;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    slong num_threads;
    flint_bitcnt_t bits;
    const fmpz_mpolyu_struct * A, * B;
    fmpz_mpolyu_struct * Abar, * Bbar;
    const fmpz_mpoly_struct * Gamma;
    slong * Adegs, * Bdegs;
    fmpz_mpolyu_t H;
    fmpz_t Hmodulus;

    nmod_zip_mpolyu_t Z;
    slong bma_target_count;
    ulong GdegboundXY;

    flint_rand_t randstate;
    mpoly_bma_interpolate_ctx_t Ictx;
    const fmpz_mpoly_ctx_struct * ctx;
    nmod_mpoly_ctx_t ctx_sp;
    fmpz_mod_mpoly_ctx_t ctx_mp;

    nmod_bma_mpoly_t GLambda_sp, AbarLambda_sp, BbarLambda_sp;
    fmpz_mod_bma_mpoly_t GLambda_mp, AbarLambda_mp, BbarLambda_mp;

    n_bpoly_t Aone_sp, Bone_sp;
    n_poly_t Gammaone_sp;

    n_polyun_t Ainc_sp, Binc_sp;
    n_poly_t Gammainc_sp;

    fmpz_mpolycu_t Aone_mp, Ainc_mp, Ared_mp;
    fmpz_mpolycu_t Bone_mp, Binc_mp, Bred_mp;
    fmpz_mpolyc_t Gammaone_mp, Gammainc_mp, Gammared_mp;

    ulong alphashift;
    fmpz_t alphashift_mp;
    mp_limb_t * alphas_sp;
    fmpz * alphas_mp;

    slong num_images_sp;
    slong evals_sp_alloc;
    _eval_sp_ret_struct * evals_sp;
    n_poly_t coeff_evals_sp;

    slong num_images_mp;
    slong evals_mp_alloc;
    _eval_mp_ret_struct * evals_mp;
    fmpz_mpolyc_t coeff_evals_mp;
}
_base_struct;

typedef _base_struct _base_t[1];

/*
    The working memory for each thread is in a _eval_{sp|mp}_worker_arg_struct.
    The function fmpz_mpolyuu_gcd_berlekamp_massey_threaded will allocate an
    array of each of length w->num_threads.
*/
typedef struct
{
    _base_struct * w;
    n_bpoly_t Aeval_sp, Beval_sp, Geval_sp, Abareval_sp, Bbareval_sp;
    n_bpoly_t Acur_sp, Bcur_sp;
    n_poly_t Gammacur_sp;
    n_poly_bpoly_stack_t Sp_sp;
    slong thread_index;
    int success;
    int cur_is_uninited;
}
_eval_sp_worker_arg_struct;

typedef struct
{
    _base_struct * w;
    fmpz_mod_mpolyn_t Aeval_mp, Beval_mp, Geval_mp, Abareval_mp, Bbareval_mp;
    fmpz_mpolycu_t Acur_mp, Bcur_mp;
    fmpz_mpolyc_t Gammacur_mp;
    slong thread_index;
    int success;
    int cur_is_uninited;
}
_eval_mp_worker_arg_struct;



/* misc helpers **************************************************************/

void fmpz_mod_mpoly_pow_skel(
    fmpz_mpolyc_t M,
    const fmpz_mpolyc_t S,
    ulong k,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
    slong i;
    fmpz_mpolyc_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        fmpz_mod_pow_ui(M->coeffs + i, S->coeffs + i, k, ctx_mp->ffinfo);
    }
}

void fmpz_mod_mpolyu_pow_skel(
    fmpz_mpolycu_t M,
    const fmpz_mpolycu_t S,
    ulong k,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
    slong i;
    fmpz_mpolycu_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        fmpz_mod_mpoly_pow_skel(M->coeffs + i, S->coeffs + i, k, ctx_mp);
    }
}

void fmpz_mpolycu_set_length(fmpz_mpolycu_t A, slong k)
{
    fmpz_mpolycu_fit_length(A, k);
    A->length = k;
}

/*
    Given the evaluations of the coefficients of A, construct its evaluation E.

    A is in ZZ[X,Y][x_0,..., x_(n-1)]
    E is in Fp[X][][Y]
*/
static void _fmpz_mpolyuu_eval_nmod_from_coeffs(
    n_bpoly_t E,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    mp_limb_t * evals)
{
    slong i;
    slong xexp, yexp;

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        mp_limb_t eval = evals[i];
        if (eval == 0)
            continue;

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & LOW_HALF_MASK;

        n_bpoly_set_coeff(E, xexp, yexp, eval);
    }
}

/* coeffs may be clobbered */
static void _fmpz_mpolyuu_eval_fmpz_mod_from_coeffs(
    fmpz_mod_mpolyn_t E,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    fmpz * coeffs)
{
    slong i;
    slong xexp, yexp;

    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_mp->minfo));

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        if (fmpz_is_zero(coeffs + i))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & (-UWORD(1) >> (FLINT_BITS - FLINT_BITS/2));

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length - 1, yexp, coeffs + i, ctx_mp->ffinfo);
        }
        else
        {
            fmpz_mod_mpolyn_fit_length(E, E->length + 1, ctx_mp);
            fmpz_mod_poly_zero(E->coeffs + E->length, ctx_mp->ffinfo);
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length, yexp, coeffs + i, ctx_mp->ffinfo);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }
}

/*
    Set the mod member of all relevent poly's. This function could be a
    nop if the poly's did not store their own ctx.
*/
void _base_args_set_mod_sp(
    _base_struct * w,
    _eval_sp_worker_arg_struct * args)
{
    return;
}

void _base_args_set_mod_mp(
    _base_struct * w,
    _eval_mp_worker_arg_struct * args)
{
    return;
}

/* Set w->num_images and fit the length of w->evals */
void _base_set_num_images_sp(_base_struct * w, slong len)
{
    slong i;

    w->num_images_sp = len;

    if (w->evals_sp_alloc < w->num_images_sp)
    {
        w->evals_sp = flint_realloc(w->evals_sp,
                                 w->num_images_sp*sizeof(_eval_sp_ret_struct));

        for (i = w->evals_sp_alloc; i < w->num_images_sp; i++)
        {
            n_bpoly_init(w->evals_sp[i].nGeval_sp);
            n_bpoly_init(w->evals_sp[i].nAbareval_sp);
            n_bpoly_init(w->evals_sp[i].nBbareval_sp);
        }

        w->evals_sp_alloc = w->num_images_sp;
    }
}

void _base_set_num_images_mp(_base_struct * w, slong len)
{
    slong i;

    w->num_images_mp = len;

    if (w->evals_mp_alloc < w->num_images_mp)
    {
        w->evals_mp = flint_realloc(w->evals_mp,
                                 w->num_images_mp*sizeof(_eval_mp_ret_struct));

        for (i = w->evals_mp_alloc; i < w->num_images_mp; i++)
        {
            fmpz_mod_mpolyn_init(w->evals_mp[i].nGeval_mp, FLINT_BITS/2, w->ctx_mp);
            fmpz_mod_mpolyn_init(w->evals_mp[i].nAbareval_mp, FLINT_BITS/2, w->ctx_mp);
            fmpz_mod_mpolyn_init(w->evals_mp[i].nBbareval_mp, FLINT_BITS/2, w->ctx_mp);
        }

        w->evals_mp_alloc = w->num_images_mp;
    }
}

/* worker functions **********************************************************/

static void nmod_mpoly_get_eval_helper_pow(
    n_poly_t Acur,
    n_poly_t Ainc,
    ulong pow,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alphas,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_limb_t * p, * q;
    const fmpz * Acoeffs = A->coeffs;
    slong j, Alen = A->length;

    n_poly_fit_length(Acur, Alen);
    Acur->length = Alen;

    n_poly_fit_length(Ainc, 2*Alen);
    Ainc->length = 2*Alen;

    p = Acur->coeffs;
    _mpoly_monomial_evals_nmod(p, A->exps, A->bits, Alen, alphas, 0,
                                                      ctx->minfo, ctx_sp->mod);
    q = Ainc->coeffs;
    for (j = 0; j < Alen; j++)
    {
        q[2*j + 0] = fmpz_fdiv_ui(Acoeffs + j, ctx_sp->mod.n);
        q[2*j + 1] = nmod_pow_ui(p[j], pow, ctx_sp->mod);
    }
}


static void _worker_skel_sp(void * varg_)
{
    _base_struct * w = (_base_struct *) varg_;
    slong Alength = w->A->length;
    slong Blength = w->B->length;
    slong i;

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= Alength + Blength)
        goto cleanup;

    if (i < Alength)
    {
        w->Ainc_sp->terms[i].exp = w->A->exps[i];
        nmod_mpoly_get_eval_helper_pow(w->Aone_sp->coeffs + i,
                          w->Ainc_sp->terms[i].coeff, w->num_threads, w->ctx_sp,
                                       w->A->coeffs + i, w->alphas_sp, w->ctx);
    }
    else
    {
        i -= Alength;
        w->Binc_sp->terms[i].exp = w->B->exps[i];
        nmod_mpoly_get_eval_helper_pow(w->Bone_sp->coeffs + i,
                          w->Binc_sp->terms[i].coeff, w->num_threads, w->ctx_sp,
                                       w->B->coeffs + i, w->alphas_sp, w->ctx);
    }

    goto get_next_index;

cleanup:

    return;
}


static void _worker_skel_mp(void * varg_)
{
    _base_struct * w = (_base_struct *) varg_;
    slong Alength = w->A->length;
    slong Blength = w->B->length;
    slong i;

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= Alength + Blength)
        goto cleanup;

    if (i < Alength)
    {
        fmpz_mod_mpoly_set_skel(w->Aone_mp->coeffs + i,
                            w->ctx_mp, w->A->coeffs + i, w->alphas_mp, w->ctx);
        fmpz_mod_mpoly_red_skel(w->Ared_mp->coeffs + i,
                                          w->A->coeffs + i, w->ctx_mp->ffinfo);
        fmpz_mod_mpoly_pow_skel(w->Ainc_mp->coeffs + i,
                            w->Aone_mp->coeffs + i, w->num_threads, w->ctx_mp);
    }
    else
    {
        i -= Alength;
        fmpz_mod_mpoly_set_skel(w->Bone_mp->coeffs + i,
                            w->ctx_mp, w->B->coeffs + i, w->alphas_mp, w->ctx);
        fmpz_mod_mpoly_red_skel(w->Bred_mp->coeffs + i,
                                          w->B->coeffs + i, w->ctx_mp->ffinfo);
        fmpz_mod_mpoly_pow_skel(w->Binc_mp->coeffs + i,
                            w->Bone_mp->coeffs + i, w->num_threads, w->ctx_mp);
    }

    goto get_next_index;

cleanup:

    return;
}


/*
    evaluation is at x_i = w->alphas[i]

    Set {Gamma|A|B}one_sp to the evaluation of the monomials in {Gamma|A|B}.
    Set {Gamma|A|B}red_sp to the reduction mod p_sp of coeffs of {Gamma|A|B}.
    Then, set {Gamma|A|B}inc to the num_threads^th power of {Gamma|A|B}one_sp.
*/
static void _set_skels_sp(
    _base_struct * w,
    _eval_sp_worker_arg_struct * args,
    const thread_pool_handle * handles)
{
    slong i;

    n_bpoly_fit_length(w->Aone_sp, w->A->length);
    w->Aone_sp->length = w->A->length;
    n_polyun_fit_length(w->Ainc_sp, w->A->length);
    w->Ainc_sp->length = w->A->length;

    n_bpoly_fit_length(w->Bone_sp, w->B->length);
    w->Bone_sp->length = w->B->length;
    n_polyun_fit_length(w->Binc_sp, w->B->length);
    w->Binc_sp->length = w->B->length;

    w->index = 0;

    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, _worker_skel_sp, w);

    nmod_mpoly_get_eval_helper_pow(w->Gammaone_sp, w->Gammainc_sp,
                    w->num_threads, w->ctx_sp, w->Gamma, w->alphas_sp, w->ctx);

    _worker_skel_sp(w);

    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    /* signal to threads that they need to initialize cur = one^(i+1) */
    for (i = 0; i < w->num_threads; i++)
    {
        args[i].thread_index = i;
        args[i].cur_is_uninited = 1;
    }
}

static void _set_skels_mp(
    _base_struct * w,
    _eval_mp_worker_arg_struct * args,
    const thread_pool_handle * handles)
{
    slong i;

    fmpz_mpolycu_set_length(w->Aone_mp, w->A->length);
    fmpz_mpolycu_set_length(w->Ared_mp, w->A->length);
    fmpz_mpolycu_set_length(w->Ainc_mp, w->A->length);
    fmpz_mpolycu_set_length(w->Bone_mp, w->B->length);
    fmpz_mpolycu_set_length(w->Bred_mp, w->B->length);
    fmpz_mpolycu_set_length(w->Binc_mp, w->B->length);

    w->index = 0;
    for (i = 1; i < w->num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                                    _worker_skel_mp, w);
    }
    fmpz_mod_mpoly_set_skel(w->Gammaone_mp, w->ctx_mp, w->Gamma, w->alphas_mp, w->ctx);
    fmpz_mod_mpoly_red_skel(w->Gammared_mp, w->Gamma, w->ctx_mp->ffinfo);
    fmpz_mod_mpoly_pow_skel(w->Gammainc_mp, w->Gammaone_mp, w->num_threads, w->ctx_mp);
    _worker_skel_mp(w);
    for (i = 1; i < w->num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i - 1]);
    }

    /* signal to threads that they need to initialize cur = one^(i+1) */
    for (i = 0; i < w->num_threads; i++)
    {
        args[i].thread_index = i;
        args[i].cur_is_uninited = 1;
    }
}


static void n_poly_mod_scalar_pow(
    n_poly_t A,
    const n_poly_t B,
    ulong pow,
    nmod_t ctx)
{
    slong i;

    n_poly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        A->coeffs[i] = nmod_pow_ui(B->coeffs[i], pow, ctx);
    
}

void n_bpoly_mod_scalar_pow(
    n_bpoly_t A,
    const n_bpoly_t B,
    ulong pow,
    nmod_t ctx)
{
    slong i;

    n_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        n_poly_mod_scalar_pow(A->coeffs + i, B->coeffs + i, pow, ctx);
}

/*
    Put image gcd's on w->evals_{sp|mp}
*/
static void _worker_eval_sp(void * varg)
{
    _eval_sp_worker_arg_struct * arg = (_eval_sp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    mp_limb_t Gammaeval_sp;
    slong i;

    i = arg->thread_index;

    FLINT_ASSERT(i < w->num_threads);
    FLINT_ASSERT(w->num_images_sp >= w->num_threads);
    FLINT_ASSERT(0 == w->num_images_sp % w->num_threads);

    if (arg->cur_is_uninited)
    {
        n_poly_mod_scalar_pow(arg->Gammacur_sp, w->Gammaone_sp, i + 1, w->ctx_sp->mod);
        n_bpoly_mod_scalar_pow(arg->Acur_sp, w->Aone_sp, i + 1, w->ctx_sp->mod);
        n_bpoly_mod_scalar_pow(arg->Bcur_sp, w->Bone_sp, i + 1, w->ctx_sp->mod);
    }
    arg->cur_is_uninited = 0;

    while (i < w->num_images_sp)
    {
        _eval_sp_ret_struct * ret;

        Gammaeval_sp = n_poly_mod_eval_step2(arg->Gammacur_sp, w->Gammainc_sp, w->ctx_sp->mod);
        nmod_mpolyuu_eval_step2(arg->Aeval_sp, arg->Acur_sp, w->Ainc_sp, w->ctx_sp);
        nmod_mpolyuu_eval_step2(arg->Beval_sp, arg->Bcur_sp, w->Binc_sp, w->ctx_sp);

        ret = w->evals_sp + i;
        i += w->num_threads;

        if (arg->Aeval_sp->length < 1 || arg->Beval_sp->length < 1 ||
            n_bpoly_bidegree(arg->Aeval_sp) != w->A->exps[0] ||
            n_bpoly_bidegree(arg->Beval_sp) != w->B->exps[0])
        {
            ret->success = 0;
            continue;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(Gammaeval_sp != 0);

        ret->success = n_bpoly_mod_gcd_brown_smprime(arg->Geval_sp,
                             arg->Abareval_sp, arg->Bbareval_sp,
                     arg->Aeval_sp, arg->Beval_sp, w->ctx_sp->mod, arg->Sp_sp);
        if (!ret->success)
            continue;

        ret->GevaldegXY = n_bpoly_bidegree(arg->Geval_sp);
        ret->success = 1;
        n_bpoly_scalar_mul_nmod(arg->Geval_sp, Gammaeval_sp, w->ctx_sp->mod);
        n_bpoly_swap(ret->nGeval_sp, arg->Geval_sp);
        n_bpoly_swap(ret->nAbareval_sp, arg->Abareval_sp);
        n_bpoly_swap(ret->nBbareval_sp, arg->Bbareval_sp);
    }
}

static void _eval_mp_worker(void * varg)
{
    _eval_mp_worker_arg_struct * arg = (_eval_mp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    fmpz_t Gammaeval_mp;
    slong i;

    fmpz_init(Gammaeval_mp);

    i = arg->thread_index;

    FLINT_ASSERT(i < w->num_threads);
    FLINT_ASSERT(w->num_images_mp >= w->num_threads);
    FLINT_ASSERT(0 == w->num_images_mp % w->num_threads);

    if (arg->cur_is_uninited)
    {
        fmpz_mod_mpoly_pow_skel(arg->Gammacur_mp, w->Gammaone_mp, i + 1, w->ctx_mp);
        fmpz_mod_mpolyu_pow_skel(arg->Acur_mp, w->Aone_mp, i + 1, w->ctx_mp);
        fmpz_mod_mpolyu_pow_skel(arg->Bcur_mp, w->Bone_mp, i + 1, w->ctx_mp);
    }
    arg->cur_is_uninited = 0;

    while (i < w->num_images_mp)
    {
        _eval_mp_ret_struct * ret;

        fmpz_mod_mpoly_use_skel_mul(Gammaeval_mp, w->Gammared_mp,
                          arg->Gammacur_mp, w->Gammainc_mp, w->ctx_mp->ffinfo);
        fmpz_mod_mpolyuu_use_skel_mul(arg->Aeval_mp, w->A, w->Ared_mp,
                                          arg->Acur_mp, w->Ainc_mp, w->ctx_mp);
        fmpz_mod_mpolyuu_use_skel_mul(arg->Beval_mp, w->B, w->Bred_mp,
                                          arg->Bcur_mp, w->Binc_mp, w->ctx_mp);

        ret = w->evals_mp + i;
        i += w->num_threads;

        if (arg->Aeval_mp->length == 0 || arg->Beval_mp->length == 0
            || fmpz_mod_mpolyn_bidegree(arg->Aeval_mp) != w->A->exps[0]
            || fmpz_mod_mpolyn_bidegree(arg->Beval_mp) != w->B->exps[0])
        {
            /* evaluation killed either lc(A) or lc(B) */
            ret->success = 0;
            continue;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval_mp));
        ret->success = fmpz_mod_mpolyn_gcd_brown_bivar(arg->Geval_mp,
                                   arg->Abareval_mp, arg->Bbareval_mp,
                                      arg->Aeval_mp, arg->Beval_mp, w->ctx_mp);
        if (!ret->success)
        {
            /* image gcd could not be computed */
            continue;
        }

        /* leading coefficient correction */
        ret->GevaldegXY = fmpz_mod_mpolyn_bidegree(arg->Geval_mp);
        ret->success = 1;
        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(arg->Geval_mp, Gammaeval_mp,
                                                                    w->ctx_mp);
        fmpz_mod_mpolyn_swap(ret->nGeval_mp, arg->Geval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_swap(ret->nAbareval_mp, arg->Abareval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_swap(ret->nBbareval_mp, arg->Bbareval_mp, w->ctx_mp);
    }

    fmpz_clear(Gammaeval_mp);
}


/*
    Evaluate the coefficients of A, B, and H at alphas
*/
static void _worker_check_eval_sp(void * varg)
{
    _base_struct * w = (_base_struct *) varg;
    const slong Alength = w->A->length;
    const slong Blength = w->B->length;
    const slong Hlength = w->H->length;
    fmpz_mpoly_struct * Acoeffs = w->A->coeffs;
    fmpz_mpoly_struct * Bcoeffs = w->B->coeffs;
    fmpz_mpoly_struct * Hcoeffs = w->H->coeffs;
    mp_limb_t * Aevals = w->coeff_evals_sp->coeffs + 0;
    mp_limb_t * Bevals = Aevals + Alength;
    mp_limb_t * Hevals = Bevals + Blength;
    slong i;

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= Alength + Blength + Hlength)
        goto cleanup;

    if (i < Alength)
    {
        Aevals[i] = fmpz_mpoly_eval_nmod(w->ctx_sp->mod,
                                  Acoeffs + i, w->alphas_sp, w->ctx);
    }
    else
    {
        i -= Alength;
        if (i < Blength)
        {
            Bevals[i] = fmpz_mpoly_eval_nmod(w->ctx_sp->mod,
                                  Bcoeffs + i, w->alphas_sp, w->ctx);
        }
        else
        {
            i -= Blength;
            Hevals[i] = fmpz_mpoly_eval_nmod(w->ctx_sp->mod,
                                  Hcoeffs + i, w->alphas_sp, w->ctx);
        }
    }

    goto get_next_index;

cleanup:

    return;
}

static void _worker_check_eval_mp(void * varg)
{
    _base_struct * w = (_base_struct *) varg;
    const slong Alength = w->A->length;
    const slong Blength = w->B->length;
    const slong Hlength = w->H->length;
    fmpz_mpoly_struct * Acoeffs = w->A->coeffs;
    fmpz_mpoly_struct * Bcoeffs = w->B->coeffs;
    fmpz_mpoly_struct * Hcoeffs = w->H->coeffs;
    fmpz * Aevals = w->coeff_evals_mp->coeffs + 0;
    fmpz * Bevals = Aevals + Alength;
    fmpz * Hevals = Bevals + Blength;
    slong i;

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= Alength + Blength + Hlength)
        goto cleanup;

    if (i < Alength)
    {
        fmpz_mpoly_eval_fmpz_mod(Aevals + i, w->ctx_mp->ffinfo,
                                            Acoeffs + i, w->alphas_mp, w->ctx);
    }
    else
    {
        i -= Alength;
        if (i < Blength)
        {
            fmpz_mpoly_eval_fmpz_mod(Bevals + i, w->ctx_mp->ffinfo,
                                            Bcoeffs + i, w->alphas_mp, w->ctx);
        }
        else
        {
            i -= Blength;
            fmpz_mpoly_eval_fmpz_mod(Hevals + i, w->ctx_mp->ffinfo,
                                            Hcoeffs + i, w->alphas_mp, w->ctx);
        }
    }

    goto get_next_index;

cleanup:

    return;
}


static void _worker_reduce_sp(void * varg)
{
    _eval_sp_worker_arg_struct * arg = (_eval_sp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    nmod_berlekamp_massey_struct * GLcoeffs = w->GLambda_sp->coeffs;
    nmod_berlekamp_massey_struct * AbarLcoeffs = w->AbarLambda_sp->coeffs;
    nmod_berlekamp_massey_struct * BbarLcoeffs = w->BbarLambda_sp->coeffs;
/*
    int Gearly = w->GLambda_sp->pointcount/2 < Gamma->length;
    int Abarearly = w->AbarLambda_sp->pointcount/2 < w->A->coeffs[0].length;
    int Bbarearly = w->BbarLambda_sp->pointcount/2 < w->B->coeffs[0].length;
*/
    slong GLlen = w->GLambda_sp->length;
    slong AbarLlen = w->AbarLambda_sp->length;
    slong BbarLlen = w->BbarLambda_sp->length;
    slong i;

next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i < GLlen)
    {
        if (!w->Gchanged && nmod_berlekamp_massey_reduce(GLcoeffs + i))
            w->Gchanged = 1;
        goto next_index;
    }

    i -= GLlen;
    if (i < AbarLlen)
    {
        if (!w->Abarchanged && nmod_berlekamp_massey_reduce(AbarLcoeffs + i))
            w->Abarchanged = 1;
        goto next_index;
    }

    i -= AbarLlen;
    if (i < BbarLlen)
    {
        if (!w->Bbarchanged && nmod_berlekamp_massey_reduce(BbarLcoeffs + i))
            w->Bbarchanged = 1;
        goto next_index;
    }

    return;
}

static void _worker_reduce_mp(void * varg)
{
    _eval_mp_worker_arg_struct * arg = (_eval_mp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    fmpz_mod_berlekamp_massey_struct * GLcoeffs = w->GLambda_mp->coeffs;
    fmpz_mod_berlekamp_massey_struct * AbarLcoeffs = w->AbarLambda_mp->coeffs;
    fmpz_mod_berlekamp_massey_struct * BbarLcoeffs = w->BbarLambda_mp->coeffs;
/*
    int Gearly = w->GLambda_mp->pointcount/2 < Gamma->length;
    int Abarearly = w->AbarLambda_mp->pointcount/2 < w->A->coeffs[0].length;
    int Bbarearly = w->BbarLambda_mp->pointcount/2 < w->B->coeffs[0].length;
*/
    slong GLlen = w->GLambda_mp->length;
    slong AbarLlen = w->AbarLambda_mp->length;
    slong BbarLlen = w->BbarLambda_mp->length;
    slong i;

next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i < GLlen)
    {
        if (!w->Gchanged && fmpz_mod_berlekamp_massey_reduce(GLcoeffs + i, w->ctx_mp->ffinfo))
            w->Gchanged = 1;
        goto next_index;
    }

    i -= GLlen;
    if (i < AbarLlen)
    {
        if (!w->Abarchanged && fmpz_mod_berlekamp_massey_reduce(AbarLcoeffs + i, w->ctx_mp->ffinfo))
            w->Abarchanged = 1;
        goto next_index;
    }

    i -= AbarLlen;
    if (i < BbarLlen)
    {
        if (!w->Bbarchanged && fmpz_mod_berlekamp_massey_reduce(BbarLcoeffs + i, w->ctx_mp->ffinfo))
            w->Bbarchanged = 1;
        goto next_index;
    }

    return;
}


static void _worker_get_mpoly_sp(void * varg)
{
    _eval_sp_worker_arg_struct * arg = (_eval_sp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    nmod_berlekamp_massey_struct * Lcoeffs;
    ulong * Lexps;
    fmpz_mpoly_struct * Hcoeffs = w->H->coeffs;
    slong length = w->H->length;
    slong i;

    if (w->which_check == 1)
    {
        Lcoeffs = w->AbarLambda_sp->coeffs;
        Lexps = w->AbarLambda_sp->exps;
        FLINT_ASSERT(length == w->AbarLambda_sp->length);

    }
    else if (w->which_check == 2)
    {
        Lcoeffs = w->BbarLambda_sp->coeffs;
        Lexps = w->BbarLambda_sp->exps;
        FLINT_ASSERT(length == w->BbarLambda_sp->length);
    }
    else
    {
        FLINT_ASSERT(w->which_check == 0);
        Lcoeffs = w->GLambda_sp->coeffs;
        Lexps = w->GLambda_sp->exps;
        FLINT_ASSERT(length == w->GLambda_sp->length);
    }

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= length)
        return;

    w->H->exps[i] = Lexps[i];
    if (!w->failed)
    {
        if (!nmod_mpoly_bma_get_fmpz_mpoly(Hcoeffs + i, w->ctx,
                        w->alphashift, Lcoeffs + i, w->Ictx, w->ctx_sp->mod) ||
            w->H->coeffs[i].length == 0)
        {
            w->failed = 1;
            return;
        }
    }

    goto get_next_index;
}


static void _worker_get_mpoly_mp(void * varg)
{
    _eval_mp_worker_arg_struct * arg = (_eval_mp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    fmpz_mod_berlekamp_massey_struct * Lcoeffs;
    ulong * Lexps;
    fmpz_mpoly_struct * Hcoeffs = w->H->coeffs;
    slong i, length = w->H->length;

    if (w->which_check == 1)
    {
        Lcoeffs = w->AbarLambda_mp->coeffs;
        Lexps = w->AbarLambda_mp->exps;
        FLINT_ASSERT(length == w->AbarLambda_mp->length);

    }
    else if (w->which_check == 2)
    {
        Lcoeffs = w->BbarLambda_mp->coeffs;
        Lexps = w->BbarLambda_mp->exps;
        FLINT_ASSERT(length == w->BbarLambda_mp->length);
    }
    else
    {
        FLINT_ASSERT(w->which_check == 0);
        Lcoeffs = w->GLambda_mp->coeffs;
        Lexps = w->GLambda_mp->exps;
        FLINT_ASSERT(length == w->GLambda_mp->length);
    }

next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= length)
        return;

    w->H->exps[i] = Lexps[i];
    if (!w->failed)
    {
        if (!fmpz_mod_bma_get_fmpz_mpoly(Hcoeffs + i, w->ctx,
                 w->alphashift_mp, Lcoeffs + i, w->Ictx, w->ctx_mp->ffinfo) ||
            w->H->coeffs[i].length == 0)
        {
            w->failed = 1;
            return;
        }
    }

    goto next_index;

}


/*
    Try to solve for coefficients in zippel interpolation.
*/
static void _worker_find_zip_coeffs(void * varg)
{
    _eval_sp_worker_arg_struct * arg = (_eval_sp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    slong i;
    nmod_poly_t T;

    nmod_poly_init_mod(T, w->ctx_sp->mod);

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= w->Z->length)
        goto cleanup;

    switch (nmod_zip_find_coeffs(w->Z->coeffs + i, T,
                                             w->Z->pointcount, w->ctx_sp->mod))
    {
        case nmod_zip_find_coeffs_no_match:
            w->zip_find_coeffs_no_match = 1; /* safe - only goes from 0 to 1 */
            break;
        case nmod_zip_find_coeffs_non_invertible:
            w->zip_find_coeffs_non_invertible = 1; /* safe - only goes from 0 to 1 */
            break;
        default:
            break;
    }

    goto get_next_index;

cleanup:

    nmod_poly_clear(T);

    return;
}

/*
    Update H by crt'ing with a zippel image.
*/
static void _worker_crt_zip_coeffs(void * varg)
{
    _eval_sp_worker_arg_struct * arg = (_eval_sp_worker_arg_struct *) varg;
    _base_struct * w = arg->w;
    int changed = 0;
    slong i, j;
    fmpz_mpoly_struct * Hc;
    nmod_zip_struct * Zc;
    fmpz_t t;

    fmpz_init(t);

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&w->mutex);
#endif
    i = w->index;
    w->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&w->mutex);
#endif

    if (i >= w->H->length)
        goto cleanup;

    Hc = w->H->coeffs + i;
    Zc = w->Z->coeffs + i;

    FLINT_ASSERT(Hc->length == Zc->mlength);
    for (j = 0; j < Hc->length; j++)
    {
        fmpz_CRT_ui(t, Hc->coeffs + j, w->Hmodulus, Zc->coeffs[j],
                                                          w->ctx_sp->mod.n, 1);
        changed |= !fmpz_equal(t, Hc->coeffs + j);
        fmpz_swap(t, Hc->coeffs + j);
    }

    goto get_next_index;

cleanup:

    if (changed)
        w->Hchanged = 1;

    fmpz_clear(t);

    return;
}


/*
    Check divisibility.
*/
typedef struct
{
    fmpz_mpolyu_struct * quo;
    const fmpz_mpolyu_struct * num, * den;
    const fmpz_mpoly_ctx_struct * ctx;
    const thread_pool_handle * handles;
    slong num_handles;
    int success;
}
_divide_arg_struct;

typedef _divide_arg_struct _divide_arg_t[1];

static void _divide_worker(void * varg)
{
    _divide_arg_struct * arg = (_divide_arg_struct *) varg;

    if (arg->num_handles > 0)
    {
        arg->success = fmpz_mpolyuu_divides_threaded_pool(arg->quo, arg->num,
                        arg->den, 2, arg->ctx, arg->handles, arg->num_handles);
    }
    else
    {
        arg->success = fmpz_mpolyuu_divides(arg->quo, arg->num,
                                                        arg->den, 2, arg->ctx);
    }
}


/*
    Try to produce an image H mod p
*/
typedef enum {
    bma_loop_good,
    insufficient_eval_points,
    ksub_probably_unlucky,
    prime_probably_unlucky,
    gcd_is_one
} bma_loop_ret_t;

static bma_loop_ret_t _bma_loop_sp(
    ulong p_sp,
    _base_struct * w,
    _eval_sp_worker_arg_struct * args,
    const thread_pool_handle * handles)
{
    slong i, j;
    int unlucky_count, consecutive_unlucky_count, point_try_count;
    mp_limb_t cur_alpha_pow_sp;

    nmod_discrete_log_pohlig_hellman_precompute_prime(w->Ictx->dlogenv_sp, p_sp);
    nmod_mpoly_ctx_set_modulus(w->ctx_sp, p_sp);
    /* unfortunate nmod_poly's store their own ctx */
    _base_args_set_mod_sp(w, args);

    nmod_bma_mpoly_reset_prime(w->GLambda_sp, w->ctx_sp->mod);
    nmod_bma_mpoly_reset_prime(w->AbarLambda_sp, w->ctx_sp->mod);
    nmod_bma_mpoly_reset_prime(w->BbarLambda_sp, w->ctx_sp->mod);
    nmod_bma_mpoly_zero(w->GLambda_sp);
    nmod_bma_mpoly_zero(w->AbarLambda_sp);
    nmod_bma_mpoly_zero(w->BbarLambda_sp);

    /*
        Set w->alphas_sp[i] to evaluation of x_i under the ksub
            and the alpha chosen as a generator of F_p*.
    */
    nmod_mpoly_bma_interpolate_alpha_powers(w->alphas_sp,
                                           1, w->Ictx, w->ctx, w->ctx_sp->mod);

    /* set skeletons for evaluation */
    _set_skels_sp(w, args, handles);

    unlucky_count = 0;
    consecutive_unlucky_count = 0;
    cur_alpha_pow_sp = 0;

next_bma_image_sp:

    /* determine how many powers to evaluate */
    j = w->bma_target_count > w->GLambda_sp->pointcount
           ? w->bma_target_count - w->GLambda_sp->pointcount
           : w->num_threads;
    j = FLINT_MAX(j, WORD(2)); /* updating with 1 point is not interesting */
    j = (j + w->num_threads - 1)/w->num_threads;
    j = FLINT_MAX(j, WORD(2));
    j = FLINT_MIN(j, WORD(20));
    j *= w->num_threads;
    _base_set_num_images_sp(w, j);

    /*
        must evaluate at powers (where i = cur_alpha_pow_sp)
        alpha ^ (i + 1), ..., alpha ^ (i + w->num_images)
    */

    cur_alpha_pow_sp += w->num_images_sp;

    if (cur_alpha_pow_sp < w->num_images_sp || cur_alpha_pow_sp >= p_sp - 1)
    {
        return insufficient_eval_points;
    }

    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                                    _worker_eval_sp, &args[i]);
    _worker_eval_sp(&args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    /* if any image gcd failed, reset lamda */
    for (i = 0; i < w->num_images_sp; i++)
    {
        if (!w->evals_sp[i].success)
        {
            nmod_bma_mpoly_zero(w->GLambda_sp);
            nmod_bma_mpoly_zero(w->AbarLambda_sp);
            nmod_bma_mpoly_zero(w->BbarLambda_sp);
            consecutive_unlucky_count = 0;
            goto next_bma_image_sp;
        }
    }

    /* if any image gcd was unlucky, take action */
    for (i = 0; i < w->num_images_sp; i++)
    {
        if (w->GdegboundXY < w->evals_sp[i].GevaldegXY)
        {
            ++consecutive_unlucky_count;
            if (consecutive_unlucky_count > 1)
                return ksub_probably_unlucky;

            ++unlucky_count;
            if (unlucky_count > 2)
                return prime_probably_unlucky;

            nmod_bma_mpoly_zero(w->GLambda_sp);
            nmod_bma_mpoly_zero(w->AbarLambda_sp);
            nmod_bma_mpoly_zero(w->BbarLambda_sp);
            for (i++; i < w->num_images_sp; i++)
            {
                if (w->GdegboundXY < w->evals_sp[i].GevaldegXY)
                    ++consecutive_unlucky_count;
                else
                    consecutive_unlucky_count = 0;
            }
            goto next_bma_image_sp;
        }
        else
        {
            consecutive_unlucky_count = 0;
        }
    }

    /* if any image gcd was revealing, reset lambda or finish */
    for (i = 0; i < w->num_images_sp; i++)
    {
        if (w->GdegboundXY > w->evals_sp[i].GevaldegXY)
        {
            /* new bound on deg_XY(G) */
            nmod_bma_mpoly_zero(w->GLambda_sp);
            nmod_bma_mpoly_zero(w->AbarLambda_sp);
            nmod_bma_mpoly_zero(w->BbarLambda_sp);
            consecutive_unlucky_count = 0;
            w->GdegboundXY = w->evals_sp[i].GevaldegXY;
            if (w->GdegboundXY == 0)
                return gcd_is_one;
            goto next_bma_image_sp;
        }
    }

    for (i = 0; i < w->num_images_sp; i++)
    {
        nmod_bma_mpoly_add_point(w->GLambda_sp, w->evals_sp[i].nGeval_sp, w->ctx_sp);
        nmod_bma_mpoly_add_point(w->AbarLambda_sp, w->evals_sp[i].nAbareval_sp, w->ctx_sp);
        nmod_bma_mpoly_add_point(w->BbarLambda_sp, w->evals_sp[i].nBbareval_sp, w->ctx_sp);
    }

    if (w->bma_target_count > w->GLambda_sp->pointcount)
    {
        goto next_bma_image_sp;
    }

    w->Gchanged = w->GLambda_sp->pointcount/2 < w->Gamma->length;
    w->Abarchanged = w->AbarLambda_sp->pointcount/2 < w->A->coeffs[0].length;
    w->Bbarchanged = w->BbarLambda_sp->pointcount/2 < w->B->coeffs[0].length;

    if (w->Gchanged && w->Abarchanged && w->Bbarchanged)
    {
        goto next_bma_image_sp;
    }

    w->index = 0;
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                                  _worker_reduce_sp, &args[i]);
    _worker_reduce_sp(&args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (!w->Gchanged)
    {
        w->which_check = 0;
        fmpz_mpolyu_fit_length(w->H, w->GLambda_sp->length, w->ctx);
        w->H->length = w->GLambda_sp->length;
    }
    else if (!w->Abarchanged)
    {
        w->which_check = 1;
        fmpz_mpolyu_fit_length(w->H, w->AbarLambda_sp->length, w->ctx);
        w->H->length = w->AbarLambda_sp->length;
    }
    else if (!w->Bbarchanged)
    {
        w->which_check = 2;
        fmpz_mpolyu_fit_length(w->H, w->BbarLambda_sp->length, w->ctx);
        w->H->length = w->BbarLambda_sp->length;
    }
    else
    {
        w->which_check = -1;
        goto next_bma_image_sp;
    }

    FLINT_ASSERT(cur_alpha_pow_sp >= w->GLambda_sp->pointcount);
    FLINT_ASSERT(w->AbarLambda_sp->pointcount >= w->GLambda_sp->pointcount);
    FLINT_ASSERT(w->BbarLambda_sp->pointcount >= w->GLambda_sp->pointcount);

    w->index = 0;
    w->failed = 0;
    w->alphashift = cur_alpha_pow_sp - w->GLambda_sp->pointcount + 1;
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                              _worker_get_mpoly_sp, &args[i]);
    _worker_get_mpoly_sp(&args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (w->failed || w->H->coeffs[0].length != 
                   (w->which_check == 1 ? w->A->coeffs + 0 :
                    w->which_check == 2 ? w->B->coeffs + 0 : w->Gamma)->length)
    {
        goto next_bma_image_sp;
    }

    FLINT_ASSERT(w->which_check != 0 || w->GdegboundXY == w->H->exps[0]);
    FLINT_ASSERT(w->which_check != 1 || w->GdegboundXY == w->A->exps[0] - w->H->exps[0]);
    FLINT_ASSERT(w->which_check != 2 || w->GdegboundXY == w->B->exps[0] - w->H->exps[0]);

    /* coeffs_evals_sp will store the the evaluations of coeffs of A, B, H */
    i = w->A->length + w->B->length + w->H->length;
    n_poly_fit_length(w->coeff_evals_sp, i);
    w->coeff_evals_sp->length = i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        mp_limb_t * Aevals, * Bevals, * Hevals;
        ulong GevaldegXY;
        mp_limb_t Gammaeval_sp;
        int success;

        for (i = 0; i < w->ctx->minfo->nvars; i++)
            w->alphas_sp[i] = n_urandint(w->randstate, w->ctx_sp->mod.n);

        /* fill in {A|B|H}evals[i] with the evaluation of {A|B|H}->coeffs + i */
        Aevals = w->coeff_evals_sp->coeffs + 0;
        Bevals = Aevals + w->A->length;
        Hevals = Bevals + w->B->length;
        w->index = 0;
        for (i = 1; i < w->num_threads; i++)
            thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                                        _worker_check_eval_sp, w);

        Gammaeval_sp = fmpz_mpoly_eval_nmod(w->ctx_sp->mod, w->Gamma,
                                                         w->alphas_sp, w->ctx);
        _worker_check_eval_sp(w);
        for (i = 1; i < w->num_threads; i++)
            thread_pool_wait(global_thread_pool, handles[i - 1]);

        _fmpz_mpolyuu_eval_nmod_from_coeffs(args->Aeval_sp, w->ctx_sp,
                                                         w->A, w->ctx, Aevals);
        _fmpz_mpolyuu_eval_nmod_from_coeffs(args->Beval_sp, w->ctx_sp,
                                                         w->B, w->ctx, Bevals);

        /* make sure that evaluation did not kill either lc(A) or lc(B) */
        if (args->Aeval_sp->length == 0 || args->Beval_sp->length == 0 ||
            n_bpoly_bidegree(args->Aeval_sp) != w->A->exps[0] ||
            n_bpoly_bidegree(args->Beval_sp) != w->B->exps[0])
        {
            continue;
        }

        /* Gamma is gcd(lc(A), lc(B)) so it evaluation should not be zero */
        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_bpoly_mod_gcd_brown_smprime(args->Geval_sp,
                       args->Abareval_sp, args->Bbareval_sp, args->Aeval_sp,
                                  args->Beval_sp, w->ctx_sp->mod, args->Sp_sp);
        if (!success)
            continue;

        n_bpoly_scalar_mul_nmod(args->Geval_sp, Gammaeval_sp, w->ctx_sp->mod);
        GevaldegXY = n_bpoly_bidegree(args->Geval_sp);

        if (w->GdegboundXY < GevaldegXY)
        {
            goto next_bma_image_sp;
        }
        else if (w->GdegboundXY > GevaldegXY)
        {
            w->GdegboundXY = GevaldegXY;
            nmod_bma_mpoly_zero(w->GLambda_sp);
            nmod_bma_mpoly_zero(w->AbarLambda_sp);
            nmod_bma_mpoly_zero(w->BbarLambda_sp);
            if (w->GdegboundXY == 0)
                return gcd_is_one;
            goto next_bma_image_sp;
        }

        if (w->which_check == 0)
        {
            _fmpz_mpolyuu_eval_nmod_from_coeffs(args->Bbareval_sp, w->ctx_sp,
                                                         w->H, w->ctx, Hevals);
            if (!n_bpoly_equal(args->Bbareval_sp, args->Geval_sp))
                goto next_bma_image_sp;
        }
        else
        {
            _fmpz_mpolyuu_eval_nmod_from_coeffs(args->Geval_sp, w->ctx_sp,
                                                         w->H, w->ctx, Hevals);
            if (!n_bpoly_equal(args->Geval_sp,
                  w->which_check == 1 ? args->Abareval_sp : args->Bbareval_sp))
            {
                goto next_bma_image_sp;
            }
        }

        /* if match, then assume H is good */
        break;
    }
    /* if no good evaluation point was found, then assume H is good */

    /* assume that H is correct modulo Hmodulus = p */
    fmpz_set_ui(w->Hmodulus, p_sp);

    return bma_loop_good;
}

static bma_loop_ret_t _bma_loop_mp(
    const fmpz_t p,
    _base_struct * w,
    _eval_mp_worker_arg_struct * args,
    const thread_pool_handle * handles)
{
    bma_loop_ret_t ret;
    slong i, j;
    int success;
    int unlucky_count, consecutive_unlucky_count, point_try_count;
    ulong GevaldegXY;
    fmpz_t pminus1;
    fmpz_t Gammaeval_mp;
    fmpz_t cur_alpha_pow_mp;

    fmpz_init(cur_alpha_pow_mp);
    fmpz_init(Gammaeval_mp);
    fmpz_init(pminus1);
    fmpz_sub_ui(pminus1, p, 1);

    fmpz_mod_ctx_set_modulus(w->ctx_mp->ffinfo, p);
    fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(w->Ictx->dlogenv, p);
    /* unfortunate fmpz_mod_poly's store their own ctx */
    _base_args_set_mod_mp(w, args);

    fmpz_mod_bma_mpoly_zero(w->GLambda_mp);
    fmpz_mod_bma_mpoly_zero(w->AbarLambda_mp);
    fmpz_mod_bma_mpoly_zero(w->BbarLambda_mp);

    /*
        Set w->alphas[i] to evaluation of x_i under the ksub
            and the alpha chosen as a generator of F_p*.
    */

    fmpz_one(cur_alpha_pow_mp);
    fmpz_mod_mpoly_bma_interpolate_alpha_powers(w->alphas_mp,
                         cur_alpha_pow_mp, w->Ictx, w->ctx, w->ctx_mp->ffinfo);

    /* set skeletons for evaluation */
    _set_skels_mp(w, args, handles);

    unlucky_count = 0;
    consecutive_unlucky_count = 0;
    fmpz_zero(cur_alpha_pow_mp);

next_bma_image_mp:

    /* determine how many powers to evaluate */
    j = w->bma_target_count > w->GLambda_mp->pointcount
           ? w->bma_target_count - w->GLambda_mp->pointcount
           : w->num_threads;
    j = FLINT_MAX(j, WORD(2)); /* updating with 1 point is not interesting */
    j = (j + w->num_threads - 1)/w->num_threads;
    j = FLINT_MAX(j, WORD(2));
    j = FLINT_MIN(j, WORD(20));
    j *= w->num_threads;
    _base_set_num_images_mp(w, j);

    /*
        must evaluate at powers (where i = cur_alpha_pow_sp)
        alpha ^ (i + 1), ..., alpha ^ (i + w->num_images)
    */

    fmpz_add_ui(cur_alpha_pow_mp, cur_alpha_pow_mp, w->num_images_mp);
    if (fmpz_cmp(cur_alpha_pow_mp, pminus1) >= 0)
    {
        ret = insufficient_eval_points;
        goto done;
    }

    for (i = 1; i < w->num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                        _eval_mp_worker, &args[i]);
    }
    _eval_mp_worker(&args[0]);
    for (i = 1; i < w->num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i - 1]);
    }

    /* if any image gcd failed, reset lamda */
    for (i = 0; i < w->num_images_mp; i++)
    {
        if (!w->evals_mp[i].success)
        {
            consecutive_unlucky_count = 0;
            fmpz_mod_bma_mpoly_zero(w->GLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->BbarLambda_mp);
            goto next_bma_image_mp;
        }
    }

    /* if any image gcd was unlucky, take action */
    for (i = 0; i < w->num_images_mp; i++)
    {
        if (w->GdegboundXY < w->evals_mp[i].GevaldegXY)
        {
            for (; i < w->num_images_mp; i++)
            {
                if (w->GdegboundXY < w->evals_mp[i].GevaldegXY)
                {
                    if (++consecutive_unlucky_count > 1)
                    {
                        ret = ksub_probably_unlucky;
                        goto done;
                    }
                    if (++unlucky_count > 2)
                    {
                        ret = prime_probably_unlucky;
                        goto done;
                    }
                }
                else
                {
                    consecutive_unlucky_count = 0;
                }
            }
            fmpz_mod_bma_mpoly_zero(w->GLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->BbarLambda_mp);
            goto next_bma_image_mp;
        }
        consecutive_unlucky_count = 0;
    }

    /* if any image gcd was revealing, reset lambda or finish */
    for (i = 0; i < w->num_images_mp; i++)
    {
        if (w->GdegboundXY > w->evals_mp[i].GevaldegXY)
        {
            /* new bound on deg_XY(G) */
            fmpz_mod_bma_mpoly_zero(w->GLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->BbarLambda_mp);
            consecutive_unlucky_count = 0;
            w->GdegboundXY = w->evals_mp[i].GevaldegXY;
            if (w->GdegboundXY == 0)
            {
                ret = gcd_is_one;
                goto done;
            }
            goto next_bma_image_mp;
        }
    }

    for (i = 0; i < w->num_images_mp; i++)
    {
        fmpz_mod_bma_mpoly_add_point(w->GLambda_mp, w->evals_mp[i].nGeval_mp, w->ctx_mp);
        fmpz_mod_bma_mpoly_add_point(w->AbarLambda_mp, w->evals_mp[i].nAbareval_mp, w->ctx_mp);
        fmpz_mod_bma_mpoly_add_point(w->BbarLambda_mp, w->evals_mp[i].nBbareval_mp, w->ctx_mp);
    }

    if (w->bma_target_count > w->GLambda_mp->pointcount)
    {
        goto next_bma_image_mp;
    }

    w->Gchanged = w->GLambda_mp->pointcount/2 < w->Gamma->length;
    w->Abarchanged = w->AbarLambda_mp->pointcount/2 < w->A->coeffs[0].length;
    w->Bbarchanged = w->BbarLambda_mp->pointcount/2 < w->B->coeffs[0].length;

    if (w->Gchanged && w->Abarchanged && w->Bbarchanged)
    {
        goto next_bma_image_mp;
    }

    w->index = 0;
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                              _worker_reduce_mp, &args[i]);
    _worker_reduce_mp(&args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (!w->Gchanged)
    {
        w->which_check = 0;
        fmpz_mpolyu_fit_length(w->H, w->GLambda_mp->length, w->ctx);
        w->H->length = w->GLambda_mp->length;
    }
    else if (!w->Abarchanged)
    {
        w->which_check = 1;
        fmpz_mpolyu_fit_length(w->H, w->AbarLambda_mp->length, w->ctx);
        w->H->length = w->AbarLambda_mp->length;
    }
    else if (!w->Bbarchanged)
    {
        w->which_check = 2;
        fmpz_mpolyu_fit_length(w->H, w->BbarLambda_mp->length, w->ctx);
        w->H->length = w->BbarLambda_mp->length;
    }
    else
    {
        w->which_check = -1;
        goto next_bma_image_mp;
    }

    FLINT_ASSERT(w->AbarLambda_mp->pointcount >= w->GLambda_mp->pointcount);
    FLINT_ASSERT(w->BbarLambda_mp->pointcount >= w->GLambda_mp->pointcount);

    /* get fmpz_mpolyu */
    w->index = 0;
    w->failed = 0;
    fmpz_sub_ui(w->alphashift_mp, cur_alpha_pow_mp, w->GLambda_mp->pointcount);
    FLINT_ASSERT(fmpz_sgn(w->alphashift_mp) >= 0);
    fmpz_add_ui(w->alphashift_mp, w->alphashift_mp, 1);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                              _worker_get_mpoly_mp, &args[i]);
    _worker_get_mpoly_mp(&args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (w->failed || w->H->coeffs[0].length != 
                   (w->which_check == 1 ? w->A->coeffs + 0 :
                    w->which_check == 2 ? w->B->coeffs + 0 : w->Gamma)->length)
    {
        goto next_bma_image_mp;
    }

    FLINT_ASSERT(w->which_check != 0 || w->GdegboundXY == w->H->exps[0]);
    FLINT_ASSERT(w->which_check != 1 || w->GdegboundXY == w->A->exps[0] - w->H->exps[0]);
    FLINT_ASSERT(w->which_check != 2 || w->GdegboundXY == w->B->exps[0] - w->H->exps[0]);

    /* coeffs_evals_mp will store the the evaluations of coeffs of A, B, H */
    i = w->A->length + w->B->length + w->H->length;
    fmpz_mpolyc_fit_length(w->coeff_evals_mp, i);
    w->coeff_evals_mp->length = i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        fmpz * Aevals, * Bevals, * Hevals;

        /* evaluate Gamma, A, and B at random point */
        for (i = 0; i < w->ctx->minfo->nvars; i++)
            fmpz_randm(w->alphas_mp + i, w->randstate,
                                      fmpz_mod_ctx_modulus(w->ctx_mp->ffinfo));

        /* fill in {A|B|H}evals[i] with the evaluation of {A|B|H}->coeffs + i */
        Aevals = w->coeff_evals_mp->coeffs + 0;
        Bevals = Aevals + w->A->length;
        Hevals = Bevals + w->B->length;
        w->index = 0;
        for (i = 1; i < w->num_threads; i++)
            thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                                     _worker_check_eval_mp, w);
        fmpz_mpoly_eval_fmpz_mod(Gammaeval_mp, w->ctx_mp->ffinfo,
                                               w->Gamma, w->alphas_mp, w->ctx);
        _worker_check_eval_mp(w);
        for (i = 1; i < w->num_threads; i++)
            thread_pool_wait(global_thread_pool, handles[i - 1]);
        _fmpz_mpolyuu_eval_fmpz_mod_from_coeffs(args->Aeval_mp, w->ctx_mp,
                                                         w->A, w->ctx, Aevals);
        _fmpz_mpolyuu_eval_fmpz_mod_from_coeffs(args->Beval_mp, w->ctx_mp,
                                                         w->B, w->ctx, Bevals);

        /* make sure that evaluation did not kill either lc(A) or lc(B) */
        if (   args->Aeval_mp->length == 0 || args->Beval_mp->length == 0 
            || fmpz_mod_mpolyn_bidegree(args->Aeval_mp) != w->A->exps[0]
            || fmpz_mod_mpolyn_bidegree(args->Beval_mp) != w->B->exps[0])
        {
            continue;
        }

        /* Gamma is gcd(lc(A), lc(B)) so it evaluation should not be zero */
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval_mp));

        success = fmpz_mod_mpolyn_gcd_brown_bivar(args->Geval_mp,
                            args->Abareval_mp, args->Bbareval_mp,
                            args->Aeval_mp, args->Beval_mp, w->ctx_mp);
        if (!success)
        {
            continue;
        }

        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(args->Geval_mp, Gammaeval_mp,
                                                                    w->ctx_mp);
        GevaldegXY = fmpz_mod_mpolyn_bidegree(args->Geval_mp);

        if (w->GdegboundXY < GevaldegXY)
        {
            /* the random evaluation point was unlucky */
            goto next_bma_image_mp;
        }
        else if (w->GdegboundXY > GevaldegXY)
        {
            /* the random evaluation point gave us a better degree bound */
            w->GdegboundXY = GevaldegXY;
            fmpz_mod_bma_mpoly_zero(w->GLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(w->BbarLambda_mp);
            if (w->GdegboundXY == 0)
            {
                ret = gcd_is_one;
                goto done;
            }
            goto next_bma_image_mp;
        }
#if 0
        /* reuse Bbareval for Heval */
        _fmpz_mpolyuu_eval_fmpz_mod_from_coeffs(args->Bbareval_mp, w->ctx_mp,
                                                         w->H, w->ctx, Hevals);
        if (!fmpz_mod_mpolyn_equal(args->Bbareval_mp, args->Geval_mp, w->ctx_mp))
        {
            goto next_bma_image_mp;
        }
#endif

        if (w->which_check == 0)
        {
            _fmpz_mpolyuu_eval_fmpz_mod_from_coeffs(args->Bbareval_mp, w->ctx_mp,
                                                         w->H, w->ctx, Hevals);
            if (!fmpz_mod_mpolyn_equal(args->Bbareval_mp, args->Geval_mp, w->ctx_mp))
                goto next_bma_image_mp;
        }
        else
        {
            _fmpz_mpolyuu_eval_fmpz_mod_from_coeffs(args->Geval_mp, w->ctx_mp,
                                                         w->H, w->ctx, Hevals);
            if (!fmpz_mod_mpolyn_equal(args->Geval_mp,
                  w->which_check == 1 ? args->Abareval_mp : args->Bbareval_mp,
                                                                    w->ctx_mp))
            {
                goto next_bma_image_mp;
            }
        }

        /* if match, then assume H is good */
        break;
    }
    /* if no good evaluation point was found, then assume H is good */

    /* assume that H is correct modulo Hmodulus = p */
    fmpz_set(w->Hmodulus, p);
    ret = bma_loop_good;

done:

    fmpz_clear(cur_alpha_pow_mp);
    fmpz_clear(Gammaeval_mp);
    fmpz_clear(pminus1);

    return ret;
}


int fmpz_mpolyuu_gcd_berlekamp_massey_threaded_pool(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t Abar,
    fmpz_mpolyu_t Bbar,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles, slong num_handles)
{
    int success, point_try_count;
    fmpz_mpoly_t Hcontent;
    _base_t w;
    _eval_sp_worker_arg_struct * eval_sp_args;
    _eval_mp_worker_arg_struct * eval_mp_args;
    slong i, j;
    fmpz_t p, subprod, cAksub, cBksub;
    mp_limb_t p_sp;
    slong zip_evals;
    ulong ABtotal_length;

    w->bits = A->bits;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(w->bits == A->bits);
    FLINT_ASSERT(w->bits == B->bits);
    FLINT_ASSERT(w->bits == G->bits);
    FLINT_ASSERT(w->bits == Abar->bits);
    FLINT_ASSERT(w->bits == Bbar->bits);
    FLINT_ASSERT(w->bits == Gamma->bits);

    /* let's initialize everything at once to avoid complicated cleanup */

    w->A = A;
    w->B = B;
    w->Gamma = Gamma;
    w->Abar = Abar;
    w->Bbar = Bbar;
    w->ctx = ctx;
    w->which_check = -1;

    ABtotal_length = 0;
    for (i = 0; i < A->length; i++)
        ABtotal_length += (A->coeffs + i)->length;
    for (i = 0; i < B->length; i++)
        ABtotal_length += (B->coeffs + i)->length;

    ABtotal_length = FLINT_MAX(ABtotal_length, 100);

    flint_randinit(w->randstate);
    fmpz_init(p);
    fmpz_init(subprod);
    fmpz_init(cAksub);
    fmpz_init(cBksub);
    fmpz_init(w->Hmodulus);
    fmpz_mpolyu_init(w->H, w->bits, ctx);
    fmpz_mpoly_init3(Hcontent, 0, w->bits, ctx);

    mpoly_bma_interpolate_ctx_init(w->Ictx, ctx->minfo->nvars);

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&w->mutex, NULL);
#endif

    w->num_threads = num_handles + 1;

    /* multiprecision workspace */

    w->evals_mp = NULL;
    w->evals_mp_alloc = 0;
    fmpz_mpolyc_init(w->coeff_evals_mp);
    fmpz_init(w->alphashift_mp);

    fmpz_set_ui(p, 2);    /* modulus no care */
    fmpz_mod_mpoly_ctx_init(w->ctx_mp, 2, ORD_LEX, p); /* modulus no care */
    fmpz_mod_bma_mpoly_init(w->GLambda_mp);
    fmpz_mod_bma_mpoly_init(w->AbarLambda_mp);
    fmpz_mod_bma_mpoly_init(w->BbarLambda_mp);

    eval_mp_args = (_eval_mp_worker_arg_struct *) flint_malloc(
                            w->num_threads*sizeof(_eval_mp_worker_arg_struct));
    for (i = 0; i < w->num_threads; i++)
    {
        eval_mp_args[i].w = w;
        fmpz_mod_mpolyn_init(eval_mp_args[i].Aeval_mp, FLINT_BITS/2, w->ctx_mp);
        fmpz_mod_mpolyn_init(eval_mp_args[i].Beval_mp, FLINT_BITS/2, w->ctx_mp);
        fmpz_mod_mpolyn_init(eval_mp_args[i].Geval_mp, FLINT_BITS/2, w->ctx_mp);
        fmpz_mod_mpolyn_init(eval_mp_args[i].Abareval_mp, FLINT_BITS/2, w->ctx_mp);
        fmpz_mod_mpolyn_init(eval_mp_args[i].Bbareval_mp, FLINT_BITS/2, w->ctx_mp);
        fmpz_mpolycu_init(eval_mp_args[i].Acur_mp);
        fmpz_mpolycu_init(eval_mp_args[i].Bcur_mp);
        fmpz_mpolyc_init(eval_mp_args[i].Gammacur_mp);
    }

    fmpz_mpolyc_init(w->Gammaone_mp);
    fmpz_mpolycu_init(w->Aone_mp);
    fmpz_mpolycu_init(w->Bone_mp);
    fmpz_mpolyc_init(w->Gammainc_mp);
    fmpz_mpolycu_init(w->Ainc_mp);
    fmpz_mpolycu_init(w->Binc_mp);
    fmpz_mpolyc_init(w->Gammared_mp);
    fmpz_mpolycu_init(w->Ared_mp);
    fmpz_mpolycu_init(w->Bred_mp);
    w->alphas_mp = _fmpz_vec_init(ctx->minfo->nvars);

    /* machine precision workspace */

    w->evals_sp = NULL;
    w->evals_sp_alloc = 0;
    n_poly_init(w->coeff_evals_sp);

    nmod_mpoly_ctx_init(w->ctx_sp, 2, ORD_LEX, 2); /* modulus no care */
    nmod_bma_mpoly_init(w->GLambda_sp);
    nmod_bma_mpoly_init(w->AbarLambda_sp);
    nmod_bma_mpoly_init(w->BbarLambda_sp);

    eval_sp_args = (_eval_sp_worker_arg_struct *) flint_malloc(
                         w->num_threads*sizeof(_eval_sp_worker_arg_struct));
    for (i = 0; i < w->num_threads; i++)
    {
        eval_sp_args[i].w = w;
        n_bpoly_init(eval_sp_args[i].Aeval_sp);
        n_bpoly_init(eval_sp_args[i].Beval_sp);
        n_bpoly_init(eval_sp_args[i].Geval_sp);
        n_bpoly_init(eval_sp_args[i].Abareval_sp);
        n_bpoly_init(eval_sp_args[i].Bbareval_sp);
        n_bpoly_init(eval_sp_args[i].Acur_sp);
        n_bpoly_init(eval_sp_args[i].Bcur_sp);
        n_poly_init(eval_sp_args[i].Gammacur_sp);
        n_poly_stack_init(eval_sp_args[i].Sp_sp->poly_stack);
        n_bpoly_stack_init(eval_sp_args[i].Sp_sp->bpoly_stack);
    }

    n_bpoly_init(w->Aone_sp);
    n_polyun_init(w->Ainc_sp);
    n_bpoly_init(w->Bone_sp);
    n_polyun_init(w->Binc_sp);
    n_poly_init(w->Gammaone_sp);
    n_poly_init(w->Gammainc_sp);

    w->alphas_sp = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));

    nmod_zip_mpolyu_init(w->Z);

    w->Adegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    w->Bdegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    w->GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
    p_sp = UWORD(1) << (FLINT_BITS - 2);

    w->H->length = 0;
    i = w->A->length + w->B->length + w->H->length;
    n_poly_fit_length(w->coeff_evals_sp, i);
    w->coeff_evals_sp->length = i;

    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        _eval_sp_worker_arg_struct * arg = eval_sp_args + 0;
        mp_limb_t * Aevals = w->coeff_evals_sp->coeffs + 0;
        mp_limb_t * Bevals = Aevals + w->A->length;

        p_sp = n_nextprime(p_sp, 1);
        nmod_mpoly_ctx_set_modulus(w->ctx_sp, p_sp);
        /* unfortunate nmod_poly's need mod set */
        _base_args_set_mod_sp(w, eval_sp_args);

        for (i = 0; i < ctx->minfo->nvars; i++)
            w->alphas_sp[i] = n_urandint(w->randstate, p_sp);

        /* fill in {A|B}evals[i] with the evaluation of {A|B}->coeffs + i */
        w->index = 0;
        for (i = 1; i < w->num_threads; i++)
            thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                                        _worker_check_eval_sp, w);
        _worker_check_eval_sp(w);
        for (i = 1; i < w->num_threads; i++)
            thread_pool_wait(global_thread_pool, handles[i - 1]);

        _fmpz_mpolyuu_eval_nmod_from_coeffs(arg->Aeval_sp, w->ctx_sp,
                                                         w->A, w->ctx, Aevals);
        _fmpz_mpolyuu_eval_nmod_from_coeffs(arg->Aeval_sp, w->ctx_sp,
                                                         w->B, w->ctx, Bevals);

        if (arg->Aeval_sp->length == 0 || arg->Beval_sp->length == 0 ||
            n_bpoly_bidegree(arg->Aeval_sp) != A->exps[0] ||
            n_bpoly_bidegree(arg->Beval_sp) != B->exps[0])
        {
            continue;
        }

        success = n_bpoly_mod_gcd_brown_smprime(arg->Geval_sp,
                                     arg->Abareval_sp, arg->Bbareval_sp,
                     arg->Aeval_sp, arg->Beval_sp, w->ctx_sp->mod, arg->Sp_sp);
        if (success)
        {
            w->GdegboundXY = n_bpoly_bidegree(arg->Geval_sp);
            break;
        }
    }

    fmpz_mpolyu_inner_degrees_si(w->Adegs, A, ctx);
    fmpz_mpolyu_inner_degrees_si(w->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
        w->Ictx->degbounds[i] = 1 + FLINT_MAX(w->Adegs[i], w->Bdegs[i]);

    /* initialization done! */

    if (w->GdegboundXY == 0)
    {
        fmpz_mpolyu_one(G, ctx);
        fmpz_mpolyu_swap(Abar, A, ctx);
        fmpz_mpolyu_swap(Bbar, B, ctx);
        success = 1;
        goto cleanup;
    }

    /* find a image_count before which we do not try to reduce */
    w->bma_target_count = FLINT_MAX(w->num_threads, 2*Gamma->length);

    for (i = 0; i < ctx->minfo->nvars; i++)
        w->Ictx->subdegs[i] = w->Ictx->degbounds[i];

    goto got_ksub;

pick_ksub:

    if (ctx->minfo->nvars > 1)
    {
        /* just increment the smallest subdegs[j] */
        j = 1;
        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            if (w->Ictx->subdegs[i] < w->Ictx->subdegs[j])
            {
                j = i;
            }
        }
        w->Ictx->subdegs[j]++;
    }

got_ksub:

    fmpz_one(subprod);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        if ((slong)(w->Ictx->subdegs[i]) <= 0)
        {
            /* ksub has overflown - absolute falure */
            success = 0;
            goto cleanup;
        }
        fmpz_mul_ui(subprod, subprod, w->Ictx->subdegs[i]);
    }

    /* see if the ksub killed either lc(A) or lc(B) */
    fmpz_mpoly_ksub_content(cAksub, A->coeffs + 0, w->Ictx->subdegs, ctx);
    fmpz_mpoly_ksub_content(cBksub, B->coeffs + 0, w->Ictx->subdegs, ctx);
    if (fmpz_is_zero(cAksub) || fmpz_is_zero(cBksub))
        goto pick_ksub;

pick_bma_prime:

    if (fmpz_cmp_ui(p, ABtotal_length) < 0)
        fmpz_set_ui(p, ABtotal_length);

    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);

    success = fmpz_next_smooth_prime(p, p);
    if (!success)
    {
        /* ran out of smooth primes - absolute falure */
        success = 0;
        goto cleanup;
    }

    /* make sure reduction does not kill either leading coeff after ksub */
    if (fmpz_divisible(cAksub, p) || fmpz_divisible(cBksub, p))
    {
        goto pick_bma_prime;
    }

    /* make sure p does not divide any coefficient of Gamma */
    for (i = 0; i < Gamma->length; i++)
    {
        if (fmpz_divisible(Gamma->coeffs + i, p))
        {
            goto pick_bma_prime;
        }
    }

    /* try to get first image mod p */
    switch (fmpz_abs_fits_ui(p)
                 ? _bma_loop_sp(fmpz_get_ui(p), w, eval_sp_args, handles)
                 : _bma_loop_mp(p, w, eval_mp_args, handles))
    {
        default:
            FLINT_ASSERT(0);
        case ksub_probably_unlucky:
            goto pick_ksub;
        case insufficient_eval_points:
        case prime_probably_unlucky:
            goto pick_bma_prime;
        case gcd_is_one:
            fmpz_mpolyu_one(G, ctx);
            fmpz_mpolyu_swap(Abar, A, ctx);
            fmpz_mpolyu_swap(Bbar, B, ctx);
            success = 1;
            goto cleanup;
        case bma_loop_good:
            break;
    }

    /* Hmodulus was supposed to be set by bma_loop */
    FLINT_ASSERT(0 <= w->which_check && w->which_check <= 2);
    FLINT_ASSERT(fmpz_equal(w->Hmodulus, p));

    FLINT_ASSERT(w->H->length > 0);
    zip_evals = w->H->coeffs[0].length;
    for (i = 1; i < w->H->length; i++)
        zip_evals = FLINT_MAX(zip_evals, w->H->coeffs[i].length);
    zip_evals += 1;
    /* batches of size num_threads so overestimate fit size */
    nmod_zip_mpolyu_fit_poly(w->Z, w->H, zip_evals + w->num_threads);

    p_sp = UWORD(1) << (FLINT_BITS - 2);

pick_zip_prime:

    if (p_sp >= UWORD_MAX_PRIME)
    {
        /* ran out of machine primes - absolute failure */
        success = 0;
        goto cleanup;
    }
    p_sp = n_nextprime(p_sp, 1);

    if (0 == fmpz_fdiv_ui(w->Hmodulus, p_sp))
        goto pick_zip_prime;

    nmod_mpoly_ctx_set_modulus(w->ctx_sp, p_sp);
    /* unfortunate nmod_poly's store their own ctx :( */
    _base_args_set_mod_sp(w, eval_sp_args);

    FLINT_ASSERT(p_sp > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        w->alphas_sp[i] = n_urandint(w->randstate, p_sp - 3) + 2;

    /* TODO thread this? */
    nmod_zip_mpolyu_set_skel(w->Z, w->ctx_sp, w->H, w->alphas_sp, ctx);

    /* set skeletons for evaluation */
    _set_skels_sp(w, eval_sp_args, handles);

next_zip_image:

    j = zip_evals > w->Z->pointcount
            ? zip_evals - w->Z->pointcount
            : w->num_threads;
    j += w->num_threads - 1;
    j /= w->num_threads;
    j *= w->num_threads;
    _base_set_num_images_sp(w, j);

    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                            _worker_eval_sp, &eval_sp_args[i]);
    _worker_eval_sp(&eval_sp_args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    for (i = 0; i < w->num_images_sp; i++)
    {
        if (!w->evals_sp[i].success ||
            w->GdegboundXY < w->evals_sp[i].GevaldegXY)
        {
            /* start over if image failed or was unlucky */
            goto pick_zip_prime;
        }

        if (w->GdegboundXY > w->evals_sp[i].GevaldegXY)
        {
            /* we have a new degree bound on deg_XY(G) */
            w->GdegboundXY = w->evals_sp[i].GevaldegXY;
            if (w->GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                fmpz_mpolyu_swap(Abar, A, ctx);
                fmpz_mpolyu_swap(Bbar, B, ctx);
                success = 1;
                goto cleanup;
            }
            goto pick_bma_prime;
        }
    }

    FLINT_ASSERT(0 <= w->which_check && w->which_check <= 2);
    for (i = 0; i < w->num_images_sp; i++)
    {
        success = nmod_zip_mpolyuu_add_point(w->Z,
                        w->which_check == 1 ? w->evals_sp[i].nAbareval_sp :
                        w->which_check == 2 ? w->evals_sp[i].nBbareval_sp :
                                              w->evals_sp[i].nGeval_sp);
        if (!success)
            goto pick_bma_prime;
    }

    if (w->Z->pointcount < zip_evals)
        goto next_zip_image;

    /* find coeffs */
    w->zip_find_coeffs_no_match = 0;
    w->zip_find_coeffs_non_invertible = 0;
    w->index = 0;
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                    _worker_find_zip_coeffs, &eval_sp_args[i]);
    _worker_find_zip_coeffs(&eval_sp_args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (w->zip_find_coeffs_no_match)
        goto pick_bma_prime;
    else if (w->zip_find_coeffs_non_invertible)
        goto pick_zip_prime;

    FLINT_ASSERT(w->bits == w->H->bits);

    /* crt */
    w->Hchanged = 0;
    w->index = 0;
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                                     _worker_crt_zip_coeffs, &eval_sp_args[i]);
    _worker_crt_zip_coeffs(&eval_sp_args[0]);
    for (i = 1; i < w->num_threads; i++)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    fmpz_mul_ui(w->Hmodulus, w->Hmodulus, w->ctx_sp->mod.n);

    if (w->Hchanged)
    {
        /* TODO if the coefficients of H are getting to large? */
        goto pick_zip_prime;
    }

    success = fmpz_mpolyu_content_mpoly_threaded_pool(Hcontent, w->H, ctx,
                                                         handles, num_handles);
    FLINT_ASSERT(Hcontent->bits == w->bits);
    if (!success)
    {
        success = 0;
        goto cleanup;
    }

    if (w->which_check == 1)
    {
        fmpz_mpolyu_divexact_mpoly(Abar, w->H, 1, Hcontent, ctx);
        if (num_handles > 0)
        {
            success = fmpz_mpolyuu_divides_threaded_pool(
                                     G, A, Abar, 2, ctx, handles, num_handles);
            success = success && fmpz_mpolyuu_divides_threaded_pool(
                                     Bbar, B, G, 2, ctx, handles, num_handles);
        }
        else
        {
            success = fmpz_mpolyuu_divides(G, A, Abar, 2, ctx) &&
                      fmpz_mpolyuu_divides(Bbar, B, G, 2, ctx);
        }
    }
    else if (w->which_check == 2)
    {
        fmpz_mpolyu_divexact_mpoly(Bbar, w->H, 1, Hcontent, ctx);
        if (num_handles > 0)
        {
            success = fmpz_mpolyuu_divides_threaded_pool(
                                     G, B, Bbar, 2, ctx, handles, num_handles);
            success = success && fmpz_mpolyuu_divides_threaded_pool(
                                     Abar, A, G, 2, ctx, handles, num_handles);
        }
        else
        {
            success = fmpz_mpolyuu_divides(G, B, Bbar, 2, ctx) &&
                      fmpz_mpolyuu_divides(Abar, A, G, 2, ctx);
        }
    }
    else
    {
        FLINT_ASSERT(w->which_check == 0);
        fmpz_mpolyu_divexact_mpoly(G, w->H, 1, Hcontent, ctx);
        if (num_handles > 0)
        {
            /*
                Set n = num_handles. For some integer 0 <= m < n,
                A/G is processed by the m + 1 threads
                    main, (handles[0]  , ..., handles[m-1])
                B/G is processed by the n - m threads
                    handles[m], (handles[m+1], ..., handles[n-1])
            */
            slong m = mpoly_divide_threads(num_handles, A->length, B->length);
            _divide_arg_t divide_arg;

            divide_arg->ctx = ctx;
            divide_arg->den = G;

            if (m <= 0)
            {
                /* process A with one thread */
                divide_arg->quo = w->Bbar;
                divide_arg->num = B;
                divide_arg->handles = handles + 1;
                divide_arg->num_handles = num_handles - 1;
                thread_pool_wake(global_thread_pool, handles[0], 0,
                                                   _divide_worker, divide_arg);
                success = fmpz_mpolyuu_divides(w->Abar, A, G, 2, ctx);
                thread_pool_wait(global_thread_pool, handles[0]);
            }
            else if (m >= num_handles - 1)
            {
                /* process B with one thread */
                divide_arg->quo = w->Abar;
                divide_arg->num = A;
                divide_arg->handles = handles + 1;
                divide_arg->num_handles = num_handles - 1;
                thread_pool_wake(global_thread_pool, handles[0], 0,
                                                   _divide_worker, divide_arg);
                success = fmpz_mpolyuu_divides(w->Bbar, B, G, 2, ctx);
                thread_pool_wait(global_thread_pool, handles[0]);
            }
            else
            {
                divide_arg->quo = w->Bbar;
                divide_arg->num = B;
                divide_arg->handles = handles + (m + 1);
                divide_arg->num_handles = num_handles - (m + 1);
                thread_pool_wake(global_thread_pool, handles[m], 0,
                                                   _divide_worker, divide_arg);
                success = fmpz_mpolyuu_divides_threaded_pool(Abar, A, G, 2,
                                                          ctx, handles + 0, m);
                thread_pool_wait(global_thread_pool, handles[m]);
            }

            success = success && divide_arg->success;
        }
        else
        {
            success = fmpz_mpolyuu_divides(Abar, A, G, 2, ctx) &&
                      fmpz_mpolyuu_divides(Bbar, B, G, 2, ctx);
        }
    }

    if (!success)
        goto pick_zip_prime;

    success = 1;

cleanup:

    fmpz_mpoly_clear(Hcontent, ctx);
    fmpz_mpolyu_clear(w->H, ctx);

    flint_free(w->Adegs);
    flint_free(w->Bdegs);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&w->mutex);
#endif

    /* the zippler */
    nmod_zip_mpolyu_clear(w->Z);

    /* machine precision workspace */

    n_poly_clear(w->coeff_evals_sp);

    flint_free(w->alphas_sp);

    for (i = 0; i < w->num_threads; i++)
    {
        n_bpoly_clear(eval_sp_args[i].Aeval_sp);
        n_bpoly_clear(eval_sp_args[i].Beval_sp);
        n_bpoly_clear(eval_sp_args[i].Geval_sp);
        n_bpoly_clear(eval_sp_args[i].Abareval_sp);
        n_bpoly_clear(eval_sp_args[i].Bbareval_sp);
        n_bpoly_clear(eval_sp_args[i].Acur_sp);
        n_bpoly_clear(eval_sp_args[i].Bcur_sp);
        n_poly_clear(eval_sp_args[i].Gammacur_sp);
        n_poly_stack_clear(eval_sp_args[i].Sp_sp->poly_stack);
        n_bpoly_stack_clear(eval_sp_args[i].Sp_sp->bpoly_stack);
    }
    flint_free(eval_sp_args);

    n_bpoly_clear(w->Aone_sp);
    n_polyun_clear(w->Ainc_sp);
    n_bpoly_clear(w->Bone_sp);
    n_polyun_clear(w->Binc_sp);
    n_poly_clear(w->Gammaone_sp);
    n_poly_clear(w->Gammainc_sp);

    for (i = 0; i < w->evals_sp_alloc; i++)
    {
        n_bpoly_clear(w->evals_sp[i].nGeval_sp);
        n_bpoly_clear(w->evals_sp[i].nAbareval_sp);
        n_bpoly_clear(w->evals_sp[i].nBbareval_sp);
    }
    if (w->evals_sp_alloc > 0)
        flint_free(w->evals_sp);

    nmod_bma_mpoly_clear(w->GLambda_sp);
    nmod_bma_mpoly_clear(w->AbarLambda_sp);
    nmod_bma_mpoly_clear(w->BbarLambda_sp);
    nmod_mpoly_ctx_clear(w->ctx_sp);

    /* multiprecision workspace */

    fmpz_clear(w->alphashift_mp);
    fmpz_mpolyc_clear(w->coeff_evals_mp);

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(w->alphas_mp + i);
    flint_free(w->alphas_mp);

    for (i = 0; i < w->num_threads; i++)
    {
        fmpz_mod_mpolyn_clear(eval_mp_args[i].Aeval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_clear(eval_mp_args[i].Beval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_clear(eval_mp_args[i].Geval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_clear(eval_mp_args[i].Abareval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_clear(eval_mp_args[i].Bbareval_mp, w->ctx_mp);
        fmpz_mpolycu_clear(eval_mp_args[i].Acur_mp);
        fmpz_mpolycu_clear(eval_mp_args[i].Bcur_mp);
        fmpz_mpolyc_clear(eval_mp_args[i].Gammacur_mp);
    }
    flint_free(eval_mp_args);

    fmpz_mpolyc_clear(w->Gammaone_mp);
    fmpz_mpolycu_clear(w->Aone_mp);
    fmpz_mpolycu_clear(w->Bone_mp);
    fmpz_mpolyc_clear(w->Gammainc_mp);
    fmpz_mpolycu_clear(w->Ainc_mp);
    fmpz_mpolycu_clear(w->Binc_mp);
    fmpz_mpolyc_clear(w->Gammared_mp);
    fmpz_mpolycu_clear(w->Ared_mp);
    fmpz_mpolycu_clear(w->Bred_mp);

    for (i = 0; i < w->evals_mp_alloc; i++)
    {
        fmpz_mod_mpolyn_clear(w->evals_mp[i].nGeval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_clear(w->evals_mp[i].nAbareval_mp, w->ctx_mp);
        fmpz_mod_mpolyn_clear(w->evals_mp[i].nBbareval_mp, w->ctx_mp);
    }
    if (w->evals_mp_alloc > 0)
        flint_free(w->evals_mp);

    fmpz_mod_bma_mpoly_clear(w->GLambda_mp, w->ctx_mp->ffinfo);
    fmpz_mod_bma_mpoly_clear(w->AbarLambda_mp, w->ctx_mp->ffinfo);
    fmpz_mod_bma_mpoly_clear(w->BbarLambda_mp, w->ctx_mp->ffinfo);
    fmpz_mod_mpoly_ctx_clear(w->ctx_mp);

    /* misc */

    mpoly_bma_interpolate_ctx_clear(w->Ictx);

    fmpz_clear(w->Hmodulus);
    fmpz_clear(cBksub);
    fmpz_clear(cAksub);
    fmpz_clear(subprod);
    fmpz_clear(p);
    flint_randclear(w->randstate);

    FLINT_ASSERT(G->bits == w->bits);
    FLINT_ASSERT(Abar->bits == w->bits);
    FLINT_ASSERT(Bbar->bits == w->bits);

    return success;
}


typedef struct
{
    const fmpz_mpoly_struct * P;
    fmpz_mpoly_struct * Pcontent;
    fmpz_mpolyu_struct * Puu;
    const slong * perm;
    const ulong * shift, * stride, * maxexps;
    const fmpz_mpoly_ctx_struct * ctx;
    const fmpz_mpoly_ctx_struct * uctx;
    const thread_pool_handle * handles;
    slong num_handles;
    int success;
}
_convertuu_arg_struct;

typedef _convertuu_arg_struct _convertuu_arg_t[1];

static void _worker_convertuu(void * varg)
{
    _convertuu_arg_struct * arg = (_convertuu_arg_struct *) varg;

    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(arg->Puu, arg->uctx, arg->P, arg->ctx,
                             arg->perm, arg->shift, arg->stride, arg->maxexps,
                                               arg->handles, arg->num_handles);

    arg->success = fmpz_mpolyu_content_mpoly_threaded_pool(arg->Pcontent,
                          arg->Puu, arg->uctx, arg->handles, arg->num_handles);
    if (arg->success)
    {
        fmpz_mpolyu_divexact_mpoly_inplace(arg->Puu, arg->Pcontent, arg->uctx);
    }
}

int fmpz_mpoly_gcd_berlekamp_massey_threaded(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, num_handles;
    thread_pool_handle * handles;
    flint_bitcnt_t wbits;
    int success = 0;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    fmpz_mpoly_t Ac, Bc, Gc, Gamma;
    slong * Adegs, * Bdegs, * perm;
    ulong * shift, * stride;
    ulong max_main_degree, max_minor_degree;
    slong thread_limit = FLINT_MIN(A->length, B->length)/16;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars < 3)
    {
        return fmpz_mpoly_gcd_zippel(G, A, B, ctx);
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->nvars >= 3);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    /* get workers */
    num_handles = flint_request_threads(&handles, thread_limit);

    /* collect degree info */
    Adegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) flint_malloc((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));

    mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

    max_main_degree = 0;
    max_minor_degree = 0;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
        FLINT_ASSERT(Adegs[i] >= 0);
        FLINT_ASSERT(Bdegs[i] >= 0);
        if (i < 2)
        {
            max_main_degree = FLINT_MAX(max_main_degree, Adegs[i]);
            max_main_degree = FLINT_MAX(max_main_degree, Bdegs[i]);
        }
        else
        {
            max_minor_degree = FLINT_MAX(max_minor_degree, Adegs[i]);
            max_minor_degree = FLINT_MAX(max_minor_degree, Bdegs[i]);
        }
    }

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 2, ORD_LEX);

    /* wbits is bits for intermediates in ZZ[x_0,x_1][x_2,...,x_(n-1)] */
    wbits = 1 + FLINT_BIT_COUNT(max_minor_degree);
    wbits = FLINT_MAX(MPOLY_MIN_BITS, wbits);
    wbits = mpoly_fix_bits(wbits, uctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mpolyu_init(Auu, wbits, uctx);
    fmpz_mpolyu_init(Buu, wbits, uctx);
    fmpz_mpolyu_init(Guu, wbits, uctx);
    fmpz_mpolyu_init(Abaruu, wbits, uctx);
    fmpz_mpolyu_init(Bbaruu, wbits, uctx);
    fmpz_mpoly_init3(Ac, 0, wbits, uctx);
    fmpz_mpoly_init3(Bc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gamma, 0, wbits, uctx);

    /* two main variables must be packed into bits = FLINT_BITS/2 */
    if (FLINT_BIT_COUNT(max_main_degree) >= FLINT_BITS/2)
    {
        success = 0;
        goto cleanup;
    }

    if (num_handles > 0)
    {
        slong s = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertuu_arg_t arg;

        FLINT_ASSERT(s >= 0);
        FLINT_ASSERT(s < num_handles);

        arg->ctx = ctx;
        arg->uctx = uctx;
        arg->P = B;
        arg->Puu = Buu;
        arg->Pcontent = Bc;
        arg->perm = perm;
        arg->shift = shift;
        arg->stride = stride;
        arg->maxexps = (const ulong *) Bdegs;
        arg->handles = handles + (s + 1);
        arg->num_handles = num_handles - (s + 1);

        thread_pool_wake(global_thread_pool, handles[s], 0, _worker_convertuu, arg);

        fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
                   perm, shift, stride, (const ulong *) Adegs, handles + 0, s);
        success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Auu, uctx,
                                                               handles + 0, s);
        if (success)
        {
            fmpz_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
        }

        thread_pool_wait(global_thread_pool, handles[s]);

        success = success && arg->success;
        if (!success)
            goto cleanup;
    }
    else
    {
        fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
                          perm, shift, stride, (const ulong *) Adegs, NULL, 0);
        fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Buu, uctx, B, ctx,
                          perm, shift, stride, (const ulong *) Bdegs, NULL, 0);

        /* remove content from A and B */
        success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Auu, uctx, NULL, 0);
        success = success &&
                  fmpz_mpolyu_content_mpoly_threaded_pool(Bc, Buu, uctx, NULL, 0);
        if (!success)
            goto cleanup;

        fmpz_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
        fmpz_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);
    }

    success = _fmpz_mpoly_gcd_threaded_pool(Gamma, wbits, Auu->coeffs + 0,
                                  Buu->coeffs + 0, uctx, handles, num_handles);
    if (!success)
        goto cleanup;

    success = fmpz_mpolyuu_gcd_berlekamp_massey_threaded_pool(Guu, Abaruu, Bbaruu,
                                  Auu, Buu, Gamma, uctx, handles, num_handles);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx,
                                                         handles, num_handles);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    fmpz_mpoly_from_mpolyuu_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                               Guu, uctx, perm, shift, stride);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

    success = 1;

cleanup:

    flint_give_back_threads(handles, num_handles);

    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    fmpz_mpolyu_clear(Auu, uctx);
    fmpz_mpolyu_clear(Buu, uctx);
    fmpz_mpolyu_clear(Guu, uctx);
    fmpz_mpolyu_clear(Abaruu, uctx);
    fmpz_mpolyu_clear(Bbaruu, uctx);
    fmpz_mpoly_clear(Ac, uctx);
    fmpz_mpoly_clear(Bc, uctx);
    fmpz_mpoly_clear(Gc, uctx);
    fmpz_mpoly_clear(Gamma, uctx);

    fmpz_mpoly_ctx_clear(uctx);

    return success;
}

