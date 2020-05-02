/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "thread_pool.h"

typedef struct
{
    volatile int gcd_is_one;
    nmod_poly_struct * gamma;
    const nmod_mpoly_ctx_struct * ctx;
    nmod_mpolyn_struct * A, * B;
    ulong num_threads;
    slong var;
    slong bound;
    const mpoly_gcd_info_struct * I;
}
_splitbase_struct;

typedef _splitbase_struct _splitbase_t[1];

typedef struct
{
    slong idx;
    _splitbase_struct * base;
    nmod_mpolyn_t G, Abar, Bbar;
    nmod_poly_t modulus;
    mp_limb_t alpha;
    slong required_images;
}
_splitworker_arg_struct;

static void _splitworker_bivar(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    nmod_poly_t modulus2, alphapow, r;
    nmod_poly_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_poly_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyn_t T;
    mp_limb_t gammaevalp, alpha, temp;
    mp_limb_t gammaevalm;
    int gstab, astab, bstab, use_stab;
    slong ldeg;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(base->A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, base->A->bits, ctx->minfo);

    FLINT_ASSERT(base->var == 1);

    nmod_poly_init(r, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), base->bound + 1));

    nmod_poly_init(Aevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Aevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalm, ctx->ffinfo->mod.n);
    nmod_mpolyn_init(T, base->A->bits, ctx);

    alpha = arg->alpha;

    use_stab = 1;
    gstab = bstab = astab = 0;

    nmod_poly_one(arg->modulus);
    while (nmod_poly_degree(arg->modulus) < arg->required_images)
    {
        /* get evaluation point */
        if (alpha <= base->num_threads)
        {
            break;
        }
        alpha -= base->num_threads;

        FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
        FLINT_ASSERT(alphapow->alloc >= 2);
        alphapow->length = 2;
        alphapow->coeffs[0] = 1;
        alphapow->coeffs[1] = alpha;

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm,
                                           base->gamma, alphapow, ctx->ffinfo);
        if (gammaevalp == 0 || gammaevalm == 0)
        {
            continue;
        }

        /* evaluation point should kill neither A nor B */
        nmod_mpolyn_interp_reduce_2sm_poly(Aevalp, Aevalm,
                                                       base->A, alphapow, ctx);
        nmod_mpolyn_interp_reduce_2sm_poly(Bevalp, Bevalm,
                                                       base->B, alphapow, ctx);
        FLINT_ASSERT(Aevalp->length > 0);
        FLINT_ASSERT(Aevalm->length > 0);
        FLINT_ASSERT(Bevalp->length > 0);
        FLINT_ASSERT(Bevalm->length > 0);

        if (use_stab && gstab)
        {
            int success;
            slong Gdeg;

            nmod_mpolyn_interp_reduce_2sm_poly(Gevalp, Gevalm,
                                                        arg->G, alphapow, ctx);
            Gdeg = ((arg->G->exps + N*0)[off]>>shift);
            success = 1;
            success = success && nmod_poly_degree(Gevalp) == Gdeg;
            success = success && nmod_poly_degree(Gevalm) == Gdeg;
            success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
            success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
            nmod_poly_divrem_basecase(Abarevalp, r, Aevalp, Gevalp);
            success = success && (r->length == 0);
            nmod_poly_divrem_basecase(Abarevalm, r, Aevalm, Gevalm);
            success = success && (r->length == 0);
            nmod_poly_divrem_basecase(Bbarevalp, r, Bevalp, Gevalp);
            success = success && (r->length == 0);
            nmod_poly_divrem_basecase(Bbarevalm, r, Bevalm, Gevalm);
            success = success && (r->length == 0);

            if (!success)
            {
                use_stab = 0;
                nmod_poly_one(arg->modulus);
                alpha = arg->alpha;
                continue;
            }

            nmod_poly_scalar_mul_nmod(Abarevalp, Abarevalp, gammaevalp);
            nmod_poly_scalar_mul_nmod(Abarevalm, Abarevalm, gammaevalm);
            nmod_poly_scalar_mul_nmod(Bbarevalp, Bbarevalp, gammaevalp);
            nmod_poly_scalar_mul_nmod(Bbarevalm, Bbarevalm, gammaevalm);
        }
        else
        {
            nmod_poly_gcd(Gevalp, Aevalp, Bevalp);
            nmod_poly_div(Abarevalp, Aevalp, Gevalp);
            nmod_poly_div(Bbarevalp, Bevalp, Gevalp);
            nmod_poly_gcd(Gevalm, Aevalm, Bevalm);
            nmod_poly_div(Abarevalm, Aevalm, Gevalm);
            nmod_poly_div(Bbarevalm, Bevalm, Gevalm);
        }

        FLINT_ASSERT(Gevalp->length > 0);
        FLINT_ASSERT(Abarevalp->length > 0);
        FLINT_ASSERT(Bbarevalp->length > 0);
        FLINT_ASSERT(Gevalm->length > 0);
        FLINT_ASSERT(Abarevalm->length > 0);
        FLINT_ASSERT(Bbarevalm->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
            break;
        }
        if (nmod_poly_degree(Gevalp) == 0 || nmod_poly_degree(Gevalm) == 0)
        {
            base->gcd_is_one = 1;
            break;
        }

        if (nmod_poly_degree(Gevalp) != nmod_poly_degree(Gevalm))
        {
            continue;
        }

        /* the Geval have matching degrees */
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            FLINT_ASSERT(arg->G->length > 0);
            if (nmod_poly_degree(Gevalp) > ((arg->G->exps + N*0)[off]>>shift))
            {
                continue;
            }
            else if (nmod_poly_degree(Gevalp) < ((arg->G->exps + N*0)[off]>>shift))
            {
                nmod_poly_one(arg->modulus);
            }
        }
        /* update interpolants */
        nmod_poly_scalar_mul_nmod(Gevalp, Gevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Gevalm, Gevalm, gammaevalm);
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(arg->modulus, alpha);
            FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(arg->modulus,
                                                  ctx->ffinfo->mod.n - alpha));
            temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
            temp = nmod_add(temp, temp, ctx->ffinfo->mod);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(arg->modulus, arg->modulus, temp);
            if (!gstab)
            {
                gstab = !nmod_mpolyn_interp_crt_2sm_poly(&ldeg, arg->G, T,
                                 Gevalp, Gevalm, arg->modulus, alphapow, ctx);
            }
            nmod_mpolyn_interp_crt_2sm_poly(&ldeg, arg->Abar, T,
                            Abarevalp, Abarevalm, arg->modulus, alphapow, ctx);
            nmod_mpolyn_interp_crt_2sm_poly(&ldeg, arg->Bbar, T,
                            Bbarevalp, Bbarevalm, arg->modulus, alphapow, ctx);
        }
        else
        {
            nmod_mpolyn_interp_lift_2sm_poly(&ldeg, arg->G,
                                                   Gevalp, Gevalm, alpha, ctx);
            nmod_mpolyn_interp_lift_2sm_poly(&ldeg, arg->Abar,
                                             Abarevalp, Abarevalm, alpha, ctx);
            nmod_mpolyn_interp_lift_2sm_poly(&ldeg, arg->Bbar,
                                             Bbarevalp, Bbarevalm, alpha, ctx);
            gstab = astab = bstab = 0;
        }
        temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus2, arg->modulus, temp);
        nmod_poly_shift_left(arg->modulus, arg->modulus, 2);
        nmod_poly_sub(arg->modulus, arg->modulus, modulus2);
    }

    nmod_poly_clear(r);
    nmod_poly_clear(modulus2);
    nmod_poly_clear(alphapow);

    nmod_poly_clear(Aevalp);
    nmod_poly_clear(Bevalp);
    nmod_poly_clear(Gevalp);
    nmod_poly_clear(Abarevalp);
    nmod_poly_clear(Bbarevalp);
    nmod_poly_clear(Aevalm);
    nmod_poly_clear(Bevalm);
    nmod_poly_clear(Gevalm);
    nmod_poly_clear(Abarevalm);
    nmod_poly_clear(Bbarevalm);
    nmod_mpolyn_clear(T, ctx);
}


static void _splitworker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    flint_bitcnt_t bits = base->A->bits;
    slong var = base->var;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong offset, shift;
    nmod_poly_t modulus2, alphapow;
    nmod_mpolyn_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_mpolyn_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyn_t T;
    mp_limb_t gammaevalp, alpha, temp;
    mp_limb_t gammaevalm;
    slong ldeg;
    int success;
    nmod_poly_stack_t Sp;

    nmod_poly_stack_init(Sp, bits, ctx);

    FLINT_ASSERT(var > 0);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, bits, ctx->minfo);

    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), base->bound + 1));

    nmod_mpolyn_init(Aevalp, bits, ctx);
    nmod_mpolyn_init(Bevalp, bits, ctx);
    nmod_mpolyn_init(Gevalp, bits, ctx);
    nmod_mpolyn_init(Abarevalp, bits, ctx);
    nmod_mpolyn_init(Bbarevalp, bits, ctx);
    nmod_mpolyn_init(Aevalm, bits, ctx);
    nmod_mpolyn_init(Bevalm, bits, ctx);
    nmod_mpolyn_init(Gevalm, bits, ctx);
    nmod_mpolyn_init(Abarevalm, bits, ctx);
    nmod_mpolyn_init(Bbarevalm, bits, ctx);
    nmod_mpolyn_init(T, bits, ctx);

    alpha = arg->alpha;

    nmod_poly_one(arg->modulus);
    while (nmod_poly_degree(arg->modulus) < arg->required_images)
    {
        /* get evaluation point */
        if (alpha <= base->num_threads)
        {
            break;
        }
        alpha -= base->num_threads;

        FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
        FLINT_ASSERT(alphapow->alloc >= 2);
        alphapow->length = 2;
        alphapow->coeffs[0] = 1;
        alphapow->coeffs[1] = alpha;

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm,
                                           base->gamma, alphapow, ctx->ffinfo);
        if (gammaevalp == 0 || gammaevalm == 0)
        {
            continue;
        }

        /* evaluation should kill neither A nor B */
        nmod_mpolyn_interp_reduce_2sm_mpolyn(Aevalp, Aevalm,
                                                  base->A, var, alphapow, ctx);
        nmod_mpolyn_interp_reduce_2sm_mpolyn(Bevalp, Bevalm,
                                                  base->B, var, alphapow, ctx);
        FLINT_ASSERT(Aevalp->length > 0);
        FLINT_ASSERT(Aevalm->length > 0);
        FLINT_ASSERT(Bevalp->length > 0);
        FLINT_ASSERT(Bevalm->length > 0);

        success = nmod_mpolyn_gcd_brown_smprime(
                    Gevalp, Abarevalp, Bbarevalp, Aevalp, Bevalp, var - 1,
                                                             ctx, base->I, Sp);
        success = success && nmod_mpolyn_gcd_brown_smprime(
                    Gevalm, Abarevalm, Bbarevalm, Aevalm, Bevalm, var - 1,
                                                             ctx, base->I, Sp);
        if (success == 0)
        {
            continue;
        }

        FLINT_ASSERT(Gevalp->length > 0);
        FLINT_ASSERT(Abarevalp->length > 0);
        FLINT_ASSERT(Bbarevalp->length > 0);
        FLINT_ASSERT(Gevalm->length > 0);
        FLINT_ASSERT(Abarevalm->length > 0);
        FLINT_ASSERT(Bbarevalm->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
            break;
        }

        if (   nmod_mpolyn_is_nonzero_nmod(Gevalp, ctx)
            || nmod_mpolyn_is_nonzero_nmod(Gevalm, ctx))
        {
            base->gcd_is_one = 1;
            break;
        }

        if (   nmod_poly_degree(Gevalp->coeffs + 0)
            != nmod_poly_degree(Gevalm->coeffs + 0))
        {
            continue;
        }
        if (!mpoly_monomial_equal(Gevalp->exps + N*0, Gevalm->exps + N*0, N))
        {
            continue;
        }

        /* the Geval have matching degrees */
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            int cmp;
            slong k;

            FLINT_ASSERT(arg->G->length > 0);

            k = nmod_poly_degree(Gevalp->coeffs + 0);
            cmp = mpoly_monomial_cmp_nomask_extra(arg->G->exps + N*0,
                                    Gevalp->exps + N*0, N, offset, k << shift);
            if (cmp < 0)
            {
                continue;
            }
            else if (cmp > 0)
            {
                nmod_poly_one(arg->modulus);
            }
        }

        /* update interpolants */
        temp = nmod_mpolyn_leadcoeff(Gevalp, ctx);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        temp = nmod_mul(gammaevalp, temp, ctx->ffinfo->mod);
        nmod_mpolyn_scalar_mul_nmod(Gevalp, temp, ctx);
        temp = nmod_mpolyn_leadcoeff(Gevalm, ctx);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        temp = nmod_mul(gammaevalm, temp, ctx->ffinfo->mod);
        nmod_mpolyn_scalar_mul_nmod(Gevalm, temp, ctx);
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(arg->modulus, alpha);
            FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(arg->modulus, ctx->ffinfo->mod.n - alpha));
            temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
            temp = nmod_add(temp, temp, ctx->ffinfo->mod);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(arg->modulus, arg->modulus, temp);
            nmod_mpolyn_interp_crt_2sm_mpolyn(&ldeg, arg->G, T,
                             Gevalp, Gevalm, var, arg->modulus, alphapow, ctx);
            nmod_mpolyn_interp_crt_2sm_mpolyn(&ldeg, arg->Abar, T,
                       Abarevalp, Abarevalm, var, arg->modulus, alphapow, ctx);
            nmod_mpolyn_interp_crt_2sm_mpolyn(&ldeg, arg->Bbar, T,
                       Bbarevalp, Bbarevalm, var, arg->modulus, alphapow, ctx);
        }
        else
        {
            nmod_mpolyn_interp_lift_2sm_mpolyn(&ldeg, arg->G,
                                              Gevalp, Gevalm, var, alpha, ctx);
            nmod_mpolyn_interp_lift_2sm_mpolyn(&ldeg, arg->Abar,
                                        Abarevalp, Abarevalm, var, alpha, ctx);
            nmod_mpolyn_interp_lift_2sm_mpolyn(&ldeg, arg->Bbar,
                                        Bbarevalp, Bbarevalm, var, alpha, ctx);
        }
        temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus2, arg->modulus, temp);
        nmod_poly_shift_left(arg->modulus, arg->modulus, 2);
        nmod_poly_sub(arg->modulus, arg->modulus, modulus2);
    }

    nmod_poly_clear(modulus2);
    nmod_poly_clear(alphapow);

    nmod_mpolyn_clear(Aevalp, ctx);
    nmod_mpolyn_clear(Bevalp, ctx);
    nmod_mpolyn_clear(Gevalp, ctx);
    nmod_mpolyn_clear(Abarevalp, ctx);
    nmod_mpolyn_clear(Bbarevalp, ctx);
    nmod_mpolyn_clear(Aevalm, ctx);
    nmod_mpolyn_clear(Bevalm, ctx);
    nmod_mpolyn_clear(Gevalm, ctx);
    nmod_mpolyn_clear(Abarevalm, ctx);
    nmod_mpolyn_clear(Bbarevalm, ctx);
    nmod_mpolyn_clear(T, ctx);

    nmod_poly_stack_clear(Sp);
}


typedef struct
{
    slong hint_start, hint_stop;
    ulong * left_exp, * right_exp;
    nmod_mpolyn_t poly;
    slong lastdeg;
    slong thread_idx;
    slong final_idx;
    int GAB; /* 0 -> G,  1 -> A,  2 -> B */
}
_joinworker_arg_struct;


/*
    A = crt(B[0], ...., B[count-1]) wrt to P
    This function takes some preallocated temp space.
*/

static void _find_edge(
    slong * start,
    slong count,
    const ulong * exp_left,
    nmod_mpolyn_struct * const * B,
    slong N)
{
    slong k;

    for (k = 0; k < count; k++)
    {
        slong Blength = B[k]->length;
        const ulong * Bexps = B[k]->exps;

        if (start[k] < Blength
            && mpoly_monomial_gt_nomask(Bexps + N*start[k], exp_left, N))
        {
            /* go right */
            do {
                start[k]++;
            } while (start[k] < Blength
                && mpoly_monomial_gt_nomask(Bexps + N*start[k], exp_left, N));
        }
        else
        {
            /* go left */
            while (start[k] > 0
                && !mpoly_monomial_gt_nomask(Bexps + N*(start[k] - 1), exp_left, N))
            {
                start[k]--;
            }
        }
    }
}

static slong _nmod_mpolyn_crt(
    const nmod_poly_multi_crt_t P,
    _joinworker_arg_struct * S,
    nmod_mpolyn_struct * const * B,
    slong count,
    nmod_poly_struct * output,
    nmod_poly_struct ** input,
    const nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong N = mpoly_words_per_exp_sp(S->poly->bits, ctx->minfo);
    slong lastdegree;
    slong Ai;
    slong j, k;
    slong * start, * stop;
    nmod_poly_t zero;
    nmod_mpolyn_t A;
    const ulong * exp_left = S->left_exp;
    const ulong * exp_right = S->right_exp;
    TMP_INIT;

    *A = *S->poly;

    TMP_START;

    start = (slong *) TMP_ALLOC(2*count*sizeof(slong));
    stop = start + count;

    for (k = 0; k < count; k++)
    {
        start[k] = exp_left ? FLINT_MIN(S->hint_start, B[k]->length) : 0;
        stop[k] = exp_right ? FLINT_MIN(S->hint_stop, B[k]->length) : B[k]->length;
    }

    if (exp_left)
        _find_edge(start, count, exp_left, B, N);

    if (exp_right)
        _find_edge(stop, count, exp_right, B, N);

#if WANT_ASSERT
    for (k = 0; k < count; k++)
    {
        FLINT_ASSERT(0 <= start[k]);
        FLINT_ASSERT(0 <= stop[k]);
        FLINT_ASSERT(start[k] <= B[k]->length);
        FLINT_ASSERT(stop[k] <= B[k]->length);
        FLINT_ASSERT(start[k] <= stop[k]);

        /* check start */
        if (!exp_left)
        {
            FLINT_ASSERT(start[k] == 0);
        }
        else
        {
            FLINT_ASSERT(start[k] == B[k]->length ||
                    !mpoly_monomial_gt_nomask(B[k]->exps + N*start[k], exp_left, N));

            FLINT_ASSERT(start[k] == 0 ||
                    mpoly_monomial_gt_nomask(B[k]->exps + N*(start[k] - 1), exp_left, N));
        }

        /* check stop */
        if (!exp_right)
        {
            FLINT_ASSERT(stop[k] == B[k]->length);
        }
        else
        {
            FLINT_ASSERT(stop[k] == B[k]->length ||
                   !mpoly_monomial_gt_nomask(B[k]->exps + N*stop[k], exp_right, N));

            FLINT_ASSERT(stop[k] == 0 ||
                    mpoly_monomial_gt_nomask(B[k]->exps + N*(stop[k] - 1), exp_right, N));
        }
    }
#endif

    nmod_poly_init_mod(zero, ctx->ffinfo->mod);

    Ai = 0;
    lastdegree = -WORD(1);
    while (1)
    {
        nmod_mpolyn_fit_length(A, Ai + 1, ctx);

        k = 0;
        do
        {
            input[k] = zero;
            if (start[k] < stop[k])
            {
                goto found_max;
            }
        } while (++k < count);

        break; /* all B[k] have been scanned completely */

    found_max:

        input[k] = B[k]->coeffs + start[k];
        mpoly_monomial_set(A->exps + N*Ai, B[k]->exps + N*start[k], N);
        start[k]++;

        for (k++; k < count; k++)
        {
            input[k] = zero;
            if (start[k] >= stop[k])
            {
                continue;
            }

            cmp = mpoly_monomial_cmp_nomask(B[k]->exps + N*start[k], A->exps + N*Ai, N);
            if (cmp == 0)
            {
                input[k] = B[k]->coeffs + start[k];
                start[k]++;
            }
            else if (cmp > 0)
            {
                /* undo previous max's */
                for (j = 0; j < k; j++)
                {
                    start[j] -= (input[j] != zero);
                    input[j] = zero;
                }
                goto found_max;
            }
        }

        _nmod_poly_multi_crt_run_p(output, P,
                                     (const nmod_poly_struct * const *) input);
        nmod_poly_swap(A->coeffs + Ai, output + 0);
        lastdegree = FLINT_MAX(lastdegree, nmod_poly_degree(A->coeffs + Ai));
        Ai += !nmod_poly_is_zero(A->coeffs + Ai);
    }
    A->length = Ai;

    nmod_poly_clear(zero);

    TMP_END;

    *S->poly = *A;

    return lastdegree;
}


typedef struct
{
    volatile int idx;
    pthread_mutex_t mutex;
    const nmod_mpoly_ctx_struct * ctx;
    nmod_poly_multi_crt_t CRT;
    nmod_mpolyn_struct ** gptrs, ** abarptrs, ** bbarptrs;
    nmod_mpolyn_struct * G, * Abar, * Bbar;
    _joinworker_arg_struct * chunks;
    slong chunks_length;
    ulong num_threads;
}
_joinbase_struct;

typedef _joinbase_struct _joinbase_t[1];

typedef struct
{
    _joinbase_struct * base;
    slong thread_idx;
}
_njoinworker_arg_struct;

/* Join the images of G, Abar, and Bbar in each of the split threads */
static void _joinworker(void * varg)
{
    _njoinworker_arg_struct * arg = (_njoinworker_arg_struct *) varg;
    _joinbase_struct * base = arg->base;
    nmod_poly_struct ** input;
    nmod_poly_struct * output;
    slong i, ls = _nmod_poly_multi_crt_local_size(base->CRT);
    TMP_INIT;

    TMP_START;

    input = (nmod_poly_struct **) TMP_ALLOC(
                               base->num_threads * sizeof(nmod_poly_struct *));
    output = (nmod_poly_struct *) TMP_ALLOC(ls*sizeof(nmod_poly_struct));
    for (i = 0; i < ls; i++)
        nmod_poly_init_mod(output + i, base->ctx->ffinfo->mod);

    while (1)
    {
        /* get exponent of either G, Abar, or Bbar to start working on */
        pthread_mutex_lock(&base->mutex);
        i = base->idx;
        base->idx = i + 1;
        pthread_mutex_unlock(&base->mutex);

        if (i >= base->chunks_length)
        {
            goto cleanup;
        }

        base->chunks[i].thread_idx = arg->thread_idx;

        if (base->chunks[i].GAB == 0)
        {
            base->chunks[i].lastdeg = _nmod_mpolyn_crt(base->CRT,
                             base->chunks + i, base->gptrs, base->num_threads,
                                                     output, input, base->ctx);
        }
        else if (base->chunks[i].GAB == 1)
        {
            base->chunks[i].lastdeg = _nmod_mpolyn_crt(base->CRT,
                          base->chunks + i, base->abarptrs, base->num_threads,
                                                     output, input, base->ctx);
        }
        else
        {
            FLINT_ASSERT(base->chunks[i].GAB == 2);

            base->chunks[i].lastdeg = _nmod_mpolyn_crt(base->CRT,
                          base->chunks + i, base->bbarptrs, base->num_threads,
                                                     output, input, base->ctx);
        }
    }

cleanup:

    for (i = 0; i < ls; i++)
        nmod_poly_clear(output + i);

    TMP_END;

    return;
}


static void _finaljoinworker(void * varg)
{
    _njoinworker_arg_struct * arg = (_njoinworker_arg_struct *) varg;
    _joinbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    flint_bitcnt_t bits = base->G->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong i, j;
    slong source_len;
    ulong * source_exps;
    nmod_poly_struct * source_coeffs;
    slong Ti;
    ulong * Texps;
    nmod_poly_struct * Tcoeffs;

    for (i = base->chunks_length - 1; i >= 0; i--)
    {
        int type = base->chunks[i].GAB;

        FLINT_ASSERT(base->chunks[i].thread_idx >= 0);

        if (base->chunks[i].thread_idx != arg->thread_idx)
            continue;

        if (type == 0)
        {
            Texps = base->G->exps;
            Tcoeffs = base->G->coeffs;
        }
        else if (type == 1)
        {
            Texps = base->Abar->exps;
            Tcoeffs = base->Abar->coeffs;
        }
        else
        {
            FLINT_ASSERT(type == 2);
            Texps = base->Bbar->exps;
            Tcoeffs = base->Bbar->coeffs;
        }

        source_len = base->chunks[i].poly->length;
        source_exps = base->chunks[i].poly->exps + N*0;
        source_coeffs = base->chunks[i].poly->coeffs + 0;

        Ti = base->chunks[i].final_idx;
        mpoly_copy_monomials(Texps + N*Ti, source_exps, source_len, N);
        for (j = 0; j < source_len; j++)
            nmod_poly_swap(Tcoeffs + Ti + j, source_coeffs + j);
    }
}

/*
    Do same as nmod_mpolyun_gcd_brown_smprime but use the threads in
        handles[0], ..., handles[num_handles - 1]
    num_handles is allowed to be zero.
*/
int nmod_mpolyn_gcd_brown_smprime_threaded_pool(
    nmod_mpolyn_t G,
    nmod_mpolyn_t Abar,
    nmod_mpolyn_t Bbar,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int divisibility_test = 0; /* 1: by G, 2: by Abar, 3: by Bbar */
    slong i, j;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong offset, shift;
    ulong num_threads;
    int success;
    ulong bound, best_est;
    slong g_stab_est, abar_stab_est, bbar_stab_est, upper_limit;
    mp_limb_t alpha;
    slong deggamma, ldegA, ldegB;
    slong ldegGs_Abars_Bbars[3];
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_struct ** mptrs;
    nmod_mpolyn_struct ** gptrs, ** abarptrs, ** bbarptrs;
    _splitworker_arg_struct * splitargs;
    _splitbase_t splitbase;
    _njoinworker_arg_struct * joinargs;
    _joinbase_t joinbase;
    nmod_poly_t t1;
    nmod_mpolyn_t T1, T2;
#if WANT_ASSERT
    nmod_mpolyn_t Aorg, Borg;
#endif

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);

#if WANT_ASSERT
    nmod_mpolyn_init(Aorg, A->bits, ctx);
    nmod_mpolyn_init(Borg, B->bits, ctx);
    nmod_mpolyn_set(Aorg, A, ctx);
    nmod_mpolyn_set(Borg, B, ctx);
#endif

    mpoly_gen_offset_shift_sp(&offset, &shift, 0, G->bits, ctx->minfo);

    nmod_poly_init_mod(t1, ctx->ffinfo->mod);
    nmod_mpolyn_init(T1, bits, ctx);
    nmod_mpolyn_init(T2, bits, ctx);

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyn_content_last(cA, A, ctx);
    nmod_mpolyn_content_last(cB, B, ctx);
    nmod_mpolyn_divexact_last(A, cA, ctx);
    nmod_mpolyn_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyn_leadcoeff_poly(A, ctx),
                         nmod_mpolyn_leadcoeff_poly(B, ctx));

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    ldegB = nmod_mpolyn_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);
    best_est = bound;

    upper_limit = mpoly_gcd_info_get_brown_upper_limit(I, var, bound);

    if (I != NULL && I->Gdeflate_deg_bounds_are_nice)
    {
        slong k = I->brown_perm[var];

        FLINT_ASSERT(var < I->mvars);
        
        g_stab_est = 2 + deggamma + I->Gdeflate_deg_bound[k];
        abar_stab_est = 2 + deggamma + I->Adeflate_deg[k] - I->Gdeflate_deg_bound[k];
        bbar_stab_est = 2 + deggamma + I->Bdeflate_deg[k] - I->Gdeflate_deg_bound[k];

        if (g_stab_est < upper_limit)
        {
            best_est = g_stab_est;
            divisibility_test = 1;
        }

        if (abar_stab_est < upper_limit
             && (divisibility_test == 0 || abar_stab_est < best_est))
        {
            best_est = abar_stab_est;
            divisibility_test = 2;
        }

        if (bbar_stab_est < upper_limit
             && (divisibility_test == 0 || bbar_stab_est < best_est))
        {
            best_est = bbar_stab_est;
            divisibility_test = 3;
        }
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);
    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    num_threads = num_handles + 1;
    gptrs = (nmod_mpolyn_struct **) flint_malloc(
                                    num_threads*sizeof(nmod_mpolyn_struct *));
    abarptrs = (nmod_mpolyn_struct **) flint_malloc(
                                    num_threads*sizeof(nmod_mpolyn_struct *));
    bbarptrs = (nmod_mpolyn_struct **) flint_malloc(
                                    num_threads*sizeof(nmod_mpolyn_struct *));
    mptrs = (nmod_poly_struct **) flint_malloc(
                                       num_threads*sizeof(nmod_poly_struct *));
    splitargs = (_splitworker_arg_struct *) flint_malloc(
                                  num_threads*sizeof(_splitworker_arg_struct));
    for (i = 0; i < num_threads; i++)
    {
        nmod_mpolyn_init(splitargs[i].G, bits, ctx);
        nmod_mpolyn_init(splitargs[i].Abar, bits, ctx);
        nmod_mpolyn_init(splitargs[i].Bbar, bits, ctx);
        nmod_poly_init(splitargs[i].modulus, ctx->ffinfo->mod.n);
    }
    splitbase->num_threads = num_threads;
    splitbase->A = A;
    splitbase->B = B;
    splitbase->ctx = ctx;
    splitbase->gamma = gamma;
    splitbase->var = var;
    splitbase->I = I;

compute_split:

    splitbase->bound = best_est;

    if (alpha <= num_threads)
    {
        success = 0;
        goto cleanup_split;
    }
    alpha -= num_threads;

    splitbase->gcd_is_one = 0;
    for (i = 0; i < num_threads; i++)
    {
        slong ri = best_est / num_threads + (i < (best_est % num_threads));
        splitargs[i].idx = i;
        splitargs[i].base = splitbase;
        splitargs[i].alpha = alpha + i;
        splitargs[i].required_images = FLINT_MAX(ri, WORD(2));
    }

    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i], 0,
                  var == 1 ? _splitworker_bivar : _splitworker, &splitargs[i]);
    }
    (var == 1 ? _splitworker_bivar : _splitworker)(&splitargs[num_threads - 1]);
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    if (splitbase->gcd_is_one)
    {
        nmod_mpolyn_one(G, ctx);
        nmod_mpolyn_swap(Abar, A);
        nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;
    }

    for (i = 0; i < num_threads; i++)
    {
        gptrs[i] = splitargs[i].G;
        abarptrs[i] = splitargs[i].Abar;
        bbarptrs[i] = splitargs[i].Bbar;
        mptrs[i] = splitargs[i].modulus;
        if (nmod_poly_degree(splitargs[i].modulus) < splitargs[i].required_images)
        {
            /* not enough evaluation points - must fail */
            success = 0;
            goto cleanup_split;
        }
        FLINT_ASSERT(gptrs[i]->length > 0);
        FLINT_ASSERT(abarptrs[i]->length > 0);
        FLINT_ASSERT(bbarptrs[i]->length > 0);
    }

    /*
        Check for consistency in the leading monomial. All args have at least
        one image, so G, Abar, Bbar are defined and nonzero for each.
    */
    for (i = 1; i < num_threads; i++)
    {
        if (!mpoly_monomial_equal(gptrs[i]->exps + N*0, gptrs[0]->exps + N*0, N))
        {
            /* very unlucky - could try again or just fail */
            goto compute_split;
        }
    }

    nmod_poly_multi_crt_init(joinbase->CRT);
    success = nmod_poly_multi_crt_precompute_p(joinbase->CRT,
                        (const nmod_poly_struct * const *) mptrs, num_threads);
    FLINT_ASSERT(success);

    joinbase->num_threads = num_threads;
    joinbase->gptrs = gptrs;
    joinbase->abarptrs = abarptrs;
    joinbase->bbarptrs = bbarptrs;
    joinbase->G = G;
    joinbase->Abar = Abar;
    joinbase->Bbar = Bbar;
    joinbase->ctx = ctx;
    pthread_mutex_init(&joinbase->mutex, NULL);

    joinargs = (_njoinworker_arg_struct *) flint_malloc(
                                   num_threads*sizeof(_njoinworker_arg_struct));

    for (i = 0; i < num_threads; i++)
    {
        joinargs[i].base = joinbase;
        joinargs[i].thread_idx = i;
    }

    joinbase->chunks_length = 3*num_threads;
    joinbase->chunks = (_joinworker_arg_struct *) flint_malloc(
                   joinbase->chunks_length*sizeof(_joinworker_arg_struct));

    FLINT_ASSERT(joinbase->chunks_length >= 3);

    for (i = 0; i < 3; i++)
    {
        nmod_mpolyn_struct * poly = (i == 0) ? gptrs[0]
                                  : (i == 1) ? abarptrs[0]
                                  :            bbarptrs[0];

        for (j = 0; j < num_threads; j++)
        {
            _joinworker_arg_struct * d = joinbase->chunks + i*num_threads + j;
            nmod_mpolyn_init(d->poly, bits, ctx);
            d->GAB = i;
            d->lastdeg = -WORD(1);
            d->thread_idx = -WORD(1);
            d->final_idx = -WORD(1);

            d->hint_start = poly->length * (j + 0) / num_threads;
            d->hint_stop  = poly->length * (j + 1) / num_threads;

            d->left_exp = NULL;
            d->right_exp = NULL;

            FLINT_ASSERT(0 <= d->hint_start);
            FLINT_ASSERT(d->hint_start <= d->hint_stop);
            FLINT_ASSERT(d->hint_stop <= poly->length);
            FLINT_ASSERT(d->hint_start < poly->length);

            if (j > 0)
            {
                FLINT_ASSERT(d->hint_start < poly->length);
                d->left_exp = poly->exps + N*d->hint_start;
            }

            if (j + 1 < num_threads)
            {
                FLINT_ASSERT(d->hint_stop < poly->length);
                d->right_exp = poly->exps + N*d->hint_stop;
            }
        }
    }

    joinbase->idx = 0;
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool,
                                 handles[i], 0, _joinworker, joinargs + i);
    }
    _joinworker(joinargs + num_threads - 1);
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    /* final trivial join */
    ldegGs_Abars_Bbars[0] = -WORD(1);
    ldegGs_Abars_Bbars[1] = -WORD(1);
    ldegGs_Abars_Bbars[2] = -WORD(1);
    {
        slong idxs[3] = {0, 0, 0};
        for (i = 0; i < joinbase->chunks_length; i++)
        {
            int type = joinbase->chunks[i].GAB;
            FLINT_ASSERT(0 <= type && type < 3);

            joinbase->chunks[i].final_idx = idxs[type];
            idxs[type] += joinbase->chunks[i].poly->length;

            ldegGs_Abars_Bbars[type] = FLINT_MAX(ldegGs_Abars_Bbars[type],
                                                 joinbase->chunks[i].lastdeg);
        }

        nmod_mpolyn_fit_length(G, idxs[0], ctx);
        nmod_mpolyn_fit_length(Abar, idxs[1], ctx);
        nmod_mpolyn_fit_length(Bbar, idxs[2], ctx);
        G->length = idxs[0];
        Abar->length = idxs[1];
        Bbar->length = idxs[2];
    }

    joinbase->idx = 0;
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool,
                                handles[i], 0, _finaljoinworker, joinargs + i);
    }
    _finaljoinworker(joinargs + num_threads - 1);
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    FLINT_ASSERT(nmod_mpolyn_is_canonical(G, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(Abar, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(Bbar, ctx));

    pthread_mutex_destroy(&joinbase->mutex);

    /* free join data */
    nmod_poly_multi_crt_clear(joinbase->CRT);
    for (i = 0; i < joinbase->chunks_length; i++)
    {
        nmod_mpolyn_clear(joinbase->chunks[i].poly, ctx);
    }

    flint_free(joinbase->chunks);
    flint_free(joinargs);

    if (divisibility_test == 1)
    {
        nmod_mpolyn_content_last(t1, G, ctx);
        nmod_mpolyn_divexact_last(G, t1, ctx);
        success =            nmod_mpolyn_divides_threaded_pool(T1, A, G,
                                                    ctx, handles, num_handles);
        success = success && nmod_mpolyn_divides_threaded_pool(T2, B, G,
                                                    ctx, handles, num_handles);
        if (success)
        {
            ulong temp;
            nmod_mpolyn_swap(T1, Abar);
            nmod_mpolyn_swap(T2, Bbar);
successful_fix_lc:
            temp = nmod_mpolyn_leadcoeff(G, ctx);
            nmod_mpolyn_scalar_mul_nmod(Abar, temp, ctx);
            nmod_mpolyn_scalar_mul_nmod(Bbar, temp, ctx);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_mpolyn_scalar_mul_nmod(G, temp, ctx);
            goto successful_put_content;
        }
    }
    else if (divisibility_test == 2)
    {
        nmod_mpolyn_content_last(t1, Abar, ctx);
        nmod_mpolyn_divexact_last(Abar, t1, ctx);
        success =            nmod_mpolyn_divides_threaded_pool(T1, A, Abar,
                                                    ctx, handles, num_handles);
        success = success && nmod_mpolyn_divides_threaded_pool(T2, B, T1,
                                                    ctx, handles, num_handles);
        if (success)
        {
            nmod_mpolyn_swap(T1, G);
            nmod_mpolyn_swap(T2, Bbar);
            goto successful_fix_lc;
        }
    }
    else if (divisibility_test == 3)
    {
        nmod_mpolyn_content_last(t1, Bbar, ctx);
        nmod_mpolyn_divexact_last(Bbar, t1, ctx);
        success =            nmod_mpolyn_divides_threaded_pool(T1, B, Bbar,
                                                    ctx, handles, num_handles);
        success = success && nmod_mpolyn_divides_threaded_pool(T2, A, T1,
                                                    ctx, handles, num_handles);
        if (success)
        {
            nmod_mpolyn_swap(T1, G);
            nmod_mpolyn_swap(T2, Abar);
            goto successful_fix_lc;
        }        
    }
    else /* divisibility test == 0 */
    {
        if (   deggamma + ldegA == ldegGs_Abars_Bbars[0] + ldegGs_Abars_Bbars[1]
            && deggamma + ldegB == ldegGs_Abars_Bbars[0] + ldegGs_Abars_Bbars[2])
        {
            goto successful;
        }
    }

    /* divisibility test failed - try again */
    best_est = bound;
    divisibility_test = 0;
    goto compute_split;

successful:

    nmod_mpolyn_content_last(t1, G, ctx);
    nmod_mpolyn_divexact_last(G, t1, ctx);
    nmod_mpolyn_divexact_last(Abar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyn_divexact_last(Bbar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyn_mul_last(G, cG, ctx);
    nmod_mpolyn_mul_last(Abar, cAbar, ctx);
    nmod_mpolyn_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup_split:

    for (i = 0; i < num_threads; i++)
    {
        nmod_mpolyn_clear(splitargs[i].G, ctx);
        nmod_mpolyn_clear(splitargs[i].Abar, ctx);
        nmod_mpolyn_clear(splitargs[i].Bbar, ctx);
        nmod_poly_clear(splitargs[i].modulus);
    }

    flint_free(gptrs);
    flint_free(abarptrs);
    flint_free(bbarptrs);
    flint_free(mptrs);
    flint_free(splitargs);

cleanup:

#if WANT_ASSERT
    if (success)
    {
        success =            nmod_mpolyn_divides(T1, Aorg, G, ctx);
        success = success && nmod_mpolyn_divides(T2, Borg, G, ctx);
        FLINT_ASSERT(success);
        FLINT_ASSERT(nmod_mpolyn_equal(T1, Abar, ctx));
        FLINT_ASSERT(nmod_mpolyn_equal(T2, Bbar, ctx));

        success = 1;
    }
    nmod_mpolyn_clear(Aorg, ctx);
    nmod_mpolyn_clear(Borg, ctx);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);
    nmod_poly_clear(gamma);

    nmod_poly_clear(t1);
    nmod_mpolyn_clear(T1, ctx);
    nmod_mpolyn_clear(T2, ctx);

    return success;
}




typedef struct
{
    nmod_mpolyn_struct * Pn;
    const nmod_mpoly_ctx_struct * nctx;
    const nmod_mpoly_struct * P;
    const nmod_mpoly_ctx_struct * ctx;
    const slong * perm;
    const ulong * shift;
    const ulong * stride;
    const thread_pool_handle * handles;
    slong num_handles;
}
_convertn_arg_struct;

typedef _convertn_arg_struct _convertn_arg_t[1];

static void _worker_convertn(void * varg)
{
    _convertn_arg_struct * arg = (_convertn_arg_struct *) varg;

    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(arg->Pn, arg->nctx,
                    arg->P, arg->ctx, arg->perm, arg->shift, arg->stride,
                                               arg->handles, arg->num_handles);
}

int nmod_mpoly_gcd_brown_threaded(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    flint_bitcnt_t ABbits;
    nmod_mpoly_ctx_t nctx;
    nmod_mpolyn_t An, Bn, Gn, Abarn, Bbarn;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit = FLINT_MIN(A->length, B->length)/16;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
        {
            nmod_mpoly_zero(G, ctx);
        }
        else
        {
            nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    if (ctx->minfo->nvars == 1)
    {
        nmod_poly_t a, b, g;
        nmod_poly_init(a, ctx->ffinfo->mod.n);
        nmod_poly_init(b, ctx->ffinfo->mod.n);
        nmod_poly_init(g, ctx->ffinfo->mod.n);
        _nmod_mpoly_to_nmod_poly_deflate(a, A, 0, shift, stride, ctx);
        _nmod_mpoly_to_nmod_poly_deflate(b, B, 0, shift, stride, ctx);
        nmod_poly_gcd(g, a, b);
        _nmod_mpoly_from_nmod_poly_inflate(G, A->bits, g, 0, shift, stride, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        success = 1;
        goto cleanup1;
    }

    FLINT_ASSERT(ctx->minfo->nvars >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    ABbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(nctx, ctx->minfo->nvars, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyn_init(An, ABbits, nctx);
    nmod_mpolyn_init(Bn, ABbits, nctx);
    nmod_mpolyn_init(Gn, ABbits, nctx);
    nmod_mpolyn_init(Abarn, ABbits, nctx);
    nmod_mpolyn_init(Bbarn, ABbits, nctx);

    num_handles = flint_request_threads(&handles, thread_limit);

    /* convert inputs */
    if (num_handles > 0)
    {
        slong m = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertn_arg_t arg;

        FLINT_ASSERT(m >= 0);
        FLINT_ASSERT(m < num_handles);

        arg->Pn = Bn;
        arg->nctx = nctx;
        arg->P = B;
        arg->ctx = ctx;
        arg->perm = perm;
        arg->shift = shift;
        arg->stride = stride;
        arg->handles = handles + (m + 1);
        arg->num_handles = num_handles - (m + 1);

        thread_pool_wake(global_thread_pool, handles[m], 0, _worker_convertn, arg);

        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                                          perm, shift, stride, handles + 0, m);

        thread_pool_wait(global_thread_pool, handles[m]);
    }
    else
    {
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx,
                                         A, ctx, perm, shift, stride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(Bn, nctx,
                                         B, ctx, perm, shift, stride, NULL, 0);
    }

    /* calculate gcd */
    success = nmod_mpolyn_gcd_brown_smprime_threaded_pool(
                          Gn, Abarn, Bbarn, An, Bn, nctx->minfo->nvars - 1,
                                             nctx, NULL, handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    if (!success)
    {
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx,
                                         A, ctx, perm, shift, stride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(Bn, nctx,
                                         B, ctx, perm, shift, stride, NULL, 0);
        success = nmod_mpolyn_gcd_brown_lgprime(Gn, Abarn, Bbarn,
                                         An, Bn, nctx->minfo->nvars - 1, nctx);
    }

    if (success)
    {
        nmod_mpoly_from_mpolyn_perm_inflate(G, ABbits, ctx,
                                                Gn, nctx, perm, shift, stride);
        nmod_mpoly_make_monic(G, G, ctx);
    }

    nmod_mpolyn_clear(An, nctx);
    nmod_mpolyn_clear(Bn, nctx);
    nmod_mpolyn_clear(Gn, nctx);
    nmod_mpolyn_clear(Abarn, nctx);
    nmod_mpolyn_clear(Bbarn, nctx);
    nmod_mpoly_ctx_clear(nctx);

cleanup1:

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}

