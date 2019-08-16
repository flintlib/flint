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

/*
    A = crt(B[0], ...., B[count-1]) wrt to P
    This functions takes some preallocated temp space.
*/
static slong _nmod_mpolyn_crt(
    const nmod_poly_multi_crt_t P,
    nmod_mpolyn_t A,
    nmod_mpolyn_struct * const * B,
    slong count,
    nmod_poly_struct * output,      /* temp space */
    nmod_poly_struct ** input,      /* temp space */
    slong * start,                  /* temp space */
    const nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdegree;
    slong Ai;
    slong j, k;
    nmod_poly_t zero;

    nmod_poly_init(zero, ctx->ffinfo->mod.n);

    /* start[k] is the next available term in B[k] */
    for (k = 0; k < count; k++)
    {
        start[k] = 0;
    }

    Ai = 0;
    lastdegree = -WORD(1);
    while (1)
    {
        nmod_mpolyn_fit_length(A, Ai + 1, ctx);

        k = 0;
        do
        {
            input[k] = zero;
            if (start[k] < B[k]->length)
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
            if (start[k] >= B[k]->length)
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
    return lastdegree;
}

/*
    Append to A the result of crt'ing the coeff of X^exp wrt P of the B's
    This functions takes some preallocated temp space.
    A and the B's are in Fp[X][x_0, ..., x_(var-1)][x_var]
*/
static slong _nmod_mpolyun_crt_exp(
    const nmod_poly_multi_crt_t P,
    nmod_mpolyun_t A,
    ulong exp,
    nmod_mpolyun_struct * const * B,
    slong count,
    nmod_mpolyn_struct ** Bcoeffs,  /* temp space */
    nmod_poly_struct ** input,      /* temp space */
    slong * start,                  /* temp space */
    nmod_poly_struct * output,      /* temp space */
    const nmod_mpoly_ctx_t ctx)
{
    slong j, k;
    slong Ai;
    slong lastdegree;
    nmod_mpolyn_t zero;

    nmod_mpolyn_init(zero, A->bits, ctx);

    for (k = 0; k < count; k++)
    {
        Bcoeffs[k] = zero;
        for (j = 0; j < B[k]->length; j++)
        {
            if (B[k]->exps[j] == exp)
            {
                Bcoeffs[k] = B[k]->coeffs + j;
                break;
            }
        }
    }

    Ai = A->length;
    nmod_mpolyun_fit_length(A, Ai + 1, ctx);
    A->exps[Ai] = exp;
    lastdegree = _nmod_mpolyn_crt(P, A->coeffs + Ai, Bcoeffs, count,
                                                    output, input, start, ctx);
    A->length += (A->coeffs + Ai)->length != 0;

    nmod_mpolyn_clear(zero, ctx);
    return lastdegree;
}


typedef struct
{
    volatile int gcd_is_one;
    nmod_poly_struct * gamma;
    const nmod_mpoly_ctx_struct * ctx;
    nmod_mpolyun_struct * A, * B;
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
    nmod_mpolyun_t G, Abar, Bbar;
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
    nmod_mpolyun_t T;
    mp_limb_t gammaevalp, alpha, temp;
    mp_limb_t gammaevalm;
    int gstab, astab, bstab, use_stab;
    slong ldeg;

    FLINT_ASSERT(base->var == 0);

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
    nmod_mpolyun_init(T, base->A->bits, ctx);

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
        nmod_mpolyun_intp_reduce_2sm_poly(Aevalp, Aevalm,
                                                       base->A, alphapow, ctx);
        nmod_mpolyun_intp_reduce_2sm_poly(Bevalp, Bevalm,
                                                       base->B, alphapow, ctx);
        FLINT_ASSERT(Aevalp->length > 0);
        FLINT_ASSERT(Aevalm->length > 0);
        FLINT_ASSERT(Bevalp->length > 0);
        FLINT_ASSERT(Bevalm->length > 0);

        if (use_stab && gstab)
        {
            int success;
            slong Gdeg;

            nmod_mpolyun_intp_reduce_2sm_poly(Gevalp, Gevalm,
                                                        arg->G, alphapow, ctx);
            Gdeg = arg->G->exps[0];
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
            if (nmod_poly_degree(Gevalp) > arg->G->exps[0])
            {
                continue;
            }
            else if (nmod_poly_degree(Gevalp) < arg->G->exps[0])
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
                gstab = !nmod_mpolyun_intp_crt_2sm_poly(&ldeg, arg->G, T,
                                 Gevalp, Gevalm, arg->modulus, alphapow, ctx);
            }
            nmod_mpolyun_intp_crt_2sm_poly(&ldeg, arg->Abar, T,
                            Abarevalp, Abarevalm, arg->modulus, alphapow, ctx);
            nmod_mpolyun_intp_crt_2sm_poly(&ldeg, arg->Bbar, T,
                            Bbarevalp, Bbarevalm, arg->modulus, alphapow, ctx);
        }
        else
        {
            nmod_mpolyun_intp_lift_2sm_poly(&ldeg, arg->G,
                                                   Gevalp, Gevalm, alpha, ctx);
            nmod_mpolyun_intp_lift_2sm_poly(&ldeg, arg->Abar,
                                             Abarevalp, Abarevalm, alpha, ctx);
            nmod_mpolyun_intp_lift_2sm_poly(&ldeg, arg->Bbar,
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
    nmod_mpolyun_clear(T, ctx);
}


static void _splitworker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    flint_bitcnt_t bits = base->A->bits;
    slong var = base->var;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong offset, shift;
    nmod_poly_t modulus2, alphapow;
    nmod_mpolyun_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_mpolyun_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
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

    nmod_mpolyun_init(Aevalp, bits, ctx);
    nmod_mpolyun_init(Bevalp, bits, ctx);
    nmod_mpolyun_init(Gevalp, bits, ctx);
    nmod_mpolyun_init(Abarevalp, bits, ctx);
    nmod_mpolyun_init(Bbarevalp, bits, ctx);
    nmod_mpolyun_init(Aevalm, bits, ctx);
    nmod_mpolyun_init(Bevalm, bits, ctx);
    nmod_mpolyun_init(Gevalm, bits, ctx);
    nmod_mpolyun_init(Abarevalm, bits, ctx);
    nmod_mpolyun_init(Bbarevalm, bits, ctx);
    nmod_mpolyun_init(T, bits, ctx);

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
        nmod_mpolyun_intp_reduce_2sm_mpolyun(Aevalp, Aevalm,
                                                  base->A, var, alphapow, ctx);
        nmod_mpolyun_intp_reduce_2sm_mpolyun(Bevalp, Bevalm,
                                                  base->B, var, alphapow, ctx);
        FLINT_ASSERT(Aevalp->length > 0);
        FLINT_ASSERT(Aevalm->length > 0);
        FLINT_ASSERT(Bevalp->length > 0);
        FLINT_ASSERT(Bevalm->length > 0);

        success = nmod_mpolyun_gcd_brown_smprime(
                    Gevalp, Abarevalp, Bbarevalp, Aevalp, Bevalp, var - 1,
                                                             ctx, base->I, Sp);
        success = success && nmod_mpolyun_gcd_brown_smprime(
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

        if (   nmod_mpolyun_is_nonzero_nmod(Gevalp, ctx)
            || nmod_mpolyun_is_nonzero_nmod(Gevalm, ctx))
        {
            base->gcd_is_one = 1;
            break;
        }

        if (Gevalp->exps[0] != Gevalm->exps[0])
        {
            continue;
        }
        if (   nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0)
            != nmod_poly_degree((Gevalm->coeffs + 0)->coeffs + 0))
        {
            continue;
        }
        if (!mpoly_monomial_equal((Gevalp->coeffs + 0)->exps + N*0,
                                  (Gevalm->coeffs + 0)->exps + N*0, N))
        {
            continue;
        }

        /* the Geval have matching degrees */
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            int cmp = 0;
            FLINT_ASSERT(arg->G->length > 0);
            if (arg->G->exps[0] != Gevalp->exps[0])
            {
                cmp = arg->G->exps[0] > Gevalp->exps[0] ? 1 : -1;
            }
            if (cmp == 0)
            {
                slong k = nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0);
                cmp = mpoly_monomial_cmp_nomask_extra(
                       (arg->G->coeffs + 0)->exps + N*0,
                       (Gevalp->coeffs + 0)->exps + N*0, N, offset, k << shift);
            }

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
        temp = nmod_mpolyn_leadcoeff(Gevalp->coeffs + 0, ctx);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        temp = nmod_mul(gammaevalp, temp, ctx->ffinfo->mod);
        nmod_mpolyun_scalar_mul_nmod(Gevalp, temp, ctx);
        temp = nmod_mpolyn_leadcoeff(Gevalm->coeffs + 0, ctx);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        temp = nmod_mul(gammaevalm, temp, ctx->ffinfo->mod);
        nmod_mpolyun_scalar_mul_nmod(Gevalm, temp, ctx);
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(arg->modulus, alpha);
            FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(arg->modulus, ctx->ffinfo->mod.n - alpha));
            temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
            temp = nmod_add(temp, temp, ctx->ffinfo->mod);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(arg->modulus, arg->modulus, temp);
            nmod_mpolyun_intp_crt_2sm_mpolyun(&ldeg, arg->G, T,
                             Gevalp, Gevalm, var, arg->modulus, alphapow, ctx);
            nmod_mpolyun_intp_crt_2sm_mpolyun(&ldeg, arg->Abar, T,
                       Abarevalp, Abarevalm, var, arg->modulus, alphapow, ctx);
            nmod_mpolyun_intp_crt_2sm_mpolyun(&ldeg, arg->Bbar, T,
                       Bbarevalp, Bbarevalm, var, arg->modulus, alphapow, ctx);
        }
        else
        {
            nmod_mpolyun_intp_lift_2sm_mpolyun(&ldeg, arg->G,
                                              Gevalp, Gevalm, var, alpha, ctx);
            nmod_mpolyun_intp_lift_2sm_mpolyun(&ldeg, arg->Abar,
                                        Abarevalp, Abarevalm, var, alpha, ctx);
            nmod_mpolyun_intp_lift_2sm_mpolyun(&ldeg, arg->Bbar,
                                        Bbarevalp, Bbarevalm, var, alpha, ctx);
        }
        temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus2, arg->modulus, temp);
        nmod_poly_shift_left(arg->modulus, arg->modulus, 2);
        nmod_poly_sub(arg->modulus, arg->modulus, modulus2);
    }

    nmod_poly_clear(modulus2);
    nmod_poly_clear(alphapow);

    nmod_mpolyun_clear(Aevalp, ctx);
    nmod_mpolyun_clear(Bevalp, ctx);
    nmod_mpolyun_clear(Gevalp, ctx);
    nmod_mpolyun_clear(Abarevalp, ctx);
    nmod_mpolyun_clear(Bbarevalp, ctx);
    nmod_mpolyun_clear(Aevalm, ctx);
    nmod_mpolyun_clear(Bevalm, ctx);
    nmod_mpolyun_clear(Gevalm, ctx);
    nmod_mpolyun_clear(Abarevalm, ctx);
    nmod_mpolyun_clear(Bbarevalm, ctx);
    nmod_mpolyun_clear(T, ctx);

    nmod_poly_stack_clear(Sp);
}



typedef struct
{
    volatile int idx;
    volatile slong G_exp, Abar_exp, Bbar_exp;
    pthread_mutex_t mutex;
    const nmod_mpoly_ctx_struct * ctx;
    nmod_poly_multi_crt_t CRT;
    nmod_mpolyun_struct ** gptrs, ** abarptrs, ** bbarptrs;
    ulong num_threads;
}
_joinbase_struct;

typedef _joinbase_struct _joinbase_t[1];

typedef struct
{
    _joinbase_struct * base;
    nmod_mpolyun_t G, Abar, Bbar;
    slong G_lastdeg, Abar_lastdeg, Bbar_lastdeg;
}
_joinworker_arg_struct;

/* Join the images of G, Abar, and Bbar in each of the split threads */
static void _joinworker(void * varg)
{
    _joinworker_arg_struct * arg = (_joinworker_arg_struct *) varg;
    _joinbase_struct * base = arg->base;
    slong t, our_G_exp, our_Abar_exp, our_Bbar_exp;
    slong count = base->num_threads;
    nmod_mpolyn_struct ** coeffs;
    nmod_poly_struct ** input;
    slong * start;
    slong j, ls;
    nmod_poly_struct * output;

    /* should already be zero, but just make sure */
    arg->G->length = 0;
    arg->Abar->length = 0;
    arg->Bbar->length = 0;

    /* some temp space */
    coeffs = (nmod_mpolyn_struct **) flint_malloc(count * sizeof(nmod_mpolyn_struct *));
    input = (nmod_poly_struct **) flint_malloc(count * sizeof(nmod_poly_struct *));
    start = (slong *) flint_malloc(count * sizeof(slong));
    ls = _nmod_poly_multi_crt_local_size(base->CRT);
    output = (nmod_poly_struct *) flint_malloc(ls*sizeof(nmod_poly_struct));
    for (j = 0; j < ls; j++)
    {
        nmod_poly_init_mod(output + j, base->ctx->ffinfo->mod);
    }

    while (1)
    {
        /* get exponent of either G, Abar, or Bbar to start working on */
        pthread_mutex_lock(&base->mutex);
        our_G_exp = base->G_exp;
        our_Abar_exp = base->Abar_exp;
        our_Bbar_exp = base->Bbar_exp;
        if (our_G_exp >= 0)
        {
            base->G_exp = our_G_exp - 1;
        }
        else if (our_Abar_exp >= 0)
        {
            base->Abar_exp = our_Abar_exp - 1;
        }
        else if (our_Bbar_exp >= 0)
        {
            base->Bbar_exp = our_Bbar_exp - 1;
        }
        pthread_mutex_unlock(&base->mutex);

        if (our_G_exp >= 0)
        {
            t = _nmod_mpolyun_crt_exp(base->CRT, arg->G, our_G_exp,
                                 base->gptrs, base->num_threads,
                                      coeffs, input, start, output, base->ctx);
            arg->G_lastdeg = FLINT_MAX(arg->G_lastdeg, t);
        }
        else if (our_Abar_exp >= 0)
        {
            t = _nmod_mpolyun_crt_exp(base->CRT, arg->Abar, our_Abar_exp,
                              base->abarptrs, base->num_threads,
                                      coeffs, input, start, output, base->ctx);
            arg->Abar_lastdeg = FLINT_MAX(arg->Abar_lastdeg, t);
        }
        else if (our_Bbar_exp >= 0)
        {
            t = _nmod_mpolyun_crt_exp(base->CRT, arg->Bbar, our_Bbar_exp,
                              base->bbarptrs, base->num_threads,
                                      coeffs, input, start, output, base->ctx);
            arg->Bbar_lastdeg = FLINT_MAX(arg->Bbar_lastdeg, t);
        }
        else
        {
            goto cleanup;
        }
    }

cleanup:

    for (j = 0; j < ls; j++)
    {
        nmod_poly_clear(output + j);
    }
    flint_free(output);
    flint_free(coeffs);
    flint_free(input);
    flint_free(start);

    return;
}

/*
    A = B[0] + ... + B[num_threads - 1]
    A and the B[i] are in Fp[X][x_0, ..., x_(var-1)][var]
    The B[i] have distinct exponents on X, so this is just a top level merge.
    The inputs B[i] are clobbered.
*/
static void _final_join(
    nmod_mpolyun_t A,
    nmod_mpolyun_struct ** B,
    slong num_threads,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Ai, total_length;
    slong * starts;
    TMP_INIT;

    TMP_START;
    starts = (slong *) TMP_ALLOC(num_threads*sizeof(slong));
    total_length = 0;
    for (i = 0; i < num_threads; i++)
    {
        starts[i] = 0;
        total_length += B[i]->length;
    }

    nmod_mpolyun_fit_length(A, total_length, ctx);
    Ai = 0;
    while (1)
    {
        slong max_pos = -WORD(1);
        slong max_exp = -WORD(1);
        for (i = 0; i < num_threads; i++)
        {
            if (starts[i] < B[i]->length
                          && (slong)(B[i]->exps[starts[i]]) > max_exp)
            {
                max_pos = i;
                max_exp = B[i]->exps[starts[i]];
            }
        }
        if (max_pos < 0)
        {
            break;
        }
        A->exps[Ai] = max_exp;
        nmod_mpolyn_swap(A->coeffs + Ai, B[max_pos]->coeffs + starts[max_pos]);
        starts[max_pos]++;
        Ai++;
    }
    A->length = Ai;
    FLINT_ASSERT(Ai == total_length);
    FLINT_ASSERT(nmod_mpolyun_is_canonical(A, ctx));

    TMP_END;
}

/*
    Do same as nmod_mpolyun_gcd_brown_smprime but use the threads in
        handles[0], ..., handles[num_handles - 1]
    num_handles is allowed to be zero.
*/
int nmod_mpolyun_gcd_brown_smprime_threaded(
    nmod_mpolyun_t G,
    nmod_mpolyun_t Abar,
    nmod_mpolyun_t Bbar,
    nmod_mpolyun_t A,
    nmod_mpolyun_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int divisibility_test = 0; /* 1: by G, 2: by Abar, 3: by Bbar */
    slong i;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong num_threads;
    int success;
    ulong bound, best_est;
    slong g_stab_est, abar_stab_est, bbar_stab_est, upper_limit;
    mp_limb_t alpha;
    slong deggamma, ldegGs, ldegAbars, ldegBbars, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    slong Gexp0, Abarexp0, Bbarexp0;
    nmod_poly_struct ** mptrs;
    nmod_mpolyun_struct ** gptrs, ** abarptrs, ** bbarptrs;
    _splitworker_arg_struct * splitargs;
    _splitbase_t splitbase;
    _joinworker_arg_struct * joinargs;
    _joinbase_t joinbase;
    nmod_poly_t t1;
    nmod_mpolyun_t T1, T2;

    nmod_poly_init_mod(t1, ctx->ffinfo->mod);
    nmod_mpolyun_init(T1, bits, ctx);
    nmod_mpolyun_init(T2, bits, ctx);

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(A, ctx),
                         nmod_mpolyun_leadcoeff_poly(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);
    best_est = bound;

    upper_limit = mpoly_gcd_info_get_brown_upper_limit(I, var + 1, bound);

    if (I != NULL && I->Gdeflate_deg_bounds_are_nice)
    {
        slong k = I->brown_perm[var + 1];

        FLINT_ASSERT(var + 1 < I->mvars);
        
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
    gptrs = (nmod_mpolyun_struct **) flint_malloc(
                                    num_threads*sizeof(nmod_mpolyun_struct *));
    abarptrs = (nmod_mpolyun_struct **) flint_malloc(
                                    num_threads*sizeof(nmod_mpolyun_struct *));
    bbarptrs = (nmod_mpolyun_struct **) flint_malloc(
                                    num_threads*sizeof(nmod_mpolyun_struct *));
    mptrs = (nmod_poly_struct **) flint_malloc(
                                       num_threads*sizeof(nmod_poly_struct *));
    splitargs = (_splitworker_arg_struct *) flint_malloc(
                                  num_threads*sizeof(_splitworker_arg_struct));
    for (i = 0; i < num_threads; i++)
    {
        nmod_mpolyun_init(splitargs[i].G, bits, ctx);
        nmod_mpolyun_init(splitargs[i].Abar, bits, ctx);
        nmod_mpolyun_init(splitargs[i].Bbar, bits, ctx);
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
        splitargs[i].required_images = FLINT_MAX(ri, WORD(1));
    }

    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i],
                  var == 0 ? _splitworker_bivar : _splitworker, &splitargs[i]);
    }
    (var == 0 ? _splitworker_bivar : _splitworker)(&splitargs[num_threads - 1]);
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    if (splitbase->gcd_is_one)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
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
    Gexp0 = gptrs[0]->exps[0];
    Abarexp0 = A->exps[0] - Gexp0;
    Bbarexp0 = B->exps[0] - Gexp0;
    for (i = 0; i < num_threads; i++)
    {
        if (gptrs[i]->exps[0] != Gexp0
                      || !mpoly_monomial_equal(
                                  (gptrs[i]->coeffs + 0)->exps + N*0,
                                  (gptrs[0]->coeffs + 0)->exps + N*0, N))
        {
            /* very unlucky - could try again or just fail */
            goto compute_split;
        }
        FLINT_ASSERT(splitargs[i].Abar->exps[0] == Abarexp0);
        FLINT_ASSERT(splitargs[i].Bbar->exps[0] == Bbarexp0);
    }

    nmod_poly_multi_crt_init(joinbase->CRT);
    success = nmod_poly_multi_crt_precompute_p(joinbase->CRT,
                        (const nmod_poly_struct * const *) mptrs, num_threads);
    FLINT_ASSERT(success);

    joinbase->num_threads = num_threads;
    joinbase->gptrs = gptrs;
    joinbase->abarptrs = abarptrs;
    joinbase->bbarptrs = bbarptrs;
    joinbase->G_exp = Gexp0;
    joinbase->Abar_exp = Abarexp0;
    joinbase->Bbar_exp = Bbarexp0;
    joinbase->ctx = ctx;
    pthread_mutex_init(&joinbase->mutex, NULL);

    joinargs = (_joinworker_arg_struct *) flint_malloc(
                                   num_threads*sizeof(_joinworker_arg_struct));

    for (i = 0; i < num_threads; i++)
    {
        joinargs[i].base = joinbase;
        joinargs[i].G_lastdeg = -WORD(1);
        joinargs[i].Abar_lastdeg = -WORD(1);
        joinargs[i].Bbar_lastdeg = -WORD(1);
        nmod_mpolyun_init(joinargs[i].G, bits, ctx);
        nmod_mpolyun_init(joinargs[i].Abar, bits, ctx);
        nmod_mpolyun_init(joinargs[i].Bbar, bits, ctx);
    }
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool,
                                        handles[i], _joinworker, joinargs + i);
    }
    _joinworker(joinargs + num_threads - 1);
    for (i = 0; i + 1 < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }
    pthread_mutex_destroy(&joinbase->mutex);

    ldegGs = ldegAbars = ldegBbars = -WORD(1);
    for (i = 0; i < num_threads; i++)
    {
        ldegGs = FLINT_MAX(ldegGs, joinargs[i].G_lastdeg);
        ldegAbars = FLINT_MAX(ldegAbars, joinargs[i].Abar_lastdeg);
        ldegBbars = FLINT_MAX(ldegBbars, joinargs[i].Bbar_lastdeg);
    }

    /* reuse gptrs, abarpts, bbarpts for final trivial join */
    for (i = 0; i < num_threads; i++)
    {
        gptrs[i] = joinargs[i].G;
        abarptrs[i] = joinargs[i].Abar;
        bbarptrs[i] = joinargs[i].Bbar;
    }
    _final_join(G, gptrs, num_threads, ctx);
    _final_join(Abar, abarptrs, num_threads, ctx);
    _final_join(Bbar, bbarptrs, num_threads, ctx);

    /* free join data */
    nmod_poly_multi_crt_clear(joinbase->CRT);
    for (i = 0; i < num_threads; i++)
    {
        nmod_mpolyun_clear(joinargs[i].G, ctx);
        nmod_mpolyun_clear(joinargs[i].Abar, ctx);
        nmod_mpolyun_clear(joinargs[i].Bbar, ctx);
    }
    flint_free(joinargs);

    if (divisibility_test == 1)
    {
        nmod_mpolyun_content_last(t1, G, ctx);
        nmod_mpolyun_divexact_last(G, t1, ctx);
        success =            nmod_mpolyun_divides(T1, A, G, ctx);
        success = success && nmod_mpolyun_divides(T2, B, G, ctx);
        if (success)
        {
            ulong temp;
            nmod_mpolyun_swap(T1, Abar);
            nmod_mpolyun_swap(T2, Bbar);
successful_fix_lc:
            temp = nmod_mpolyun_leadcoeff(G, ctx);
            nmod_mpolyun_scalar_mul_nmod(Abar, temp, ctx);
            nmod_mpolyun_scalar_mul_nmod(Bbar, temp, ctx);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_mpolyun_scalar_mul_nmod(G, temp, ctx);
            goto successful_put_content;
        }
    }
    else if (divisibility_test == 2)
    {
        nmod_mpolyun_content_last(t1, Abar, ctx);
        nmod_mpolyun_divexact_last(Abar, t1, ctx);
        success =            nmod_mpolyun_divides(T1, A, Abar, ctx);
        success = success && nmod_mpolyun_divides(T2, B, T1, ctx);
        if (success)
        {
            nmod_mpolyun_swap(T1, G);
            nmod_mpolyun_swap(T2, Bbar);
            goto successful_fix_lc;
        }
    }
    else if (divisibility_test == 3)
    {
        nmod_mpolyun_content_last(t1, Bbar, ctx);
        nmod_mpolyun_divexact_last(Bbar, t1, ctx);
        success =            nmod_mpolyun_divides(T1, B, Bbar, ctx);
        success = success && nmod_mpolyun_divides(T2, A, T1, ctx);
        if (success)
        {
            nmod_mpolyun_swap(T1, G);
            nmod_mpolyun_swap(T2, Abar);
            goto successful_fix_lc;
        }        
    }
    else /* divisibility test == 0 */
    {
        if (   deggamma + ldegA == ldegGs + ldegAbars
            && deggamma + ldegB == ldegGs + ldegBbars)
        {
            goto successful;
        }
    }

    /* divisibility test failed - try again */
    best_est = bound;
    divisibility_test = 0;
    goto compute_split;

successful:

    nmod_mpolyun_content_last(t1, G, ctx);
    nmod_mpolyun_divexact_last(G, t1, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup_split:

    for (i = 0; i < num_threads; i++)
    {
        nmod_mpolyun_clear(splitargs[i].G, ctx);
        nmod_mpolyun_clear(splitargs[i].Abar, ctx);
        nmod_mpolyun_clear(splitargs[i].Bbar, ctx);
        nmod_poly_clear(splitargs[i].modulus);
    }

    flint_free(gptrs);
    flint_free(abarptrs);
    flint_free(bbarptrs);
    flint_free(mptrs);
    flint_free(splitargs);

cleanup:

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);
    nmod_poly_clear(gamma);

    nmod_poly_clear(t1);
    nmod_mpolyun_clear(T1, ctx);
    nmod_mpolyun_clear(T2, ctx);

    return success;
}


typedef struct
{
    nmod_mpolyun_struct * Pn;
    const nmod_mpoly_ctx_struct * uctx;
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

    nmod_mpoly_to_mpolyun_perm_deflate(arg->Pn, arg->uctx,
                    arg->P, arg->ctx, arg->perm, arg->shift, arg->stride,
                                               arg->handles, arg->num_handles);
}


int nmod_mpoly_gcd_brown_threaded(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    slong thread_limit)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    flint_bitcnt_t ABbits;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyun_t An, Bn, Gn, Abarn, Bbarn;
    thread_pool_handle * handles;
    slong num_handles;

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

    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyun_init(An, ABbits, uctx);
    nmod_mpolyun_init(Bn, ABbits, uctx);
    nmod_mpolyun_init(Gn, ABbits, uctx);
    nmod_mpolyun_init(Abarn, ABbits, uctx);
    nmod_mpolyun_init(Bbarn, ABbits, uctx);

    handles = NULL;
    num_handles = 0;
    if (global_thread_pool_initialized)
    {
        slong max_num_handles = thread_pool_get_size(global_thread_pool);
        max_num_handles = FLINT_MIN(thread_limit - 1, max_num_handles);
        if (max_num_handles > 0)
        {
            handles = (thread_pool_handle *) flint_malloc(
                               max_num_handles*sizeof(thread_pool_handle));
            num_handles = thread_pool_request(global_thread_pool,
                                                 handles, max_num_handles);
        }
    }

    /* convert inputs */
    if (num_handles > 0)
    {
        slong m = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertn_arg_t arg;

        FLINT_ASSERT(m >= 0);
        FLINT_ASSERT(m < num_handles);

        arg->Pn = Bn;
        arg->uctx = uctx;
        arg->P = B;
        arg->ctx = ctx;
        arg->perm = perm;
        arg->shift = shift;
        arg->stride = stride;
        arg->handles = handles + (m + 1);
        arg->num_handles = num_handles - (m + 1);

        thread_pool_wake(global_thread_pool, handles[m], _worker_convertn, arg);

        nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx,
                                          perm, shift, stride, handles + 0, m);

        thread_pool_wait(global_thread_pool, handles[m]);
    }
    else
    {
        nmod_mpoly_to_mpolyun_perm_deflate(An, uctx,
                                         A, ctx, perm, shift, stride, NULL, 0);
        nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx,
                                         B, ctx, perm, shift, stride, NULL, 0);
    }

    /* calculate gcd */
    success = nmod_mpolyun_gcd_brown_smprime_threaded(
                          Gn, Abarn, Bbarn, An, Bn, uctx->minfo->nvars - 1,
                                             uctx, NULL, handles, num_handles);

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_give_back(global_thread_pool, handles[i]);
    }

    if (handles)
        flint_free(handles);

    if (!success)
    {
        nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx,
                                                 perm, shift, stride, NULL, 0);
        nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx, B, ctx,
                                                 perm, shift, stride, NULL, 0);
        success = nmod_mpolyun_gcd_brown_lgprime(Gn, Abarn, Bbarn,
                                         An, Bn, uctx->minfo->nvars - 1, uctx);
    }

    if (success)
    {
        nmod_mpoly_from_mpolyun_perm_inflate(G, ABbits, ctx,
                                                Gn, uctx, perm, shift, stride);
        nmod_mpoly_make_monic(G, G, ctx);
    }

    nmod_mpolyun_clear(An, uctx);
    nmod_mpolyun_clear(Bn, uctx);
    nmod_mpolyun_clear(Gn, uctx);
    nmod_mpolyun_clear(Abarn, uctx);
    nmod_mpolyun_clear(Bbarn, uctx);
    nmod_mpoly_ctx_clear(uctx);

cleanup1:

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}
