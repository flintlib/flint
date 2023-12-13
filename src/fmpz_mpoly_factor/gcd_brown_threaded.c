/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpq.h"
#include "fmpq_vec.h"

typedef struct
{
    volatile int gcd_is_one;
    volatile mp_limb_t p;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    fmpz_t gamma;
    const fmpz_mpoly_ctx_struct * ctx;
    fmpz_mpoly_struct * A, * B;
    ulong num_threads;
    slong var;
    const mpoly_gcd_info_struct * I;
}
_splitbase_struct;

typedef _splitbase_struct _splitbase_t[1];

typedef struct
{
    slong idx;
    _splitbase_struct * base;
    fmpz_mpoly_t G, Abar, Bbar;
    fmpz_t modulus;
    slong image_count;
    slong required_images;
    thread_pool_handle master_handle;
    slong num_handles;
    thread_pool_handle * worker_handles;

    nmod_mpoly_ctx_t pctx;
    nmod_mpolyn_t Ap, Bp, Gp, Abarp, Bbarp;
    fmpz_mpoly_t T, T1, T2;
}
_splitworker_arg_struct;



/* worker for reducing polynomial over ZZ to polynomial over Fp */
static void _reduce_Bp_worker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    fmpz_mpoly_interp_reduce_p_mpolyn(arg->Bp, arg->pctx, arg->base->B,
                                                               arg->base->ctx);
}

/* workers for crt'ing polynomials */
static void _join_Abar_worker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    if (!fmpz_is_one(arg->modulus))
        fmpz_mpoly_interp_crt_p_mpolyn(arg->Abar, arg->T1, arg->base->ctx,
                                          arg->modulus, arg->Abarp, arg->pctx);
    else
        fmpz_mpoly_interp_lift_p_mpolyn(arg->Abar, arg->base->ctx,
                                                        arg->Abarp, arg->pctx);
}

static void _join_Bbar_worker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    if (!fmpz_is_one(arg->modulus))
        fmpz_mpoly_interp_crt_p_mpolyn(arg->Bbar, arg->T2, arg->base->ctx,
                                          arg->modulus, arg->Bbarp, arg->pctx);
    else
        fmpz_mpoly_interp_lift_p_mpolyn(arg->Bbar, arg->base->ctx,
                                                        arg->Bbarp, arg->pctx);
}


static void _splitworker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const fmpz_mpoly_ctx_struct * ctx = base->ctx;
    flint_bitcnt_t bits = base->A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong offset, shift;
    int success;
    mp_limb_t p, gammared;
    nmod_poly_stack_t Sp;

    mpoly_gen_offset_shift_sp(&offset, &shift,
                                      ctx->minfo->nvars - 1, bits, ctx->minfo);

    fmpz_one(arg->modulus);
    arg->image_count = 0;

    nmod_mpoly_ctx_init(arg->pctx, ctx->minfo->nvars, ORD_LEX, 2);
    nmod_poly_stack_init(Sp, bits, arg->pctx);
    nmod_mpolyn_init(arg->Ap, bits, arg->pctx);
    nmod_mpolyn_init(arg->Bp, bits, arg->pctx);
    nmod_mpolyn_init(arg->Gp, bits, arg->pctx);
    nmod_mpolyn_init(arg->Abarp, bits, arg->pctx);
    nmod_mpolyn_init(arg->Bbarp, bits, arg->pctx);
    fmpz_mpoly_init3(arg->T, 0, bits, ctx);
    fmpz_mpoly_init3(arg->T1, 0, bits, ctx);
    fmpz_mpoly_init3(arg->T2, 0, bits, ctx);

    while (arg->image_count < arg->required_images)
    {
        /* get prime */
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
	p = base->p;
        if (p >= UWORD_MAX_PRIME)
        {
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(&base->mutex);
#endif
            break;
        }
        p = n_nextprime(base->p, 1);
        base->p = p;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif

        /* make sure reduction does not kill both lc(A) and lc(B) */
        gammared = fmpz_fdiv_ui(base->gamma, p);
        if (gammared == 0)
        {
            continue;
        }

        nmod_mpoly_ctx_set_modulus(arg->pctx, p);

        /* the unfortunate nmod poly's store their own context :( */
        nmod_poly_stack_set_ctx(Sp, arg->pctx);
        nmod_mpolyn_set_mod(arg->Ap, arg->pctx->mod);
        nmod_mpolyn_set_mod(arg->Bp, arg->pctx->mod);
        nmod_mpolyn_set_mod(arg->Gp, arg->pctx->mod);
        nmod_mpolyn_set_mod(arg->Abarp, arg->pctx->mod);
        nmod_mpolyn_set_mod(arg->Bbarp, arg->pctx->mod);

        /* reduce to Fp and calculate an image gcd */
        if (arg->num_handles > 0)
        {
            thread_pool_wake(global_thread_pool, arg->worker_handles[0], 0,
                                                       _reduce_Bp_worker, arg);

            fmpz_mpoly_interp_reduce_p_mpolyn(arg->Ap, arg->pctx, base->A, ctx);

            thread_pool_wait(global_thread_pool, arg->worker_handles[0]);

            success = nmod_mpolyn_gcd_brown_smprime_threaded_pool(
                                  arg->Gp, arg->Abarp, arg->Bbarp,
                   arg->Ap, arg->Bp, ctx->minfo->nvars - 1, arg->pctx, base->I,
                                       arg->worker_handles, arg->num_handles);
        }
        else
        {
            /* reduction should kill neither A nor B */
            fmpz_mpoly_interp_reduce_p_mpolyn(arg->Ap, arg->pctx, base->A, ctx);
            fmpz_mpoly_interp_reduce_p_mpolyn(arg->Bp, arg->pctx, base->B, ctx);
            FLINT_ASSERT(arg->Ap->length > 0);
            FLINT_ASSERT(arg->Bp->length > 0);
            success = nmod_mpolyn_gcd_brown_smprime(
                            arg->Gp, arg->Abarp, arg->Bbarp, arg->Ap, arg->Bp,
                                ctx->minfo->nvars - 1, arg->pctx, base->I, Sp);
        }

        if (!success)
        {
            continue;
        }

        FLINT_ASSERT(arg->Gp->length > 0);
        FLINT_ASSERT(arg->Abarp->length > 0);
        FLINT_ASSERT(arg->Bbarp->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
            break;
        }

        if (nmod_mpolyn_is_nonzero_nmod(arg->Gp, arg->pctx))
        {
            base->gcd_is_one = 1;
            break;
        }

        if (!fmpz_is_one(arg->modulus))
        {
            int cmp = 0;
            slong k;

            FLINT_ASSERT(arg->G->length > 0);

            k = n_poly_degree(arg->Gp->coeffs + 0);
            cmp = mpoly_monomial_cmp_nomask_extra(arg->G->exps + N*0,
                                   arg->Gp->exps + N*0, N, offset, k << shift);
            if (cmp < 0)
            {
                continue;
            }
            else if (cmp > 0)
            {
                fmpz_one(arg->modulus);
                arg->image_count = 0;
            }
        }

        FLINT_ASSERT(1 == nmod_mpolyn_leadcoeff(arg->Gp, arg->pctx));
        nmod_mpolyn_scalar_mul_nmod(arg->Gp, gammared, arg->pctx);

        /* crt image gcd */
        if (arg->num_handles > 0)
        {
            thread_pool_wake(global_thread_pool, arg->worker_handles[0], 0,
                                                       _join_Abar_worker, arg);
            if (arg->num_handles > 1)
            {
                thread_pool_wake(global_thread_pool, arg->worker_handles[1], 0,
                                                       _join_Bbar_worker, arg);
            }
            else
            {
                _join_Bbar_worker(arg);
            }

            if (!fmpz_is_one(arg->modulus))
                fmpz_mpoly_interp_crt_p_mpolyn(arg->G, arg->T, ctx, arg->modulus,
                                                           arg->Gp, arg->pctx);
            else
                fmpz_mpoly_interp_lift_p_mpolyn(arg->G, ctx, arg->Gp, arg->pctx);

            thread_pool_wait(global_thread_pool, arg->worker_handles[0]);
            if (arg->num_handles > 1)
                thread_pool_wait(global_thread_pool, arg->worker_handles[1]);
        }
        else
        {
            if (!fmpz_is_one(arg->modulus))
            {
                fmpz_mpoly_interp_crt_p_mpolyn(arg->G, arg->T, ctx,
                                             arg->modulus, arg->Gp, arg->pctx);
                fmpz_mpoly_interp_crt_p_mpolyn(arg->Abar, arg->T, ctx,
                                          arg->modulus, arg->Abarp, arg->pctx);
                fmpz_mpoly_interp_crt_p_mpolyn(arg->Bbar, arg->T, ctx,
                                          arg->modulus, arg->Bbarp, arg->pctx);
            }
            else
            {
                fmpz_mpoly_interp_lift_p_mpolyn(arg->G, ctx,
                                                           arg->Gp, arg->pctx);
                fmpz_mpoly_interp_lift_p_mpolyn(arg->Abar, ctx,
                                                        arg->Abarp, arg->pctx);
                fmpz_mpoly_interp_lift_p_mpolyn(arg->Bbar, ctx,
                                                        arg->Bbarp, arg->pctx);
            }
        }

        fmpz_mul_ui(arg->modulus, arg->modulus, p);
        arg->image_count++;
    }

    fmpz_mpoly_clear(arg->T, ctx);
    fmpz_mpoly_clear(arg->T1, ctx);
    fmpz_mpoly_clear(arg->T2, ctx);

    nmod_mpolyn_clear(arg->Ap, arg->pctx);
    nmod_mpolyn_clear(arg->Bp, arg->pctx);
    nmod_mpolyn_clear(arg->Gp, arg->pctx);
    nmod_mpolyn_clear(arg->Abarp, arg->pctx);
    nmod_mpolyn_clear(arg->Bbarp, arg->pctx);
    nmod_poly_stack_clear(Sp);
    nmod_mpoly_ctx_clear(arg->pctx);
}

typedef struct
{
    slong hint_start, hint_stop;
    ulong * left_exp, * right_exp;
    fmpz_mpoly_t poly;
    fmpz_t maxcoeff, sumcoeff;
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
    fmpz_mpoly_struct * const * B,
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

/* NOTE: Here `input' has to be an array of size `count'.  It will be
 * overwritten, and its initial element values do not matter.  It will contain
 * dangling pointers after returning, that is, the values `input[ix]' will be
 * invalid after returning.  Hence, we silence warnings about dangling pointers. */
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wdangling-pointer"
#endif
static slong _fmpz_mpoly_crt(
    const fmpz_multi_CRT_t P,
    _joinworker_arg_struct * S,
    fmpz_mpoly_struct * const * B,
    slong count,
    fmpz * output,
    fmpz ** input,
    const fmpz_mpoly_ctx_t ctx)
{
    int cmp;
    slong N = mpoly_words_per_exp_sp(S->poly->bits, ctx->minfo);
    slong lastdegree;
    slong Ai;
    slong j, k;
    slong * start, * stop;
    fmpz * ins;
    fmpz_t zero, max, sum;
    fmpz_mpoly_t A;
    const ulong * exp_left = S->left_exp;
    const ulong * exp_right = S->right_exp;
    TMP_INIT;

    FLINT_ASSERT(count > 0);
    FLINT_ASSERT(count == P->moduli_count);

    *A = *S->poly;

    TMP_START;

    ins = TMP_ARRAY_ALLOC(count, fmpz);
    start = TMP_ARRAY_ALLOC(2*count, slong);
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

#ifdef FLINT_WANT_ASSERT
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

    fmpz_init(zero);
    fmpz_init(max);
    fmpz_init(sum);

    Ai = 0;
    lastdegree = -WORD(1);
    while (1)
    {
        fmpz_mpoly_fit_length(A, Ai + 1, ctx);

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

        k = 0;
        do {
            ins[k] = *input[k];
        } while (++k < count);

        _fmpz_multi_CRT_precomp(output, P, ins, 1);
        fmpz_swap(A->coeffs + Ai, output + 0);
        cmp = fmpz_sgn(A->coeffs + Ai);
        if (cmp != 0)
        {
            if (fmpz_cmpabs(max, A->coeffs + Ai) < 0)
                fmpz_abs(max, A->coeffs + Ai);
            if (cmp > 0)
                fmpz_add(sum, sum, A->coeffs + Ai);
            else
                fmpz_sub(sum, sum, A->coeffs + Ai);
            Ai++;
        }
    }
    A->length = Ai;

    fmpz_swap(S->maxcoeff, max);
    fmpz_swap(S->sumcoeff, sum);

    fmpz_clear(zero);
    fmpz_clear(max);
    fmpz_clear(sum);

    TMP_END;

    *S->poly = *A;

    return lastdegree;
}
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif


typedef struct
{
    volatile int idx;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    const fmpz_mpoly_ctx_struct * ctx;
    fmpz_multi_CRT_t CRT;
    fmpz_mpoly_struct ** gptrs, ** abarptrs, ** bbarptrs;
    fmpz_mpoly_struct * G, * Abar, * Bbar;
    _joinworker_arg_struct * chunks;
    slong chunks_length;
    ulong num_images;
}
_joinbase_struct;

typedef _joinbase_struct _joinbase_t[1];

typedef struct
{
    _joinbase_struct * base;
    slong thread_idx;
}
_njoinworker_arg_struct;

static void _joinworker(void * varg)
{
    _njoinworker_arg_struct * arg = (_njoinworker_arg_struct *) varg;
    _joinbase_struct * base = arg->base;
    fmpz ** input;
    fmpz * output;
    slong i, ls = base->CRT->localsize;
    TMP_INIT;

    TMP_START;

    input = (fmpz **) TMP_ALLOC(base->num_images * sizeof(fmpz *));
    output = (fmpz *) TMP_ALLOC(ls*sizeof(fmpz));
    for (i = 0; i < ls; i++)
        fmpz_init(output + i);

    while (1)
    {
        /* get exponent of either G, Abar, or Bbar to start working on */
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
	i = base->idx;
        base->idx = i + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif

        if (i >= base->chunks_length)
        {
            goto cleanup;
        }

        base->chunks[i].thread_idx = arg->thread_idx;

        if (base->chunks[i].GAB == 0)
        {
            _fmpz_mpoly_crt(base->CRT, base->chunks + i, base->gptrs,
                                  base->num_images, output, input, base->ctx);
        }
        else if (base->chunks[i].GAB == 1)
        {
            _fmpz_mpoly_crt(base->CRT, base->chunks + i, base->abarptrs,
                                  base->num_images, output, input, base->ctx);
        }
        else
        {
            FLINT_ASSERT(base->chunks[i].GAB == 2);

            _fmpz_mpoly_crt(base->CRT, base->chunks + i, base->bbarptrs,
                                  base->num_images, output, input, base->ctx);
        }
    }

cleanup:

    for (i = 0; i < ls; i++)
        fmpz_clear(output + i);

    TMP_END;

    return;
}


static void _finaljoinworker(void * varg)
{
    _njoinworker_arg_struct * arg = (_njoinworker_arg_struct *) varg;
    _joinbase_struct * base = arg->base;
    const fmpz_mpoly_ctx_struct * ctx = base->ctx;
    flint_bitcnt_t bits = base->G->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong i, j;
    slong source_len;
    ulong * source_exps;
    fmpz * source_coeffs;
    slong Ti;
    ulong * Texps;
    fmpz * Tcoeffs;

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
            fmpz_swap(Tcoeffs + Ti + j, source_coeffs + j);
    }
}

/*
    Set 1 <= l <= min(n, m) and fractions v + 0, ..., v + l - 1
    More comments in usage

    example operation on input n = 10, m = 16

    gcd is 2:
        5/8, 5/8
    split first 5/8:
        2/3, 5/8, 3/5,
    split first 5/8:
        2/3, 2/3, 3/5, 3/5,
    split first 3/5:
        2/3, 2/3, 2/3, 3/5, 1/2
    split first 3/5:
        2/3, 2/3, 2/3, 2/3, 1/2, 1/2

    The maximum fraction is now 0.666 which is not much bigger than n/m
*/
static slong _divide_master_threads(fmpq * v, slong n, slong m)
{
    slong l, i;
    double score_threashold;
    fmpq_t left, right;

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(m > 0);

    fmpq_init(left);
    fmpq_init(right);

    score_threashold = (double)(n)/(double)(m);
    score_threashold *= 1.1;

    /* initial choice for v */
    l = n_gcd(n, m);
    for (i = 0; i < l; i++)
    {
        fmpq_set_si(v + i, n, m);
    }

    i = 0;
    while (i < l)
    {
        if (fmpz_cmp_ui(fmpq_denref(v + i), 2) >= 0)
        {
            fmpq_farey_neighbors(left, right, v + i, fmpq_denref(v + i));

            FLINT_ASSERT(fmpz_get_ui(fmpq_denref(v + i))
                            ==  fmpz_get_ui(fmpq_denref(left))
                              + fmpz_get_ui(fmpq_denref(right)));

            if (fmpq_sgn(left) > 0 && fmpq_get_d(right) < score_threashold)
            {
                /* delete v + i, add left and right */
                FLINT_ASSERT(l < m);
                fmpq_set(v + i, right);
                fmpq_set(v + l, left);
                l++;
                continue;
            }
        }
        i++;
    }

    fmpq_clear(left);
    fmpq_clear(right);

    return l;
}

/* inputs A and B are modified */
int fmpz_mpolyl_gcd_brown_threaded_pool(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    fmpz_mpoly_t A,
    fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, j, k;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong num_threads = num_handles + 1;
    slong num_master_threads;
    slong num_images;
    int success;
    fmpz_t bound, modulus, temp, temp2;
    fmpz maxcoeff_Gs_Abars_Bbars[3];
    fmpz sumcoeff_Gs_Abars_Bbars[3];
    fmpz_t cA, cB, cG, cAbar, cBbar;
    fmpz * moduli;
    fmpz_mpoly_struct ** gptrs, ** abarptrs, ** bbarptrs;
    fmpq * qvec;
    _splitworker_arg_struct * splitargs;
    _splitbase_t splitbase;
    _njoinworker_arg_struct * joinargs;
    _joinbase_t joinbase;

    fmpz_init(bound);
    fmpz_init(modulus);
    fmpz_init(temp);
    fmpz_init(temp2);
    for (i = 0; i < 3; i++)
    {
        fmpz_init(maxcoeff_Gs_Abars_Bbars + i);
        fmpz_init(sumcoeff_Gs_Abars_Bbars + i);
    }

    /* compute contents of G, Abar, Bbar, A, B */
    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);
    fmpz_init(cAbar);
    fmpz_init(cBbar);
    _fmpz_vec_content(cA, A->coeffs, A->length);
    _fmpz_vec_content(cB, B->coeffs, B->length);
    fmpz_gcd(cG, cA, cB);
    fmpz_divexact(cAbar, cA, cG);
    fmpz_divexact(cBbar, cB, cG);

    /* remove content from inputs */
    fmpz_mpoly_scalar_divexact_fmpz(A, A, cA, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(B, B, cB, ctx);

    /* init split info */
    qvec = _fmpq_vec_init(num_threads);

    fmpz_init(splitbase->gamma);
    fmpz_gcd(splitbase->gamma, A->coeffs + 0, B->coeffs + 0);

    /*
        if compute_split is jumped back to, there could be as many as
        num_threads + 1 images that need to be joined.
    */
    gptrs = FLINT_ARRAY_ALLOC(num_threads + 1, fmpz_mpoly_struct *);
    abarptrs = FLINT_ARRAY_ALLOC(num_threads + 1, fmpz_mpoly_struct *);
    bbarptrs =FLINT_ARRAY_ALLOC(num_threads + 1, fmpz_mpoly_struct *);
    moduli = FLINT_ARRAY_ALLOC(num_threads + 1, fmpz); /* shallow copies */
    splitargs = FLINT_ARRAY_ALLOC(num_threads, _splitworker_arg_struct);
    for (i = 0; i < num_threads; i++)
    {
        fmpz_mpoly_init3(splitargs[i].G, 0, bits, ctx);
        fmpz_mpoly_init3(splitargs[i].Abar, 0, bits, ctx);
        fmpz_mpoly_init3(splitargs[i].Bbar, 0, bits, ctx);
        fmpz_init(splitargs[i].modulus);
        splitargs[i].worker_handles = (thread_pool_handle *) flint_malloc(
                                       num_threads*sizeof(thread_pool_handle));
    }

    splitbase->num_threads = num_threads;
    splitbase->A = A;
    splitbase->B = B;
    splitbase->ctx = ctx;
    splitbase->p = UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX);
    splitbase->I = I;

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&splitbase->mutex, NULL);
#endif

    /* initial bound on target modulus */
    fmpz_mpoly_height(bound, A, ctx);
    fmpz_mpoly_height(temp, B, ctx);

    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, splitbase->gamma);
    fmpz_add(bound, bound, bound);

    /* no images yet */
    fmpz_one(modulus);

compute_split:

    splitbase->gcd_is_one = 0;
    fmpz_cdiv_q(temp, bound, modulus);
    fmpz_add_ui(temp, temp, 2);

    /*
        n := fmpz_clog_ui(temp, splitbase->p) is the number of images we need.
        We have m := num_threads threads available for work.
        We can calculate images in parallel (parallel on p), but the
        calculation of each image can itself also be done in parallel.

        An integer l := num_master_threads with 1 <= l <= min(n, m) is selected
        and fractions a_i/b_i, 0 <= i < l are also selected with sum_i(a_i) = n
        and sum_i(b_i) = m. Then l jobs are spawned, each doing a_i images
        where each image uses b_i threads.
    */
    num_master_threads = _divide_master_threads(qvec,
                                fmpz_clog_ui(temp, splitbase->p), num_threads);
    FLINT_ASSERT(num_master_threads > 0);

    k = 0;
    for (i = 0; i < num_master_threads; i++)
    {
        splitargs[i].idx = i;
        splitargs[i].base = splitbase;
        FLINT_ASSERT(fmpz_fits_si(fmpq_numref(qvec + i)));
        FLINT_ASSERT(fmpz_fits_si(fmpq_denref(qvec + i)));
        splitargs[i].required_images = fmpz_get_si(fmpq_numref(qvec + i));
        splitargs[i].num_handles = fmpz_get_si(fmpq_denref(qvec + i)) - 1;
        FLINT_ASSERT(splitargs[i].required_images > 0);
        FLINT_ASSERT(splitargs[i].num_handles >= 0);

        if (i == 0)
        {
            splitargs[i].master_handle = -1;
        }
        else
        {
            splitargs[i].master_handle = handles[k++];
        }
        FLINT_ASSERT(splitargs[i].num_handles <= num_handles);
        for (j = 0; j < splitargs[i].num_handles; j++)
        {
            splitargs[i].worker_handles[j] = handles[k++];
        }
    }
    /* all handles should have been distributed */
    FLINT_ASSERT(k == num_handles);

    for (i = 1; i < num_master_threads; i++)
    {
        thread_pool_wake(global_thread_pool, splitargs[i].master_handle, 0,
                                                 _splitworker, &splitargs[i]);
    }
    _splitworker(&splitargs[0]);
    for (i = 1; i < num_master_threads; i++)
    {
        thread_pool_wait(global_thread_pool, splitargs[i].master_handle);
    }

    if (splitbase->gcd_is_one)
    {
        fmpz_mpoly_one(G, ctx);
        fmpz_mpoly_swap(Abar, A, ctx);
        fmpz_mpoly_swap(Bbar, B, ctx);
        goto successful_put_content;
    }

    /* check each thread reached its goal */
    for (i = 0; i < num_master_threads; i++)
    {
        if (splitargs[i].image_count < splitargs[i].required_images)
        {
            /* ran out of rational primes - must fail */
            success = 0;
            goto cleanup_split;
        }
    }

    /* find images to join */
    num_images = 0;
    if (!fmpz_is_one(modulus))
    {
        gptrs[num_images] = G;
        abarptrs[num_images] = Abar;
        bbarptrs[num_images] = Bbar;
        moduli[num_images] = *modulus;
        num_images++;
    }

    i = 0;
    if (num_images == 0)
    {
        gptrs[num_images] = splitargs[i].G;
        abarptrs[num_images] = splitargs[i].Abar;
        bbarptrs[num_images] = splitargs[i].Bbar;
        moduli[num_images] = *splitargs[i].modulus;
        num_images++;
        i++;
    }

    FLINT_ASSERT(num_images <= num_master_threads + 1);

    for (; i < num_master_threads; i++)
    {
        int cmp = 0;
        cmp = mpoly_monomial_cmp_nomask(gptrs[0]->exps + N*0,
                                        splitargs[i].G->exps + N*0, N);

        if (cmp < 0)
        {
            /* splitarg[i] was unlucky - ignore it */
        }
        else
        {
            if (cmp > 0)
            {
                /* splitarg[0], ..., splitarg[i - 1] were unlucky */
                num_images = 0;
            }
            gptrs[num_images] = splitargs[i].G;
            abarptrs[num_images] = splitargs[i].Abar;
            bbarptrs[num_images] = splitargs[i].Bbar;
            moduli[num_images] = *splitargs[i].modulus;
            num_images++;
        }
        FLINT_ASSERT(num_images <= num_master_threads + 1);
    }

    /* now must join ptrs[0], ..., ptrs[num_images-1] where num_images > 0 */
    FLINT_ASSERT(num_images > 0);
    fmpz_multi_CRT_init(joinbase->CRT);
    success = fmpz_multi_CRT_precompute(joinbase->CRT, moduli, num_images);
    FLINT_ASSERT(success);

    joinbase->num_images = num_images;
    joinbase->gptrs = gptrs;
    joinbase->abarptrs = abarptrs;
    joinbase->bbarptrs = bbarptrs;
    joinbase->G = G;
    joinbase->Abar = Abar;
    joinbase->Bbar = Bbar;
    joinbase->ctx = ctx;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&joinbase->mutex, NULL);
#endif

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
        fmpz_mpoly_struct * poly = (i == 0) ? gptrs[0]
                                 : (i == 1) ? abarptrs[0]
                                 :            bbarptrs[0];

        for (j = 0; j < num_threads; j++)
        {
            _joinworker_arg_struct * d = joinbase->chunks + i*num_threads + j;
            fmpz_mpoly_init3(d->poly, 0, bits, ctx);
            fmpz_init(d->maxcoeff);
            fmpz_init(d->sumcoeff);
            d->GAB = i;
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
    for (i = 0; i < 3; i++)
    {
        fmpz_zero(maxcoeff_Gs_Abars_Bbars + i);
        fmpz_zero(sumcoeff_Gs_Abars_Bbars + i);
    }

    {
        slong idxs[3] = {0, 0, 0};
        for (i = 0; i < joinbase->chunks_length; i++)
        {
            int type = joinbase->chunks[i].GAB;
            FLINT_ASSERT(0 <= type && type < 3);

            joinbase->chunks[i].final_idx = idxs[type];
            idxs[type] += joinbase->chunks[i].poly->length;

            fmpz_add(sumcoeff_Gs_Abars_Bbars + type,
                     sumcoeff_Gs_Abars_Bbars + type,
                     joinbase->chunks[i].sumcoeff);

            if (fmpz_cmp(maxcoeff_Gs_Abars_Bbars + type,
                         joinbase->chunks[i].maxcoeff) < 0)
            {
                fmpz_set(maxcoeff_Gs_Abars_Bbars + type,
                             joinbase->chunks[i].maxcoeff);
            }
        }

        fmpz_mpoly_fit_length(G, idxs[0], ctx);
        fmpz_mpoly_fit_length(Abar, idxs[1], ctx);
        fmpz_mpoly_fit_length(Bbar, idxs[2], ctx);
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

    FLINT_ASSERT(fmpz_mpoly_is_canonical(G, ctx));
    FLINT_ASSERT(fmpz_mpoly_is_canonical(Abar, ctx));
    FLINT_ASSERT(fmpz_mpoly_is_canonical(Bbar, ctx));

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&joinbase->mutex);
#endif

    /* free join data */
    fmpz_multi_CRT_clear(joinbase->CRT);
    for (i = 0; i < joinbase->chunks_length; i++)
    {
        fmpz_mpoly_clear(joinbase->chunks[i].poly, ctx);
        fmpz_clear(joinbase->chunks[i].maxcoeff);
        fmpz_clear(joinbase->chunks[i].sumcoeff);
    }

    flint_free(joinbase->chunks);
    flint_free(joinargs);

    /* update modulus - modulus could be one of the moduli[i] */
    fmpz_one(temp);
    for (i = 0; i < num_images; i++)
        fmpz_mul(temp, temp, moduli + i);
    fmpz_swap(modulus, temp);

    for (i = 1; i < 3; i++)
    {
        fmpz_mul(temp, maxcoeff_Gs_Abars_Bbars + 0,
                       sumcoeff_Gs_Abars_Bbars + i);
        fmpz_mul(temp2, sumcoeff_Gs_Abars_Bbars + 0,
                        maxcoeff_Gs_Abars_Bbars + i);
        if (fmpz_cmp(temp, temp2) > 0)
            fmpz_swap(temp, temp2);
        fmpz_mul_2exp(temp, temp, 1);
        if (fmpz_cmp(temp, modulus) >= 0)
        {
            fmpz_mul_2exp(bound, modulus, 2*FLINT_BITS);
            goto compute_split;
        }
    }

    FLINT_ASSERT(fmpz_equal(splitbase->gamma, G->coeffs + 0));

    _fmpz_vec_content(temp, G->coeffs, G->length);
    fmpz_mpoly_scalar_divexact_fmpz(G, G, temp, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(Abar, Abar, G->coeffs + 0, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(Bbar, Bbar, G->coeffs + 0, ctx);

successful_put_content:

    fmpz_mpoly_scalar_mul_fmpz(G, G, cG, ctx);
    fmpz_mpoly_scalar_mul_fmpz(Abar, Abar, cAbar, ctx);
    fmpz_mpoly_scalar_mul_fmpz(Bbar, Bbar, cBbar, ctx);

    success = 1;

cleanup_split:

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&splitbase->mutex);
#endif
    fmpz_clear(splitbase->gamma);

    for (i = 0; i < num_threads; i++)
    {
        fmpz_mpoly_clear(splitargs[i].G, ctx);
        fmpz_mpoly_clear(splitargs[i].Abar, ctx);
        fmpz_mpoly_clear(splitargs[i].Bbar, ctx);
        fmpz_clear(splitargs[i].modulus);
        flint_free(splitargs[i].worker_handles);
    }

    flint_free(gptrs);
    flint_free(abarptrs);
    flint_free(bbarptrs);
    flint_free(moduli); /* they were just shallow copies */
    flint_free(splitargs);

    _fmpq_vec_clear(qvec, num_threads);

    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);
    fmpz_clear(cAbar);
    fmpz_clear(cBbar);

    fmpz_clear(bound);
    fmpz_clear(modulus);
    fmpz_clear(temp);
    fmpz_clear(temp2);
    for (i = 0; i < 3; i++)
    {
        fmpz_clear(maxcoeff_Gs_Abars_Bbars + i);
        fmpz_clear(sumcoeff_Gs_Abars_Bbars + i);
    }

    return success;
}

