/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/* solve z^2+z=c */
static int _quadratic_root_const(
    mp_limb_t * z,
    const mp_limb_t * c,
    const fq_nmod_ctx_t fqctx)
{
    slong i, d = fq_nmod_ctx_degree(fqctx);
    mp_limb_t * t, * p, * u, * cp, * ut, * up, * ct;
    int success;
    TMP_INIT;
#if FLINT_WANT_ASSERT
    mp_limb_t * c_org = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    _n_fq_set(c_org, c, d);
#endif

    TMP_START;
    i = FLINT_MAX(N_FQ_REDUCE_ITCH, N_FQ_MUL_INV_ITCH);
    t = (mp_limb_t *) TMP_ALLOC((i + 7)*d*sizeof(mp_limb_t));
    p = t + d*i;
    u = p + d*2;
    ut = u + d;
    up = ut + d;
    ct = up + d;
    cp = ct + d;

    /* find u with trace 1 */
    for (i = 0; i < d; i++)
    {
        p[2*i] = 0;
        p[2*i + 1] = nmod_poly_get_coeff_ui(fqctx->modulus, 2*i + 1);
    }
    _n_fq_reduce2(ut, p, fqctx, t);
    FLINT_ASSERT(!_n_fq_is_zero(ut, d));
    _n_fq_inv(u, ut, fqctx, t);

    _n_fq_set(cp, c, d);
    _n_fq_zero(ut, d);
    _n_fq_set(up, u, d);
    _n_fq_set(ct, c, d);
    _n_fq_zero(z, d);
    for (i = 0; i < d - 1; i++)
    {
        _n_fq_add(ut, ut, up, d, fqctx->mod);
        _n_fq_mul(cp, cp, cp, fqctx, t);
        _n_fq_mul(p, cp, ut, fqctx, t);
        _n_fq_add(z, z, p, d, fqctx->mod);
        _n_fq_add(ct, ct, cp, d, fqctx->mod);
        _n_fq_mul(up, up, up, fqctx, t);
    }

    success = _n_fq_is_zero(ct, d);
#if FLINT_WANT_ASSERT
    if (success)
    {
        _n_fq_add(ut, ut, up, d, fqctx->mod);
        FLINT_ASSERT(_n_fq_is_one(ut, d));
        _n_fq_mul(p, z, z, fqctx, t);
        _n_fq_add(p, p, z, d, fqctx->mod);
        FLINT_ASSERT(_n_fq_equal(p, c_org, d));
    }
    flint_free(c_org);
#endif

    TMP_END;

    return success;
}

/*
    solve x^2+A*x=B
    if x = Q + x' for the candidate solution terms Q,
        x'^2+Ax'=B-Q^2-A*Q

      B has heap idx (-1, j)
    Q^2 has heap idx (-2, j)
    A*Q has heap idx (i, j) for 0 <= i < Alen
*/
static int _fq_nmod_mpoly_quadratic_root_heap(
    fq_nmod_mpoly_t Q,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    slong bits,
    slong N,
    const ulong * cmpmask,
    const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    slong i, j, Qlen, Qs, As;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    mp_limb_t * t, * c, * lcAinv;
    int mcmp;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    t = (mp_limb_t *) TMP_ALLOC(8*d*sizeof(mp_limb_t));
    c = t + 6*d;
    lcAinv = c + d;
    _n_fq_inv(lcAinv, Acoeffs + d*0, fqctx, t);

    /* alloc array of heap nodes which can be chained together */
    next_loc = Alen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Alen + 3)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC((Alen + 2)*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*(Alen + 2)*sizeof(slong));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC((Alen + 2)*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC((Alen + 2)*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < Alen + 2; i++)
        exp_list[i] = exps + i*N;

    mask = (bits <= FLINT_BITS) ? mpoly_overflow_mask_sp(bits) : 0;

    Qs = 1;
    As = Alen;
    mcmp = 1;
    Qlen = 0;

    /* insert (-1, 0, Bexps[0]) into heap */
    x = chain + Alen + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, Bexps + N*0, N);

    while (heap_len > 1)
    {
        FLINT_ASSERT(heap_len - 1 <= Alen + 2);

        _fq_nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc, d,
                                  &Qexps, &Q->exps_alloc, N, Qlen + 1);

        /* exp can overflow, but divisibility & halving check their answers */
        mpoly_monomial_set(exp, heap[1].exp, N);

        _nmod_vec_zero(t, 6*d);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -UWORD(1))
                {
                    FLINT_ASSERT(x->j < Blen);
                    _nmod_vec_add(t, t, Bcoeffs + d*x->j, d, fqctx->mod);
                }
                else
                {
                    const mp_limb_t * s = (x->i == -UWORD(2)) ?
                                           Qcoeffs + d*x->j : Acoeffs + d*x->i;
                    FLINT_ASSERT(x->j < Qlen);
                    FLINT_ASSERT(x->i == -UWORD(2) || x->i < Alen);
                    _n_fq_madd2(t, s, Qcoeffs + d*x->j, fqctx, t + 2*d);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        _n_fq_reduce2(Qcoeffs + d*Qlen, t, fqctx, t + 2*d);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -UWORD(1))
            {
                if (j + 1 < Blen)
                {
                    x = chain + Alen;
                    x->i = i;
                    x->j = j + 1;

                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], Bexps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);

                    FLINT_ASSERT(exp_next <= Alen + 2);
                }
            }
            else if (i == -UWORD(2))
            {
                if (j + 1 < Qlen)
                {
                    x = chain + Alen + 1;
                    x->i = i;
                    x->j = j + 1;

                    x->next = NULL;
                    mpoly_monomial_add_mp(exp_list[exp_next], Qexps + N*x->j,
                                                            Qexps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);

                    FLINT_ASSERT(exp_next <= Alen + 2);
                    FLINT_ASSERT(heap_len - 1 <= Alen + 2);
                }
                else
                {
                    FLINT_ASSERT(j + 1 == Qlen);
                    FLINT_ASSERT(Qs == 0);
                    Qs = 1;
                }
            }
            else
            {
                FLINT_ASSERT(0 <= i && i < Alen);
                if (j + 1 < Qlen)
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;

                    x->next = NULL;
                    mpoly_monomial_add_mp(exp_list[exp_next], Aexps + N*x->i,
                                                            Qexps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);

                    FLINT_ASSERT(exp_next <= Alen + 2);
                    FLINT_ASSERT(heap_len - 1 <= Alen + 2);
                }
                else
                {
                    FLINT_ASSERT(j + 1 == Qlen);
                    As++;
                    FLINT_ASSERT(As <= Alen);
                }
            }
        }

        if (_n_fq_is_zero(Qcoeffs + d*Qlen, d))
            continue;

        /*
            mcmp > 0: The last written Qexp is > lm(A)
            mcmp = 0:                          = lm(A)
            mcmp < 0:                          < lm(A)

            must find an m such that m^2 + lt(A)*m = acc*exp
        */

        if (mcmp <= 0)
            goto try_less;

        if (bits <= FLINT_BITS ?
                !mpoly_monomial_halves(Qexps + N*Qlen, exp, N, mask) :
                !mpoly_monomial_halves_mp(Qexps + N*Qlen, exp, N, bits))
        {
            goto try_less;
        }

        if (mpoly_monomial_gt(Qexps + N*Qlen, Aexps + N*0, N, cmpmask))
        {
            for (j = 1; j < d; j++)
                _n_fq_mul(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen,
                                            Qcoeffs + d*Qlen, fqctx, t);
            goto mfound;
        }

        /* try equal */

        if (!mpoly_monomial_equal(Qexps + Qlen*N, Aexps + N*0, N))
            goto try_less;

        if (d < 2)
            goto try_less;

        _n_fq_mul(c, Qcoeffs + d*Qlen, lcAinv, fqctx, t);
        _n_fq_mul(c, c, lcAinv, fqctx, t);
        if (_quadratic_root_const(c, c, fqctx))
        {
            _n_fq_mul(Qcoeffs + d*Qlen, c, Acoeffs + d*0, fqctx, t);
            mcmp = 0;
            goto mfound;
        }

    try_less:

        mcmp = -1;

        if (bits <= FLINT_BITS ?
            !mpoly_monomial_divides(Qexps + Qlen*N, exp, Aexps + N*0, N, mask) :
            !mpoly_monomial_divides_mp(Qexps + Qlen*N, exp, Aexps + N*0, N, bits))
        {
            goto no_solution;
        }

        if (!mpoly_monomial_lt(Qexps + N*Qlen, Aexps + N*0, N, cmpmask))
            goto no_solution;

        _n_fq_mul(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, lcAinv, fqctx, t);

    mfound:

        /*
            verify heap consistency
            (i >= 0, j) should be in the heap iff i >= As
            (-2, j) should be in the heap iff Qs = 0
        */
        FLINT_ASSERT(Qs == 0 || Qs == 1);
        FLINT_ASSERT(As <= Alen);
    #if FLINT_WANT_ASSERT
        {
            slong Asleft = Alen, Qsleft = 1;
            for (i = 1; i < heap_len; i++)
            {
                mpoly_heap_t * x = (mpoly_heap_t *) heap[i].next;
                do {
                    if (x->i == -UWORD(1))
                    {
                        continue;
                    }
                    else if (x->i == -UWORD(2))
                    {
                        Qsleft--;
                    }
                    else
                    {
                        FLINT_ASSERT(x->i >= As);
                        Asleft--;
                    }
                } while ((x = x->next) != NULL);
            }
            FLINT_ASSERT(Asleft == As);
            FLINT_ASSERT(Qsleft == Qs);
        }
    #endif

        FLINT_ASSERT(mcmp < 0 || Qs == 1);
        if ((mcmp >= 0) < Qs)
        {
            /* the new Q^2 term did not cancel exp */
            x = chain + Alen + 1;
            x->i = -UWORD(2);
            x->j = Qlen;
            x->next = NULL;
            mpoly_monomial_add_mp(exp_list[exp_next], Qexps + N*x->j,
                                                      Qexps + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                     &next_loc, &heap_len, N, cmpmask);

            FLINT_ASSERT(exp_next <= Alen + 2);
            FLINT_ASSERT(heap_len - 1 <= Alen + 2);
        }
        Qs = FLINT_MIN(Qs, (mcmp >= 0));

        for (i = (mcmp <= 0); i < As; i++)
        {
            /* the new Q*A[i] term did not not cancel exp */
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            mpoly_monomial_add_mp(exp_list[exp_next], Aexps + N*x->i,
                                                     Qexps + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);

            FLINT_ASSERT(exp_next <= Alen + 2);
            FLINT_ASSERT(heap_len - 1 <= Alen + 2);
        }
        As = FLINT_MIN(As, (mcmp <= 0));

        FLINT_ASSERT(!_n_fq_is_zero(Qcoeffs + d*Qlen, d));

        Qlen++;
    }

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    TMP_END;

    return 1;

no_solution:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    TMP_END;

    return 0;
}


int fq_nmod_mpoly_quadratic_root(
    fq_nmod_mpoly_t Q,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t Qbits;
    ulong * cmpmask;
    ulong * Aexps = A->exps, * Bexps = B->exps;
    int success, freeAexps = 0, freeBexps = 0;
    TMP_INIT;

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        return fq_nmod_mpoly_sqrt(Q, B, ctx);
    }

    if (ctx->fqctx->mod.n != 2)
    {
        mp_limb_t mhalf = (ctx->fqctx->mod.n - 1)/2;
        fq_nmod_mpoly_t t1, t2;
        fq_nmod_t c;

        fq_nmod_mpoly_init(t1, ctx);
        fq_nmod_mpoly_init(t2, ctx);
        fq_nmod_init(c, ctx->fqctx);

        fq_nmod_mpoly_mul(t1, A, A, ctx);
        fq_nmod_set_ui(c, nmod_mul(mhalf, mhalf, ctx->fqctx->mod), ctx->fqctx);
        fq_nmod_mpoly_scalar_addmul_fq_nmod(t2, B, t1, c, ctx);
        success = fq_nmod_mpoly_sqrt(t1, t2, ctx);
        if (success)
        {
            fq_nmod_set_ui(c, mhalf, ctx->fqctx);
            fq_nmod_mpoly_scalar_addmul_fq_nmod(Q, t1, A, c, ctx);
        }

        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_clear(t1, ctx);
        fq_nmod_mpoly_clear(t2, ctx);

        return success;
    }

    TMP_START;

    Qbits = FLINT_MAX(A->bits, B->bits);
    Qbits = mpoly_fix_bits(Qbits, ctx->minfo);
    N = mpoly_words_per_exp(Qbits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Qbits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    if (Qbits > A->bits)
    {
        freeAexps = 1;
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, Qbits, A->exps, A->bits,
                                                        A->length, ctx->minfo);
    }

    if (Qbits > B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Qbits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    if (Q == A || Q == B)
    {
        fq_nmod_mpoly_t T;
        fq_nmod_mpoly_init3(T, B->length/A->length + 1, Qbits, ctx);
        success = _fq_nmod_mpoly_quadratic_root_heap(T,
                                                A->coeffs, Aexps, A->length,
                                                B->coeffs, Bexps, B->length,
                                                Qbits, N, cmpmask, ctx->fqctx);
        fq_nmod_mpoly_swap(T, Q, ctx);
        fq_nmod_mpoly_clear(T, ctx);
    }
    else
    {
        fq_nmod_mpoly_fit_length_reset_bits(Q, B->length/A->length + 1, Qbits, ctx);
        success = _fq_nmod_mpoly_quadratic_root_heap(Q,
                                                A->coeffs, Aexps, A->length,
                                                B->coeffs, Bexps, B->length,
                                                Qbits, N, cmpmask, ctx->fqctx);
    }

    if (freeAexps)
        flint_free(Aexps);

    if (freeBexps)
        flint_free(Bexps);

    TMP_END;

    return success;
}

