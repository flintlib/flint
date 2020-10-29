/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void _fq_nmod_mpoly_mul_johnson1(
    fq_nmod_mpoly_t A,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    const mp_limb_t * Ccoeffs, const ulong * Cexps, slong Clen,
    ulong maskhi,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong * hind;
    ulong exp;
    mp_limb_t * t;
    int lazy_size = _n_fq_dot_lazy_size(Blen, ctx);
    mp_limb_t * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    slong Acoeffs_alloc = A->coeffs_alloc;
    slong Aexps_alloc = A->exps_alloc;
    slong Alen;
    TMP_INIT;

    TMP_START;

    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    t = (mp_limb_t *) TMP_ALLOC(6*d*sizeof(mp_limb_t));

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    /* put (0, 0, exp2[0] + exp3[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    HEAP_ASSIGN(heap[1], Bexps[0] + Cexps[0], x);
    hind[0] = 2*1 + 0;

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fq_nmod_mpoly_fit_length(&Acoeffs, &Acoeffs_alloc, d,
                                  &Aexps, &Aexps_alloc, 1, Alen + 1);

        Aexps[Alen] = exp;

        _nmod_vec_zero(t, 6*d);

        switch (lazy_size)
        {
        case 1:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy1(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            _n_fq_reduce2_lazy1(t, d, ctx->mod);
            break;

        case 2:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy2(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            _n_fq_reduce2_lazy2(t, d, ctx->mod);
            break;

        case 3:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy3(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            _n_fq_reduce2_lazy3(t, d, ctx->mod);
            break;

        default:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2(t, Bcoeffs + d*x->i,
                                   Ccoeffs + d*x->j, ctx, t + 2*d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            break;
        }

        _n_fq_reduce2(Acoeffs + d*Alen, t, ctx, t + 2*d);
        Alen += !_n_fq_is_zero(Acoeffs + d*Alen, d);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if ((i + 1 < Blen) &&
                (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexps[x->i] + Cexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }

            /* should we go up? */
            if ((j + 1 < Clen) &&
                ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1))
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexps[x->i] + Cexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;
    A->coeffs_alloc = Acoeffs_alloc;
    A->exps_alloc = Aexps_alloc;
    A->length = Alen;

    TMP_END;
}


void _fq_nmod_mpoly_mul_johnson(
    fq_nmod_mpoly_t A,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    const mp_limb_t * Ccoeffs, const ulong * Cexps, slong Clen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    mp_limb_t * t;
    int lazy_size = _n_fq_dot_lazy_size(Blen, ctx);
    mp_limb_t * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    slong Alen;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);
    FLINT_ASSERT(A->bits == bits);

    if (N == 1)
    {
        _fq_nmod_mpoly_mul_johnson1(A, Bcoeffs, Bexps, Blen,
                                       Ccoeffs, Cexps, Clen, cmpmask[0], ctx);
        return;
    }

    TMP_START;

    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    t = (mp_limb_t *) TMP_ALLOC(6*d*sizeof(mp_limb_t));

    for (i = 0; i < Blen; i++)
    {
        exp_list[i] = exps + i*N;
        hind[i] = 1;
    }

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put (0, 0, exp2[0] + exp3[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, Bexps + N*0, Cexps + N*0, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, Bexps + N*0, Cexps + N*0, N);

    hind[0] = 2*1 + 0;

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fq_nmod_mpoly_fit_length(&Acoeffs, &A->coeffs_alloc, d,
                                  &Aexps, &A->exps_alloc, N, Alen + 1);

        mpoly_monomial_set(Aexps + N*Alen, exp, N);

        _nmod_vec_zero(t, 6*d);

        switch (lazy_size)
        {
        case 1:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy1(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            _n_fq_reduce2_lazy1(t, d, ctx->mod);
            break;

        case 2:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy2(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            _n_fq_reduce2_lazy2(t, d, ctx->mod);
            break;

        case 3:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2_lazy3(t, Bcoeffs + d*x->i, Ccoeffs + d*x->j, d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            _n_fq_reduce2_lazy3(t, d, ctx->mod);
            break;

        default:
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    _n_fq_madd2(t, Bcoeffs + d*x->i,
                                   Ccoeffs + d*x->j, ctx, t + 2*d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            break;
        }

        _n_fq_reduce2(Acoeffs + d*Alen, t, ctx, t + 2*d);
        Alen += !_n_fq_is_zero(Acoeffs + d*Alen, d);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if ((i + 1 < Blen) &&
                (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexps + x->i*N,
                                                           Cexps + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexps + x->i*N,
                                                              Cexps + x->j*N, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if ((j + 1 < Clen) &&
                ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1))
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexps + x->i*N,
                                                           Cexps + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexps + x->i*N,
                                                              Cexps + x->j*N, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;
    A->length = Alen;

    TMP_END;
}

void fq_nmod_mpoly_mul_johnson(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_t C,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t Abits;
    fmpz * Bmaxfields, * Cmaxfields;
    ulong * cmpmask;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    int freeBexps = 0, freeCexps = 0;
    fq_nmod_mpoly_struct * P, T[1];
    TMP_INIT;

    if (B->length < 1 || C->length < 1)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    Bmaxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    Cmaxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(Bmaxfields + i);
        fmpz_init(Cmaxfields + i);
    }
    mpoly_max_fields_fmpz(Bmaxfields, Bexps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(Cmaxfields, Cexps, C->length, C->bits, ctx->minfo);
    _fmpz_vec_add(Bmaxfields, Bmaxfields, Cmaxfields, ctx->minfo->nfields);

    Abits = 1 + _fmpz_vec_max_bits(Bmaxfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(Abits, B->bits);
    Abits = FLINT_MAX(Abits, C->bits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(Bmaxfields + i);
        fmpz_clear(Cmaxfields + i);
    }

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    if (Abits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, B->length, ctx->minfo);
    }

    if (Abits != C->bits)
    {
        freeCexps = 1;
        Cexps = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits, C->length, ctx->minfo);
    }

    if (A == B || A == C)
    {
        fq_nmod_mpoly_init(T, ctx);
        P = T;
    }
    else
    {
        P = A;
    }

    fq_nmod_mpoly_fit_length_reset_bits(P, B->length + C->length, Abits, ctx);

    if (B->length > C->length)
    {
        _fq_nmod_mpoly_mul_johnson(P, C->coeffs, Cexps, C->length,
                   B->coeffs, Bexps, B->length, Abits, N, cmpmask, ctx->fqctx);
    }
    else
    {
        _fq_nmod_mpoly_mul_johnson(P, B->coeffs, Bexps, B->length,
                   C->coeffs, Cexps, C->length, Abits, N, cmpmask, ctx->fqctx);
    }

    if (A == B || A == C)
    {
        fq_nmod_mpoly_swap(A, T, ctx);
        fq_nmod_mpoly_clear(T, ctx);
    }

    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;
}
