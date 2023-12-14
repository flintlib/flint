/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mpoly.h"

#ifdef fmpz_mod_ctx_get_modulus_mpz_read_only
# undef fmpz_mod_ctx_get_modulus_mpz_read_only
#endif

static inline void
fmpz_mod_ctx_get_modulus_mpz_read_only(mpz_t m, const fmpz_mod_ctx_t ctx)
{
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    if (COEFF_IS_MPZ(*p))
    {
        *m = *COEFF_TO_PTR(*p);
    }
    else
    {
        m->_mp_size = 1;
        m->_mp_alloc = 1;
        m->_mp_d = (mp_ptr) p;
    }
}

void _fmpz_mod_mpoly_mul_johnson1(
    fmpz_mod_mpoly_t A,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    const fmpz * Ccoeffs, const ulong * Cexps, slong Clen,
    ulong cmpmask,
    const fmpz_mod_ctx_t ctx)
{
    slong n = fmpz_size(fmpz_mod_ctx_modulus(ctx));
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong * hind;
    ulong exp;
    fmpz * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    slong Alen;
    mpz_t t, acc, modulus;
    mp_limb_t * Bcoeffs_packed = NULL;
    mp_limb_t * Ccoeffs_packed = NULL;
    TMP_INIT;

    TMP_START;

    mpz_init(t);
    mpz_init(acc);
    fmpz_mod_ctx_get_modulus_mpz_read_only(modulus, ctx);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    if (Blen > 8*n)
    {
        Bcoeffs_packed = FLINT_ARRAY_ALLOC(n*(Blen + Clen), mp_limb_t);
        Ccoeffs_packed = Bcoeffs_packed + n*Blen;
        for (i = 0; i < Blen; i++)
            fmpz_get_ui_array(Bcoeffs_packed + n*i, n, Bcoeffs + i);
        for (i = 0; i < Clen; i++)
            fmpz_get_ui_array(Ccoeffs_packed + n*i, n, Ccoeffs + i);
    }

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

        _fmpz_mod_mpoly_fit_length(&Acoeffs, &A->coeffs_alloc,
                                   &Aexps, &A->exps_alloc, 1, Alen + 1);
        Aexps[Alen] = exp;

        if (Bcoeffs_packed)
        {
            mp_limb_t * acc_d, * t_d;
            slong acc_len;

            FLINT_MPZ_REALLOC(acc, 2*n+1);
            FLINT_MPZ_REALLOC(t, 2*n);
            acc_d = acc->_mp_d;
            t_d = t->_mp_d;

            flint_mpn_zero(acc_d, 2*n+1);
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                    hind[x->i] |= WORD(1);
                    flint_mpn_mul_n(t_d, Bcoeffs_packed + n*x->i,
                                   Ccoeffs_packed + n*x->j, n);
                    acc_d[2*n] += mpn_add_n(acc_d, acc_d, t_d, 2*n);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            acc_len = 2*n+1;
            MPN_NORM(acc_d, acc_len);
            acc->_mp_size = acc_len;
        }
        else
        {
            mpz_set_ui(acc, 0);
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do {
                    fmpz Bi = Bcoeffs[x->i];
                    fmpz Cj = Ccoeffs[x->j];

                    *store++ = x->i;
                    *store++ = x->j;

                    hind[x->i] |= WORD(1);

                    if (COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Cj))
                    {
                        mpz_addmul(acc, COEFF_TO_PTR(Bi), COEFF_TO_PTR(Cj));
                    }
                    else if (COEFF_IS_MPZ(Bi) && !COEFF_IS_MPZ(Cj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Bi), Cj);
                    }
                    else if (!COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Cj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Cj), Bi);
                    }
                    else
                    {
                        ulong pp1, pp0;
                        umul_ppmm(pp1, pp0, Bi, Cj);
                        flint_mpz_add_uiui(acc, acc, pp1, pp0);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        mpz_tdiv_qr(t, _fmpz_promote(Acoeffs + Alen), acc, modulus);
        _fmpz_demote_val(Acoeffs + Alen);
        Alen += !fmpz_is_zero(Acoeffs + Alen);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if ((i + 1 < Blen) &&
                (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexps[x->i] + Cexps[x->j], x,
                                                 &next_loc, &heap_len, cmpmask);
            }

            /* should we go up? */
            if ((j + 1 < Clen) &&
                ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1)))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexps[x->i] + Cexps[x->j], x,
                                                 &next_loc, &heap_len, cmpmask);
            }
        }
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;
    A->length = Alen;

    mpz_clear(t);
    mpz_clear(acc);
    flint_free(Bcoeffs_packed);

    TMP_END;
}


void _fmpz_mod_mpoly_mul_johnson(
    fmpz_mod_mpoly_t A,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    const fmpz * Ccoeffs, const ulong * Cexps, slong Clen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    const fmpz_mod_ctx_t ctx)
{
    slong n = fmpz_size(fmpz_mod_ctx_modulus(ctx));
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    ulong * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    fmpz * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    slong Alen;
    mpz_t t, acc, modulus;
    mp_limb_t * Bcoeffs_packed = NULL;
    mp_limb_t * Ccoeffs_packed = NULL;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);
    FLINT_ASSERT(A->bits == bits);

    if (N == 1)
    {
        _fmpz_mod_mpoly_mul_johnson1(A, Bcoeffs, Bexps, Blen,
                                        Ccoeffs, Cexps, Clen, cmpmask[0], ctx);
        return;
    }

    TMP_START;

    mpz_init(t);
    mpz_init(acc);
    fmpz_mod_ctx_get_modulus_mpz_read_only(modulus, ctx);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));

    for (i = 0; i < Blen; i++)
    {
        exp_list[i] = exps + i*N;
        hind[i] = 1;
    }

    if (Blen > 8*n)
    {
        Bcoeffs_packed = FLINT_ARRAY_ALLOC(n*(Blen + Clen), mp_limb_t);
        Ccoeffs_packed = Bcoeffs_packed + n*Blen;
        for (i = 0; i < Blen; i++)
            fmpz_get_ui_array(Bcoeffs_packed + n*i, n, Bcoeffs + i);
        for (i = 0; i < Clen; i++)
            fmpz_get_ui_array(Ccoeffs_packed + n*i, n, Ccoeffs + i);
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

    mpoly_monomial_add_mp(heap[1].exp, Bexps + N*0, Cexps + N*0, N);

    hind[0] = 2*1 + 0;

    Alen = 0;
    while (heap_len > 1)
    {
        _fmpz_mod_mpoly_fit_length(&Acoeffs, &A->coeffs_alloc,
                                   &Aexps, &A->exps_alloc, N, Alen + 1);

        mpoly_monomial_set(Aexps + N*Alen, heap[1].exp, N);

        if (Bcoeffs_packed)
        {
            mp_limb_t * acc_d, * t_d;
            slong acc_len;

            FLINT_MPZ_REALLOC(acc, 2*n+1);
            FLINT_MPZ_REALLOC(t, 2*n);
            acc_d = acc->_mp_d;
            t_d = t->_mp_d;

            flint_mpn_zero(acc_d, 2*n+1);
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                    hind[x->i] |= WORD(1);
                    flint_mpn_mul_n(t_d, Bcoeffs_packed + n*x->i,
                                   Ccoeffs_packed + n*x->j, n);
                    acc_d[2*n] += mpn_add_n(acc_d, acc_d, t_d, 2*n);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 &&
                     mpoly_monomial_equal(heap[1].exp, Aexps + N*Alen, N));
            acc_len = 2*n+1;
            MPN_NORM(acc_d, acc_len);
            acc->_mp_size = acc_len;
        }
        else
        {
            mpz_set_ui(acc, 0);
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    fmpz Bi, Cj;

                    *store++ = x->i;
                    *store++ = x->j;

                    Bi = Bcoeffs[x->i];
                    Cj = Ccoeffs[x->j];

                    hind[x->i] |= WORD(1);

                    if (COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Cj))
                    {
                        mpz_addmul(acc, COEFF_TO_PTR(Bi), COEFF_TO_PTR(Cj));
                    }
                    else if (COEFF_IS_MPZ(Bi) && !COEFF_IS_MPZ(Cj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Bi), Cj);
                    }
                    else if (!COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Cj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Cj), Bi);
                    }
                    else
                    {
                        ulong pp1, pp0;
                        umul_ppmm(pp1, pp0, Bi, Cj);
                        flint_mpz_add_uiui(acc, acc, pp1, pp0);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 &&
                     mpoly_monomial_equal(heap[1].exp, Aexps + N*Alen, N));
        }

        mpz_tdiv_qr(t, _fmpz_promote(Acoeffs + Alen), acc, modulus);
        _fmpz_demote_val(Acoeffs + Alen);
        Alen += !fmpz_is_zero(Acoeffs + Alen);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if ((i + 1 < Blen) &&
                (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexps + N*x->i,
                                                          Cexps + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if ((j + 1 < Clen) &&
                ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1)))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexps + N*x->i,
                                                          Cexps + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;
    A->length = Alen;

    mpz_clear(t);
    mpz_clear(acc);
    flint_free(Bcoeffs_packed);

    TMP_END;
}

void _fmpz_mod_mpoly_mul_johnson_maxfields(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B, fmpz * maxBfields,
    const fmpz_mod_mpoly_t C, fmpz * maxCfields,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t Abits;
    ulong * cmpmask;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    int freeBexps = 0, freeCexps = 0;
    fmpz_mod_mpoly_struct * P, T[1];
    TMP_INIT;

    FLINT_ASSERT(B->length > 0 && C->length > 0);

    TMP_START;

    _fmpz_vec_add(maxBfields, maxBfields, maxCfields, ctx->minfo->nfields);

    Abits = 1 + _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(Abits, B->bits);
    Abits = FLINT_MAX(Abits, C->bits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);

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
        fmpz_mod_mpoly_init(T, ctx);
        P = T;
    }
    else
    {
        P = A;
    }

    fmpz_mod_mpoly_fit_length_reset_bits(P, B->length + C->length, Abits, ctx);

    if (B->length > C->length)
    {
        _fmpz_mod_mpoly_mul_johnson(P, C->coeffs, Cexps, C->length,
                  B->coeffs, Bexps, B->length, Abits, N, cmpmask, ctx->ffinfo);
    }
    else
    {
        _fmpz_mod_mpoly_mul_johnson(P, B->coeffs, Bexps, B->length,
                  C->coeffs, Cexps, C->length, Abits, N, cmpmask, ctx->ffinfo);
    }

    if (A == B || A == C)
    {
        fmpz_mod_mpoly_swap(A, T, ctx);
        fmpz_mod_mpoly_clear(T, ctx);
    }

    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;
}

void fmpz_mod_mpoly_mul_johnson(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_t C,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length < 1 || C->length < 1)
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    maxBfields = TMP_ARRAY_ALLOC(2*ctx->minfo->nfields, fmpz);
    maxCfields = maxBfields + ctx->minfo->nfields;
    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    _fmpz_mod_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);

    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;
}

