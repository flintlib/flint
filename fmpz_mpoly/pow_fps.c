/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017,2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "long_extras.h"
#include "fmpz_mod_mpoly.h"

/*
Dense algorithm (reference) for the k^th power:

f = f[0] + f[1]*x + ... + f[d]*x^d
g = f^k
  = g[0] + g[1]*x + ... + g[k*d]*x^(k*d)

where i*g[i]*f[0] = sum_{0 <= j <= min(i, d)} ((k + 1)*j - i)*f[j]*g[i - j]

g_[k*d] = f[d]^k
for i from k*d - 1 to 0
    c = 0
    for j from 1 to min(d,k*d-i)
        c += (i + j - k*(d - j))*g[i + j]*f[d - j]
    end
    g[i] = c/((k*d - i)*f[d])
end


Dense algorithm (used here) for the k^th power and (k - 1)^st power:

g comes out as f^(k-1)
a comes out as f^k

g[(k - 1)*d] = f[d]^(k-1)
a[k*d] = f[d]^k
for i from k*d - 1 to 0
    s = c = 0
    for k from max(0, d + i - d*k) to min(d - 1, i)
        s += g[i - j]*f[j]
        if i >= d
            c += (i - k*j)*g[i - j]*f[j]
        end
    end
    c = c/((d*k - i)*f[d])
    s += c*f[d]
    if i >= d
        g[i - d] = c
    end
    a[i] = s
end

For the sparse version, the only thing to note is that the i >= d condition
is implemented as _divisibility_ of monomials and not merely a comparison.
*/

static slong _fmpz_mpoly_pow_fps1(
    fmpz_mpoly_t A,
    const fmpz * Fcoeffs, const ulong * Fexps, slong Flen,
    ulong k,
    ulong cmpmask,
    ulong ofmask)
{
    slong i, j;
    slong next_loc, Qlen = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    ulong exp;
    slong * hind;
    slong * Q;
    mpoly_heap_t * x;
    fmpz * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    slong Alen;
    fmpz * Gcoeffs = NULL;
    ulong * Gexps;
    slong Galloc, Glen, Gdemote = 0;
    fmpz_t t1, temp1;
    fmpz * S, * C;
    int divides;
    TMP_INIT;

    TMP_START;

    next_loc = Flen + 4;   /* something bigger than heap can ever be */
    heap = TMP_ARRAY_ALLOC(Flen + 1, mpoly_heap1_s);
    chain = TMP_ARRAY_ALLOC(Flen, mpoly_heap_t);
    Q = TMP_ARRAY_ALLOC(3*Flen, slong);
    hind = Q + 2*Flen; /* flagged heap indices */

    for (i = 0; i < Flen; i++)
        hind[i] = 1;

    fmpz_init(t1);
    fmpz_init(temp1);

    _fmpz_mpoly_fit_length(&Acoeffs, &Aexps, &A->alloc, 2, 1);

    Galloc = k*(Flen - 1) + 2;
    Gexps = FLINT_ARRAY_ALLOC(Galloc, ulong);
    Gcoeffs = (fmpz *) flint_calloc(Galloc, sizeof(fmpz));

    Gexps[0] = Fexps[0]*(k - 1);
    Aexps[0] = Fexps[0]*k;
    fmpz_pow_ui(Gcoeffs + 0, Fcoeffs + 0, k - 1);
    fmpz_mul(Acoeffs + 0, Gcoeffs + 0, Fcoeffs + 0);
    Glen = 1;
    Alen = 1;

    x = chain + 1;
    x->i = 1;
    x->j = 0;
    x->next = NULL;
    hind[1] = 2*1 + 0;
    HEAP_ASSIGN(heap[1], Fexps[1] + Gexps[0], x);

    while (heap_len > 1)
    {
        exp = heap[1].exp;
        Aexps[Alen] = exp;
        S = Acoeffs + Alen;
        C = Gcoeffs + Glen;

        fmpz_zero(C);
        fmpz_zero(S);
        Qlen = 0;

        Gexps[Glen] = exp - Fexps[0];

        divides = !mpoly_monomial_overflows1(Gexps[Glen], ofmask);

        do {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
            do {
                fmpz fi;
                slong dd, ddd;

                Q[Qlen++] = i = x->i;
                Q[Qlen++] = j = x->j;
                hind[i] |= 1;

                fi = Fcoeffs[i];

                FLINT_ASSERT(j >= Gdemote);

                dd = (k - 1)*Fexps[i] - Gexps[j];

                if (divides)
                {
                    if (!COEFF_IS_MPZ(fi) && !z_mul_checked(&ddd, dd, fi))
                    {
                        fmpz_addmul_si(S, Gcoeffs + j, fi);
                        fmpz_addmul_si(C, Gcoeffs + j, ddd);
                    }
                    else
                    {
                        fmpz_mul(t1, Fcoeffs + i, Gcoeffs + j);
                        fmpz_add(S, S, t1);
                        fmpz_addmul_si(C, t1, dd);
                    }
                }
                else
                {
                    fmpz_addmul(S, Fcoeffs + x->i, Gcoeffs + x->j);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        FLINT_ASSERT(Qlen <= 2*Flen);
      
        while (Qlen > 0)
        {
            /* take node from store */
            j = Q[--Qlen];
            i = Q[--Qlen];

            /* should we go right? */
            if (i + 1 < Flen && hind[i + 1] == 2*j + 1)
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*j + 2;

                FLINT_ASSERT(j >= Gdemote);
                _mpoly_heap_insert1(heap, Fexps[i + 1] + Gexps[j], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            /* should we go up? */
            if (j + 1 < Glen && hind[i] < 2*j + 4)
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*j + 4;

                FLINT_ASSERT(j + 1 >= Gdemote);
                _mpoly_heap_insert1(heap, Fexps[i] + Gexps[j + 1], x,
                                                &next_loc, &heap_len, cmpmask);
            }
        }

        /* cleanup unused Gcoeffs */
        j = hind[Flen - 1]/2 - 1;
        for ( ; Gdemote < j; Gdemote++)
            fmpz_clear(Gcoeffs + Gdemote);

        if (!fmpz_is_zero(C))
        {
            if (fmpz_is_one(Fcoeffs + 0))
            {
                fmpz_divexact_si(C, C, exp - k*Fexps[0]);
                fmpz_add(S, S, C);
            }
            else
            {
                fmpz_divexact_si(temp1, C, exp - k*Fexps[0]);
                fmpz_add(S, S, temp1);
                fmpz_divexact(C, temp1, Fcoeffs + 0);
            }

            if ((hind[1] & 1) != 0)
            {
                x = chain + 1;
                x->i = 1;
                x->j = Glen;
                x->next = NULL;

                hind[x->i] = 2*(Glen + 1) + 0;

                FLINT_ASSERT(Glen >= Gdemote);
                _mpoly_heap_insert1(heap, Fexps[1] + Gexps[Glen], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            Glen++;

            if (Glen >= Galloc)
            {
                Gexps = FLINT_ARRAY_REALLOC(Gexps, 2*Galloc, ulong);
                Gcoeffs = FLINT_ARRAY_REALLOC(Gcoeffs, 2*Galloc, fmpz);
                flint_mpn_zero(Gcoeffs + Galloc, Galloc);
                Galloc *= 2;
            }
        }

        Alen += !fmpz_is_zero(S);
        _fmpz_mpoly_fit_length(&Acoeffs, &Aexps, &A->alloc, Alen + 1, 1);
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;

    fmpz_clear(t1);
    fmpz_clear(temp1);

    for ( ; Gdemote < Galloc; Gdemote++)
        fmpz_clear(Gcoeffs + Gdemote);
    flint_free(Gcoeffs);
    flint_free(Gexps);

    TMP_END;

    return Alen;
}


static slong _fmpz_mpoly_pow_fps(
    fmpz_mpoly_t A,
    const fmpz * Fcoeffs, const ulong * Fexps, slong Flen,
    ulong k,
    slong N,
    const ulong * cmpmask)
{
    flint_bitcnt_t bits = A->bits;
    ulong ofmask = (bits > FLINT_BITS) ? 0 : mpoly_overflow_mask_sp(bits);
    slong i, j, exp_next;
    slong next_loc, heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    slong Qlen = 0;
    slong * hind;
    mpoly_heap_t * x;
    fmpz * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    slong Alen;
    fmpz * Gcoeffs;
    slong Galloc, Glen, Gdemote = 0;
    ulong * Gexps;
    ulong * fik, * exps;
    ulong ** exp_list;
    ulong * temp2;
    fmpz_t t1, t2, C;
    int divides;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_pow_fps1(A, Fcoeffs, Fexps, Flen, k, cmpmask[0], ofmask);

    TMP_START;

    next_loc = Flen + 4;   /* something bigger than heap can ever be */
    heap = TMP_ARRAY_ALLOC(Flen + 1, mpoly_heap_s);
    chain = TMP_ARRAY_ALLOC(Flen, mpoly_heap_t);
    exp_list = TMP_ARRAY_ALLOC((Flen + 1), ulong *);
    exps = TMP_ARRAY_ALLOC(N*(Flen + 1) + N*Flen + N, ulong);
    fik = exps + N*(Flen + 1);
    temp2 = fik + N*Flen;
    Q = TMP_ARRAY_ALLOC(3*Flen, slong);
    hind = Q + 2*Flen; /* flagged heap indices */

    for (i = 0; i < Flen; i++)
        mpoly_monomial_mul_ui_mp(fik + i*N, Fexps + i*N, N, k - 1);

    for (i = 0; i < Flen; i++)
        hind[i] = 1;

    exp_next = 0;
    for (i = 0; i < Flen + 1; i++)
        exp_list[i] = exps + N*i;

    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(C);

    _fmpz_mpoly_fit_length(&Acoeffs, &Aexps, &A->alloc, 2, N);

    Galloc = (k - 1)*Flen + 2;
    Gexps = (ulong *) flint_malloc(Galloc*sizeof(ulong)*N);
    Gcoeffs = (fmpz *) flint_calloc(Galloc, sizeof(fmpz));

    mpoly_monomial_mul_ui_mp(Gexps + 0, Fexps + 0, N, k - 1);
    mpoly_monomial_mul_ui_mp(Aexps + 0, Fexps + 0, N, k);
    fmpz_pow_ui(Gcoeffs + 0, Fcoeffs + 0, k - 1);
    fmpz_mul(Acoeffs + 0, Gcoeffs + 0, Fcoeffs + 0);
    Glen = 1;
    Alen = 1;

    hind[1] = 2*1 + 0;
    x = chain + 1;
    x->i = 1;
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_add_mp(heap[1].exp, Fexps + N, Gexps + 0, N);

    while (heap_len > 1)
    {
        fmpz * const S = Acoeffs + Alen;

        mpoly_monomial_set(Aexps + N*Alen, heap[1].exp, N);
        mpoly_monomial_sub_mp(Gexps + N*Glen, Aexps + N*Alen, Fexps + 0, N);

        if (bits > FLINT_BITS)
            divides = !mpoly_monomial_overflows_mp(Gexps + N*Glen, N, bits);
        else
            divides = !mpoly_monomial_overflows(Gexps + N*Glen, N, ofmask);

        fmpz_zero(C);
        fmpz_zero(S);
        Qlen = 0;

        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, Aexps + N*Alen, N))
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            do {
                Q[Qlen++] = i = x->i;
                Q[Qlen++] = j = x->j;
                hind[i] |= 1;

                FLINT_ASSERT(j >= Gdemote);

                fmpz_mul(t1, Fcoeffs + i, Gcoeffs + j);
                fmpz_add(S, S, t1);
                if (divides)
                {
                    mpn_sub_n(temp2, fik + N*i, Gexps + N*j, N);
                    fmpz_set_signed_ui_array(t2, temp2, N);
                    fmpz_addmul(C, t1, t2);
                }
            } while ((x = x->next) != NULL);
        }

        FLINT_ASSERT(Qlen <= 2*Flen);

        while (Qlen > 0)
        {
            /* take node from store */
            j = Q[--Qlen];
            i = Q[--Qlen];

            /* should we go right? */
            if (i + 1 < Flen && hind[i + 1] == 2*j + 1)
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*j + 2;

                FLINT_ASSERT(j >= Gdemote);
                mpoly_monomial_add_mp(exp_list[exp_next], Fexps + N*(i + 1),
                                                          Gexps + N*j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if (j + 1 < Glen && hind[i] < 2*j + 4)
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*j + 4;

                FLINT_ASSERT(j + 1 >= Gdemote);
                mpoly_monomial_add_mp(exp_list[exp_next], Fexps + N*i,
                                                          Gexps + N*(j + 1), N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }

        /* cleanup unused Gcoeffs */
        j = hind[Flen - 1]/2 - 1;
        for ( ; Gdemote < j; Gdemote++)
            fmpz_clear(Gcoeffs + Gdemote);

        if (!fmpz_is_zero(C))
        {
            mpoly_monomial_mul_ui_mp(temp2, Fexps + 0, N, k);
            mpn_sub_n(temp2, Aexps + N*Alen, temp2, N);
            fmpz_set_signed_ui_array(t2, temp2, N);

            if (fmpz_is_one(Fcoeffs + 0))
            {
                fmpz_divexact(Gcoeffs + Glen, C, t2);
                fmpz_add(S, S, Gcoeffs + Glen);
            }
            else
            {
                fmpz_divexact(t1, C, t2);
                fmpz_add(S, S, t1);
                fmpz_divexact(Gcoeffs + Glen, t1, Fcoeffs + 0);
            }

            if ((hind[1] & 1) != 0)
            {
                x = chain + 1;
                x->i = 1;
                x->j = Glen;
                x->next = NULL;

                hind[x->i] = 2*(Glen + 1) + 0;

                FLINT_ASSERT(Glen >= Gdemote);
                mpoly_monomial_add_mp(exp_list[exp_next], Fexps + N,
                                                          Gexps + N*Glen, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            Glen++;

            if (Glen >= Galloc)
            {
                Gexps = FLINT_ARRAY_REALLOC(Gexps, 2*N*Galloc, ulong);
                Gcoeffs = FLINT_ARRAY_REALLOC(Gcoeffs, 2*Galloc, fmpz);
                flint_mpn_zero(Gcoeffs + Galloc, Galloc);
                Galloc *= 2;
            }
        }

        Alen += !fmpz_is_zero(Acoeffs + Alen);
        _fmpz_mpoly_fit_length(&Acoeffs, &Aexps, &A->alloc, Alen + 1, N);
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;
   
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(C);

    for ( ; Gdemote < Galloc; Gdemote++)
        fmpz_clear(Gcoeffs + Gdemote);
    flint_free(Gcoeffs);
    flint_free(Gexps);

    TMP_END;

    return Alen;
}

void fmpz_mpoly_pow_fps(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           ulong k, const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, len;
    fmpz * maxBfields;
    flint_bitcnt_t Abits;
    ulong * cmpmask;
    ulong * Bexps;
    int freeBexps;
    TMP_INIT;

    FLINT_ASSERT(k >= 2);
    FLINT_ASSERT(B->length > 0);

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    _fmpz_vec_scalar_mul_ui(maxBfields, maxBfields, ctx->minfo->nfields, k);

    Abits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(MPOLY_MIN_BITS, Abits + 1);
    Abits = FLINT_MAX(Abits, B->bits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);
    N = mpoly_words_per_exp(Abits, ctx->minfo);

    if (B->length == 1)
    {
        /* powering a monomial */
        fmpz_mpoly_fit_length_reset_bits(A, 1, Abits, ctx);

        if (B->bits == Abits && B != A)
            mpoly_monomial_mul_ui_mp(A->exps, B->exps, N, k);
        else
            mpoly_pack_vec_fmpz(A->exps, maxBfields, Abits,
                                                       ctx->minfo->nfields, 1);
        fmpz_pow_ui(A->coeffs + 0, B->coeffs + 0, k);
        len = 1;
        goto cleanup;
    }

    freeBexps = 0;
    Bexps = B->exps;
    if (Abits > B->bits)
    {
       freeBexps = 1;
       Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
       mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    if (A == B)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, k*(B->length - 1) + 1, Abits, ctx);
        len = _fmpz_mpoly_pow_fps(T, B->coeffs, Bexps, B->length, k, N, cmpmask);
        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, k*(B->length - 1) + 1, Abits, ctx);
        len = _fmpz_mpoly_pow_fps(A, B->coeffs, Bexps, B->length, k, N, cmpmask);
    }

    if (freeBexps)
        flint_free(Bexps);

cleanup:

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    _fmpz_mpoly_set_length(A, len, ctx);

    TMP_END;
}
