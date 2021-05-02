/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017 Daniel Schultz

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
Dense algorithm for the k^th power

f = f_0 + f_1*x + .... + f_d*x^d
g = f^k
  = g_0 + g_1*x +      + g_{k*d}*x^(k*d)

where i*g_i*f_0 = sum_{0 <= j <= min(i,d)} ((k+1)*j-i)*f_j*g_{i-j}

g_{k*d} = f_d^k
for i from k*d-1 to 0
    c = 0
    for 1 <= j <= min(d,k*d-i)
        c += (i+j-(d-j)*k)*g_{i+j}*f_{d-j} [*x^(i+d)]
    end
    g_i = c/((k*d-i)*f_d [*x^d])
*/


/*
Sparse algorithm for k^th power:

output g = (f1 + f2 + ... + ft)^k:

g = g1^k
H = {2,1,f2*g1}
while #H > 0 && H1 >= f1
    M = H1
    C = 0
    Q = {};
    while #H > 0 && H1 == M
        (i,j) = popmax(H)
        C += (exp(gj)-k*exp(fi))*coeff(fi)*coeff(gj)
        Q += (i,j)
    end
    for (i,j) in Q
        ...
    end
    if C != 0
        g += C/((exp(g1)-M)*coeff(f1))*x^(M-exp(f1))
        ...
    end
end
*/

/*  
Sparse algorithm for k^th power using the (k-1)^st power:

output h = (f1 + f2 + ... + ft)^k:

h = 0;
g = g1^(k-1)
H = {2,1,f2*g1}
while #H > 0
    M = H1
    C = 0
    S = 0
    Q = {};
    while #H > 0 && H1 == M
        (i,j) = popmax(H)
        S += coeff(fi)*coeff(gj)
        if M >= exp(f1)
            C += (exp(gj)-k*exp(fi))*coeff(fi)*coeff(gj)
        end
        Q += (i,j)
    end
    for (i,j) in Q
        ...
    end
    if C != 0
        C /= (exp(g1)-M)*coeff(f1)
        S += C*coeff(f1)
        g += C*x^(M-exp(f1))
        ...
    end
    h += S*x^M
end
*/


slong _fmpz_mpoly_pow_fps1(
    fmpz_mpoly_t A,
    const fmpz * fcoeffs, const ulong * fexps, slong flen,
    ulong k,
    ulong cmpmask)
{
    slong N = 1;
    flint_bitcnt_t bits = A->bits;
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
    slong galloc, glen;
    fmpz_t t1, temp1;
    fmpz * S, * C;
    slong gdemote = 0;
    ulong ofmask = mpoly_overflow_mask_sp(bits);
    TMP_INIT;

    TMP_START;

    next_loc = flen + 4;   /* something bigger than heap can ever be */
    heap = TMP_ARRAY_ALLOC(flen + 1, mpoly_heap1_s);
    chain = TMP_ARRAY_ALLOC(flen, mpoly_heap_t);
    Q = TMP_ARRAY_ALLOC(3*flen, slong);
    hind = Q + 2*flen; /* space for heap indices */    
    for (i = 0; i < flen; i++)
        hind[i] = 1;

    fmpz_init(t1);
    fmpz_init(temp1);

    _fmpz_mpoly_fit_length(&Acoeffs, &Aexps, &A->alloc, 2, N);

    galloc = k*(flen - 1) + 1;
    Gexps = FLINT_ARRAY_ALLOC(galloc, ulong);
    glen = 1;
    Alen = 1;

    Gexps[0] = fexps[0]*(k - 1);

    Aexps[0] = fexps[0]*k;

    Gcoeffs = (fmpz *) flint_calloc(galloc, sizeof(fmpz));
    fmpz_pow_ui(Gcoeffs + 0, fcoeffs + 0, k - 1);
    fmpz_mul(Acoeffs + 0, Gcoeffs + 0, fcoeffs + 0);

    x = chain + 1;
    x->i = 1;
    x->j = 0;
    x->next = NULL;

    HEAP_ASSIGN(heap[1], fexps[1] + Gexps[0], x);

    hind[1] = 2*1 + 0;

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fmpz_mpoly_fit_length(&Acoeffs, &Aexps, &A->alloc, Alen + 1, N);

        if (glen >= galloc)
        {
            Gexps = FLINT_ARRAY_REALLOC(Gexps, 2*galloc, ulong);
            Gcoeffs = FLINT_ARRAY_REALLOC(Gcoeffs, 2*galloc, fmpz);
            flint_mpn_zero(Gcoeffs + galloc, galloc);
            galloc *= 2;
        }

        Aexps[Alen] = exp;
        S = Acoeffs + Alen;
        C = Gcoeffs + glen;

        fmpz_zero(C);
        fmpz_zero(S);

        Qlen = 0;

        Gexps[glen] = exp - fexps[0];

        if ((Gexps[glen] & ofmask) == 0)
        {
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do {
                    fmpz fi;
                    slong dd, ddd;
                    hind[x->i] |= 1;
                    Q[Qlen++] = i = x->i;
                    Q[Qlen++] = j = x->j;
                    fi = fcoeffs[i];

                    FLINT_ASSERT(j >= gdemote);

                    dd = (k - 1)*fexps[i] - Gexps[j];

                    if (!COEFF_IS_MPZ(fi) && !z_mul_checked(&ddd, dd, fi))
                    {
                        fmpz_addmul_si(S, Gcoeffs + j, fi);
                        fmpz_addmul_si(C, Gcoeffs + j, ddd);
                    }
                    else
                    {
                        fmpz_mul(t1, fcoeffs + i, Gcoeffs + j);
                        fmpz_add(S, S, t1);
                        fmpz_addmul_si(C, t1, dd);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }
        else
        {
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do {
                    hind[x->i] |= 1;
                    Q[Qlen++] = i = x->i;
                    Q[Qlen++] = j = x->j;
                    FLINT_ASSERT(j >= gdemote);
                    fmpz_addmul(S, fcoeffs + x->i, Gcoeffs + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        FLINT_ASSERT(Qlen <= 2*flen);
      
        while (Qlen > 0)
        {
            /* take node from store */
            j = Q[--Qlen];
            i = Q[--Qlen];

            /* should we go right? */
            if (i + 1 < flen && hind[i + 1] == 2*j + 1)
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*j + 2;
                _mpoly_heap_insert1(heap, fexps[i + 1] + Gexps[j], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            /* should we go up */
            if (j + 1 < glen && hind[i] < 2*j + 4)
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*j + 4;
                _mpoly_heap_insert1(heap, fexps[i] + Gexps[j + 1], x,
                                                &next_loc, &heap_len, cmpmask);
            }
        }

        j = hind[flen - 1]/2 - 1;
        for ( ; gdemote < j; gdemote++)
            fmpz_clear(Gcoeffs + gdemote);

        if (!fmpz_is_zero(C))
        {
            if (fmpz_is_one(fcoeffs + 0))
            {
                fmpz_divexact_si(C, C, exp - k*fexps[0]);
                fmpz_add(S, S, C);
            }
            else
            {
                fmpz_divexact_si(temp1, C, exp - k*fexps[0]);
                fmpz_add(S, S, temp1);
                fmpz_divexact(C, temp1, fcoeffs + 0);
            }

            if ((hind[1] & 1) != 0)
            {
                x = chain + 1;

                x->i = 1;
                x->j = glen;
                x->next = NULL;
                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, fexps[1] + Gexps[glen], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            glen += 1;
        }

        Alen += !fmpz_is_zero(S);
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;

    fmpz_clear(t1);
    fmpz_clear(temp1);

    for ( ; gdemote < galloc; gdemote++)
        fmpz_clear(Gcoeffs + gdemote);
    flint_free(Gcoeffs);
    flint_free(Gexps);

    TMP_END;

    return Alen;
}


slong _fmpz_mpoly_pow_fps(
    fmpz_mpoly_t A,
    const fmpz * Fcoeffs, const ulong * Fexps, slong Flen,
    ulong k,
    slong N,
    const ulong * cmpmask)
{
    flint_bitcnt_t bits = A->bits;
    ulong ofmask = (bits > FLINT_BITS) ? 0 : mpoly_overflow_mask_sp(bits);
    slong i, j, exp_next;
    slong next_loc;
    slong heap_len; /* heap zero index unused */
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
    slong Galloc, Glen;
    ulong * Gexps;
    ulong * fik, * exps;
    ulong ** exp_list;
    ulong * temp2;
    fmpz_t t1, t2, C;
    int divides;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_pow_fps1(A, Fcoeffs, Fexps, Flen, k, cmpmask[0]);

    TMP_START;

    next_loc = Flen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Flen + 1)*sizeof(mpoly_heap_s));
    chain = TMP_ARRAY_ALLOC(Flen, mpoly_heap_t);
    exps = (ulong *) TMP_ALLOC((Flen + 1)*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC((Flen + 1)*sizeof(ulong *));
    fik = (ulong *) TMP_ALLOC(N*Flen*sizeof(ulong));
    temp2 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    Q = TMP_ARRAY_ALLOC(3*Flen, slong);
    hind = Q + 2*Flen; /* space for heap indices */

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

    Galloc = 2*Flen + 3;
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
    heap_len = 2;

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

        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, Aexps + N*Alen, N))
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            do {
                Q[Qlen++] = i = x->i;
                Q[Qlen++] = j = x->j;
                hind[i] |= 1;

                fmpz_mul(t1, Fcoeffs + x->i, Gcoeffs + x->j);
                fmpz_add(S, S, t1);
                if (divides)
                {
                    mpn_sub_n(temp2, fik + x->i*N, Gexps + x->j*N, N);
                    fmpz_set_signed_ui_array(t2, temp2, N);
                    fmpz_addmul(C, t1, t2);
                }
            } while ((x = x->next) != NULL);
        }
      
        while (Qlen > 0)
        {
            /* take node from store */
            j = Q[--Qlen];
            i = Q[--Qlen];

            if (i + 1 < Flen && hind[i + 1] == 2*j + 1)
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*j + 2;

                mpoly_monomial_add_mp(exp_list[exp_next], Fexps + N*(i + 1),
                                                          Gexps + N*j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            if (j + 1 < Glen && hind[i] < 2*j + 4)
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*j + 4;

                mpoly_monomial_add_mp(exp_list[exp_next], Fexps + i*N,
                                                          Gexps + N*(j + 1), N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }

        if (!fmpz_is_zero(C))
        {
            mpoly_monomial_mul_ui_mp(temp2, Fexps + 0, N, k);
            mpn_sub_n(temp2, Aexps + N*Alen, temp2, N);
            fmpz_set_signed_ui_array(t2, temp2, N);

            fmpz_divexact(t1, C, t2);
            fmpz_add(S, S, t1);
            fmpz_divexact(Gcoeffs + Glen, t1, Fcoeffs + 0);

            if ((hind[1] & 1) != 0)
            {
                x = chain + 1;

                x->i = 1;
                x->j = Glen;
                x->next = NULL;

                hind[x->i] = 2*(Glen + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Fexps + N,
                                                          Gexps + N*Glen, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            Glen++;

            if (Glen >= Galloc)
            {
                Gexps = (ulong *) flint_realloc(Gexps, 2*N*sizeof(ulong)*Galloc);
                Gcoeffs = (fmpz *) flint_realloc(Gcoeffs, 2*sizeof(fmpz)*Galloc);
                flint_mpn_zero(Gcoeffs + Galloc, Galloc);
                Galloc *= 2;
            }
        }

        Alen += !fmpz_is_zero(Acoeffs + Alen);
        _fmpz_mpoly_fit_length(&Acoeffs, &Aexps, &A->alloc, Alen+1, N);
    }

    A->coeffs = Acoeffs;
    A->exps = Aexps;
   
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(C);

    for (i = 0; i < Galloc; i++)
        fmpz_clear(Gcoeffs + i);
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
