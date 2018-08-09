/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpolyun_init(nmod_mpolyun_t A, mp_bitcnt_t bits,
                                                    const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void nmod_mpolyun_clear(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpolyn_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void nmod_mpolyun_swap(nmod_mpolyun_t A, nmod_mpolyun_t B)
{
   nmod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyun_zero(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        nmod_mpolyn_clear(A->coeffs + i, ctx);
        nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
    }
    A->length = 0;
}

void nmod_mpolyun_print_pretty(const nmod_mpolyun_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        FLINT_ASSERT((poly->coeffs + i)->bits == poly->bits);
        nmod_mpolyn_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void nmod_mpolyun_fit_length(nmod_mpolyun_t A, slong length,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpolyn_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpolyn_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
        }
        A->alloc = new_alloc;
    }
}

void nmod_mpolyun_shift_right(nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void nmod_mpolyun_shift_left(nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

slong nmod_mpolyun_lastdeg(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            deg = FLINT_MAX(deg, nmod_poly_degree((A->coeffs + i)->coeffs + j));
        }
    }
    FLINT_ASSERT(deg >= 0);
    return deg;
}

void nmod_mpolyun_set(nmod_mpolyun_t A, const nmod_mpolyun_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_set(Acoeff + i, Bcoeff + i, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


/* take the last variable of B out */
void nmod_mpoly_cvtto_mpolyn(nmod_mpolyn_t A, nmod_mpoly_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    ulong * oneexp;
    slong offset;
    slong shift;
    ulong mask;
    slong N;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    oneexp = TMP_ALLOC(N*sizeof(ulong));
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var,
                                                       N, B->bits, ctx->minfo);

    nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    k = 0;
    nmod_mpolyn_fit_length(A, k + 1, ctx);
    for (i = 0; i < B->length; i++)
    {
        ulong c = (B->exps[N*i + offset] >> shift) & mask;
        mpoly_monomial_msub(A->exps + N*k, B->exps + N*i, c, oneexp, N);

        if (k > 0 && mpoly_monomial_equal(A->exps + N*k, A->exps + N*(k - 1), N))
        {
            nmod_poly_set_coeff_ui(A->coeffs + k - 1, c, B->coeffs[i]);
        } else
        {
            nmod_poly_zero(A->coeffs + k);
            nmod_poly_set_coeff_ui(A->coeffs + k, c, B->coeffs[i]);
            k++;
            nmod_mpolyn_fit_length(A, k + 1, ctx);
        }
    }

    nmod_mpolyn_set_length(A, k, ctx);
    TMP_END;
}

void nmod_mpolyu_cvtto_mpolyun(nmod_mpolyun_t A, nmod_mpolyu_t B,
                                           slong k, const nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff;
    nmod_mpoly_struct * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpoly_cvtto_mpolyn(Acoeff + i, Bcoeff + i, k, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;  
}




/* put the last variable of B back into A */
void nmod_mpoly_cvtfrom_mpolyn(nmod_mpoly_t A, nmod_mpolyn_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong k;
    slong N;
    ulong * oneexp;
    slong offset;
    slong shift;
    TMP_INIT;

    FLINT_ASSERT(B->bits == A->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    oneexp = TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var,
                                                       N, B->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, B->length, ctx);

    k = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = (B->coeffs + i)->length - 1; j >= 0; j--)
        {
            mp_limb_t c = (B->coeffs + i)->coeffs[j];
            if (c != UWORD(0))
            {
                nmod_mpoly_fit_length(A, k + 1, ctx);
                A->coeffs[k] = c;
                mpoly_monomial_madd(A->exps + N*k, B->exps + N*i, j, oneexp, N);                
                k++;
            }
        }
    }

    A->length = k;
    TMP_END;
}

void nmod_mpolyu_cvtfrom_mpolyun(nmod_mpolyu_t A, nmod_mpolyun_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpolyu_fit_length(A, B->length, ctx);

    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_cvtfrom_mpolyn(A->coeffs + i, B->coeffs + i, var, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}



void nmod_mpolyn_set_mpoly(nmod_mpolyn_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;
    nmod_poly_struct * Acoeff;
    mp_limb_t * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    Blen = B->length;
    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        nmod_poly_zero(Acoeff + i);
        nmod_poly_set_coeff_ui(Acoeff + i, 0, Bcoeff[i]);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_poly_clear(Acoeff + i);
        nmod_poly_init(Acoeff + i, ctx->ffinfo->mod.n);
    }
    A->length = Blen;
}

void nmod_mpolyun_set_mpolyu(nmod_mpolyun_t A, const nmod_mpolyu_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);
    nmod_mpolyun_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        nmod_mpolyn_set_mpoly(A->coeffs + i, B->coeffs + i, ctx);

        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);

    }
    A->length = B->length;
}



void nmod_mpolyun_mul_poly(nmod_mpolyun_t A, const nmod_mpolyun_t B,
                              const nmod_poly_t c, const nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_mul_poly(Acoeff + i, Bcoeff + i, c, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


void nmod_mpolyun_content_last(nmod_poly_t a, nmod_mpolyun_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    nmod_poly_zero(a);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            nmod_poly_gcd(a, a, (B->coeffs + i)->coeffs + j);
            if (nmod_poly_degree(a) == 0)
                break;
        }
    }
}

void nmod_mpolyun_divexact_last(nmod_mpolyun_t A, nmod_poly_t b,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    nmod_poly_t r;
    nmod_poly_init(r, ctx->ffinfo->mod.n);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            nmod_poly_divrem((A->coeffs + i)->coeffs + j, r,
                             (A->coeffs + i)->coeffs + j, b);
            FLINT_ASSERT(nmod_poly_is_zero(r));
        }
    }
    nmod_poly_clear(r);
}


/* Convert Ap to A using the lowest degree representative */
void nmod_mpolyn_set_fq_nmod_mpoly(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx,
                            fq_nmod_mpoly_t Ap, const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;
    slong N;

    FLINT_ASSERT(Ap->bits == A->bits);
    nmod_mpolyn_fit_length(A, Ap->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    for (i = 0; i < Ap->length; i++) {
        mpoly_monomial_set(A->exps + N*i, Ap->exps + N*i, N);
        nmod_poly_set(A->coeffs + i, Ap->coeffs + i);
    }
    A->length = Ap->length;
}

void nmod_mpolyun_set_fq_nmod_mpolyu(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx,
                           fq_nmod_mpolyu_t Ap, const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;

    FLINT_ASSERT(Ap->bits == A->bits);
    nmod_mpolyun_fit_length(A, Ap->length, ctx);
    for (i = 0; i < Ap->length; i++)
    {
        A->exps[i] = Ap->exps[i];
        nmod_mpolyn_set_fq_nmod_mpoly(A->coeffs + i, ctx, Ap->coeffs + i, ctxp);
    }
    A->length = Ap->length;
}


/*
    Update H so that it does not change mod m, and is now A mod p
    It is asserted that the monomials in H and A match
*/
int nmod_mpolyn_CRT_fq_nmod_mpoly(slong * lastdeg,
                          nmod_mpolyn_t H, const nmod_mpoly_ctx_t ctx,
                          nmod_poly_t m, fq_nmod_t inv_m_eval,
                             fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;
    slong N;
    int changed = 0;
    fq_nmod_t u, v;
    nmod_poly_t w;

    fq_nmod_init(u, ctxp->fqctx);
    fq_nmod_init(v, ctxp->fqctx);
    nmod_poly_init(w, ctxp->fqctx->modulus->mod.n);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        nmod_poly_rem(u, H->coeffs + i, ctxp->fqctx->modulus);
        fq_nmod_sub(v, A->coeffs + i, u, ctxp->fqctx);
        if (!fq_nmod_is_zero(v, ctxp->fqctx))
        {
            changed = 1;
            fq_nmod_mul(u, v, inv_m_eval, ctxp->fqctx);
            nmod_poly_mul(w, u, m);
            nmod_poly_add(H->coeffs + i, H->coeffs + i, w);
        }

        lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(H->coeffs + i));

        FLINT_ASSERT(nmod_poly_degree(H->coeffs + i) < nmod_poly_degree(m)
                                                     + nmod_poly_degree(ctxp->fqctx->modulus));
    }

    fq_nmod_clear(u, ctxp->fqctx);
    fq_nmod_clear(v, ctxp->fqctx);
    nmod_poly_clear(w);

    return changed;
}

int nmod_mpolyun_CRT_fq_nmod_mpolyu(slong * lastdeg,
                        nmod_mpolyun_t H, const nmod_mpoly_ctx_t ctx,
             nmod_poly_t m, fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;
    int changed = 0;
    fq_nmod_t inv_m_eval;

    lastdeg[0] = -WORD(1);

    fq_nmod_init(inv_m_eval, ctxp->fqctx);
    nmod_poly_rem(inv_m_eval, m, ctxp->fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ctxp->fqctx);

    FLINT_ASSERT(H->bits == A->bits);
    FLINT_ASSERT(H->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(H->exps[i] == A->exps[i]);
        changed |= nmod_mpolyn_CRT_fq_nmod_mpoly(lastdeg, H->coeffs + i, ctx, m,
                                              inv_m_eval, A->coeffs + i, ctxp);
    }
    H->length = A->length;
    fq_nmod_clear(inv_m_eval, ctxp->fqctx);
    return changed;
}


/* reduce B via evaluation of the last variable */
void nmod_mpolyn_redto_fq_nmod_mpoly(fq_nmod_mpoly_t A, nmod_mpolyn_t B,
                   const fq_nmod_mpoly_ctx_t ffctx, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ffctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ffctx->minfo->nvars == ctx->minfo->nvars);

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    k = 0;
    fq_nmod_mpoly_fit_length(A, k + 1, ffctx);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_fit_length(A, k + 1, ffctx);
        mpoly_monomial_set(A->exps + N*k, B->exps + N*i, N);
        nmod_poly_rem(A->coeffs + k, B->coeffs + i, ffctx->fqctx->modulus);
        k += !fq_nmod_is_zero(A->coeffs + k, ffctx->fqctx);
    }

    fq_nmod_mpoly_set_length(A, k, ffctx);
}

void nmod_mpolyun_redto_fq_nmod_mpolyu(fq_nmod_mpolyu_t A, nmod_mpolyun_t B,
                   const fq_nmod_mpoly_ctx_t ffctx, const nmod_mpoly_ctx_t ctx)
{
    slong i, k, Blen;
    fq_nmod_mpoly_struct * Acoeff;
    nmod_mpolyn_struct * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyu_fit_length(A, Blen, ffctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    k = 0;
    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_redto_fq_nmod_mpoly(Acoeff + k, Bcoeff + i, ffctx, ctx);
        Aexp[k] = Bexp[i];
        k += !fq_nmod_mpoly_is_zero(Acoeff + i, ffctx);
    }

    for (i = k; i < A->length; i++)
    {
        fq_nmod_mpoly_clear(Acoeff + i, ffctx);
        fq_nmod_mpoly_init(Acoeff + i, ffctx);
        fq_nmod_mpoly_fit_bits(Acoeff + i, A->bits, ffctx);
        (Acoeff + i)->bits = A->bits;
    }
    A->length = k;  
}


/* evaluate A at lastvar = alpha */
void nmod_mpolyn_eval_last(nmod_mpoly_t B, nmod_mpolyn_t A, mp_limb_t alpha,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;
    slong k;
    FLINT_ASSERT(B->bits == A->bits);

    nmod_mpoly_fit_length(B, A->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(B->exps + N*k, A->exps + N*i, N);
        B->coeffs[k] = nmod_poly_evaluate_nmod(A->coeffs + i, alpha);
        if (B->coeffs[k] != UWORD(0))
        {
            k++;
        }
    }
    B->length = k;
}

void nmod_mpolyun_eval_last(nmod_mpolyu_t B, nmod_mpolyun_t A, mp_limb_t alpha,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;

    FLINT_ASSERT(B->bits == A->bits);

    nmod_mpolyu_fit_length(B, A->length, ctx);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        B->exps[k] = A->exps[i];
        nmod_mpolyn_eval_last(B->coeffs + k, A->coeffs + i, alpha, ctx);
        if (!nmod_mpoly_is_zero(B->coeffs + k, ctx))
        {
            k++;
        }
    }
    B->length = k;
}



/*
    F = F + modulus*(A - F(alpha))
    no assumptions about matching monomials
*/
int
nmod_mpolyn_addinterp(slong * lastdeg,
            nmod_mpolyn_t F, nmod_mpolyn_t T, nmod_mpoly_t A,
            nmod_poly_t modulus,  mp_limb_t alpha,  const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    mp_limb_t v;
    mp_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    mp_limb_t * Acoeff = A->coeffs;
    nmod_poly_struct * Fcoeff = F->coeffs;
    nmod_poly_struct * Tcoeff;
    nmod_poly_t tp;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(modulus));

            /* F term ok, A term missing */
            v = nmod_poly_evaluate_nmod(Fcoeff + i, alpha);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_sub(Tcoeff + k, Fcoeff + i, tp);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (Acoeff[j] != UWORD(0))
            {
                changed = 1;
                nmod_poly_zero(Tcoeff + k);
                nmod_poly_scalar_mul_nmod(Tcoeff + k, modulus, Acoeff[j]);
                lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(modulus));

            /* F term ok, A term ok */
            v = nmod_poly_evaluate_nmod(Fcoeff + i, alpha);
            v = nmod_sub(Acoeff[j], v, ctx->ffinfo->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_add(Tcoeff + k, Fcoeff + i, tp);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    nmod_poly_clear(tp);
    return changed;
}

int nmod_mpolyun_addinterp(slong * lastdeg,
             nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_mpolyu_t A,
              nmod_poly_t modulus, mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    nmod_mpolyn_t S;
    nmod_mpolyn_struct * Tcoeff;
    nmod_mpolyn_struct * Fcoeff;
    nmod_mpoly_struct  * Acoeff;
    nmod_mpoly_t zero;

    lastdeg[0] = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    nmod_mpolyn_init(S, F->bits, ctx);

    Flen = F->length;
    Alen = A->length;
    nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    nmod_mpoly_init(zero, ctx);
    nmod_mpoly_fit_bits(zero, A->bits, ctx);
    zero->bits = A->bits;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_addinterp(lastdeg, Tcoeff + k,
                                                 S, zero, modulus, alpha, ctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= nmod_mpolyn_addinterp(lastdeg, Tcoeff + k,
                                           S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_addinterp(lastdeg, Tcoeff + k,
                                           S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    /* demote remaining coefficients */
    for (i = k; i < T->length; i++)
    {
        nmod_mpolyn_clear(Tcoeff + i, ctx);
        nmod_mpolyn_init(Tcoeff + i, T->bits, ctx);
    }
    T->length = k;

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    nmod_mpolyn_clear(S, ctx);
    nmod_mpoly_clear(zero, ctx);
    return changed;    
}


/*
    F = F + modulus*((A - F(alpha))/(modulus(alpha)))
    no assumptions about matching monomials
*/
int
nmod_mpolyn_addinterp_fq_nmod_mpoly(slong * lastdeg,
                         nmod_mpolyn_t F, nmod_mpolyn_t T, nmod_poly_t m,
                         const nmod_mpoly_ctx_t ctx, fq_nmod_mpoly_t A,
                         fq_nmod_t inv_m_eval, const fq_nmod_mpoly_ctx_t ffctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t u, v;
    nmod_poly_t w;
    mp_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fq_nmod_struct * Acoeff = A->coeffs;
    nmod_poly_struct * Fcoeff = F->coeffs;
    nmod_poly_struct * Tcoeff;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(u, ffctx->fqctx);
    fq_nmod_init(v, ffctx->fqctx);
    nmod_poly_init(w, ffctx->fqctx->modulus->mod.n);

    nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(m));

            /* F term ok, A term missing */
            nmod_poly_rem(v, Fcoeff + i, ffctx->fqctx->modulus);
            if (!fq_nmod_is_zero(v, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ffctx->fqctx);
                nmod_poly_mul(w, u, m);
                nmod_poly_sub(Tcoeff + k, Fcoeff + i, w);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!fq_nmod_is_zero(Acoeff + j, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, Acoeff + j, inv_m_eval, ffctx->fqctx);
                nmod_poly_mul(Tcoeff + k, m, u);
                lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(m));

            /* F term ok, A term ok */
            nmod_poly_rem(u, Fcoeff + i, ffctx->fqctx->modulus);
            fq_nmod_sub(v, Acoeff + j, u, ffctx->fqctx);
            if (!fq_nmod_is_zero(v, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ffctx->fqctx);
                nmod_poly_mul(w, m, u);
                nmod_poly_add(Tcoeff + k, Fcoeff + i, w);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, ffctx->fqctx);
    fq_nmod_clear(v, ffctx->fqctx);
    nmod_poly_clear(w);

    return changed;
}


int nmod_mpolyun_addinterp_fq_nmod_mpolyu(slong * lastdeg,
                            nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_poly_t m,
                            const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t ffctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    nmod_mpolyn_t S;
    nmod_mpolyn_struct * Tcoeff;
    nmod_mpolyn_struct * Fcoeff;
    fq_nmod_mpoly_struct  * Acoeff;
    fq_nmod_mpoly_t zero;
    fq_nmod_t inv_m_eval;

    lastdeg[0] = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    nmod_mpolyn_init(S, F->bits, ctx);

    Flen = F->length;
    Alen = A->length;
    nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    fq_nmod_mpoly_init(zero, ffctx);
    fq_nmod_mpoly_fit_bits(zero, A->bits, ffctx);
    zero->bits = A->bits;

    fq_nmod_init(inv_m_eval, ffctx->fqctx);
    nmod_poly_rem(inv_m_eval, m, ffctx->fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ffctx->fqctx);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_addinterp_fq_nmod_mpoly(lastdeg, Tcoeff + k,
                                           S, m, ctx, zero, inv_m_eval, ffctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= nmod_mpolyn_addinterp_fq_nmod_mpoly(lastdeg, Tcoeff + k,
                                     S, m, ctx, Acoeff + j, inv_m_eval, ffctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_addinterp_fq_nmod_mpoly(lastdeg, Tcoeff + k,
                                     S, m, ctx, Acoeff + j, inv_m_eval, ffctx);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
            j++;
        } else 
        {

            FLINT_ASSERT(0);
        }
    }

    /* demote remaining coefficients */
    for (i = k; i < T->length; i++)
    {
        nmod_mpolyn_clear(Tcoeff + i, ctx);
        nmod_mpolyn_init(Tcoeff + i, T->bits, ctx);
    }
    T->length = k;

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    fq_nmod_clear(inv_m_eval, ffctx->fqctx);

    nmod_mpolyn_clear(S, ctx);
    fq_nmod_mpoly_clear(zero, ffctx);
    return changed;    
}
