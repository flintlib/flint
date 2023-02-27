/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"


void nmod_mpolyun_init(
    nmod_mpolyun_t A,
    flint_bitcnt_t bits,
    const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void nmod_mpolyun_clear(
    nmod_mpolyun_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpolyn_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

int nmod_mpolyun_is_canonical(
    const nmod_mpolyun_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (A->length > A->alloc)
    {
        return 0;
    }

    for (i = 0; i < A->length; i++)
    {
        if (!nmod_mpolyn_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }
    }

    return 1;
}

void nmod_mpolyun_print_pretty(
    const nmod_mpolyun_t poly,
    const char ** x,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        nmod_mpolyn_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
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

int nmod_mpolyn_is_nonzero_nmod(const nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length != WORD(1))
    {
        return 0;
    }

    if (n_poly_degree(A->coeffs + 0) != 0)
    {
        return 0;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    return mpoly_monomial_is_zero(A->exps + N*0, N);
}

int nmod_mpolyun_is_nonzero_nmod(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
    {
        return 0;
    }

    return nmod_mpolyn_is_nonzero_nmod(A->coeffs + 0, ctx);
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
            deg = FLINT_MAX(deg, n_poly_degree((A->coeffs + i)->coeffs + j));
        }
    }
    FLINT_ASSERT(deg >= 0);
    return deg;
}

slong nmod_mpolyn_lastdeg(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        deg = FLINT_MAX(deg, n_poly_degree(A->coeffs + i));
    }

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

void nmod_mpolyn_one(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    n_poly_struct * Acoeff;
    ulong * Aexp;
    slong N;

    nmod_mpolyn_fit_length(A, 1, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    n_poly_one(Acoeff + 0);
    mpoly_monomial_zero(Aexp + N*0, N);

    A->length = 1;
}


void nmod_mpolyun_one(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(ctx->mod.n > 1);
    nmod_mpolyun_fit_length(A, 1, ctx);
    nmod_mpolyn_one(A->coeffs + 0, ctx);
    A->exps[0] = 0;
    A->length = 1;
}

void nmod_mpolyn_set_mod(nmod_mpolyn_t A, const nmod_t mod)
{
}

void nmod_mpolyun_set_mod(nmod_mpolyun_t A, const nmod_t mod)
{
}

void nmod_mpolyn_scalar_mul_nmod(
    nmod_mpolyn_t A,
    mp_limb_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (c == 1)
        return;

    for (i = 0; i < A->length; i++)
        _n_poly_mod_scalar_mul_nmod_inplace(A->coeffs + i, c, ctx->mod);
}

void nmod_mpolyun_scalar_mul_nmod(
    nmod_mpolyun_t A,
    mp_limb_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    FLINT_ASSERT(c != 0);
    for (i = 0; i < A->length; i++)
    {
        nmod_mpolyn_scalar_mul_nmod(A->coeffs + i, c, ctx);
    }
}


void nmod_mpolyn_mul_last(
    nmod_mpolyn_t A,
    n_poly_t b,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    n_poly_t t;

    FLINT_ASSERT(!n_poly_is_zero(b));

    if (n_poly_is_one(b))
        return;

    n_poly_init(t);

    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_mul(t, A->coeffs + i, b, ctx->mod);
        n_poly_swap(t, A->coeffs + i);
    }

    n_poly_clear(t);
}

/*
    A *= b

    A is in R[X][x_0,..., x_(v-1)][x_v]
    b is in R[x_v]
*/
void nmod_mpolyun_mul_last(
    nmod_mpolyun_t A,
    n_poly_t b,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    n_poly_t t;

    FLINT_ASSERT(!n_poly_is_zero(b));

    if (n_poly_is_one(b))
        return;

    n_poly_init(t);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            n_poly_mod_mul(t, (A->coeffs + i)->coeffs + j, b, ctx->mod);
            n_poly_swap(t, (A->coeffs + i)->coeffs + j);
        }
    }

    n_poly_clear(t);
}


int nmod_mpolyn_equal(
    const nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (!mpoly_monomial_equal(A->exps + N*i, B->exps + N*i, N))
        {
            return 0;
        }
        if (!n_poly_equal(A->coeffs + i, B->coeffs + i))
        {
            return 0;
        }
    }
    return 1;
}

int nmod_mpolyun_equal(
    const nmod_mpolyun_t A,
    const nmod_mpolyun_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
        {
            return 0;
        }
        if (!nmod_mpolyn_equal(A->coeffs + i, B->coeffs + i, ctx))
        {
            return 0;
        }
    }
    return 1;
}



/* if the coefficient doesn't exist, a new one is created (and set to zero) */
n_poly_struct * _nmod_mpolyn_get_coeff(nmod_mpolyn_t A,
                                      ulong * pow, const nmod_mpoly_ctx_t uctx)
{
    slong i, j, a, b;
    n_poly_struct * xk;
    slong N = mpoly_words_per_exp_sp(A->bits, uctx->minfo);
    int cmp;

    a = 0;
    b = A->length;

    if (b == 0 || mpoly_monomial_gt_nomask(pow, A->exps + N*0, N))
    {
        i = 0;
        goto create_new;
    }

    if (mpoly_monomial_equal(pow , A->exps + N*(b - 1), N))
    {
        return A->coeffs + b - 1;
    }

try_again:

    if (b - a < 4)
    {
        for (i = a; i < b && (cmp = mpoly_monomial_cmp_nomask(A->exps + N*i, pow, N)) >= 0; i++)
        {
            if (cmp == 0) 
            {
                return A->coeffs + i;
            }
        }
        goto create_new;
    }
    else
    {
        i = a + (b - a)/2;
        cmp = mpoly_monomial_cmp_nomask(A->exps + N*i, pow, N);
        if (cmp == 0) 
        {
            return A->coeffs + i;
        }
        else if (cmp > 0)
        {
            a = i;
        }
        else
        {
            b = i;
        }
        goto try_again;
    }

create_new:

    nmod_mpolyn_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        mpoly_monomial_set(A->exps + N*j, A->exps + N*(j - 1), N);
        n_poly_swap(A->coeffs + j, A->coeffs + j - 1);
    }

    mpoly_monomial_set(A->exps + N*i, pow, N);
    A->length++;
    xk = A->coeffs + i;
    xk->length = 0;

    return xk;
}



/* if the coefficient doesn't exist, a new one is created (and set to zero) */
nmod_mpolyn_struct * _nmod_mpolyun_get_coeff(nmod_mpolyun_t A,
                                        ulong pow, const nmod_mpoly_ctx_t uctx)
{
    slong i, j, a, b;
    nmod_mpolyn_struct * xk;

    a = 0;
    b = A->length;

    if (b == 0 || pow > A->exps[0])
    {
        i = 0;
        goto create_new;
    }

    if (pow == A->exps[b - 1])
    {
        return A->coeffs + b - 1;
    }

try_again:

    if (b - a < 8)
    {
        for (i = a; i < b && A->exps[i] >= pow; i++)
        {
            if (A->exps[i] == pow) 
            {
                return A->coeffs + i;
            }
        }
        goto create_new;
    }
    else
    {
        i = a + (b - a)/2;
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
        else if (A->exps[i] > pow)
        {
            a = i;
        }
        else
        {
            b = i;
        }
        goto try_again;
    }

create_new:

    nmod_mpolyun_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        nmod_mpolyn_swap(A->coeffs + j, A->coeffs + j - 1);
    }
    
    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}



void nmod_mpoly_to_mpolyun_perm_deflate_bivar(
    nmod_mpolyun_t A,
    const nmod_mpoly_t B,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const nmod_mpoly_ctx_t uctx,
    const nmod_mpoly_ctx_t ctx)
{
    slong j;
    slong NB, NA;
    slong p0 = perm[0], p1 = perm[1];
    ulong shift0 = shift[p0], shift1 = shift[p1];
    ulong stride0 = stride[p0], stride1 = stride[p1];
    ulong Bexp0, Bexp1;
    slong Boff0, Bshift0, Boff1, Bshift1;
    ulong mask;
    nmod_mpolyn_struct * Ac;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(1 == uctx->minfo->nvars);
    NA = mpoly_words_per_exp_sp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&Boff0, &Bshift0, p0, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boff1, &Bshift1, p1, B->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);

    for (j = 0; j < B->length; j++)
    {
        Bexp0 = ((B->exps + NB*j)[Boff0] >> Bshift0) & mask;
        Bexp1 = ((B->exps + NB*j)[Boff1] >> Bshift1) & mask;

        Ac = _nmod_mpolyun_get_coeff(A, stride0 == 1 ? (Bexp0 - shift0)
                                           : (Bexp0 - shift0) / stride0, uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        if (Ac->length == 0)
        {
            nmod_mpolyn_fit_length(Ac, 1, uctx);
            n_poly_zero(Ac->coeffs + 0);
        }
        Ac->length = 1;

        n_poly_set_coeff(Ac->coeffs + 0, stride1 == 1 ? (Bexp1 - shift1)
                                   : (Bexp1 - shift1) / stride1, B->coeffs[j]);

        mpoly_monomial_zero(Ac->exps + NA*0, NA);
    }
}


/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k < m + 1
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses k = 0
    the coefficients of A use variables k = 1 ... m
    the variable corresponding to k = m is moved to dense storage.
*/
void nmod_mpoly_to_mpolyun_perm_deflate(
    nmod_mpolyun_t A,
    const nmod_mpoly_ctx_t uctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong j, k, l;
    slong NA = mpoly_words_per_exp_sp(A->bits, uctx->minfo);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    ulong * Bexps;
    ulong * texp;
    slong * offs, * shifts;
    nmod_mpolyn_struct * Ac;
    n_poly_struct * Acc;
    TMP_INIT;

    A->length = 0;

    if (m == 1)
    {
        nmod_mpoly_to_mpolyun_perm_deflate_bivar(A, B,
                                               perm, shift, stride, uctx, ctx);
        return;
    }

    if (m > 2)
    {
        nmod_mpolyu_t Au;
        nmod_mpolyu_init(Au, A->bits, uctx);
        nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, B, ctx,
                                    perm, shift, stride, handles, num_handles);
        nmod_mpolyu_cvtto_mpolyun(A, Au, m - 1, uctx);
        nmod_mpolyu_clear(Au, uctx);
        return;
    }

    TMP_START;
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    texp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));

    offs   = (slong *) TMP_ALLOC(m*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(m*sizeof(slong));
    for (k = 0; k < m; k++)
    {
        mpoly_gen_offset_shift_sp(offs + k, shifts + k, k, A->bits, uctx->minfo);
    }

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        l = perm[0];
        Ac = _nmod_mpolyun_get_coeff(A, stride[l] == 1 ? (Bexps[l] - shift[l]) : (Bexps[l] - shift[l]) / stride[l], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        mpoly_monomial_zero(texp, NA);
        for (k = 1; k < m; k++)
        {
            l = perm[k];
            texp[offs[k - 1]] += (stride[l] == 1 ? (Bexps[l] - shift[l]) : (Bexps[l] - shift[l]) / stride[l]) << shifts[k - 1];
        }

        Acc = _nmod_mpolyn_get_coeff(Ac, texp, uctx);
        l = perm[m];
        n_poly_set_coeff(Acc, stride[l] == 1 ? (Bexps[l] - shift[l]) : (Bexps[l] - shift[l]) / stride[l], B->coeffs[j]);
    }

    TMP_END;
}


void nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(
    nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t nctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong j, k, l;
    slong NA = mpoly_words_per_exp_sp(A->bits, nctx->minfo);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong n = ctx->minfo->nvars;
    slong m = nctx->minfo->nvars;
    ulong * Bexps;
    slong * offs, * shifts;
    nmod_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(m <= n);

    TMP_START;
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    offs   = (slong *) TMP_ALLOC(m*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(m*sizeof(slong));
    for (k = 0; k < m; k++)
    {
        mpoly_gen_offset_shift_sp(offs + k, shifts + k, k, A->bits, nctx->minfo);
    }

    nmod_mpoly_init3(T, B->length, A->bits, nctx);
    T->length = B->length;
    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        T->coeffs[j] = B->coeffs[j];
        mpoly_monomial_zero(T->exps + NA*j, NA);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            (T->exps + NA*j)[offs[k]] += ((Bexps[l] - shift[l]) / stride[l]) << shifts[k];
        }
    }

    nmod_mpoly_sort_terms(T, nctx);

    nmod_mpoly_cvtto_mpolyn(A, T, nctx->minfo->nvars - 1, nctx);

    nmod_mpoly_clear(T, nctx);

    TMP_END;
}



/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k <= m
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void nmod_mpoly_from_mpolyun_perm_inflate(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t ctx,
    const nmod_mpolyun_t B,
    const nmod_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, h, k, l;
    slong NA, NB;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    ulong * uexps;
    ulong * Aexps, * tAexp, * tAgexp;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    tAexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    tAgexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    mpoly_gen_monomial_sp(tAgexp, perm[m], Abits, ctx->minfo);
    for (i = 0; i < NA; i++)
    {
        tAgexp[i] *= stride[perm[m]];
    }

    nmod_mpoly_fit_length_reset_bits(A, 0, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        nmod_mpolyn_struct * Bc = B->coeffs + i;
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            mpoly_get_monomial_ui(uexps + 1, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i];
            FLINT_ASSERT(uexps[m] == 0);
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 1; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }

            mpoly_set_monomial_ui(tAexp, Aexps, Abits, ctx->minfo);

            l = perm[m];
            h = (Bc->coeffs + j)->length;
            _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                                   &Aexp, &A->exps_alloc, NA, Alen + h);
            for (h--; h >= 0; h--)
            {
                mp_limb_t c = (Bc->coeffs + j)->coeffs[h];
                if (c == 0)
                    continue;
                mpoly_monomial_madd(Aexp + NA*Alen, tAexp, h, tAgexp, NA);
                Acoeff[Alen] = c;
                Alen++;
            }
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    _nmod_mpoly_set_length(A, Alen, ctx);

    nmod_mpoly_sort_terms(A, ctx);
    TMP_END;
}

void nmod_mpoly_from_mpolyn_perm_inflate(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t ctx,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t nctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = nctx->minfo->nvars;
    slong i, h, k, l;
    slong NA, NB;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    ulong * Bexps;
    ulong * Aexps, * tAexp, * tAgexp;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);
    TMP_START;

    Bexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, nctx->minfo);

    tAexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    tAgexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    mpoly_gen_monomial_sp(tAgexp, perm[m - 1], Abits, ctx->minfo);
    for (i = 0; i < NA; i++)
        tAgexp[i] *= stride[perm[m - 1]];

    nmod_mpoly_fit_length_reset_bits(A, 0, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, nctx->minfo);
        FLINT_ASSERT(Bexps[m - 1] == 0);
        for (l = 0; l < n; l++)
        {
            Aexps[l] = shift[l];
        }
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[l] += stride[l]*Bexps[k];
        }

        mpoly_set_monomial_ui(tAexp, Aexps, Abits, ctx->minfo);

        h = (B->coeffs + i)->length;
        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, NA, Alen + h);
        for (h--; h >= 0; h--)
        {
            mp_limb_t c = (B->coeffs + i)->coeffs[h];
            if (c == 0)
                continue;
            mpoly_monomial_madd(Aexp + NA*Alen, tAexp, h, tAgexp, NA);
            Acoeff[Alen] = c;
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    _nmod_mpoly_set_length(A, Alen, ctx);

    nmod_mpoly_sort_terms(A, ctx);
    TMP_END;
}


/* take the last variable of B out */
void nmod_mpoly_cvtto_mpolyn(nmod_mpolyn_t A, const nmod_mpoly_t B,
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

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift, var,
                                                          B->bits, ctx->minfo);

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
            n_poly_set_coeff(A->coeffs + k - 1, c, B->coeffs[i]);
        }
        else
        {
            n_poly_zero(A->coeffs + k);
            n_poly_set_coeff(A->coeffs + k, c, B->coeffs[i]);
            k++;
            nmod_mpolyn_fit_length(A, k + 1, ctx);
        }
    }

    nmod_mpolyn_set_length(A, k, ctx);
    TMP_END;
}

void nmod_mpolyu_cvtto_mpolyun(nmod_mpolyun_t A, const nmod_mpolyu_t B,
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

    A->length = Blen;  
}




/* put the last variable of B back into A */
void nmod_mpoly_cvtfrom_mpolyn(
    nmod_mpoly_t A,
    const nmod_mpolyn_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    ulong * genexp;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(genexp, var, B->bits, ctx->minfo);

    nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

    k = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->coeffs[i].length - 1; j >= 0; j--)
        {
            mp_limb_t c = B->coeffs[i].coeffs[j];
            if (c == 0)
                continue;

            _nmod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                                   &A->exps, &A->exps_alloc, N, k + 1);
            A->coeffs[k] = c;
            mpoly_monomial_madd(A->exps + N*k, B->exps + N*i, j, genexp, N);                
            k++;
        }
    }

    A->length = k;
    TMP_END;
}

void nmod_mpolyu_cvtfrom_mpolyun(nmod_mpolyu_t A, const nmod_mpolyun_t B,
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

void nmod_mpolyun_mul_poly(nmod_mpolyun_t A, const nmod_mpolyun_t B,
                                  const n_poly_t c, const nmod_mpoly_ctx_t ctx)
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


void nmod_mpolyun_content_last(
    n_poly_t a,
    nmod_mpolyun_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    n_poly_zero(a);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            n_poly_mod_gcd(a, a, (B->coeffs + i)->coeffs + j, ctx->mod);
            if (n_poly_degree(a) == 0)
                break;
        }
    }
}

void nmod_mpolyn_content_last(
    n_poly_t a,
    nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    n_poly_zero(a);
    for (i = 0; i < B->length; i++)
    {
        n_poly_mod_gcd(a, a, B->coeffs + i, ctx->mod);
        if (n_poly_degree(a) == 0)
            break;
    }
}

void nmod_mpolyun_divexact_last(
    nmod_mpolyun_t A,
    n_poly_t b,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    n_poly_t r;

    if (n_poly_is_one(b))
        return;

    n_poly_init(r);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            n_poly_mod_divrem((A->coeffs + i)->coeffs + j, r,
                              (A->coeffs + i)->coeffs + j, b, ctx->mod);
            FLINT_ASSERT(n_poly_is_zero(r));
        }
    }

    n_poly_clear(r);
}

void nmod_mpolyn_divexact_last(
    nmod_mpolyn_t A,
    n_poly_t b,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    n_poly_t r;

    if (n_poly_is_one(b))
        return;

    n_poly_init(r);

    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_divrem(A->coeffs + i, r, A->coeffs + i, b, ctx->mod);
        FLINT_ASSERT(n_poly_is_zero(r));
    }

    n_poly_clear(r);
}

