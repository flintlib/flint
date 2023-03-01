/*
    Copyright (C) 2018, 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* evaluate B(xbar) at xbar = C */
int fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A,
                     const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    slong i;
    fmpz_mat_t M;

    FLINT_ASSERT(A != B);

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctxAC);
        return 1;
    }

    fmpz_mat_init(M, ctxAC->minfo->nfields + 1, ctxB->minfo->nfields);
    fmpz_mat_zero(M);

    for (i = 0; i < ctxB->minfo->nvars; i++)
    {
        if (C[i]->length > 1)
            goto matrix_no_good;

        if (C[i]->length == 0)
        {
            mpoly_compose_mat_fill_column(M, NULL, 0, i,
                                                    ctxB->minfo, ctxAC->minfo);
        }
        else
        {
            if (!fmpz_is_one(C[i]->coeffs + 0))
                goto matrix_no_good;

            mpoly_compose_mat_fill_column(M, C[i]->exps, C[i]->bits, i,
                                                    ctxB->minfo, ctxAC->minfo);
        }
    }

    _fmpz_mpoly_compose_mat(A, B, M, ctxB, ctxAC);

    fmpz_mat_clear(M);

    return 1;

matrix_no_good:

    fmpz_mat_clear(M);

    for (i = 0; i < ctxB->minfo->nvars; i++)
    {
        if (C[i]->length > 1)
        {
            return fmpz_mpoly_compose_fmpz_mpoly_horner(A, B, C, ctxB, ctxAC);
        }
    }

    return fmpz_mpoly_compose_fmpz_mpoly_geobucket(A, B, C, ctxB, ctxAC);
}

/* evaluate B(x_1,...,x_n) at x_i = y_c[i], y_j are vars of ctxAC */
void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A,
                             const fmpz_mpoly_t B, const slong * c,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    fmpz_mat_t M;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctxAC);
        return;
    }

    fmpz_mat_init(M, ctxAC->minfo->nfields + 1, ctxB->minfo->nfields);
    mpoly_compose_mat_gen(M, c, ctxB->minfo, ctxAC->minfo);

    if (A == B)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init(T, ctxAC);
        _fmpz_mpoly_compose_mat(T, B, M, ctxB, ctxAC);
        fmpz_mpoly_swap(A, T, ctxAC);
        fmpz_mpoly_clear(T, ctxAC);
    }
    else
    {
        _fmpz_mpoly_compose_mat(A, B, M, ctxB, ctxAC);
    }

    fmpz_mat_clear(M);

    return;
}

/* evaluate B(xbar) at xbar = C */
int fmpz_mpoly_compose_fmpz_mpoly_geobucket(fmpz_mpoly_t A,
                  const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    int success = 1;
    slong i, j;
    slong Blen = B->length;
    const fmpz * Bcoeff = B->coeffs;
    const ulong * Bexp = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong BN = mpoly_words_per_exp(Bbits, ctxB->minfo);
    fmpz_mpoly_t U, V, W;
    fmpz_mpoly_geobucket_t T;
    fmpz * e;

    fmpz_mpoly_init(U, ctxAC);
    fmpz_mpoly_init(V, ctxAC);
    fmpz_mpoly_init(W, ctxAC);
    fmpz_mpoly_geobucket_init(T, ctxAC);
    e = _fmpz_vec_init(ctxB->minfo->nvars);

    for (i = 0; success && i < Blen; i++)
    {
        fmpz_mpoly_set_fmpz(U, Bcoeff + i, ctxAC);
        mpoly_get_monomial_ffmpz(e, Bexp + BN*i, Bbits, ctxB->minfo);
        for (j = 0; j < ctxB->minfo->nvars; j++)
        {
            success = success && fmpz_mpoly_pow_fmpz(V, C[j], e + j, ctxAC);
            fmpz_mpoly_mul(W, U, V, ctxAC);
            fmpz_mpoly_swap(U, W, ctxAC);
        }
        fmpz_mpoly_geobucket_add(T, U, ctxAC);
    }

    if (success)
        fmpz_mpoly_geobucket_empty(A, T, ctxAC);

    fmpz_mpoly_clear(U, ctxAC);
    fmpz_mpoly_clear(V, ctxAC);
    fmpz_mpoly_clear(W, ctxAC);
    fmpz_mpoly_geobucket_clear(T, ctxAC);
    _fmpz_vec_clear(e, ctxB->minfo->nvars);

    return success;
}

/*
The conversion to Horner form can be stated as recursive. However, the call
stack has depth proportial to the length of the input polynomial in the worst
case. Therefore, we must convert it to an iterative algorithm.

The procedure is

HornerForm(f):

    if f is simple to evaluate

        return eval(f)

    else
        choose a variable v and the smallest non zero exponent e appearing
            in the terms of f

        write f = q * v^e + r  where r is independent of the variable v

        return  HornerForm(q) * v^e + HornerForm(r)
*/

typedef struct
{
    slong f;
    slong r;
    slong v_var;
    fmpz_t v_exp;   /* will be managed as stack grows / shrinks */
    int ret;
} stack_entry_struct;

typedef stack_entry_struct stack_entry_t[1];


/* A = A * X^pow */
static int _fmpz_mpoly_pmul(fmpz_mpoly_t A, const fmpz_mpoly_t X,
                  const fmpz_t pow, fmpz_mpoly_t T, const fmpz_mpoly_ctx_t ctx)
{
    ulong p;
    FLINT_ASSERT(fmpz_sgn(pow) > 0);

    if (!fmpz_fits_si(pow))
    {
        if (!fmpz_mpoly_pow_fmpz(T, X, pow, ctx))
        {
            fmpz_mpoly_zero(A, ctx);
            return 0;
        }

        fmpz_mpoly_mul(A, A, T, ctx);
        return 1;
    }

    p = fmpz_get_ui(pow);

    if (X->length <= WORD(2) || A->length/p < X->length)
    {
        if (!fmpz_mpoly_pow_ui(T, X, p, ctx))
        {
            fmpz_mpoly_zero(A, ctx);
            return 0;
        }

        fmpz_mpoly_mul(A, A, T, ctx);
    }
    else
    {
        while (p >= 1)
        {
            fmpz_mpoly_mul(T, A, X, ctx);
            fmpz_mpoly_swap(A, T, ctx);
            p--;
        }
    }

    return 1;
}

/* evaluate B(xbar) at xbar = C */
int fmpz_mpoly_compose_fmpz_mpoly_horner(fmpz_mpoly_t A,
                  const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    int success = 1;
    int ret;
    slong nvars = ctxB->minfo->nvars;
    slong i, j, k, cur, next, f, r, f_prev, r_prev, v;
    slong sp, rp;
    stack_entry_struct * stack;
    fmpz_mpoly_struct * regs;
    fmpz_mpoly_t temp;
    slong * rtypes;
    ulong totalcounts, maxcounts;
    ulong * counts;
    slong Blen = B->length;
    slong * Blist;
    const fmpz * Bcoeff = B->coeffs;
    const ulong * Bexp = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong BN = mpoly_words_per_exp(Bbits, ctxB->minfo);
    fmpz * Buexp;
    fmpz * mdegs;
    fmpz_t score, tz;
    TMP_INIT;

    if (Blen < 1)
    {
        fmpz_mpoly_zero(A, ctxAC);
        return 1;
    }

    if (nvars < 1)
    {
       FLINT_ASSERT(Blen == 1);
       fmpz_mpoly_set_fmpz(A, B->coeffs + 0, ctxAC);
       return 1;
    }

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    fmpz_init(score);
    fmpz_init(tz);

    /* unpack B exponents */
    Buexp = _fmpz_vec_init(nvars*Blen);
    for (i = 0; i < Blen; i++)
        mpoly_get_monomial_ffmpz(Buexp + nvars*i, Bexp + BN*i, Bbits, ctxB->minfo);

    counts = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    mdegs = _fmpz_vec_init(nvars);

    /* stack */
    sp = -WORD(1); /* start with empty stack */
    stack = (stack_entry_struct *) TMP_ALLOC(nvars*(Blen + 1)*sizeof(stack_entry_struct));
    Blist = (slong *) TMP_ALLOC(Blen*sizeof(slong));

    /* registers of polynomials */
    rp = 0;
    rtypes = (slong *) TMP_ALLOC((nvars + 1)*sizeof(slong));
    regs   = (fmpz_mpoly_struct *) TMP_ALLOC(nvars*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < nvars; i++)
        fmpz_mpoly_init(regs + i, ctxAC);
    fmpz_mpoly_init(temp, ctxAC);

    /* polynomials will be stored as link lists */
    for (i = 0; i + 1 < Blen; i++)
        Blist[i] = i + 1;
    Blist[i] = -WORD(1);

    sp++;
    fmpz_init((stack + sp)->v_exp);
    (stack + sp)->ret = 0;
    (stack + sp)->f = 0;

HornerForm:

    f = (stack + sp)->f;

    FLINT_ASSERT(f != -WORD(1)); /* f is not supposed to be zero */

    /* obtain a count of the number of terms containing each variable */
    for (i = 0; i < nvars; i++)
    {
        counts[i] = 0;
        fmpz_set_si(mdegs + i, -WORD(1));
    }

    for (j = f; j != -WORD(1); j = Blist[j])
    {
        for (i = 0; i < nvars; i++)
        {
            if (!fmpz_is_zero(Buexp + nvars*j + i ))
            {
                counts[i]++;
                if (fmpz_sgn(mdegs + i) < 0
                    || fmpz_cmp(mdegs + i, Buexp + nvars*j + i) > 0)
                {
                    fmpz_set(mdegs + i, Buexp + nvars*j + i);
                }
            }
        }
    }

    totalcounts = 0;
    maxcounts = 0;
    v = -WORD(1);
    for (i = 0; i < nvars; i++)
    {
        maxcounts = FLINT_MAX(maxcounts, counts[i]);
        totalcounts += counts[i];
        if (counts[i] != 0)
            v = i;
    }

    /* handle simple cases */
    if (totalcounts == 0)
    {
        FLINT_ASSERT(Blist[f] == -WORD(1));    /* f should have had only one term */
        rtypes[rp] = f;
        goto HornerFormReturn;
    }
    else if (totalcounts == 1)
    {
        FLINT_ASSERT(!fmpz_is_zero(Buexp + nvars*f + v)); /* this term should not be a scalar */
        if (!fmpz_mpoly_pow_fmpz(regs + rp, C[v], Buexp + nvars*f + v, ctxAC))
        {
            success = 0;
        }
        fmpz_mpoly_scalar_mul_fmpz(regs + rp, regs + rp, Bcoeff + f, ctxAC);

        if (Blist[f] != -WORD(1)) /* if f has a second term */
        {
            /* this term should be a scalar */
            FLINT_ASSERT(fmpz_is_zero(Buexp + nvars*Blist[f] + v));
            fmpz_mpoly_add_fmpz(regs + rp, regs + rp,  Bcoeff + Blist[f], ctxAC);
        }   

        rtypes[rp] = -WORD(1);

        goto HornerFormReturn;
    }

    /* pick best power to pull out */
    k = 0;
    if (maxcounts == 1)
    {
        fmpz_set_si(score, -WORD(1));
        for (i = 0; i < nvars; i++)
        {
            if (counts[i] == 1 && (fmpz_sgn(score) < 0
                                   || fmpz_cmp(mdegs + i, score) < 0))
            {
                FLINT_ASSERT(fmpz_sgn(mdegs + i) > 0);
                fmpz_set(score, mdegs + i);
                k = i;
            }
        }
    }
    else
    {
        fmpz_zero(score);
        for (i = 0; i < nvars; i++)
        {
            if (counts[i] > 1)
            {
                FLINT_ASSERT(fmpz_sgn(mdegs + i) > 0);
                fmpz_mul_ui(tz, mdegs + i, counts[i] - 1);
                if (fmpz_cmp(tz, score) > 0)
                {
                    fmpz_swap(score, tz);
                    k = i;
                }
            }
        }
    }

    /* set variable power v */
    (stack + sp)->v_var = k;
    fmpz_set((stack + sp)->v_exp, mdegs + k);

    /* scan f and split into q and v with f = q*v + r then set f = q */
    r = -WORD(1);
    cur = f;
    f_prev = -WORD(1);
    r_prev = -WORD(1);
    while (cur != -WORD(1))
    {
        next = Blist[cur];
        if (fmpz_is_zero(Buexp + nvars*cur + k))
        {
            if (f_prev == -WORD(1))
                f = Blist[cur];
            else
                Blist[f_prev] = Blist[cur];

            if (r_prev == -WORD(1))
                r = cur;
            else
                Blist[r_prev] = cur;

            Blist[cur] = -WORD(1);
            r_prev = cur;
        }
        else
        {
            /* mdegs[k] should be minimum non zero exponent */
            fmpz_sub(Buexp + nvars*cur + k, Buexp + nvars*cur + k, mdegs + k);
            FLINT_ASSERT(fmpz_sgn(Buexp + nvars*cur + k) >= 0);
            f_prev = cur;
        }
        cur = next;
    }
    (stack + sp)->r = r;

    /* convert the quotient */
    sp++;
    fmpz_init((stack + sp)->v_exp);
    (stack + sp)->ret = 1;
    (stack + sp)->f = f;
    goto HornerForm;

HornerForm1:

    /* convert the remainder */
    r = (stack + sp)->r;
    if (r != -WORD(1))
    {
        /* remainder is non zero */
        rp++;
        FLINT_ASSERT(0 <= rp && rp <= nvars);
        sp++;
        fmpz_init((stack + sp)->v_exp);
        (stack + sp)->ret = 2;
        (stack + sp)->f = r;
        goto HornerForm;

HornerForm2:

        if (rtypes[rp - 1] == -WORD(1) && rtypes[rp] == -WORD(1))
        {
            /* both quotient and remainder are polynomials */
            if (!_fmpz_mpoly_pmul(regs + rp - 1, C[(stack + sp)->v_var],
                                             (stack + sp)->v_exp, temp, ctxAC))
            {
                success = 0;
            }
            fmpz_mpoly_add(temp, regs + rp - 1, regs + rp, ctxAC);
            fmpz_mpoly_swap(temp, regs + rp - 1, ctxAC);
        }
        else if (rtypes[rp - 1] == -WORD(1) && rtypes[rp] != -WORD(1))
        {
            /* quotient is a polynomial, remainder is a scalar */
            if (!_fmpz_mpoly_pmul(regs + rp - 1, C[(stack + sp)->v_var],
                                             (stack + sp)->v_exp, temp, ctxAC))
            {
                success = 0;
            }
            fmpz_mpoly_add_fmpz(regs + rp - 1, regs + rp - 1,
                                                   Bcoeff + rtypes[rp], ctxAC);
        }
        else if (rtypes[rp - 1] != -WORD(1) && rtypes[rp] == -WORD(1))
        {
            /* quotient is a scalar, remainder is a polynomial */
            if (!fmpz_mpoly_pow_fmpz(temp, C[(stack + sp)->v_var],
                                                   (stack + sp)->v_exp, ctxAC))
            {
                success = 0;
            }
            fmpz_mpoly_scalar_mul_fmpz(temp, temp, Bcoeff + rtypes[rp - 1], ctxAC);
            fmpz_mpoly_add(regs + rp - 1, temp, regs + rp, ctxAC);
        }
        else
        {
            /* quotient is a scalar, remainder is a scalar */
            FLINT_ASSERT(0);    /* this should have been handled by simple case */
        }
        rp--;
        FLINT_ASSERT(0 <= rp && rp <= nvars);
    }
    else
    {
        /* remainder is zero */
        FLINT_ASSERT(rtypes[rp] == -WORD(1)); /* quotient is not a scalar */

        /* quotient is a polynomial */
        if (!_fmpz_mpoly_pmul(regs + rp, C[(stack + sp)->v_var],
                                            (stack + sp)->v_exp, temp, ctxAC))
        {
            success = 0;
        }
    }

    rtypes[rp] = -WORD(1);

HornerFormReturn:

    if (!success)
    {
        while (sp >= 0)
        {
            fmpz_clear((stack + sp)->v_exp);
            sp--;
        }

        goto cleanup;
    }

    ret = (stack + sp)->ret;
    fmpz_clear((stack + sp)->v_exp);
    sp--;
    if (ret == 1) goto HornerForm1;
    if (ret == 2) goto HornerForm2;

    FLINT_ASSERT(rp == 0);
    FLINT_ASSERT(sp == -WORD(1));

    if (rtypes[rp] == -WORD(1))
    {
        fmpz_mpoly_swap(A, regs + rp, ctxAC);
    }
    else
    {
        fmpz_mpoly_set_fmpz(A, Bcoeff + rtypes[rp], ctxAC);
    }

cleanup:

    for (i = 0; i < nvars; i++)
        fmpz_mpoly_clear(regs + i, ctxAC);
    fmpz_mpoly_clear(temp, ctxAC);

    fmpz_clear(score);
    fmpz_clear(tz);
    _fmpz_vec_clear(mdegs, nvars);
    _fmpz_vec_clear(Buexp, nvars*Blen);

    TMP_END;

    return success;
}

static int _fmpz_poly_pow_fmpz_is_not_feasible(const fmpz_poly_t b, const fmpz_t e)
{
    if (b->length < 2)
    {
        if (b->length == 1)
        {
            return _fmpz_pow_fmpz_is_not_feasible(fmpz_bits(b->coeffs + 0), e);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return fmpz_cmp_ui(e, limit/(ulong)(b->length)) >= 0;
    }
}

static int _fmpz_poly_pow_ui_is_not_feasible(const fmpz_poly_t b, ulong e)
{
    if (b->length < 2)
    {
        if (b->length == 1)
        {
            return _fmpz_pow_ui_is_not_feasible(fmpz_bits(b->coeffs + 0), e);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return e >= limit/(ulong)(b->length);
    }
}

int _fmpz_mpoly_compose_fmpz_poly_sp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, shift, off;
    slong Blen = B->length;
    const fmpz * Bcoeff = B->coeffs;
    const ulong * Bexp = B->exps;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    fmpz_poly_struct * powers;
    fmpz_poly_t t, t2;
    TMP_INIT;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    degrees = TMP_ARRAY_ALLOC(nvars, slong);
    mpoly_degrees_si(degrees, Bexp, Blen, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_poly_pow_ui_is_not_feasible(C[i], degrees[i]))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += FLINT_BIT_COUNT(degrees[i]);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpz_poly_struct);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits = FLINT_BIT_COUNT(degrees[i]);

        mpoly_gen_offset_shift_sp(&off, &shift, i, bits, ctx->minfo);
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off;
            masks[k] = UWORD(1) << (shift + j);
            fmpz_poly_init(powers + k);
            if (j == 0)
                fmpz_poly_set(powers + k, C[i]);
            else
                fmpz_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpz_poly_zero(A);
    fmpz_poly_init(t);
    fmpz_poly_init(t2);
    for (i = 0; i < Blen; i++)
    {
        fmpz_poly_set_fmpz(t, Bcoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fmpz_poly_mul(t2, t, powers + k);
                fmpz_poly_swap(t, t2);
            }
        }
        fmpz_poly_add(A, A, t);
    }
    fmpz_poly_clear(t);
    fmpz_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        fmpz_poly_clear(powers + k);

cleanup_degrees:

    TMP_END;

    return success;
}


int _fmpz_mpoly_compose_fmpz_poly_mp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    ulong l;
    slong i, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, off;
    slong Blen = B->length;
    fmpz * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    fmpz_poly_struct * powers;
    fmpz_poly_t t, t2;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);

    TMP_START;

    degrees = TMP_ARRAY_ALLOC(nvars, fmpz);
    for (i = 0; i < nvars; i++)
        fmpz_init(degrees + i);

    mpoly_degrees_ffmpz(degrees, Bexp, Blen, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_poly_pow_fmpz_is_not_feasible(C[i], degrees + i))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += fmpz_bits(degrees + i);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpz_poly_struct);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits = fmpz_bits(degrees + i);

        off = mpoly_gen_offset_mp(i, bits, ctx->minfo);

        for (l = 0; l < varibits; l++)
        {
            offs[k] = off + (l / FLINT_BITS);
            masks[k] = UWORD(1) << (l % FLINT_BITS);
            fmpz_poly_init(powers + k);
            if (l == 0)
                fmpz_poly_set(powers + k, C[i]);
            else
                fmpz_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpz_poly_zero(A);
    fmpz_poly_init(t);
    fmpz_poly_init(t2);
    for (i = 0; i < Blen; i++)
    {
        fmpz_poly_set_fmpz(t, Bcoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fmpz_poly_mul(t2, t, powers + k);
                fmpz_poly_swap(t, t2);
            }
        }
        fmpz_poly_add(A, A, t);
    }
    fmpz_poly_clear(t);
    fmpz_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        fmpz_poly_clear(powers + k);

cleanup_degrees:

    for (i = 0; i < nvars; i++)
        fmpz_clear(degrees + i);

    TMP_END;

    return success;
}


int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fmpz_poly_zero(A);
        return 1;
    }

    if (B->bits <= FLINT_BITS)
    {
        return _fmpz_mpoly_compose_fmpz_poly_sp(A, B, C, ctx);
    }
    else
    {
        return _fmpz_mpoly_compose_fmpz_poly_mp(A, B, C, ctx);
    }
}

/* essentially exps(A) = M*exps(B) */
void _fmpz_mpoly_compose_mat(fmpz_mpoly_t A,
                            const fmpz_mpoly_t B, const fmpz_mat_t M,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    slong i;
    fmpz * u, * v;
    flint_bitcnt_t vbits;
    slong Blen = B->length;
    flint_bitcnt_t Bbits = B->bits;
    slong BN = mpoly_words_per_exp(Bbits, ctxB->minfo);
    const ulong * Bexp = B->exps;
    const fmpz * Bcoeffs = B->coeffs;
    slong Alen_old = A->length;
    slong AN;

    FLINT_ASSERT(A != B);

    FLINT_ASSERT(fmpz_mat_nrows(M) == ctxAC->minfo->nfields + 1);
    FLINT_ASSERT(fmpz_mat_ncols(M) == ctxB->minfo->nfields);

    u = _fmpz_vec_init(ctxB->minfo->nfields);
    v = _fmpz_vec_init(ctxAC->minfo->nfields + 1);

    fmpz_mpoly_fit_length(A, Blen, ctxAC);
    A->length = 0;
    fmpz_mpoly_fit_bits(A, MPOLY_MIN_BITS, ctxAC);
    A->bits = MPOLY_MIN_BITS;
    for (i = 0; i < Blen; i++)
    {
        mpoly_unpack_vec_fmpz(u, Bexp + BN*i, Bbits, ctxB->minfo->nfields, 1);
        fmpz_mat_mul_vec(v, M, u);
        if (!fmpz_is_zero(v + ctxAC->minfo->nfields))
            continue;
        vbits = _fmpz_vec_max_bits(v, ctxAC->minfo->nfields);
        FLINT_ASSERT(vbits >= 0);
        fmpz_mpoly_fit_bits(A, mpoly_fix_bits(vbits + 1, ctxAC->minfo), ctxAC);
        fmpz_set(A->coeffs + A->length, Bcoeffs + i);
        AN = mpoly_words_per_exp(A->bits, ctxAC->minfo);
        mpoly_pack_vec_fmpz(A->exps + AN*A->length, v, A->bits, ctxAC->minfo->nfields, 1);
        A->length++;
    }

    while (--Alen_old >= A->length)
        _fmpz_demote(A->coeffs + Alen_old);

    _fmpz_vec_clear(u, ctxB->minfo->nfields);
    _fmpz_vec_clear(v, ctxAC->minfo->nfields + 1);

    fmpz_mpoly_sort_terms(A, ctxAC);
    fmpz_mpoly_combine_like_terms(A, ctxAC);
    return;
}
