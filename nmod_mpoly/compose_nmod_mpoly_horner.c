/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


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
static int _nmod_mpoly_pmul(nmod_mpoly_t A, const nmod_mpoly_t X,
                  const fmpz_t pow, nmod_mpoly_t T, const nmod_mpoly_ctx_t ctx)
{
    ulong p;
    FLINT_ASSERT(fmpz_sgn(pow) > 0);

    if (!fmpz_fits_si(pow))
    {
        if (!nmod_mpoly_pow_fmpz(T, X, pow, ctx))
        {
            nmod_mpoly_zero(A, ctx);
            return 0;
        }

        nmod_mpoly_mul(A, A, T, ctx);
        return 1;
    }

    p = fmpz_get_ui(pow);

    if (X->length <= WORD(2) || A->length/p < X->length)
    {
        if (!nmod_mpoly_pow_ui(T, X, p, ctx))
        {
            nmod_mpoly_zero(A, ctx);
            return 0;
        }

        nmod_mpoly_mul(A, A, T, ctx);
    }
    else
    {
        while (p >= 1)
        {
            nmod_mpoly_mul(T, A, X, ctx);
            nmod_mpoly_swap(A, T, ctx);
            p--;
        }
    }

    return 1;
}

/*
    evaluate B(xbar) at xbar = C,
*/
int nmod_mpoly_compose_nmod_mpoly_horner(nmod_mpoly_t A,
                   const nmod_mpoly_t B, nmod_mpoly_struct * const * C,
                     const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
{
    int success = 1;
    int ret;
    slong nvars = ctxB->minfo->nvars;
    slong i, j, k, cur, next, f, r, f_prev, r_prev, v;
    slong sp, rp;
    stack_entry_struct * stack;
    nmod_mpoly_struct * regs;
    nmod_mpoly_t temp;
    slong * rtypes;
    ulong totalcounts, maxcounts;
    ulong * counts;
    slong Blen = B->length;
    slong * Blist;
    const mp_limb_t * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong BN = mpoly_words_per_exp(Bbits, ctxB->minfo);
    fmpz * Buexp;
    fmpz * mdegs;
    fmpz_t score, tz;
    TMP_INIT;

    if (Blen < 1)
    {
        nmod_mpoly_zero(A, ctxAC);
        return 1;
    }

    if (nvars < 1)
    {
       FLINT_ASSERT(Blen == 1);
       nmod_mpoly_set_ui(A, B->coeffs[0], ctxAC);
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
    regs   = (nmod_mpoly_struct *) TMP_ALLOC(nvars*sizeof(nmod_mpoly_struct));
    for (i = 0; i < nvars; i++)
        nmod_mpoly_init(regs + i, ctxAC);
    nmod_mpoly_init(temp, ctxAC);

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
        FLINT_ASSERT(Blist[f] == -WORD(1)); /* f should have had only one term */
        rtypes[rp] = f;
        goto HornerFormReturn;
    }
    else if (totalcounts == 1)
    {
        FLINT_ASSERT(!fmpz_is_zero(Buexp + nvars*f + v)); /* this term should not be a scalar */
        if (!nmod_mpoly_pow_fmpz(regs + rp, C[v], Buexp + nvars*f + v, ctxAC))
        {
            success = 0;
        }
        nmod_mpoly_scalar_mul_ui(regs + rp, regs + rp, Bcoeff[f], ctxAC);

        if (Blist[f] != -WORD(1)) /* if f has a second term */
        {
            /* this term should be a scalar */
            FLINT_ASSERT(fmpz_is_zero(Buexp + nvars*Blist[f] + v));
            nmod_mpoly_add_ui(regs + rp, regs + rp,  Bcoeff[Blist[f]], ctxAC);
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
            if (!_nmod_mpoly_pmul(regs + rp - 1, C[(stack + sp)->v_var],
                                             (stack + sp)->v_exp, temp, ctxAC))
            {
                success = 0;
            }
            nmod_mpoly_add(temp, regs + rp - 1, regs + rp, ctxAC);
            nmod_mpoly_swap(temp, regs + rp - 1, ctxAC);
        }
        else if (rtypes[rp - 1] == -WORD(1) && rtypes[rp] != -WORD(1))
        {
            /* quotient is a polynomial, remainder is a scalar */
            if (!_nmod_mpoly_pmul(regs + rp - 1, C[(stack + sp)->v_var],
                                             (stack + sp)->v_exp, temp, ctxAC))
            {
                success = 0;
            }
            nmod_mpoly_add_ui(regs + rp - 1, regs + rp - 1, Bcoeff[rtypes[rp]], ctxAC);
        }
        else if (rtypes[rp - 1] != -WORD(1) && rtypes[rp] == -WORD(1))
        {
            /* quotient is a scalar, remainder is a polynomial */
            if (!nmod_mpoly_pow_fmpz(temp, C[(stack + sp)->v_var],
                                                   (stack + sp)->v_exp, ctxAC))
            {
                success = 0;
            }
            nmod_mpoly_scalar_mul_ui(temp, temp, Bcoeff[rtypes[rp - 1]], ctxAC);
            nmod_mpoly_add(regs + rp - 1, temp, regs + rp, ctxAC);
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
        FLINT_ASSERT(rtypes[rp] == -WORD(1)); /* quotient is a scalar */

        /* quotient is a polynomial */
        if (!_nmod_mpoly_pmul(regs + rp, C[(stack + sp)->v_var],
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
        nmod_mpoly_swap(A, regs + rp, ctxAC);
    }
    else
    {
        nmod_mpoly_set_ui(A, Bcoeff[rtypes[rp]], ctxAC);
    }

cleanup:

    for (i = 0; i < nvars; i++)
        nmod_mpoly_clear(regs + i, ctxAC);
    nmod_mpoly_clear(temp, ctxAC);

    fmpz_clear(score);
    fmpz_clear(tz);
    _fmpz_vec_clear(mdegs, nvars);
    _fmpz_vec_clear(Buexp, nvars*Blen);

    TMP_END;

    return success;
}

