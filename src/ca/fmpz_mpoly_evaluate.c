/*
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "ca.h"

/*
The conversion to Horner form can be stated as recursive. However, the call
stack has depth proportial to the length of the input polynomial in the worst
case. Therefore, we must convert it to an iterative algorithm.

The proceedure is

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
/* todo: cache squares, ...? */
static void _ca_pmul(ca_t A, const ca_t X,
                  const fmpz_t pow, ca_t T, ca_ctx_t ctx)
{
    if (fmpz_is_one(pow))
    {
        ca_mul(A, A, X, ctx);
    }
    else
    {
        ca_pow_fmpz(T, X, pow, ctx);
        ca_mul(A, A, T, ctx);
    }
}

void
ca_fmpz_mpoly_evaluate_horner(ca_t A, const fmpz_mpoly_t B, ca_srcptr C, const fmpz_mpoly_ctx_t ctxB, ca_ctx_t ctx)
{
    int ret;
    slong nvars = ctxB->minfo->nvars;
    slong i, j, k, cur, next, f, r, f_prev, r_prev, v;
    slong sp, rp;
    stack_entry_struct * stack;
    ca_struct * regs;
    ca_t temp;
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

    if (Blen == 0)
    {
        ca_zero(A, ctx);
        return;
    }

    if (Blen == 1 && fmpz_mpoly_is_fmpz(B, ctxB))
    {
        ca_set_fmpz(A, B->coeffs, ctx);
        return;
    }

    /* flint_printf("========================== HORNER %wd ==========================\n", Blen); */

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

    /* registers of qqbars */
    rp = 0;
    rtypes = (slong *) TMP_ALLOC((nvars + 1)*sizeof(slong));
    regs   = (ca_struct *) TMP_ALLOC(nvars*sizeof(ca_struct));
    for (i = 0; i < nvars; i++)
        ca_init(regs + i, ctx);
    ca_init(temp, ctx);

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

        ca_pow_fmpz(regs + rp, C + v, Buexp + nvars*f + v, ctx);
        ca_mul_fmpz(regs + rp, regs + rp, Bcoeff + f, ctx);

        if (Blist[f] != -WORD(1)) /* if f has a second term */
        {
            /* this term should be a scalar */
            FLINT_ASSERT(fmpz_is_zero(Buexp + nvars*Blist[f] + v));
            ca_add_fmpz(regs + rp, regs + rp,  Bcoeff + Blist[f], ctx);
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
            _ca_pmul(regs + rp - 1, C + (stack + sp)->v_var, (stack + sp)->v_exp, temp, ctx);
            ca_add(temp, regs + rp - 1, regs + rp, ctx);
            ca_swap(temp, regs + rp - 1, ctx);
        }
        else if (rtypes[rp - 1] == -WORD(1) && rtypes[rp] != -WORD(1))
        {
            /* quotient is a polynomial, remainder is a scalar */
            _ca_pmul(regs + rp - 1, C + (stack + sp)->v_var, (stack + sp)->v_exp, temp, ctx);
            ca_add_fmpz(regs + rp - 1, regs + rp - 1, Bcoeff + rtypes[rp], ctx);
        }
        else if (rtypes[rp - 1] != -WORD(1) && rtypes[rp] == -WORD(1))
        {
            /* quotient is a scalar, remainder is a polynomial */
            ca_pow_fmpz(temp, C + (stack + sp)->v_var, (stack + sp)->v_exp, ctx);
            ca_mul_fmpz(temp, temp, Bcoeff + rtypes[rp - 1], ctx);
            ca_add(regs + rp - 1, temp, regs + rp, ctx);
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
        _ca_pmul(regs + rp, C + (stack + sp)->v_var, (stack + sp)->v_exp, temp, ctx);
    }

    rtypes[rp] = -WORD(1);

HornerFormReturn:

    ret = (stack + sp)->ret;
    fmpz_clear((stack + sp)->v_exp);
    sp--;
    if (ret == 1) goto HornerForm1;
    if (ret == 2) goto HornerForm2;

    FLINT_ASSERT(rp == 0);
    FLINT_ASSERT(sp == -WORD(1));

    if (rtypes[rp] == -WORD(1))
    {
        ca_swap(A, regs + rp, ctx);
    }
    else
    {
        ca_set_fmpz(A, Bcoeff + rtypes[rp], ctx);
    }

    for (i = 0; i < nvars; i++)
        ca_clear(regs + i, ctx);
    ca_clear(temp, ctx);

    fmpz_clear(score);
    fmpz_clear(tz);
    _fmpz_vec_clear(mdegs, nvars);
    _fmpz_vec_clear(Buexp, nvars*Blen);

    TMP_END;
}

/* todo: accept fmpz exponents */
void
ca_evaluate_fmpz_mpoly_iter(ca_t res, const fmpz_mpoly_t pol, ca_srcptr x, const fmpz_mpoly_ctx_t ctx, ca_ctx_t cactx)
{
    slong i, j, len, nvars;
    ca_t s, t, u;
    ulong * exp;

    len = fmpz_mpoly_length(pol, ctx);

    if (len == 0)
    {
        ca_zero(res, cactx);
        return;
    }

    if (len == 1 && fmpz_mpoly_is_fmpz(pol, ctx))
    {
        ca_set_fmpz(res, pol->coeffs, cactx);
        return;
    }

    nvars = ctx->minfo->nvars;
    exp = flint_malloc(sizeof(ulong) * nvars);

    ca_init(s, cactx);
    ca_init(t, cactx);
    ca_init(u, cactx);

    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_get_term_exp_ui(exp, pol, i, ctx);

        ca_one(t, cactx);

        for (j = 0; j < nvars; j++)
        {
            if (exp[j] == 1)
            {
                ca_mul(t, t, x + j, cactx);
            }
            else if (exp[j] >= 2)
            {
                ca_pow_ui(u, x + j, exp[j], cactx);
                ca_mul(t, t, u, cactx);
            }
        }

        ca_mul_fmpz(t, t, pol->coeffs + i, cactx);
        ca_add(s, s, t, cactx);
    }

    ca_swap(res, s, cactx);

    flint_free(exp);

    ca_clear(s, cactx);
    ca_clear(t, cactx);
    ca_clear(u, cactx);
}

void
ca_fmpz_mpoly_evaluate(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t ctx, ca_ctx_t cactx)
{
    ca_fmpz_mpoly_evaluate_horner(res, f, x, ctx, cactx);
}
