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
#include "gr.h"
#include "gr_generic.h"

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
/* todo: cache squares, ...? */
static int _gr_pmul(gr_ptr A, gr_srcptr X, const fmpz_t pow, gr_ptr T, gr_ctx_t ctx)
{
    if (fmpz_is_one(pow))
    {
        return gr_mul(A, A, X, ctx);
    }
    else
    {
        int status;
        status = gr_pow_fmpz(T, X, pow, ctx);
        status |= gr_mul(A, A, T, ctx);
        return status;
    }
}

int
gr_fmpz_mpoly_evaluate_horner(gr_ptr A, const fmpz_mpoly_t B, gr_srcptr C, const fmpz_mpoly_ctx_t ctxB, gr_ctx_t ctx)
{
    int ret;
    slong nvars = ctxB->minfo->nvars;
    slong i, j, k, cur, next, f, r, f_prev, r_prev, v;
    slong sp, rp;
    stack_entry_struct * stack;
    gr_ptr regs;
    gr_ptr temp;
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
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    TMP_INIT;

    if (Blen == 0)
    {
        return gr_zero(A, ctx);
    }

    if (Blen == 1 && fmpz_mpoly_is_fmpz(B, ctxB))
    {
        return gr_set_fmpz(A, B->coeffs, ctx);
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

    /* registers of grs */
    rp = 0;
    rtypes = (slong *) TMP_ALLOC((nvars + 1)*sizeof(slong));
    GR_TMP_INIT_VEC(regs, nvars, ctx);
    GR_TMP_INIT(temp, ctx);

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

        status |= gr_pow_fmpz(GR_ENTRY(regs, rp, sz), GR_ENTRY(C, v, sz), Buexp + nvars*f + v, ctx);
        status |= gr_mul_fmpz(GR_ENTRY(regs, rp, sz), GR_ENTRY(regs, rp, sz), Bcoeff + f, ctx);

        if (Blist[f] != -WORD(1)) /* if f has a second term */
        {
            /* this term should be a scalar */
            FLINT_ASSERT(fmpz_is_zero(Buexp + nvars*Blist[f] + v));
            status |= gr_add_fmpz(GR_ENTRY(regs, rp, sz), GR_ENTRY(regs, rp, sz),  Bcoeff + Blist[f], ctx);
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
            status |= _gr_pmul(GR_ENTRY(regs, rp - 1, sz), GR_ENTRY(C, (stack + sp)->v_var, sz), (stack + sp)->v_exp, temp, ctx);
            status |= gr_add(temp, GR_ENTRY(regs, rp - 1, sz), GR_ENTRY(regs, rp, sz), ctx);
            gr_swap(temp, GR_ENTRY(regs, rp - 1, sz), ctx);
        }
        else if (rtypes[rp - 1] == -WORD(1) && rtypes[rp] != -WORD(1))
        {
            /* quotient is a polynomial, remainder is a scalar */
            status |= _gr_pmul(GR_ENTRY(regs, rp - 1, sz), GR_ENTRY(C, (stack + sp)->v_var, sz), (stack + sp)->v_exp, temp, ctx);
            status |= gr_add_fmpz(GR_ENTRY(regs, rp - 1, sz), GR_ENTRY(regs, rp - 1, sz), Bcoeff + rtypes[rp], ctx);
        }
        else if (rtypes[rp - 1] != -WORD(1) && rtypes[rp] == -WORD(1))
        {
            /* quotient is a scalar, remainder is a polynomial */
            status |= gr_pow_fmpz(temp, GR_ENTRY(C, (stack + sp)->v_var, sz), (stack + sp)->v_exp, ctx);
            status |= gr_mul_fmpz(temp, temp, Bcoeff + rtypes[rp - 1], ctx);
            status |= gr_add(GR_ENTRY(regs, rp - 1, sz), temp, GR_ENTRY(regs, rp, sz), ctx);
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
        status |= _gr_pmul(GR_ENTRY(regs, rp, sz), GR_ENTRY(C, (stack + sp)->v_var, sz), (stack + sp)->v_exp, temp, ctx);
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
        gr_swap(A, GR_ENTRY(regs, rp, sz), ctx);
    }
    else
    {
        status |= gr_set_fmpz(A, Bcoeff + rtypes[rp], ctx);
    }

    GR_TMP_CLEAR_VEC(regs, nvars, ctx);
    GR_TMP_CLEAR(temp, ctx);

    fmpz_clear(score);
    fmpz_clear(tz);
    _fmpz_vec_clear(mdegs, nvars);
    _fmpz_vec_clear(Buexp, nvars*Blen);

    TMP_END;

    return status;
}

/* todo: accept fmpz exponents */
int
gr_evaluate_fmpz_mpoly_iter(gr_ptr res, const fmpz_mpoly_t pol, gr_srcptr x, const fmpz_mpoly_ctx_t ctx, gr_ctx_t cactx)
{
    slong i, j, len, nvars;
    gr_ptr s, t, u;
    ulong * exp;
    int status = GR_SUCCESS;
    slong sz = cactx->sizeof_elem;

    len = fmpz_mpoly_length(pol, ctx);

    if (len == 0)
        return gr_zero(res, cactx);

    if (len == 1 && fmpz_mpoly_is_fmpz(pol, ctx))
        return gr_set_fmpz(res, pol->coeffs, cactx);

    nvars = ctx->minfo->nvars;
    exp = flint_malloc(sizeof(ulong) * nvars);

    GR_TMP_INIT3(s, t, u, cactx);

    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_get_term_exp_ui(exp, pol, i, ctx);

        status |= gr_one(t, cactx);

        for (j = 0; j < nvars; j++)
        {
            if (exp[j] == 1)
            {
                status |= gr_mul(t, t, GR_ENTRY(x, j, sz), cactx);
            }
            else if (exp[j] >= 2)
            {
                status |= gr_pow_ui(u, GR_ENTRY(x, j, sz), exp[j], cactx);
                status |= gr_mul(t, t, u, cactx);
            }
        }

        status |= gr_mul_fmpz(t, t, pol->coeffs + i, cactx);
        status |= gr_add(s, s, t, cactx);
    }

    gr_swap(res, s, cactx);

    flint_free(exp);

    GR_TMP_CLEAR3(s, t, u, cactx);

    return status;
}

int
gr_fmpz_mpoly_evaluate(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t ctx, gr_ctx_t cactx)
{
    return gr_fmpz_mpoly_evaluate_horner(res, f, x, ctx, cactx);
}
