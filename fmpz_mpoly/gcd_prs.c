/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "assert.h"


/* gcd of two monomials */
void fmpz_mpoly_gcd_monomial(fmpz_mpoly_t poly1, const fmpz_mpoly_t polyA,
                          const fmpz_mpoly_t polyB, const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, bits;
    ulong * texpA, * texpB, * exps;
    ulong mask;
    fmpz_t igcd;
    TMP_INIT;

    assert(polyA->length == 1);
    assert(polyB->length == 1);
    if (polyA->bits > FLINT_BITS || polyB->bits > FLINT_BITS)
        flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_gcd_monomial");

    TMP_START;
    texpA = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    texpB = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    exps = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));

    fmpz_init(igcd);
    fmpz_gcd(igcd, polyA->coeffs + 0, polyB->coeffs + 0);

    bits = FLINT_MAX(polyA->bits, polyB->bits);
    N = mpoly_words_per_exp(bits, ctx->minfo);
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    mpoly_repack_monomials(texpA, bits, polyA->exps, polyA->bits, 1, ctx->minfo);
    mpoly_repack_monomials(texpB, bits, polyB->exps, polyB->bits, 1, ctx->minfo);
    mpoly_monomial_min(texpA, texpA, texpB, bits, N, mask);
    mpoly_get_monomial_ui(exps, texpA, bits, ctx->minfo);
    
    fmpz_mpoly_fit_length(poly1, 1, ctx);
    fmpz_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;
    mpoly_set_monomial_ui(poly1->exps + N*0, exps, bits, ctx->minfo);
    fmpz_set(poly1->coeffs + 0, igcd);
    _fmpz_mpoly_set_length(poly1, 1, ctx);

    fmpz_clear(igcd);

    TMP_END;
}



int _fmpz_mpoly_gcd_prs(fmpz_mpoly_t poly1, const fmpz_mpoly_t polyA,
                          const fmpz_mpoly_t polyB, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    slong shift, off, bits, N;
    slong i, m, c, d, v, k, var, nvars = ctx->minfo->nvars;
    ulong mask;
    slong * a_degs, * b_degs, * a_leads, * b_leads;
    fmpz_mpoly_t ac, bc, gc, gabc, g;
    fmpz_mpoly_univar_t ax, bx, gx;
    TMP_INIT;

    TMP_START;

    fmpz_mpoly_init(ac, ctx);
    fmpz_mpoly_init(bc, ctx);
    fmpz_mpoly_init(gc, ctx);
    fmpz_mpoly_init(gabc, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_univar_init(ax, ctx);
    fmpz_mpoly_univar_init(bx, ctx);
    fmpz_mpoly_univar_init(gx, ctx);

    if (polyA->length == 0)
    {
        fmpz_mpoly_set(poly1, polyB, ctx);
        goto done;
    } else if (polyB->length == 0)
    {
        fmpz_mpoly_set(poly1, polyA, ctx);
        goto done;
    }

    if (polyA->length == 1)
    {
        fmpz_mpoly_term_content(bc, polyB, ctx);
        fmpz_mpoly_gcd_monomial(poly1, polyA, bc, ctx);
        goto done;
    } else if (polyB->length == 1)
    {
        fmpz_mpoly_term_content(ac, polyA, ctx);
        fmpz_mpoly_gcd_monomial(poly1, polyB, ac, ctx);
        goto done;
    }

    a_degs = (slong *)TMP_ALLOC(nvars*sizeof(slong));
    b_degs = (slong *)TMP_ALLOC(nvars*sizeof(slong));
    a_leads = (slong *)TMP_ALLOC(nvars*sizeof(slong));
    b_leads = (slong *)TMP_ALLOC(nvars*sizeof(slong));
    fmpz_mpoly_degrees_si(a_degs, polyA, ctx);
    fmpz_mpoly_degrees_si(b_degs, polyB, ctx);

    for (v = 0; v < nvars; v++)
    {
        if (a_degs[v] == 0 && b_degs[v] != 0)
        {
            if (!fmpz_mpoly_to_univar(bx, polyB, v, ctx))
                goto failed;
            fmpz_mpoly_set(poly1, polyA, ctx);
            for (i = 0; i < bx->length; i++)
                _fmpz_mpoly_gcd_prs(poly1, poly1, bx->coeffs + i, ctx);
            goto done;
        }
        if (b_degs[v] == 0 && a_degs[v] != 0)
        {
            if (!fmpz_mpoly_to_univar(ax, polyA, v, ctx))
                goto failed;
            fmpz_mpoly_set(poly1, polyB, ctx);
            for (i = 0; i < ax->length; i++)
                _fmpz_mpoly_gcd_prs(poly1, poly1, ax->coeffs + i, ctx);
            goto done;
        }
    }

    for (v = 0; v < nvars; v++)
    {
        bits = polyA->bits;
        mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        N = mpoly_words_per_exp(bits, ctx->minfo);
        mpoly_gen_offset_shift(&off, &shift, v, N, bits, ctx->minfo);

        if (a_degs[v] != 0)
        {
            a_leads[v] = 0;
            for (i = 0; i < polyA->length; i++)
                a_leads[v] += ((polyA->exps[N*i + off] >> shift) & mask) == a_degs[v];
        }

        bits = polyB->bits;
        mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        N = mpoly_words_per_exp(bits, ctx->minfo);
        mpoly_gen_offset_shift(&off, &shift, v, N, bits, ctx->minfo);

        if (b_degs[v] != 0)
        {
            b_leads[v] = 0;
            for (i = 0; i < polyB->length; i++)
                b_leads[v] += ((polyB->exps[N*i + off] >> shift) & mask) == b_degs[v];
        }
    }

    k = 0;
    m = 1000000;
    for (v = 0; v < nvars; v++)
    {
        if (a_degs[v] != 0)
        {
            if (a_degs[v] >= b_degs[v])
            {
                c = FLINT_MAX(FLINT_BIT_COUNT(b_leads[v] - 1)*a_degs[v], 1);
                c *= b_degs[v];
            } else
            {
                c = FLINT_MAX(FLINT_BIT_COUNT(a_leads[v] - 1)*b_degs[v], 1);
                c *= a_degs[v];
            }
            if (c < m)
            {
                m = c;
                k = v;
            }
        }
    }
    var = k;

    /* check for univariate case */
    d = 0;
    for (v = 0; v < nvars; v++)
    {
        if (v != k)
            d |= a_degs[v] | b_degs[v];
    }
    if (d == 0)
    {
        slong e1, e2, e3;
        fmpz_poly_t u1, u2, u3;
        fmpz_poly_init(u1);
        fmpz_poly_init(u2);
        fmpz_poly_init(u3);
        fmpz_mpoly_to_fmpz_poly(u2, &e2, polyA, k, ctx);
        fmpz_mpoly_to_fmpz_poly(u3, &e3, polyB, k, ctx);
        fmpz_poly_gcd(u1, u2, u3);
        e1 = FLINT_MIN(e2, e3);
        fmpz_mpoly_from_fmpz_poly(poly1, u1, e1, k, ctx);
        fmpz_poly_clear(u1);
        fmpz_poly_clear(u2);
        fmpz_poly_clear(u3);
        goto done;
    }

    if (!fmpz_mpoly_to_univar(ax, polyA, var, ctx))
        goto failed;
    if (!fmpz_mpoly_to_univar(bx, polyB, var, ctx))
        goto failed;

    gx->var = var;

    /* divide out content in ax and bx */
    fmpz_mpoly_set(ac, ax->coeffs + 0, ctx);
    for (i = 1; i < ax->length; i++)
    {
        _fmpz_mpoly_gcd_prs(ac, ac, ax->coeffs + i, ctx);
    }
    fmpz_mpoly_set(bc, bx->coeffs + 0, ctx);
    for (i = 1; i < bx->length; i++)
    {
        _fmpz_mpoly_gcd_prs(bc, bc, bx->coeffs + i, ctx);
    }
    for (i = 0; i < ax->length; i++)
    {
        fmpz_mpoly_divides_monagan_pearce(ax->coeffs + i, ax->coeffs + i, ac, ctx);
    }
    for (i = 0; i < bx->length; i++)
    {
        fmpz_mpoly_divides_monagan_pearce(bx->coeffs + i, bx->coeffs + i, bc, ctx);
    }

    /* compute pseudo gcd of contentless polynomials */
    if (ax->exps[0] == 0 || bx->exps[0] == 0)
    {
        fmpz_mpoly_univar_fit_length(gx, 1, ctx);
        fmpz_mpoly_set_ui(gx->coeffs + 0, WORD(1), ctx);
        gx->exps[0] = 0;
        gx->length = 1;
    } else {
        if (ax->exps[0] >= bx->exps[0])
        {
            _fmpz_mpoly_univar_pgcd(gx, ax, bx, ctx);
        } else {
            _fmpz_mpoly_univar_pgcd(gx, bx, ax, ctx);
        }
        assert(gx->length > 0);
    }

    /* try to divide out easy content from gcd */
    assert(gx->length > 0);
    if ((gx->coeffs + 0)->length != 1
        && (gx->coeffs + gx->length - 1)->length != 1)
    {
        if ((ax->coeffs + 0)->length == 1 || (bx->coeffs + 0)->length == 1)
        {
            fmpz_mpoly_term_content(gc, gx->coeffs + 0, ctx);
            fmpz_mpoly_divides_monagan_pearce(gabc, gx->coeffs + 0, gc, ctx);
            for (i = 0; i < gx->length; i++)
            {
                if (!fmpz_mpoly_divides_monagan_pearce(gx->coeffs + i, gx->coeffs + i, gabc, ctx))
                    assert(0 && "not lead div");
            }
        } else if ((ax->coeffs + ax->length - 1)->length == 1
                || (bx->coeffs + bx->length - 1)->length == 1)
        {
            fmpz_mpoly_term_content(gc, gx->coeffs + gx->length-1, ctx);
            fmpz_mpoly_divides_monagan_pearce(gabc, gx->coeffs + gx->length-1, gc, ctx);
            for (i = 0; i < gx->length; i++)
            {
                if (!fmpz_mpoly_divides_monagan_pearce(gx->coeffs + i, gx->coeffs + i, gabc, ctx))
                    assert(0 && "not trail div");
            }
        } else
        {
            _fmpz_mpoly_gcd_prs(gc, ax->coeffs + 0, bx->coeffs + 0, ctx);
            if (gc->length == 1)
            {
                fmpz_mpoly_term_content(gc, gx->coeffs + 0, ctx);
                fmpz_mpoly_divides_monagan_pearce(gabc, gx->coeffs + 0, gc, ctx);
                for (i = 0; i < gx->length; i++)
                {
                    if (!fmpz_mpoly_divides_monagan_pearce(gx->coeffs + i, gx->coeffs + i, gabc, ctx))
                        assert(0 && "not lead div");
                }                
            }
        }
    }

    /* divide out monomial content from gcd */
    fmpz_mpoly_term_content(gc, gx->coeffs + 0, ctx);
    for (i = 1; i < gx->length; i++)
    {
        fmpz_mpoly_term_content(gabc, gx->coeffs + i, ctx);
        fmpz_mpoly_gcd_monomial(gc, gc, gabc, ctx);
    }
    for (i = 0; i < gx->length; i++)
    {
        fmpz_mpoly_divides_monagan_pearce(gx->coeffs + i, gx->coeffs + i, gc, ctx);
    }

    /* divide out content from gcd */
    fmpz_mpoly_set(gc, gx->coeffs + 0, ctx);
    for (i = 1; i < gx->length; i++)
    {
        _fmpz_mpoly_gcd_prs(gc, gc, gx->coeffs + i, ctx);
    }
    for (i = 0; i < gx->length; i++)
    {
        fmpz_mpoly_divides_monagan_pearce(gx->coeffs + i, gx->coeffs + i, gc, ctx);
    }

    /* put back real content */
    fmpz_mpoly_from_univar(g, gx, ctx);
    _fmpz_mpoly_gcd_prs(gabc, ac, bc, ctx);
    fmpz_mpoly_mul_johnson(poly1, g, gabc, ctx);

done:

    TMP_END;

    fmpz_mpoly_clear(ac, ctx);
    fmpz_mpoly_clear(bc, ctx);
    fmpz_mpoly_clear(gc, ctx);
    fmpz_mpoly_clear(gabc, ctx);
    fmpz_mpoly_clear(g, ctx);

    fmpz_mpoly_univar_clear(ax, ctx);
    fmpz_mpoly_univar_clear(bx, ctx);
    fmpz_mpoly_univar_clear(gx, ctx);

    return success;

failed:
    success = 0;
    fmpz_mpoly_zero(poly1, ctx);
    goto done;
}


int fmpz_mpoly_gcd_prs(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = _fmpz_mpoly_gcd_prs(poly1, poly2, poly3, ctx);

    if (success && (poly1->length > 0) && (fmpz_sgn(poly1->coeffs + 0) < 0))
        fmpz_mpoly_neg(poly1, poly1, ctx);

    return success;
}


