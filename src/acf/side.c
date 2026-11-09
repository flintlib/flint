/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz.h"
#include "acf.h"
#include "impl.h"

/*
    Deciding a rounding exactly: a value V near a rounding boundary beta
    (a representable number, or a midpoint for ARF_RND_NEAR) rounds like
    any number on the same side of beta and closer to it than the next
    boundary, so that the side sgn(V - beta) decides it.  For the
    algebraic functions the side is the sign of a polynomial in the exact
    inputs and beta, a sum of exact terms, which _acf_sign_sum finds
    without forming the sum when its terms lie far apart.
*/

void
_acf_terms_init(_acf_terms_t T)
{
    slong i;
    T->n = 0;
    for (i = 0; i < ACF_TERMS_MAX; i++)
        arf_init(T->t + i);
}

void
_acf_terms_clear(_acf_terms_t T)
{
    slong i;
    for (i = 0; i < ACF_TERMS_MAX; i++)
        arf_clear(T->t + i);
}

/* T += c x1 x2 ... xk (exactly; x NULL-terminated, c a small integer) */
void
_acf_terms_add_monomial(_acf_terms_t T, slong c, arf_srcptr x1, arf_srcptr x2, arf_srcptr x3, arf_srcptr x4)
{
    arf_ptr t;

    if (T->n >= ACF_TERMS_MAX)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    t = T->t + T->n;
    arf_set_si(t, c);
    if (x1 != NULL) arf_mul(t, t, x1, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (x2 != NULL) arf_mul(t, t, x2, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (x3 != NULL) arf_mul(t, t, x3, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (x4 != NULL) arf_mul(t, t, x4, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (!arf_is_zero(t))
        T->n++;
}

/* R = c X Y (all products of the terms; R distinct from X and Y) */
void
_acf_terms_mul(_acf_terms_t R, const _acf_terms_t X, const _acf_terms_t Y, slong c)
{
    slong i, j;

    R->n = 0;
    for (i = 0; i < X->n; i++)
        for (j = 0; j < Y->n; j++)
            _acf_terms_add_monomial(R, c, X->t + i, Y->t + j, NULL, NULL);
}

/* R += X */
void
_acf_terms_append(_acf_terms_t R, const _acf_terms_t X)
{
    slong i;

    for (i = 0; i < X->n; i++)
        _acf_terms_add_monomial(R, 1, X->t + i, NULL, NULL, NULL);
}

static int
_acf_cmp_exp_desc(const void * x, const void * y)
{
    arf_srcptr u = *(arf_srcptr const *) x, v = *(arf_srcptr const *) y;
    return -fmpz_cmp(ARF_EXPREF(u), ARF_EXPREF(v));
}

/* The sign of the sum of the n exact terms, of any exponents: the terms
   in decreasing magnitude are summed exactly until the remaining ones,
   each below 2^E(t_i), cannot reach half the partial sum S (r of them
   sum to less than 2^(E(t_i) + bits(r)) <= 2^(E(S) - 2) < |S| / 2), so
   that no exact sum spans the gap to a negligible term. */
int
_acf_sign_sum(const _acf_terms_t T)
{
    arf_srcptr * order;
    arf_t S;
    fmpz_t e;
    slong i, n = T->n;
    int sign;

    if (n == 0)
        return 0;

    order = flint_malloc(n * sizeof(arf_srcptr));
    for (i = 0; i < n; i++)
        order[i] = T->t + i;
    qsort(order, n, sizeof(arf_srcptr), _acf_cmp_exp_desc);

    arf_init(S);
    fmpz_init(e);

    for (i = 0; i < n; i++)
    {
        if (!arf_is_zero(S))
        {
            /* E(t_i) + bits(n - i) + 2 <= E(S): the rest is negligible */
            fmpz_add_ui(e, ARF_EXPREF(order[i]), FLINT_BIT_COUNT(n - i) + 2);
            if (fmpz_cmp(e, ARF_EXPREF(S)) <= 0)
                break;
        }
        arf_add(S, S, order[i], ARF_PREC_EXACT, ARF_RND_DOWN);
    }

    sign = arf_sgn(S);

    arf_clear(S);
    fmpz_clear(e);
    flint_free(order);
    return sign;
}

/* the rounding boundary nearest to an approximation y (nonzero, within
   far less than an ulp at prec + 1 bits of the value): the nearest
   number of prec bits for the directed modes, of prec + 1 bits for
   ARF_RND_NEAR (a midpoint, or a number of prec bits that the rounding
   to nearest keeps either way) */
void
_acf_round_boundary(arf_t beta, const arf_t y, slong prec, arf_rnd_t rnd)
{
    FLINT_ASSERT(!arf_is_special(y));
    arf_set_round(beta, y, (rnd == ARF_RND_NEAR) ? prec + 1 : prec, ARF_RND_NEAR);
}

/* the rounding to prec bits in the strict mode rnd of a value V with
   sgn(V - beta) = side, V closer to beta than to any other number of
   prec + 1 bits: beta itself for side 0, otherwise beta moved by far
   less than that spacing (2^(E - prec - 8), below a quarter of the
   spacing even beneath a power of two), with the inexact flag */
int
_acf_round_side(arf_t res, const arf_t beta, int side, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int inexact;

    FLINT_ASSERT(!arf_is_special(beta));
    FLINT_ASSERT(rnd < ARF_RND_FAST);

    if (side == 0)
        return arf_set_round(res, beta, prec, rnd);

    arf_init(t);
    arf_one(t);
    arf_mul_2exp_fmpz(t, t, ARF_EXPREF(beta));
    arf_mul_2exp_si(t, t, -prec - 8);
    if (side < 0)
        arf_neg(t, t);
    arf_add(t, t, beta, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_set_round(res, t, prec, rnd);
    arf_clear(t);
    inexact = 1;
    return inexact;
}
