/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_compose_divconquer) (
    TEMPLATE(T, struct) * rop,
    const TEMPLATE(T, struct) * op1, slong len1,
    const TEMPLATE(T, struct) * op2, slong len2,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, k, n;
    slong *hlen, alloc, powlen;
    TEMPLATE(T, struct) * v, **h, *pow, *temp;

    if (len1 <= 2 || len2 <= 1)
    {
        if (len1 == 1)
            TEMPLATE(T, set) (rop, op1, ctx);
        else if (len2 == 1)
            _TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (rop, op1, len1, op2,
                                                      ctx);
        else                    /* len1 == 2 */
            _TEMPLATE(T, poly_compose_horner) (rop, op1, len1, op2, len2, ctx);
        return;
    }

    /* Initialisation */

    hlen = (slong *) flint_malloc(((len1 + 1) / 2) * sizeof(slong));

    k = FLINT_CLOG2(len1) - 1;

    hlen[0] = hlen[1] = ((1 << k) - 1) * (len2 - 1) + 1;
    for (i = k - 1; i > 0; i--)
    {
        slong hi = (len1 + (1 << i) - 1) / (1 << i);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = ((1 << i) - 1) * (len2 - 1) + 1;
    }
    powlen = (1 << k) * (len2 - 1) + 1;

    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    v = _TEMPLATE(T, vec_init) (alloc + 2 * powlen, ctx);
    h = (TEMPLATE(T, struct) **) flint_malloc(((len1 + 1) / 2) *
                                              sizeof(TEMPLATE(T, struct) *));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
    {
        h[i + 1] = h[i] + hlen[i];
        hlen[i] = 0;
    }
    hlen[(len1 - 1) / 2] = 0;
    pow = v + alloc;
    temp = pow + powlen;

    /* Let's start the actual work */

    for (i = 0, j = 0; i < len1 / 2; i++, j += 2)
    {
        if (!TEMPLATE(T, is_zero) (op1 + (j + 1), ctx))
        {
            _TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (h[i], op2, len2,
                                                        op1 + j + 1, ctx);
            TEMPLATE(T, add) (h[i], h[i], op1 + j, ctx);
            hlen[i] = len2;
        }
        else if (!TEMPLATE(T, is_zero) (op1 + j, ctx))
        {
            TEMPLATE(T, set) (h[i], op1 + j, ctx);
            hlen[i] = 1;
        }
    }
    if ((len1 & WORD(1)))
    {
        if (!TEMPLATE(T, is_zero) (op1 + j, ctx))
        {
            TEMPLATE(T, set) (h[i], op1 + j, ctx);
            hlen[i] = 1;
        }
    }

    _TEMPLATE(T, poly_sqr) (pow, op2, len2, ctx);
    powlen = 2 * len2 - 1;

    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            slong templen = powlen + hlen[1] - 1;
            _TEMPLATE(T, poly_mul) (temp, pow, powlen, h[1], hlen[1], ctx);
            _TEMPLATE(T, poly_add) (h[0], temp, templen, h[0], hlen[0], ctx);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }

        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2 * i + 1] > 0)
            {
                _TEMPLATE(T, poly_mul) (h[i], pow, powlen, h[2 * i + 1],
                                        hlen[2 * i + 1], ctx);
                hlen[i] = hlen[2 * i + 1] + powlen - 1;
            }
            else
                hlen[i] = 0;
            _TEMPLATE(T, poly_add) (h[i], h[i], hlen[i], h[2 * i], hlen[2 * i],
                                    ctx);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2 * i]);
        }
        if ((n & WORD(1)))
        {
            _TEMPLATE(T, poly_set) (h[i], h[2 * i], hlen[2 * i], ctx);
            hlen[i] = hlen[2 * i];
        }

        _TEMPLATE(T, poly_sqr) (temp, pow, powlen, ctx);
        powlen += powlen - 1;
        {
            TEMPLATE(T, struct) * t = pow;
            pow = temp;
            temp = t;
        }
    }

    _TEMPLATE(T, poly_mul) (rop, pow, powlen, h[1], hlen[1], ctx);
    _TEMPLATE(T, poly_add) (rop, rop, hlen[0], h[0], hlen[0], ctx);

    _TEMPLATE(T, vec_clear) (v, alloc + 2 * powlen, ctx);
    flint_free(h);
    flint_free(hlen);
}

void
TEMPLATE(T, poly_compose_divconquer) (TEMPLATE(T, poly_t) rop,
                                      const TEMPLATE(T, poly_t) op1,
                                      const TEMPLATE(T, poly_t) op2,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong lenr = (len1 - 1) * (len2 - 1) + 1;

    if (len1 == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
    }
    else if (len1 == 1 || len2 == 0)
    {
        TEMPLATE(T, TEMPLATE(poly_set, T)) (rop, op1->coeffs + 0, ctx);
    }
    else if (rop != op1 && rop != op2)
    {
        TEMPLATE(T, poly_fit_length) (rop, lenr, ctx);
        _TEMPLATE(T, poly_compose_divconquer) (rop->coeffs, op1->coeffs, len1,
                                               op2->coeffs, len2, ctx);
        _TEMPLATE(T, poly_set_length) (rop, lenr, ctx);
        _TEMPLATE(T, poly_normalise) (rop, ctx);
    }
    else
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, lenr, ctx);
        _TEMPLATE(T, poly_compose_divconquer) (t->coeffs, op1->coeffs, len1,
                                               op2->coeffs, len2, ctx);
        _TEMPLATE(T, poly_set_length) (t, lenr, ctx);
        _TEMPLATE(T, poly_normalise) (t, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
}


#endif
