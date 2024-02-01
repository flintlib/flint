/*
    Copyright (C) 2011, 2020 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"

void
_fmpq_add_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2)
{
    ulong pp, qq, rr, ss;
    mp_limb_t hi, lo;
    int s1, s2;

    if (q1 == q2)
    {
        p1 += p2;

        if (q1 != 1)
        {
            ulong g = n_gcd(FLINT_ABS(p1), q1);

            if (g != 1)
            {
                p1 /= (slong) g;
                q1 /= g;
            }
        }

        fmpz_set_si(rnum, p1);
        _fmpz_demote(rden);
        *rden = q1;
        return;
    }

    if (p1 == 0)
    {
        _fmpz_demote(rnum);
        _fmpz_demote(rden);
        *rnum = p2;
        *rden = q2;
        return;
    }

    if (p2 == 0)
    {
        _fmpz_demote(rnum);
        _fmpz_demote(rden);
        *rnum = p1;
        *rden = q1;
        return;
    }

    qq = q1;
    ss = q2;

    if (p1 >= 0)
    {
        s1 = 0;
        pp = p1;
    }
    else
    {
        s1 = 1;
        pp = -p1;
    }

    if (p2 >= 0)
    {
        s2 = 0;
        rr = p2;
    }
    else
    {
        s2 = 1;
        rr = -p2;
    }

    if (ss == 1)
    {
        umul_ppmm(hi, lo, rr, qq);
        if (s1 == s2)
        {
            add_ssaaaa(hi, lo, hi, lo, 0, pp);
        }
        else
        {
            if (hi || (lo >= pp))
            {
                sub_ddmmss(hi, lo, hi, lo, 0, pp);
                s1 = s2;
            }
            else
            {
                hi = 0;
                lo = pp - lo;
            }
        }

        _fmpz_demote(rden);
        *rden = qq;
    }
    else if (qq == 1)
    {
        umul_ppmm(hi, lo, pp, ss);
        if (s1 == s2)
        {
            add_ssaaaa(hi, lo, hi, lo, 0, rr);
        }
        else
        {
            if (hi || (lo >= rr))
            {
                sub_ddmmss(hi, lo, hi, lo, 0, rr);
            }
            else
            {
                hi = 0;
                lo = rr - lo;
                s1 = s2;
            }
        }

        _fmpz_demote(rden);
        *rden = ss;
    }
    else
    {
        ulong g, h, a, b, t, u, v, denhi, denlo;

        g = n_gcd(qq, ss);

        if (g == 1)
        {
            umul_ppmm(hi, lo, pp, ss);
            umul_ppmm(t, u, qq, rr);
            if (s1 == s2)
            {
                add_ssaaaa(hi, lo, hi, lo, t, u);
            }
            else
            {
                if (hi > t || (hi == t && lo >= u))
                {
                    sub_ddmmss(hi, lo, hi, lo, t, u);
                }
                else
                {
                    sub_ddmmss(hi, lo, t, u, hi, lo);
                    s1 = s2;
                }
            }

            umul_ppmm(denhi, denlo, qq, ss);
        }
        else
        {
            a = qq / g;
            b = ss / g;

            umul_ppmm(hi, lo, pp, b);
            umul_ppmm(t, u, rr, a);
            if (s1 == s2)
            {
                add_ssaaaa(hi, lo, hi, lo, t, u);
            }
            else
            {
                if (hi > t || (hi == t && lo >= u))
                {
                    sub_ddmmss(hi, lo, hi, lo, t, u);
                }
                else
                {
                    sub_ddmmss(hi, lo, t, u, hi, lo);
                    s1 = s2;
                }
            }

            if (hi == 0)
            {
                h = n_gcd(lo, g);
            }
            else
            {
                udiv_qrnnd(t, h, hi % g, lo, g);
                h = n_gcd(h, g);
            }

            if (h != 1)
            {
                if (hi == 0)
                {
                    lo /= h;
                }
                else
                {
                    t = hi / h;
                    u = hi - t * h;
                    udiv_qrnnd(v, u, u, lo, h);
                    hi = t;
                    lo = v;
                }

                qq /= h;
            }

            umul_ppmm(denhi, denlo, qq, b);
        }

        fmpz_set_uiui(rden, denhi, denlo);
    }

    if (s1)
        fmpz_neg_uiui(rnum, hi, lo);
    else
        fmpz_set_uiui(rnum, hi, lo);
}

void
_fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            const fmpz_t r, const fmpz_t s)
{
    fmpz_t g, a, b, t, u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && !COEFF_IS_MPZ(*r) && !COEFF_IS_MPZ(*s))
    {
        _fmpq_add_small(rnum, rden, *p, *q, *r, *s);
        return;
    }

    /* Same denominator */
    if (fmpz_equal(q, s))
    {
        fmpz_add(rnum, p, r);

        /* Both are integers */
        if (fmpz_is_one(q))
        {
            fmpz_set(rden, q);
        }
        else
        {
            fmpz_init(g);
            fmpz_gcd(g, rnum, q);

            if (fmpz_is_one(g))
            {
                fmpz_set(rden, q);
            }
            else
            {
                fmpz_divexact(rnum, rnum, g);
                fmpz_divexact(rden, q, g);
            }
            fmpz_clear(g);
        }
        return;
    }

    /* p/q is an integer */
    if (fmpz_is_one(q))
    {
        fmpz_init(t);
        fmpz_mul(t, p, s);
        fmpz_add(rnum, t, r);
        fmpz_set(rden, s);
        fmpz_clear(t);
        return;
    }

    /* r/s is an integer */
    if (fmpz_is_one(s))
    {
        fmpz_init(t);
        fmpz_mul(t, r, q);
        fmpz_add(rnum, t, p);
        fmpz_set(rden, q);
        fmpz_clear(t);
        return;
    }

    /*
    We want to compute p/q + r/s where the inputs are already
    in canonical form.

    If q and s are coprime, then (p*s + q*r, q*s) is in canonical form.

    Otherwise, let g = gcd(q, s) with q = g*a, s = g*b. Then the sum
    is given by ((p*b + r*a) / (a*b)) / g.

    As above, (p*b + r*a) / (a*b) is in canonical form, and g has
    no common factor with a*b. Thus we only need to reduce (p*b + r*a, g).
    If the gcd is 1, the reduced denominator is g*a*b = q*b.
    */
    fmpz_init(g);
    fmpz_gcd(g, q, s);

    if (fmpz_is_one(g))
    {
        fmpz_init(t);
        fmpz_init(u);

        fmpz_mul(t, p, s);
        fmpz_mul(u, q, r);
        fmpz_add(rnum, t, u);
        fmpz_mul(rden, q, s);

        fmpz_clear(t);
        fmpz_clear(u);
    }
    else
    {
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(t);
        fmpz_init(u);

        fmpz_divexact(a, q, g);
        fmpz_divexact(b, s, g);

        fmpz_mul(t, p, b);
        fmpz_mul(u, r, a);
        fmpz_add(rnum, t, u);

        fmpz_gcd(t, rnum, g);

        if (fmpz_is_one(t))
        {
            fmpz_mul(rden, q, b);
        }
        else
        {
            fmpz_divexact(rnum, rnum, t);
            fmpz_divexact(g, q, t);
            fmpz_mul(rden, g, b);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(t);
        fmpz_clear(u);
    }

    fmpz_clear(g);
}

void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_add(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1),
              fmpq_numref(op2), fmpq_denref(op2));
}

void
_fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            const fmpz_t r)
{
    fmpz_t u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && !COEFF_IS_MPZ(*r))
    {
        _fmpq_add_small(rnum, rden, *p, *q, *r, 1);
        return;
    }

    /* both are integers */
    if (fmpz_is_one(q))
    {
        fmpz_add(rnum, p, r);

        fmpz_set(rden, q);

        return;
    }

    /*
    We want to compute p/q + r/1 where the inputs are already
    in canonical form.

    Note (p + q*r, q) is in canonical form.

    */

    fmpz_init(u);

    fmpz_mul(u, q, r);
    fmpz_add(rnum, p, u);
    fmpz_set(rden, q);

    fmpz_clear(u);
}

void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
{
    _fmpq_add_fmpz(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}

void
_fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            slong r)
{
    fmpz_t u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && r >= COEFF_MIN && r <= COEFF_MAX)
    {
        _fmpq_add_small(rnum, rden, *p, *q, r, 1);
        return;
    }

    /* both are integers */
    if (fmpz_is_one(q))
    {
        if (r >= 0)
           fmpz_add_ui(rnum, p, r);
        else
           fmpz_sub_ui(rnum, p, -r);

        fmpz_set(rden, q);

        return;
    }

    /*
    We want to compute p/q + r/1 where the inputs are already
    in canonical form.

    Note (p + q*r, q) is in canonical form.

    */

    fmpz_init(u);

    fmpz_mul_si(u, q, r);
    fmpz_add(rnum, p, u);
    fmpz_set(rden, q);

    fmpz_clear(u);
}

void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c)
{
    _fmpq_add_si(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}

void
_fmpq_add_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            ulong r)
{
    fmpz_t u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && r <= COEFF_MAX)
    {
        _fmpq_add_small(rnum, rden, *p, *q, r, 1);
        return;
    }

    /* both are integers */
    if (fmpz_is_one(q))
    {
        fmpz_add_ui(rnum, p, r);
        fmpz_set(rden, q);
        return;
    }

    /*
    We want to compute p/q + r/1 where the inputs are already
    in canonical form.

    Note (p + q*r, q) is in canonical form.

    */

    fmpz_init(u);

    fmpz_mul_ui(u, q, r);
    fmpz_add(rnum, p, u);
    fmpz_set(rden, q);

    fmpz_clear(u);
}

void fmpq_add_ui(fmpq_t res, const fmpq_t op1, ulong c)
{
    _fmpq_add_ui(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}
