/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

