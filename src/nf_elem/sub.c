/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf_elem.h"

void
_nf_elem_sub_lf(nf_elem_t a, const nf_elem_t b,
                const nf_elem_t c, const nf_t nf, int can)
{
    const fmpz *const p = LNF_ELEM_NUMREF(b);
    const fmpz *const q = LNF_ELEM_DENREF(b);
    const fmpz *const r = LNF_ELEM_NUMREF(c);
    const fmpz *const s = LNF_ELEM_DENREF(c);
    fmpz *const rnum = LNF_ELEM_NUMREF(a);
    fmpz *const rden = LNF_ELEM_DENREF(a);
    fmpz_t t;

    if (can)
        _fmpq_sub(rnum, rden, p, q, r, s);
    else
    {
        /* Same denominator */
        if (fmpz_equal(q, s))
        {
            fmpz_sub(rnum, p, r);
            fmpz_set(rden, q);

            return;
        }

        /* p/q is an integer */
        if (fmpz_is_one(q))
        {
            fmpz_init(t);

            fmpz_mul(t, p, s);
            fmpz_sub(rnum, t, r);
            fmpz_set(rden, s);

            fmpz_clear(t);

            return;
        }

        /* r/s is an integer */
        if (fmpz_is_one(s))
        {
            fmpz_init(t);

            fmpz_mul(t, r, q);
            fmpz_sub(rnum, t, p);
            fmpz_set(rden, q);

            fmpz_clear(t);

            return;
        }

        /*
           We want to compute p/q - r/s which is (p*s - q*r, q*s).
         */

        fmpz_init(t);

        fmpz_mul(t, q, r);
        fmpz_mul(rnum, p, s);
        fmpz_sub(rnum, rnum, t);
        fmpz_mul(rden, q, s);

        fmpz_clear(t);
    }
}

void
_nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b,
                const nf_elem_t c, const nf_t nf, int can)
{
    fmpz_t d;

    const fmpz *const bnum = QNF_ELEM_NUMREF(b);
    const fmpz *const bden = QNF_ELEM_DENREF(b);

    const fmpz *const cnum = QNF_ELEM_NUMREF(c);
    const fmpz *const cden = QNF_ELEM_DENREF(c);

    fmpz *const anum = QNF_ELEM_NUMREF(a);
    fmpz *const aden = QNF_ELEM_DENREF(a);

    fmpz_init(d);
    fmpz_one(d);

    if (fmpz_equal(bden, cden))
    {
        fmpz_sub(anum, bnum, cnum);
        fmpz_sub(anum + 1, bnum + 1, cnum + 1);
        fmpz_sub(anum + 2, bnum + 2, cnum + 2);
        fmpz_set(aden, bden);

        if (can && !fmpz_is_one(aden))
        {
            fmpz_gcd3(d, anum + 0, anum + 1, anum + 2);

            if (!fmpz_is_one(d))
            {
                fmpz_gcd(d, d, aden);

                if (!fmpz_is_one(d))
                {
                    fmpz_divexact(anum, anum, d);
                    fmpz_divexact(anum + 1, anum + 1, d);
                    fmpz_divexact(anum + 2, anum + 2, d);
                    fmpz_divexact(aden, aden, d);
                }
            }
        }

        fmpz_clear(d);

        return;
    }

    if (!fmpz_is_one(bden) && !fmpz_is_one(cden))
        fmpz_gcd(d, bden, cden);

    if (fmpz_is_one(d))
    {
        fmpz_mul(anum, bnum, cden);
        fmpz_mul(anum + 1, bnum + 1, cden);
        fmpz_mul(anum + 2, bnum + 2, cden);
        fmpz_submul(anum, cnum, bden);
        fmpz_submul(anum + 1, cnum + 1, bden);
        fmpz_submul(anum + 2, cnum + 2, bden);
        fmpz_mul(aden, bden, cden);
    }
    else
    {
        fmpz_t bden1;
        fmpz_t cden1;

        fmpz_init(bden1);
        fmpz_init(cden1);

        fmpz_divexact(bden1, bden, d);
        fmpz_divexact(cden1, cden, d);

        fmpz_mul(anum, bnum, cden1);
        fmpz_mul(anum + 1, bnum + 1, cden1);
        fmpz_mul(anum + 2, bnum + 2, cden1);
        fmpz_submul(anum, cnum, bden1);
        fmpz_submul(anum + 1, cnum + 1, bden1);
        fmpz_submul(anum + 2, cnum + 2, bden1);

        if (fmpz_is_zero(anum) && fmpz_is_zero(anum + 1)
            && fmpz_is_zero(anum + 2))
            fmpz_one(aden);
        else
        {
            if (can)
            {
                fmpz_t e;

                fmpz_init(e);
                fmpz_gcd3(e, anum + 0, anum + 1, anum + 2);

                if (!fmpz_is_one(e))
                    fmpz_gcd(e, e, d);

                if (fmpz_is_one(e))
                    fmpz_mul(aden, bden, cden1);
                else
                {
                    fmpz_divexact(anum, anum, e);
                    fmpz_divexact(anum + 1, anum + 1, e);
                    fmpz_divexact(anum + 2, anum + 2, e);
                    fmpz_divexact(bden1, bden, e);
                    fmpz_mul(aden, bden1, cden1);
                }

                fmpz_clear(e);
            }
            else
                fmpz_mul(aden, bden, cden1);
        }

        fmpz_clear(bden1);
        fmpz_clear(cden1);
    }

    fmpz_clear(d);
}

void
nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b,
               const nf_elem_t c, const nf_t nf)
{
    if (a == c)
    {
        nf_elem_t t;

        nf_elem_init(t, nf);

        _nf_elem_sub_qf(t, b, c, nf, 1);
        nf_elem_swap(t, a, nf);

        nf_elem_clear(t, nf);
    }
    else
        _nf_elem_sub_qf(a, b, c, nf, 1);
}

void
_nf_elem_sub(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        _nf_elem_sub_lf(a, b, c, nf, 0);
    else if (nf->flag & NF_QUADRATIC)
        _nf_elem_sub_qf(a, b, c, nf, 0);
    else
        fmpq_poly_sub_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 0);
}

void
nf_elem_sub(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        _nf_elem_sub_lf(a, b, c, nf, 1);
    else if (nf->flag & NF_QUADRATIC)
        nf_elem_sub_qf(a, b, c, nf);
    else
        fmpq_poly_sub_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 1);
}

void
nf_elem_fmpq_sub(nf_elem_t a, const fmpq_t c, const nf_elem_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);
        const fmpz *const den2 = LNF_ELEM_DENREF(b);
        const fmpz *const num2 = LNF_ELEM_NUMREF(b);

        _fmpq_sub(num, den, fmpq_numref(c), fmpq_denref(c), num2, den2);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);
        const fmpz *const den2 = QNF_ELEM_DENREF(b);
        const fmpz *const num2 = QNF_ELEM_NUMREF(b);
        slong len = 2;

        while (len != 0 && fmpz_is_zero(num2 + len - 1))
            len--;

        if (len == 0)
        {
            fmpz_set(num, fmpq_numref(c));
            fmpz_set(den, fmpq_denref(c));
        }
        else if (len == 1)
            _fmpq_sub(num, den, fmpq_numref(c), fmpq_denref(c), num2, den2);
        else
        {
            /* fast path */
            if (fmpz_equal(fmpq_denref(c), den2))
            {
                fmpz_sub(num, fmpq_numref(c), num2);
                fmpz_neg(num + 1, num2 + 1);
                fmpz_set(den, den2);
            }
            else                /* slow path */
            {
                fmpz_t d1, d2, g;

                fmpz_init(d1);
                fmpz_init(d2);
                fmpz_init(g);

                nf_elem_set(a, b, nf);

                fmpz_gcd(g, fmpq_denref(c), den);
                fmpz_divexact(d1, fmpq_denref(c), g);
                fmpz_divexact(d2, den, g);

                fmpz_mul(num + 1, num + 1, d1);
                fmpz_mul(num, num, d1);
                fmpz_mul(den, den, d1);

                fmpz_submul(num, d2, fmpq_numref(c));
                fmpz_neg(num, num);
                fmpz_neg(num + 1, num + 1);

                fmpz_clear(g);
                fmpz_clear(d1);
                fmpz_clear(d2);
            }

            _fmpq_poly_canonicalise(num, den, 2);
        }
    }
    else
    {
        fmpq_poly_fmpq_sub(NF_ELEM(a), c, NF_ELEM(b));
    }
}

void
nf_elem_sub_fmpq(nf_elem_t a, const nf_elem_t b, const fmpq_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);
        const fmpz *const den2 = LNF_ELEM_DENREF(b);
        const fmpz *const num2 = LNF_ELEM_NUMREF(b);

        _fmpq_sub(num, den, num2, den2, fmpq_numref(c), fmpq_denref(c));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);
        const fmpz *const den2 = QNF_ELEM_DENREF(b);
        const fmpz *const num2 = QNF_ELEM_NUMREF(b);
        slong len = 2;

        nf_elem_set(a, b, nf);
        while (len != 0 && fmpz_is_zero(num2 + len - 1))
            len--;

        if (len == 0)
        {
            fmpz_neg(num, fmpq_numref(c));
            fmpz_set(den, fmpq_denref(c));
        }
        else if (len == 1)
            _fmpq_sub(num, den, num2, den2, fmpq_numref(c), fmpq_denref(c));
        else
        {
            /* fast path */
            if (fmpz_equal(fmpq_denref(c), den2))
            {
                fmpz_sub(num, num2, fmpq_numref(c));
                fmpz_set(den, den2);
            }
            else                /* slow path */
            {
                fmpz_t d1, d2, g;

                fmpz_init(d1);
                fmpz_init(d2);
                fmpz_init(g);


                fmpz_gcd(g, fmpq_denref(c), den);
                fmpz_divexact(d1, fmpq_denref(c), g);
                fmpz_divexact(d2, den, g);

                fmpz_mul(num + 1, num + 1, d1);
                fmpz_mul(num, num, d1);
                fmpz_mul(den, den, d1);

                fmpz_submul(num, d2, fmpq_numref(c));

                fmpz_clear(g);
                fmpz_clear(d1);
                fmpz_clear(d2);
            }

            _fmpq_poly_canonicalise(num, den, 2);
        }
    }
    else
    {
        fmpq_poly_sub_fmpq(NF_ELEM(a), NF_ELEM(b), c);
    }
}

void
nf_elem_fmpz_sub(nf_elem_t a, const fmpz_t c, const nf_elem_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);
        const fmpz *const den2 = LNF_ELEM_DENREF(b);
        const fmpz *const num2 = LNF_ELEM_NUMREF(b);

        _fmpq_sub_fmpz(num, den, num2, den2, c);
        fmpz_neg(num, num);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);

        nf_elem_neg(a, b, nf);
        fmpz_addmul(num, den, c);
        _fmpq_poly_canonicalise(num, den, 2);
    }
    else
    {
        fmpq_poly_fmpz_sub(NF_ELEM(a), c, NF_ELEM(b));
    }
}

void
nf_elem_sub_fmpz(nf_elem_t a, const nf_elem_t b, const fmpz_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);
        const fmpz *const den2 = LNF_ELEM_DENREF(b);
        const fmpz *const num2 = LNF_ELEM_NUMREF(b);

        _fmpq_sub_fmpz(num, den, num2, den2, c);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);
        slong len = 2;

        nf_elem_set(a, b, nf);

        while (len != 0 && fmpz_is_zero(num + len - 1))
            len--;

        fmpz_submul(num, den, c);
        _fmpq_poly_canonicalise(num, den, len);
    }
    else
    {
        fmpq_poly_sub_fmpz(NF_ELEM(a), NF_ELEM(b), c);
    }
}

void
nf_elem_sub_si(nf_elem_t a, const nf_elem_t b, slong c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);

        nf_elem_set(a, b, nf);

        if (c >= 0)
            fmpz_submul_ui(num, den, c);
        else
            fmpz_addmul_ui(num, den, -c);
        _fmpq_canonicalise(num, den);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);
        slong len = 2;

        nf_elem_set(a, b, nf);

        while (len != 0 && fmpz_is_zero(num + len - 1))
            len--;

        if (c >= 0)
            fmpz_submul_ui(num, den, c);
        else
            fmpz_addmul_ui(num, den, -c);
        _fmpq_poly_canonicalise(num, den, len);
    }
    else
    {
        fmpq_poly_sub_si(NF_ELEM(a), NF_ELEM(b), c);
    }
}

void
nf_elem_si_sub(nf_elem_t a, slong c, const nf_elem_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);

        nf_elem_neg(a, b, nf);

        if (c >= 0)
            fmpz_addmul_ui(num, den, c);
        else
            fmpz_submul_ui(num, den, -c);
        _fmpq_canonicalise(num, den);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);

        nf_elem_neg(a, b, nf);

        if (c >= 0)
            fmpz_addmul_ui(num, den, c);
        else
            fmpz_submul_ui(num, den, -c);

        _fmpq_poly_canonicalise(num, den, 2);
    }
    else
    {
        fmpq_poly_si_sub(NF_ELEM(a), c, NF_ELEM(b));
    }
}
