/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "nf_elem.h"

int
_nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        slong d, bits1, bits2;
        int res = 1;

        const fmpz *const anum = LNF_ELEM_NUMREF(a);
        const fmpz *const bnum = LNF_ELEM_NUMREF(b);
        const fmpz *const aden = LNF_ELEM_DENREF(a);
        const fmpz *const bden = LNF_ELEM_DENREF(b);

        fmpz_t t1, t2;

        if (fmpz_equal(aden, bden))
            return fmpz_equal(anum, bnum);

        d = fmpz_bits(aden) - fmpz_bits(bden) + 1;

        bits1 = fmpz_bits(anum);
        bits2 = fmpz_bits(bnum);
        if (!(bits1 == 0 && bits2 == 0) && (ulong) (bits1 - bits2 + d) > 2)
            return 0;

        fmpz_init(t1);
        fmpz_init(t2);

        fmpz_mul(t1, anum, bden);
        fmpz_mul(t2, bnum, aden);

        if (!fmpz_equal(t1, t2))
            res = 0;

        fmpz_clear(t1);
        fmpz_clear(t2);

        return res;
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        slong d, bits1, bits2;
        int res = 1;

        const fmpz *const anum = QNF_ELEM_NUMREF(a);
        const fmpz *const bnum = QNF_ELEM_NUMREF(b);
        const fmpz *const aden = QNF_ELEM_DENREF(a);
        const fmpz *const bden = QNF_ELEM_DENREF(b);

        fmpz_t t1, t2;

        if (fmpz_equal(aden, bden))
            return fmpz_equal(anum, bnum) && fmpz_equal(anum + 1, bnum + 1);

        d = fmpz_bits(aden) - fmpz_bits(bden) + 1;

        bits1 = fmpz_bits(anum + 1);
        bits2 = fmpz_bits(bnum + 1);
        if (!(bits1 == 0 && bits2 == 0) && (ulong) (bits1 - bits2 + d) > 2)
            return 0;

        bits1 = fmpz_bits(anum);
        bits2 = fmpz_bits(bnum);
        if (!(bits1 == 0 && bits2 == 0) && (ulong) (bits1 - bits2 + d) > 2)
            return 0;

        fmpz_init(t1);
        fmpz_init(t2);

        fmpz_mul(t1, anum, bden);
        fmpz_mul(t2, bnum, aden);

        if (!fmpz_equal(t1, t2))
        {
            res = 0;
            goto cleanup;
        }

        fmpz_mul(t1, anum + 1, bden);
        fmpz_mul(t2, bnum + 1, aden);

        if (!fmpz_equal(t1, t2))
        {
            res = 0;
            goto cleanup;
        }

      cleanup:

        fmpz_clear(t1);
        fmpz_clear(t2);

        return res;
    }
    else
    {
        const slong len1 = NF_ELEM(a)->length;
        const slong len2 = NF_ELEM(b)->length;

        if (len1 != len2)
            return 0;

        if (fmpz_equal
            (fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b))))
            return _fmpz_vec_equal(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1);
        else
        {
            slong i;
            slong d = fmpz_bits(fmpq_poly_denref(NF_ELEM(b)))
                - fmpz_bits(fmpq_poly_denref(NF_ELEM(a))) + 1;
            fmpz *p1 = NF_ELEM_NUMREF(a);
            fmpz *p2 = NF_ELEM_NUMREF(b);
            fmpz_t gcd, den1, den2;
            fmpz *t1, *t2;
            int equal;

            for (i = 0; i < len1; i++)
            {
                slong b1 = fmpz_bits(p1 + i);
                slong b2 = fmpz_bits(p2 + i);
                if (!(b1 == 0 && b2 == 0) && (ulong) (b1 - b2 + d) > 2)
                    return 0;
            }

            fmpz_init(gcd);
            fmpz_init(den1);
            fmpz_init(den2);

            /* TODO: possibly only compute GCD if it will save time */
            fmpz_gcd(gcd, fmpq_poly_denref(NF_ELEM(a)),
                     fmpq_poly_denref(NF_ELEM(b)));
            fmpz_divexact(den1, fmpq_poly_denref(NF_ELEM(a)), gcd);
            fmpz_divexact(den2, fmpq_poly_denref(NF_ELEM(b)), gcd);

            t1 = _fmpz_vec_init(len1);
            t2 = _fmpz_vec_init(len1);

            _fmpz_vec_scalar_mul_fmpz(t1, p1, len1, den2);
            _fmpz_vec_scalar_mul_fmpz(t2, p2, len2, den1);

            equal = _fmpz_vec_equal(t1, t2, len1);

            fmpz_clear(gcd);
            fmpz_clear(den1);
            fmpz_clear(den2);

            _fmpz_vec_clear(t1, len1);
            _fmpz_vec_clear(t2, len1);

            return equal;
        }
    }
}

int
nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        if (!fmpz_equal(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b)))
            return 0;

        if (!fmpz_equal(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b)))
            return 0;

        return 1;
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        if (!fmpz_equal(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b)))
            return 0;

        if (!fmpz_equal(QNF_ELEM_NUMREF(a), QNF_ELEM_NUMREF(b)))
            return 0;

        if (!fmpz_equal(QNF_ELEM_NUMREF(a) + 1, QNF_ELEM_NUMREF(b) + 1))
            return 0;

        return 1;
    }
    else
    {
        const slong len1 = NF_ELEM(a)->length;
        const slong len2 = NF_ELEM(b)->length;

        if (len1 != len2)
            return 0;

        if (fmpz_equal
            (fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b))))
            return _fmpz_vec_equal(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1);
        else
            return 0;
    }
}

int
nf_elem_equal_si(const nf_elem_t a, const slong b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_is_one(LNF_ELEM_DENREF(a)) &&
            fmpz_equal_si(LNF_ELEM_NUMREF(a), b);
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
            fmpz_is_one(QNF_ELEM_DENREF(a)) &&
            fmpz_equal_si(QNF_ELEM_NUMREF(a), b);
    else
    {
        if (b == 0)
            return fmpq_poly_is_zero(NF_ELEM(a));
        else
            return NF_ELEM(a)->length == 1 &&
                fmpz_is_one(NF_ELEM_DENREF(a)) &&
                fmpz_equal_si(NF_ELEM_NUMREF(a), b);
    }
}

int
nf_elem_equal_ui(const nf_elem_t a, const ulong b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_is_one(LNF_ELEM_DENREF(a)) &&
            fmpz_equal_ui(LNF_ELEM_NUMREF(a), b);
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
            fmpz_is_one(QNF_ELEM_DENREF(a)) &&
            fmpz_equal_ui(QNF_ELEM_NUMREF(a), b);
    else
    {
        if (b == 0)
            return fmpq_poly_is_zero(NF_ELEM(a));
        else
            return NF_ELEM(a)->length == 1 &&
                fmpz_is_one(NF_ELEM_DENREF(a)) &&
                fmpz_equal_ui(NF_ELEM_NUMREF(a), b);
    }
}

int
nf_elem_equal_fmpz(const nf_elem_t a, const fmpz_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_is_one(LNF_ELEM_DENREF(a))
            && fmpz_equal(LNF_ELEM_NUMREF(a), b);
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
            fmpz_is_one(QNF_ELEM_DENREF(a)) &&
            fmpz_equal(QNF_ELEM_NUMREF(a), b);
    else
    {
        if (NF_ELEM(a)->length == 0)
            return fmpz_is_zero(b);
        else if (NF_ELEM(a)->length == 1)
            return fmpz_is_one(NF_ELEM_DENREF(a)) &&
                fmpz_equal(NF_ELEM_NUMREF(a), b);
        else
            return 0;
    }
}

int
nf_elem_equal_fmpq(const nf_elem_t a, const fmpq_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_equal(LNF_ELEM_NUMREF(a), fmpq_numref(b)) &&
            fmpz_equal(LNF_ELEM_DENREF(a), fmpq_denref(b));
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
            fmpz_equal(QNF_ELEM_NUMREF(a), fmpq_numref(b)) &&
            fmpz_equal(QNF_ELEM_DENREF(a), fmpq_denref(b));
    else
    {
        if (NF_ELEM(a)->length == 0)
            return fmpq_is_zero(b);
        else if (NF_ELEM(a)->length == 1)
            return fmpz_equal(NF_ELEM_NUMREF(a), fmpq_numref(b)) &&
                fmpz_equal(NF_ELEM_DENREF(a), fmpq_denref(b));
        else
            return 0;
    }
}
