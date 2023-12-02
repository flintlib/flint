/*
    Copyright 1991, 1993-1996, 2000, 2001, 2005 Free Software Foundation, Inc.
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

/* these functions were adapted from similar functions in an old version of GMP */
void _mpz_tdiv_qr_preinvn(mpz_ptr q, mpz_ptr r,
                          mpz_srcptr a, mpz_srcptr d, const fmpz_preinvn_t inv)
{
    slong size1 = a->_mp_size, size2 = d->_mp_size;
    ulong usize1 = FLINT_ABS(size1);
    ulong usize2 = FLINT_ABS(size2);
    ulong qsize = usize1 - usize2 + 1;
    int nm = (inv->norm != 0);
    TMP_INIT;

    mp_ptr qp, rp, ap, dp, tp, sp;

    if (r->_mp_alloc < usize1 + nm)
        mpz_realloc2(r, (usize1 + nm)*FLINT_BITS);

    if (usize1 < usize2) /* special case preinv code can't deal with */
    {
        mpz_set(r, a); /* remainder equals numerator */
        q->_mp_size = 0; /* quotient is zero */

        return;
    }

    if (q->_mp_alloc < qsize + nm)
        mpz_realloc2(q, (qsize + nm)*FLINT_BITS);

    dp = d->_mp_d;
    ap = a->_mp_d;
    qp = q->_mp_d;
    rp = r->_mp_d;

    TMP_START;
    if ((r == d || q == d) && !nm) /* we have alias with d */
    {
        tp = TMP_ALLOC(usize2*FLINT_BITS);
        mpn_copyi(tp, dp, usize2);
        dp = tp;
    }

    if (r == a || q == a) /* we have alias with a */
    {
        tp = TMP_ALLOC(usize1*FLINT_BITS);
        mpn_copyi(tp, ap, usize1);
        ap = tp;
    }

    /*
       TODO: speedup mpir's mullow and mulhigh and use in
       flint_mpn_divrem_preinvn so we can remove this first
       case here
    */
    if (usize2 == 2 || (usize2 > 15 && usize2 < 120))
        mpn_tdiv_qr(qp, rp, 0, ap, usize1, dp, usize2);
    else {
        if (nm) {
            tp = TMP_ALLOC(usize2*FLINT_BITS);
            mpn_lshift(tp, dp, usize2, inv->norm);
            dp = tp;

            rp[usize1] = mpn_lshift(rp, ap, usize1, inv->norm);
            if (rp[usize1] != 0) usize1++, qsize++;
            sp = rp;
        } else
            sp = ap;

        qp[qsize - 1] = flint_mpn_divrem_preinvn(qp, rp, sp, usize1, dp, usize2, inv->dinv);

        if (nm)
            mpn_rshift(rp, rp, usize2, inv->norm);
    }

    qsize -= (qp[qsize - 1] == 0);
    MPN_NORM(rp, usize2);

    q->_mp_size = ((size1 ^ size2) < 0 ? -qsize : qsize);
    r->_mp_size = (size1 < 0 ? -usize2 : usize2);

    TMP_END;
}

void _mpz_fdiv_qr_preinvn(mpz_ptr q, mpz_ptr r,
                          mpz_srcptr a, mpz_srcptr d, const fmpz_preinvn_t inv)
{
    slong size1 = a->_mp_size;
    slong size2 = d->_mp_size;
    ulong usize2 = FLINT_ABS(size2);
    mpz_t t;
    TMP_INIT;

    TMP_START;
    if (q == d || r == d) /* we need d later, so make sure it doesn't alias */
    {
        t->_mp_d = TMP_ALLOC(usize2*FLINT_BITS);
        t->_mp_size = d->_mp_size;
        t->_mp_alloc = d->_mp_alloc;
        mpn_copyi(t->_mp_d, d->_mp_d, usize2);
        d = t;
    }

    _mpz_tdiv_qr_preinvn(q, r, a, d, inv);

    if ((size1 ^ size2) < 0 && r->_mp_size != 0)
    {
        flint_mpz_sub_ui(q, q, 1);
        mpz_add(r, r, d);
    }

    TMP_END;
}

void
fmpz_fdiv_qr_preinvn(fmpz_t f, fmpz_t s, const fmpz_t g,
                         const fmpz_t h, const fmpz_preinvn_t inv)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_fdiv_q). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
            fmpz_fdiv_qr(f, s, g, h);
        else                    /* h is large and g is small */
        {
            if (c1 == WORD(0))
            {
                fmpz_set_ui(f, WORD(0)); /* g is zero */
                fmpz_set_si(s, c1);
            }
            else if ((c1 < WORD(0) && fmpz_sgn(h) < 0) || (c1 > WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
            {
                fmpz_zero(f);   /* quotient is positive, round down to zero */
                fmpz_set_si(s, c1);
            }
            else
            {
                fmpz_add(s, g, h);
                fmpz_set_si(f, WORD(-1));    /* quotient is negative, round down to minus one */
            }
        }
    }
    else /* g is large */
    {
        __mpz_struct * mf, * ms;

        if (!COEFF_IS_MPZ(c2))  /* h is small */
            fmpz_fdiv_qr(f, s, g, h);
        else
        {
            _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
            ms = _fmpz_promote(s);
            mf  = COEFF_TO_PTR(*f);

            _mpz_fdiv_qr_preinvn(mf, ms, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2), inv);

            _fmpz_demote_val(f);    /* division by h may result in small value */
            _fmpz_demote_val(s);    /* division by h may result in small value */
        }
    }
}

/* todo: avoid the temporary variable */
void
fmpz_fdiv_r_preinvn(fmpz_t f, const fmpz_t g,
                         const fmpz_t h, const fmpz_preinvn_t inv)
{
    fmpz_t q;
    fmpz_init(q);
    fmpz_fdiv_qr_preinvn(q, f, g, h, inv);
    fmpz_clear(q);
}
