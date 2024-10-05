/*
    Copyright 1991, 1993-1996, 2000, 2001, 2005 Free Software Foundation, Inc.
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "fmpz.h"

void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        fmpz_set_si(f, d>>FLINT_MIN(exp, SMALL_FMPZ_BITCOUNT_MAX));
    }
    else  /*g is large */
    {
        mpz_ptr mf = _fmpz_promote(f);  /* g is already large */
        mpz_fdiv_q_2exp(mf, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}

void
fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
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
        {
            fmpz q = c1 / c2;       /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            if (r && (c2 ^ r) < WORD(0))
                --q;

            fmpz_set_si(f, q);
        }
        else                    /* h is large and g is small */
        {
            if ((c1 > WORD(0) && fmpz_sgn(h) < 0) || (c1 < WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
                fmpz_set_si(f, WORD(-1));   /* quotient is negative, round down to minus one */
            else
                fmpz_zero(f);
        }
    }
    else                        /* g is large */
    {
        mpz_ptr mf;

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            mf = _fmpz_promote(f);

            if (c2 > 0)         /* h > 0 */
            {
                flint_mpz_fdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_cdiv_q_ui(mf, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_fdiv_q_newton(f, g, h);
            }
            else
            {
                mf = _fmpz_promote(f);
                mpz_fdiv_q(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
                _fmpz_demote_val(f);  /* division by h may result in small value */
            }
        }
    }
}

void
fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
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
        {
            fmpz q = c1 / c2;   /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            if ((c2 > WORD(0) && r < WORD(0)) || (c2 < WORD(0) && r > WORD(0)))
            {
                q--;            /* q cannot overflow as remainder implies |c2| != 1 */
                r += c2;
            }

            fmpz_set_si(f, q);
            fmpz_set_si(s, r);
        }
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
    else                        /* g is large */
    {
        mpz_ptr mf, ms;

		if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
            ms = _fmpz_promote(s);
            mf  = COEFF_TO_PTR(*f);

            if (c2 > 0)         /* h > 0 */
            {
                flint_mpz_fdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_cdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
            _fmpz_demote_val(s);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_fdiv_qr_newton(f, s, g, h);
            }
            else
            {
                _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
                ms = _fmpz_promote(s);
                mf  = COEFF_TO_PTR(*f);

                mpz_fdiv_qr(mf, ms, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));

                _fmpz_demote_val(f);    /* division by h may result in small value */
                _fmpz_demote_val(s);    /* division by h may result in small value */
            }
        }
    }
}

/* these functions were adapted from similar functions in an old version of GMP */
static void _mpz_tdiv_qr_preinvn(mpz_ptr q, mpz_ptr r,
                          mpz_srcptr a, mpz_srcptr d, const fmpz_preinvn_t inv)
{
    slong size1 = a->_mp_size, size2 = d->_mp_size;
    ulong usize1 = FLINT_ABS(size1);
    ulong usize2 = FLINT_ABS(size2);
    ulong qsize = usize1 - usize2 + 1;
    int nm = (inv->norm != 0);
    TMP_INIT;

    nn_ptr qp, rp, ap, dp, tp, sp;

    rp = FLINT_MPZ_REALLOC(r, usize1 + nm);

    if (usize1 < usize2) /* special case preinv code can't deal with */
    {
        mpz_set(r, a); /* remainder equals numerator */
        q->_mp_size = 0; /* quotient is zero */

        return;
    }

    dp = d->_mp_d;
    ap = a->_mp_d;
    qp = FLINT_MPZ_REALLOC(q, qsize + nm);

    TMP_START;
    if ((r == d || q == d) && !nm) /* we have alias with d */
    {
        tp = TMP_ALLOC(usize2*sizeof(ulong));
        mpn_copyi(tp, dp, usize2);
        dp = tp;
    }

    if (r == a || q == a) /* we have alias with a */
    {
        tp = TMP_ALLOC(usize1*sizeof(ulong));
        mpn_copyi(tp, ap, usize1);
        ap = tp;
    }

    /*
       TODO: speedup flint_mpn_divrem_preinvn so we can remove this first
       case here
    */
    slong rn = size2;
    slong qn = size1 - size2 + 1;

    if (rn <= 3 || (rn >= 15 && qn == 1) || (rn >= 60 && qn <= 40))
        mpn_tdiv_qr(qp, rp, 0, ap, usize1, dp, usize2);
    else {
        if (nm) {
            tp = TMP_ALLOC(usize2*sizeof(ulong));
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

static void _mpz_fdiv_qr_preinvn(mpz_ptr q, mpz_ptr r,
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
        t->_mp_d = TMP_ALLOC(usize2*sizeof(ulong));
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
        mpz_ptr mf, ms;

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

void
fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
{
    fmpz c1 = *g;
    slong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_fdiv_q_si). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz q = c1 / c2;       /* compute C quotient */
        fmpz r = c1 - c2 * q;   /* compute remainder */

        if (r && ((c1 ^ c2) < 0))
            --q;

        fmpz_set_si(f, q);
    }
    else                        /* g is large */
    {
        mpz_ptr mf = _fmpz_promote(f);

        if (c2 > 0)
        {
            flint_mpz_fdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        }
        else
        {
            flint_mpz_cdiv_q_ui(mf, COEFF_TO_PTR(c1), -(ulong) c2);
            mpz_neg(mf, mf);
        }
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}

void
fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_fdiv_q_ui). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 > 0)
        {
            fmpz_set_ui(f, c1 / c2);
        }
        else
        {
            ulong q = ((ulong) -c1) / c2;
            ulong r = ((ulong) -c1) - c2 * q;

            if (r)
                ++q;

            fmpz_set_si(f, - (slong) q);
        }
    }
    else                        /* g is large */
    {
        mpz_ptr mf = _fmpz_promote(f);

        flint_mpz_fdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}

void fmpz_fdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        if (d >= 0)
        {
            fmpz_set_ui(f, exp < (SMALL_FMPZ_BITCOUNT_MAX) ? d & ((WORD(1) << exp) - 1) : d);
        }
        else
        {
            if (exp <= SMALL_FMPZ_BITCOUNT_MAX)
            {
                fmpz_set_ui(f, d & ((WORD(1) << exp) - 1));
            }
            else
            {
                mpz_ptr mf = _fmpz_promote(f);

                flint_mpz_set_ui(mf, 1);
                mpz_mul_2exp(mf, mf, exp);
                flint_mpz_sub_ui(mf, mf, -d);
            }
        }
    }
    else  /*g is large */
    {
        mpz_ptr mf = _fmpz_promote(f);  /* g is already large */
        mpz_fdiv_r_2exp(mf, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}

void
fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_fdiv_r). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz q = c1 / c2;   /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            if ((c2 > WORD(0) && r < WORD(0)) || (c2 < WORD(0) && r > WORD(0)))
                r += c2;

            fmpz_set_si(f, r);
        }
        else                    /* h is large and g is small */
        {
            if (c1 == WORD(0))
            {
                fmpz_set_si(f, c1);
            }
            else if ((c1 < WORD(0) && fmpz_sgn(h) < 0) || (c1 > WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
            {
                fmpz_set_si(f, c1);
            }
            else
            {
                fmpz_add(f, g, h);
            }
        }
    }
    else                        /* g is large */
    {
        mpz_ptr mf;

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            /* todo: should not create an mpz */
            mf = _fmpz_promote(f);

            if (c2 > 0)         /* h > 0 */
            {
                flint_mpz_fdiv_r_ui(mf, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_cdiv_r_ui(mf, COEFF_TO_PTR(c1), -c2);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_fdiv_r_newton(f, g, h);
            }
            else
            {
                mf = _fmpz_promote(f);
                mpz_fdiv_r(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
                _fmpz_demote_val(f);    /* division by h may result in small value */
            }
        }
    }
}

ulong
fmpz_fdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong r;

    if (h == UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_fdiv_ui). Division by 0.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 < WORD(0))
        {
            r = h - (-c1 % h);  /* C doesn't correctly handle negative mods */
            if (r == h)
                r = 0;
        }
        else
            r = c1 % h;

        return r;
    }
    else                        /* g is large */
    {
        return flint_mpz_fdiv_ui(COEFF_TO_PTR(c1), h);
    }
}
