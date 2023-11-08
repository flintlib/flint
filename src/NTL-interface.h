/*
    Copyright (C) 2007 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_NTL_INT_H
#define FLINT_NTL_INT_H

#include <NTL/ZZX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEX.h>

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_poly.h"

NTL_CLIENT

/*
   Converts an NTL ZZ to an fmpz_t.

   Assumes the fmpz_t has already been allocated to have sufficient space.
*/
inline void fmpz_set_ZZ(fmpz_t rop, const ZZ& op)
{
    const _ntl_gbigint x = op.rep;

    if (!x)
        fmpz_zero(rop);
    else
    {
        const mp_size_t lw = op.size();
        const mp_limb_t *xp = ((mp_limb_t *) (((slong *) (x)) + 2));

        if (lw == 0)
            fmpz_zero(rop);
        else if (lw == 1)
            fmpz_set_ui(rop, xp[0]);
        else
        {
            __mpz_struct *mf = _fmpz_promote(rop);

            mpz_import(mf, lw, -1, sizeof(mp_limb_t), 0, 0, xp);
        }

        if (op < WORD(0))
            fmpz_neg(rop, rop);
    }
}

/*
   Converts an fmpz_t to an NTL ZZ. Allocation is automatically handled.
 */
inline void fmpz_get_ZZ(NTL_NNS ZZ& rop, const fmpz_t op)
{
   mp_limb_t *xp;
   _ntl_gbigint *x = &rop.rep;
   slong lw = fmpz_size(op);
   fmpz c = *op;

   if (lw == 0)
   {
      if (*x) ((slong *)(*x))[1] = 0;  // size
      return;
   }

   _ntl_gsetlength(x, lw);
   xp = ((mp_limb_t *) (((slong *) (*x)) + 2));  // data

   if (COEFF_IS_MPZ(c))
   {
      __mpz_struct * m = COEFF_TO_PTR(c);
      mpn_copyi(xp, m->_mp_d, lw);
   } else
   {
      if (c < WORD(0))
         xp[0] = -c;
      else
         xp[0] = c;
   }

   if (fmpz_sgn(op) < 0) ((slong *)(*x))[1] = -lw;
   else ((slong *)(*x))[1] = lw;
}

/*
   Converts an NTL ZZ_p to an fmpz_t.

   Assumes the fmpz_t has already been allocated to have sufficient space.
*/
inline void fmpz_set_ZZ_p(fmpz_t rop, const NTL_NNS ZZ_p& op)
{
    fmpz_set_ZZ(rop, rep(op));
}

/*
   Converts an fmpz_t to an NTL ZZ_p. Allocation is automatically handled.
 */
inline void fmpz_get_ZZ_p(NTL_NNS ZZ_p& rop, const fmpz_t op)
{
    ZZ a;
    fmpz_get_ZZ(a, op);
    conv(rop, a);
}

/*
   Converts an NTL zz_p to an fmpz_t.
*/
inline void fmpz_set_zz_p(fmpz_t rop, const NTL_NNS zz_p& op)
{
    fmpz_set_si(rop, rep(op));
}

/*
   Converts an fmpz_t to an NTL zz_p.
 */
void fmpz_get_zz_p(NTL_NNS zz_p& rop, const fmpz_t op)
{
    conv(rop, fmpz_get_si(op));
}

/*
  Converts an fmpz_poly_t to an NTL ZZX.
*/
inline void fmpz_poly_get_ZZX(NTL_NNS ZZX& rop, const fmpz_poly_t op)
{
    const slong len = op->length;

    if (len == 0)
    {
        rop = 0;
    }
    else
    {
        slong i;
        ZZ *ap;

        rop.rep.SetLength(len);

        for (i = 0, ap = rop.rep.elts(); i < len; i++, ap++)
        {
            fmpz_get_ZZ(*ap, op->coeffs + i);
        }
    }
}

/*
  Converts an NTL ZZX to an fmpz_poly_t.
*/
inline void fmpz_poly_set_ZZX(fmpz_poly_t rop, const NTL_NNS ZZX& op)
{
    const slong len = deg(op) + 1;

    if (len == 0)
    {
        fmpz_poly_zero(rop);
    }
    else
    {
        slong i;
        const ZZ *ap;

        fmpz_poly_fit_length(rop, len);
        _fmpz_poly_set_length(rop, len);

        for (i = 0, ap = op.rep.elts(); i < len; i++, ap++)
        {
            fmpz_set_ZZ(rop->coeffs + i, *ap);
        }
    }
}
/*
  Converts an fmpz_mod_poly_t to an NTL ZZ_pX.
*/
inline void fmpz_mod_poly_get_ZZ_pX(NTL_NNS ZZ_pX& rop,
                             const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx)
{
    const slong len = op->length;

    if (len == 0)
    {
        rop = 0;
    }
    else
    {
        slong i;
        ZZ_p *ap;

        rop.rep.SetLength(len);

        for (i = 0, ap = rop.rep.elts(); i < len; i++, ap++)
        {
            fmpz_get_ZZ_p(*ap, op->coeffs + i);
        }
    }
}

/*
  Converts an NTL ZZ_pX to an fmpz_poly_t.
*/
inline void fmpz_mod_poly_set_ZZ_pX(fmpz_mod_poly_t rop,
                              const NTL_NNS ZZ_pX& op, const fmpz_mod_ctx_t ctx)
{
    const slong len = deg(op) + 1;

    if (len == 0)
    {
        fmpz_mod_poly_zero(rop, ctx);
    }
    else
    {
        slong i;
        const ZZ_p *ap;

        fmpz_mod_poly_fit_length(rop, len, ctx);
        _fmpz_mod_poly_set_length(rop, len);

        for (i = 0, ap = op.rep.elts(); i < len; i++, ap++)
        {
            fmpz_set_ZZ_p(rop->coeffs + i, *ap);
        }
    }
}

/*
  Converts an fq_t to an NTL ZZ_pE.
*/
inline void fq_get_ZZ_pE(NTL_NNS ZZ_pE& rop, const fq_t op, const fq_ctx_t ctx)
{
    ZZ_pX p;

    const slong len = op->length;

    if (len == 0)
    {
        rop = 0;
    }
    else
    {
        slong i;
        ZZ_p *ap;

        p.rep.SetLength(len);

        for (i = 0, ap = p.rep.elts(); i < len; i++, ap++)
        {
            fmpz_get_ZZ_p(*ap, op->coeffs + i);
        }
        conv(rop, p);
    }
}

/*
  Converts an NTL ZZ_pE to an fq_t.
*/
inline void fq_set_ZZ_pE(fq_t rop, const NTL_NNS ZZ_pE& op, const fq_ctx_t ctx)
{
    const slong len = deg(rep(op)) + 1;

    if (len == 0)
    {
        fq_zero(rop, ctx);
    }
    else
    {
        slong i;
        const ZZ_p *ap;

        fmpz_poly_fit_length(rop, len);
        _fmpz_poly_set_length(rop, len);

        for (i = 0, ap = rep(op).rep.elts(); i < len; i++, ap++)
        {
            fmpz_set_ZZ_p(rop->coeffs + i, *ap);
        }
        _fmpz_poly_normalise(rop);
    }
}

/*
  Converts an fq_poly_t to an NTL ZZ_pEX.
*/
inline void fq_poly_get_ZZ_pEX(NTL_NNS ZZ_pEX& rop, const fq_poly_t op,
                                                             const fq_ctx_t ctx)
{
    const slong len = op->length;

    if (len == 0)
    {
        rop = 0;
    }
    else
    {
        slong i;
        ZZ_pE *ap;

        rop.rep.SetLength(len);

        for (i = 0, ap = rop.rep.elts(); i < len; i++, ap++)
        {
            fq_get_ZZ_pE(*ap, op->coeffs + i, ctx);
        }
    }
}

/*
  Converts an NTL ZZ_pEX to an fq_poly_t.
*/
inline void fq_poly_set_ZZ_pEX(fq_poly_t rop, const NTL_NNS ZZ_pEX& op,
                                                             const fq_ctx_t ctx)
{
    const slong len = deg(op) + 1;

    if (len == 0)
    {
        fq_poly_zero(rop, ctx);
    }
    else
    {
        slong i;
        const ZZ_pE *ap;

        fq_poly_fit_length(rop, len, ctx);
        _fq_poly_set_length(rop, len, ctx);

        for (i = 0, ap = op.rep.elts(); i < len; i++, ap++)
        {
            fq_set_ZZ_pE(rop->coeffs + i, *ap, ctx);
        }
        _fq_poly_normalise(rop, ctx);
    }
}

/*
  Converts an fmpz_mod_poly_t to an NTL zz_pX.
*/
inline void fmpz_mod_poly_get_zz_pX(NTL_NNS zz_pX& rop,
                             const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx)
{
    const slong len = op->length;

    if (len == 0)
    {
        rop = 0;
    }
    else
    {
        slong i;
        zz_p *ap;

        rop.rep.SetLength(len);

        for (i = 0, ap = rop.rep.elts(); i < len; i++, ap++)
        {
            fmpz_get_zz_p(*ap, op->coeffs + i);
        }
    }
}

/*
  Converts an NTL zz_pX to an fmpz_poly_t.
*/
inline void fmpz_mod_poly_set_zz_pX(fmpz_mod_poly_t rop,
                              const NTL_NNS zz_pX& op, const fmpz_mod_ctx_t ctx)
{
    const slong len = deg(op) + 1;

    if (len == 0)
    {
        fmpz_mod_poly_zero(rop, ctx);
    }
    else
    {
        slong i;
        const zz_p *ap;

        fmpz_mod_poly_fit_length(rop, len, ctx);
        _fmpz_mod_poly_set_length(rop, len);

        for (i = 0, ap = op.rep.elts(); i < len; i++, ap++)
        {
            fmpz_set_zz_p(rop->coeffs + i, *ap);
        }
    }
}

/*
  Converts an fq_t to an NTL zz_pE.
*/
inline void fq_get_zz_pE(NTL_NNS zz_pE& rop, const fq_t op, const fq_ctx_t ctx)
{
    zz_pX p;

    const slong len = op->length;

    if (len == 0)
    {
        rop = 0;
    }
    else
    {
        slong i;
        zz_p *ap;

        p.rep.SetLength(len);

        for (i = 0, ap = p.rep.elts(); i < len; i++, ap++)
        {
            fmpz_get_zz_p(*ap, op->coeffs + i);
        }
        conv(rop, p);
    }
}

/*
  Converts an NTL zz_pE to an fq_t.
*/
inline void fq_set_zz_pE(fq_t rop, const NTL_NNS zz_pE& op, const fq_ctx_t ctx)
{
    const slong len = deg(rep(op)) + 1;

    if (len == 0)
    {
        fq_zero(rop, ctx);
    }
    else
    {
        slong i;
        const zz_p *ap;

        fmpz_poly_fit_length(rop, len);
        _fmpz_poly_set_length(rop, len);

        for (i = 0, ap = rep(op).rep.elts(); i < len; i++, ap++)
        {
            fmpz_set_zz_p(rop->coeffs + i, *ap);
        }
        _fmpz_poly_normalise(rop);
    }
}

/*
  Converts an fq_poly_t to an NTL zz_pEX.
*/
inline void fq_poly_get_zz_pEX(NTL_NNS zz_pEX& rop, const fq_poly_t op,
                                                             const fq_ctx_t ctx)
{
    const slong len = op->length;

    if (len == 0)
    {
        rop = 0;
    }
    else
    {
        slong i;
        zz_pE *ap;

        rop.rep.SetLength(len);

        for (i = 0, ap = rop.rep.elts(); i < len; i++, ap++)
        {
            fq_get_zz_pE(*ap, op->coeffs + i, ctx);
        }
    }
}

/*
  Converts an NTL zz_pEX to an fq_poly_t.
*/
inline void fq_poly_set_zz_pEX(fq_poly_t rop, const NTL_NNS zz_pEX& op,
                                                             const fq_ctx_t ctx)
{
    const slong len = deg(op) + 1;

    if (len == 0)
    {
        fq_poly_zero(rop, ctx);
    }
    else
    {
        slong i;
        const zz_pE *ap;

        fq_poly_fit_length(rop, len, ctx);
        _fq_poly_set_length(rop, len, ctx);

        for (i = 0, ap = op.rep.elts(); i < len; i++, ap++)
        {
            fq_set_zz_pE(rop->coeffs + i, *ap, ctx);
        }
        _fq_poly_normalise(rop, ctx);
    }
}

#endif
