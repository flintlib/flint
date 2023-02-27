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

#include <cstdio>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ctools.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "NTL-interface.h"

#define ZZ_SIZE(p) (((slong *) (p))[1])
#define ZZ_DATA(p) ((mp_limb_t *) (((slong *) (p)) + 2))

NTL_CLIENT

static void fmpz_set_limbs(fmpz_t f, mp_srcptr x, mp_size_t limbs)
{
    if (limbs == 0)
        fmpz_zero(f);
    else if (limbs == 1)
        fmpz_set_ui(f, x[0]);
    else
    {
        __mpz_struct *mf = _fmpz_promote(f);

        mpz_import(mf, limbs, -1, sizeof(mp_limb_t), 0, 0, x);
    }
}

void fmpz_set_ZZ(fmpz_t rop, const ZZ& op)
{
    const _ntl_gbigint x = op.rep;

    if (!x) 
        fmpz_zero(rop);
    else
    {
        const mp_size_t lw = op.size();
        const mp_limb_t *xp = ZZ_DATA(x);

        fmpz_set_limbs(rop, xp, lw);

        if (op < WORD(0))
            fmpz_neg(rop, rop);
    }
}

void fmpz_set_ZZ_p(fmpz_t rop, const ZZ_p& op) 
{
    fmpz_set_ZZ(rop, rep(op));
}

void fmpz_set_zz_p(fmpz_t rop, const zz_p& op) 
{
    fmpz_set_si(rop, rep(op));
}

void fmpz_get_ZZ(ZZ& rop, const fmpz_t op)
{
   mp_limb_t *xp;
   _ntl_gbigint *x = &rop.rep;
   slong lw = fmpz_size(op);
   fmpz c = *op;

   if (lw == 0) 
   {
      if (*x) ZZ_SIZE(*x) = 0;
      return;
   }

   _ntl_gsetlength(x, lw); 
   xp = ZZ_DATA(*x);

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
   
   if (fmpz_sgn(op) < 0) ZZ_SIZE(*x) = -lw;
   else ZZ_SIZE(*x) = lw;
}

void fmpz_get_ZZ_p(ZZ_p& rop, const fmpz_t op) 
{
    ZZ a;
    fmpz_get_ZZ(a, op);
    conv(rop, a);
}

void fmpz_get_zz_p(zz_p& rop, const fmpz_t op) 
{
    conv(rop, fmpz_get_si(op));
}

void fmpz_poly_get_ZZX(ZZX& rop, const fmpz_poly_t op)
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

void fmpz_poly_set_ZZX(fmpz_poly_t rop, const ZZX& op)
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


void fmpz_mod_poly_get_ZZ_pX(ZZ_pX& rop, const fmpz_mod_poly_t op,
                                                      const fmpz_mod_ctx_t ctx)
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

void fmpz_mod_poly_set_ZZ_pX(fmpz_mod_poly_t rop, const ZZ_pX& op,
                                                      const fmpz_mod_ctx_t ctx)
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

void fmpz_mod_poly_get_zz_pX(zz_pX& rop, const fmpz_mod_poly_t op,
                                                      const fmpz_mod_ctx_t ctx)
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

void fmpz_mod_poly_set_zz_pX(fmpz_mod_poly_t rop, const zz_pX& op,
                                                      const fmpz_mod_ctx_t ctx)
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

void fq_get_ZZ_pE(ZZ_pE& rop, const fq_t op, const fq_ctx_t ctx)
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

void fq_set_ZZ_pE(fq_t rop, const ZZ_pE& op, const fq_ctx_t ctx)
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


void fq_poly_get_ZZ_pEX(ZZ_pEX& rop, const fq_poly_t op, const fq_ctx_t ctx)
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

void fq_poly_set_ZZ_pEX(fq_poly_t rop, const ZZ_pEX& op, const fq_ctx_t ctx)
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

/* ----------------------------------------- */
void fq_get_zz_pE(zz_pE& rop, const fq_t op, const fq_ctx_t ctx)
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

void fq_set_zz_pE(fq_t rop, const zz_pE& op, const fq_ctx_t ctx)
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


void fq_poly_get_zz_pEX(zz_pEX& rop, const fq_poly_t op, const fq_ctx_t ctx)
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

void fq_poly_set_zz_pEX(fq_poly_t rop, const zz_pEX& op, const fq_ctx_t ctx)
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


#undef ZZ_SIZE
#undef ZZ_DATA

