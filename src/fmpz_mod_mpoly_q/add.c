/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mod_mpoly_q.h"

void
_fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(x_num, ctx))
    {
        if (fmpz_mod_is_one(y_den->coeffs, ctx->ffinfo))
        {
            fmpz_mod_mpoly_set(res_num, y_num, ctx);
            fmpz_mod_mpoly_set(res_den, y_den, ctx);
            return;
        }
        else
        {
            fmpz_t g;
            fmpz_mod_mpoly_t t;

            fmpz_init(g);
            fmpz_mod_mpoly_init(t, ctx);

            fmpz_mod_mpoly_gcd_assert_successful(t, y_num, y_den,ctx);

            fmpz_mod_set_fmpz(g, y_den->coeffs, ctx->ffinfo);
            fmpz_mod_inv(g, g, ctx->ffinfo);

            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, y_num, g, ctx);
            fmpz_mod_mpoly_make_monic(res_den, y_den, ctx);

            fmpz_mod_mpoly_gcd_assert_successful(t, res_num, res_den, ctx);

            if (!fmpz_mod_mpoly_is_one(t, ctx))
            {
                fmpz_mod_mpoly_divexact(res_num, res_num, t, ctx);
                fmpz_mod_mpoly_divexact(res_den, res_den, t, ctx);

            }
            
            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_clear(g);
        }

        return;
    }

    if (fmpz_mod_mpoly_is_zero(y_num, ctx))
    {
        if (fmpz_mod_is_one(x_den->coeffs, ctx->ffinfo))
        {
            fmpz_mod_mpoly_set(res_num, x_num, ctx);
            fmpz_mod_mpoly_set(res_den, x_den, ctx);
            return;
        }
        else
        {
            fmpz_t g;
            fmpz_mod_mpoly_t t;

            fmpz_init(g);
            fmpz_mod_mpoly_init(t, ctx);

            fmpz_mod_mpoly_gcd_assert_successful(t, x_num, x_den,ctx);

            fmpz_mod_set_fmpz(g, x_den->coeffs, ctx->ffinfo);
            fmpz_mod_inv(g, g, ctx->ffinfo);

            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, x_num, g, ctx);

            fmpz_mod_mpoly_make_monic(res_den, x_den, ctx);
            
            if (!fmpz_mod_mpoly_is_one(t, ctx))
            {
                fmpz_mod_mpoly_divexact(res_num, res_num, t, ctx);
                fmpz_mod_mpoly_divexact(res_den, res_den, t, ctx);

            }

            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_clear(g);
        }

        return;
    }

    if (fmpz_mod_mpoly_equal(x_den, y_den, ctx))
    {
        fmpz_mod_mpoly_add(res_num, x_num, y_num, ctx);

        if (fmpz_mod_mpoly_is_one(x_den, ctx) || fmpz_mod_mpoly_is_zero(res_num, ctx))
        {
            fmpz_mod_mpoly_one(res_den, ctx);
            return;
        }
        else
        {
            fmpz_t g;
            fmpz_mod_mpoly_t t;

            fmpz_init(g);
            fmpz_mod_mpoly_init(t, ctx);

            fmpz_mod_set_fmpz(g, x_den->coeffs, ctx->ffinfo);
            fmpz_mod_inv(g, g, ctx->ffinfo);

            fmpz_mod_mpoly_gcd_assert_successful(t, res_num, x_den, ctx);
            fmpz_mod_mpoly_set(res_den, x_den, ctx);
            
            if (fmpz_mod_mpoly_is_one(t, ctx))
            {
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g, ctx);
                fmpz_mod_mpoly_make_monic(res_den, res_den, ctx);
                
                fmpz_mod_mpoly_clear(t, ctx);
                fmpz_clear(g);                
                return;
            }
            else
            {
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g, ctx);
                fmpz_mod_mpoly_make_monic(res_den, res_den, ctx);
                fmpz_mod_mpoly_divides(res_num, res_num, t, ctx);
                fmpz_mod_mpoly_divides(res_den, res_den, t, ctx);

                fmpz_mod_mpoly_clear(t, ctx);
                fmpz_clear(g);                 
                return;
            }
        }
    }

    if (fmpz_mod_mpoly_is_one(x_den, ctx))
    {  
        fmpz_t g;
        fmpz_mod_mpoly_t t, u;

        fmpz_init(g);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);


        fmpz_mod_set_fmpz(g, y_den->coeffs, ctx->ffinfo);
        fmpz_mod_inv(g, g, ctx->ffinfo);
        
        fmpz_mod_mpoly_mul(t, x_num, y_den, ctx);
        fmpz_mod_mpoly_add(res_num, t, y_num, ctx);

        if (fmpz_mod_mpoly_is_zero(res_num,ctx))
        {
            fmpz_mod_mpoly_clear(u, ctx);
            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_clear(g);
            return;
        }

        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g, ctx);
        fmpz_mod_mpoly_make_monic(res_den, y_den, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(u, res_num, res_den, ctx);
        fmpz_mod_mpoly_divexact(res_num, res_num, u, ctx);
        fmpz_mod_mpoly_divexact(res_den, res_den, u, ctx);

        fmpz_mod_mpoly_clear(u, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_clear(g);
        return;
    }

    if (fmpz_mod_mpoly_is_one(y_den, ctx))
    {
        fmpz_t g;
        fmpz_mod_mpoly_t t, u;

        fmpz_init(g);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);

        fmpz_mod_set_fmpz(g, x_den->coeffs, ctx->ffinfo);
        fmpz_mod_inv(g, g, ctx->ffinfo);
        
        fmpz_mod_mpoly_mul(t, y_num, x_den, ctx);
        fmpz_mod_mpoly_add(res_num, t, x_num, ctx);

        if (fmpz_mod_mpoly_is_zero(res_num,ctx))
        {
            fmpz_mod_mpoly_clear(u, ctx);
            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_clear(g);
            return;
        }

        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g, ctx);
        fmpz_mod_mpoly_make_monic(res_den, x_den, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(u, res_num, res_den, ctx);
        fmpz_mod_mpoly_divexact(res_num, res_num, u, ctx);
        fmpz_mod_mpoly_divexact(res_den, res_den, u, ctx);

        fmpz_mod_mpoly_clear(u, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_clear(g);
        return;
    }

    {
        fmpz_mod_mpoly_t g, t, u;

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(g, x_den, y_den, ctx);

        if (fmpz_mod_mpoly_is_one(g, ctx))
        {   
            fmpz_t g_coeff;
            
            fmpz_init(g_coeff);

            fmpz_mod_mpoly_mul(t, x_num, y_den, ctx);
            fmpz_mod_mpoly_mul(u, y_num, x_den, ctx);
            fmpz_mod_mpoly_add(res_num, t, u, ctx);
            fmpz_mod_mpoly_mul(res_den, x_den, y_den, ctx);
                        
            fmpz_mod_set_fmpz(g_coeff,res_den->coeffs,ctx->ffinfo);
            fmpz_mod_inv(g_coeff, g_coeff, ctx->ffinfo);
            
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g_coeff, ctx);
            fmpz_mod_mpoly_make_monic(res_den, res_den, ctx);

            fmpz_clear(g_coeff);

        }
        else
        {
            fmpz_mod_mpoly_t a, b;
            fmpz_t g_coeff;
            
            fmpz_init(g_coeff);

            fmpz_mod_mpoly_init(a, ctx);
            fmpz_mod_mpoly_init(b, ctx);

            _fmpz_mod_mpoly_q_mpoly_divexact(a, x_den, g, ctx);
            _fmpz_mod_mpoly_q_mpoly_divexact(b, y_den, g, ctx);

            fmpz_mod_mpoly_mul(t, x_num, b, ctx);
            fmpz_mod_mpoly_mul(u, y_num, a, ctx);
            fmpz_mod_mpoly_add(res_num, t, u, ctx);

            fmpz_mod_mpoly_gcd_assert_successful(t, res_num, g, ctx);

            if (fmpz_mod_mpoly_is_one(t, ctx))
            {
                fmpz_mod_mpoly_mul(res_den, x_den, b, ctx);
            }
            else
            {
                _fmpz_mod_mpoly_q_mpoly_divexact(res_num, res_num, t, ctx);
                _fmpz_mod_mpoly_q_mpoly_divexact(g, x_den, t, ctx);
                fmpz_mod_mpoly_mul(res_den, g, b, ctx);
            }

            fmpz_mod_set_fmpz(g_coeff,res_den->coeffs,ctx->ffinfo);
            fmpz_mod_inv(g_coeff, g_coeff, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g_coeff, ctx);
            fmpz_mod_mpoly_make_monic(res_den, res_den, ctx);

            fmpz_mod_mpoly_clear(a, ctx);
            fmpz_mod_mpoly_clear(b, ctx);
            fmpz_clear(g_coeff);
        }

        fmpz_mod_mpoly_gcd_assert_successful(g, res_num, res_den, ctx);
        fmpz_mod_mpoly_divexact(res_num, res_num, g, ctx);
        fmpz_mod_mpoly_divexact(res_den, res_den, g, ctx);

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(u, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        return;
    }
}

void
_fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t yy, yy_num, yy_den;

    fmpz_init(yy);
    fmpz_init(yy_num);
    fmpz_init(yy_den);
    fmpz_mod_set_fmpz(yy_num, y_num, ctx->ffinfo);
    fmpz_mod_set_fmpz(yy_den, y_den, ctx->ffinfo);
    
    fmpz_mod_inv(yy_den, yy_den, ctx->ffinfo);
    fmpz_mod_mul(yy, yy_num, yy_den, ctx->ffinfo);

    if (fmpz_mod_mpoly_is_zero(x_num, ctx))
    {
        fmpz_mod_mpoly_set_fmpz(res_num, yy, ctx);
        fmpz_mod_mpoly_one(res_den, ctx);
        
        fmpz_clear(yy_num);
        fmpz_clear(yy_den);
        fmpz_clear(yy);
        return;
    }

    if (fmpz_is_zero(yy))
    {
        fmpz_mod_mpoly_set(res_num, x_num, ctx);
        fmpz_mod_mpoly_set(res_den, x_den, ctx);

        if (!fmpz_mod_mpoly_is_one(res_den, ctx))
        {
            fmpz_t g;
            fmpz_mod_mpoly_t t;

            fmpz_init(g);
            fmpz_mod_mpoly_init(t, ctx);

            fmpz_mod_set_fmpz(g,res_den->coeffs, ctx->ffinfo);
            fmpz_mod_mpoly_gcd_assert_successful(t, res_num, res_den, ctx);

            fmpz_mod_inv(g, g, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g, ctx);
            fmpz_mod_mpoly_divexact(res_den, res_den, t, ctx);

            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_clear(g);
        }

        fmpz_clear(yy_num);
        fmpz_clear(yy_den);
        fmpz_clear(yy);
        return;
    }

    if (fmpz_mod_mpoly_is_one(x_den, ctx))
    {
        fmpz_mod_mpoly_add_fmpz(res_num, x_num, yy, ctx);
        fmpz_mod_mpoly_one(res_den, ctx);

        fmpz_clear(yy_num);
        fmpz_clear(yy_den);
        fmpz_clear(yy);
        return;
    }

    fmpz_t g;
    fmpz_mod_mpoly_t t, u;

    fmpz_init(g);
    fmpz_mod_mpoly_init(t, ctx);
    fmpz_mod_mpoly_init(u, ctx);

    fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(u, x_den, yy, ctx);
    fmpz_mod_mpoly_add(res_num, x_num, u, ctx);

    if (fmpz_mod_mpoly_is_zero(res_num, ctx))
    {
        fmpz_mod_mpoly_one(res_den, ctx);
        
        fmpz_clear(g);
        fmpz_clear(yy);
        fmpz_clear(yy_num);
        fmpz_clear(yy_den);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(u, ctx);
        return;
    }
    
    fmpz_mod_mpoly_set(res_den, x_den, ctx);

    fmpz_mod_set_fmpz(g,res_den->coeffs, ctx->ffinfo);
    fmpz_mod_mpoly_gcd_assert_successful(t, res_num, res_den, ctx);

    fmpz_mod_inv(g, g, ctx->ffinfo);
    fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g, ctx);
    fmpz_mod_mpoly_make_monic(res_den, res_den, ctx);
    fmpz_mod_mpoly_divexact(res_num, res_num, t, ctx);
    fmpz_mod_mpoly_divexact(res_den, res_den, t, ctx);

    fmpz_clear(g);
    fmpz_clear(yy);
    fmpz_clear(yy_num);
    fmpz_clear(yy_den);
    fmpz_mod_mpoly_clear(t, ctx);
    fmpz_mod_mpoly_clear(u, ctx);
}

void
fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                fmpz_mod_mpoly_q_numref(y), fmpz_mod_mpoly_q_denref(y),
                ctx);
}

void
fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                fmpq_numref(y), fmpq_denref(y),
                ctx);
}

void
fmpz_mod_mpoly_q_add_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t one;
    *one = 1;
    _fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                y, one,
                ctx);
}
