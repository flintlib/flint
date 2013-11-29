/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
 
******************************************************************************/

#ifndef FQ_H
#define FQ_H

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

/* Data types and context ****************************************************/

#ifdef __cplusplus
extern "C" {
#endif

typedef fmpz_poly_t fq_t;
typedef fmpz_poly_struct fq_struct;

typedef struct
{
    fmpz p;

    int sparse_modulus;

    fmpz *a;
    slong *j;
    slong len;

    fmpz_mod_poly_t modulus;
    fmpz_mod_poly_t inv;

    char *var;
}
fq_ctx_struct;

typedef fq_ctx_struct fq_ctx_t[1];

void fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var);

int _fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_ctx_init_modulus(fq_ctx_t ctx,
                         fmpz_mod_poly_t modulus,
                         const char *var);

void fq_ctx_randtest(fq_ctx_t ctx, flint_rand_t state);

void fq_ctx_clear(fq_ctx_t ctx);

static __inline__ slong fq_ctx_degree(const fq_ctx_t ctx)
{
    return ctx->modulus->length - 1;
}

#define fq_ctx_prime(ctx)  (&((ctx)->p))

static __inline__ void fq_ctx_order(fmpz_t f, const fq_ctx_t ctx)
{
    fmpz_set(f, fq_ctx_prime(ctx));
    fmpz_pow_ui(f, f, fq_ctx_degree(ctx));
}

static __inline__ int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx)
{
    int r;

    r = flint_fprintf(file, "p = ");
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, fq_ctx_prime(ctx));
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "\nd = %wd\n", fq_ctx_degree(ctx));
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "f(X) = ");
    if (r <= 0)
        return r;

    r = fmpz_mod_poly_fprint_pretty(file, ctx->modulus, "X");
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "\n");

    return r;
}

static __inline__ void fq_ctx_print(const fq_ctx_t ctx)
{
    fq_ctx_fprint(stdout, ctx);
}


/* Memory managment  *********************************************************/

static __inline__ void fq_init(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_init(rop);
}

static __inline__ void fq_init2(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_init2(rop, fq_ctx_degree(ctx));
}

static __inline__ void fq_clear(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_clear(rop);
}

static __inline__ 
void _fq_sparse_reduce(fmpz *R, slong lenR, const fq_ctx_t ctx)
{
    const slong d = ctx->j[ctx->len - 1];

    FMPZ_VEC_NORM(R, lenR);

    if (lenR > d)
    {
        slong i, k;

        for (i = lenR - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                fmpz_submul(R + ctx->j[k] + i - d, R + i, ctx->a + k);
            }
            fmpz_zero(R + i);
        }
    }

    _fmpz_vec_scalar_mod_fmpz(R, R, FLINT_MIN(d, lenR), fq_ctx_prime(ctx));
}

static __inline__ void _fq_dense_reduce(fmpz* R, slong lenR, const fq_ctx_t ctx)
{
    fmpz  *q, *r;

    if (lenR < ctx->modulus->length)
    {
        _fmpz_vec_scalar_mod_fmpz(R, R, lenR, fq_ctx_prime(ctx));
        return;
    }
    
    q = _fmpz_vec_init(lenR - ctx->modulus->length + 1);
    r = _fmpz_vec_init(ctx->modulus->length - 1);

    _fmpz_mod_poly_divrem_newton_n_preinv(q, r, R, lenR, 
                                        ctx->modulus->coeffs, ctx->modulus->length,
                                        ctx->inv->coeffs, ctx->inv->length,
                                        fq_ctx_prime(ctx));

    _fmpz_vec_set(R, r, ctx->modulus->length - 1);
    _fmpz_vec_clear(q, lenR - ctx->modulus->length + 1);
    _fmpz_vec_clear(r, ctx->modulus->length - 1);

}


static __inline__ void _fq_reduce(fmpz* R, slong lenR, const fq_ctx_t ctx)
{
    if (ctx->sparse_modulus)
        _fq_sparse_reduce(R, lenR, ctx);
    else
        _fq_dense_reduce(R, lenR, ctx);    
}

static __inline__ void fq_reduce(fq_t rop, const fq_ctx_t ctx)
{
    _fq_reduce(rop->coeffs, rop->length, ctx);
    rop->length = FLINT_MIN(rop->length, ctx->modulus->length - 1);
    _fmpz_poly_normalise(rop);
}

/* Basic arithmetic **********************************************************/

void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

void fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void fq_neg(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx);

void fq_mul_si(fq_t rop, const fq_t op, slong x, const fq_ctx_t ctx);

void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx);

void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx);

void fq_inv(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void _fq_pow(fmpz *rop, const fmpz *op, slong len, const fmpz_t e,
             const fq_ctx_t ctx);

void fq_pow(fq_t rop, const fq_t op1, const fmpz_t e, const fq_ctx_t ctx);

void fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx);

void fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

/* Randomisation *************************************************************/

void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

void fq_randtest_dense(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

/* Comparison ****************************************************************/

static __inline__ int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
{
    return fmpz_poly_equal(op1, op2);
}

static __inline__ int fq_is_zero(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_is_zero(op);
}

static __inline__ int fq_is_one(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_is_one(op);
}

/* Assignments and conversions ***********************************************/

static __inline__ void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx)
{
    fmpz_poly_set(rop, op);
}

static __inline__ void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx)
{
    fmpz_poly_set_fmpz(rop, x);
    fq_reduce(rop, ctx);
}

static __inline__ void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx)
{
    fmpz_poly_set_ui(rop, x);
    fq_reduce(rop, ctx);
}

static __inline__ void fq_swap(fq_t op1, fq_t op2, const fq_ctx_t ctx)
{
    fmpz_poly_swap(op1, op2);
}

static __inline__ void fq_zero(fq_t rop,  const fq_ctx_t ctx)
{
    fmpz_poly_zero(rop);
}

static __inline__ void fq_one(fq_t rop,  const fq_ctx_t ctx)
{
    fmpz_poly_one(rop);
}

static __inline__ void fq_gen(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_zero(rop);
    fmpz_poly_set_coeff_ui(rop, 0, 0);
    fmpz_poly_set_coeff_ui(rop, 1, 1);
}

/* Output ********************************************************************/

static __inline__
int fq_fprint(FILE * file, const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_fprint(file, op);
}

static __inline__
void fq_print(const fq_t op, const fq_ctx_t ctx)
{
    fmpz_poly_print(op);
}

static __inline__ 
int fq_fprint_pretty(FILE * file, const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_fprint_pretty(file, op, ctx->var);
}

static __inline__ 
int fq_print_pretty(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_print_pretty(op, ctx->var);
}

char *
fq_get_str(const fq_t op, const fq_ctx_t ctx);

char *
fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx);

/* Special functions *********************************************************/

void _fq_trace(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx);

void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);

void _fq_frobenius(fmpz *rop, const fmpz *op, slong len, slong e, 
                   const fq_ctx_t ctx);

void fq_frobenius(fq_t rop, const fq_t op, slong e, const fq_ctx_t ctx);

void _fq_norm(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx);

void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);


/* Bit packing ******************************************************/

void
fq_bit_pack(fmpz_t f, const fq_t op, mp_bitcnt_t bit_size,
            const fq_ctx_t ctx);

void
fq_bit_unpack(fq_t rop, const fmpz_t f, mp_bitcnt_t bit_size,
              const fq_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

