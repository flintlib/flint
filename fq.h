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

    Copyright (C) 2011 Sebastian Pancratz, 2012 Andres Goens
 
******************************************************************************/


#ifndef FQ_H
#define FQ_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "padic.h"
#include "padic_poly.h"
#include "qadic.h"

/* Data types and context ****************************************************/

/*
#define fq_t padic_poly_t
#define fq_struct padic_poly_struct
*/

typedef qadic_t fq_t;
typedef qadic_struct fq_struct;

void fq_ctx_init_conway(qadic_ctx_t ctx, 
                           const fmpz_t p, long d, const char *var, 
			enum padic_print_mode mode); 

void fq_ctx_clear(fq_ctx_t ctx);

static __inline__ long fq_ctx_degree(const fq_ctx_t ctx)
{
    return ctx->j[ctx->len - 1];
}

static __inline__ void 
fq_ctx_print(const qadic_ctx_t ctx)
{
    printf("p = "), fmpz_print((&ctx->pctx)->p), printf("\n");
    printf("d = %ld\n", ctx->j[ctx->len - 1]);
}

