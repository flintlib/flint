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

    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#ifdef T

#ifdef __cplusplus
extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef struct
{
    TEMPLATE(T, poly_struct) *poly;
    slong *exp;
    slong num;
    slong alloc;
} TEMPLATE(T, poly_factor_struct);

typedef TEMPLATE(T, poly_factor_struct) TEMPLATE(T, poly_factor_t)[1];


FLINT_DLL void TEMPLATE(T, poly_factor_init)(TEMPLATE(T, poly_factor_t) fac,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_clear)(TEMPLATE(T, poly_factor_t) fac,
                               const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_realloc)(TEMPLATE(T, poly_factor_t) fac,
                                 slong alloc,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_fit_length)(TEMPLATE(T, poly_factor_t) fac,
                                    slong len,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_set)(TEMPLATE(T, poly_factor_t) res,
                             const TEMPLATE(T, poly_factor_t) fac,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_insert)(TEMPLATE(T, poly_factor_t) fac,
                                const TEMPLATE(T, poly_t) poly,
                                slong exp, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_print)(const TEMPLATE(T, poly_factor_t) fac,
                               const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_print_pretty)(const TEMPLATE(T, poly_factor_t) fac,
                                      const char * var,
                                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_concat)(TEMPLATE(T, poly_factor_t) res,
                                     const TEMPLATE(T, poly_factor_t) fac,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_pow)(TEMPLATE(T, poly_factor_t) fac, slong exp,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int _TEMPLATE(T, poly_is_squarefree)(const TEMPLATE(T, struct) * f, slong len,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_is_squarefree)(const TEMPLATE(T, poly_t) f,
                                const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_squarefree)(TEMPLATE(T, poly_factor_t) res,
                                    const TEMPLATE(T, poly_t) f,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_is_irreducible)(const TEMPLATE(T, poly_t) f,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_is_irreducible_ddf)(const TEMPLATE(T, poly_t) f,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_is_irreducible_ben_or)(const TEMPLATE(T, poly_t) f,
                                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_distinct_deg)(TEMPLATE(T, poly_factor_t) res,
                                      const TEMPLATE(T, poly_t) poly,
                                      slong * const *degs,
                                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_factor_equal_deg_prob)(TEMPLATE(T, poly_t) factor,
                                        flint_rand_t state,
                                        const TEMPLATE(T, poly_t) pol, slong d,
                                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_equal_deg)(TEMPLATE(T, poly_factor_t) factors,
                                   const TEMPLATE(T, poly_t) pol,
                                   slong d, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_cantor_zassenhaus)(TEMPLATE(T, poly_factor_t) res,
                                           const TEMPLATE(T, poly_t) f,
                                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_kaltofen_shoup)(TEMPLATE(T, poly_factor_t) res,
                                        const TEMPLATE(T, poly_t) poly,
                                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_berlekamp)(TEMPLATE(T, poly_factor_t) factors,
                                   const TEMPLATE(T, poly_t) f,
                                   const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_with_berlekamp)(TEMPLATE(T, poly_factor_t) result,
                                        TEMPLATE(T, t) leading_coeff,
                                        const TEMPLATE(T, poly_t) input,
                                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_with_cantor_zassenhaus)(
    TEMPLATE(T, poly_factor_t) result,
    TEMPLATE(T, t) leading_coeff,
    const TEMPLATE(T, poly_t) input,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor_with_kaltofen_shoup)(TEMPLATE(T, poly_factor_t) result,
                                             TEMPLATE(T, t) leading_coeff,
                                             const TEMPLATE(T, poly_t) input,
                                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_factor)(TEMPLATE(T, poly_factor_t) result,
                         TEMPLATE(T, t) leading_coeff,
                         const TEMPLATE(T, poly_t) input,
                         const TEMPLATE(T, ctx_t) ctx);


FLINT_DLL void TEMPLATE(T, poly_iterated_frobenius_preinv)(TEMPLATE(T, poly_t)* rop,
                                            slong n,
                                            const TEMPLATE(T, poly_t) v,
                                            const TEMPLATE(T, poly_t) vinv,
                                            const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
}
#endif

#endif
