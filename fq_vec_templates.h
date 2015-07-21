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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#ifdef T

#include "flint.h"
#include "templates.h"
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Memory management  *******************************************************/

FLINT_DLL TEMPLATE(T, struct) * _TEMPLATE(T, vec_init)(slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_clear)(TEMPLATE(T, struct) * vec, slong len,
                             const TEMPLATE(T, ctx_t) ctx);

/*  Randomisation  ***********************************************************/
FLINT_DLL void _TEMPLATE(T, vec_randtest)(TEMPLATE(T, struct) * f, 
                           flint_rand_t state, 
                           slong len, 
                           const TEMPLATE(T, ctx_t) ctx);

/*  Norms  *******************************************************************/

/*  Input and output  ********************************************************/
FLINT_DLL int _TEMPLATE(T, vec_fprint)(FILE * file,
                             const TEMPLATE(T, struct) * vec,
                             slong len,
                             const TEMPLATE(T, ctx_t) ctx);

FQ_VEC_TEMPLATES_INLINE
int _TEMPLATE(T, vec_print)(const TEMPLATE(T, struct) * vec, slong len,
                            const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, vec_fprint)(stdout, vec, len, ctx);
}

/*  Conversions  *************************************************************/

/*  Assignment and basic manipulation  ***************************************/
FLINT_DLL void _TEMPLATE(T, vec_set)(TEMPLATE(T, struct) * v,
                      const TEMPLATE(T, struct) * f, 
                      slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_swap)(TEMPLATE(T, struct) * vec1,
                       TEMPLATE(T, struct) * vec2,
                       slong len2, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_zero)(TEMPLATE(T, struct) * v, 
                       slong len,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_neg)(TEMPLATE(T, struct) * vec1,
                      const TEMPLATE(T, struct) * vec2,
                      slong len2,
                      const TEMPLATE(T, ctx_t) ctx);

/*  Comparison  **************************************************************/
FLINT_DLL int _TEMPLATE(T, vec_is_zero)(const TEMPLATE(T, struct) * vec, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int _TEMPLATE(T, vec_equal)(const TEMPLATE(T, struct) * vec1,
                        const TEMPLATE(T, struct) * vec2, slong len,
                        const TEMPLATE(T, ctx_t) ctx);

/*  Sorting  *****************************************************************/

/*  Addition  ****************************************************************/

FLINT_DLL void _TEMPLATE(T, vec_add)(TEMPLATE(T, struct) * res,
                           const TEMPLATE(T, struct) * vec1, 
                           const TEMPLATE(T, struct) * vec2, 
                           slong len2,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_sub)(TEMPLATE(T, struct) * res,
                           const TEMPLATE(T, struct) * vec1, 
                           const TEMPLATE(T, struct) * vec2, 
                           slong len2,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T))(TEMPLATE(T, struct) * poly1, 
                                             const TEMPLATE(T, struct) * poly2,
                                             slong len2,
                                             const TEMPLATE(T, t) x,
                                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, TEMPLATE(vec_scalar_submul, T))(TEMPLATE(T, struct) * poly1, 
                                             const TEMPLATE(T, struct) * poly2,
                                             slong len2,
                                             const TEMPLATE(T, t) x,
                                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE3(T, vec_scalar_mul, T) (TEMPLATE(T, struct) * poly1,
                                  const TEMPLATE(T, struct) * poly2,
                                  slong len2, const TEMPLATE(T, t) x,
                                  const TEMPLATE(T, ctx_t) ctx);
     
/* ****************************************************************************/
FLINT_DLL void _TEMPLATE(T, vec_dot)(TEMPLATE(T, t) res,
                      const TEMPLATE(T, struct) * vec1,
                      const TEMPLATE(T, struct) * vec2,
                      slong len2,
                      const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
 }
#endif

#endif
