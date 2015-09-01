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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#ifndef D_VEC_H
#define D_VEC_H

#ifdef D_VEC_INLINES_C
#define D_VEC_INLINE FLINT_DLL
#else
#define D_VEC_INLINE static __inline__
#endif

#include <math.h>
#include "double_extras.h"
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Memory management  *******************************************************/

FLINT_DLL double * _d_vec_init(slong len);

FLINT_DLL void _d_vec_clear(double * vec);

/*  Randomisation  ***********************************************************/

FLINT_DLL void _d_vec_randtest(double * f, flint_rand_t state, 
                        slong len, slong minexp, slong maxexp);

/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void _d_vec_set(double * vec1, const double * vec2, slong len2);

FLINT_DLL void _d_vec_zero(double * vec, slong len);

/*  Comparison  **************************************************************/

FLINT_DLL int _d_vec_equal(const double * vec1, const double * vec2, slong len);

FLINT_DLL int _d_vec_approx_equal(const double * vec1, const double * vec2, slong len, double eps);

FLINT_DLL int _d_vec_is_zero(const double * vec, slong len);

FLINT_DLL int _d_vec_is_approx_zero(const double * vec, slong len, double eps);

/*  Addition  ****************************************************************/

FLINT_DLL void _d_vec_add(double * res, const double * vec1, const double * vec2, slong len2);

FLINT_DLL void _d_vec_sub(double * res, const double * vec1, const double * vec2, slong len2);

/*  Dot product and norm  **************************************/

FLINT_DLL double _d_vec_dot(const double * vec1, const double * vec2, slong len2);

FLINT_DLL double _d_vec_norm(const double * vec, slong len);

FLINT_DLL double _d_vec_dot_heuristic(const double * vec1, const double * vec2, slong len2, double * err);

FLINT_DLL double _d_vec_dot_thrice(const double * vec1, const double * vec2, slong len2, double * err);

#ifdef __cplusplus
}
#endif

#endif

