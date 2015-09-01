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
#ifndef FQ_ZECH_POLY_H
#define FQ_ZECH_POLY_H

#ifdef FQ_ZECH_POLY_INLINES_C
#define FQ_ZECH_POLY_INLINE FLINT_DLL
#define FQ_POLY_TEMPLATES_INLINE FLINT_DLL
#else
#define FQ_ZECH_POLY_INLINE static __inline__
#define FQ_POLY_TEMPLATES_INLINE static __inline__
#endif

#include "fq_zech_mat.h"

#define FQ_ZECH_POLY_DIVREM_DIVCONQUER_CUTOFF  16
#define FQ_ZECH_COMPOSE_MOD_LENH_CUTOFF 6
#define FQ_ZECH_COMPOSE_MOD_PREINV_LENH_CUTOFF 6
#define FQ_ZECH_SQR_CLASSICAL_CUTOFF 100
#define FQ_ZECH_MUL_CLASSICAL_CUTOFF 90
#define FQ_ZECH_MULLOW_CLASSICAL_CUTOFF 90

#define FQ_ZECH_POLY_HGCD_CUTOFF 35
#define FQ_ZECH_POLY_GCD_CUTOFF 96
#define FQ_ZECH_POLY_SMALL_GCD_CUTOFF 96

#ifdef T
#undef T
#endif

#define T fq_zech
#define CAP_T FQ_ZECH
#include "fq_poly_templates.h"
#undef CAP_T
#undef T

#include "fq_zech_poly_factor.h"

#endif
