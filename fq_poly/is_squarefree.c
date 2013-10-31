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
#include "ulong_extras.h"
#include "fq_poly.h"

int
_fq_poly_is_squarefree(const fq_struct * f, slong len, const fq_ctx_t ctx)
{
    fq_struct *fd, *g;
    fq_t invfd;
    slong dlen;
    int res;

    if (len <= 2)
        return len != 0;

    fd = _fq_vec_init(2 * (len - 1), ctx);
    g = fd + len - 1;

    _fq_poly_derivative(fd, f, len, ctx);
    dlen = len - 1;
    FQ_VEC_NORM(fd, dlen, ctx);

    if (dlen)
    {
        fq_init(invfd, ctx);
        fq_inv(invfd, fd + (dlen - 1), ctx);
        res = (_fq_poly_gcd(g, f, len, fd, dlen, invfd, ctx) == 1);
        fq_clear(invfd, ctx);
    }
    else
        res = 0;                /* gcd(f, 0) = f, and len(f) > 2 */

    _fq_vec_clear(fd, 2 * (len - 1), ctx);
    return res;
}

int
fq_poly_is_squarefree(const fq_poly_t f, const fq_ctx_t ctx)
{
    return _fq_poly_is_squarefree(f->coeffs, f->length, ctx);
}
