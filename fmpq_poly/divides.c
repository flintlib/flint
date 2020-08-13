/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

int _fmpq_poly_divides(fmpz * qpoly, fmpz_t qden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1, 
                    const fmpz * poly2, const fmpz_t den2, slong len2)
{
    fmpz_t c1, c2, n;
    fmpz * poly2_pp, * poly1_pp;
    int divides;

    fmpz_init(c1);
    fmpz_init(c2);
    _fmpz_poly_content(c1, poly1, len1);
    _fmpz_poly_content(c2, poly2, len2);
    
    if (!fmpz_is_one(c1))
    {
	poly1_pp = _fmpz_vec_init(len1);
        _fmpz_vec_scalar_divexact_fmpz(poly1_pp, poly1, len1, c1);
    } else
	poly1_pp = (fmpz *) poly1;

    if (!fmpz_is_one(c2))
    {
        poly2_pp = _fmpz_vec_init(len2);
        _fmpz_vec_scalar_divexact_fmpz(poly2_pp, poly2, len2, c2);
    } else
        poly2_pp = (fmpz *) poly2;    

    divides = _fmpz_poly_divides(qpoly, poly1_pp, len1, poly2_pp, len2);

    if (divides)
    {
        fmpz_init(n);

        _fmpq_mul(n, qden, c1, den1, den2, c2);
        _fmpz_vec_scalar_mul_fmpz(qpoly, qpoly, len1 - len2 + 1, n);

        fmpz_clear(n);
    } else
        fmpz_set_ui(qden, 1);

    if (!fmpz_is_one(c1))
        _fmpz_vec_clear(poly1_pp, len1);

    if (!fmpz_is_one(c2))
        _fmpz_vec_clear(poly2_pp, len2);

    fmpz_clear(c1);
    fmpz_clear(c2);

    return divides;
}

int fmpq_poly_divides(fmpq_poly_t q, const fmpq_poly_t poly1, 
                                    const fmpq_poly_t poly2)
{
    fmpz * qpoly;
    fmpz_t qden;
    slong len1, len2;
    int divides;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len2 == 0)
    {
        if (len1 == 0)
        {
            fmpq_poly_zero(q);
	    return 1;
        } else
	    return 0;
    } 

    if (fmpq_poly_is_zero(poly1))
    {
        fmpq_poly_zero(q);
	return 1;
    }

    if (len2 > len1)
        return 0;

    fmpq_poly_fit_length(q, len1 - len2 + 1);

    if (q == poly1 || q == poly2)
    {
        qpoly = _fmpz_vec_init(len1 - len2 + 1);
	fmpz_init(qden);

	divides = _fmpq_poly_divides(qpoly, qden,
             poly1->coeffs, poly1->den, len1, poly2->coeffs, poly2->den, len2);
        
	_fmpz_vec_set(q->coeffs, qpoly, len1 - len2 + 1);
	fmpz_set(q->den, qden);

	fmpz_clear(qden);
        _fmpz_vec_clear(qpoly, len1 - len2 + 1);
    } else
	divides = _fmpq_poly_divides(q->coeffs, q->den,
            poly1->coeffs, poly1->den, len1, poly2->coeffs, poly2->den, len2);

    _fmpq_poly_set_length(q, len1 - len2 + 1);
    _fmpq_poly_normalise(q);

    return divides;
}

