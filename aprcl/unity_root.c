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

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"
#include "fmpz_poly.h"

void unity_init(unity_root element, ulong n)
{
    element->power = n;
    fmpz_poly_init2(element->poly, n);
    element->poly->length = n;
}

void unity_clear(unity_root element)
{
    fmpz_poly_clear(element->poly);
}

void unity_print(unity_root element)
{
    fmpz_poly_print(element->poly); flint_printf("\n");
}

void unity_nth_root(unity_root element, ulong n)
{
    fmpz_poly_set_coeff_ui(element->poly, n, 1);
}

void unity_nth_root_inc(unity_root element, ulong n)
{
    fmpz_add_ui(element->poly->coeffs + n, element->poly->coeffs + n, 1);
}

void unity_nth_root_dec(unity_root element, ulong n)
{
    fmpz_sub_ui(element->poly->coeffs + n, element->poly->coeffs + n, 1);
}

void unity_roots_add(unity_root res, const unity_root element1, const unity_root element2)
{
    res->power = element1->power;
    fmpz_poly_add(res->poly, element1->poly, element2->poly);
}

void unity_roots_mul(unity_root res, const unity_root element1, const unity_root element2)
{
    fmpz_poly_t temp, divider;
    fmpz_poly_init(temp);
    fmpz_poly_mul_classical(temp, element1->poly, element2->poly);
    fmpz_poly_init(divider);
    fmpz_poly_set_coeff_si(divider, 0, -1);
    fmpz_poly_set_coeff_ui(divider, element1->power, 1);
    fmpz_poly_rem(res->poly, temp, divider);
    fmpz_poly_clear(temp);
    fmpz_poly_clear(divider);
    res->power = element1->power;
}

void unity_roots_mul_sub(unity_root res, const unity_root element1, const unity_root element2)
{
    int i, power;
    fmpz_poly_mul_classical(res->poly, element1->poly, element2->poly);

    power = element1->power;
    for (i = res->poly->length - 1; i > power; i--)
    {
        slong coeff = fmpz_get_si(res->poly->coeffs + i);
        if (coeff != 0)
        {
            slong old_coeff = fmpz_get_si(res->poly->coeffs + i - power);
            old_coeff += coeff;
            fmpz_set_si(res->poly->coeffs + i, 0);
            fmpz_set_si(res->poly->coeffs + i - power, old_coeff);
        }
    }
    _fmpz_poly_normalise(res->poly);
}

void unity_roots_reduce_cyclotomic(unity_root res, ulong p)
{
    int i;
    ulong power;
    fmpz_poly_t divider;
    fmpz_poly_init(divider);

    power = res->power / p;
    for (i = 0; i < p; i++)
        fmpz_poly_set_coeff_ui(divider, i * power, 1);

    fmpz_poly_rem(res->poly, res->poly, divider);
    fmpz_poly_clear(divider);
}

