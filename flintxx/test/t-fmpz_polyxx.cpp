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

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

#include <iostream>
#include <sstream>
#include <string>

#include "fmpz_polyxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpz_polyxx p(10);
    tassert(p.length() == 0);
}

void
test_manipulation()
{
    fmpz_polyxx p, q;
    p.set_coeff(5, 17);
    tassert(p.degree() == 5);
    q.set_coeff(5, 17u);
    p.set_coeff(0, fmpzxx(1));
    q.get_coeff(0) = 1;
    tassert(p == q);

    tassert(p.lead() == 17);
    tassert(p.length() == 6);

    p.zero_coeffs(0, 6);
    tassert(p.is_zero() && !p.is_one() && !p.is_unit());
    p.set_coeff(0, 1);
    tassert(p.is_one() && p.is_unit());
    p.set_coeff(0, -1);
    tassert(p.is_unit());
}

void
test_assignment()
{
    // TODO
}

void
test_conversion()
{
    // TODO
}

void
test_arithmetic()
{
    // TODO
}

// Won't compile if the expression is not done using addmul
template<class T>
bool is_ternary(const T&)
{
    return T::ev_traits_t::temp_rule_t::TERNARY_OP_MARKER + 1;
}

// test stuff which we should get automatically - addmul, references etc
void
test_extras()
{
    // TODO
}

void
test_functions()
{
    // TODO
}

int
main()
{
    std::cout << "fmpz_polyxx....";

    test_init();
    test_manipulation();
    test_assignment();
    test_conversion();
    test_arithmetic();
    test_functions();
    test_extras();

    std::cout << "PASS" << std::endl;
    return 0;
}

