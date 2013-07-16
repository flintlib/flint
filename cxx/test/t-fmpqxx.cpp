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

#include "cxx/fmpqxx.h"
#include "cxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpqxx a, b(fmpzxx(6), fmpzxx(10)), c((unsigned short)6, 10u);

    tassert(a.num() == 0 && a.den() == 1);
    tassert(b.num() == 3 && b.den() == 5);
    tassert(c.num() == 3 && c.den() == 5);

    a.num() = -2;
    a.den() = 4;
    tassert(!a.is_canonical());
    a.canonicalise();
    tassert(a.num() == -1 && a.den() == 2);
}

void
test_assignment()
{
    fmpqxx a;
    fmpqxx b(fmpzxx("100000000000000000000"), fmpzxx("100000000000000000001"));
    fmpqxx c(b);
    a = b;
    tassert(a == b && c == b);
    tassert(a.num()._fmpz()[0] != b.num()._fmpz()[0]);
    tassert(a.den()._fmpz()[0] != b.den()._fmpz()[0]);
    tassert(c.num()._fmpz()[0] != b.num()._fmpz()[0]);
    tassert(c.den()._fmpz()[0] != b.den()._fmpz()[0]);

    a = 7;
    tassert(a.num() == 7 && a.den() == 1);
    a = (unsigned short)8;
    tassert(a.num() == 8 && a.den() == 1);
}

void
test_conversion()
{
    fmpqxx a(3, 5u);
    tassert(a.to_string() == "3/5");
    tassert(a.to_string(5) == "3/10");

    std::ostringstream oss;
    oss << a;
    tassert(oss.str() == "3/5");
}

void test_order()
{
    fmpqxx a(0, 1u), b(1, 1u);
    tassert(a < b);
    tassert(a <= a);
    tassert(b > a);
    tassert(b >= b);
    tassert(a == a);
    tassert(a != b);
}

void
test_arithmetic()
{
    fmpqxx a(3, 5u), b(2, 7u);
    fmpzxx c(2);

    tassert(a + b == fmpqxx(3*7 + 2*5, 35u));
    tassert(a * b == fmpqxx(6, 35u));
    tassert(a - b == fmpqxx(3*7 - 5*2, 35u));
    tassert(a / b == fmpqxx(3*7, 10u));

    tassert((a+a) * (c+c) == fmpqxx(24, 5u));
    tassert(c * a == fmpqxx(6, 5u));
    tassert(a / c == fmpqxx(3, 10u));

    tassert(-a == fmpqxx(-3, 5u));
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
    fmpqxx a(3, 5u), b(2, 7u), c(3, 1u);
    tassert(is_ternary((a+a) - b*c));
    tassert(is_ternary(b*c + (a+a)));
    tassert((a+a) - b*c == fmpqxx(6*7 - 5*6, 5u*7u));
    tassert(b*c + (a+a) == fmpqxx(6*7 + 5*6, 5u*7u));

    fmpqxx_ref ar(a);
    fmpqxx_srcref asr(a);
    tassert(a == ar && ar == asr);
    ar = 3;
    tassert(a == c && asr == c);
}

int
main()
{
    std::cout << "fmpqxx....";

    test_init();
    test_assignment();
    test_conversion();
    test_order();
    test_arithmetic();
    test_extras();

    std::cout << "PASS" << std::endl;
    return 0;
}

