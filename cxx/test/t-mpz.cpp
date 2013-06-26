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
#include <stdexcept>
#include <stdio.h>

#include "cxx.h"
#include "cxx/test/helpers.h"

using namespace flint;

void
test_printing()
{
    mpz a(31);

    tassert(a.to_string() == "31");
    tassert(a.to_string(2) == "11111");

    std::ostringstream oss;
    oss << a << '\n' << std::oct << a << '\n'
        << std::hex << a << '\n' << std::dec << a;
    tassert(oss.str() == "31\n37\n1f\n31");

    mpz b(-15);
    tassert((a + b).to_string() == "16");
}

void
test_order()
{
    mpz a(0);
    mpz b(1);
    mpz c(0);
    mpz d(-1);

#define TO(xzero, zero, one, mone) \
    tassert(xzero == zero); \
    tassert(xzero != one); \
    tassert(!(xzero == one)); \
    tassert(!(xzero != zero)); \
    tassert(xzero < one); \
    tassert(xzero <= one); \
    tassert(xzero <= zero); \
    tassert(!(xzero < zero)); \
    tassert(!(xzero <= mone)); \
    tassert(xzero > mone); \
    tassert(xzero >= mone); \
    tassert(xzero >= zero); \
    tassert(!(xzero > zero)); \
    tassert(!(xzero >= one));

    TO(a, c, b, d);
    TO(a, 0, 1, -1);
    TO(a, 0u, 1u, -1);
    TO(a, (signed short)(0), (signed short)(1), (signed short)(-1));
    TO(a, (unsigned short)(0), (unsigned short)(1), -1);
    TO(a, 0l, 1l, -1l);
    TO(0, c, b, d);
    TO(0u, c, b, d);
    TO(0l, c, b, d);
    TO((short)0, c, b, d);
    TO((unsigned short)0, c, b, d);
}

void
test_conversion()
{
    mpz a(4);
    tassert(a.to<slong>() == 4);
    tassert(a.to<ulong>() == 4);
    tassert(a.to<double>() == 4.0);
    // NB: to_string is tested in test_printing
}

void
test_initialisation_assignment()
{
    mpz a(4), b(4l), c(4u), d("4");
    mpz e(a);
    mpz f, g, h, i;
    f = 4;
    g = 4l;
    h = 4u;
    i = "4";
    tassert(a == b && a == c&& a == d && a == e && a == f && a == g && a == h
            && a == i);

    // test deep copying of (f)mpz with more than one digit
    a = "100000000000000000000";
    b = a;
    tassert(a._data()[0] != b._data()[0]);
    mpz j(a);
    tassert(a._data()[0] != j._data()[0]);
    tassert(a == b && a == j);

    // just here to test our assumptions on data format
    tassert(c._data()[0] == d._data()[0]);
}

void
test_arithmetic()
{
#define TAC(seven, three) \
    tassert(seven + three == 10); \
    tassert(seven * three == 21);

#define TA(seven, three) \
    TAC(seven, three); \
    tassert(seven - three == 4); \
    tassert(seven / three == 2); \
    tassert(seven % three == 1)

    TA(mpz(7), mpz(3));
    TA(mpz(7), 3u);
    TAC(7ul, mpz(3));

    // test signed builtins (only div and mul)
    tassert(-7 * mpz(3) == -21);
    tassert(mpz(7) * (-3l) == -21);
    tassert(mpz(21) / -3 == -7);

    // test composite arithmetic
    mpz a(3), b(7);
    tassert(3*(a + b) - (b + (a - 4u)) + ((-(a - b)) % (b / 2)) == 25);

    // test unary minus
    tassert(-a == -3);
}

int
main()
{
    std::cout << "mpz....";

    test_printing();
    test_order();
    test_conversion();
    test_initialisation_assignment();
    test_arithmetic();

    // TODO test counts of allocated temporaries

    std::cout << "PASS" << std::endl;
}
