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
    //tassert(a*b + a*a + b*b == 21+9+49); // XXX

    // test unary minus
    tassert(-a == -3);

    // test assignment arithmetic
#define TAA(op, res) \
    { \
        mpz tmp1(10), tmp2(10), tmp3(10); \
        mpz three(3); \
        tmp1 op three; \
        tmp2 op 3u; \
        tmp3 op three*1; \
        tassert(tmp1 == res); \
        tassert(tmp2 == res); \
        tassert(tmp3 == res); \
    }
    TAA(+=, 13);
    TAA(*=, 30);
    TAA(/=, 3);
    TAA(%=, 1);
}

void
test_functions()
{
    mpz a(2);
    mpz b(16);

    tassert(pow(a, 4u) == 16);
    tassert(root(b, 4) == 2);
    tassert(root(b, (unsigned short)4) == 2);
    tassert(sqrt(b) == 4);
    tassert(sqrt(a) == 1);
    tassert(rfac(a, 3u) == 2*3*4);
    tassert(a + fac(4u) == 2 + 4*3*2);
    tassert(bin(4u, 2u) == 6);

    // check namespace std
    tassert(std::pow(a, 4u) == 16);
    tassert(std::sqrt(b) == 4);

    // check a composite expression
    tassert(2u + sqrt(a + a) == 4);

    // check immediate functions
    tassert(divisible(b, a + a));
    tassert(divisible(b, a + a + 2u) == false);
    tassert(divisible(a + a, (unsigned short)2));
    tassert(divisible(-a + b, 3) == false);
}

template<class T>
void assert_is_mpz(const T&)
{
    tassert(traits::is_mpz<T>::val);
}
struct newtype { };
void
test_traits()
{
    tassert(traits::is_mpz<mpz>::val == true);
    tassert(traits::is_mpz<int>::val == false);
    assert_is_mpz(mpz(1) + mpz(2));
    tassert((traits::is_mpz<mpz_expression<
                operations::immediate, newtype> >::val == false));
}

template<class T>
unsigned count_temporaries(const T&)
{
    return T::ev_traits_t::rule_t::temporaries_t::len;
}

template<class T>
unsigned count_temporaries2(const T&)
{
    return T::ev_traits_t::temp_rule_t::temporaries_t::len
         // this term is always zero, but causes compiler error
         // if we are not actually in the ternary case
         + T::ev_traits_t::temp_rule_t::TERNARY_OP_MARKER;
}

void
test_temporaries()
{
    mpz a, b, c;
    tassert(count_temporaries(a + b) == 0);
    tassert(count_temporaries(a + b + c + a + b + c) == 1);
    tassert(count_temporaries(((a / c) + (b % a)) / ((b + c) + (c / a))) == 3);
    tassert(count_temporaries((a/b) + (a/c) + (b/c) + (c/b)) == 2);
    tassert(count_temporaries2((a*b) + (a*c) + (b*c) + (c*b)) == 1);
}

void
test_ternary()
{
    mpz b(2), c(3), d(4);

#define T0 fac(4u)
#define T1 (b + b + b)
#define T2 (T1 + T1)
#define T3 (T2 + T2)
#define T4 (T3 + T3)

#define TT3(m1, m2, m3, ntemps) \
    tassert(count_temporaries2(m1 + m2*m3) == ntemps); \
    tassert(b + (m1 + m2*m3) == 2 + m1.to<long>() + m2.to<long>()*m3.to<long>()); \
    tassert(count_temporaries2(m1 + m3*m2) == ntemps); \
    tassert(b + (m1 + m3*m2) == 2 + m1.to<long>() + m2.to<long>()*m3.to<long>())
#define TT(m1, m2, ntemps) TT3(m1, m2, d, ntemps)

    TT(T0, c, 1);
    TT(T1, c, 1);
    TT(T2, c, 2);
    TT(T3, c, 3);

    TT(T0, T0, 2);
    TT(T0, T1, 2);
    TT(T0, T2, 2);
    TT(T0, T3, 3);

    TT(T1, T0, 2);
    TT(T1, T1, 2);
    TT(T1, T2, 2);
    TT(T1, T3, 3);

    TT(T2, T0, 2);
    TT(T2, T1, 2);
    TT(T2, T2, 3);
    TT(T2, T3, 3);

    TT(T3, T0, 3);
    TT(T3, T1, 3);
    TT(T3, T2, 3);
    TT(T3, T3, 4);

    // NB: TT3 is symmetric in m2 and m3
#define TT6(m1, m2, m3, ntemps) \
    TT3(m1, m2, m3, ntemps); \
    TT3(m2, m1, m3, ntemps); \
    TT3(m3, m1, m2, ntemps);

    TT6(fac(2u), fac(3u), fac(4u), 3);
    TT6(T1, T2, T3, 3);
    TT6(T1, T2, T4, 4);
    TT6(T1, (d+d+d) /* T1' */, T4, 4);
    TT6(T0, fac(2u), T2, 3);
    TT6(T0, T1, (d+d+d), 3);
    TT6(T1, T3, (T2 + T1) /* T3' */, 3);
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
    test_functions();
    test_traits();
    test_temporaries();
    test_ternary();

    // TODO test that certain things *don't* compile?
    // TODO test enable_all_mpz

    std::cout << "PASS" << std::endl;
}
