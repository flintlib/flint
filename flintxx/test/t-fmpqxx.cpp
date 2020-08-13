/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <sstream>
#include <string>

#include "fmpqxx.h"
#include "fmpz_vecxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpqxx a, b = fmpqxx::frac(6, fmpzxx(10)),
           c = fmpqxx::frac((unsigned short)6, 10u);

    tassert(a.num() == 0 && a.den() == 1);
    tassert(b.num() == 3 && b.den() == 5);
    tassert(c.num() == 3 && c.den() == 5);

    a.num() = -2;
    a.den() = 4;
    tassert(!a.is_canonical());
    a.canonicalise();
    tassert(a.num() == -1 && a.den() == 2);

    tassert(fmpqxx::frac("4", 2) == fmpqxx::integer(-fmpzxx(-2)));

    tassert(fmpqxx::zero().is_zero() && fmpqxx::one().is_one());
}

void
test_assignment()
{
    fmpqxx a;
    fmpqxx b = fmpqxx::frac("100000000000000000000", "100000000000000000001");
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
    fmpqxx a = fmpqxx::frac(3, 5);
    tassert(a.to_string() == "3/5");
    tassert(a.to_string(5) == "3/10");
    tassert(0.59 < a.to<double>() && a.to<double>() < 0.61);

    std::ostringstream oss;
    oss << a;
    tassert(oss.str() == "3/5");

    tassert((a + a).num() == 6);
    tassert((a * a).den() == 25);
}

void test_order()
{
    fmpqxx a(0, 1u), b(1, 1u);
    fmpzxx c(2), d(-2);

    tassert(a < b);
    tassert(a <= a);
    tassert(b > a);
    tassert(b >= b);
    tassert(a == a);
    tassert(a != b);

    tassert(c > a);
    tassert(c >= a + b);
    tassert(c != a);
    tassert(c > a);
    tassert(a <= c);
    tassert(a + b >= d);
    tassert(a - b != d);
    tassert(a - b > d);
    tassert(d <= b);

    tassert(b < 2u);
    tassert(a + b < c);
    tassert(a + b < 2u);
    tassert(a + b < c.to<ulong>());
    tassert(a + b < 2);
    tassert(a + b < c.to<slong>());
    tassert(1u <= a + b);
    tassert(1 <= b);
    tassert(b == 1u);
    tassert(a + b == 1);
    tassert(2 != b);
    tassert(2u != b);
}

void
test_arithmetic()
{
    fmpqxx a(3, 5u), b(2, 7u);
    fmpzxx c(2);

    tassert(a + b == fmpqxx::frac(3*7 + 2*5, 35u));
    tassert(a + c == fmpqxx::frac(13, 5u));
    tassert(a + c.to<ulong>() == a + c);
    tassert(a + c.to<slong>() == a + c);
    tassert(a + 2u            == a + c);
    tassert(a + 2             == a + c);
    tassert(c + a             == a + c);
    tassert(c.to<slong>() + a == a + c);
    tassert(c.to<ulong>() + a == a + c);
    tassert(2u + a            == a + c);
    tassert(2 + a             == a + c);

    tassert(a * b == fmpqxx::frac(6, 35u));
    tassert(a * c             == a + a);
    tassert(a * c.to<ulong>() == a + a);
    tassert(a * c.to<slong>() == a + a);
    tassert(c * a             == a + a);
    tassert(c.to<ulong>() * a == a + a);
    tassert(c.to<slong>() * a == a + a);
    tassert(a * 2u            == a + a);
    tassert(a * 2             == a + a);
    tassert(2u * a            == a + a);
    tassert(2 * a             == a + a);

    tassert(a - b == fmpqxx::frac(3*7 - 5*2, 35u));
    tassert(a - c == fmpqxx::frac(-7, 5u));
    tassert(a - c.to<ulong>() == a - c);
    tassert(a - c.to<slong>() == a - c);
    tassert(a - 2u            == a - c);
    tassert(a - 2             == a - c);

    tassert(a / b == fmpqxx::frac(3*7, 10u));

    tassert((a+a) * (c+c) == fmpqxx::frac(24, 5u));
    tassert(c * a == fmpqxx::frac(6, 5u));
    tassert(a / c == fmpqxx::frac(3, 10u));

    tassert(-a == fmpqxx::frac(-3, 5u));

    tassert(((a << 5) >> 4) == fmpzxx(2)*a);
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
    tassert((a+a) - b*c == fmpqxx::frac(6*7 - 5*6, 5u*7u));
    tassert(b*c + (a+a) == fmpqxx::frac(6*7 + 5*6, 5u*7u));

    fmpqxx_ref ar(a);
    fmpqxx_srcref asr(a);
    tassert(a == ar && ar == asr);
    ar = 3;
    tassert(a == c && asr == c);

    tassert((-a) + a == fmpqxx::frac(0, 0u));

    tassert(a.pow(3) == pow(a, 3));
    tassert(a.height() == height(a));
}

void
test_functions()
{
    fmpqxx a(-3, 5u);

    // test lazy functions
    tassert(abs(a) == -a);
    tassert(height(a) == 5);
    tassert(a % fmpzxx(7) == 5);
    assert_exception((a % fmpzxx(5)).evaluate());

    // test immediate functions
    tassert(height_bits(a) == 3);
    tassert((inv(a)*a).is_one());
    tassert(sgn(a) == -1 && sgn(-a) == 1 && sgn(fmpqxx::frac(0, 0u)) == 0);

    // test member functions
    const fmpqxx zero(0, 0u), one(1, 1u);
    tassert(zero.is_zero() && !zero.is_one());
    tassert(!one.is_zero() && one.is_one());
    tassert(pow(a, -3) == inv(a*a*a));
    tassert(zero.next_minimal().next_minimal().next_minimal() == fmpqxx::frac(2, 1u));
    tassert(zero.next_signed_minimal().next_signed_minimal() == fmpqxx::frac(-1, 1u));
    tassert(zero.next_calkin_wilf().next_calkin_wilf() == fmpqxx::frac(1, 2u));
    tassert(zero.next_signed_calkin_wilf().next_signed_calkin_wilf()
            == fmpqxx::frac(-1, 1u));

    // test static member functions
    frandxx rand;
    tassert(abs(fmpqxx::randbits(rand, 5).den()) <= 31);
    // NB: rand stuff comes from a single macro, no need for further testing
    tassert(a == fmpqxx::reconstruct(
                a % fmpzxx(41), fmpzxx(41), fmpzxx(3), fmpzxx(5)));
    assert_exception(fmpqxx::reconstruct(
                a % fmpzxx(7), fmpzxx(7), fmpzxx(1), fmpzxx(1)).evaluate());
    tassert(a == fmpqxx::reconstruct(a % fmpzxx(71), fmpzxx(71)));
    assert_exception(fmpqxx::reconstruct(a % fmpzxx(7), fmpzxx(7)).evaluate());

    // test partial fractions
    fmpz_vecxx v(5);
    fmpqxx tmp(7, 5u);
    fmpqxx rem;
    tassert(tmp == fmpqxx::from_cfrac(v, get_cfrac(v, rem, tmp)));
    tassert(rem.is_zero());
    tassert(3 <= tmp.cfrac_bound());

    // test swap
    a = 1;
    fmpqxx b(zero);
    swap(a, b);
    tassert(a.is_zero() && b.is_one());

    tassert_fprint(fmpqxx::frac(7, 29), "7/29");
}

void
test_vector()
{
    fmpq_vecxx v1(10), v2(10), v3(1);
    tassert(v1 == v2);
    tassert(v1 != v3);
    v1[0] = fmpqxx::frac(1, 1u);
    tassert(v1 != v2);
    v2[0] = v1[0];
    tassert(v1 == v2);
}

void
test_unified_access()
{
    fmpqxx a = fmpqxx::frac(1, 2);
    const fmpqxx& b = a;
    tassert(b.num() == 1);
    const fmpqxx_ref c(a);
    c.num() = 3;
    tassert(c.num() == 3);
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
    test_functions();
    test_extras();
    test_vector();
    test_unified_access();

    std::cout << "PASS" << std::endl;
    return 0;
}

