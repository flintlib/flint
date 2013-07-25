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

#include "fmpq_polyxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpq_polyxx p(10);
    tassert(p.length() == 0);
}

void
test_manipulation()
{
    fmpq_polyxx p, q;
    p.set_coeff(5, 17);
    p.denref() = 2;
    tassert(p.degree() == 5);
    q.set_coeff(5, fmpqxx(17, 2u));
    tassert((q + fmpq_polyxx()).get_coeff(5) == fmpqxx(17, 2u));
    p.set_coeff(0, fmpzxx(1));
    q.get_coeff_numref(0) = 2;
    tassert(p == q);
    tassert(q.is_canonical());

    tassert(p.length() == 6);

    p.realloc(0);
    tassert(p.is_zero() && !p.is_one());
    p.set_coeff(0, 1);
    tassert(p.is_one());
}

void
test_assignment()
{
    fmpq_polyxx p, q;
    p = 1;
    tassert(p.is_one());
    q = 0ul;
    tassert(q.is_zero());
    tassert(p != q);
    p = q;
    tassert(p == q);

    p = "4  0 0 0 1";
    q.set_coeff(3, 1);
    tassert(p == q);

    // TODO XXX this does not always fail?
    //assert_exception(p = "4  1 2");
    assert_exception(p = "2 1 2");

    fmpz_polyxx pz;
    pz = 1;
    p = pz;
    tassert(p.is_one());

    p = fmpzxx(0);
    tassert(p.is_zero());
    p = fmpqxx(1, 1ul);
    tassert(p.is_one());
}

void
test_conversion()
{
    fmpq_polyxx p;
    p.set_coeff(3, fmpqxx(1, 2u));
    tassert(p.to_string() == "4  0 0 0 1/2");
    tassert(p.pretty("x") == "1/2*x^3");
}

void
test_arithmetic()
{
    fmpq_polyxx p, q;
    p = 1;
    q = 2;
    tassert(p < q);
    p.set_coeff(1, 1);
    tassert(q < p);

    p = -2;
    tassert(p == -q);
    q = fmpqxx(-1, 2u);
    tassert(inv(p) == q);

    p = 1;
    q = "4  0 0 0 1";
    tassert((p + q).to_string() == "4  1 0 0 1");
    tassert((p - q).to_string() == "4  1 0 0 -1");
    tassert((-p).to_string() == "1  -1");

    fmpzxx two(2);
    fmpqxx twoq(2, 1u);
    tassert(two * q == 2 * q && 2u * q == q * 2 && twoq * q == 2 * q
            && (two * q).to_string() == "4  0 0 0 2");
    q *= 2;
    tassert(q / two == q / 2u && q / 2 == q / two && q / twoq == q / 2u &&
            (q / two).to_string() == "4  0 0 0 1");
    // q == "4  0 0 0 2"

    q.set_coeff(1, 17); // q == "4 0 17 0 2"

    p = "3  1 0 1";
    tassert((p*q).to_string() == "6  0 17 0 19 0 2");

    tassert((p*q) / p == q);
    tassert(p + q % q == p);

    fmpqxx five(5, 1u), one(1, 1u);
    tassert(p(one + one) == five);
    tassert(p(fmpzxx(2)) == five);
    q = "3  0 0 1";
    tassert(p(q).to_string() == "5  1 0 0 0 1");
    tassert(p(q) == compose(p, q));
    tassert(p(one) == evaluate(p, one));
    tassert(p(fmpzxx(1)) == evaluate(p, fmpzxx(1)));
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
    // TODO addmul
}

void
test_functions()
{
    fmpq_polyxx f, g, res, x;
    x = "2  0 1";

    f.set_coeff(5, 1);
    f.set_coeff(0, 1);
    tassert(fmpq_polyxx::get_slice(f, 0, 5).is_one());
    f.truncate(3);
    tassert(f.is_one());

    f = "4  1 2 3 4";
    g = "5  5 4 3 2 1";
    res = f*g; res.truncate(4);
    tassert(res == mullow(f, g, 4));

    tassert(pow(f, 5u) == f*f*f*f*f);
    tassert(poly_shift_left(f, 5) == f*pow(x, 5u));
    tassert(f == poly_shift_right(poly_shift_left(f, 5), 5));

    fmpq_polyxx Q, R;
    divrem(Q, R, g, f);
    tassert(Q*f + R == g);

    f = "2  1 -1";
    tassert(inv_series(f, 10).to_string() == "10  1 1 1 1 1 1 1 1 1 1");
    tassert(inv_series_newton(f, 20) == inv_series(f, 20));
    tassert(div_series(g, f, 10) * f % pow(x, 10u) == g);

    f = "4  2 0 0 1";
    g = "5  1 2 3 4 5";
    tassert(resultant(f, g) == fmpqxx(1797, 1u));
    tassert(gcd(f, g).is_one());
    tassert(lcm(f, g) == f*g / 5);

    fmpq_polyxx r, s;
    tassert(xgcd(r, s, f, g).is_one());
    tassert((r*f + s*g).is_one());

    res = "6  0 1 1 1 1 1";
    tassert(derivative(res) == g);
    tassert(integral(g) == res);

    tassert(pow(sqrt_series(g, 10), 2u) % pow(x, 7u) == g);
    tassert((pow(invsqrt_series(g, 10), 2u)*g % pow(x, 7u)).is_one());

    // test static functions
    frandxx state;
    tassert(fmpq_polyxx::randtest(state, 4, 10).length() <= 4);
    tassert(fmpq_polyxx::randtest_unsigned(state, 4, 10).get_coeff_numref(0) >= 0);
    tassert(fmpq_polyxx::randtest_not_zero(state, 4, 10).is_zero() == false);
}

void
test_transcendental_series()
{
    fmpq_polyxx x, xp1;
    x = "2  0 1";
    xp1 = "2  1 1";
    tassert(log_series(xp1, 5).to_string() == "5  0 1 -1/2 1/3 -1/4");
    tassert(atan_series(x, 5).to_string() == "4  0 1 0 -1/3");
    tassert(atanh_series(x, 5).to_string() == "4  0 1 0 1/3");
    tassert(asin_series(x, 5).to_string() == "4  0 1 0 1/6");
    tassert(asinh_series(x, 5).to_string() == "4  0 1 0 -1/6");

    tassert(log_series(exp_series(x, 10), 10) % pow(x, 10u) == x);
    tassert(sin_series(asin_series(x, 10), 10) % pow(x, 10u) == x);
    tassert(tan_series(atan_series(x, 10), 10) % pow(x, 10u) == x);
    tassert(sinh_series(asinh_series(x, 10), 10) % pow(x, 10u) == x);
    tassert(((pow(cos_series(x, 10), 2u)
                + pow(sin_series(x, 10), 2u)) % pow(x, 10u)).is_one());
    tassert(((pow(cosh_series(x, 10), 2u)
                - pow(sinh_series(x, 10), 2u)) % pow(x, 10u)).is_one());
}

int
main()
{
    std::cout << "fmpq_polyxx....";

    test_init();
    test_manipulation();
    test_assignment();
    test_conversion();
    test_arithmetic();
    test_functions();
    test_transcendental_series();
    test_extras();

    std::cout << "PASS" << std::endl;
    return 0;
}

