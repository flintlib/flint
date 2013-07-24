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
    fmpz_polyxx p, q;
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
}

void
test_conversion()
{
    fmpz_polyxx p;
    p.set_coeff(3, 1);
    tassert(p.to_string() == "4  0 0 0 1");
    tassert(p.pretty("x") == "x^3");
}

void
test_arithmetic()
{
    fmpz_polyxx p, q;
    p = 1;
    q = "4  0 0 0 1";
    tassert((p + q).to_string() == "4  1 0 0 1");
    tassert((p - q).to_string() == "4  1 0 0 -1");
    tassert((-p).to_string() == "1  -1");

    fmpzxx two(2);
    tassert(two * q == 2 * q && 2u * q == q * 2
            && (two * q).to_string() == "4  0 0 0 2");
    q *= 2;
    tassert(q / two == q / 2u && q / 2 == q / two &&
            (q / two).to_string() == "4  0 0 0 1");
    // q == "4  0 0 0 2"

    q.set_coeff(1, 17); // q == "4 0 17 0 2"
    tassert((q % fmpzxx(5)).to_string() == "4  0 2 0 2");

    p = "3  1 0 1";
    tassert((p*q).to_string() == "6  0 17 0 19 0 2");

    tassert((p*q) / p == q);
    tassert(p + q % q == p);

    tassert(p(fmpzxx(1) + fmpzxx(1)) == 5);
    q = "3  0 0 1";
    tassert(p(q).to_string() == "5  1 0 0 0 1");
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

ulong pow(ulong base, ulong exp)
{
    ulong res = 1;
    while(exp-- > 0)
        res *= base;
    return base;
}

void
test_functions()
{
    // test swap
    fmpz_polyxx p, q;
    p = 1; q = 0;
    swap(p, q);
    tassert(p.is_zero() && q.is_one());
    // p = 0, q = 1

    fmpz_polyxx f, g, xp1;
    f = "4  2 0 0 1"; // f = x^3 + 2
    g = "5  1 2 3 4 5";
    xp1 = "2  1 1";
    fmpzxx two(2);

    // test lazy functions
    tassert(reverse(q, 4u).pretty("x") == "x^3");
    tassert(mul_2exp(f, 3u) == f * 8);
    tassert(f == fdiv_2exp(mul_2exp(f, 3u), 3u));
    tassert(tdiv(-f, two) == tdiv(-f, 2) && tdiv(-f, 2) == tdiv(-f, 2u)
            && tdiv(-f, two).to_string() == "4  -1 0 0 0");
    tassert(f == divexact(2*f, two) && f == divexact(2*f, 2u)
            && f == divexact(2*f, 2));
    tassert(-f == tdiv_2exp(-8*f - q, 3u));
    tassert(smod(5*f, fmpzxx(3)).to_string() == "4  1 0 0 -1");
    tassert(f == poly_bit_unpack(poly_bit_pack(f, 10u), 10u));
    tassert(f == poly_bit_unpack_unsigned(poly_bit_pack(f, 10u), 10u));

    tassert(mul_classical(f, g) == f*g);
    tassert(mul_karatsuba(f, g) == f*g);
    tassert(mul_SS(f, g) == f*g);
    tassert(mul_KS(f, g) == f*g);

    fmpz_polyxx res; res = f*g; res.truncate(3);
    tassert(mullow_classical(f, g, 3) == res);
    tassert(mullow_karatsuba_n(f, g, 3) == res);
    tassert(mullow_KS(f, g, 3) == res);
    tassert(mullow_SS(f, g, 3) == res);
    tassert(mullow(f, g, 3) == res);

    res = f*g; res.zero_coeffs(0, 5);
    tassert(mulhigh_classical(f, g, 5) == res);
    tassert(mulhigh_karatsuba_n(f, g, 6) == res);
    res = mulhigh_n(f, g, 5);
    res.zero_coeffs(0, 5);
    tassert(res == mulhigh_classical(f, g, 5));

    tassert(mulmid_classical(g, f).to_string() == "2  9 12");

    tassert(sqr(f) == f*f);
    tassert(sqr_KS(f) == f*f);
    tassert(sqr_karatsuba(f) == f*f);
    tassert(sqr_classical(f) == f*f);

    res = sqr(f);res.truncate(4);
    tassert(sqrlow(f, 4) == res);
    tassert(sqrlow_KS(f, 4) == res);
    tassert(sqrlow_karatsuba_n(f, 4) == res);
    tassert(sqrlow_classical(f, 4) == res);

    tassert(pow(f, 3u) == f*f*f);
    tassert(pow_multinomial(f, 3u) == f*f*f);
    tassert(pow_binexp(f, 3u) == f*f*f);
    tassert(pow_addchains(f, 3u) == f*f*f);
    res = pow(f, 10u);
    res.truncate(10);
    tassert(pow_trunc(f, 10u, 10) == res);
    fmpz_polyxx binomial;
    binomial = "2  1 2";
    tassert(pow_binomial(binomial, 3u) == binomial*binomial*binomial);

    fmpz_polyxx x; x = "2  0 1";
    tassert(poly_shift_left(f, 5) == f*pow(x, 5u));
    tassert(poly_shift_right(poly_shift_left(f, 5), 5) == f);

    tassert(height(g) == 5);
    tassert(poly_2norm(g) == 7);

    tassert(resultant(f, g) == 1797);
    tassert(gcd(f, g).is_one());
    tassert(gcd_subresultant(f, g).is_one());
    tassert(gcd_heuristic(f, g).is_one());
    tassert(gcd_modular(f, g).is_one());

    fmpz_polyxx r, s;
    res = 1797;
    tassert(xgcd(r, s, f, g) == 1797);
    tassert(r*f + s*g == res);
    tassert(xgcd_modular(r, s, f, g) == 1797);
    tassert(r*f + s*g == res);

    tassert(lcm(f, g) == f*g);

    tassert(content(2*g) == 2);
    tassert(primitive_part(2*g) == g);

    tassert(div_basecase(f*g, g) == f);
    tassert(div_divconquer(f*g, g) == f);
    tassert(rem_basecase(f + g, g) == f);
    res = 1;
    tassert(div_root(f*(x - res), fmpzxx(1)) == f);

    tassert(inv_series(xp1, 5).to_string() == "5  1 -1 1 -1 1");
    tassert(inv_series(xp1, 10) == inv_series_newton(xp1, 10));

    tassert(derivative(xp1).is_one());

    tassert((compose(xp1, f) - f).is_one());
    tassert((compose_divconquer(xp1, f) - f).is_one());
    tassert((compose_horner(xp1, f) - f).is_one());
    tassert(evaluate(xp1, fmpzxx(1)) == 2);
    tassert(evaluate_horner(xp1, fmpzxx(1)) == 2);
    tassert(evaluate_divconquer(xp1, fmpzxx(1)) == 2);
    tassert(evaluate_mod(xp1, 1, 10) == 2);

    // test immediate functions
    p.set_coeff(3, 1);
    p.truncate(2);
    tassert(p.is_zero());
    tassert(f.max_limbs() == 1);
    tassert(g.max_bits() == 3);

    r = 0; s = 0;
    divrem(r, s, g, f);
    tassert(r*f + s == g);
    tassert(r.to_string() == "2  4 5");
    r = 0; s = 0;
    divrem_basecase(r, s, g, f);
    tassert(r*f + s == g);
    tassert(r.to_string() == "2  4 5");
    r = 0; s = 0;
    divrem_divconquer(r, s, g, f);
    tassert(r*f + s == g);
    tassert(r.to_string() == "2  4 5");

    tassert(poly_divides(res, f*g, g));
    tassert(res == f);

    tassert(div_series(f, xp1, 10) * xp1 % pow(x, 10u) == f);

    f *= 2;
    ulong d = 0;
    pseudo_divrem(r, s, d, g, f);
    tassert(r*f + s == g*pow(2, d));
    r = 0; r = 0; d = 0;
    pseudo_divrem_basecase(r, s, d, g, f);
    tassert(r*f + s == g*pow(2, d));
    r = 0; r = 0; d = 0;
    pseudo_divrem_divconquer(r, s, d, g, f);
    tassert(r*f + s == g*pow(2, d));

    tassert(pseudo_div(d, g, f) == r);
    tassert(pseudo_rem(d, g, f) == s);

    r = 0; r = 0;
    pseudo_divrem_cohen(r, s, g, f);
    tassert(r*f + s == g*4);
    tassert(pseudo_rem_cohen(g, f).to_string() == "3  -28 -32 12");

    // test static functions
    frandxx state;
    tassert(fmpz_polyxx::randtest(state, 4, 10).length() == 4);
    tassert(fmpz_polyxx::randtest_unsigned(state, 4, 10).get_coeff(0) >= 0);
    tassert(fmpz_polyxx::randtest_not_zero(state, 4, 10).is_zero() == false);
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

