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

#include "fmpz_mod_polyxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpz_mod_polyxx p(fmpzxx(2003));
    tassert(p.length() == 0);
    tassert(p.modulus() == 2003);
}

void
test_manipulation()
{
    fmpzxx M(1031);
    fmpz_mod_polyxx p(M), q(M);
    p.set_coeff(5, 17u + M);
    tassert(p.degree() == 5);
    q.set_coeff(5, fmpzxx(16) + fmpzxx(1));
    tassert((q + fmpz_mod_polyxx(M)).get_coeff(5) == 17);
    p.set_coeff(0, 1);
    tassert(p != q);
    p.set_coeff(0, 0);
    tassert(p == q);

    tassert((p + p).lead() == 2*p.lead());

    q.lead() = 0;
    q._normalise();
    tassert(q.is_zero());

    p.realloc(0);
    tassert(p.is_zero());
}

void
test_assignment()
{
    fmpzxx M(31);
    fmpz_mod_polyxx p(M), q(M);
    p.set_coeff(0, 1);
    tassert(p != q);
    p = q;
    tassert(p == q);
}

void
test_conversion()
{
    fmpzxx M(1031);
    fmpz_mod_polyxx p(M);

    p = 4u + 1031;
    tassert(p.length() == 1 && p.get_coeff(0) == 4);

    p = fmpzxx(5) + M;
    tassert(p.length() == 1 && p.get_coeff(0) == 5);

    frandxx rand;
    fmpz_polyxx P = fmpz_polyxx::randtest(rand, 10, 20);
    p = P;
    for(slong i = 0;i < P.length();++i)
        tassert(P.get_coeff(i) % M == p.get_coeff(i));
    fmpz_polyxx Pp = p.to<fmpz_polyxx>();
    for(slong i = 0;i < P.length();++i)
        tassert(P.get_coeff(i) % M == Pp.get_coeff(i));
}

void
test_arithmetic()
{
    fmpzxx M(1031);
    fmpz_mod_polyxx g(M), h(M);
    g.set_coeff(0, 17); h.set_coeff(0, 15u + M);
    tassert((g + h).get_coeff(0) == 15 + 17);

    frandxx state;
    g.set_randtest(state, 10);
    h.set_randtest(state, 10);

    tassert(((-g) + g).is_zero());
    tassert(g - h == g + (-h));

    tassert(g*fmpzxx(3) == g + g + g);
    tassert(g.make_monic() == g*g.lead().invmod(M));

    fmpz_mod_polyxx f(M);f = 15u;
    tassert(f*g == fmpzxx(15)*g);

    f = h*g;f.truncate(7);
    tassert(f == mullow(h, g, 7));

    f = h / g;
    tassert(f*g + (h % g) == h);
    tassert(((h*g) % h).is_zero());

    f.set_randtest(state, 10);
    tassert(h.mulmod(g, f) == ((h*g) % f));

    fmpz_mod_polyxx X(M);X.set_coeff(1, 1);
    fmpz_mod_polyxx one(M);one.set_coeff(0, 1);
    f = X*X + one;
    fmpzxx x(7);
    tassert(evaluate(f, x) == x*x + 1u);
    tassert(f(x) == evaluate(f, x));

    fmpz_mod_polyxx seven(M);
    seven.set_coeff(0, x);
    tassert(compose(f, seven).get_coeff(0) == f(x));
    tassert(f(seven).length() == 1);
}

void
test_functions()
{
    fmpzxx M(1031);
    fmpz_mod_polyxx g(M), res(M);

    g.set_coeff(5, 15);
    //tassert(g.max_bits() == 4);

    g.truncate(3);
    tassert(g.is_zero());

    g.set_coeff(15, 1);
    fmpz_mod_polyxx one(M);one = 1u;
    tassert(g.poly_shift_right(15) == one);
    tassert(g.poly_shift_right(15).poly_shift_left(15) == g);

    frandxx rand;
    g.set_randtest(rand, 15);
    tassert(g.length() <= 15);
    g.set_randtest_irreducible(rand, 15);
    tassert(g.length() <= 15);
    g.set_randtest_not_zero(rand, 15);
    tassert(g.length() <= 15 && !g.is_zero());

    g.set_coeff(15, 1);
    g.zero_coeffs(14, 15);
    tassert(g.get_coeff(14) == 0);

    //tassert(g == fmpz_mod_polyxx::bit_unpack(g.bit_pack(5u), 5u, ctx));

    // multiplication, division, modulo tested in arithmetic

    tassert(g.pow(3u) == g*g*g);

    res = g.pow(15u);res.truncate(12);
    tassert(res == g.pow_trunc(15u, 12));
    tassert(res == g.pow_trunc_binexp(15u, 12));

    fmpz_mod_polyxx f(M);f.set_randtest(rand, 10);
    res = g.pow(10u) % f;
    tassert(res == g.powmod_binexp(10u, f));
    tassert(res == g.powmod_binexp(fmpzxx(10), f));

    fmpz_mod_polyxx tmp(M);
    ltupleref(res, tmp) = f.gcdinv(g);
    tassert(res == gcd(f, g) && tmp*f % g == res);

    g.set_randtest_irreducible(rand, 5);
    tassert(f.invmod(g)*f % g == one);
    assert_exception((f*g).invmod(g).evaluate());

    res = g*f;
    res.remove(f);
    tassert(res == g);

    fmpz_polyxx lift;
    lift = "5  1 1 1 1 1";
    res = lift;
    tassert(res.derivative().to<fmpz_polyxx>().to_string() == "4  1 2 3 4");
    //tassert(g.integral().derivative() == g);

    tassert(f.divrem(g) == ltuple(f / g, f % g));
    tassert(f.divrem_basecase(g) == f.divrem(g));
    tassert(f.divrem_divconquer(g) == f.divrem(g));
    tassert(f.divrem_f(g) == ltuple(1, f / g, f % g));

    tassert(f.div_basecase(g) == f / g);
    tassert(f.rem_basecase(g) == f % g);

    f.set_coeff(0, 17);
    res = f*f.inv_series_newton(15);res.truncate(15);
    tassert(res == one);

    tassert(f(g) == f.compose_divconquer(g));
    tassert(f(g) == f.compose_horner(g));

    fmpz_mod_polyxx h(M);
    h.set_randtest(rand, 15);
    tassert(f.compose_mod(g, h) == f(g) % h);
    tassert(f.compose_mod(g, h) == f.compose_mod_horner(g, h));
    tassert(f.compose_mod(g, h) == f.compose_mod_brent_kung(g, h));

    h.set_randtest_irreducible(rand, 12);
    tassert(h.gcd(f) == one);
    tassert(f.gcd_euclidean(f) == f.make_monic());
    tassert(f.gcd_f(g) == ltuple(1, f.gcd(g)));
    tassert(f.gcd_euclidean_f(g) == ltuple(1, f.gcd(g)));

    fmpz_mod_polyxx R(M), S(M);
    ltupleref(res, R, S) = f.xgcd(g);
    tassert(res == R*f + S*g && res == gcd(f, g));
    tassert(f.xgcd(g) == f.xgcd_euclidean(g));

}

// test stuff which we should get automatically - addmul, references etc
void
test_extras()
{
    // TODO
}

int
main()
{
    std::cout << "fmpz_mod_polyxx....";

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

