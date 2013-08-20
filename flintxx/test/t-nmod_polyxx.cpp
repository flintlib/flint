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
#include "fmpz_polyxx.h"
#include "nmod_polyxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    nmod_polyxx p(10);
    tassert(p.length() == 0);
    tassert(p.modulus() == 10);

    nmodxx_ctx_srcref ctx = p.estimate_ctx();
    tassert(p == nmod_polyxx::from_ground(0 % ctx));
    p.set_coeff(0, 1);
    tassert(p == nmod_polyxx::from_ground((1 % ctx) + (0 % ctx)));
}

void
test_manipulation()
{
    mp_limb_t M = 31;
    nmod_polyxx p(M), q(M);
    nmodxx_ctx_srcref ctx = p.estimate_ctx();
    p.set_coeff(5, 17 + M);
    tassert(p.degree() == 5);
    q.set_coeff(5, nmodxx::red(17, ctx));
    tassert((q + nmod_polyxx(M)).get_coeff(5) ==
            nmodxx::red(17, ctx));
    p.set_coeff(0, nmodxx::red(1, ctx));
    tassert(p != q);
    p.set_coeff(0, 0);
    tassert(p == q);

    tassert(p.length() == 6);

    p.realloc(0);
    tassert(p.is_zero() && !p.is_one());
    p.set_coeff(0, 1);
    tassert(p.is_one());
}

void
test_assignment()
{
    mp_limb_t M = 31;
    nmod_polyxx p(M), q(M);
    p.set_coeff(0, 1);
    tassert(p != q);
    p = q;
    tassert(p == q);

    p = "4 31  0 0 0 1";
    q.set_coeff(3, 1);
    tassert(p == q);

    // TODO XXX this does not always fail?
    //assert_exception(p = "2 1 2");
    assert_exception(p = "2  x 2");
}

void
test_conversion()
{
    nmod_polyxx p(31);
    p.set_coeff(3, 1);
    tassert(p.to_string() == "4 31  0 0 0 1");
}

void
test_arithmetic()
{
    mp_limb_t M = 31;
    nmod_polyxx g(M), h(M);
    nmodxx_ctx_srcref ctx = g.estimate_ctx();
    g.set_coeff(0, 17); h.set_coeff(0, 15);
    tassert((g + h).get_coeff(0) == nmodxx::red(15 + 17, ctx));

    frandxx state;
    g.set_randtest(state, 10);
    h.set_randtest(state, 10);

    tassert(((-g) + g).is_zero());
    tassert(g - h == g + (-h));

    tassert(g*nmodxx::red(3, ctx) == g + g + g);
    tassert(g.make_monic() == g*inv(g.get_coeff(g.degree())));

    nmod_polyxx f(M);f.set_coeff(0, 15);
    tassert(f*g == nmodxx::red(15, ctx)*g);
    tassert(h.mul_classical(g) == h.mul_KS(g) && h.mul_KS(g) == h*g);

    f = h*g;f.truncate(7);
    tassert(f == mullow(h, g, 7));
    tassert(f == h.mullow_KS(g, 7));
    tassert(f == h.mullow_classical(g, 7));

    f = (h*g).shift_right(7);
    tassert(f == h.mulhigh(g, 7).shift_right(7));
    tassert(f == h.mulhigh_classical(g, 7).shift_right(7));

    f = h / g;
    tassert(f*g + (h % g) == h);
    tassert(((h*g) % h).is_zero());

    f.set_randtest(state, 10);
    tassert(h.mulmod(g, f) == ((h*g) % f));


    f = "3 31  1 0 1";
    nmodxx x = nmodxx::red(7, ctx);
    tassert(evaluate(f, x) == x*x + nmodxx::red(1, ctx));
    f.realloc(0);f.set_coeff(31, 1);
    tassert(evaluate(f, x) == x);
    tassert(f(x) == x);

    nmod_polyxx seven(M);
    seven.set_coeff(0, x);
    tassert(compose(f, seven).get_coeff(0) == f(x));
    tassert(f(seven).length() == 1);

    nmod_vecxx vec1(2, ctx), vec2(2, ctx);
    vec1[0] = nmodxx::red(7, ctx); vec1[1] = nmodxx::red(15, ctx);
    vec2[0] = f(vec1[0]); vec2[1] = f(vec1[1]);
    tassert(f(vec1) == vec2);
}

void
test_functions()
{
    mp_limb_t M = 31;
    nmod_polyxx g(M);
    nmodxx_ctx_srcref ctx = g.estimate_ctx();

    g.set_coeff(5, 15);
    tassert(g.max_bits() == 4);

    g.truncate(3);
    tassert(g.is_zero());

    g.set_coeff(15, 1);
    tassert(g.shift_right(15).is_one());
    tassert(g.shift_right(15).shift_left(15) == g);

    frandxx rand;
    g.set_randtest(rand, 15);
    tassert(g.length() <= 15);
    g.set_randtest_irreducible(rand, 15);
    tassert(g.length() <= 15);
    tassert(g.is_squarefree());
    tassert(g.is_irreducible());

    tassert(g == nmod_polyxx::bit_unpack(g.bit_pack(5u), 5u, ctx));

    // multiplication, division, modulo tested in arithmetic

    tassert(g.pow(3u) == g*g*g);
    tassert(g.pow(5u) == g.pow_binexp(5u));

    nmod_polyxx res(g.pow(15u));res.truncate(12);
    tassert(res == g.pow_trunc(15u, 12));
    tassert(res == g.pow_trunc_binexp(15u, 12));

    nmod_polyxx f(M);f.set_randtest(rand, 10);
    res = g.pow(10u) % f;
    tassert(res == g.powmod_binexp(10u, f));

    res = "5 31  1 1 1 1 1";
    tassert(res.derivative().to_string() == "4 31  1 2 3 4");
    tassert(g.integral().derivative() == g);

    tassert(f.divrem(g) == ltuple(f / g, f % g));
    tassert(f.divrem_basecase(g) == f.divrem(g));
    tassert(f.divrem_divconquer(g) == f.divrem(g));

    tassert(f.div_basecase(g) == f / g);
    tassert(f.div_divconquer(g) == f / g);

    tassert(f.rem_basecase(g) == f % g);

    f.set_coeff(0, 17); // non-zero mod 31, so a unit
    res = f*f.inv_series(15);res.truncate(15);
    tassert(res.is_one());
    tassert(f.inv_series(15) == f.inv_series_basecase(15));
    tassert(f.inv_series(15) == f.inv_series_newton(15));

    res = g * f.inv_series(15);res.truncate(15);
    tassert(g.div_series(f, 15) == res);

    f.set_coeff(f.degree(), 12); // unit
    tassert(g.div_newton(f) == g / f);
    tassert(g.divrem_newton(f) == g.divrem(f));
    tassert(g.divrem(f) == g.divrem_newton21_preinv(f,
                f.reverse(f.length()).inv_series(f.length())));
    tassert(g /f == g.div_newton21_preinv(f,
                f.reverse(f.length()).inv_series(f.length())));

    res = "2 31  5 1";
    tassert(f.div_root(-nmodxx::red(5, ctx)) == f / res);

    nmod_vecxx v(10, ctx);
    _nmod_vec_randtest(v._array(), rand._data(), v.size(), ctx._nmod());
    tassert(f.evaluate_fast(v) == f(v));
    tassert(f.evaluate_iter(v) == f(v));

    nmod_vecxx xs(10, ctx);
    for(unsigned i = 0;i < xs.size();++i)
        xs[i] = nmodxx::red(i, ctx);
    res = nmod_polyxx::interpolate(xs, v);
    tassert(res.degree() < xs.size());
    for(unsigned i = 0;i < xs.size();++i)
        tassert(res(xs[i]) == v[i]);
    tassert(nmod_polyxx::interpolate_fast(xs, v) == res);
    tassert(nmod_polyxx::interpolate_newton(xs, v) == res);
    tassert(nmod_polyxx::interpolate_barycentric(xs, v) == res);

    tassert(f(g) == f.compose_divconquer(g));
    tassert(f(g) == f.compose_horner(g));

    res = "2 31  7 1";
    tassert(f(res) == f.taylor_shift(nmodxx::red(7, ctx)));
    tassert(f(res) == f.taylor_shift_horner(nmodxx::red(7, ctx)));
    tassert(f(res) == f.taylor_shift_convolution(nmodxx::red(7, ctx)));

    nmod_polyxx h(M);
    h.set_randtest(rand, 15);
    tassert(f.compose_mod(g, h) == f(g) % h);
    tassert(f.compose_mod(g, h) == f.compose_mod_horner(g, h));
    tassert(f.compose_mod(g, h) == f.compose_mod_brent_kung(g, h));

    h.set_randtest_irreducible(rand, 12);
    tassert(h.gcd(f).is_one());
    tassert(f.gcd_euclidean(f) == f.make_monic());
    tassert(f.gcd_hgcd(g) == f.gcd(g));

    nmod_polyxx R(M), S(M);
    ltupleref(res, R, S) = f.xgcd(g);
    tassert(res == R*f + S*g && res == gcd(f, g));
    tassert(f.xgcd(g) == f.xgcd_hgcd(g));
    tassert(f.xgcd(g) == f.xgcd_euclidean(g));

    fmpz_polyxx lift1 = fmpz_polyxx::randtest(rand, 10, 6);
    fmpz_polyxx lift2 = fmpz_polyxx::randtest(rand, 10, 6);
    lift1.lead() = 1;
    lift2.lead() = 1;
    f = lift1.reduce(ctx);
    for(unsigned i = 0;i < f.length();++i)
        tassert(f.get_coeff(i) == nmodxx::red(lift1.get_coeff(i), ctx));

    g = lift2.reduce(ctx);
    tassert(f.resultant(g) == nmodxx::red(lift1.resultant(lift2), ctx));
    tassert(f.resultant(g) == f.resultant_euclidean(g));

    g.set_coeff(0, 0);
    res = f(g); res.truncate(15);
    tassert(f.compose_series(g, 15) == res);
    tassert(f.compose_series_horner(g, 15) == res);
    tassert(f.compose_series_brent_kung(g, 15) == res);
    tassert(f.compose_series_divconquer(g, 15) == res);

    res = "2 31  0 1";
    g.set_coeff(1, 17); // unit
    tassert(g.compose_series(g.revert_series(15), 15) == res);
    tassert(g.revert_series_newton(15) == g.revert_series(15));
    tassert(g.revert_series_lagrange(15) == g.revert_series(15));
    tassert(g.revert_series_lagrange_fast(15) == g.revert_series(15));

    f.set_coeff(0, 1);
    tassert(f.sqrt_series(15).pow_trunc(2u, 15) == f);
    tassert(f.invsqrt_series(15).pow_trunc(2u, 15) == f.inv_series(15));

    tassert((f*f).sqrt() == f);
    res = "1 31  1";
    assert_exception((f*f + res).sqrt().evaluate());

    f = nmod_polyxx::product_roots(xs);
    tassert(f.degree() == xs.size());
    for(unsigned i = 0;i < xs.size();++i)
        tassert(f(nmodxx::red(i, ctx)).to<mp_limb_t>() == 0);

    res = "2 31  0 1";
    tassert(f.inflate(5u) == f(res.pow(5u)));
    tassert(f.inflate(5u).deflate(5u) == f);
    tassert(f.inflate(5u).deflation() >= 5);
    tassert(f.deflate(f.deflation()).deflation() == 1);
}

void
test_transcendental_functions()
{
    frandxx state;
    mp_limb_t M = 1031; // prime
    nmod_polyxx f(M);
    nmodxx_ctx_srcref ctx = f.estimate_ctx();
    fmpq_polyxx lift = fmpq_polyxx::randtest(state, 10, 9);
    lift.set_coeff(0, 0);
    f = lift.reduce(ctx);
    for(unsigned i = 0;i < f.length();++i)
        tassert(f.get_coeff(i) == nmodxx::red(lift.get_coeff(i), ctx));

    tassert(f.exp_series(15) == lift.exp_series(15).reduce(ctx));
    tassert(f.atan_series(15) == lift.atan_series(15).reduce(ctx));
    tassert(f.atanh_series(15) == lift.atanh_series(15).reduce(ctx));
    tassert(f.asin_series(15) == lift.asin_series(15).reduce(ctx));
    tassert(f.asinh_series(15) == lift.asinh_series(15).reduce(ctx));
    tassert(f.sin_series(15) == lift.sin_series(15).reduce(ctx));
    tassert(f.cos_series(15) == lift.cos_series(15).reduce(ctx));
    tassert(f.tan_series(15) == lift.tan_series(15).reduce(ctx));
    tassert(f.sinh_series(15) == lift.sinh_series(15).reduce(ctx));
    tassert(f.cosh_series(15) == lift.cosh_series(15).reduce(ctx));
    tassert(f.tanh_series(15) == lift.tanh_series(15).reduce(ctx));

    tassert(f.exp_series_basecase(15) == f.exp_series(15));

    f.set_coeff(0, 1); lift.set_coeff(0, 1);
    tassert(f.log_series(15) == lift.log_series(15).reduce(ctx));

    f.realloc(0);
    nmodxx a = nmodxx::red(7, ctx);
    f.set_coeff(5, a);
    tassert(f.exp_series(15) == exp_series_monomial(a, 5u, 15));
    f.set_coeff(0, 1);
    tassert(f.log_series(15) == log_series_monomial(a, 5u, 15));
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
    std::cout << "nmod_polyxx....";

    test_init();
    test_manipulation();
    test_assignment();
    test_conversion();
    test_arithmetic();
    test_functions();
    test_transcendental_functions();
    test_extras();

    std::cout << "PASS" << std::endl;
    return 0;
}

