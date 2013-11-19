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

#ifndef FMPQ_POLYXX_H
#define FMPQ_POLYXX_H

#include <cstdlib>
#include <string>

#include "flint.h"
#include "fmpq_poly.h"

#include "fmpqxx.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/traits.h"

// TODO exhibit this as a specialisation of a generic poly<fmpzxx>
// TODO lazy version of numref?

namespace flint {
// function "declarations"
FLINT_DEFINE_BINOP(fmpq_polyxx_get_coeff)
FLINT_DEFINE_BINOP(fmpq_polyxx_interpolate)
FLINT_DEFINE_UNOP(fmpq_polyxx_den)

// TODO move to stdmath?

// TODO move to stdmath?
FLINT_DEFINE_BINOP(rescale)

namespace detail {
template<class Poly>
struct fmpq_poly_traits
{
    typedef FLINT_UNOP_BUILD_RETTYPE(fmpq_polyxx_den, fmpzxx, Poly) coeff_ref_t;
    typedef coeff_ref_t coeff_srcref_t;

    static coeff_srcref_t den(const Poly& p) {return fmpq_polyxx_den(p);}
};
} // detail

template<class Operation, class Data>
class fmpq_polyxx_expression
    : public expression<derived_wrapper<fmpq_polyxx_expression>,
                            Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpq_polyxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpq_poly_traits<fmpq_polyxx_expression> poly_traits_t;
    typedef typename poly_traits_t::coeff_ref_t coeff_ref_t;
    typedef typename poly_traits_t::coeff_srcref_t coeff_srcref_t;

    FLINTXX_DEFINE_BASICS(fmpq_polyxx_expression)
    FLINTXX_DEFINE_CTORS(fmpq_polyxx_expression)
    FLINTXX_DEFINE_C_REF(fmpq_polyxx_expression, fmpq_poly_struct, _poly)

    // static methods which only make sense with fmpq_polyxx
    static fmpq_polyxx_expression randtest(frandxx& state, slong len,
            mp_bitcnt_t bits)
    {
        fmpq_polyxx_expression res;
        fmpq_poly_randtest(res._poly(), state._data(), len, bits);
        return res;
    }
    static fmpq_polyxx_expression randtest_unsigned(frandxx& state, slong len,
            mp_bitcnt_t bits)
    {
        fmpq_polyxx_expression res;
        fmpq_poly_randtest_unsigned(res._poly(), state._data(), len, bits);
        return res;
    }
    static fmpq_polyxx_expression randtest_not_zero(frandxx& state, slong len,
            mp_bitcnt_t bits)
    {
        fmpq_polyxx_expression res;
        fmpq_poly_randtest_not_zero(res._poly(), state._data(), len, bits);
        return res;
    }

    // TODO make lazy
    // TODO perhaps as member function (non-static?)
    template<class Poly>
    static fmpq_polyxx_expression get_slice(const Poly& p, slong i, slong j)
    {
        fmpq_polyxx_expression res;
        fmpq_poly_get_slice(res._poly(), p.evaluate()._poly(), i, j);
        return res;
    }

    template<class Fmpq_vec1, class Fmpq_vec2>
    static FLINT_BINOP_ENABLE_RETTYPE(fmpq_polyxx_interpolate,
        Fmpq_vec1, Fmpq_vec2)
    interpolate(const Fmpq_vec1& xs, const Fmpq_vec2& ys)
    {
        return fmpq_polyxx_interpolate(xs, ys);
    }

    static fmpq_polyxx_expression zero() {return fmpq_polyxx_expression();}
    static fmpq_polyxx_expression one()
    {
        fmpq_polyxx_expression res;
        res.set_one();
        return res;
    }

    // These only make sense with immediates
    void realloc(slong alloc) {fmpq_poly_realloc(_poly(), alloc);}
    void fit_length(slong len) {fmpq_poly_fit_length(_poly(), len);}
    void _normalise() {_fmpq_poly_normalise(_poly());}
    void _set_length(slong len) {_fmpq_poly_set_length(_poly(), len);}
    void canonicalise() {fmpq_poly_canonicalise(_poly());}
    bool is_canonical() const {return fmpq_poly_is_canonical(_poly());}
    void set_zero() {fmpq_poly_zero(_poly());}
    void set_one() {fmpq_poly_one(_poly());}

    coeff_ref_t get_coeff_numref(slong n)
    {
        return coeff_ref_t::make(fmpq_poly_numref(_poly()) + n);
    }
    coeff_srcref_t get_coeff_numref(slong n) const
    {
        return coeff_srcref_t::make(fmpq_poly_numref(_poly()) + n);
    }
    coeff_ref_t den() {return poly_traits_t::den(*this);}
    coeff_srcref_t den() const {return poly_traits_t::den(*this);}

    // These only make sense with target immediates
    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type
    set_coeff(slong n, const Fmpz& x)
    {
        fmpq_poly_set_coeff_fmpz(_poly(), n, x.evaluate()._fmpz());
    }
    template<class Fmpq>
    typename mp::enable_if<traits::is_fmpqxx<Fmpq> >::type
    set_coeff(slong n, const Fmpq& x)
    {
        fmpq_poly_set_coeff_fmpq(_poly(), n, x.evaluate()._fmpq());
    }
    template<class T>
    typename mp::enable_if<traits::is_signed_integer<T> >::type
    set_coeff(slong n, T x)
    {
        fmpq_poly_set_coeff_si(_poly(), n, x);
    }
    template<class T>
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type
    set_coeff(slong n, T x)
    {
        fmpq_poly_set_coeff_ui(_poly(), n, x);
    }

    void truncate(slong n) {fmpq_poly_truncate(_poly(), n);}

    // These cause evaluation
    slong length() const {return fmpq_poly_length(this->evaluate()._poly());}
    slong degree() const {return fmpq_poly_degree(this->evaluate()._poly());}
    bool is_one() const {return fmpq_poly_is_one(this->evaluate()._poly());}
    bool is_zero() const {return fmpq_poly_is_zero(this->evaluate()._poly());}
    bool is_monic() const {return fmpq_poly_is_monic(this->evaluate()._poly());}
    bool is_squarefree() const
        {return fmpq_poly_is_squarefree(this->evaluate()._poly());}

    std::string pretty(const char* x) const
    {
        char* str = fmpq_poly_get_str_pretty(this->evaluate()._poly(), x);
        std::string res(str);
        flint_free(str);
        return res;
    }

    // lazy members
    FLINTXX_DEFINE_MEMBER_BINOP_(operator(), compeval)

    template<class Slong> // NB: template to instantiate lazily
    FLINT_BINOP_ENABLE_RETTYPE(
            fmpq_polyxx_get_coeff, fmpq_polyxx_expression, Slong)
    get_coeff(const Slong& n) const
    {
        return fmpq_polyxx_get_coeff(*this, n);
    }

    FLINTXX_DEFINE_MEMBER_3OP(compose_series)
    FLINTXX_DEFINE_MEMBER_3OP(compose_series_brent_kung)
    FLINTXX_DEFINE_MEMBER_3OP(compose_series_horner)
    FLINTXX_DEFINE_MEMBER_3OP(div_series)
    FLINTXX_DEFINE_MEMBER_3OP(mullow)
    FLINTXX_DEFINE_MEMBER_BINOP(asinh_series)
    FLINTXX_DEFINE_MEMBER_BINOP(asin_series)
    FLINTXX_DEFINE_MEMBER_BINOP(atanh_series)
    FLINTXX_DEFINE_MEMBER_BINOP(atan_series)
    FLINTXX_DEFINE_MEMBER_BINOP(cosh_series)
    FLINTXX_DEFINE_MEMBER_BINOP(cos_series)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem)
    FLINTXX_DEFINE_MEMBER_BINOP(exp_series)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd)
    FLINTXX_DEFINE_MEMBER_BINOP(inv_series)
    FLINTXX_DEFINE_MEMBER_BINOP(inv_series_newton)
    FLINTXX_DEFINE_MEMBER_BINOP(invsqrt_series)
    FLINTXX_DEFINE_MEMBER_BINOP(lcm)
    FLINTXX_DEFINE_MEMBER_BINOP(log_series)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_BINOP_(rescale, rescale)
    FLINTXX_DEFINE_MEMBER_BINOP(resultant)
    FLINTXX_DEFINE_MEMBER_BINOP(reverse)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series_lagrange)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series_lagrange_fast)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series_newton)
    FLINTXX_DEFINE_MEMBER_BINOP_(shift_left, shift_left)
    FLINTXX_DEFINE_MEMBER_BINOP_(shift_right, shift_right)
    FLINTXX_DEFINE_MEMBER_BINOP(sinh_series)
    FLINTXX_DEFINE_MEMBER_BINOP(sin_series)
    FLINTXX_DEFINE_MEMBER_BINOP(sqrt_series)
    FLINTXX_DEFINE_MEMBER_BINOP(tanh_series)
    FLINTXX_DEFINE_MEMBER_BINOP(tan_series)
    FLINTXX_DEFINE_MEMBER_BINOP(xgcd)
    FLINTXX_DEFINE_MEMBER_UNOP(derivative)
    FLINTXX_DEFINE_MEMBER_UNOP(integral)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_UNOP(make_monic)
    FLINTXX_DEFINE_MEMBER_UNOP(primitive_part)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpqxx, content)
};

namespace detail {
struct fmpq_poly_data;
}

typedef fmpq_polyxx_expression<operations::immediate, detail::fmpq_poly_data>
           fmpq_polyxx;
typedef fmpq_polyxx_expression<operations::immediate,
            flint_classes::ref_data<fmpq_polyxx, fmpq_poly_struct> >
           fmpq_polyxx_ref;
typedef fmpq_polyxx_expression<operations::immediate,
            flint_classes::srcref_data<
                fmpq_polyxx, fmpq_polyxx_ref, fmpq_poly_struct> >
           fmpq_polyxx_srcref;

namespace detail {
template<>
struct fmpq_poly_traits<fmpq_polyxx_srcref>
{
    typedef fmpzxx_srcref coeff_ref_t;
    typedef fmpzxx_srcref coeff_srcref_t;
    template<class P> static coeff_srcref_t den(const P& p)
        {return coeff_srcref_t::make(fmpq_poly_denref(p._poly()));}
};
template<>
struct fmpq_poly_traits<fmpq_polyxx_ref>
{
    typedef fmpzxx_ref coeff_ref_t;
    typedef fmpzxx_ref coeff_srcref_t;
    template<class P>
    static coeff_ref_t den(P p)
        {return coeff_ref_t::make(fmpq_poly_denref(p._poly()));}
};
template<>
struct fmpq_poly_traits<fmpq_polyxx>
{
    typedef fmpzxx_ref coeff_ref_t;
    typedef fmpzxx_srcref coeff_srcref_t;
    template<class P>
    static coeff_ref_t den(P& p)
        {return coeff_ref_t::make(fmpq_poly_denref(p._poly()));}
    template<class P>
    static coeff_srcref_t den(const P& p)
        {return coeff_srcref_t::make(fmpq_poly_denref(p._poly()));}
};

struct fmpq_poly_data
{
    fmpq_poly_t inner;
    typedef fmpq_poly_t& data_ref_t;
    typedef const fmpq_poly_t& data_srcref_t;

    fmpq_poly_data() {fmpq_poly_init(inner);}
    ~fmpq_poly_data() {fmpq_poly_clear(inner);}

    fmpq_poly_data(const fmpq_poly_data& o)
    {
        fmpq_poly_init(inner);
        fmpq_poly_set(inner, o.inner);
    }

    fmpq_poly_data(fmpq_polyxx_srcref r)
    {
        fmpq_poly_init(inner);
        fmpq_poly_set(inner, r._poly());
    }

    fmpq_poly_data(slong alloc)
    {
        fmpq_poly_init2(inner, alloc);
    }

    fmpq_poly_data(const char* str)
    {
        fmpq_poly_init(inner);
        execution_check(!fmpq_poly_set_str(inner, str),
                "construct from string", "fmpq_polyxx");
    }
};
} // detail

namespace traits {
template<class T> struct is_fmpq_polyxx : mp::or_<
     traits::is_T_expr<T, fmpq_polyxx>,
     flint_classes::is_source<fmpq_polyxx, T> > { };
} // traits
namespace mp {
template<class T1, class T2 = void, class T3 = void, class T4 = void>
struct all_fmpq_polyxx : mp::and_<all_fmpq_polyxx<T1>, all_fmpq_polyxx<T2, T3, T4> > { };
template<class T>
struct all_fmpq_polyxx<T, void, void, void> : traits::is_fmpq_polyxx<T> { };

template<class Out, class T1, class T2 = void, class T3 = void, class T4 = void>
struct enable_all_fmpq_polyxx
    : mp::enable_if<all_fmpq_polyxx<T1, T2, T3, T4>, Out> { };
} // mp
} // flint

// here to deal with circular dependencies...
#include "fmpz_polyxx.h"

namespace flint {
namespace rules {
#define FMPQ_POLYXX_COND_S FLINTXX_COND_S(fmpq_polyxx)
#define FMPQ_POLYXX_COND_T FLINTXX_COND_T(fmpq_polyxx)

FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_POLYXX_COND_T, FMPQ_POLYXX_COND_S,
        fmpq_poly_set(to._poly(), from._poly()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_POLYXX_COND_T,
        traits::is_signed_integer, fmpq_poly_set_si(to._poly(), from))
FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_POLYXX_COND_T,
        traits::is_unsigned_integer, fmpq_poly_set_ui(to._poly(), from))
FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_POLYXX_COND_T, FMPZXX_COND_S,
        fmpq_poly_set_fmpz(to._poly(), from._fmpz()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_POLYXX_COND_T, FMPQXX_COND_S,
        fmpq_poly_set_fmpq(to._poly(), from._fmpq()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_POLYXX_COND_T, FMPZ_POLYXX_COND_S,
        fmpq_poly_set_fmpz_poly(to._poly(), from._poly()))

FLINTXX_DEFINE_ASSIGN_STR(fmpq_polyxx, execution_check(
            !fmpq_poly_set_str(to._poly(), from), "assign string", "fmpq_polyxx"))

FLINTXX_DEFINE_TO_STR(fmpq_polyxx, fmpq_poly_get_str(from._poly()))
FLINTXX_DEFINE_SWAP(fmpq_polyxx, fmpq_poly_swap(e1._poly(), e2._poly()))

FLINT_DEFINE_PRINT_COND(FMPQ_POLYXX_COND_S, fmpq_poly_fprint(to, from._poly()))
FLINT_DEFINE_PRINT_PRETTY_COND_2(FMPQ_POLYXX_COND_S, const char*,
        fmpq_poly_fprint_pretty(to, from._poly(), extra))
FLINT_DEFINE_READ_COND(FMPQ_POLYXX_COND_T, fmpq_poly_fread(from, to._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(reverse_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_reverse(to._poly(), e1._poly(), e2))

FLINTXX_DEFINE_CMP(fmpq_polyxx, fmpq_poly_cmp(e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpq_polyxx, FMPQ_POLYXX_COND_S,
        fmpq_poly_neg(to._poly(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(inv_op, fmpq_polyxx, FMPQ_POLYXX_COND_S,
        fmpq_poly_inv(to._poly(), from._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(fmpq_polyxx_get_coeff_op, fmpqxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_get_coeff_fmpq(to._fmpq(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_add(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_sub(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpq_poly_scalar_mul_fmpz(to._poly(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQXX_COND_S,
        fmpq_poly_scalar_mul_fmpq(to._poly(), e1._poly(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::is_signed_integer,
        fmpq_poly_scalar_mul_si(to._poly(), e1._poly(), e2))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpq_poly_scalar_mul_ui(to._poly(), e1._poly(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND2(divided_by, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQXX_COND_S,
        fmpq_poly_scalar_div_fmpq(to._poly(), e1._poly(), e2._fmpq()))
FLINT_DEFINE_CBINARY_EXPR_COND2(divided_by, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpq_poly_scalar_div_fmpz(to._poly(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpq_poly_scalar_div_ui(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::is_signed_integer,
        fmpq_poly_scalar_div_si(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_mul(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpq_poly_pow(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(shift_left_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_shift_left(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(shift_right_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_shift_right(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_div(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_rem(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(inv_series_newton_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_inv_series_newton(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(inv_series_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_inv_series(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(gcd_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_gcd(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(lcm_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_lcm(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(resultant_op, fmpqxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_resultant(to._fmpq(), e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(derivative_op, fmpq_polyxx, FMPQ_POLYXX_COND_S,
        fmpq_poly_derivative(to._poly(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(integral_op, fmpq_polyxx, FMPQ_POLYXX_COND_S,
        fmpq_poly_integral(to._poly(), from._poly()))

#define FMPQ_POLYXX_DEFINE_SERIES(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_series_op, fmpq_polyxx, \
        FMPQ_POLYXX_COND_S, traits::fits_into_slong, \
        fmpq_poly_##name##_series(to._poly(), e1._poly(), e2))
FMPQ_POLYXX_DEFINE_SERIES(sqrt)
FMPQ_POLYXX_DEFINE_SERIES(invsqrt)
FMPQ_POLYXX_DEFINE_SERIES(exp)
FMPQ_POLYXX_DEFINE_SERIES(log)
FMPQ_POLYXX_DEFINE_SERIES(atan)
FMPQ_POLYXX_DEFINE_SERIES(atanh)
FMPQ_POLYXX_DEFINE_SERIES(asin)
FMPQ_POLYXX_DEFINE_SERIES(asinh)
FMPQ_POLYXX_DEFINE_SERIES(tan)
FMPQ_POLYXX_DEFINE_SERIES(sin)
FMPQ_POLYXX_DEFINE_SERIES(cos)
FMPQ_POLYXX_DEFINE_SERIES(sinh)
FMPQ_POLYXX_DEFINE_SERIES(cosh)
FMPQ_POLYXX_DEFINE_SERIES(tanh)

FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_op, fmpqxx,
        FMPQ_POLYXX_COND_S, FMPQXX_COND_S,
        fmpq_poly_evaluate_fmpq(to._fmpq(), e1._poly(), e2._fmpq()))
FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_op, fmpqxx,
        FMPQ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpq_poly_evaluate_fmpz(to._fmpq(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(compose_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
        fmpq_poly_compose(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(fmpq_polyxx_interpolate_op, fmpq_polyxx,
        FMPZ_VECXX_COND_S, FMPZ_VECXX_COND_S,
        fmpq_poly_interpolate_fmpz_vec(to._poly(), e1._data().array,
            e2._data().array, e2.size()))

FLINT_DEFINE_BINARY_EXPR_COND2(rescale_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, FMPQXX_COND_S,
        fmpq_poly_rescale(to._poly(), e1._poly(), e2._fmpq()))

FLINT_DEFINE_BINARY_EXPR_COND2(revert_series_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_revert_series(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(revert_series_lagrange_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_revert_series_lagrange(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(revert_series_lagrange_fast_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_revert_series_lagrange_fast(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(revert_series_newton_op, fmpq_polyxx,
        FMPQ_POLYXX_COND_S, traits::fits_into_slong,
        fmpq_poly_revert_series_newton(to._poly(), e1._poly(), e2))

FLINT_DEFINE_UNARY_EXPR_COND(content_op, fmpqxx, FMPQ_POLYXX_COND_S,
        fmpq_poly_content(to._fmpq(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(primitive_part_op, fmpq_polyxx, FMPQ_POLYXX_COND_S,
        fmpq_poly_primitive_part(to._poly(), from._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(make_monic_op, fmpq_polyxx, FMPQ_POLYXX_COND_S,
        fmpq_poly_make_monic(to._poly(), from._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(fmpq_polyxx_den_op, fmpzxx, FMPQ_POLYXX_COND_S,
        fmpz_set(to._fmpz(), fmpq_poly_denref(from._poly())))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<
    fmpq_polyxx, fmpq_polyxx, fmpq_polyxx>::type>::type fmpq_polyxx_triple;
typedef make_ltuple<mp::make_tuple<
    fmpq_polyxx, fmpq_polyxx>::type>::type fmpq_polyxx_pair;
} // rdetail

FLINT_DEFINE_BINARY_EXPR_COND2(xgcd_op, rdetail::fmpq_polyxx_triple,
    FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
    fmpq_poly_xgcd(to.template get<0>()._poly(), to.template get<1>()._poly(),
        to.template get<2>()._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(divrem_op, rdetail::fmpq_polyxx_pair,
    FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S,
    fmpq_poly_divrem(to.template get<0>()._poly(), to.template get<1>()._poly(),
        e1._poly(), e2._poly()))

#define FMPQ_POLYXX_DEFINE_SERIES_FUNC(name) \
FLINT_DEFINE_THREEARY_EXPR_COND3(name##_op, fmpq_polyxx, \
    FMPQ_POLYXX_COND_S, FMPQ_POLYXX_COND_S, traits::fits_into_slong, \
    fmpq_poly_##name(to._poly(), e1._poly(), e2._poly(), e3))
FMPQ_POLYXX_DEFINE_SERIES_FUNC(mullow)
FMPQ_POLYXX_DEFINE_SERIES_FUNC(div_series)
FMPQ_POLYXX_DEFINE_SERIES_FUNC(compose_series)
FMPQ_POLYXX_DEFINE_SERIES_FUNC(compose_series_brent_kung)
FMPQ_POLYXX_DEFINE_SERIES_FUNC(compose_series_horner)
} // rules

// NB: fmpq_poly addmul is just done by hand currently, no need to wrap that ..

} // flint

#endif
