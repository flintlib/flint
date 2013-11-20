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

#ifndef FMPZ_POLY_QXX_H
#define FMPZ_POLY_QXX_H FMPZ_POLY_QXX_H

#include <cstdlib>
#include <string>

#include "flint.h"
#include "fmpz_poly_q.h"

#include "fmpz_polyxx.h"
#include "fmpqxx.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"

namespace flint {
FLINT_DEFINE_UNOP(fmpz_poly_qxx_num)
FLINT_DEFINE_UNOP(fmpz_poly_qxx_den)

namespace detail {
template<class Polyq>
struct fmpz_poly_q_traits
{
    typedef FLINT_UNOP_BUILD_RETTYPE(
                fmpz_poly_qxx_num, fmpz_polyxx, Polyq) num_ref_t;
    typedef FLINT_UNOP_BUILD_RETTYPE(
                fmpz_poly_qxx_den, fmpz_polyxx, Polyq) den_ref_t;
    typedef num_ref_t num_srcref_t;
    typedef den_ref_t den_srcref_t;
    static num_ref_t num(const Polyq& p) {return fmpz_poly_qxx_num(p);}
    static den_ref_t den(const Polyq& p) {return fmpz_poly_qxx_den(p);}
};
} // detail

template<class Operation, class Data>
class fmpz_poly_qxx_expression
    : public expression<derived_wrapper<fmpz_poly_qxx_expression>,
                            Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpz_poly_qxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpz_poly_q_traits<fmpz_poly_qxx_expression> poly_traits_t;
    typedef typename poly_traits_t::num_ref_t num_ref_t;
    typedef typename poly_traits_t::num_srcref_t num_srcref_t;
    typedef typename poly_traits_t::den_ref_t den_ref_t;
    typedef typename poly_traits_t::den_srcref_t den_srcref_t;

    FLINTXX_DEFINE_BASICS(fmpz_poly_qxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_poly_qxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_poly_qxx_expression, fmpz_poly_q_struct, _polyq)

    // static methods which only make sense with fmpq_polyxx
    static fmpz_poly_qxx_expression randtest(frandxx& state,
            slong len1, mp_bitcnt_t bits1, slong len2, mp_bitcnt_t bits2)
    {
        fmpz_poly_qxx_expression res;
        fmpz_poly_q_randtest(res._polyq(), state._data(),
                len1, bits1, len2, bits2);
        return res;
    }
    static fmpz_poly_qxx_expression randtest_not_zero(frandxx& state,
            slong len1, mp_bitcnt_t bits1, slong len2, mp_bitcnt_t bits2)
    {
        fmpz_poly_qxx_expression res;
        fmpz_poly_q_randtest_not_zero(res._polyq(), state._data(),
                len1, bits1, len2, bits2);
        return res;
    }

    static fmpz_poly_qxx_expression zero(){return fmpz_poly_qxx_expression();}
    static fmpz_poly_qxx_expression one()
    {
        fmpz_poly_qxx_expression res;
        res.set_one();
        return res;
    }

    // Numerator and denominator access
    num_ref_t num() {return poly_traits_t::num(*this);}
    num_srcref_t num() const {return poly_traits_t::num(*this);}
    den_ref_t den() {return poly_traits_t::den(*this);}
    den_srcref_t den() const {return poly_traits_t::den(*this);}
    bool is_canonical() const {return fmpz_poly_q_is_canonical(_polyq());}

    // these only make sense with target immediates
    void canonicalise() {fmpz_poly_q_canonicalise(_polyq());}
    void set_zero() {fmpz_poly_q_zero(_polyq());}
    void set_one() {fmpz_poly_q_one(_polyq());}

    // These cause evaluation
    std::string pretty(const char* x) const
    {
        char* str = fmpz_poly_q_get_str_pretty(this->evaluate()._polyq(), x);
        std::string res(str);
        flint_free(str);
        return res;
    }
    bool is_one() const {return fmpz_poly_q_is_one(this->evaluate()._polyq());}
    bool is_zero() const {return fmpz_poly_q_is_zero(this->evaluate()._polyq());}

    // forwarded lazy member functions
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_UNOP(derivative)
};

namespace detail {
struct fmpz_poly_q_data;
}

typedef fmpz_poly_qxx_expression<operations::immediate, detail::fmpz_poly_q_data>
           fmpz_poly_qxx;
typedef fmpz_poly_qxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_poly_qxx, fmpz_poly_q_struct> >
           fmpz_poly_qxx_ref;
typedef fmpz_poly_qxx_expression<operations::immediate,
            flint_classes::srcref_data<
                fmpz_poly_qxx, fmpz_poly_qxx_ref, fmpz_poly_q_struct> >
           fmpz_poly_qxx_srcref;

namespace detail {
template<>
struct fmpz_poly_q_traits<fmpz_poly_qxx_srcref>
{
    typedef fmpz_polyxx_srcref num_ref_t;
    typedef fmpz_polyxx_srcref num_srcref_t;
    typedef fmpz_polyxx_srcref den_ref_t;
    typedef fmpz_polyxx_srcref den_srcref_t;
    template<class P> static num_srcref_t num(const P& p)
        {return num_srcref_t::make(fmpz_poly_q_numref(p._polyq()));}
    template<class P> static den_srcref_t den(const P& p)
        {return num_srcref_t::make(fmpz_poly_q_denref(p._polyq()));}
};
template<>
struct fmpz_poly_q_traits<fmpz_poly_qxx_ref>
{
    typedef fmpz_polyxx_ref num_ref_t;
    typedef fmpz_polyxx_ref den_ref_t;
    typedef fmpz_polyxx_ref num_srcref_t;
    typedef fmpz_polyxx_ref den_srcref_t;
    template<class P> static num_ref_t num(P p)
        {return num_ref_t::make(fmpz_poly_q_numref(p._polyq()));}
    template<class P> static den_ref_t den(P p)
        {return num_ref_t::make(fmpz_poly_q_denref(p._polyq()));}
};
template<>
struct fmpz_poly_q_traits<fmpz_poly_qxx>
{
    typedef fmpz_polyxx_ref num_ref_t;
    typedef fmpz_polyxx_ref den_ref_t;
    typedef fmpz_polyxx_srcref num_srcref_t;
    typedef fmpz_polyxx_srcref den_srcref_t;
    template<class P> static num_ref_t num(P& p)
        {return num_ref_t::make(fmpz_poly_q_numref(p._polyq()));}
    template<class P> static den_ref_t den(P& p)
        {return num_ref_t::make(fmpz_poly_q_denref(p._polyq()));}
    template<class P> static num_srcref_t num(const P& p)
        {return num_srcref_t::make(fmpz_poly_q_numref(p._polyq()));}
    template<class P> static den_srcref_t den(const P& p)
        {return num_srcref_t::make(fmpz_poly_q_denref(p._polyq()));}
};

struct fmpz_poly_q_data
{
    fmpz_poly_q_t inner;
    typedef fmpz_poly_q_t& data_ref_t;
    typedef const fmpz_poly_q_t& data_srcref_t;

    fmpz_poly_q_data() {fmpz_poly_q_init(inner);}
    ~fmpz_poly_q_data() {fmpz_poly_q_clear(inner);}

    fmpz_poly_q_data(const fmpz_poly_q_data& o)
    {
        fmpz_poly_q_init(inner);
        fmpz_poly_q_set(inner, o.inner);
    }

    fmpz_poly_q_data(fmpz_poly_qxx_srcref r)
    {
        fmpz_poly_q_init(inner);
        fmpz_poly_q_set(inner, r._polyq());
    }

    fmpz_poly_q_data(const char* str)
    {
        fmpz_poly_q_init(inner);
        execution_check(!fmpz_poly_q_set_str(inner, str),
                "construct from string", "fmpz_poly_qxx");
    }
};
} // detail

namespace rules {
#define FMPZ_POLY_QXX_COND_S FLINTXX_COND_S(fmpz_poly_qxx)
#define FMPZ_POLY_QXX_COND_T FLINTXX_COND_T(fmpz_poly_qxx)

FLINTXX_DEFINE_TO_STR(fmpz_poly_qxx, fmpz_poly_q_get_str(from._polyq()))
FLINTXX_DEFINE_SWAP(fmpz_poly_qxx, fmpz_poly_q_swap(e1._polyq(), e2._polyq()))

FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_POLY_QXX_COND_T, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_set(to._polyq(), from._polyq()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_POLY_QXX_COND_T,
        traits::fits_into_slong, fmpz_poly_q_set_si(to._polyq(), from))

FLINTXX_DEFINE_ASSIGN_STR(fmpz_poly_qxx, execution_check(
            !fmpz_poly_q_set_str(to._polyq(), from),
            "assign string", "fmpz_poly_qxx"))

FLINT_DEFINE_PRINT_COND(FMPZ_POLY_QXX_COND_S, fmpz_poly_q_print(from._polyq()))
FLINT_DEFINE_PRINT_PRETTY_COND_2(FMPZ_POLY_QXX_COND_S, const char*,
        fmpz_poly_q_print_pretty(from._polyq(), extra))

FLINTXX_DEFINE_EQUALS(fmpz_poly_qxx, fmpz_poly_q_equal(e1._polyq(), e2._polyq()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpz_poly_qxx, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_neg(to._polyq(), from._polyq()))
FLINT_DEFINE_UNARY_EXPR_COND(inv_op, fmpz_poly_qxx, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_inv(to._polyq(), from._polyq()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_add(to._polyq(), e1._polyq(), e2._polyq()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_sub(to._polyq(), e1._polyq(), e2._polyq()))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, traits::fits_into_slong,
        fmpz_poly_q_scalar_mul_si(to._polyq(), e1._polyq(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, traits::fits_into_slong,
        fmpz_poly_q_scalar_div_si(to._polyq(), e1._polyq(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_mul(to._polyq(), e1._polyq(), e2._polyq()))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_div(to._polyq(), e1._polyq(), e2._polyq()))

FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, traits::is_unsigned_integer,
        fmpz_poly_q_pow(to._polyq(), e1._polyq(), e2))

FLINT_DEFINE_UNARY_EXPR_COND(derivative_op, fmpz_poly_qxx, FMPZ_POLY_QXX_COND_S,
        fmpz_poly_q_derivative(to._polyq(), from._polyq()))

FLINT_DEFINE_UNARY_EXPR_COND(fmpz_poly_qxx_num_op, fmpz_polyxx,
        FMPZ_POLY_QXX_COND_S,
        fmpz_poly_set(to._poly(), fmpz_poly_q_numref(from._polyq())))
FLINT_DEFINE_UNARY_EXPR_COND(fmpz_poly_qxx_den_op, fmpz_polyxx,
        FMPZ_POLY_QXX_COND_S,
        fmpz_poly_set(to._poly(), fmpz_poly_q_denref(from._polyq())))

// XXX these should really be in fmpz_poly_q.h ...
#if 0
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPZXX_COND_S,
        to.numref() = e1.numref()*e2;to.denref() = e1.denref();to.canonicalise())
FLINT_DEFINE_CBINARY_EXPR_COND2(divided_by, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPZXX_COND_S,
        to.denref() = e1.denref()*e2;to.numref() = e1.numref();to.canonicalise())
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPQXX_COND_S,
        to.numref() = e1.numref()*e2.num();to.denref() = e1.denref()*e2.den();
        to.canonicalise())
FLINT_DEFINE_CBINARY_EXPR_COND2(divided_by, fmpz_poly_qxx,
        FMPZ_POLY_QXX_COND_S, FMPQXX_COND_S,
        to.numref() = e1.numref()*e2.den();to.denref() = e1.denref()*e2.num();
        to.canonicalise())
#endif
} // rules
// NB: addmul is not actually optimised currently
} // flint

#endif
