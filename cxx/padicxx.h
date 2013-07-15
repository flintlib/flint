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

#ifndef CXX_PADICXX_H
#define CXX_PADICXX_H CXX_PADICXX_H

#include <algorithm> // std::max
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "padic.h"

#include "cxx/expression.h"
#include "cxx/stdmath.h"
#include "cxx/traits.h"
#include "cxx/tuple.h"

#include "fmpzxx.h"
#include "fmpqxx.h"

namespace flint {
class padicxx_exception
    : public std::domain_error
{
public:
    padicxx_exception(const std::string& what)
        : std::domain_error("padic computation failed: " + what) {}
};

namespace detail {
template<class T>
struct padicxx_max_prec
{
    // XXX is this a good idea?
    static slong get(const T&) {return 0;}
};

void padic_check(bool worked, const std::string& where)
{
    if(!worked)
        throw padicxx_exception(where);
}
} // detail

class padicxx_ctx
{
private:
    mutable padic_ctx_t ctx;

public:
    // NB: you must not modify user-visible state of ctx through a constant
    // instance of padicxx_ctx
    padic_ctx_t& _ctx() const {return ctx;}

    // XXX these two are not actually exposed in the C api ...
    fmpzxx_ref get_p() {return fmpzxx_ref::make(ctx[0].p);}
    fmpzxx_srcref get_p() const {return fmpzxx_srcref::make(ctx[0].p);}

    // TODO more constructors? Should we wrap padic_print_mode?
    padicxx_ctx(fmpzxx_srcref p, long min, long max, padic_print_mode mode)
    {
        padic_ctx_init(ctx, p._fmpz(), min, max, mode);
    }

    ~padicxx_ctx() {padic_ctx_clear(ctx);}
};

template<class Operation, class Data>
class padicxx_expression
    : public expression<derived_wrapper<padicxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::padicxx_expression>,
              Operation, Data> base_t;
    typedef typename base_t::evaluated_t evaluated_t;

    template<class T>
    explicit padicxx_expression(const T& t)
        : base_t(t) {}

    template<class T, class U>
    padicxx_expression(const T& t, const U& u)
        : base_t(t, u) {}

    template<class T>
    padicxx_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

    // These only make sense with immediates
    padic_t& _padic() {return this->_data().p;}
    const padic_t& _padic() const {return this->_data().p;}
    const padicxx_ctx& get_ctx() const {return this->_data().ctx;}
    fmpzxx_ref unit() {return fmpzxx_ref::make(padic_unit(_padic()));}
    fmpzxx_srcref unit() const {return fmpzxx_srcref::make(padic_unit(_padic()));}
    slong val() const {return padic_val(_padic());}
    slong _prec() const {return padic_prec(_padic());}
    padic_ctx_t& _ctx() const {return get_ctx()._ctx();}
    // TODO reduce? canonicalise?

    // Compute the maximal precision of all subexpressions
    slong prec() const
    {
        return detail::padicxx_max_prec<padicxx_expression>::get(*this);
    }

    // The above method get_ctx() only works on immediates, i.e. instances of
    // padicxx. This next one here works on any composite expression which
    // contains at least one instance of padicxx. Returns the context of one of
    // those immediate subexpressions.
    const padicxx_ctx& estimate_ctx() const;

    // Create a temporary. The context will be estimated, and the precision
    // will be the maximum of all subexpressions.
    // TODO should we incorporate target precision, if available?
    evaluated_t create_temporary() const
    {
        return evaluated_t(estimate_ctx(), prec());
    }

protected:
    explicit padicxx_expression(const Data& d) : base_t(d) {}

    template<class D, class O, class Da>
    friend class expression;
};

namespace detail {
struct padic_data
{
    const padicxx_ctx& ctx;
    padic_t p;

    padic_data(const padicxx_ctx& c)
        : ctx(c)
    {
        padic_init(p);
    }

    padic_data(const padicxx_ctx& c, long N)
        : ctx(c)
    {
        padic_init2(p, N);
    }

    padic_data(const padic_data& o)
        : ctx(o.ctx)
    {
        padic_init2(p, padic_prec(o.p));
        padic_set(p, o.p, ctx._ctx());
    }

    // TODO more constructors? (e.g. unit, val?)

    ~padic_data() {padic_clear(p);}
};
} // detail

typedef padicxx_expression<operations::immediate, detail::padic_data> padicxx;

template<class Operation, class Data>
inline const padicxx_ctx&
padicxx_expression<Operation, Data>::estimate_ctx() const
{
    return tools::find_subexpr_T<padicxx>(*this).get_ctx();
}

namespace detail {
template<>
struct padicxx_max_prec<padicxx>
{
    static slong get(const padicxx& p) {return p._prec();}
};

template<class Data>
struct padicxx_max_prec_h;
template<class Head, class Tail>
struct padicxx_max_prec_h<tuple<Head, Tail> >
{
    static slong get(const tuple<Head, Tail>& t)
    {
        slong p1 = padicxx_max_prec_h<Tail>::get(t.tail);
        slong p2 =
            padicxx_max_prec<typename traits::basetype<Head>::type>::get(t.head);
        return std::max(p1, p2);
    }
};
template<>
struct padicxx_max_prec_h<empty_tuple>
{
    static slong get(empty_tuple) {return 0;}
};

template<class Op, class Data>
struct padicxx_max_prec<padicxx_expression<Op, Data> >
{
    static slong get(const padicxx_expression<Op, Data>& e)
    {
        return padicxx_max_prec_h<Data>::get(e._data());
    }
};
} // detail

namespace rules {
FLINT_DEFINE_DOIT(assignment, padicxx, padicxx,
        padic_set(to._padic(), from._padic(), to._ctx()))
FLINT_DEFINE_DOIT_COND(assignment, padicxx, traits::is_signed_integer, 
        padic_set_si(to._padic(), from, to._ctx()))
FLINT_DEFINE_DOIT_COND(assignment, padicxx, traits::is_unsigned_integer, 
        padic_set_ui(to._padic(), from, to._ctx()))
FLINT_DEFINE_DOIT_COND(assignment, padicxx, FMPZXX_COND_S,
        padic_set_fmpz(to._padic(), from._fmpz(), to._ctx()))
FLINT_DEFINE_DOIT(assignment, padicxx, fmpqxx,
        padic_set_fmpq(to._padic(), from._fmpq(), to._ctx()))

template<>
struct conversion<fmpzxx, padicxx>
{
    static fmpzxx get(const padicxx& from)
    {
        fmpzxx res;
        padic_get_fmpz(res._fmpz(), from._padic(), from._ctx());
        return res;
    }
};

template<>
struct conversion<fmpqxx, padicxx>
{
    static fmpqxx get(const padicxx& from)
    {
        fmpqxx res;
        padic_get_fmpq(res._fmpq(), from._padic(), from._ctx());
        return res;
    }
};

template<>
struct to_string<padicxx>
{
    static std::string get(const padicxx& p, int base /* ignored! */)
    {
        char* str = padic_get_str(0, p._padic(), p._ctx());
        std::string res(str);
        std::free(str);
        return res;
    }
};

template<>
struct equals<padicxx, padicxx>
{
    static bool get(const padicxx& p1, const padicxx& p2)
    {
        return padic_equal(p1._padic(), p2._padic());
    }
};

FLINT_DEFINE_CBINARY_EXPR(plus, padicxx,
        padic_add(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR(minus, padicxx,
        padic_sub(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_CBINARY_EXPR(times, padicxx,
        padic_mul(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR(divided_by, padicxx,
        padic_div(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND(shift, padicxx,
        traits::fits_into_slong,
        padic_shift(to._padic(), e1._padic(), e2, to._ctx()))

FLINT_DEFINE_UNARY_EXPR(negate, padicxx,
        padic_neg(to._padic(), from._padic(), to._ctx()))

// lazy functions
FLINT_DEFINE_UNARY_EXPR(sqrt_op, padicxx,
        detail::padic_check(
            padic_sqrt(to._padic(), from._padic(), to._ctx()), "sqrt"))
FLINT_DEFINE_BINARY_EXPR_COND(pow_op, padicxx, traits::fits_into_slong,
        padic_pow_si(to._padic(), e1._padic(), e2, to._ctx()))
FLINT_DEFINE_UNARY_EXPR(exp_op, padicxx,
        detail::padic_check(
            padic_exp(to._padic(), from._padic(), to._ctx()), "exp"))
FLINT_DEFINE_UNARY_EXPR(log_op, padicxx,
        detail::padic_check(
            padic_log(to._padic(), from._padic(), to._ctx()), "log"))
// TODO some more
} // rules
// TODO non-lazy functions
} // flint

#endif
