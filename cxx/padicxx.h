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
#include "cxx/flint_classes.h"
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

    FLINTXX_DEFINE_BASICS(padicxx_expression)
    FLINTXX_DEFINE_CTORS(padicxx_expression)
    FLINTXX_DEFINE_C_REF(padicxx_expression, padic_struct, _padic)

public:
    // These only make sense with immediates
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
};

namespace detail {
struct padic_data;
}

typedef padicxx_expression<operations::immediate, detail::padic_data> padicxx;
typedef padicxx_expression<operations::immediate,
            flint_classes::ref_data<padicxx, padic_struct> > padicxx_ref;
typedef padicxx_expression<operations::immediate,
            flint_classes::srcref_data<padicxx, padicxx_ref, padic_struct> > padicxx_srcref;

namespace flint_classes {
template<class Padic>
struct ref_data<Padic, padic_struct>
{
    typedef void IS_REF_OR_CREF;
    typedef Padic wrapped_t;

    typedef padic_struct* data_ref_t;
    typedef const padic_struct* data_srcref_t;

    padic_struct* inner;
    const padicxx_ctx& ctx;

    ref_data(Padic& o) : inner(o._data().inner), ctx(o._data().ctx) {}

    static ref_data make(padic_struct* f, const padicxx_ctx& ctx)
    {
        return ref_data(f, ctx);
    }

private:
    ref_data(padic_struct* fp, const padicxx_ctx& c) : inner(fp), ctx(c) {}
};

template<class Padic, class Ref>
struct srcref_data<Padic, Ref, padic_struct>
{
    typedef void IS_REF_OR_CREF;
    typedef Padic wrapped_t;

    typedef const padic_struct* data_ref_t;
    typedef const padic_struct* data_srcref_t;

    const padic_struct* inner;
    const padicxx_ctx& ctx;

    srcref_data(const Padic& o) : inner(o._data().inner), ctx(o._data().ctx) {}
    srcref_data(Ref o) : inner(o._data().inner) {}

    static srcref_data make(const padic_struct* f, const padicxx_ctx& ctx)
    {
        return srcref_data(f, ctx);
    }

private:
    srcref_data(const padic_struct* fp, const padicxx_ctx& c) : inner(fp), ctx(c) {}
};
} // flint_classes

namespace detail {
struct padic_data
{
    typedef padic_t& data_ref_t;
    typedef const padic_t& data_srcref_t;

    const padicxx_ctx& ctx;
    padic_t inner;

    padic_data(const padicxx_ctx& c)
        : ctx(c)
    {
        padic_init(inner);
    }

    padic_data(const padicxx_ctx& c, long N)
        : ctx(c)
    {
        padic_init2(inner, N);
    }

    padic_data(const padic_data& o)
        : ctx(o.ctx)
    {
        padic_init2(inner, padic_prec(o.inner));
        padic_set(inner, o.inner, ctx._ctx());
    }

    // TODO more constructors? (e.g. unit, val?)

    ~padic_data() {padic_clear(inner);}
};
} // detail

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
#define PADICXX_COND_S FLINTXX_COND_S(padicxx)
#define PADICXX_COND_T FLINTXX_COND_T(padicxx)

FLINT_DEFINE_DOIT_COND2(assignment, PADICXX_COND_T, PADICXX_COND_S,
        padic_set(to._padic(), from._padic(), to._ctx()))
FLINT_DEFINE_DOIT_COND2(assignment, PADICXX_COND_T, traits::is_signed_integer, 
        padic_set_si(to._padic(), from, to._ctx()))
FLINT_DEFINE_DOIT_COND2(assignment, PADICXX_COND_T, traits::is_unsigned_integer, 
        padic_set_ui(to._padic(), from, to._ctx()))
FLINT_DEFINE_DOIT_COND2(assignment, PADICXX_COND_T, FMPZXX_COND_S,
        padic_set_fmpz(to._padic(), from._fmpz(), to._ctx()))
FLINT_DEFINE_DOIT_COND2(assignment, PADICXX_COND_T, FMPQXX_COND_S,
        padic_set_fmpq(to._padic(), from._fmpq(), to._ctx()))

FLINTXX_DEFINE_CONVERSION_TMP(fmpzxx, padicxx,
        padic_get_fmpz(to._fmpz(), from._padic(), from._ctx()))
FLINTXX_DEFINE_CONVERSION_TMP(fmpqxx, padicxx,
        padic_get_fmpq(to._fmpq(), from._padic(), from._ctx()))

FLINTXX_DEFINE_TO_STR(padicxx, padic_get_str(0, from._padic(), from._ctx()))

FLINTXX_DEFINE_EQUALS(padicxx, padic_equal(e1._padic(), e2._padic()))


FLINT_DEFINE_CBINARY_EXPR_COND2(plus, padicxx, PADICXX_COND_S, PADICXX_COND_S,
        padic_add(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, padicxx, PADICXX_COND_S, PADICXX_COND_S,
        padic_sub(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, padicxx, PADICXX_COND_S, PADICXX_COND_S,
        padic_mul(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, padicxx, PADICXX_COND_S, PADICXX_COND_S,
        padic_div(to._padic(), e1._padic(), e2._padic(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND2(shift, padicxx, PADICXX_COND_S,
        traits::fits_into_slong,
        padic_shift(to._padic(), e1._padic(), e2, to._ctx()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, padicxx, PADICXX_COND_S,
        padic_neg(to._padic(), from._padic(), to._ctx()))

// lazy functions
FLINT_DEFINE_UNARY_EXPR_COND(sqrt_op, padicxx, PADICXX_COND_S,
        detail::padic_check(
            padic_sqrt(to._padic(), from._padic(), to._ctx()), "sqrt"))
FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, padicxx, PADICXX_COND_S,
        traits::fits_into_slong,
        padic_pow_si(to._padic(), e1._padic(), e2, to._ctx()))
FLINT_DEFINE_UNARY_EXPR_COND(exp_op, padicxx, PADICXX_COND_S,
        detail::padic_check(
            padic_exp(to._padic(), from._padic(), to._ctx()), "exp"))
FLINT_DEFINE_UNARY_EXPR_COND(log_op, padicxx, PADICXX_COND_S,
        detail::padic_check(
            padic_log(to._padic(), from._padic(), to._ctx()), "log"))
// TODO some more
} // rules
// TODO non-lazy functions
} // flint

#endif
