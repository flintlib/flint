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

#include "padic.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/flint_exception.h"
#include "flintxx/frandxx.h"
#include "flintxx/stdmath.h"
#include "flintxx/traits.h"
#include "flintxx/tuple.h"

#include "fmpzxx.h"
#include "fmpqxx.h"

// TODO check codegen ...
// TODO padic_output_prec does not work on non-padic expressions,
//      is that a problem?

namespace flint {
// function "declarations"
FLINT_DEFINE_UNOP(exp_rectangular)
FLINT_DEFINE_UNOP(exp_balanced)
FLINT_DEFINE_UNOP(log_rectangular)
FLINT_DEFINE_UNOP(log_balanced)
FLINT_DEFINE_UNOP(log_satoh)
FLINT_DEFINE_UNOP(teichmuller)
FLINT_DEFINE_BINOP(padic_val_fac)

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

    padic_print_mode mode() const {return _ctx()->mode;}
    padic_print_mode& mode() {return _ctx()->mode;}

    // TODO more constructors? Should we wrap padic_print_mode?
    padicxx_ctx(fmpzxx_srcref p, slong min, slong max, padic_print_mode mode)
    {
        padic_ctx_init(ctx, p._fmpz(), min, max, mode);
    }

    ~padicxx_ctx() {padic_ctx_clear(ctx);}
};

class padicxx_ctx_srcref
{
private:
    mutable padic_ctx_struct* ctx;

    padicxx_ctx_srcref(padic_ctx_struct* c) : ctx(c) {}

public:
    // NB: you must not modify user-visible state of ctx through a constant
    // instance of padicxx_ctx
    padic_ctx_struct* _ctx() const {return ctx;}

    // XXX these two are not actually exposed in the C api ...
    fmpzxx_ref get_p() {return fmpzxx_ref::make(ctx[0].p);}
    fmpzxx_srcref get_p() const {return fmpzxx_srcref::make(ctx[0].p);}

    padic_print_mode mode() const {return _ctx()->mode;}

    padicxx_ctx_srcref(padicxx_ctx& c)
        : ctx(c._ctx()) {}

    static padicxx_ctx_srcref make(padic_ctx_struct* c)
        {return padicxx_ctx_srcref(c);}
};

namespace traits {
template<class T> struct has_padicxx_ctx : mp::false_ { };

template<class T> struct is_padic_expr
    : has_padicxx_ctx<typename tools::evaluation_helper<T>::type> { };
} // traits
namespace detail {
struct has_padicxx_ctx_predicate
{
    template<class T> struct type : traits::has_padicxx_ctx<T> { };
};

template<class T, class Enable = void>
struct padicxx_max_prec;

template<class T>
struct padicxx_max_prec<T,
    typename mp::enable_if<traits::has_padicxx_ctx<T> >::type>
{
    static slong get(const T& p) {return p.prec();}
};

template<class T>
struct padicxx_max_prec<T,
    typename mp::disable_if<traits::is_padic_expr<T> >::type>
{
    static slong get(const T&) {return 0;}
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

template<class T>
struct padicxx_max_prec<T, typename mp::enable_if<mp::and_<
    traits::is_padic_expr<T>,
    mp::not_<traits::is_immediate<T> > > >::type>
{
    static slong get(const T& e)
        {return padicxx_max_prec_h<typename T::data_t>::get(e._data());}
};
} // detail
namespace tools {
template<class Expr>
padicxx_ctx_srcref find_padicxx_ctx(const Expr& e)
{
    return tools::find_subexpr<detail::has_padicxx_ctx_predicate>(e).get_ctx();
}

template<class Expr>
slong padic_output_prec(const Expr& e)
{
    return detail::padicxx_max_prec<Expr>::get(e);
}
} // tools

FLINT_DEFINE_UNOP(padicxx_unit)

namespace detail {
template<class Padic>
struct padic_traits
{
    typedef FLINT_UNOP_BUILD_RETTYPE(padicxx_unit, fmpzxx, Padic)
        unit_srcref_t;

    typedef slong prec_ref_t;
    typedef slong prec_srcref_t;

    typedef slong val_ref_t;
    typedef slong val_srcref_t;

    static slong prec(const Padic& p) {return tools::padic_output_prec(p);}
    static slong val(const Padic& p) {return padic_val(p.evaluate()._padic());}

    static unit_srcref_t unit(const Padic& p) {return padicxx_unit(p);}
};
} //detail

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
    typedef detail::padic_traits<padicxx_expression> traits_t;

    // These only make sense with immediates
    void reduce() {padic_reduce(_padic(), _ctx());}
    void set_zero() {padic_zero(_padic());}
    void set_one() {padic_one(_padic());}

#define PADICXX_DEFINE_CTX \
    padicxx_ctx_srcref get_ctx() const {return this->_data().ctx;} \
    padic_ctx_struct* _ctx() const {return get_ctx()._ctx();}
    PADICXX_DEFINE_CTX

    static padicxx_expression zero(padicxx_ctx_srcref ctx)
        {return padicxx_expression(ctx);}
    static padicxx_expression zero(padicxx_ctx_srcref ctx, slong N)
        {return padicxx_expression(ctx, N);}
    static padicxx_expression one(padicxx_ctx_srcref ctx)
    {
        padicxx_expression res(ctx);
        res.set_one();
        return res;
    }
    static padicxx_expression one(padicxx_ctx_srcref ctx, slong N)
    {
        padicxx_expression res(ctx, N);
        res.set_one();
        return res;
    }

    template<class T>
    static padicxx_expression from_QQ(const T& q, padicxx_ctx_srcref ctx,
            typename mp::enable_if<mp::or_<traits::is_fmpqxx<T>,
                traits::is_fmpzxx<T>, traits::is_integer<T> > >::type* = 0)
    {
        padicxx_expression res(ctx);
        res = q;
        return res;
    }
    template<class T>
    static padicxx_expression from_QQ(const T& q, padicxx_ctx_srcref ctx,
            slong N,
            typename mp::enable_if<mp::or_<traits::is_fmpqxx<T>,
                traits::is_fmpzxx<T>, traits::is_integer<T> > >::type* = 0)
    {
        padicxx_expression res(ctx, N);
        res = q;
        return res;
    }
    // TODO more?

    // The above method get_ctx() only works on immediates, i.e. instances of
    // padicxx. This next one here works on any composite expression which
    // contains at least one instance of padicxx. Returns the context of one of
    // those immediate subexpressions.
#define PADICXX_DEFINE_ESTIMATE_CTX \
    padicxx_ctx_srcref estimate_ctx() const \
    { \
        return tools::find_padicxx_ctx(*this); \
    }
    PADICXX_DEFINE_ESTIMATE_CTX

    // Create a temporary. The context will be estimated, and the precision
    // will be the maximum of all subexpressions.
    evaluated_t create_temporary() const
    {
        return evaluated_t(estimate_ctx(), prec());
    }

    // static methods which only make sense with padicxx
    static padicxx_expression randtest(frandxx& state,
            padicxx_ctx_srcref ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        padicxx_expression res(ctx, prec);
        padic_randtest(res._padic(), state._data(), ctx._ctx());
        return res;
    }
    static padicxx_expression randtest_not_zero(frandxx& state,
            padicxx_ctx_srcref ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        padicxx_expression res(ctx, prec);
        padic_randtest_not_zero(res._padic(), state._data(), ctx._ctx());
        return res;
    }
    static padicxx_expression randtest_int(frandxx& state,
            padicxx_ctx_srcref ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        padicxx_expression res(ctx, prec);
        padic_randtest_int(res._padic(), state._data(), ctx._ctx());
        return res;
    }

#define PADICXX_DEFINE_TON \
    typename flint_classes::to_srcref<typename base_t::derived_t>::type \
    toN(slong N) const \
    { \
        return flint_classes::to_srcref<typename base_t::derived_t>::type::make( \
                this->_data().inner, get_ctx(), N); \
    }
    PADICXX_DEFINE_TON

    typename traits_t::unit_srcref_t unit() const
        {return traits_t::unit(*this);}

    // Compute the maximal precision of all subexpressions
#define PADICXX_DEFINE_PREC \
    typename traits_t::prec_ref_t prec() {return traits_t::prec(*this);} \
    typename traits_t::prec_srcref_t prec() const \
        {return traits_t::prec(*this);}
    PADICXX_DEFINE_PREC

#define PADICXX_DEFINE_VAL \
    typename traits_t::val_ref_t val() {return traits_t::val(*this);} \
    typename traits_t::val_srcref_t val() const \
        {return traits_t::val(*this);}
    PADICXX_DEFINE_VAL

    // these cause evaluation
    bool is_zero() const {return padic_is_zero(this->evaluate()._padic());}
    bool is_one() const {return padic_is_one(this->evaluate()._padic());}

    // forwarding of lazy functions
    FLINTXX_DEFINE_MEMBER_UNOP(exp)
    FLINTXX_DEFINE_MEMBER_UNOP(exp_balanced)
    FLINTXX_DEFINE_MEMBER_UNOP(exp_rectangular)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_UNOP(log)
    FLINTXX_DEFINE_MEMBER_UNOP(log_balanced)
    FLINTXX_DEFINE_MEMBER_UNOP(log_satoh)
    FLINTXX_DEFINE_MEMBER_UNOP(sqrt)
    FLINTXX_DEFINE_MEMBER_UNOP(teichmuller)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
};

#define PADICXX_DEFINE_STD \
    PADICXX_DEFINE_CTX \
    PADICXX_DEFINE_ESTIMATE_CTX \
    PADICXX_DEFINE_TON \
    PADICXX_DEFINE_PREC \
    PADICXX_DEFINE_VAL

namespace detail {
struct padic_data;
}

typedef padicxx_expression<operations::immediate, detail::padic_data> padicxx;
typedef padicxx_expression<operations::immediate,
            flint_classes::ref_data<padicxx, padic_struct> > padicxx_ref;
typedef padicxx_expression<operations::immediate, flint_classes::srcref_data<
    padicxx, padicxx_ref, padic_struct> > padicxx_srcref;

namespace traits {
template<> struct has_padicxx_ctx<padicxx> : mp::true_ { };
template<> struct has_padicxx_ctx<padicxx_ref> : mp::true_ { };
template<> struct has_padicxx_ctx<padicxx_srcref> : mp::true_ { };

template<class T> struct is_padicxx : flint_classes::is_Base<padicxx, T> { };
} // traits

namespace detail {
template<>
struct padic_traits<padicxx_srcref>
{
    typedef fmpzxx_srcref unit_srcref_t;
    template<class Poly>
    static fmpzxx_srcref unit(const Poly& p)
        {return fmpzxx_srcref::make(padic_unit(p._padic()));}

    typedef slong prec_ref_t;
    typedef slong prec_srcref_t;

    typedef slong val_ref_t;
    typedef slong val_srcref_t;

    template<class P>
    static slong prec(P p) {return p._data().N;}
    template<class P>
    static slong val(P p) {return padic_val(p._padic());}
};

template<>
struct padic_traits<padicxx_ref>
    : padic_traits<padicxx_srcref>
{
    typedef slong& prec_ref_t;
    typedef slong& val_ref_t;

    template<class P>
    static slong& prec(P& p) {return padic_prec(p._padic());}
    template<class P>
    static slong prec(const P& p) {return padic_prec(p._padic());}

    template<class P>
    static slong& val(P& p) {return padic_val(p._padic());}
    template<class P>
    static slong val(const P& p) {return padic_val(p._padic());}
};
template<>
struct padic_traits<padicxx>
    : padic_traits<padicxx_ref> { };

template<class T, class U>
struct faketemplate : mp::enable_if<mp::equal_types<T, U> > { };
} // detail

// NB: usually, the "padicname" parameter would not be necessary. We would
// leave that as a template, identifying the type by structname alone, and
// conveniently delay all instantiations. Unfortunately qadic and padic_poly
// have the same structname, so we cannot do this.
// Instead we pass in the class name explicitly, and delay relevant functions
// by hand...
#define PADICXX_DEFINE_REF_STRUCTS_(padicname, structname, precname, ctxtype) \
namespace flint_classes {                                                     \
template<>                                                                    \
struct ref_data<padicname, structname>                                        \
{                                                                             \
    typedef void IS_REF_OR_CREF;                                              \
    typedef padicname wrapped_t;                                              \
                                                                              \
    typedef structname* data_ref_t;                                           \
    typedef const structname* data_srcref_t;                                  \
                                                                              \
    structname* inner;                                                        \
    ctxtype ctx;                                                              \
                                                                              \
    template<class T>                                                         \
    ref_data(T& o, typename detail::faketemplate<T, padicname>::type* = 0)    \
        : inner(o._data().inner), ctx(o._data().ctx) {}                       \
                                                                              \
    static ref_data make(structname* f, ctxtype ctx)                          \
    {                                                                         \
        return ref_data(f, ctx);                                              \
    }                                                                         \
                                                                              \
private:                                                                      \
    ref_data(structname* fp, ctxtype c) : inner(fp), ctx(c) {}                \
};                                                                            \
                                                                              \
template<class Ref>                                                           \
struct srcref_data<padicname, Ref, structname>                                \
{                                                                             \
    typedef void IS_REF_OR_CREF;                                              \
    typedef padicname wrapped_t;                                              \
                                                                              \
    typedef const structname* data_ref_t;                                     \
    typedef const structname* data_srcref_t;                                  \
                                                                              \
    const structname* inner;                                                  \
    ctxtype ctx;                                                              \
    slong N;                                                                  \
                                                                              \
    template<class T>                                                         \
    srcref_data(const T& o,                                                   \
            typename detail::faketemplate<T, padicname>::type* = 0)           \
        : inner(o._data().inner), ctx(o._data().ctx), N(o.prec()) {}          \
    template<class T>                                                         \
    srcref_data(T o, typename detail::faketemplate<T, Ref>::type* = 0)        \
        : inner(o._data().inner), ctx(o._data().ctx), N(o.prec()) {}          \
                                                                              \
    static srcref_data make(const structname* f, ctxtype ctx)                 \
    {                                                                         \
        return srcref_data(f, ctx);                                           \
    }                                                                         \
    static srcref_data make(const structname* f, ctxtype ctx,                 \
            slong N)                                                          \
    {                                                                         \
        return srcref_data(f, ctx, N);                                        \
    }                                                                         \
                                                                              \
private:                                                                      \
    srcref_data(const structname* fp, ctxtype c)                   \
        : inner(fp), ctx(c), N(precname(fp)) {}                               \
    srcref_data(const structname* fp, ctxtype c, slong n)          \
        : inner(fp), ctx(c), N(n) {}                                          \
};                                                                            \
} /* flint_classes */
#define PADICXX_DEFINE_REF_STRUCTS(padicname, structname, precname) \
    PADICXX_DEFINE_REF_STRUCTS_(padicname, structname, precname, padicxx_ctx_srcref)
PADICXX_DEFINE_REF_STRUCTS(padicxx, padic_struct, padic_prec)

namespace detail {
struct padic_data
{
    typedef padic_t& data_ref_t;
    typedef const padic_t& data_srcref_t;

    padicxx_ctx_srcref ctx;
    padic_t inner;

    padic_data(padicxx_ctx_srcref c)
        : ctx(c)
    {
        padic_init(inner);
    }

    padic_data(padicxx_ctx_srcref c, slong N)
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

    ~padic_data() {padic_clear(inner);}

    padic_data(padicxx_srcref c)
        : ctx(c.get_ctx())
    {
        padic_init2(inner, c.prec());
        padic_set(inner, c._padic(), ctx._ctx());
    }
};
} // detail

#define PADICXX_COND_S FLINTXX_COND_S(padicxx)
#define PADICXX_COND_T FLINTXX_COND_T(padicxx)

namespace rules {

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

FLINTXX_DEFINE_SWAP(padicxx, padic_swap(e1._padic(), e2._padic()))

FLINTXX_DEFINE_EQUALS(padicxx, padic_equal(e1._padic(), e2._padic()))

FLINT_DEFINE_UNARY_EXPR_COND(padicxx_unit_op, fmpzxx, PADICXX_COND_S,
        fmpz_set(to._fmpz(), padic_unit(from._padic())))

FLINT_DEFINE_PRINT_COND(PADICXX_COND_S,
        padic_fprint(to, from._padic(), from._ctx()))


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
        execution_check(
            padic_sqrt(to._padic(), from._padic(), to._ctx()), "sqrt", "padic"))
FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, padicxx, PADICXX_COND_S,
        traits::fits_into_slong,
        padic_pow_si(to._padic(), e1._padic(), e2, to._ctx()))
FLINT_DEFINE_UNARY_EXPR_COND(exp_op, padicxx, PADICXX_COND_S,
        execution_check(
            padic_exp(to._padic(), from._padic(), to._ctx()), "exp", "padic"))
FLINT_DEFINE_UNARY_EXPR_COND(exp_balanced_op, padicxx, PADICXX_COND_S,
        execution_check(padic_exp_balanced(
                to._padic(), from._padic(), to._ctx()), "exp_balanced", "padic"))
FLINT_DEFINE_UNARY_EXPR_COND(exp_rectangular_op, padicxx, PADICXX_COND_S,
        execution_check(padic_exp_rectangular(
                to._padic(), from._padic(), to._ctx()),
            "exp_rectangular", "padic"))
FLINT_DEFINE_UNARY_EXPR_COND(log_op, padicxx, PADICXX_COND_S,
        execution_check(
            padic_log(to._padic(), from._padic(), to._ctx()), "log", "padic"))
FLINT_DEFINE_UNARY_EXPR_COND(log_rectangular_op, padicxx, PADICXX_COND_S,
        execution_check(padic_log_rectangular(
                to._padic(), from._padic(), to._ctx()),
            "log_rectangular", "padic"))
FLINT_DEFINE_UNARY_EXPR_COND(log_balanced_op, padicxx, PADICXX_COND_S,
        execution_check(padic_log_balanced(
                to._padic(), from._padic(), to._ctx()), "log_balanced", "padic"))
FLINT_DEFINE_UNARY_EXPR_COND(log_satoh_op, padicxx, PADICXX_COND_S,
        execution_check(padic_log_satoh(to._padic(), from._padic(), to._ctx()),
            "log_satoh", "padic"))
FLINT_DEFINE_UNARY_EXPR_COND(inv_op, padicxx, PADICXX_COND_S,
            padic_inv(to._padic(), from._padic(), to._ctx()))
FLINT_DEFINE_UNARY_EXPR_COND(teichmuller_op, padicxx, PADICXX_COND_S,
            padic_teichmuller(to._padic(), from._padic(), to._ctx()))

FLINT_DEFINE_BINARY_EXPR_COND2(padic_val_fac_op, fmpzxx,
        FMPZXX_COND_S, FMPZXX_COND_S,
        ::padic_val_fac(to._fmpz(), e1._fmpz(), e2._fmpz()))
} // rules

// immediate version of padic_val_fac
template<class Fmpz, class T>
inline typename mp::enable_if<mp::and_<
        traits::is_unsigned_integer<T>, traits::is_fmpzxx<Fmpz> >,
    ulong>::type
padic_val_fac(T n, const Fmpz& p)
{
    return padic_val_fac_ui(n, p._fmpz());
}
} // flint

#endif
