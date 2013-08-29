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

#ifndef PADIC_POLYXX_H
#define PADIC_POLYXX_H

#include "padic_poly.h"

#include "padicxx.h"

namespace flint {
FLINT_DEFINE_BINOP(padic_polyxx_get_coeff)

namespace detail {
template<class Padic>
struct padic_poly_traits
{
    typedef slong prec_ref_t;
    typedef slong prec_srcref_t;

    typedef slong val_ref_t;
    typedef slong val_srcref_t;

    static slong prec(const Padic& p) {return tools::padic_output_prec(p);}
    static slong val(const Padic& p) {return padic_poly_val(p.evaluate()._padic());}
};
} //detail

template<class Operation, class Data>
class padic_polyxx_expression
    : public expression<derived_wrapper<padic_polyxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::padic_polyxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(padic_polyxx_expression)
    FLINTXX_DEFINE_CTORS(padic_polyxx_expression)
    FLINTXX_DEFINE_C_REF(padic_polyxx_expression, padic_poly_struct, _poly)

public:
    typedef detail::padic_poly_traits<padic_polyxx_expression> traits_t;

    PADICXX_DEFINE_STD

    static padic_polyxx_expression zero(const padicxx_ctx& ctx)
        {return padic_polyxx_expression(ctx);}
    static padic_polyxx_expression zero(const padicxx_ctx& ctx, long N)
        {return padic_polyxx_expression(ctx, N);}
    static padic_polyxx_expression one(const padicxx_ctx& ctx)
    {
        padic_polyxx_expression res(ctx);
        res.set_one();
        return res;
    }
    static padic_polyxx_expression one(const padicxx_ctx& ctx, long N)
    {
        padic_polyxx_expression res(ctx, N);
        res.set_one();
        return res;
    }

    template<class T>
    static padic_polyxx_expression from_QQ(const T& q, const padicxx_ctx& ctx,
            typename mp::enable_if<mp::or_<traits::is_fmpqxx<T>,
                traits::is_fmpzxx<T>, traits::is_integer<T> > >::type* = 0)
    {
        padic_polyxx_expression res(ctx);
        res = q;
        return res;
    }
    template<class T>
    static padic_polyxx_expression from_QQ(const T& q, const padicxx_ctx& ctx,
            slong N,
            typename mp::enable_if<mp::or_<traits::is_fmpqxx<T>,
                traits::is_fmpzxx<T>, traits::is_integer<T> > >::type* = 0)
    {
        padic_polyxx_expression res(ctx, N);
        res = q;
        return res;
    }
    // TODO more?

    // Create a temporary. The context will be estimated, and the precision
    // will be the maximum of all subexpressions.
    evaluated_t create_temporary() const
    {
        return evaluated_t(estimate_ctx(), prec());
    }

#if 0
    // static methods which only make sense with padicxx
    static padic_polyxx_expression randtest(frandxx& state,
            const padicxx_ctx& ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        padic_polyxx_expression res(ctx, prec);
        padic_randtest(res._padic(), state._data(), ctx._ctx());
        return res;
    }
    static padic_polyxx_expression randtest_not_zero(frandxx& state,
            const padicxx_ctx& ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        padic_polyxx_expression res(ctx, prec);
        padic_randtest_not_zero(res._padic(), state._data(), ctx._ctx());
        return res;
    }
    static padic_polyxx_expression randtest_int(frandxx& state,
            const padicxx_ctx& ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        padic_polyxx_expression res(ctx, prec);
        padic_randtest_int(res._padic(), state._data(), ctx._ctx());
        return res;
    }
#endif

    // These only make sense with immediates
    void reduce() {padic_poly_reduce(_poly(), _ctx());}
    void realloc(slong alloc) {padic_poly_realloc(_poly(), alloc);}
    void fit_length(slong len) {padic_poly_fit_length(_poly(), len);}
    void _normalise() {_padic_poly_normalise(_poly());}
    void _set_length(slong len) {_padic_poly_set_length(_poly(), len);}
    void set_zero() {padic_poly_zero(_poly());}
    void set_one() {padic_poly_one(_poly());}
    void truncate(slong n) {fmpz_poly_truncate(_poly(), n);}

    template<class Padic>
    void set_coeff(slong n, const Padic& p,
            typename mp::enable_if<traits::is_padicxx<Padic> >::type* = 0)
        {padic_poly_set_coeff_padic(_poly(), n, p._padic());}

    // these cause evaluation
    bool is_zero() const {return padic_poly_is_zero(this->evaluate()._padic());}
    bool is_one() const {return padic_poly_is_one(this->evaluate()._padic());}
    slong length() const {return padic_poly_length(this->evaluate()._poly());}
    slong degree() const {return padic_poly_degree(this->evaluate()._poly());}

    // forwarding of lazy functions
    FLINTXX_DEFINE_MEMBER_BINOP_(get_coeff, padic_polyxx_get_coeff)
};

namespace detail {
struct padic_poly_data;
}

typedef padic_polyxx_expression<operations::immediate,
            detail::padic_poly_data> padic_polyxx;
typedef padic_polyxx_expression<operations::immediate,
            flint_classes::ref_data<padic_polyxx, padic_poly_struct> >
                padic_polyxx_ref;
typedef padic_polyxx_expression<operations::immediate, flint_classes::srcref_data<
    padic_polyxx, padic_polyxx_ref, padic_poly_struct> > padic_polyxx_srcref;

namespace traits {
template<> struct has_padicxx_ctx<padic_polyxx> : mp::true_ { };
template<> struct has_padicxx_ctx<padic_polyxx_ref> : mp::true_ { };
template<> struct has_padicxx_ctx<padic_polyxx_srcref> : mp::true_ { };
} // traits

namespace detail {
template<>
struct padic_poly_traits<padic_polyxx_srcref>
{
    typedef slong prec_ref_t;
    typedef slong prec_srcref_t;

    typedef slong val_ref_t;
    typedef slong val_srcref_t;

    template<class P>
    static slong prec(P p) {return p._data().N;}
    template<class P>
    static slong val(P p) {return padic_poly_val(p._padic_poly());}
};

template<>
struct padic_poly_traits<padic_polyxx_ref>
    : padic_poly_traits<padic_polyxx_srcref>
{
    typedef slong& prec_ref_t;
    typedef slong& val_ref_t;

    template<class P>
    static slong& prec(P& p) {return padic_poly_prec(p._padic_poly());}
    template<class P>
    static slong prec(const P& p) {return padic_poly_prec(p._padic_poly());}

    template<class P>
    static slong& val(P& p) {return padic_poly_val(p._padic_poly());}
    template<class P>
    static slong val(const P& p) {return padic_poly_val(p._padic_poly());}
};
template<>
struct padic_poly_traits<padic_polyxx>
    : padic_poly_traits<padic_polyxx_ref> { };
} // detail

PADICXX_DEFINE_REF_STRUCTS(padic_poly_struct, padic_poly_prec)

namespace detail {
struct padic_poly_data
{
    typedef padic_poly_t& data_ref_t;
    typedef const padic_poly_t& data_srcref_t;

    const padicxx_ctx& ctx;
    padic_poly_t inner;

    padic_poly_data(const padicxx_ctx& c)
        : ctx(c)
    {
        padic_poly_init(inner);
    }

    padic_poly_data(const padicxx_ctx& c, slong N, slong alloc = 0)
        : ctx(c)
    {
        padic_poly_init2(inner, alloc, N);
    }

    padic_poly_data(const padic_poly_data& o)
        : ctx(o.ctx)
    {
        padic_poly_init2(inner, padic_poly_length(o.inner),
                padic_poly_prec(o.inner));
        padic_poly_set(inner, o.inner, ctx._ctx());
    }

    ~padic_poly_data() {padic_poly_clear(inner);}

    padic_poly_data(padic_polyxx_srcref c)
        : ctx(c.get_ctx())
    {
        padic_poly_init2(inner, c.length(), c.prec());
        padic_poly_set(inner, c._poly(), ctx._ctx());
    }
};
} // detail

#define PADIC_POLYXX_COND_S FLINTXX_COND_S(padic_polyxx)
#define PADIC_POLYXX_COND_T FLINTXX_COND_T(padic_polyxx)

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, PADIC_POLYXX_COND_T, PADIC_POLYXX_COND_S,
        padic_poly_set(to._poly(), from._poly(), to._ctx()))

FLINTXX_DEFINE_SWAP(padic_polyxx, padic_poly_swap(e1._poly(), e2._poly()))

FLINTXX_DEFINE_EQUALS(padic_polyxx, padic_poly_equal(e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(padic_polyxx_get_coeff_op, padicxx,
        PADIC_POLYXX_COND_S, traits::fits_into_slong,
        padic_poly_get_coeff_padic(to._padic(), e1._poly(), e2, to._ctx()))
} // rules

} // flint

#endif
