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

#ifndef PADIC_MATXX_H
#define PADIC_MATXX_H

#include "padic_mat.h"

#include "padicxx.h"
#include "fmpq_matxx.h"

#include "flintxx/matrix.h"

// TODO input and output

namespace flint {
FLINT_DEFINE_THREEARY(padic_matxx_get_entry)

namespace detail {
template<class Padic>
struct padic_mat_traits
    : matrices::generic_traits<Padic>
{
    typedef slong prec_ref_t;
    typedef slong prec_srcref_t;

    typedef slong val_ref_t;
    typedef slong val_srcref_t;

    static slong prec(const Padic& p) {return tools::padic_output_prec(p);}
    static slong val(const Padic& p) {return padic_mat_val(p.evaluate()._mat());}
};
} //detail

template<class Operation, class Data>
class padic_matxx_expression
    : public expression<derived_wrapper<padic_matxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::padic_matxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(padic_matxx_expression)
    FLINTXX_DEFINE_CTORS(padic_matxx_expression)
    FLINTXX_DEFINE_C_REF(padic_matxx_expression, padic_mat_struct, _mat)

public:
    typedef detail::padic_mat_traits<padic_matxx_expression> traits_t;

    PADICXX_DEFINE_STD

    static padic_matxx_expression zero(padicxx_ctx_srcref ctx,
            slong rows, slong cols)
        {return padic_matxx_expression(ctx, rows, cols);}
    static padic_matxx_expression zero(padicxx_ctx_srcref ctx,
            slong rows, slong cols, slong N)
        {return padic_matxx_expression(ctx, rows, cols, N);}

    static padic_matxx_expression one(padicxx_ctx_srcref ctx,
            slong rows, slong cols)
    {
        padic_matxx_expression res(ctx, rows, cols);
        res.set_one();
        return res;
    }
    static padic_matxx_expression one(padicxx_ctx_srcref ctx,
            slong rows, slong cols, slong N)
    {
        padic_matxx_expression res(ctx, rows, cols, N);
        res.set_one();
        return res;
    }

    template<class T>
    static padic_matxx_expression from_QQ(const T& q, padicxx_ctx_srcref ctx,
            typename mp::enable_if<traits::is_fmpq_matxx<T> >::type* = 0)
    {
        padic_matxx_expression res(ctx, q.rows(), q.cols());
        res = q;
        return res;
    }
    template<class T>
    static padic_matxx_expression from_QQ(const T& q, padicxx_ctx_srcref ctx,
            slong N,
            typename mp::enable_if<traits::is_fmpq_matxx<T> >::type* = 0)
    {
        padic_matxx_expression res(ctx, q.rows(), q.cols(), N);
        res = q;
        return res;
    }

    template<class Expr>
    static evaluated_t create_temporary_rowscols(
            const Expr& e, slong rows, slong cols)
    {
        return evaluated_t(tools::find_padicxx_ctx(e),
                rows, cols, tools::padic_output_prec(e));
    }
    FLINTXX_DEFINE_MATRIX_METHODS(traits_t)
        
    // static methods which only make sense with padicxx
    static padic_matxx_expression randtest(slong rows, slong cols,
            frandxx& state,
            padicxx_ctx_srcref ctx, slong prec = PADIC_DEFAULT_PREC)
    {
        padic_matxx_expression res(ctx, rows, cols, prec);
        padic_mat_randtest(res._mat(), state._data(), ctx._ctx());
        return res;
    }

    // These only make sense with immediates
    void reduce() {padic_mat_reduce(_mat(), _ctx());}
    void set_zero() {padic_mat_zero(_mat());}
    void set_one() {padic_mat_one(_mat());}
    void truncate(slong n) {fmpz_poly_truncate(_mat(), n);}
    void canonicalise() {padic_mat_canonicalise(_mat());}
    bool is_canonical() const {return padic_mat_is_canonical(_mat());}
    bool is_reduced() const {return padic_mat_is_reduced(_mat());}

    template<class Padic>
    void set_entry(slong i, slong j, const Padic& p, 
            typename mp::enable_if<traits::is_padicxx<Padic> >::type* = 0)
    {
        padic_mat_set_entry_padic(_mat(), i, j, p.evaluate()._padic(), _ctx());
    }

    // these cause evaluation
    bool is_zero() const {return padic_mat_is_zero(this->evaluate()._mat());}
    bool is_empty() const {return padic_mat_is_empty(this->evaluate()._mat());}
    bool is_square() const {return padic_mat_is_square(this->evaluate()._mat());}

    // forwarding of lazy functions
    FLINTXX_DEFINE_MEMBER_UNOP(transpose)
    FLINTXX_DEFINE_MEMBER_3OP_(get_entry, padic_matxx_get_entry)
};

namespace detail {
struct padic_mat_data;
}

typedef padic_matxx_expression<operations::immediate,
            detail::padic_mat_data> padic_matxx;
typedef padic_matxx_expression<operations::immediate,
            flint_classes::ref_data<padic_matxx, padic_mat_struct> >
                padic_matxx_ref;
typedef padic_matxx_expression<operations::immediate, flint_classes::srcref_data<
    padic_matxx, padic_matxx_ref, padic_mat_struct> > padic_matxx_srcref;

template<>
struct matrix_traits<padic_matxx>
{
    template<class M> static slong rows(const M& m)
    {
        return padic_mat_nrows(m._mat());
    }
    template<class M> static slong cols(const M& m)
    {
        return padic_mat_ncols(m._mat());
    }

    template<class M> static fmpzxx_srcref at(const M& m, slong i, slong j)
    {
        return fmpzxx_srcref::make(padic_mat_entry(m._mat(), i, j));
    }
    template<class M> static fmpzxx_ref at(M& m, slong i, slong j)
    {
        return fmpzxx_ref::make(padic_mat_entry(m._mat(), i, j));
    }
};

namespace traits {
template<> struct has_padicxx_ctx<padic_matxx> : mp::true_ { };
template<> struct has_padicxx_ctx<padic_matxx_ref> : mp::true_ { };
template<> struct has_padicxx_ctx<padic_matxx_srcref> : mp::true_ { };
} // traits

namespace detail {
struct padic_matxx_srcref_traits_no_std_matrix
{
    typedef slong prec_ref_t;
    typedef slong prec_srcref_t;

    typedef slong val_ref_t;
    typedef slong val_srcref_t;

    template<class P>
    static slong prec(P p) {return p._data().N;}
    template<class P>
    static slong val(P p) {return padic_mat_val(p._mat());}
};
struct padic_matxx_ref_traits_no_std_matrix
    : padic_matxx_srcref_traits_no_std_matrix
{
    typedef slong& prec_ref_t;
    typedef slong& val_ref_t;

    template<class P>
    static slong& prec(P& p) {return padic_mat_prec(p._mat());}
    template<class P>
    static slong prec(const P& p) {return padic_mat_prec(p._mat());}

    template<class P>
    static slong& val(P& p) {return padic_mat_val(p._mat());}
    template<class P>
    static slong val(const P& p) {return padic_mat_val(p._mat());}
};
template<>
struct padic_mat_traits<padic_matxx_srcref>
    : matrices::generic_traits_srcref<fmpzxx_srcref>,
      padic_matxx_srcref_traits_no_std_matrix { };
template<>
struct padic_mat_traits<padic_matxx_ref>
    : matrices::generic_traits_ref<fmpzxx_ref>, 
      padic_matxx_ref_traits_no_std_matrix { };
template<>
struct padic_mat_traits<padic_matxx>
    : matrices::generic_traits_nonref<fmpzxx_ref, fmpzxx_srcref>, 
      padic_matxx_ref_traits_no_std_matrix { };
} // detail

PADICXX_DEFINE_REF_STRUCTS(padic_matxx, padic_mat_struct, padic_mat_prec)

namespace detail {
struct padic_mat_data
{
    typedef padic_mat_t& data_ref_t;
    typedef const padic_mat_t& data_srcref_t;

    padicxx_ctx_srcref ctx;
    padic_mat_t inner;

    padic_mat_data(padicxx_ctx_srcref c, slong rows, slong cols)
        : ctx(c)
    {
        padic_mat_init(inner, rows, cols);
    }

    padic_mat_data(padicxx_ctx_srcref c, slong rows, slong cols, slong N)
        : ctx(c)
    {
        padic_mat_init2(inner, rows, cols, N);
    }

    padic_mat_data(const padic_mat_data& o)
        : ctx(o.ctx)
    {
        padic_mat_init2(inner, padic_mat_nrows(o.inner),
                padic_mat_ncols(o.inner), padic_mat_prec(o.inner));
        padic_mat_set(inner, o.inner, ctx._ctx());
    }

    ~padic_mat_data() {padic_mat_clear(inner);}

    padic_mat_data(padic_matxx_srcref c)
        : ctx(c.get_ctx())
    {
        padic_mat_init2(inner, c.rows(), c.cols(), c.prec());
        padic_mat_set(inner, c._mat(), ctx._ctx());
    }
};
} // detail

// matrix temporary stuff
FLINTXX_DEFINE_TEMPORARY_RULES(padic_matxx)

#define PADIC_MATXX_COND_S FLINTXX_COND_S(padic_matxx)
#define PADIC_MATXX_COND_T FLINTXX_COND_T(padic_matxx)

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, PADIC_MATXX_COND_T, PADIC_MATXX_COND_S,
        padic_mat_set(to._mat(), from._mat(), to._ctx()))
FLINT_DEFINE_DOIT_COND2(assignment, PADIC_MATXX_COND_T, FMPQ_MATXX_COND_S,
        padic_mat_set_fmpq_mat(to._mat(), from._mat(), to._ctx()))

FLINTXX_DEFINE_SWAP(padic_matxx, padic_mat_swap(e1._mat(), e2._mat()))

FLINTXX_DEFINE_EQUALS(padic_matxx, padic_mat_equal(e1._mat(), e2._mat()))

FLINT_DEFINE_PRINT_COND(PADIC_MATXX_COND_S,
        padic_mat_fprint(to, from._mat(), from._ctx()))
FLINT_DEFINE_PRINT_PRETTY_COND(PADIC_MATXX_COND_S,
        padic_mat_fprint_pretty(to, from._mat(), from._ctx()))

template<class T>
struct conversion<fmpq_matxx, T,
    typename mp::enable_if< PADIC_MATXX_COND_S<T> >::type>
{
    static fmpq_matxx get(const T& from)
    {
        fmpq_matxx to(from.rows(), from.cols());
        padic_mat_get_fmpq_mat(to._mat(), from._mat(), from._ctx());
        return to;
    }
};

FLINT_DEFINE_THREEARY_EXPR_COND3(padic_matxx_get_entry_op, padicxx,
        PADIC_MATXX_COND_S, traits::fits_into_slong, traits::fits_into_slong,
        padic_mat_get_entry_padic(to._padic(), e1._mat(), e2, e3, to._ctx()))

FLINT_DEFINE_UNARY_EXPR_COND(transpose_op, padic_matxx, PADIC_MATXX_COND_S,
        padic_mat_transpose(to._mat(), from._mat()))

FLINT_DEFINE_CBINARY_EXPR_COND2(plus, padic_matxx,
        PADIC_MATXX_COND_S, PADIC_MATXX_COND_S,
        padic_mat_add(to._mat(), e1._mat(), e2._mat(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, padic_matxx,
        PADIC_MATXX_COND_S, PADIC_MATXX_COND_S,
        padic_mat_sub(to._mat(), e1._mat(), e2._mat(), to._ctx()))
FLINT_DEFINE_UNARY_EXPR_COND(negate, padic_matxx, PADIC_MATXX_COND_S,
        padic_mat_neg(to._mat(), from._mat(), to._ctx()))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, padic_matxx,
        PADIC_MATXX_COND_S, PADICXX_COND_S,
        padic_mat_scalar_mul_padic(to._mat(), e1._mat(), e2._padic(), to._ctx()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, padic_matxx,
        PADIC_MATXX_COND_S, FMPZXX_COND_S,
        padic_mat_scalar_mul_fmpz(to._mat(), e1._mat(), e2._fmpz(), to._ctx()))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, padic_matxx,
        PADIC_MATXX_COND_S, FMPZXX_COND_S,
        padic_mat_scalar_div_fmpz(to._mat(), e1._mat(), e2._fmpz(), to._ctx()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, padic_matxx,
        PADIC_MATXX_COND_S, PADIC_MATXX_COND_S,
        padic_mat_mul(to._mat(), e1._mat(), e2._mat(), to._ctx()))
} // rules
}

#endif
