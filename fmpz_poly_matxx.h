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

#ifndef FMPZ_POLY_MATXX_H
#define FMPZ_POLY_MATXX_H FMPZ_POLY_MATXX_H

#include "fmpz_poly_mat.h"

#include "fmpz_matxx.h"
#include "fmpz_polyxx.h"
#include "permxx.h"

#include "flintxx/matrix.h"

namespace flint {
FLINT_DEFINE_UNOP(prod)

namespace detail {
template<class Mat>
struct fmpz_poly_matxx_traits : matrices::generic_traits<Mat> { };
} // detail

template<class Operation, class Data>
class fmpz_poly_matxx_expression
    : public expression<derived_wrapper<fmpz_poly_matxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpz_poly_matxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpz_poly_matxx_traits<fmpz_poly_matxx_expression> traits_t;

    FLINTXX_DEFINE_BASICS(fmpz_poly_matxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_poly_matxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_poly_matxx_expression, fmpz_poly_mat_struct, _mat)

    template<class Expr>
    static evaluated_t create_temporary_rowscols(
            const Expr&, slong rows, slong cols)
    {
        return evaluated_t(rows, cols);
    }
    FLINTXX_DEFINE_MATRIX_METHODS(traits_t)

    // static functions for fmpz_poly_matxx
    template<class Fmpz_matxx>
    static fmpz_poly_matxx_expression from_ground(const Fmpz_matxx& f)
    {
        return _from_ground(f.evaluate());
    }
    template<class Fmpz_matxx>
    static fmpz_poly_matxx_expression _from_ground(const Fmpz_matxx& f)
    {
        fmpz_poly_matxx_expression res(f.rows(), f.cols());
        for(slong i = 0;i < f.rows();++i)
            for(slong j = 0;j < f.rows();++j)
                res.at(i, j).set_coeff(0, f.at(i, j));
        return res;
    }

    static fmpz_poly_matxx_expression randtest(slong rows, slong cols,
            frandxx& state, slong len, mp_bitcnt_t bits)
    {
        fmpz_poly_matxx_expression res(rows, cols);
        res.set_randtest(state, len, bits);
        return res;
    }
    static fmpz_poly_matxx_expression randtest_unsigned(slong rows, slong cols,
            frandxx& state, slong len, mp_bitcnt_t bits)
    {
        fmpz_poly_matxx_expression res(rows, cols);
        res.set_randtest_unsigned(state, len, bits);
        return res;
    }
    static fmpz_poly_matxx_expression randtest_sparse(slong rows, slong cols,
            frandxx& state, slong len, mp_bitcnt_t bits, float density)
    {
        fmpz_poly_matxx_expression res(rows, cols);
        res.set_randtest_sparse(state, len, bits, density);
        return res;
    }

    static fmpz_poly_matxx_expression zero(slong rows, slong cols)
        {return fmpz_poly_matxx_expression(rows, cols);}
    static fmpz_poly_matxx_expression one(slong rows, slong cols)
    {
        fmpz_poly_matxx_expression res(rows, cols);
        res.set_one();
        return res;
    }

    // these only make sense with targets
    void set_randtest(frandxx& state, slong len, mp_bitcnt_t bits)
        {fmpz_poly_mat_randtest(_mat(), state._data(), len, bits);}
    void set_randtest_unsigned(frandxx& state, slong len, mp_bitcnt_t bits)
        {fmpz_poly_mat_randtest_unsigned(_mat(), state._data(), len, bits);}
    void set_randtest_sparse(frandxx& state, slong len, mp_bitcnt_t bits,
            float density)
        {fmpz_poly_mat_randtest_sparse(_mat(), state._data(), len, bits, density);}
    void truncate(slong len) {fmpz_poly_mat_truncate(_mat(), len);}
    void set_zero()
        {fmpz_poly_mat_zero(_mat());}
    void set_one()
        {fmpz_poly_mat_one(_mat());}

    // these cause evaluation
    slong rank() const {return fmpz_poly_mat_rank(this->evaluate()._mat());}
    bool is_zero() const
        {return fmpz_poly_mat_is_zero(this->evaluate()._mat());}
    bool is_one() const
        {return fmpz_poly_mat_is_one(this->evaluate()._mat());}
    bool is_empty() const
        {return fmpz_poly_mat_is_empty(this->evaluate()._mat());}
    bool is_square() const
        {return fmpz_poly_mat_is_square(this->evaluate()._mat());}
    slong max_length() const
        {return fmpz_poly_mat_max_length(this->evaluate()._mat());}
    slong max_bits() const
        {return fmpz_poly_mat_max_bits(this->evaluate()._mat());}
    slong find_pivot_any(slong start, slong end, slong c) const
    {
        return fmpz_poly_mat_find_pivot_any(
                this->evaluate()._mat(), start, end, c);
    }
    slong find_pivot_partial(slong start, slong end, slong c) const
    {
        return fmpz_poly_mat_find_pivot_partial(
                this->evaluate()._mat(), start, end, c);
    }

    // forwarded lazy ops
    FLINTXX_DEFINE_MEMBER_BINOP_(operator(), compeval)

    FLINTXX_DEFINE_MEMBER_3OP(mullow)
    FLINTXX_DEFINE_MEMBER_3OP(pow_trunc)

    FLINTXX_DEFINE_MEMBER_BINOP(solve)
    FLINTXX_DEFINE_MEMBER_BINOP(solve_fflu)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_KS)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_classical)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_BINOP(sqrlow)

    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpz_polyxx, det)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpz_polyxx, det_fflu)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpz_polyxx, det_interpolate)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpz_polyxx, trace)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_classical)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_KS)
    FLINTXX_DEFINE_MEMBER_UNOP(transpose)

    //FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(nullspace) // TODO
    //FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(???, inv) // TODO
    //FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(???, rref) // TODO

    FLINTXX_DEFINE_MEMBER_FFLU
};

namespace detail {
struct fmpz_poly_mat_data;
} // detail

typedef fmpz_poly_matxx_expression<operations::immediate,
            detail::fmpz_poly_mat_data> fmpz_poly_matxx;
typedef fmpz_poly_matxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_poly_matxx,
            fmpz_poly_mat_struct> > fmpz_poly_matxx_ref;
typedef fmpz_poly_matxx_expression<operations::immediate,
            flint_classes::srcref_data<
                fmpz_poly_matxx, fmpz_poly_matxx_ref,
                fmpz_poly_mat_struct> > fmpz_poly_matxx_srcref;

template<>
struct matrix_traits<fmpz_poly_matxx>
{
    template<class M> static slong rows(const M& m)
    {
        return fmpz_poly_mat_nrows(m._mat());
    }
    template<class M> static slong cols(const M& m)
    {
        return fmpz_poly_mat_ncols(m._mat());
    }

    template<class M> static fmpz_polyxx_srcref at(const M& m, slong i, slong j)
    {
        return fmpz_polyxx_srcref::make(fmpz_poly_mat_entry(m._mat(), i, j));
    }
    template<class M> static fmpz_polyxx_ref at(M& m, slong i, slong j)
    {
        return fmpz_polyxx_ref::make(fmpz_poly_mat_entry(m._mat(), i, j));
    }
};

namespace detail {
template<>
struct fmpz_poly_matxx_traits<fmpz_poly_matxx_srcref>
    : matrices::generic_traits_srcref<fmpz_polyxx_srcref> { };
template<>
struct fmpz_poly_matxx_traits<fmpz_poly_matxx_ref>
    : matrices::generic_traits_ref<fmpz_polyxx_ref> { };
template<> struct fmpz_poly_matxx_traits<fmpz_poly_matxx>
    : matrices::generic_traits_nonref<fmpz_polyxx_ref, fmpz_polyxx_srcref> { };

struct fmpz_poly_mat_data
{
    typedef fmpz_poly_mat_t& data_ref_t;
    typedef const fmpz_poly_mat_t& data_srcref_t;

    fmpz_poly_mat_t inner;

    fmpz_poly_mat_data(slong m, slong n)
    {
        fmpz_poly_mat_init(inner, m, n);
    }

    fmpz_poly_mat_data(const fmpz_poly_mat_data& o)
    {
        fmpz_poly_mat_init_set(inner, o.inner);
    }

    fmpz_poly_mat_data(fmpz_poly_matxx_srcref o)
    {
        fmpz_poly_mat_init_set(inner, o._data().inner);
    }

    ~fmpz_poly_mat_data() {fmpz_poly_mat_clear(inner);}
};
} // detail

#define FMPZ_POLY_MATXX_COND_S FLINTXX_COND_S(fmpz_poly_matxx)
#define FMPZ_POLY_MATXX_COND_T FLINTXX_COND_T(fmpz_poly_matxx)

namespace matrices {
template<>
struct outsize<operations::mul_KS_op>
    : outsize<operations::times> { };
template<>
struct outsize<operations::mullow_op>
    : outsize<operations::times> { };
} // matrices

FLINTXX_DEFINE_TEMPORARY_RULES(fmpz_poly_matxx)

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_POLY_MATXX_COND_T, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_set(to._mat(), from._mat()))

FLINTXX_DEFINE_SWAP(fmpz_poly_matxx, fmpz_poly_mat_swap(e1._mat(), e2._mat()))

FLINTXX_DEFINE_EQUALS(fmpz_poly_matxx, fmpz_poly_mat_equal(e1._mat(), e2._mat()))

FLINT_DEFINE_PRINT_PRETTY_COND_2(FMPZ_POLY_MATXX_COND_S, const char*,
        (fmpz_poly_mat_print(from._mat(), extra), 1))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_mat_scalar_mul_fmpz(to._mat(), e1._mat(), e2._fmpz()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_mat_scalar_mul_fmpz_poly(to._mat(), e1._mat(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_add(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_sub(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpz_poly_matxx, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_neg(to._mat(), from._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(transpose_op, fmpz_poly_matxx, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_transpose(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(trace_op, fmpz_polyxx, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_trace(to._poly(), from._mat()))

FLINT_DEFINE_THREEARY_EXPR_COND3(mat_at_op, fmpz_polyxx,
        FMPZ_POLY_MATXX_COND_S, traits::fits_into_slong, traits::fits_into_slong,
        fmpz_poly_set(to._poly(), fmpz_poly_mat_entry(e1._mat(), e2, e3)))

FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_op, fmpz_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_mat_evaluate_fmpz(to._mat(), e1._mat(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(mul_classical_op, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_mul_classical(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(mul_KS_op, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_mul_KS(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_THREEARY_EXPR_COND3(mullow_op, fmpz_poly_matxx,
    FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S, traits::fits_into_slong,
    fmpz_poly_mat_mullow(to._mat(), e1._mat(), e2._mat(), e3))

FLINT_DEFINE_UNARY_EXPR_COND(sqr_op, fmpz_poly_matxx, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_sqr(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_classical_op, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_sqr_classical(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_KS_op, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_sqr_KS(to._mat(), from._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(sqrlow_op, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, traits::fits_into_slong,
        fmpz_poly_mat_sqrlow(to._mat(), e1._mat(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpz_poly_matxx,
        FMPZ_POLY_MATXX_COND_S, traits::is_unsigned_integer,
        fmpz_poly_mat_pow(to._mat(), e1._mat(), e2))

FLINT_DEFINE_THREEARY_EXPR_COND3(pow_trunc_op, fmpz_poly_matxx,
    FMPZ_POLY_MATXX_COND_S, traits::is_unsigned_integer, traits::fits_into_slong,
    fmpz_poly_mat_pow_trunc(to._mat(), e1._mat(), e2, e3))

FLINT_DEFINE_UNARY_EXPR_COND(det_op, fmpz_polyxx, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_det(to._poly(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(det_fflu_op, fmpz_polyxx, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_det_fflu(to._poly(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(det_interpolate_op, fmpz_polyxx,
        FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_det_interpolate(to._poly(), from._mat()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<bool, fmpz_poly_matxx, fmpz_polyxx>::type >::type
    fmpz_poly_mat_inv_rt;

typedef make_ltuple<mp::make_tuple<slong, fmpz_poly_matxx>::type >::type
    fmpz_poly_mat_nullspace_rt;
}
FLINT_DEFINE_UNARY_EXPR_COND(inv_op, rdetail::fmpz_poly_mat_inv_rt,
        FMPZ_POLY_MATXX_COND_S,
        to.template get<0>() = fmpz_poly_mat_inv(to.template get<1>()._mat(),
            to.template get<2>()._poly(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(nullspace_op, rdetail::fmpz_poly_mat_nullspace_rt,
        FMPZ_POLY_MATXX_COND_S, to.template get<0>() = fmpz_poly_mat_nullspace(
            to.template get<1>()._mat(), from._mat()))

FLINT_DEFINE_BINARY_EXPR_COND2(solve_op, rdetail::fmpz_poly_mat_inv_rt,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        to.template get<0>() = fmpz_poly_mat_solve(to.template get<1>()._mat(),
            to.template get<2>()._poly(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(solve_fflu_op, rdetail::fmpz_poly_mat_inv_rt,
        FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        to.template get<0>() = fmpz_poly_mat_solve_fflu(
            to.template get<1>()._mat(),
            to.template get<2>()._poly(), e1._mat(), e2._mat()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<slong, fmpz_poly_matxx, fmpz_polyxx>::type>::type
    fmpz_poly_matxx_fflu_rt;
} // rdetail

FLINT_DEFINE_THREEARY_EXPR_COND3(fflu_op, rdetail::fmpz_poly_matxx_fflu_rt,
        FMPZ_POLY_MATXX_COND_S, traits::is_maybe_perm, tools::is_bool,
        to.template get<0>() = fmpz_poly_mat_fflu(to.template get<1>()._mat(),
            to.template get<2>()._poly(), maybe_perm_data(e2), e1._mat(), e3))

FLINT_DEFINE_UNARY_EXPR_COND(rref_op, rdetail::fmpz_poly_matxx_fflu_rt,
        FMPZ_POLY_MATXX_COND_S,
        to.template get<0>() = fmpz_poly_mat_rref(to.template get<1>()._mat(),
            to.template get<2>()._poly(), from._mat()))

FLINT_DEFINE_THREEARY_EXPR_COND3(solve_fflu_precomp_op, fmpz_poly_matxx,
        traits::is_permxx, FMPZ_POLY_MATXX_COND_S, FMPZ_POLY_MATXX_COND_S,
        fmpz_poly_mat_solve_fflu_precomp(to._mat(), e1._data(),
            e2._mat(), e3._mat()))
} // rules


////////////////////////////////////////////////////////////////////////////
// fmpz_poly_mat_vecxx and prod
////////////////////////////////////////////////////////////////////////////
namespace detail {
struct fmpz_poly_mat_vector_data
{
    slong size;
    fmpz_poly_mat_t* array;

    fmpz_poly_mat_vector_data(slong n, slong rows, slong cols)
        : size(n)
    {
        array = new fmpz_poly_mat_t[n];
        for(slong i = 0; i < n; ++i)
            fmpz_poly_mat_init(array[i], rows, cols);
    }

    ~fmpz_poly_mat_vector_data()
    {
        for(slong i = 0; i < size; ++i)
            fmpz_poly_mat_clear(array[i]);
        delete[] array;
    }

    fmpz_poly_mat_vector_data(const fmpz_poly_mat_vector_data& o)
        : size(o.size)
    {
        array = new fmpz_poly_mat_t[size];
        for(slong i = 0; i < size; ++i)
            fmpz_poly_mat_init_set(array[i], o.array[i]);
    }

    fmpz_poly_matxx_ref at(slong i)
        {return fmpz_poly_matxx_ref::make(array[i]);}
    fmpz_poly_matxx_srcref at(slong i) const
        {return fmpz_poly_matxx_srcref::make(array[i]);}

    bool equals(const fmpz_poly_mat_vector_data& o) const
    {
        if(size != o.size)
            return false;
        for(slong i = 0; i < size; ++i)
            if(!fmpz_poly_mat_equal(array[i], o.array[i]))
                return false;
        return true;
    }
};
} // detail
// TODO temporary allocation

typedef vector_expression<
    detail::wrapped_vector_traits<fmpz_poly_matxx, slong,
        fmpz_poly_matxx_ref, fmpz_poly_matxx_srcref, fmpz_poly_mat_t>,
    operations::immediate,
    detail::fmpz_poly_mat_vector_data> fmpz_poly_mat_vecxx;
// TODO references

template<>
struct enable_vector_rules<fmpz_poly_mat_vecxx> : mp::false_ { };

namespace matrices {
template<>
struct outsize<operations::prod_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head[0].rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head[0].cols();
    }
};
}

namespace rules {
// TODO hack to make code look like references are implemented
template<class T> struct FMPZ_POLY_MAT_VECXX_COND_S
    : mp::equal_types<T, fmpz_poly_mat_vecxx> { };
#define FMPZ_POLY_MAT_VECXX_COND_T FMPZ_POLY_MAT_VECXX_COND_S

// TODO references
FLINT_DEFINE_GET(equals, bool, fmpz_poly_mat_vecxx, e1._data().equals(e2._data()))

FLINT_DEFINE_UNARY_EXPR_COND(prod_op, fmpz_poly_matxx,
        FMPZ_POLY_MAT_VECXX_COND_S,
        fmpz_poly_mat_prod(to._mat(), (fmpz_poly_mat_t * const) from._array(), from.size()))
} // rules

} // flint

#endif
