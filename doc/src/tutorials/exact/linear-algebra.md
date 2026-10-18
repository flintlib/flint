---
title: Exact linear algebra
summary: Eigenvalues of the Hilbert matrix as exact algebraic numbers with verification of the trace and determinant, and an exact discrete Fourier transform with roots of unity; the example programs hilbert_matrix_ca and dft.
order: 3
---

Matrices and vectors of exact numbers combine the field arithmetic of
[](ca) with the generic linear algebra of [](gr_mat) (see
[](algorithms/exact/linear-algebra)). This tutorial computes the
eigenvalues of a Hilbert matrix and a discrete Fourier transform
exactly, following the example programs `hilbert_matrix_ca` and `dft`
in the `examples` directory.

## Eigenvalues of the Hilbert matrix {#hilbert}

The Hilbert matrix $H_n = (1/(i + j - 1))_{i,j=1}^n$ is a classical
example of an ill-conditioned matrix: its smallest eigenvalue decreases
exponentially with $n$. Its eigenvalues are the roots of the
characteristic polynomial, which has rational coefficients and, for
small $n$, is irreducible. `ca_mat_eigenvalues` computes the
characteristic polynomial, factors it and returns the roots with their
multiplicities as `ca_t` numbers; each eigenvalue becomes the generator
of a number field of degree $n$. The program checks the identities

$$
\lambda_1 + \cdots + \lambda_n = \operatorname{tr}(H_n), \qquad
\lambda_1 \cdots \lambda_n = \det(H_n).
$$

The sum and the product of the conjugate eigenvalues mix $n$ generators
whose only relations are that they are the roots of one polynomial.
[](ca) recognises such sums and products of all roots of a polynomial
with Vieta's formulas, for degrees up to the option
`CA_OPT_VIETA_LIMIT` (6 by default), and otherwise falls back on the
zero test, which evaluates the numbers as algebraic numbers in
[](qqbar).

## Discrete Fourier transform {#dft}

The second part computes the discrete Fourier transform

$$
X_k = \sum_{n=0}^{N-1} \omega^{-kn} x_n, \qquad \omega = e^{2 \pi i / N},
$$

of the vector $x_n = \sqrt{n + 2}$ for $N = 5$ by naive summation, and
then checks that the inverse transform recovers $x$ exactly. The root of
unity $\omega = e^{2\pi i/5}$ is recognised as an algebraic number and
represented by its minimal polynomial $a^4 - a^3 + a^2 - a + 1$, so all
computations take place in a number field with five algebraic
generators and the differences reduce to 0 automatically.

```c example name=tutorial_exact_linear_algebra
#include <flint/calcium.h>
#include <flint/ca.h>
#include <flint/ca_vec.h>
#include <flint/ca_mat.h>

int main(void)
{
    ca_ctx_t ctx;
    ca_mat_t H;
    ca_vec_t eig;
    ulong mult[5];
    ca_t s, p, t;
    ca_ptr x, X, w;
    slong i, k, n, N = 5;

    ca_ctx_init(ctx);
    ca_init(s, ctx);
    ca_init(p, ctx);
    ca_init(t, ctx);

    /* eigenvalues of the 5 x 5 Hilbert matrix */
    ca_mat_init(H, 5, 5, ctx);
    ca_vec_init(eig, 0, ctx);
    ca_mat_hilbert(H, ctx);
    if (ca_mat_eigenvalues(eig, mult, H, ctx))
    {
        ca_zero(s, ctx);
        ca_one(p, ctx);
        for (i = 0; i < ca_vec_length(eig, ctx); i++)
        {
            flint_printf("lambda_%wd = ", i);
            ca_printn(ca_vec_entry_ptr(eig, i), 10, ctx);
            flint_printf("   (multiplicity %wu)\n", mult[i]);
            ca_add(s, s, ca_vec_entry_ptr(eig, i), ctx);
            ca_mul(p, p, ca_vec_entry_ptr(eig, i), ctx);
        }
        ca_mat_trace(t, H, ctx);
        flint_printf("sum = trace: "); truth_print(ca_check_equal(s, t, ctx));
        flint_printf("\n");
        ca_mat_det(t, H, ctx);
        flint_printf("product = det = "); ca_print(t, ctx);
        flint_printf(": "); truth_print(ca_check_equal(p, t, ctx));
        flint_printf("\n");
    }

    /* x - IDFT(DFT(x)) for x_k = sqrt(k + 2), length N = 5 */
    x = _ca_vec_init(N, ctx);
    X = _ca_vec_init(N, ctx);
    w = _ca_vec_init(N, ctx);
    for (k = 0; k < N; k++)
    {
        ca_sqrt_ui(x + k, k + 2, ctx);
        /* w_k = exp(2 pi i k / N) */
        ca_pi_i(w + k, ctx);
        ca_mul_ui(w + k, w + k, 2 * k, ctx);
        ca_div_ui(w + k, w + k, N, ctx);
        ca_exp(w + k, w + k, ctx);
    }
    for (k = 0; k < N; k++)
    {
        ca_zero(X + k, ctx);
        for (n = 0; n < N; n++)
        {
            ca_mul(t, x + n, w + (N - (k * n) % N) % N, ctx);
            ca_add(X + k, X + k, t, ctx);
        }
    }
    flint_printf("DFT(x)_1 = "); ca_print(X + 1, ctx); flint_printf("\n");
    for (k = 0; k < N; k++)
    {
        ca_zero(s, ctx);
        for (n = 0; n < N; n++)
        {
            ca_mul(t, X + n, w + (k * n) % N, ctx);
            ca_add(s, s, t, ctx);
        }
        ca_div_ui(s, s, N, ctx);
        ca_sub(s, x + k, s, ctx);
        flint_printf("x_%wd - IDFT(DFT(x))_%wd = ", k, k);
        ca_print(s, ctx);
        flint_printf("\n");
    }

    _ca_vec_clear(x, N, ctx);
    _ca_vec_clear(X, N, ctx);
    _ca_vec_clear(w, N, ctx);
    ca_mat_clear(H, ctx);
    ca_vec_clear(eig, ctx);
    ca_clear(s, ctx);
    ca_clear(p, ctx);
    ca_clear(t, ctx);
    ca_ctx_clear(ctx);
    return 0;
}
```

```output
lambda_0 = 1.567050691   (multiplicity 1)
lambda_1 = 0.2085342186   (multiplicity 1)
lambda_2 = 0.01140749162   (multiplicity 1)
lambda_3 = 0.0003058980402   (multiplicity 1)
lambda_4 = 3.287928772e-6   (multiplicity 1)
sum = trace: T_TRUE
product = det = 3.74930e-12 {1/266716800000}: T_TRUE
DFT(x)_1 = -0.720670 + 0.821082*I {a^3*c-a^3*d-a^2*c+a^2*d*e+a*c-2*a-c+e where a = 0.809017 + 0.587785*I [a^4-a^3+a^2-a+1=0], b = 2.44949 [b^2-6=0], c = 2.23607 [c^2-5=0], d = 1.73205 [d^2-3=0], e = 1.41421 [e^2-2=0]}
x_0 - IDFT(DFT(x))_0 = 0
x_1 - IDFT(DFT(x))_1 = 0
x_2 - IDFT(DFT(x))_2 = 0
x_3 - IDFT(DFT(x))_3 = 0
x_4 - IDFT(DFT(x))_4 = 0
```

## The example programs {#programs}

`hilbert_matrix_ca n` performs the computation of the first part for
$H_n$ and prints the trace and the determinant. It raises
`CA_OPT_QQBAR_DEG_LIMIT` to 10000, since the verification may need
algebraic numbers of high degree. With `-novieta` it disables Vieta's
formulas (`CA_OPT_VIETA_LIMIT` = 0): the sum of the eigenvalues is
then printed as the formal sum $a + b + c + d + e$ of five generators
with the same minimal polynomial, and `ca_check_equal` still proves that
it equals the trace. With `-qqbar` it computes the eigenvalues with
`qqbar_eigenvalues_fmpq_mat` and adds and multiplies them as
[](qqbar) numbers; the intermediate sums have degrees up to
$\binom{n}{\lfloor n/2 \rfloor}$, so this is much slower for large $n$.

`dft N` computes the transform for several input sequences selected by
`-input i`: integers, square roots, logarithms, roots of unity and
reciprocals such as $1/(1 + \pi \sqrt{n + 2})$. For logarithms and $\pi$, the fields
contain transcendental generators and the zero test does the final
verification. The option `-timing` selects whether fields are cached
between repeated runs, which shows the cost of constructing fields and
finding relations compared with the arithmetic itself.
