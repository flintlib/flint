---
title: Evaluating Dirichlet L-functions
functions: acb_dirichlet_l acb_dirichlet_l_jet acb_dirichlet_root_number
---

`acb_dirichlet_l` evaluates $L(s, \chi)$ at any complex $s$ for a
Dirichlet character $\chi$ modulo $q$, `acb_dirichlet_l_jet` its Taylor
expansion, and `acb_dirichlet_root_number` the root number of the
functional equation. With `G` and `chi` set to `NULL` (or $q = 1$) these
compute the Riemann zeta function. The methods are compared on
[](algorithms/special-functions/zeta#l-functions).

## Choice of method {#method}

`acb_dirichlet_l` uses, in this order:

1. for an exact integer $s$, `acb_dirichlet_l_si`, which uses closed forms
   with generalized Bernoulli numbers when $s$ has the same parity as $\chi$,
   and otherwise the Euler product, a generalization of Borwein's
   algorithm (`acb_dirichlet_l_ui_borwein_bsplit`) or the approximate
   functional equation;
2. for an exact dyadic $s$ such as $1/2$ (more generally an exact rational
   number of small height), the approximate functional equation
   (`acb_dirichlet_l_fmpq_afe_jet`) when `_acb_dirichlet_l_fmpq_use_afe`
   estimates it to be faster, i.e. at high precision or for large
   conductor;
3. for $q = 1$, `acb_dirichlet_zeta`;
4. for primitive $\chi$ and $\operatorname{Re}(s) < -1/2$ (or $s$ near 0
   and $\chi$ even), the functional equation
   $L(s, \chi) = \epsilon(\chi) (\pi/q)^{s - 1/2} \frac{\Gamma((1-s+\delta)/2)}{\Gamma((s+\delta)/2)} \overline{L(1-\bar s, \chi)}$,
   which moves the evaluation to the right half-plane;
5. for $\operatorname{Re}(s) \ge 8 + p/(2 \log p)$ at precision $p$, the
   Euler product `acb_dirichlet_l_euler_product`;
6. otherwise the decomposition
   $L(s,\chi) = q^{-s}\sum_{k=1}^q \chi(k) \zeta(s, k/q)$
   (`acb_dirichlet_l_hurwitz`), where the $\phi(q)$ Hurwitz zeta values
   come from a precomputed grid of Taylor expansions in $a$
   (`acb_dirichlet_hurwitz_precomp_init_num`) when this is cheaper than
   evaluating each by Euler–Maclaurin summation.

The cost of step 6 grows linearly with $q$, so these functions are meant
for moderate moduli. For many characters with the same modulus,
`acb_dirichlet_l_vec_hurwitz` computes all $\phi(q)$ values with one set
of Hurwitz zeta values and a DFT. For rational $s$ that cannot be represented exactly
by an `acb_t` (such as $1/3$), `acb_dirichlet_l_fmpq` makes the same
choices.

For exact rational $s$ of small height, `acb_dirichlet_l_jet` uses the
approximate functional equation when `_acb_dirichlet_l_jet_use_afe`
estimates it to be faster; otherwise it sums the Hurwitz zeta expansions
$\zeta(s + x, k/q)$ (deflated at $s = 1$, where the poles cancel for
non-principal $\chi$) weighted by the character values, for any length;
for principal characters, `deflate` removes the pole
$\phi(q)/(q(s-1))$.

`acb_dirichlet_root_number` requires a primitive character. For
composite moduli it computes $\epsilon(\chi) = G_q(\chi)/(i^\delta\sqrt q)$
from the Gauss sum, which is a product of local Gauss sums; for prime
power moduli it uses the theta series,
$\epsilon(\chi) = \Theta_q(\chi, 1) / \overline{\Theta_q(\chi, 1)}$,
which needs $O(\sqrt{pq})$ terms at precision $p$.

## Example {#example}

$L(1, \chi_4) = \pi/4$ for the nontrivial character modulo 4, and the
value, derivative and root number at $s = 1/2$ of the character modulo 13
with Conrey number 2 (a generator, so the character has order 12).

```c example name=acb_dirichlet_l
#include <flint/acb.h>
#include <flint/acb_dirichlet.h>

int main(void)
{
    dirichlet_group_t G;
    dirichlet_char_t chi;
    acb_t s, eps;
    acb_ptr L;
    slong prec = 128;

    acb_init(s);
    acb_init(eps);
    L = _acb_vec_init(2);

    dirichlet_group_init(G, 4);
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, 3);
    acb_one(s);
    acb_dirichlet_l(L, s, G, chi, prec);
    flint_printf("L(1, chi_4)  = "); acb_printn(L, 30, 0); flint_printf("\n");
    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);

    dirichlet_group_init(G, 13);
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, 2);
    acb_set_d(s, 0.5);
    acb_dirichlet_l_jet(L, s, G, chi, 0, 2, prec);
    flint_printf("L(1/2, chi)  = "); acb_printn(L, 20, 0); flint_printf("\n");
    flint_printf("L'(1/2, chi) = "); acb_printn(L + 1, 20, 0); flint_printf("\n");
    acb_dirichlet_root_number(eps, G, chi, prec);
    flint_printf("epsilon      = "); acb_printn(eps, 20, 0); flint_printf("\n");
    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);

    acb_clear(s);
    acb_clear(eps);
    _acb_vec_clear(L, 2);
    return 0;
}
```
```output
L(1, chi_4)  = [0.785398163397448309615660845820 +/- 1.25e-31]
L(1/2, chi)  = [1.3468137049240637212 +/- 3.29e-20] + [0.75459682860501534409 +/- 4.17e-21]*I
L'(1/2, chi) = [-0.038182523834367601160 +/- 3.96e-22] + [-0.45987153676368484477 +/- 3.13e-21]*I
epsilon      = [0.52216611241417884070 +/- 2.28e-21] + [0.85284380225587800538 +/- 3.36e-21]*I
```
