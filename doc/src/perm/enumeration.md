---
title: Enumerating permutations
functions: _perm_next_lex _perm_next_heap
---

`_perm_next_lex` and `_perm_next_heap` step through all permutations of
$\{0, \ldots, n - 1\}$ in place, one call per permutation. A loop of the
form

```
do { /* use p */ } while (_perm_next_lex(p, n));
```

started at the identity visits each of the $n!$ permutations exactly
once.

## Lexicographic order {#lex}

`_perm_next_lex` uses the classical successor rule [@Knu1997]: find the
largest $k$ with $p_k < p_{k+1}$, swap $p_k$ with the last entry larger
than it, and reverse the suffix after position $k$. It needs no state
besides the permutation, so the enumeration can start from any
permutation and visits the permutations that follow it in lexicographic
order. A call costs $O(n)$ in the worst case and $O(1)$ on average. The
function returns 0, leaving the array unchanged, at the last permutation
$(n - 1, \ldots, 1, 0)$.

## Heap's order {#heap}

`_perm_next_heap` follows the iterative form of Heap's algorithm
[@Hea1963]. Consecutive permutations differ by a single transposition, so
their parities alternate, and quantities that change little under a
transposition (a determinant expansion, a cost function) can be updated
instead of recomputed. Each call does $O(1)$ work on average. The state
is kept in the array `stack` of length $n$, which must be zero at the
start; the starting permutation may be arbitrary. After the last
permutation the function returns 0, leaving the permutation unchanged and
`stack` zero again. The final permutation is not the reverse of the
first.

## Example {#example}

```c example name=perm_enumerate
#include <flint/flint.h>
#include <flint/perm.h>

int main(void)
{
    slong p[3], stack[3] = {0, 0, 0};
    slong q[6], st[6] = {0, 0, 0, 0, 0, 0};
    slong count = 0;
    int parity = 0, alternates = 1;

    _perm_one(p, 3);
    flint_printf("lexicographic:");
    do
        flint_printf(" %{slong*}", p, (slong) 3);
    while (_perm_next_lex(p, 3));
    flint_printf("\n");

    _perm_one(p, 3);
    flint_printf("Heap:         ");
    do
        flint_printf(" %{slong*}", p, (slong) 3);
    while (_perm_next_heap(p, 3, stack));
    flint_printf("\n");

    /* Heap's order: 6! permutations with alternating parity */
    _perm_one(q, 6);
    do
    {
        count++;
        if (_perm_parity(q, 6) != parity)
            alternates = 0;
        parity ^= 1;
    }
    while (_perm_next_heap(q, 6, st));
    flint_printf("n = 6: %wd permutations, parities alternate: %d\n",
                 count, alternates);

    return 0;
}
```

```output
lexicographic: [0, 1, 2] [0, 2, 1] [1, 0, 2] [1, 2, 0] [2, 0, 1] [2, 1, 0]
Heap:          [0, 1, 2] [1, 0, 2] [2, 0, 1] [0, 2, 1] [1, 2, 0] [2, 1, 0]
n = 6: 720 permutations, parities alternate: 1
```
