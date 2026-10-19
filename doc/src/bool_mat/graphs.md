---
title: Walks, cycles and components
functions: bool_mat_transitive_closure bool_mat_get_strongly_connected_components bool_mat_nilpotency_degree bool_mat_all_pairs_longest_walk
---

A square Boolean matrix $A$ of order $n$ is the adjacency matrix of a
directed graph on the vertices $0, \ldots, n - 1$, and the Boolean
product of matrices composes walks: $(A^k)_{ij} = 1$ exactly when there
is a walk of length $k$ from $i$ to $j$. The four functions described here
compute properties of all powers of $A$ at once from the structure of the
graph, without forming the powers.

- `bool_mat_transitive_closure` gives $\sum_{k \ge 1} A^k$, the
  reachability relation, by Warshall's algorithm [@War1962] in $O(n^3)$
  operations.
- `bool_mat_get_strongly_connected_components` partitions the vertices
  into strongly connected components (maximal sets of vertices that can
  all reach each other) with Tarjan's depth-first search [@Tar1972], in
  $O(n^2)$ operations since the graph is given by a dense matrix. The
  components are numbered so that edges between different components go
  from higher to lower numbers.
- `bool_mat_nilpotency_degree` finds the least $k$ with $A^k = 0$. Such
  a $k$ exists exactly when the graph has no cycle, and it is then one
  more than the number of edges of a longest path. The function finds a
  topological order by depth-first search and computes longest paths
  along it.
- `bool_mat_all_pairs_longest_walk` refines this to each pair of
  vertices: the largest $k$ with $(A^k)_{ij} = 1$, or $-1$ if there is
  none, or $-2$ if there is none largest. A pair $(i, j)$ has arbitrarily
  long walks when some walk from $i$ to $j$ passes through a component
  containing a cycle. The function works on the condensation of the graph
  (the acyclic graph of components).

The last function serves to bound the error of truncated matrix power
series: if $M$ is a real or complex matrix with the sparsity pattern
$A$, then $(M^k)_{ij} = 0$ whenever $(A^k)_{ij} = 0$, so the entry
$(i, j)$ of $\sum_k a_k M^k$ is a finite sum when $B_{ij} \ge 0$.
`arb_mat_exp` and `acb_mat_exp` use this to compute exponentials of
nilpotent and partly nilpotent matrices with exact zero error bounds where
possible.

## Example {#example}

The graph with edges $0 \to 1$, $1 \to 2$, $2 \to 1$ and $1 \to 3$ has
the components $\{3\}$, $\{1, 2\}$ and $\{0\}$, numbered 0, 1 and 2.
Every walk into the cycle $1 \to 2 \to 1$ can be made arbitrarily long.
Removing the edge $2 \to 1$ leaves an acyclic graph whose longest path,
$0 \to 1 \to 3$ (or $0 \to 1 \to 2$), has two edges, so $A^3 = 0$.

```c example name=bool_mat_graphs
#include <flint/flint.h>
#include <flint/bool_mat.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    bool_mat_t A, C;
    fmpz_mat_t W;
    slong p[4], k, w;

    bool_mat_init(A, 4, 4);
    bool_mat_init(C, 4, 4);
    fmpz_mat_init(W, 4, 4);
    bool_mat_set_entry(A, 0, 1, 1);
    bool_mat_set_entry(A, 1, 2, 1);
    bool_mat_set_entry(A, 2, 1, 1);
    bool_mat_set_entry(A, 1, 3, 1);

    flint_printf("A =\n");
    bool_mat_print(A);

    bool_mat_transitive_closure(C, A);
    flint_printf("transitive closure =\n");
    bool_mat_print(C);

    k = bool_mat_get_strongly_connected_components(p, A);
    flint_printf("%wd components: %{slong*}\n", k, p, (slong) 4);

    flint_printf("nilpotent: %d\n", bool_mat_is_nilpotent(A));
    w = bool_mat_all_pairs_longest_walk(W, A);
    flint_printf("longest walks (return value %wd):\n", w);
    fmpz_mat_print_pretty(W);
    flint_printf("\n");

    /* remove the cycle 1 -> 2 -> 1 */
    bool_mat_set_entry(A, 2, 1, 0);
    flint_printf("without 2 -> 1: nilpotency degree %wd\n",
                 bool_mat_nilpotency_degree(A));
    w = bool_mat_all_pairs_longest_walk(W, A);
    flint_printf("longest walks (return value %wd):\n", w);
    fmpz_mat_print_pretty(W);
    flint_printf("\n");

    bool_mat_clear(A);
    bool_mat_clear(C);
    fmpz_mat_clear(W);
    return 0;
}
```

```output
A =
[0, 1, 0, 0]
[0, 0, 1, 1]
[0, 1, 0, 0]
[0, 0, 0, 0]
transitive closure =
[0, 1, 1, 1]
[0, 1, 1, 1]
[0, 1, 1, 1]
[0, 0, 0, 0]
3 components: [2, 1, 1, 0]
nilpotent: 0
longest walks (return value -2):
[[0 -2 -2 -2]
[-1 -2 -2 -2]
[-1 -2 -2 -2]
[-1 -1 -1 0]
]
without 2 -> 1: nilpotency degree 3
longest walks (return value 2):
[[0 1 2 2]
[-1 0 1 1]
[-1 -1 0 -1]
[-1 -1 -1 0]
]
```

The diagonal entries of the walk matrix are 0 for vertices off any cycle
(the empty walk) and $-2$ for vertices on a cycle.
