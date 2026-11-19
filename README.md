# The whisker threshold for WARM graphs with equal rates

Yitzchak Shmalo

In the WARM model of van der Hofstad, Holmes, Kuznetsov and Ruszel, each vertex of a finite graph fires at Poisson
times and reinforces one incident edge, chosen with probability proportional to the alpha-th power of its count.
With equal firing rates, we show that for every alpha > 4/3 the infinitely reinforced edges almost surely form
whiskers (trees of diameter at most three), and for 4/3 < alpha <= 2 stars. The constant 4/3 cannot be lowered. This
improves the thresholds 25 of Holmes and Kleptsyn and 17/4 of Kriebel.

- `main.tex`, `main.pdf`: the paper.
- `check_computations.py`: exact SymPy checks of the Bernstein coefficients, the factorization in Lemma 4.1, the
  derivative in Lemma 3.3, the remark at alpha = 9, the star equilibria and the example at alpha = 3/2.
  Run with `python check_computations.py` (needs SymPy).
