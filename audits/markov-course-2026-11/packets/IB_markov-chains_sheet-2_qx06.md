# IB/markov-chains/sheet-2/qx06

## Source locator
```
archive: private/raw/sheets/IB/markov-chains/2012-13/sheet-2.pdf
pages: [2]
original_question: 14
url: https://www.statslab.cam.ac.uk/~rrw1/markov/exmarkov2.pdf
```

## Question (verbatim question.tex)
```latex
We wish to simulate a random variable $X$ for which we know $P(X=i)\propto\alpha_i$, $i\in I$, but we do not know the constant of proportionality because $\sum_{i\in I}\alpha_i$ is too difficult calculate. Let $P=(p_{ij})$ be the (known) transition matrix of some irreducible Markov chain on $I$. Define
\[
q_{ij}=\min\{p_{ij},\alpha_jp_{ji}/\alpha_i\},\ i\ne j,\quad\text{and}\quad q_{ii}=1-\sum_{j\ne i}q_{ij}.
\]
Prove that $Q=(q_{ij})$ is the transition matrix of a reversible Markov chain with invariant distribution proportional to $\alpha$. How can this help with the problem of simulating $X$?
```
