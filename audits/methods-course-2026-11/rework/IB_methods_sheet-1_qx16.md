---
id: IB/methods/sheet-1/qx16
slug: IB_methods_sheet-1_qx16
source_archive: private/raw/sheets/IB/methods/2001-02/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against Question 10 on page 2 of `2001-02/sheet-1.pdf` (pre-rendered 200 DPI PNG `IB_methods_sheet-1_qx16_p2.png`). Also confirmed that page 2 is the final page of the sheet (ending with Question 12 and the sheet footnote; no page 3 exists).

- **Result:** `transcription: ok`.
- **Every word, symbol, sign, bracket, limit, and formula matches exactly**, including:
  - $\int_a^b f(x)\delta[\phi(x)]dx = \frac{f(c)}{|\phi'(c)|}$.
  - $\int_{-\infty}^{+\infty}|x|\delta(x^2-a^2)dx = 1$.
- **Purely typographic note (not a discrepancy):** the printed heading *Delta function properties.* is underlined as well as italicised, whereas `question.tex` uses `\emph{Delta function properties.}`.

---

## Full Solution

Throughout, we assume $\phi$ is continuously differentiable ($C^1$) on $[a, b]$ with non-vanishing derivative at its zeros (simple zeros), and $f$ is a continuous test function at the zeros of $\phi$. (Each change-of-variables manipulation below is made completely rigorous by replacing $\delta$ with any nascent delta sequence $\delta_\epsilon(u) = \epsilon^{-1}\eta(u/\epsilon)$ with $\int_{-\infty}^\infty \eta(s)\,ds = 1$ and taking $\epsilon \to 0^+$.)

### 1. $\phi(x)$ monotone increasing on $[a, b]$ with a simple zero at $x = c \in (a, b)$

Since $\phi$ is strictly monotone increasing and $C^1$ on $[a, b]$ with $\phi(c) = 0$ and $\phi'(c) > 0$ (so $x = c$ is a simple zero), $\phi$ has a $C^1$ inverse $x = \phi^{-1}(y)$ on a neighbourhood of $y = 0$ (indeed on $[\phi(a), \phi(b)]$ if $\phi' > 0$ on $[a, b]$; away from $c$, $\phi(x) \ne 0$ so $\delta[\phi(x)] = 0$). Moreover,
$$\phi(a) < \phi(c) = 0 < \phi(b), \qquad \phi'(c) > 0 \implies \phi'(c) = |\phi'(c)|.$$

Make the change of variables $y = \phi(x)$, so $dy = \phi'(x)\,dx$ and $dx = \dfrac{dy}{\phi'\bigl(\phi^{-1}(y)\bigr)}$:
$$\int_a^b f(x)\delta[\phi(x)]\,dx = \int_{\phi(a)}^{\phi(b)} \frac{f\bigl(\phi^{-1}(y)\bigr)}{\phi'\bigl(\phi^{-1}(y)\bigr)}\,\delta(y)\,dy.$$
Since $\phi(a) < 0 < \phi(b)$ and the integrand factor $\dfrac{f(\phi^{-1}(y))}{\phi'(\phi^{-1}(y))}$ is continuous at $y = 0$ with $\phi^{-1}(0) = c$, the sifting property of the Dirac delta function ($\int_\alpha^\beta g(y)\delta(y)\,dy = g(0)$ whenever $\alpha < 0 < \beta$) gives
$$\int_a^b f(x)\delta[\phi(x)]\,dx = \frac{f\bigl(\phi^{-1}(0)\bigr)}{\phi'\bigl(\phi^{-1}(0)\bigr)} = \frac{f(c)}{\phi'(c)} = \frac{f(c)}{|\phi'(c)|}.$$

### 2. $\phi(x)$ monotone decreasing on $[a, b]$ with a simple zero at $x = c \in (a, b)$

Now $\phi$ is strictly monotone decreasing on $[a, b]$ with $\phi(c) = 0$ and $\phi'(c) < 0$. Hence
$$\phi(b) < \phi(c) = 0 < \phi(a), \qquad \phi'(c) = -|\phi'(c)|.$$
Making the same substitution $y = \phi(x)$, $dy = \phi'(x)\,dx$:
$$\int_a^b f(x)\delta[\phi(x)]\,dx = \int_{\phi(a)}^{\phi(b)} \frac{f\bigl(\phi^{-1}(y)\bigr)}{\phi'\bigl(\phi^{-1}(y)\bigr)}\,\delta(y)\,dy.$$
Because $\phi(b) < \phi(a)$, we reverse the limits of integration to place the smaller limit at the bottom, which introduces a minus sign:
$$\int_a^b f(x)\delta[\phi(x)]\,dx = -\int_{\phi(b)}^{\phi(a)} \frac{f\bigl(\phi^{-1}(y)\bigr)}{\phi'\bigl(\phi^{-1}(y)\bigr)}\,\delta(y)\,dy.$$
Now $\phi(b) < 0 < \phi(a)$, so the sifting property at $y = 0$ gives
$$\int_a^b f(x)\delta[\phi(x)]\,dx = -\frac{f(c)}{\phi'(c)} = \frac{f(c)}{-\phi'(c)} = \frac{f(c)}{|\phi'(c)|},$$
since $\phi'(c) < 0 \implies -\phi'(c) = |\phi'(c)|$. Thus the exact same formula holds whether $\phi$ is monotone increasing or monotone decreasing.

### 3. General $\phi(x)$ with simple zeros

Suppose $\phi \in C^1([a, b])$ has simple zeros $\{x_i\}$ lying in the open interval $(a, b)$ — meaning $\phi(x_i) = 0$ and $\phi'(x_i) \ne 0$ for each $i$ (and no zeros at the endpoints $a, b$; for integrals over $(-\infty, \infty)$, isolated simple zeros).

Because $\phi'$ is continuous and $\phi'(x_i) \ne 0$, around each zero $x_i$ there is a small closed interval $[a_i, b_i] \subset (a, b)$ containing $x_i$ in its interior such that:
1. the intervals $[a_i, b_i]$ are pairwise disjoint;
2. on each $[a_i, b_i]$, $\phi'$ keeps a constant non-zero sign, so $\phi$ is strictly monotone (increasing if $\phi'(x_i) > 0$, decreasing if $\phi'(x_i) < 0$) and has $x_i$ as its unique zero in $[a_i, b_i]$;
3. on the complementary set $[a, b] \setminus \bigcup_i [a_i, b_i]$, $\phi(x) \ne 0$, so $\delta[\phi(x)] = 0$ there and that region contributes $0$ to the integral.

Splitting the integral across these disjoint subintervals and applying the monotone formula to each $[a_i, b_i]$:
$$\int_a^b f(x)\delta[\phi(x)]\,dx = \sum_{i} \int_{a_i}^{b_i} f(x)\delta[\phi(x)]\,dx = \sum_{i:\,\phi(x_i)=0,\,x_i\in(a,b)} \frac{f(x_i)}{|\phi'(x_i)|}.$$
Equivalently, as a distributional identity (holding when integrated against any continuous test function $f$):
$$\delta[\phi(x)] = \sum_{i:\,\phi(x_i)=0} \frac{\delta(x - x_i)}{|\phi'(x_i)|}.$$

### 4. Evaluation of $\int_{-\infty}^{+\infty} |x|\delta(x^2 - a^2)\,dx$

Assume $a \in \mathbb{R} \setminus \{0\}$ (so that the zeros of $\phi(x) = x^2 - a^2$ are simple). Then:
- The zeros of $\phi(x) = x^2 - a^2$ are $x_1 = a$ and $x_2 = -a$ (two distinct real zeros).
- The derivative is $\phi'(x) = 2x$, so
  $$\phi'(a) = 2a \ne 0, \qquad \phi'(-a) = -2a \ne 0 \implies |\phi'(\pm a)| = 2|a| > 0.$$
  Both zeros are simple, so by the general formula:
  $$\delta(x^2 - a^2) = \frac{\delta(x - a) + \delta(x + a)}{2|a|}.$$
- Since $f(x) = |x|$ is continuous (in fact smooth) at $x = \pm a$ for $a \ne 0$, applying the general formula yields
  $$\int_{-\infty}^{+\infty} |x|\delta(x^2 - a^2)\,dx = \frac{|a|}{|\phi'(a)|} + \frac{|-a|}{|\phi'(-a)|} = \frac{|a|}{2|a|} + \frac{|a|}{2|a|} = \frac{1}{2} + \frac{1}{2} = 1.$$

*(Alternative direct check for $a \ne 0$:* by evenness of the integrand, $\int_{-\infty}^{+\infty} |x|\delta(x^2 - a^2)\,dx = 2\int_0^\infty x\delta(x^2 - a^2)\,dx$; substituting $u = x^2 - a^2$ with $du = 2x\,dx$ gives $\int_{-a^2}^{+\infty} \delta(u)\,du = 1$ since $-a^2 < 0 < +\infty$.)

*Remark on the degenerate case $a = 0$:* When $a = 0$, $\phi(x) = x^2$ has a **double** (non-simple) zero at $x = 0$ ($\phi'(0) = 0$), so the simple-zero formula does not apply. Indeed, with a nascent delta sequence $\delta_\epsilon(u) = \epsilon^{-1}\eta(u/\epsilon)$, substituting $u = x^2$ gives
$$\int_{-\infty}^{+\infty} |x|\delta_\epsilon(x^2)\,dx = 2\int_0^\infty x\delta_\epsilon(x^2)\,dx = \int_0^\infty \delta_\epsilon(u)\,du = \int_0^\infty \eta(s)\,ds,$$
which depends on the choice of regularisation $\eta$ (e.g., $\frac{1}{2}$ for any symmetric nascent delta, and other values for asymmetric kernels) rather than equalling $1$. Thus $a \ne 0$ is an essential implicit hypothesis.

---

## Final Answers

- **Monotone $\phi$ (increasing or decreasing) with simple zero $c \in (a, b)$:**
  $$\int_a^b f(x)\delta[\phi(x)]\,dx = \frac{f(c)}{|\phi'(c)|}.$$
- **General $\phi$ with simple zeros $\{x_i\}$:**
  $$\int_a^b f(x)\delta[\phi(x)]\,dx = \sum_{\substack{i:\,\phi(x_i)=0 \\ x_i \in (a,b)}} \frac{f(x_i)}{|\phi'(x_i)|}, \qquad \text{i.e.,} \qquad \delta[\phi(x)] = \sum_{i:\,\phi(x_i)=0} \frac{\delta(x - x_i)}{|\phi'(x_i)|}.$$
- **Application ($a \ne 0$):**
  $$\delta(x^2 - a^2) = \frac{\delta(x - a) + \delta(x + a)}{2|a|} \implies \int_{-\infty}^{+\infty} |x|\delta(x^2 - a^2)\,dx = \frac{|a|}{2|a|} + \frac{|-a|}{2|a|} = 1.$$

---

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx16_check.py` (exit code 0).

```text
Both symmetric and asymmetric nascent delta kernels integrate to 1.
Monotone increasing (Gaussian, eps=1e-05): num=0.6248459577, exact=0.6248459576, err=9.04e-12
Monotone increasing (Asymmetric, eps=1e-05): num=0.6248451727, exact=0.6248459576, err=7.85e-07
Monotone decreasing (Gaussian, eps=1e-05): num=0.5495737569, exact=0.5495737569, err=2.32e-11
Monotone decreasing (Asymmetric, eps=1e-05): num=0.5495754743, exact=0.5495737569, err=1.72e-06
General multi-zero (Gaussian, eps=1e-05): num=2.0000000006, exact=2.0000000000, err=5.98e-10
General multi-zero (Asymmetric, eps=1e-05): num=1.9999875157, exact=2.0000000000, err=1.25e-05
Target integral a=+1.70 (Gaussian, eps=0.0001): val=1.0000000000, err=6.24e-13
Target integral a=+1.70 (Asymmetric, eps=0.0001): val=1.0000000000, err=3.53e-13
Target integral a=-0.80 (Gaussian, eps=0.0001): val=1.0000000000, err=2.30e-14
Target integral a=-0.80 (Asymmetric, eps=0.0001): val=1.0000000000, err=5.24e-14
Target integral a=+0.35 (Gaussian, eps=0.0001): val=1.0000000000, err=0.00e+00
Target integral a=+0.35 (Asymmetric, eps=0.0001): val=1.0000000000, err=6.00e-15
Degenerate a=0 check: Gaussian -> 0.500000 (=0.5), Asymmetric -> 0.750950 (=0.750950) => ill-defined at a=0!
ALL CHECKS PASSED for IB/methods/sheet-1/qx16.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** high.
- **Phase B watchpoints:**
  1. In the monotone-decreasing proof, verify that BOTH sign flips are shown explicitly: reversing the limits $\int_{\phi(a)}^{\phi(b)} = -\int_{\phi(b)}^{\phi(a)}$ (since $\phi(b) < 0 < \phi(a)$) AND $-\phi'(c) = |\phi'(c)|$ (since $\phi'(c) < 0$).
  2. In the general formula, verify the hypothesis that zeros $x_i$ are **simple** ($\phi'(x_i) \ne 0$) and strictly inside $(a, b)$ (not at endpoints), and that the domain is partitioned into small disjoint intervals on each of which $\phi$ is strictly monotone.
  3. In the final integral $\int_{-\infty}^{+\infty}|x|\delta(x^2-a^2)\,dx = 1$, watch for:
     - using $|\phi'(\pm a)| = 2|a|$ (NOT $2a$, which fails for $a < 0$);
     - noting that $a \ne 0$ is needed so the zeros are simple (for $a = 0$, $x^2$ has a double zero and the integral is regularisation-dependent, giving $1/2$ for symmetric nascent deltas).
