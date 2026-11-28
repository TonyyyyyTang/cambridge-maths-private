---
id: IB/methods/sheet-1/qx10
slug: IB_methods_sheet-1_qx10
source_archive: private/raw/sheets/IB/methods/2001-02/sheet-1.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against Question 1 on page 1 of `2001-02/sheet-1.pdf` (pre-rendered 200 DPI PNG `IB_methods_sheet-1_qx10_p1.png`).

- **Result:** `transcription: ok`.
- **Every word, symbol, sign, domain, exponent, and formula matches exactly**, including:
  - Interval $-1 \le x < 1$ (half-open) in part (i), and $0 \le x \le \pi$ (closed) in part (ii).
  - Target formula $f(x) = \frac{8}{15} + \frac{48}{\pi^4}\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n^4}\cos n\pi x$.
  - Exact wording in part (ii): "(a) a Fourier sine series, and (b) a cosine series, each having period $2\pi$", range "$-6\pi$ to $6\pi$", and "$g(x)=x$ and $h(x)=|x|$ each in the range $(-\pi,\pi)$".
- **Purely typographic notes (not discrepancies):** the printed heading *Fourier coefficients (full- and half-range series).* is underlined as well as italicised, whereas `question.tex` renders it with `\emph{...}` and formats `(i)`/`(ii)` using an `enumerate` list.

---

## Full Solution

### Part (i)

#### 1. Computation of the Fourier series

We are given $f(x) = (x^2 - 1)^2 = x^4 - 2x^2 + 1$ on $[-1, 1)$, extended periodically with period $T = 2$ (half-period $L = 1$). The Fourier series on $[-1, 1]$ takes the form
$$f(x) \sim \frac{a_0}{2} + \sum_{n=1}^\infty \bigl(a_n \cos(n\pi x) + b_n \sin(n\pi x)\bigr),$$
where
$$a_n = \int_{-1}^{1} f(x)\cos(n\pi x)\,dx \quad (n \ge 0), \qquad b_n = \int_{-1}^{1} f(x)\sin(n\pi x)\,dx \quad (n \ge 1).$$
Since $f(-x) = ((-x)^2 - 1)^2 = f(x)$ is even on $(-1, 1)$ (the single endpoint $x = -1$ does not affect the integrals), $f(x)\sin(n\pi x)$ is odd and hence
$$b_n = 0 \qquad \text{for all } n \ge 1.$$

For $n = 0$:
$$a_0 = 2\int_0^1 (x^4 - 2x^2 + 1)\,dx = 2\left[\frac{x^5}{5} - \frac{2x^3}{3} + x\right]_0^1 = 2\left(\frac{1}{5} - \frac{2}{3} + 1\right) = \frac{16}{15},$$
so the constant term is $\frac{a_0}{2} = \frac{8}{15}$.

For $n \ge 1$, we compute $a_n = 2\int_0^1 f(x)\cos(n\pi x)\,dx$ by repeated integration by parts. First record the derivatives of $f(x) = (x^2-1)^2$ and their endpoint values:
$$\begin{aligned}
f(x) &= (x^2 - 1)^2, & f(0) &= 1, & f(1) &= 0,\\
f'(x) &= 4x(x^2 - 1) = 4x^3 - 4x, & f'(0) &= 0, & f'(1) &= 0,\\
f''(x) &= 12x^2 - 4, & f''(0) &= -4, & f''(1) &= 8,\\
f'''(x) &= 24x, & f'''(0) &= 0, & f'''(1) &= 24,\\
f^{(4)}(x) &= 24.
\end{aligned}$$

Integrating by parts four times on $[0, 1]$:
$$\begin{aligned}
\int_0^1 f(x)\cos(n\pi x)\,dx
&= \underbrace{\left[f(x)\frac{\sin(n\pi x)}{n\pi}\right]_0^1}_{=\,0} - \frac{1}{n\pi}\int_0^1 f'(x)\sin(n\pi x)\,dx \\
&= \underbrace{\left[f'(x)\frac{\cos(n\pi x)}{(n\pi)^2}\right]_0^1}_{=\,0} - \frac{1}{(n\pi)^2}\int_0^1 f''(x)\cos(n\pi x)\,dx \\
&= \underbrace{\left[-f''(x)\frac{\sin(n\pi x)}{(n\pi)^3}\right]_0^1}_{=\,0} + \frac{1}{(n\pi)^3}\int_0^1 f'''(x)\sin(n\pi x)\,dx \\
&= \left[-f'''(x)\frac{\cos(n\pi x)}{(n\pi)^4}\right]_0^1 + \frac{1}{(n\pi)^4}\int_0^1 24\cos(n\pi x)\,dx \\
&= -\frac{24\cos(n\pi) - 0}{n^4\pi^4} + \frac{24}{n^4\pi^4}\underbrace{\left[\frac{\sin(n\pi x)}{n\pi}\right]_0^1}_{=\,0}
= \frac{24(-1)^{n+1}}{n^4\pi^4}.
\end{aligned}$$
Multiplying by $2$ yields
$$a_n = \frac{48(-1)^{n+1}}{n^4\pi^4} \qquad (n \ge 1).$$
Since the $2$-periodic extension $\tilde{f}$ of $f$ is continuous everywhere and piecewise smooth (indeed $C^2(\mathbb{R})$), by Dirichlet's theorem (or the Weierstrass $M$-test since $|a_n| = O(n^{-4})$) the Fourier series converges uniformly to $\tilde{f}(x)$ for all $x \in \mathbb{R}$, and therefore
$$f(x) = \frac{8}{15} + \frac{48}{\pi^4}\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n^4}\cos(n\pi x)$$
for all $x \in [-1, 1]$ (including $x = \pm 1$).

#### 2. Sketch, differentiability, and convergence rate

- **Sketch of $f(x)$ and its $2$-periodic extension $\tilde{f}(x)$:**
  On $[-1, 1]$, $f(x) = (x^2-1)^2 \ge 0$ is even and bell-shaped: it has a maximum of $1$ at $x = 0$ ($f'(0) = 0$, $f''(0) = -4 < 0$), points of inflection at $x = \pm 1/\sqrt{3}$ where $f(\pm 1/\sqrt{3}) = 4/9$, and meets the $x$-axis tangentially at $x = \pm 1$ ($f(\pm 1) = 0$, $f'(\pm 1) = 0$, $f''(\pm 1) = 8 > 0$). Extending with period $2$ produces a smooth periodic train of identical symmetric arches peaking at $y = 1$ at even integers $x = 2k$ and touching the $x$-axis with zero slope at odd integers $x = 2k+1$ ($k \in \mathbb{Z}$).
- **Differentiability:**
  Inside every open interval $(2k-1, 2k+1)$, $\tilde{f}$ is a polynomial and hence $C^\infty$. At the periodic join points $x = 2k+1$ (matching $x = 1^-$ with $x = -1^+$):
  $$\begin{aligned}
  \tilde{f}(1^-) &= f(1) = 0 = f(-1) = \tilde{f}(1^+),\\
  \tilde{f}'(1^-) &= f'(1) = 0 = f'(-1) = \tilde{f}'(1^+),\\
  \tilde{f}''(1^-) &= f''(1) = 8 = f''(-1) = \tilde{f}''(1^+),\\
  \tilde{f}'''(1^-) &= f'''(1) = +24 \neq -24 = f'''(-1) = \tilde{f}'''(1^+).
  \end{aligned}$$
  Thus $\tilde{f}$ is **twice continuously differentiable on all of $\mathbb{R}$** ($\tilde{f} \in C^2(\mathbb{R})$): $\tilde{f}$, $\tilde{f}'$, and $\tilde{f}''$ are all continuous everywhere, while the third derivative $\tilde{f}'''$ is piecewise smooth with finite jump discontinuities of $\Delta\tilde{f}''' = -48$ at the odd integers $x = 2k+1$.
- **Convergence rate:**
  In general, if a periodic function is $C^{m-1}(\mathbb{R})$ and its $m$-th derivative has finite jump discontinuities (with piecewise smooth $(m+1)$-th derivative), integrating by parts $m+1$ times shows that the Fourier coefficients decay precisely as $O(n^{-(m+1)})$. Here $m = 3$ (the first discontinuity is in $\tilde{f}'''$), so $a_n = O(n^{-4})$. Consequently, the Fourier series converges **very rapidly**, absolutely and uniformly on $\mathbb{R}$, with truncation error bounded by
  $$|f(x) - S_N(x)| \le \frac{48}{\pi^4}\sum_{n=N+1}^\infty \frac{1}{n^4} = O(N^{-3}).$$
  Furthermore, the series may be differentiated term-by-term up to twice to obtain uniformly convergent Fourier series for $\tilde{f}'$ ($O(n^{-3})$) and $\tilde{f}''$ ($O(n^{-2})$), and a third time to obtain a pointwise convergent series for $\tilde{f}'''$ away from the odd integers ($O(n^{-1})$).

---

### Part (ii)

Let $f(x) = x^2$ for $0 \le x \le \pi$.

#### (a) Fourier sine series of period $2\pi$

Extending $f$ as an odd $2\pi$-periodic function gives $f_{\text{odd}}(x) = x|x|$ on $(-\pi, \pi)$ (i.e., $x^2$ on $(0, \pi)$ and $-x^2$ on $(-\pi, 0)$), repeated with period $2\pi$. Its sine series is $S_a(x) = \sum_{n=1}^\infty b_n \sin(nx)$ with
$$b_n = \frac{2}{\pi}\int_0^\pi x^2 \sin(nx)\,dx.$$
Integrating by parts twice:
$$\begin{aligned}
\int_0^\pi x^2\sin(nx)\,dx
&= \left[-x^2\frac{\cos(nx)}{n}\right]_0^\pi + \frac{2}{n}\int_0^\pi x\cos(nx)\,dx \\
&= \frac{\pi^2(-1)^{n+1}}{n} + \frac{2}{n}\left(\left[x\frac{\sin(nx)}{n}\right]_0^\pi - \frac{1}{n}\int_0^\pi \sin(nx)\,dx\right) \\
&= \frac{\pi^2(-1)^{n+1}}{n} + \frac{2}{n^2}\left[\frac{\cos(nx)}{n}\right]_0^\pi
= \frac{\pi^2(-1)^{n+1}}{n} + \frac{2((-1)^n - 1)}{n^3}.
\end{aligned}$$
Multiplying by $\frac{2}{\pi}$ gives
$$b_n = \frac{2\pi(-1)^{n+1}}{n} - \frac{4\bigl(1 - (-1)^n\bigr)}{\pi n^3} = \begin{cases} \dfrac{2\pi}{n} - \dfrac{8}{\pi n^3}, & n \text{ odd}, \\[8pt] -\dfrac{2\pi}{n}, & n \text{ even}. \end{cases}$$
Hence
$$S_a(x) = \sum_{n=1}^\infty \left[\frac{2\pi(-1)^{n+1}}{n} - \frac{4\bigl(1 - (-1)^n\bigr)}{\pi n^3}\right]\sin(nx) = \sum_{n=1}^\infty \frac{2\pi(-1)^{n+1}}{n}\sin(nx) - \frac{8}{\pi}\sum_{m=1}^\infty \frac{\sin\bigl((2m-1)x\bigr)}{(2m-1)^3}.$$
By Dirichlet's theorem, $S_a(x) = x^2$ for $0 \le x < \pi$, whereas at $x = \pi$ (where $f_{\text{odd}}$ jumps from $+\pi^2$ to $-\pi^2$), $S_a(\pi) = 0 = \frac{1}{2}(\pi^2 - \pi^2)$.

#### (b) Fourier cosine series of period $2\pi$

Extending $f$ as an even $2\pi$-periodic function gives $f_{\text{even}}(x) = x^2$ on $[-\pi, \pi]$, repeated with period $2\pi$. Its cosine series is $S_b(x) = \frac{a_0}{2} + \sum_{n=1}^\infty a_n \cos(nx)$ with
$$a_0 = \frac{2}{\pi}\int_0^\pi x^2\,dx = \frac{2\pi^2}{3} \implies \frac{a_0}{2} = \frac{\pi^2}{3},$$
and for $n \ge 1$, integrating by parts twice:
$$\begin{aligned}
a_n = \frac{2}{\pi}\int_0^\pi x^2\cos(nx)\,dx
&= \frac{2}{\pi}\left(\left[x^2\frac{\sin(nx)}{n}\right]_0^\pi - \frac{2}{n}\int_0^\pi x\sin(nx)\,dx\right) \\
&= -\frac{4}{\pi n}\left(\left[-x\frac{\cos(nx)}{n}\right]_0^\pi + \frac{1}{n}\int_0^\pi \cos(nx)\,dx\right) \\
&= \frac{4(-1)^n}{n^2}.
\end{aligned}$$
Hence
$$S_b(x) = \frac{\pi^2}{3} + 4\sum_{n=1}^\infty \frac{(-1)^n}{n^2}\cos(nx).$$
Because $f_{\text{even}}$ is continuous everywhere (matching $\pi^2$ at $x = \pm\pi$) and piecewise smooth, $S_b(x)$ converges uniformly on $\mathbb{R}$ and equals $x^2$ on the entire closed interval $0 \le x \le \pi$.

#### Sketches on $[-6\pi, 6\pi]$

- **Function represented by (a) ($S_a(x)$):**
  On each open interval $((2k-1)\pi, (2k+1)\pi)$ ($k \in \mathbb{Z}$), setting $u = x - 2k\pi \in (-\pi, \pi)$, the graph is the odd S-shaped curve $y = u|u|$ (passing with zero slope and an inflection point through $(2k\pi, 0)$, rising to $+\pi^2$ as $x \to (2k+1)\pi^-$ and falling to $-\pi^2$ as $x \to (2k-1)\pi^+$). At every odd multiple of $\pi$, $x = \pm\pi, \pm 3\pi, \pm 5\pi$, the function has jump discontinuities of $-2\pi^2$ and the series evaluates to the midpoint value $0$ (isolated dots on the $x$-axis). Across $[-6\pi, 6\pi]$, this consists of:
  - the right half-branch $(x+6\pi)^2$ on $[-6\pi, -5\pi)$ rising from $0$ to $\pi^2$;
  - five complete odd S-shaped branches centred at $x = -4\pi, -2\pi, 0, 2\pi, 4\pi$;
  - the left half-branch $-(x-6\pi)^2$ on $(5\pi, 6\pi]$ rising from $-\pi^2$ to $0$;
  - isolated points at $(\pm\pi, 0), (\pm 3\pi, 0), (\pm 5\pi, 0)$.
- **Function represented by (b) ($S_b(x)$):**
  Continuous, even, $2\pi$-periodic chain of upward-opening parabolas $y = (x - 2k\pi)^2$ on $[(2k-1)\pi, (2k+1)\pi]$, touching the $x$-axis smoothly (slope $0$) at $x = 0, \pm 2\pi, \pm 4\pi, \pm 6\pi$ and meeting at sharp peaks (corners) of height $\pi^2$ at $x = \pm\pi, \pm 3\pi, \pm 5\pi$ (where the derivative jumps from $+2\pi$ to $-2\pi$). Across $[-6\pi, 6\pi]$, this comprises a half-parabola on $[-6\pi, -5\pi]$, five complete U-shaped parabolic arcs centred at $-4\pi, -2\pi, 0, 2\pi, 4\pi$, and a half-parabola on $[5\pi, 6\pi]$.

#### Term-by-term differentiation and relation to $g(x) = x$ and $h(x) = |x|$ on $(-\pi, \pi)$

First compute the Fourier series of $g(x) = x$ and $h(x) = |x|$ on $(-\pi, \pi)$:

1. **$g(x) = x$ on $(-\pi, \pi)$** (odd, so $a_n^{(g)} = 0$):
   $$b_n^{(g)} = \frac{2}{\pi}\int_0^\pi x\sin(nx)\,dx = \frac{2}{\pi}\left(\left[-x\frac{\cos(nx)}{n}\right]_0^\pi + \frac{1}{n}\int_0^\pi \cos(nx)\,dx\right) = \frac{2(-1)^{n+1}}{n},$$
   so
   $$g(x) = x = 2\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n}\sin(nx) \qquad \text{for } x \in (-\pi, \pi).$$

2. **$h(x) = |x|$ on $(-\pi, \pi)$** (even, so $b_n^{(h)} = 0$):
   $$a_0^{(h)} = \frac{2}{\pi}\int_0^\pi x\,dx = \pi \implies \frac{a_0^{(h)}}{2} = \frac{\pi}{2},$$
   $$a_n^{(h)} = \frac{2}{\pi}\int_0^\pi x\cos(nx)\,dx = \frac{2}{\pi}\left(\left[x\frac{\sin(nx)}{n}\right]_0^\pi - \frac{1}{n}\int_0^\pi \sin(nx)\,dx\right) = -\frac{2\bigl(1 - (-1)^n\bigr)}{\pi n^2},$$
   so
   $$h(x) = |x| = \frac{\pi}{2} - \frac{2}{\pi}\sum_{n=1}^\infty \frac{1 - (-1)^n}{n^2}\cos(nx) = \frac{\pi}{2} - \frac{4}{\pi}\sum_{m=1}^\infty \frac{\cos\bigl((2m-1)x\bigr)}{(2m-1)^2} \qquad \text{for } x \in [-\pi, \pi].$$

Now differentiate (b) and (a) term-by-term:

- **Differentiating series (b) term-by-term:**
  $$S_b'(x) = -4\sum_{n=1}^\infty \frac{(-1)^n}{n}\sin(nx) = 4\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n}\sin(nx) = 2\left(2\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n}\sin(nx)\right),$$
  which is **precisely the Fourier series for $2g(x) = 2x$** on $(-\pi, \pi)$ (and equals $0$ at $x = \pm\pi$). Term-by-term differentiation succeeds here because the periodic extension $f_{\text{even}}(x)$ is **continuous everywhere on $\mathbb{R}$** and piecewise smooth, so integrating by parts the Fourier coefficients of $f_{\text{even}}'$ introduces no boundary jump terms, giving $b_n(f_{\text{even}}') = -n a_n(f_{\text{even}})$.

- **Differentiating series (a) term-by-term:**
  $$S_a'(x) = \sum_{n=1}^\infty \left[2\pi(-1)^{n+1} - \frac{4\bigl(1 - (-1)^n\bigr)}{\pi n^2}\right]\cos(nx) = 2\pi\sum_{n=1}^\infty (-1)^{n+1}\cos(nx) - \frac{8}{\pi}\sum_{m=1}^\infty \frac{\cos\bigl((2m-1)x\bigr)}{(2m-1)^2}.$$
  On $(-\pi, \pi)$, the classical derivative of $f_{\text{odd}}(x) = x|x|$ is $2|x| = 2h(x)$, whose Fourier series is
  $$2h(x) = \pi - \frac{8}{\pi}\sum_{m=1}^\infty \frac{\cos\bigl((2m-1)x\bigr)}{(2m-1)^2}.$$
  Comparing $S_a'(x)$ with $2h(x)$:
  1. **Classically**, term-by-term differentiation of (a) is **invalid**: the series $S_a'(x)$ **diverges for every $x \in \mathbb{R}$** because its general term contains $2\pi(-1)^{n+1}\cos(nx)$, which does not tend to $0$ as $n \to \infty$. Moreover, $S_a'(x)$ lacks the constant term $\pi$ of $2h(x)$. The failure occurs because $f_{\text{odd}}$ has **jump discontinuities of size $-2\pi^2$** at $x = (2k+1)\pi$, so its Fourier coefficients decay only as $O(n^{-1})$ and differentiation multiplies by $n$.
  2. **Distributionally (generalised functions)**, $S_a'(x)$ is closely and explicitly related to $2h(x)$: rewriting $S_a'(x)$ gives
     $$S_a'(x) = \underbrace{\left[\pi - \frac{8}{\pi}\sum_{m=1}^\infty \frac{\cos\bigl((2m-1)x\bigr)}{(2m-1)^2}\right]}_{=\,2h_{\text{per}}(x)} - 2\pi^2\underbrace{\left[\frac{1}{2\pi} + \frac{1}{\pi}\sum_{n=1}^\infty (-1)^n \cos(nx)\right]}_{=\,\sum_{k=-\infty}^\infty \delta(x - (2k+1)\pi)}.$$
     Here $\frac{1}{2\pi} + \frac{1}{\pi}\sum_{n=1}^\infty (-1)^n \cos(nx)$ is the Fourier series of the $2\pi$-periodic Dirac comb $\sum_{k=-\infty}^\infty \delta\bigl(x - (2k+1)\pi\bigr)$ (since $\frac{1}{\pi}\int_0^{2\pi}\delta(x-\pi)\cos(nx)\,dx = \frac{(-1)^n}{\pi}$ and $\frac{a_0}{2} = \frac{1}{2\pi}$). Thus, in the sense of distributions,
     $$S_a'(x) = 2h_{\text{per}}(x) - 2\pi^2\sum_{k=-\infty}^\infty \delta\bigl(x - (2k+1)\pi\bigr),$$
     which is the exact distributional derivative of $f_{\text{odd}}$: the classical derivative $2|x|$ on $(-\pi, \pi)$ together with the delta-function spikes of weight $-2\pi^2$ produced by the downward jumps from $+\pi^2$ to $-\pi^2$ at every odd multiple of $\pi$.

---

## Final Answers

- **(i)**
  - Fourier series: $f(x) = \dfrac{8}{15} + \dfrac{48}{\pi^4}\displaystyle\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n^4}\cos(n\pi x)$ on $[-1, 1]$ (uniformly convergent everywhere on $\mathbb{R}$ to the $2$-periodic extension).
  - Differentiability & convergence: the $2$-periodic extension $\tilde{f}$ is $C^2(\mathbb{R})$ ($\tilde{f}, \tilde{f}', \tilde{f}''$ continuous everywhere, with $\tilde{f}''(\pm 1) = 8$), while $\tilde{f}'''$ has jump discontinuities of $-48$ at $x = 2k+1$; consequently the Fourier coefficients decay as $O(n^{-4})$ and the series converges rapidly, absolutely, and uniformly with tail $O(N^{-3})$.
- **(ii)**
  - **(a) Sine series:**
    $$S_a(x) = \sum_{n=1}^\infty \left[\frac{2\pi(-1)^{n+1}}{n} - \frac{4\bigl(1-(-1)^n\bigr)}{\pi n^3}\right]\sin(nx)$$
    (converges to $x^2$ on $[0, \pi)$ and to $0$ at $x = \pi$).
  - **(b) Cosine series:**
    $$S_b(x) = \frac{\pi^2}{3} + 4\sum_{n=1}^\infty \frac{(-1)^n}{n^2}\cos(nx)$$
    (converges uniformly to $x^2$ on $[0, \pi]$).
  - **Fourier series of $g(x) = x$ and $h(x) = |x|$ on $(-\pi, \pi)$:**
    $$g(x) = 2\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n}\sin(nx), \qquad h(x) = \frac{\pi}{2} - \frac{4}{\pi}\sum_{m=1}^\infty \frac{\cos\bigl((2m-1)x\bigr)}{(2m-1)^2}.$$
  - **Term-by-term differentiation:**
    - $S_b'(x) = 4\sum_{n=1}^\infty \frac{(-1)^{n+1}}{n}\sin(nx) = 2g(x)$ on $(-\pi, \pi)$ (valid because $f_{\text{even}}$ is continuous and piecewise smooth).
    - $S_a'(x) = \sum_{n=1}^\infty \left[2\pi(-1)^{n+1} - \frac{4(1-(-1)^n)}{\pi n^2}\right]\cos(nx)$ diverges everywhere classically (because $f_{\text{odd}}$ has jumps of $-2\pi^2$ at odd multiples of $\pi$), and distributionally equals $2h_{\text{per}}(x) - 2\pi^2\sum_{k\in\mathbb{Z}}\delta\bigl(x - (2k+1)\pi\bigr)$.

---

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx10_check.py` (exit code 0).

```text
=== Part (i): f(x) = (x^2 - 1)^2 on [-1, 1) ===
Exact a_0/2 = 8/15
Fourier coeffs n=1..12 verified: max |a_n err| = 1.69e-14, max |b_n| = 7.77e-17
Join checks at x=+-1: f=0, f'=0, f''=8 (C^2 continuous), jump in f''' = -48.0
Series convergence (N=2000): max err f=2.05e-11, f'=6.11e-10, f''=2.43e-03
=== Part (ii): f(x) = x^2 on [0, pi], sine/cosine series & differentiation ===
Half-range coeffs n=1..15 verified: max |b_n err|=4.46e-11, max |a_n err|=5.71e-12
Extensions on [-6pi, 6pi] (N=1000): max err S_a=4.32e-03 (0 at odd pi), max err S_b=4.00e-03
Term-by-term S_b'(x) coeffs match 2*g(x) identically for n=1..49.
S_a'(0.7) terms for n=100..119 reach 6.222 (does not tend to 0 -> classically divergent).
Distributional pairing <S_a', exp(cos x)>: |LHS - RHS| = 1.25e-12
Distributional pairing <S_a', 1/(2 + cos x)>: |LHS - RHS| = 3.94e-13
Distributional pairing <S_a', cos(3x) + 0.5*sin(2x)>: |LHS - RHS| = 7.36e-12
ALL CHECKS PASSED for IB/methods/sheet-1/qx10.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** high.
- **Phase B watchpoints:**
  1. In part (i), check that the solution explicitly identifies $\tilde{f} \in C^2(\mathbb{R})$ (both $\tilde{f}'$ and $\tilde{f}''$ are continuous across $x = \pm 1$, since $f''(\pm 1) = 8$) and that the first jump is in $\tilde{f}'''$ ($+24 \to -24$), which is what explains the $O(n^{-4})$ decay rate (not merely "$C^1$").
  2. In part (ii)(a), check the sign of the $n^{-3}$ term in $b_n$: $+\frac{4((-1)^n - 1)}{\pi n^3} = -\frac{4(1-(-1)^n)}{\pi n^3}$ (so for odd $n$, $b_n = \frac{2\pi}{n} - \frac{8}{\pi n^3}$, NOT $+$), and verify the midpoint value $S_a(\pm\pi) = 0$ in the sketch/convergence statement.
  3. In part (ii) differentiation, check that the solution explains BOTH why $S_b'$ equals $2g(x)$ (continuity of $f_{\text{even}}$) AND why $S_a'$ fails classically (divergence + missing constant term $\pi$ from jumps at odd multiples of $\pi$), ideally identifying the exact distributional relation $S_a'(x) = 2h_{\text{per}}(x) - 2\pi^2\sum_{k\in\mathbb{Z}}\delta(x-(2k+1)\pi)$.
