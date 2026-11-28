---
id: IB/methods/sheet-3/qx24
slug: IB_methods_sheet-3_qx24
source_archive: private/raw/sheets/IB/methods/2016-17/sheet-3.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-3_qx24_p2.png` (Michaelmas 2016, 1B Methods, David Skinner, Sheet 3, Question 6). Question 6 sits in the middle of page 2 and finishes well above Question 7.

- **Status:** `ok` (exact match for every word, symbol, paragraph break, `\widetilde\delta_n(k)`, integral limits, and punctuation).
- **Convention note from printed sheet:** Question 7 immediately below on page 2 states that $f(x) = e^{-x}$ ($x > 0$), $0$ ($x < 0$) has Fourier transform $\tilde{f}(k) = (1 - ik)/(1 + k^2)$, which confirms that this course uses the standard IB Methods convention
  \[
  \tilde{f}(k) = \mathcal{F}[f](k) = \int_{-\infty}^{\infty} f(x)\,e^{-ikx}\,dx, \qquad
  f(x) = \frac{1}{2\pi}\int_{-\infty}^{\infty} \tilde{f}(k)\,e^{ikx}\,dk.
  \]

## Full Solution

### 1. Fourier transform convention, differentiation properties, and shift property

We adopt the IB Methods convention (confirmed by Question 7 on the same sheet)
\[
\tilde{h}(k) = \mathcal{F}[h](k) = \int_{-\infty}^{\infty} h(x)\,e^{-ikx}\,dx, \qquad
h(x) = \frac{1}{2\pi}\int_{-\infty}^{\infty} \tilde{h}(k)\,e^{ikx}\,dk.
\]
We assume $n > 0$ and $\mu \in \mathbb{R}$ throughout (as implied by $\delta_n(x) = (n/\sqrt{\pi})f(x)$ and taking $n \to \infty$ below; for general real $n \ne 0$, simply replace $n$ by $|n|$ in the prefactor of $\tilde{f}(k)$, while $n = 0$ gives the constant function $f(x) \equiv 1$, which is not in $L^1(\mathbb{R})$ and only has the distributional transform $2\pi\delta(k)$).

First we record the shift and differentiation properties for a smooth rapidly decaying function $g(x)$ (such as $g(x) = e^{-n^2 x^2}$, which decays faster than any power of $|x|$ as $|x| \to \infty$):

1. **Shift property:** If $f(x) = g(x - \mu)$, substituting $u = x - \mu$ gives
   \[
   \tilde{f}(k) = \int_{-\infty}^{\infty} g(x - \mu)\,e^{-ikx}\,dx
   = \int_{-\infty}^{\infty} g(u)\,e^{-ik(u + \mu)}\,du
   = e^{-ik\mu}\,\tilde{g}(k).
   \]
2. **Spatial differentiation:** Integrating by parts and noting $g(x) \to 0$ as $|x| \to \infty$,
   \[
   \mathcal{F}[g'](k) = \int_{-\infty}^{\infty} g'(x)\,e^{-ikx}\,dx
   = \Bigl[g(x)e^{-ikx}\Bigr]_{-\infty}^{\infty} + ik\int_{-\infty}^{\infty} g(x)\,e^{-ikx}\,dx
   = ik\,\tilde{g}(k).
   \]
3. **Frequency differentiation:** Since $|x\,g(x)| \in L^1(\mathbb{R})$, differentiating under the integral sign is justified by dominated convergence:
   \[
   \frac{d\tilde{g}}{dk} = \int_{-\infty}^{\infty} (-ix)\,g(x)\,e^{-ikx}\,dx = -i\,\mathcal{F}[x\,g(x)](k)
   \implies
   \mathcal{F}[x\,g(x)](k) = i\,\frac{d\tilde{g}}{dk}.
   \]

### 2. Fourier transform of $f(x) = \exp[-n^2(x - \mu)^2]$

First consider the unshifted Gaussian $g(x) = e^{-n^2 x^2}$ ($n > 0$). It satisfies the first-order linear ODE
\[
g'(x) = -2n^2 x\,g(x).
\]
Taking the Fourier transform of both sides and using the two differentiation properties above gives
\[
ik\,\tilde{g}(k) = -2n^2\left(i\,\frac{d\tilde{g}}{dk}\right) = -2in^2\frac{d\tilde{g}}{dk}.
\]
Dividing both sides by $-2in^2$ yields a separable first-order ODE for $\tilde{g}(k)$:
\[
\frac{d\tilde{g}}{dk} = -\frac{k}{2n^2}\,\tilde{g}(k).
\]
Solving via the integrating factor $e^{k^2/(4n^2)}$,
\[
\frac{d}{dk}\!\left(\tilde{g}(k)\,e^{k^2/(4n^2)}\right) = 0
\implies
\tilde{g}(k) = \tilde{g}(0)\,\exp\!\left(-\frac{k^2}{4n^2}\right).
\]

To determine the integration constant $\tilde{g}(0)$, substitute $u = nx$ ($du = n\,dx$, using $n > 0$):
\[
\tilde{g}(0) = \int_{-\infty}^{\infty} e^{-n^2 x^2}\,dx = \frac{1}{n}\int_{-\infty}^{\infty} e^{-u^2}\,du = \frac{\sqrt{\pi}}{n},
\]
where the standard Gaussian integral $I = \int_{-\infty}^{\infty} e^{-u^2}\,du > 0$ satisfies
\[
I^2 = \iint_{\mathbb{R}^2} e^{-(u^2 + v^2)}\,du\,dv
= \int_0^{2\pi}\int_0^{\infty} e^{-r^2}\,r\,dr\,d\theta
= 2\pi\left[-\frac{1}{2}e^{-r^2}\right]_0^{\infty} = \pi
\implies I = \sqrt{\pi}.
\]
(For general real $n \ne 0$, substituting $u = |n|x$ yields $\tilde{g}(0) = \sqrt{\pi}/|n|$.)

Hence, for $n > 0$,
\[
\tilde{g}(k) = \frac{\sqrt{\pi}}{n}\,\exp\!\left(-\frac{k^2}{4n^2}\right).
\]
Now $f(x) = \exp[-n^2(x - \mu)^2] = g(x - \mu)$. Applying the shift property gives
\[
\tilde{f}(k) = e^{-ik\mu}\,\tilde{g}(k) = \frac{\sqrt{\pi}}{n}\,\exp\!\left(-ik\mu - \frac{k^2}{4n^2}\right)
\qquad \left(\text{or } \frac{\sqrt{\pi}}{|n|}\exp\!\left(-ik\mu - \frac{k^2}{4n^2}\right) \text{ for } n \in \mathbb{R}\setminus\{0\}\right).
\]

### 3. The normalised Gaussian $\delta_n(x)$, its transform $\widetilde{\delta}_n(k)$, and sketches

Now set $\mu = 0$ and take $n > 0$:
\[
\delta_n(x) = \frac{n}{\sqrt{\pi}}\,e^{-n^2 x^2}.
\]
By linearity of the Fourier transform,
\[
\widetilde{\delta}_n(k) = \frac{n}{\sqrt{\pi}}\,\tilde{g}(k)
= \frac{n}{\sqrt{\pi}} \cdot \frac{\sqrt{\pi}}{n}\,\exp\!\left(-\frac{k^2}{4n^2}\right)
= \exp\!\left(-\frac{k^2}{4n^2}\right).
\]

#### Key features for sketching $\delta_n(x)$ and $\widetilde{\delta}_n(k)$

1. **$\delta_n(x) = \dfrac{n}{\sqrt{\pi}}e^{-n^2 x^2}$ (even, strictly positive Gaussian bell curve):**
   - Maximum at $x = 0$ with peak height $\delta_n(0) = n/\sqrt{\pi}$.
   - Characteristic width $\sim 1/n$ (specifically, $1/e$-half-width $1/n$, inflection points at $x = \pm 1/(\sqrt{2}\,n)$, standard deviation $\sigma_x = 1/(\sqrt{2}\,n)$).
   - Area under the curve is fixed at $1$ for every $n > 0$.
   - **Small $n$ ($0 < n \ll 1$):** a very **low, broad** bell curve peaking at $n/\sqrt{\pi} \ll 1$ and decaying slowly over $|x| \sim 1/n \gg 1$.
   - **Large $n$ ($n \gg 1$):** a very **tall, narrow** spike peaking at $n/\sqrt{\pi} \gg 1$ and concentrated within $|x| \sim 1/n \ll 1$, practically zero away from $x = 0$.

2. **$\widetilde{\delta}_n(k) = \exp\!\left(-\dfrac{k^2}{4n^2}\right)$ (even, strictly positive real Gaussian bell curve):**
   - Maximum at $k = 0$ with peak height $\widetilde{\delta}_n(0) = 1$, **independent of $n$** (reflecting $\widetilde{\delta}_n(0) = \int_{-\infty}^\infty \delta_n(x)\,dx = 1$).
   - Characteristic width $\sim 2n$ ($1/e$-half-width $2n$, inflection points at $k = \pm \sqrt{2}\,n$, profile standard deviation $\sigma_k = \sqrt{2}\,n$). Notice the reciprocal scaling $\sigma_x \sigma_k = 1$ (Fourier uncertainty relation).
   - **Small $n$ ($0 < n \ll 1$):** a **narrow** peak of height $1$ at $k = 0$, dropping rapidly to $0$ over $|k| \sim 2n \ll 1$.
   - **Large $n$ ($n \gg 1$):** a **very wide, flat-topped** bell curve of height $1$ at $k = 0$, staying close to the constant line $\widetilde{\delta}_n(k) = 1$ across $|k| \ll 2n$ before decaying to $0$ only for $|k| \gg 2n$.

### 4. Evaluation of $\int_{-\infty}^{\infty} \delta_n(x)\,dx$

Using the substitution $u = nx$ ($n > 0$)—or immediately from $\int_{-\infty}^{\infty}\delta_n(x)\,dx = \widetilde{\delta}_n(0)$:
\[
\int_{-\infty}^{\infty} \delta_n(x)\,dx
= \frac{n}{\sqrt{\pi}}\int_{-\infty}^{\infty} e^{-n^2 x^2}\,dx
= \frac{n}{\sqrt{\pi}} \cdot \frac{\sqrt{\pi}}{n} = 1.
\]
(Notice this value $1$ is completely independent of $n > 0$.)

### 5. Behaviour as $n \to \infty$

As $n \to \infty$, $\delta_n(x)$ forms a **nascent delta function** (a Gaussian approximation to the Dirac delta distribution $\delta(x)$):

- **Pointwise behaviour in $x$-space:** For every $x \ne 0$, $\delta_n(x) = (n/\sqrt{\pi})e^{-n^2 x^2} \to 0$ exponentially fast as $n \to \infty$, whereas at $x = 0$, $\delta_n(0) = n/\sqrt{\pi} \to +\infty$, while maintaining constant unit area $\int_{-\infty}^\infty \delta_n(x)\,dx = 1$ and $\delta_n(x) > 0$.
- **Distributional limit (sifting property):** For any bounded continuous test function $\phi(x)$ (in particular any Schwartz or $C_c^\infty$ test function),
  \[
  \lim_{n \to \infty} \int_{-\infty}^{\infty} \delta_n(x)\,\phi(x)\,dx = \phi(0).
  \]
  *Proof:* Given $\varepsilon > 0$, continuity at $0$ gives $\eta > 0$ such that $|\phi(x) - \phi(0)| < \varepsilon/2$ for $|x| < \eta$. Since $\int_{-\infty}^\infty \delta_n(x)\,dx = 1$ and $\delta_n > 0$,
  \begin{align*}
  \left|\int_{-\infty}^{\infty} \delta_n(x)\,\phi(x)\,dx - \phi(0)\right|
  &= \left|\int_{-\infty}^{\infty} \delta_n(x)\bigl(\phi(x) - \phi(0)\bigr)\,dx\right| \\
  &\le \frac{\varepsilon}{2}\int_{|x| < \eta} \delta_n(x)\,dx + 2\|\phi\|_\infty \int_{|x| \ge \eta} \delta_n(x)\,dx \\
  &\le \frac{\varepsilon}{2} + \frac{4\|\phi\|_\infty}{\sqrt{\pi}}\int_{n\eta}^{\infty} e^{-u^2}\,du.
  \end{align*}
  Because $\int_{n\eta}^\infty e^{-u^2}\,du \to 0$ as $n \to \infty$, the right-hand side is $< \varepsilon$ for all sufficiently large $n$. Thus $\delta_n \to \delta$ in the sense of distributions.
- **Behaviour in $k$-space:** For every fixed $k \in \mathbb{R}$ (and uniformly on every bounded interval $[-K, K]$),
  \[
  \lim_{n \to \infty} \widetilde{\delta}_n(k) = \lim_{n \to \infty} \exp\!\left(-\frac{k^2}{4n^2}\right) = 1.
  \]
  This mirrors the reciprocal width scaling of Fourier transforms: as $\delta_n(x)$ becomes infinitely localised at $x = 0$ (approaching $\delta(x)$), its Fourier transform $\widetilde{\delta}_n(k)$ spreads uniformly over all frequencies and converges to the constant function $1 = \tilde{\delta}(k)$, the Fourier transform of the Dirac delta function (consistent with the formal inverse relation $\delta(x) = \frac{1}{2\pi}\int_{-\infty}^\infty e^{ikx}\,dk$).

## Final Answers

- **Fourier transform of $f(x) = \exp[-n^2(x - \mu)^2]$** (under $\tilde{f}(k) = \int_{-\infty}^\infty f(x)e^{-ikx}\,dx$, for $n > 0$):
  \[
  \tilde{f}(k) = \frac{\sqrt{\pi}}{n}\,\exp\!\left(-ik\mu - \frac{k^2}{4n^2}\right)
  \qquad \left(\text{general } n \in \mathbb{R}\setminus\{0\}:\ \frac{\sqrt{\pi}}{|n|}\exp\!\left(-ik\mu - \frac{k^2}{4n^2}\right)\right).
  \]
- **Fourier transform of $\delta_n(x) = (n/\sqrt{\pi})e^{-n^2 x^2}$ ($n > 0$):** $\widetilde{\delta}_n(k) = \exp\!\left(-\dfrac{k^2}{4n^2}\right)$.
- **Integral:** $\displaystyle \int_{-\infty}^{\infty} \delta_n(x)\,dx = 1$ (for all $n > 0$).
- **Sketches:** $\delta_n(x)$ is a unit-area Gaussian centred at $0$ of peak $n/\sqrt{\pi}$ and width $\sim 1/n$ (low/wide for small $n$, tall/narrow for large $n$); $\widetilde{\delta}_n(k)$ is a Gaussian centred at $0$ of constant peak height $1$ and width $\sim 2n$ (narrow for small $n$, broad/flat near $1$ for large $n$).
- **Limit as $n \to \infty$:** $\delta_n(x) \to \delta(x)$ in the sense of distributions (sifting property $\int \delta_n(x)\phi(x)\,dx \to \phi(0)$), while $\widetilde{\delta}_n(k) \to 1 = \tilde{\delta}(k)$ pointwise (and uniformly on bounded $k$-intervals).

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-3_qx24_check.py`

```text
[PASS] Gaussian FT verified across positive/negative n, mu, k (max err=3.34e-15)
[PASS] Fourier inversion formula verified (max err=4.62e-16)
[PASS] Differentiation ODE dg_tilde/dk = -(k/(2n^2)) g_tilde verified (max err=4.92e-11)
[PASS] Area int delta_n(x) dx = 1 (max err=7.77e-16) and delta_n_tilde(k) = exp(-k^2/(4n^2)) (max err=2.11e-15)
[PASS] Distributional limit delta_n -> delta verified on 3 test functions and max |delta_1000_tilde(k) - 1| on [-10,10] = 2.50e-05
ALL CHECKS PASSED for IB_methods_sheet-3_qx24.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. **Fourier transform convention:** Confirm the stored solution states its convention explicitly and uses it consistently! Under $\tilde{f}(k) = \int f(x)e^{-ikx}\,dx$ (confirmed by Q7 on the same Michaelmas 2016 sheet), the shift factor is $e^{-ik\mu}$ and $\widetilde{\delta}_n(k) = e^{-k^2/(4n^2)} \to 1$. (Under the opposite $e^{+ikx}$ sign convention, the shift factor becomes $e^{+ik\mu}$; under the unitary $(2\pi)^{-1/2}$ convention, an extra $(2\pi)^{-1/2}$ appears in $\tilde{f}(k)$ and $\widetilde{\delta}_n(k) \to (2\pi)^{-1/2}$.)
  2. **Method requirement:** The question explicitly asks to use **differentiation and shift properties** (deriving the ODE $d\tilde{g}/dk = -(k/(2n^2))\tilde{g}$ and then shifting by $\mu$), NOT merely completing the square by contour shifting! Check that the existing solution actually uses differentiation and shift properties.
  3. **Sign of $n$:** Check whether the solution assumes $n > 0$ explicitly (or uses $|n|$ for general real $n \ne 0$).
  4. **Peak height of $\widetilde{\delta}_n(k)$ in sketches:** A common sketch mistake is varying the peak height of $\widetilde{\delta}_n(0)$ with $n$; it must remain fixed at $\widetilde{\delta}_n(0) = 1$ (since $\int \delta_n(x)\,dx = 1$), while only its width ($\propto n$) grows!
