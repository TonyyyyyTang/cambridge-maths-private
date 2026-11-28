---
id: IB/methods/sheet-1/q10
slug: IB_methods_sheet-1_q10
source_archive: private/raw/sheets/IB/methods/2025-26/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against the 200 DPI render `IB_methods_sheet-1_q10_p2.png` (middle of page 2, Question 10, followed by "Additional problems" and Question 11).

- **Result:** `ok` — exact mathematical and verbal match for every equation, derivative sign, index ($m$, $\alpha$, $\beta$, $j_{mk}$, $j_{m\ell}$, $\delta_{k\ell}$), domain ($r \in (0, 1)$), condition ($\alpha \ne \beta$, $J_m(j_{mk}) = 0$, $k = 1, 2, \ldots$), and factor ($\frac{1}{2}$).
- **Purely typographic differences:** In the printed PDF, the final integral identity $\int_0^1 J_m(j_{mk}r)J_m(j_{m\ell}r)r\,\mathrm{d}r = \frac{1}{2}[J_m'(j_{mk})]^2\delta_{k\ell}$ is typeset inline and the sentence ends with `...` (without a trailing period after the ellipsis), whereas `question.tex` places that equation in display mode and adds a trailing full stop after `\ldots`. No mathematical or verbal content is affected.

## Full Solution

### 1. Proof of the Wronskian derivative identity

We are given that $y_\alpha(r) = J_m(\alpha r)$ satisfies
\[
-\frac{\mathrm{d}}{\mathrm{d}r}\!\left(r\frac{\mathrm{d}y_\alpha}{\mathrm{d}r}\right) + \frac{m^2}{r}y_\alpha = \alpha^2 r y_\alpha, \qquad r \in (0, 1),
\]
and analogously $y_\beta(r) = J_m(\beta r)$ satisfies
\[
-\frac{\mathrm{d}}{\mathrm{d}r}\!\left(r\frac{\mathrm{d}y_\beta}{\mathrm{d}r}\right) + \frac{m^2}{r}y_\beta = \beta^2 r y_\beta, \qquad r \in (0, 1).
\]
Multiplying the first equation by $y_\beta(r)$ and the second equation by $y_\alpha(r)$, and subtracting the second from the first, the $\frac{m^2}{r} y_\alpha y_\beta$ terms cancel identically:
\[
y_\alpha(r)\frac{\mathrm{d}}{\mathrm{d}r}\!\left(r\frac{\mathrm{d}y_\beta}{\mathrm{d}r}\right) - y_\beta(r)\frac{\mathrm{d}}{\mathrm{d}r}\!\left(r\frac{\mathrm{d}y_\alpha}{\mathrm{d}r}\right) = (\alpha^2 - \beta^2) r y_\alpha(r) y_\beta(r).
\]
Now expand the derivative of $r(y_\alpha y_\beta' - y_\beta y_\alpha')$ using the product rule:
\[
\big[r(y_\alpha y_\beta' - y_\beta y_\alpha')\big]'
= \big(y_\alpha (r y_\beta') - y_\beta (r y_\alpha')\big)'
= y_\alpha' (r y_\beta') + y_\alpha (r y_\beta')' - y_\beta' (r y_\alpha') - y_\beta (r y_\alpha')'
= y_\alpha (r y_\beta')' - y_\beta (r y_\alpha')',
\]
where the middle terms $r y_\alpha' y_\beta'$ and $-r y_\beta' y_\alpha'$ cancel. Comparing the two equations yields
\[
\big[r(y_\alpha y_\beta' - y_\beta y_\alpha')\big]' = (\alpha^2 - \beta^2) r y_\alpha y_\beta, \qquad r \in (0, 1).
\]

---

### 2. Deduction of the integral for $\alpha \neq \beta$

Integrate both sides over $[\varepsilon, 1]$ for $0 < \varepsilon < 1$:
\[
(\alpha^2 - \beta^2)\int_\varepsilon^1 J_m(\alpha r) J_m(\beta r) r\,\mathrm{d}r
= \Big[r\big(y_\alpha(r) y_\beta'(r) - y_\beta(r) y_\alpha'(r)\big)\Big]_\varepsilon^1.
\]
By the chain rule,
\[
y_\alpha'(r) = \alpha J_m'(\alpha r), \qquad y_\beta'(r) = \beta J_m'(\beta r).
\]
At the upper limit $r = 1$:
\[
1\cdot\big(y_\alpha(1) y_\beta'(1) - y_\beta(1) y_\alpha'(1)\big)
= \beta J_m(\alpha) J_m'(\beta) - \alpha J_m(\beta) J_m'(\alpha).
\]
At the lower limit $r = \varepsilon \to 0^+$:
- For non-negative integer $m \in \mathbb{Z}_{\ge 0}$ (the standard Bessel functions of integer order, and likewise for $m \in \mathbb{Z}_{<0}$ since $J_{-n} = (-1)^n J_n$), $J_m$ is entire (analytic at $0$), so $y_\alpha$, $y_\beta$, $y_\alpha'$, $y_\beta'$ all have finite limits as $r \to 0^+$. Hence $\varepsilon\big(y_\alpha(\varepsilon) y_\beta'(\varepsilon) - y_\beta(\varepsilon) y_\alpha'(\varepsilon)\big) \to 0$ as $\varepsilon \to 0^+$, and the integrand $J_m(\alpha r) J_m(\beta r) r$ is continuous on $[0, 1]$.
- (More generally, for any real $m > -1$, the Frobenius series $J_m(z) = \frac{(z/2)^m}{\Gamma(m+1)}\big(1 - \frac{z^2}{4(m+1)} + O(z^4)\big)$ shows that $y_\alpha(r) y_\beta'(r) - y_\beta(r) y_\alpha'(r) = O(r^{2m+1})$ because the leading $r^{2m-1}$ terms cancel identically, so $r(y_\alpha y_\beta' - y_\beta y_\alpha') = O(r^{2m+2}) \to 0$ as $r \to 0^+$, while $r y_\alpha y_\beta = O(r^{2m+1})$ is integrable at $0$.)

Taking $\varepsilon \to 0^+$ therefore gives
\[
(\alpha^2 - \beta^2)\int_0^1 J_m(\alpha r) J_m(\beta r) r\,\mathrm{d}r
= \beta J_m(\alpha) J_m'(\beta) - \alpha J_m(\beta) J_m'(\alpha).
\]
Dividing by $\alpha^2 - \beta^2$ (assuming $\alpha^2 \neq \beta^2$, which holds for distinct $\alpha, \beta > 0$, or more generally $\alpha \neq \pm\beta$) yields
\[
\int_0^1 J_m(\alpha r) J_m(\beta r) r\,\mathrm{d}r
= \frac{\beta J_m(\alpha) J_m'(\beta) - \alpha J_m(\beta) J_m'(\alpha)}{\alpha^2 - \beta^2}, \qquad \alpha \neq \beta.
\]

---

### 3. Orthogonality and normalisation for the positive zeros $j_{mk}$

Let $0 < j_{m1} < j_{m2} < \cdots$ denote the positive zeros of $J_m(z)$, so $J_m(j_{mk}) = 0$ and $j_{mk} > 0$ for $k = 1, 2, \ldots$. (Note that for $m \ge 1$, $z = 0$ is also a zero of $J_m$, which is excluded from $\{j_{mk}\}_{k \ge 1}$: at $z = 0$ the integral vanishes identically, whereas for $m = 1$, $\frac{1}{2}[J_1'(0)]^2 = \frac{1}{8} \neq 0$.)

- **Case $k \neq \ell$:** Since $j_{mk}$ and $j_{m\ell}$ are distinct positive zeros, $j_{mk} \neq j_{m\ell}$ and $j_{mk}^2 - j_{m\ell}^2 \neq 0$. Setting $\alpha = j_{mk}$ and $\beta = j_{m\ell}$ in the integral formula above, and using $J_m(j_{mk}) = J_m(j_{m\ell}) = 0$, we obtain
  \[
  \int_0^1 J_m(j_{mk} r) J_m(j_{m\ell} r) r\,\mathrm{d}r
  = \frac{j_{m\ell}\cdot 0\cdot J_m'(j_{m\ell}) - j_{mk}\cdot 0\cdot J_m'(j_{mk})}{j_{mk}^2 - j_{m\ell}^2}
  = 0.
  \]

- **Case $k = \ell$:** Fix $\alpha = j_{mk} > 0$ (so $J_m(\alpha) = J_m(j_{mk}) = 0$) and let $\beta$ vary near $j_{mk}$ with $\beta \neq \pm j_{mk}$. The first term in the numerator vanishes identically, leaving
  \[
  \int_0^1 J_m(j_{mk} r) J_m(\beta r) r\,\mathrm{d}r
  = \frac{-j_{mk} J_m(\beta) J_m'(j_{mk})}{j_{mk}^2 - \beta^2}
  = \frac{j_{mk} J_m'(j_{mk})}{\beta + j_{mk}} \cdot \frac{J_m(\beta) - J_m(j_{mk})}{\beta - j_{mk}}.
  \]
  Now take the limit $\beta \to j_{mk}$. On the left-hand side, for $m \ge 0$ the integrand $(r, \beta) \mapsto J_m(j_{mk} r) J_m(\beta r) r$ is continuous on the compact rectangle $[0, 1] \times [\frac{1}{2}j_{mk}, \frac{3}{2}j_{mk}]$ (and for general $m > -1$ is dominated by an integrable multiple of $r^{2m+1}$), so we may pass the limit inside the integral. On the right-hand side, $J_m$ is differentiable at $j_{mk} > 0$, so $\frac{J_m(\beta) - J_m(j_{mk})}{\beta - j_{mk}} \to J_m'(j_{mk})$ and $\frac{j_{mk}}{\beta + j_{mk}} \to \frac{j_{mk}}{2 j_{mk}} = \frac{1}{2}$ (using $j_{mk} > 0$). (Equivalently, applying L'Hôpital's rule with respect to $\beta$ to $\frac{-j_{mk} J_m(\beta) J_m'(j_{mk})}{j_{mk}^2 - \beta^2}$ gives $\left.\frac{-j_{mk} J_m'(\beta) J_m'(j_{mk})}{-2\beta}\right|_{\beta = j_{mk}} = \frac{1}{2}[J_m'(j_{mk})]^2$.) Hence
  \[
  \int_0^1 \big[J_m(j_{mk} r)\big]^2 r\,\mathrm{d}r = \frac{1}{2}\big[J_m'(j_{mk})\big]^2.
  \]
  (Note also that $J_m'(j_{mk}) \neq 0$ since every positive zero of $J_m$ is simple: at the ordinary point $z = j_{mk} > 0$ of Bessel's second-order linear ODE, $J_m(j_{mk}) = J_m'(j_{mk}) = 0$ would force $J_m \equiv 0$ by uniqueness.)

Combining the two cases $k \neq \ell$ and $k = \ell$ gives
\[
\int_0^1 J_m(j_{mk} r) J_m(j_{m\ell} r) r\,\mathrm{d}r
= \frac{1}{2}\big[J_m'(j_{mk})\big]^2 \delta_{k\ell}, \qquad k, \ell = 1, 2, \ldots
\]

## Final Answers

All three required identities are proved:
1. $\big[r(y_\alpha y_\beta' - y_\beta y_\alpha')\big]' = (\alpha^2 - \beta^2) r y_\alpha y_\beta$ on $(0, 1)$.
2. $\displaystyle\int_0^1 J_m(\alpha r) J_m(\beta r) r\,\mathrm{d}r = \frac{\beta J_m(\alpha) J_m'(\beta) - \alpha J_m(\beta) J_m'(\alpha)}{\alpha^2 - \beta^2}$ for $\alpha \neq \pm\beta$ (in particular $\alpha \neq \beta$ with $\alpha, \beta > 0$).
3. $\displaystyle\int_0^1 J_m(j_{mk} r) J_m(j_{m\ell} r) r\,\mathrm{d}r = \frac{1}{2}\big[J_m'(j_{mk})\big]^2 \delta_{k\ell}$ for the positive zeros $j_{mk}$ ($k = 1, 2, \ldots$) of $J_m$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_q10_check.py` (exit code 0).

```text
ODE & Wronskian identity passed: max_ode_err=1.08e-14, max_wron_err=2.72e-09.
General alpha != beta integral check passed (max_err = 3.17e-13).
  m=0: first 5 positive zeros = [2.404826, 5.520078, 8.653728, 11.791534, 14.930918]
  m=1: first 5 positive zeros = [3.831706, 7.015587, 10.173468, 13.323692, 16.47063]
  m=2: first 5 positive zeros = [5.135622, 8.417244, 11.619841, 14.795952, 17.959819]
  m=3: first 5 positive zeros = [6.380162, 9.761023, 13.015201, 16.223466, 19.409415]
Zero-orthogonality & norm checks passed for m=0,1,2,3 (max_err = 5.67e-13).
Confirmed edge case: j_{mk} must be positive zeros (at z=0 for m=1, J_1(0)=0 but (1/2)[J_1'(0)]^2 = 1/8 != 0).
ALL CHECKS PASSED FOR IB/methods/sheet-1/q10.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. **Lower boundary term at $r \to 0^+$:** Do not assert $r(y_\alpha y_\beta' - y_\beta y_\alpha')|_{r=0} = 0$ without checking regularity/boundedness of $y_\alpha, y_\alpha'$ at $r = 0$ (for $m \in \mathbb{Z}_{\ge 0}$, $J_m$ is analytic at $0$; for general real $m > -1$, the leading $r^{2m-1}$ terms cancel leaving $O(r^{2m+2}) \to 0$).
  2. **Sign of the Wronskian bracket and numerator:** Because the ODE has a minus sign in $-(r y_\alpha')'$, subtracting the $\beta$ equation from the $\alpha$ equation gives $+y_\alpha(r y_\beta')' - y_\beta(r y_\alpha')' = (\alpha^2 - \beta^2) r y_\alpha y_\beta$, which yields $+\beta J_m(\alpha) J_m'(\beta) - \alpha J_m(\beta) J_m'(\alpha)$ over $\alpha^2 - \beta^2$. Watch for sign flips!
  3. **Positive zeros $j_{mk} > 0$ and $\alpha \neq -\beta$:** Cancelling $j_{mk}/(2j_{mk}) = 1/2$ uses $j_{mk} > 0$ (for $m \ge 1$, $z = 0$ is a zero of $J_m$ where the formula fails for $m = 1$), and dividing by $\alpha^2 - \beta^2$ formally requires $\alpha^2 \neq \beta^2$ (automatic for distinct positive $\alpha, \beta$).
  4. **Chain rule factor in $y_\alpha'(1) = \alpha J_m'(\alpha)$:** Keep track of the distinction between $y_\alpha'(r) = \frac{\mathrm{d}}{\mathrm{d}r} J_m(\alpha r) = \alpha J_m'(\alpha r)$ and $J_m'(z)$.
