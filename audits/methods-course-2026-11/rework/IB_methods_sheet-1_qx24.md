---
id: IB/methods/sheet-1/qx24
slug: IB_methods_sheet-1_qx24
source_archive: private/raw/sheets/IB/methods/2007-08/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-1_qx24_p2.png` (page 2 of `2007-08/sheet-1.pdf`, Question 6). Page 3 does not exist (`IB_methods_sheet-1_qx27_p3.png` absent; Question 7 follows immediately on page 2).

- **Question 6 body:** Exact verbatim match to the printed sheet, including lowercase `fourier transforms`, punctuation (periods after each of the four equations), and wording.
- **Purely typographic differences:**
  1. Printed heading is underlined italic *Functions defined on $|x| < c$.*; stored `question.tex` uses `\emph{Functions defined on $|x|<c$.}`.
  2. The four equations $f(x)=1.$, $f(x)=e^{iax}.$, $f(x)=\sin(ax).$, $f(x)=\cos(ax).$ are printed as indented left-aligned lines rather than a centred `gathered` block.
- **Important standalone context note (section preamble on page 2):** Immediately above Question 6 on page 2, the sheet prints a section header and Fourier transform convention that applies to Questions 6–11:
  > **Fourier transforms**
  > As in the lectures, for a function $f(x)$, $\tilde{f}(k)$ is the Fourier transform of $f(x)$, defined as
  > $$\tilde{f}(k) = \int_{-\infty}^{\infty} f(x) e^{-ikx}\,dx.$$
  Because `question.tex` omits this section-level preamble, Question 6 taken standalone does not specify the Fourier transform convention (sign in the exponent and $2\pi$ normalisation). The solution below uses the sheet's exact convention $\tilde{f}(k) = \int_{-\infty}^{\infty} f(x) e^{-ikx}\,dx$.

## Full Solution

Using the convention stated at the top of the Fourier transforms section of the sheet, the Fourier transform of $f(x)$ is
\[
\tilde{f}(k) = \int_{-\infty}^{\infty} f(x) e^{-ikx}\,dx.
\]
Each function below vanishes for $|x| > c$, so the integral restricts to $(-c, c)$ (the values at the endpoints $x = \pm c$ do not affect the integral). We assume $c > 0$ (if $c \le 0$ the open interval $|x| < c$ is empty and all four transforms vanish identically), and take $a \in \mathbb{R}$ as a real parameter (noting that the algebraic expressions for $k \neq 0, \pm a$ remain valid for complex $a$ as well).

### 1. $f_1(x) = 1$ for $|x| < c$, $0$ elsewhere

\[
\tilde{f}_1(k) = \int_{-c}^{c} e^{-ikx}\,dx.
\]
- For $k \neq 0$:
  \[
  \tilde{f}_1(k) = \left[\frac{e^{-ikx}}{-ik}\right]_{-c}^{c}
  = \frac{e^{-ikc} - e^{ikc}}{-ik}
  = \frac{e^{ikc} - e^{-ikc}}{ik}
  = \frac{2\sin(kc)}{k}.
  \]
- For $k = 0$:
  \[
  \tilde{f}_1(0) = \int_{-c}^{c} 1\,dx = 2c = \lim_{k \to 0} \frac{2\sin(kc)}{k}.
  \]

Thus $\tilde{f}_1(k)$ is continuous (indeed entire) for all $k \in \mathbb{R}$, with
\[
\tilde{f}_1(k) = \begin{cases}
\dfrac{2\sin(kc)}{k}, & k \neq 0, \\[6pt]
2c, & k = 0.
\end{cases}
\]
As a consistency check, $f_1$ is real and even, and $\tilde{f}_1$ is correspondingly real and even.

### 2. $f_2(x) = e^{iax}$ for $|x| < c$, $0$ elsewhere

Observing that $f_2(x) = e^{iax} f_1(x)$, we apply the frequency-shift property (or integrate directly):
\[
\tilde{f}_2(k) = \int_{-c}^{c} e^{iax} e^{-ikx}\,dx = \int_{-c}^{c} e^{-i(k-a)x}\,dx = \tilde{f}_1(k - a).
\]
- For $k \neq a$:
  \[
  \tilde{f}_2(k) = \left[\frac{e^{-i(k-a)x}}{-i(k-a)}\right]_{-c}^{c} = \frac{2\sin\bigl((k-a)c\bigr)}{k - a}.
  \]
- For $k = a$:
  \[
  \tilde{f}_2(a) = \int_{-c}^{c} 1\,dx = 2c = \lim_{k \to a} \frac{2\sin\bigl((k-a)c\bigr)}{k - a}.
  \]

Hence
\[
\tilde{f}_2(k) = \begin{cases}
\dfrac{2\sin\bigl((k-a)c\bigr)}{k - a}, & k \neq a, \\[6pt]
2c, & k = a.
\end{cases}
\]
As a check, for real $a$, $f_2(-x) = e^{-iax} = f_2(x)^*$ is Hermitian-symmetric, so $\tilde{f}_2(k)$ is real for all real $k$.

### 3. $f_3(x) = \sin(ax)$ for $|x| < c$, $0$ elsewhere

Writing $\sin(ax) = \dfrac{e^{iax} - e^{-iax}}{2i}$ and using linearity together with result (2) for $+a$ and $-a$:
\[
\tilde{f}_3(k) = \frac{1}{2i}\bigl(\tilde{f}_1(k - a) - \tilde{f}_1(k + a)\bigr)
= -\frac{i}{2}\bigl(\tilde{f}_1(k - a) - \tilde{f}_1(k + a)\bigr).
\]
- For $k \neq \pm a$:
  \[
  \tilde{f}_3(k) = -i\left[\frac{\sin\bigl((k-a)c\bigr)}{k - a} - \frac{\sin\bigl((k+a)c\bigr)}{k + a}\right]
  = i\left[\frac{\sin\bigl((k+a)c\bigr)}{k + a} - \frac{\sin\bigl((k-a)c\bigr)}{k - a}\right].
  \]
  Combining over the common denominator $k^2 - a^2$ using $\sin\bigl((k \pm a)c\bigr) = \sin(kc)\cos(ac) \pm \cos(kc)\sin(ac)$ gives the equivalent rational-trigonometric form (for $k \neq \pm a$):
  \[
  \tilde{f}_3(k) = \frac{2i\bigl(k\cos(kc)\sin(ac) - a\sin(kc)\cos(ac)\bigr)}{k^2 - a^2}.
  \]
- Special cases:
  - If $a = 0$, then $f_3(x) \equiv 0$, so $\tilde{f}_3(k) = 0$ for all $k$.
  - If $a \neq 0$ and $k = a$: using $\tilde{f}_1(0) = 2c$ and $\tilde{f}_1(2a) = \frac{2\sin(2ac)}{2a} = \frac{\sin(2ac)}{a}$,
    \[
    \tilde{f}_3(a) = -i\left(c - \frac{\sin(2ac)}{2a}\right).
    \]
  - If $a \neq 0$ and $k = -a$:
    \[
    \tilde{f}_3(-a) = i\left(c - \frac{\sin(2ac)}{2a}\right) = -\tilde{f}_3(a).
    \]

As a check, $f_3$ is real and odd (for real $a$), and $\tilde{f}_3$ is purely imaginary and odd in $k$.

### 4. $f_4(x) = \cos(ax)$ for $|x| < c$, $0$ elsewhere

Writing $\cos(ax) = \dfrac{e^{iax} + e^{-iax}}{2}$ and using linearity together with result (2) for $+a$ and $-a$:
\[
\tilde{f}_4(k) = \frac{1}{2}\bigl(\tilde{f}_1(k - a) + \tilde{f}_1(k + a)\bigr).
\]
- For $k \neq \pm a$:
  \[
  \tilde{f}_4(k) = \frac{\sin\bigl((k-a)c\bigr)}{k - a} + \frac{\sin\bigl((k+a)c\bigr)}{k + a}.
  \]
  Combining over the common denominator $k^2 - a^2$ yields the equivalent form (for $k \neq \pm a$):
  \[
  \tilde{f}_4(k) = \frac{2\bigl(k\sin(kc)\cos(ac) - a\cos(kc)\sin(ac)\bigr)}{k^2 - a^2}.
  \]
- Special cases:
  - If $a = 0$, then $f_4 = f_1$, so $\tilde{f}_4(k) = \frac{2\sin(kc)}{k}$ for $k \neq 0$ and $\tilde{f}_4(0) = 2c$.
  - If $a \neq 0$ and $k = \pm a$:
    \[
    \tilde{f}_4(\pm a) = c + \frac{\sin(2ac)}{2a}.
    \]

As a check, $f_4$ is real and even (for real $a$), $\tilde{f}_4$ is real and even in $k$, and Euler's relation $\tilde{f}_4(k) + i\tilde{f}_3(k) = \tilde{f}_2(k)$ holds identically for all $k$.

## Final Answers

With $\tilde{f}(k) = \int_{-\infty}^{\infty} f(x)e^{-ikx}\,dx$ and $c > 0$ (all functions zero for $|x| > c$, and removable singularities filled by continuity):

1. $f(x) = 1$:
   \[
   \tilde{f}(k) = \frac{2\sin(kc)}{k} \quad (k \neq 0), \qquad \tilde{f}(0) = 2c.
   \]
2. $f(x) = e^{iax}$:
   \[
   \tilde{f}(k) = \frac{2\sin\bigl((k-a)c\bigr)}{k - a} \quad (k \neq a), \qquad \tilde{f}(a) = 2c.
   \]
3. $f(x) = \sin(ax)$:
   \[
   \tilde{f}(k) = -i\left[\frac{\sin\bigl((k-a)c\bigr)}{k - a} - \frac{\sin\bigl((k+a)c\bigr)}{k + a}\right] = \frac{2i\bigl(k\cos(kc)\sin(ac) - a\sin(kc)\cos(ac)\bigr)}{k^2 - a^2} \quad (k \neq \pm a),
   \]
   with $\tilde{f}(\pm a) = \mp i\left(c - \dfrac{\sin(2ac)}{2a}\right)$ for $a \neq 0$ (and $\tilde{f}(k) \equiv 0$ if $a = 0$).
4. $f(x) = \cos(ax)$:
   \[
   \tilde{f}(k) = \frac{\sin\bigl((k-a)c\bigr)}{k - a} + \frac{\sin\bigl((k+a)c\bigr)}{k + a} = \frac{2\bigl(k\sin(kc)\cos(ac) - a\cos(kc)\sin(ac)\bigr)}{k^2 - a^2} \quad (k \neq \pm a),
   \]
   with $\tilde{f}(\pm a) = c + \dfrac{\sin(2ac)}{2a}$ for $a \neq 0$ (and $\tilde{f}(k) = \tilde{f}_1(k)$ if $a = 0$).

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx24_check.py` (executed with `python3`, exit code 0).

```text
IB_methods_sheet-1_qx24_check.py: ALL CHECKS PASSED
  Max quadrature error (all 4 transforms, all params/special k): 1.698e-12
  Max split-vs-combined error (sin, cos):                        1.110e-15
  Max Euler identity error (FT(cos)+i*FT(sin)-FT(e^{iax})):     8.882e-16
  Parseval truncation errors (K=600): [0.001062, 0.001061, 0.000683, 0.000378]
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Phase B watchpoints:**
  1. **Fourier transform convention:** Because `question.tex` omits the sheet's section preamble ($\tilde{f}(k) = \int_{-\infty}^\infty f(x)e^{-ikx}\,dx$), verify that the repository solution uses this exact convention (not $1/\sqrt{2\pi}$ or $e^{+ikx}$), and check whether the sign in $\tilde{f}_3(k)$ is $-i\bigl[\frac{\sin((k-a)c)}{k-a} - \frac{\sin((k+a)c)}{k+a}\bigr]$ (a sign error on $\frac{1}{2i} = -\frac{i}{2}$ is a very common slip!).
  2. **Removable singularities & edge cases:** Check whether the solution explicitly addresses $k = 0$ in (1), $k = a$ in (2), $k = \pm a$ (and $a = 0$) in (3)–(4), and $c > 0$.
  3. **Equivalent algebraic forms:** Both the split sinc form and the combined $(k^2 - a^2)^{-1}$ rational form are exact; watch signs carefully when comparing.
