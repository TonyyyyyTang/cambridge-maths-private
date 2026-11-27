# IB/methods/sheet-2/qx25

## Source locator
```
archive: private/raw/sheets/IB/methods/2004-05/sheet-2.pdf
pages: [2]
original_question: 8
url: https://web.archive.org/web/20050129090433id_/http://www.damtp.cam.ac.uk:80/user/examples/B8b.pdf
```

## Question figures you may open
- none

## Question (verbatim question.tex)
```latex
\emph{Bessel’s equation.}
\begin{enumerate}
\item[(i)] Show that the eigenvalues of the Sturm-Liouville problem
\[
\frac{d}{dx}\left(x\frac{du}{dx}\right)+\lambda xu=0,\qquad 0<x<1,
\]
with $u(x)$ bounded as $x\to0$ and $u(1)=0$, are $\lambda=j_n^2$ $(n=1,2,\ldots)$, where the $j_n$ are the zeros of the Bessel function $J_0(z)$, arranged in ascending order. [Note: Bessel’s equation of order zero is $\frac{d}{dz}(z\frac{dy}{dz})+zy=0$, $(z>0)$, which you may assume has one solution $J_0(z)$ defined as
\[
J_0(z)=\sum_{m=0}^\infty(-1)^m\frac{z^{2m}}{2^{2m}(m!)^2}
\]
and a second solution that is the sum of a regular function and $J_0(z)\log z$.]
\item[(ii)] Using integration by parts on the differential equations for $J_0(\alpha x)$ and $J_0(\beta x)$, show that
\[
\int_0^1 J_0(\alpha x)J_0(\beta x)x\,dx
=\frac{\beta J_0(\alpha)J_0'(\beta)-\alpha J_0(\beta)J_0'(\alpha)}{\alpha^2-\beta^2}
\qquad(\beta\ne\alpha)
\]
and deduce that
\[
\int_0^1 J_0(j_nx)J_0(j_mx)x\,dx=0,\qquad(n\ne m),
\]
\[
\int_0^1[J_0(j_nx)]^2x\,dx=\frac12[J_0'(j_n)]^2.
\]
[\emph{Hint:} Here, consider $\beta=j_n+\epsilon$ as $\epsilon\to0$.]
\item[(iii)] Assume that the inhomogeneous equation
\[
\frac{d}{dx}\left(x\frac{du}{dx}\right)+\widetilde\lambda xu=xf(x),
\]
where $\widetilde\lambda$ is not an eigenvalue, has a unique solution such that $u(x)$ is bounded as $x\to0$ and $u(1)=0$. Assuming also the completeness of the eigenfunctions $J_0(j_nx)$, obtain the eigenfunction expansion of $u$.
\end{enumerate}
```
