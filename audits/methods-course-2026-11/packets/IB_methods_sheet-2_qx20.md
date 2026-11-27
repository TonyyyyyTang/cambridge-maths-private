# IB/methods/sheet-2/qx20

## Source locator
```
archive: private/raw/sheets/IB/methods/2007-08/sheet-2.pdf
pages: [1]
original_question: 4
url: https://web.archive.org/web/20081001012634id_/http://www.damtp.cam.ac.uk:80/user/examples/B8b.pdf
```

## Question figures you may open
- none

## Question (verbatim question.tex)
```latex
\emph{Higher order initial value problem$^*$.} Show that the Green’s function for the initial value problem $({}'\equiv\frac{d}{dt})$
\[
y''''+k^2y''=f(t),\qquad y(0)=y'(0)=y''(0)=y'''(0)=0,
\]
is given by
\[
G(t,\tau)=\begin{cases}0,&0<t<\tau,\\k^{-2}(t-\tau)-k^{-3}\sin k(t-\tau),&t>\tau.\end{cases}
\]
Hence write down the solution for $f(t)=e^{-t}$ and verify that it satisfies the equation and the initial conditions.

[\emph{Hint:} Make life easy by noting $G(\tau,\tau)=0$ for an IVP Green’s function and so use the time invariance of the equation to take $G(t,\tau)=\mathrm{fn}(t-\tau)$ for $t>\tau$.]
```
