# IB/methods/sheet-2/qx10

## Source locator
```
archive: private/raw/sheets/IB/methods/2026-27/sheet-2.pdf
pages: [2]
original_question: 10
url: https://www.damtp.cam.ac.uk/user/examples/B8b.pdf
```

## Question figures you may open
- none

## Question (verbatim question.tex)
```latex
Let $\Omega\subset\mathbf R^3$ be a bounded domain and $(\mathbf x,t)\in\Omega\times(0,\infty)$. We will be concerned with the following initial-boundary value problems for the heat and wave equations, respectively:
\[
(A)\quad\begin{cases}\varphi_t-\kappa\Delta\varphi=0,&\text{in }\Omega\times(0,\infty)\\
\varphi=f,&\text{on }\Omega\times\{t=0\}\\
\varphi=0,&\text{on }\partial\Omega\times[0,\infty)\end{cases}
\]
\[
(B)\quad\begin{cases}\varphi_{tt}-c^2\Delta\varphi=0,&\text{in }\Omega\times(0,\infty)\\
\varphi=g,&\text{on }\Omega\times\{t=0\}\\
\varphi_t=h,&\text{on }\Omega\times\{t=0\}\\
\varphi=0,&\text{on }\partial\Omega\times[0,\infty)\end{cases}
\]
You may \emph{assume} the following: there is a collection $\{(\psi_n,\lambda_n)\}_{n=1}^\infty$ of real eigenfunction-eigenvalue pairs such that (a) $-\Delta\psi_n=\lambda_n\psi_n$ in $\Omega$; (b) $\psi_n=0$ on $\partial\Omega$; (c) each eigenvalue has finite multiplicity; (d) $\{\psi_n\}$ are complete on $\Omega$. The latter means for $f:\Omega\to\mathbf R$ satisfying $f=0$ on $\partial\Omega$ we can write $f=\sum_n\alpha_n\psi_n$ for some $\{\alpha_n\}$.
\begin{enumerate}
\item[(i)] Show that $\lambda_n>0$ for each $n$ and $\int_\Omega\psi_n\psi_m\,dV=0$ for $\lambda_n\ne\lambda_m$.
\item[(ii)] Explain why we can assume, without loss of generality, that $\int_\Omega\psi_n\psi_m\,dV=0$ for $n\ne m$.
\item[(iii)] Using separation of variables, show that the solution to $(A)$ is given by
\[
\varphi(\mathbf x,t)=\sum_{n=1}^\infty\alpha_n e^{-\lambda_n\kappa t}\psi_n(\mathbf x)\quad\text{where}\quad\alpha_n=\frac{\int_\Omega f\psi_n\,dV}{\int_\Omega\psi_n^2\,dV}.
\]
Explain why this might be formally interpreted as $\varphi(\mathbf x,t)=e^{\kappa t\Delta}\varphi(\mathbf x,0)$ where $e^{\kappa t\Delta}=\sum_{p=0}^\infty\frac{(\kappa t)^p}{p!}\Delta^p$.
\item[(iv)] Solve $(B)$, again using separation of variables. Relate your answer to the formal expression
\[
\varphi(\mathbf x,t)=\frac{\sin\left(ct\sqrt{-\Delta}\right)}{c\sqrt{-\Delta}}\varphi_t(\mathbf x,0)+\cos\left(ct\sqrt{-\Delta}\right)\varphi(\mathbf x,0).
\]
\end{enumerate}
```
