# Jacobian of a Variable Transformation

The `rgeneric` model in INLA requires defining the function which
returns the log prior density for internal scale *θ*. So, the Jacobian
needs to be added to the log prior density.

The precision matrix *Q* of the multivariate Leroux model is specified
as:

*Q* = *Λ* ⊗ {(1−*λ*)*I* + *λ*(*D*−*W*)}

where *Λ* is defined as:

$$
\Lambda =
\begin{pmatrix}
\tau_1 & \frac{\rho}{\sqrt{\tau_1\tau_2}} \\
\frac{\rho}{\sqrt{\tau_1\tau_2}} & \tau_2
\end{pmatrix}
$$

and *I* is the identity matrix, *W* is the weight matrix, *D* is the
diagonal matrix with row-wise sums, *τ* is the precision parameter, *ρ*
is the correlation parameter, and *λ* is the spatial dependence
parameter.

Hyperparameters were transformed to the internal scale *θ* as:
*θ* = log (*τ*)
$$
\theta = \mathrm{logit}(\frac{\rho+1}{2})
$$
*θ* = logit(*λ*).

Each Jacobian is shown below.

## Precision *τ*

Assuming that the prior distribution of the precision *τ* follows a
Gamma distribution with shape parameter *a* and rate parameter *b*, the
probability density function of *τ* is given by:

$$
f(\tau; a, b) = \text{Gamma}(a, b) = \frac{b^a}{\Gamma(a)} \tau^{a-1} e^{-b\tau}, \quad \text{for } \tau \> 0.
$$
Considering the variable transformation *θ* = log (*τ*), the log prior
needs to be provided on the internal scale (*θ*) instead of the
parameter scale (*τ*). Since *τ* = *e*<sup>*θ*</sup>, the Jacobian is
given by:

$$
\frac{\mathrm{d}\tau}{\mathrm{d}\theta} = \frac{\mathrm{d}}{\mathrm{d}\theta}(e^\theta) = e^\theta.
$$

The density function on the internal scale is:

$$
f(\theta; a, b) = f(\tau; a, b) \left\|\frac{d\tau}{d\theta}\right\| = \text{Gamma}(a, b) \cdot e^\theta.
$$

Therefore, the log prior density on the parameter scale has
log (*e*<sup>*θ*</sup>) = *θ* added to “log.prior”.

## Correlation *ρ*

Correlation parameter was transformed as
$\theta = \mathrm{logit}(\frac{\rho+1}{2})$. Since
$\rho = \frac{e^\theta-1}{1+e^\theta}$, the Jacobian is given by:

$$
\begin{align\*}
\left\|\frac{\mathrm{d}\rho}{\mathrm{d}\theta}\right\| &= 
\frac{\mathrm{d}}{\mathrm{d}\theta} \left( \frac{e^{\theta}-1}{1 + e^{\theta}} \right) \\
&=  \frac{ e^{\theta}(1+e^{\theta})-(e^\theta-1)e^{\theta} }{ (1 + e^{\theta})^2 } \\
&=  \frac{ 2e^{\theta} }{ (1 + e^{\theta})^2 } \\
\end{align\*}
$$
Then, below is added to “log.prior”:

$$
\log\left(\frac{ 2e^{\theta} }{ (1 + e^{\theta})^2 }\right) = 
\log(2)+\log(e^\theta)-2\log(1+e^\theta) \text{.}
$$

## Spatial dependence *λ*

Spatial dependence paramter *λ* was transformed as *θ* = logit(*λ*).
Since $\lambda = \frac{1}{1+e^-\theta}$, the Jacobian is given by:

$$
\begin{align\*}
\left\|\frac{ \mathrm{d}\lambda }{ \mathrm{d}\theta}\right\| &= 
\frac{ \mathrm{d} }{ \mathrm{d}\theta } \left( \frac{ 1 }{ 1+e^{-\theta} } \right) \\
&= \frac{ e^{-\theta} }{ (1+e^{-\theta})^2 } \\
&= \frac{ 1 }{ 1+e^{-\theta} } \frac {e^{-\theta}}{1+e^{-\theta}}     \\
&= \lambda(1-\lambda) \\
\end{align\*}
$$
Then, below is added to “log.prior”:

log (*λ*(1−*λ*)) = log (*λ*) + log (1−*λ*).
