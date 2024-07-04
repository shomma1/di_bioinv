------------------------------------------------------------------------

title: “Jacobian of a Variable Transformation” output: github_document

------------------------------------------------------------------------

# Jacobian of a Variable Transformation

The `rgeneric` model in INLA requires defining the function which
returns the log prior density for internal scale $\theta$. So, the
Jacobian needs to be added to the log prior density.

The precision matrix $Q$ of the multivariate Leroux model is specified
as:

$$
Q=\Lambda\otimes\{(1-\lambda)I+\lambda(D-W)\}
$$

where $\Lambda$ is defined as:

$$
\Lambda =
\begin{pmatrix}
\tau_1 & \frac{\rho}{\sqrt{\tau_1\tau_2}} \\
\frac{\rho}{\sqrt{\tau_1\tau_2}} & \tau_2
\end{pmatrix}
$$

and $I$ is the identity matrix, $W$ is the weight matrix, $D$ is the
diagonal matrix with row-wise sums, $\tau$ is the precision parameter,
$\rho$ is the correlation parameter, and $\lambda$ is the spatial
dependence parameter.

Hyperparameters were transformed to the internal scale $\theta$ as:

$$
\theta = \log { ( \tau ) }
$$

$$
\theta = \mathrm{logit}(\frac{\rho+1}{2})
$$

$$
\theta = \mathrm{logit}(\lambda) \text{.}
$$

Each Jacobian is shown below.

## Precision $\tau$

Assuming that the prior distribution of the precision $\tau$ follows a
Gamma distribution with shape parameter $a$ and rate parameter $b$, the
probability density function of $\tau$ is given by:

$$
f(\tau; a, b) = \text{Gamma}(a, b) = \frac{b^a}{\Gamma(a)} \tau^{a-1} e^{-b\tau}, \quad \text{for } \tau > 0.
$$

Considering the variable transformation $\theta = \log(\tau)$, the log
prior needs to be provided on the internal scale ($\theta$) instead of
the parameter scale ($\tau$). Since $\tau = e^\theta$, the Jacobian is
given by:

$$
\frac{\mathrm{d}\tau}{\mathrm{d}\theta} = \frac{\mathrm{d}}{\mathrm{d}\theta}(e^\theta) = e^\theta.
$$

The density function on the internal scale is:

$$
f(\theta; a, b) = f(\tau; a, b) \left|\frac{d\tau}{d\theta}\right| = \text{Gamma}(a, b) \cdot e^\theta.
$$

Therefore, the log prior density on the parameter scale has
$\log(e^\theta) = \theta$ added to “log.prior”.

## Correlation $\rho$

Correlation parameter was transformed as
$\theta = \mathrm{logit}(\frac{\rho+1}{2})$. Since
$\rho = \frac{e^\theta-1}{1+e^\theta}$, the Jacobian is given by:

$$
\begin{align*}
\left|\frac{\mathrm{d}\rho}{\mathrm{d}\theta}\right| &= 
\frac{\mathrm{d}}{\mathrm{d}\theta} \left( \frac{e^{\theta}-1}{1 + e^{\theta}} \right) \\
&=  \frac{ e^{\theta}(1+e^{\theta})-(e^\theta-1)e^{\theta} }{ (1 + e^{\theta})^2 } \\
&=  \frac{ 2e^{\theta} }{ (1 + e^{\theta})^2 } \\
\end{align*}
$$

Then, below is added to “log.prior”:

$$
\log\left(\frac{ 2e^{\theta} }{ (1 + e^{\theta})^2 }\right) = 
\log(2)+\log(e^\theta)-2\log(1+e^\theta) \text{.}
$$

## Spatial dependence $\lambda$

Spatial dependence paramter $\lambda$ was transformed as
$\theta = \mathrm{logit}(\lambda)$. Since
$\lambda = \frac{1}{1+e^{-\theta}}$, the Jacobian is given by:

$$
\begin{align*}
\left|\frac{ \mathrm{d}\lambda }{ \mathrm{d}\theta}\right| &= 
\frac{ \mathrm{d} }{ \mathrm{d}\theta } \left( \frac{ 1 }{ 1+e^{-\theta} } \right) \\
&= \frac{ e^{-\theta} }{ (1+e^{-\theta})^2 } \\
&= \frac{ 1 }{ 1+e^{-\theta} } \frac {e^{-\theta}}{1+e^{-\theta}}     \\
&= \lambda(1-\lambda) \\
\end{align*} \text{.}
$$

Then, below is added to “log.prior”:

$$
\log\left(\lambda(1-\lambda)\right) = 
\log(\lambda)+\log(1-\lambda) \text{.}
$$
