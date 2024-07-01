---
title: |

  ENEL445: Engineering Optimization

  Source Geo-Location

  Final Report
---

Submitted by:

Daegan Rose-Love (33480267)

Ben Mangin (92699294)

Date: 31/05/2024

# Table of Contents {#table-of-contents .TOC-Heading .unnumbered}

[1. Background [3](#background)](#background)

[1.1. Brief [3](#brief)](#brief)

[1.2. Math Background [3](#math-background)](#math-background)

[1.3 Assumptions [4](#assumptions)](#assumptions)

[2. FDOA-Only Objective Function Formulation
[5](#fdoa-only-objective-function-formulation)](#fdoa-only-objective-function-formulation)

[3. Grid Search [6](#grid-search)](#grid-search)

[4. Find FDOA With Unknown Frequency
[8](#find-fdoa-with-unknown-frequency)](#find-fdoa-with-unknown-frequency)

[5. Gradient Derivation
[10](#gradient-derivation)](#gradient-derivation)

[6. Derivative Verification
[12](#derivative-verification)](#derivative-verification)

[7. Gradient Based algorithm
[14](#gradient-based-algorithm)](#gradient-based-algorithm)

[8. Monty Carlo Comparison
[17](#monty-carlo-comparison)](#monty-carlo-comparison)

[Appendix A: Levenberg-Marquardt Update Derivation
[19](#appendix-a-levenberg-marquardt-update-derivation)](#appendix-a-levenberg-marquardt-update-derivation)

[Appendix B: Algorithm and Grid Search Code.
[21](#appendix-b-algorithm-and-grid-search-code.)](#appendix-b-algorithm-and-grid-search-code.)

# Background

## Brief

Radio source geo-location involves locating the coordinates of where a
particular signal is broadcast by using arbitrary measurements. This
project uses four satellites to locate a ground-based signal using the
frequency-difference of arrival (FDOA) measurements. Source geo-location
has a range of applications, from military to civil, including
surveillance, navigation, and search and rescue.

## Math Background

When calculating the Doppler shift of the signal the coordinates need to
be in the Earth-centred Earth Fixed (ECEF) coordinate system as shown in
vector (1) as x, y, and z. The position and velocity of the satellites
are already in the ECEF system which is shown in vectors (1.2) and (1.3)
but the grid points in the search area (grid points) are geodetic
positions as shown in vector (1.1) where B is latitude, L is longitude,
and H is altitude.

$$\begin{array}{r}
{\overrightarrow{u}}^{0} = \ \begin{bmatrix}
x^{0} & y^{0} & z^{0} \\
\end{bmatrix}^{T}\ \#(1) \\
\end{array}$$

$$\begin{array}{r}
{\overrightarrow{p}}^{0} = \begin{bmatrix}
B^{0} & L^{0} & H^{0} \\
\end{bmatrix}^{T}\ \#(1.1) \\
\end{array}$$

$$\begin{array}{r}
{\overrightarrow{s}}_{i} = \ \begin{bmatrix}
s_{x,i} & s_{y,i} & s_{z,i} \\
\end{bmatrix}^{T}\ \#(1.2) \\
\end{array}$$

$$\begin{array}{r}
{\overrightarrow{\dot{s}}}_{i} = \begin{bmatrix}
{\dot{s}}_{x,i} & {\dot{s}}_{y,i} & {\dot{s}}_{z,i} \\
\end{bmatrix}^{T}\ \#(1.3) \\
\end{array}$$

To convert the coordinates the geodetic coordinates, need to be
converted into radians. Once the coordinates were converted to radians
the equatorial radius is found as shown in equation (2).
$\begin{array}{r}
R_{E}\left( B^{0} \right) = \ \frac{R_{0}}{\sqrt{(1 - e^{2}\sin^{2}\left( B^{0} \right)}}\ \#(2) \\
\end{array}\begin{array}{r}
R_{E}\left( B^{0} \right) = \ \frac{R_{0}}{\sqrt{(1 - e^{2}\sin^{2}\left( B^{0} \right)}}\ \#(\ SEQ\ Equation\ \backslash*\ ARABIC\ 1) \\
\end{array}$

This uses the Earth\'s equatorial radius $R_{0} = 6378.137\ km$ and the
Earth\'s eccentricity $e = 0.081819198425$. The resultant can then be
used in the following equations which will give the converted x, y, and
z coordinates.

$$\begin{array}{r}
x^{0} = \left( R_{E}\left( B^{0} \right) + H \right)\cos\left( B^{0} \right)\cos\left( L^{0} \right)\#(2.1) \\
\end{array}$$

$$\begin{array}{r}
y^{0} = \left( R_{E}\left( B^{0} \right) + H \right)\cos\left( B^{0} \right)\sin\left( L^{0} \right)\#(2.2) \\
\end{array}$$

$$\begin{array}{r}
z^{0} = \left\lbrack \left( 1 + e^{2} \right)R_{E}\left( B^{0} \right) + H \right\rbrack\sin\left( B^{0} \right)\#(2.3) \\
\end{array}$$

The equation shown below in (3) calculates the frequency of the signal
received by the satellite. This is not the same frequency as the
transmitted frequency. The transmitted signal will be subjected to
Doppler shift and channel noise. $\begin{array}{r}
f_{i} = \ f_{0} + \frac{1}{\lambda}\frac{\left( {\overrightarrow{u}}^{0} - {\overrightarrow{s}}_{i} \right){\overrightarrow{\dot{s}}}_{i}}{\left\| {\overrightarrow{u}}^{0} - {\overrightarrow{s}}_{i} \right\|} + n_{i}\ \#(3) \\
\end{array}$

The equation uses the frequency of the source and the position vector of
the point that is being measured, the position vector of the satellite,
the velocity vector of the satellite, and additive noise. As the
frequency of the source is unknown this equation can't be used. To get
around this the frequency of the received signal can be compared to the
frequency of the first satellite as shown in the following equation
(3.1). This will result in a vector of frequencies as shown below in
vector (3.2).

$$\begin{array}{r}
\overrightarrow{f} = f_{j} - f_{1}\  = \ \frac{1}{\lambda}\frac{\left( {\overrightarrow{u}}^{0} - {\overrightarrow{s}}_{j} \right){\overrightarrow{\dot{s}}}_{j}}{\left\| {\overrightarrow{u}}^{0} - {\overrightarrow{s}}_{j} \right\|} - \frac{1}{\lambda}\frac{\left( {\overrightarrow{u}}^{0} - {\overrightarrow{s}}_{1} \right){\overrightarrow{\dot{s}}}_{1}}{\left\| {\overrightarrow{u}}^{0} - {\overrightarrow{s}}_{1} \right\|} + (n_{j} - n_{1})\ ,\ for\ j = 2,3,4\ \#(3.1) \\
\end{array}$$

$$\begin{array}{r}
 \rightarrow \overrightarrow{f} = \ \begin{bmatrix}
f_{21} & f_{31} & f_{41} \\
\end{bmatrix}^{T}\ \#(3.2) \\
\end{array}$$

The following vector (3.3) is the noise vector used to add noise to each
frequency vector at each coordinate. Each n value is a random value
chosen within a Gaussian function.

$$\begin{array}{r}
 \rightarrow \overrightarrow{n} = \begin{bmatrix}
n_{2} - n_{1} \\
n_{3} - n_{1} \\
n_{4} - n_{1} \\
\end{bmatrix}\ \#(3.3) \\
\end{array}$$

When creating a contour plot of the objective function the frequency at
each grid point needs to be evaluated. This will be done as shown in the
vector below (4). Where each value is the Doppler shift a satellite
receives minus the Doppler shift of satellite one.

$$\begin{array}{r}
\overrightarrow{g}\left( B^{0},\ L^{0} \right) = \ \begin{bmatrix}
g\left( B^{0},\ L^{0},\ s_{2},{\dot{s}}_{2} \right) - g\left( B^{0},L^{0},s_{1},{\dot{s}}_{1} \right) \\
g\left( B^{0},\ L^{0},\ s_{3},{\dot{s}}_{3} \right) - g\left( B^{0},L^{0},s_{1},{\dot{s}}_{1} \right) \\
g\left( B^{0},\ L^{0},\ s_{4},{\dot{s}}_{4} \right) - g\left( B^{0},L^{0},s_{1},{\dot{s}}_{1} \right) \\
\end{bmatrix}\ \#(4) \\
\end{array}$$

The following equation (5) allows us to find the Q matrix for our
weighting function used to find theta.

$$\begin{array}{r}
Cov\left( x_{i},\ x_{j} \right) = E\left\lbrack \left( x_{i} - µ_{i} \right)\left( x_{j} - µ_{j} \right) \right\rbrack = E\left\lbrack \Delta\overrightarrow{f}\Delta{\overrightarrow{f}}^{T} \right\rbrack\#(5) \\
\end{array}$$

The function for a multivariate Gaussian probability function (PDF) is
shown below.

$$\begin{array}{r}
\mathcal{N}(x;\mu,\ \Sigma) = \frac{1}{\sqrt{2\pi\Sigma}}\exp\left( - \frac{(x - \mu)^{T}\Sigma^{- 1}(x - \mu)}{2} \right)\ \#(6) \\
\end{array}$$

X is the known measurement vector, $\mu$ is the mean vector, and the
$\Sigma$ is the covariance matrix composed of $\sigma^{2}Q$.

## 1.3 Assumptions  {#assumptions .unnumbered}

The assumptions that were made for this simulation are as follows:

-   Altitude is zero and does not change.

-   Measurements do not account for special relativity.

-   The source frequency is known.

-   The location of the source is known for the grid-search and
    Monte-Carlo simulation.

# FDOA-Only Objective Function Formulation

For convenience, set the noise vector found in (3.3) to be:

$$\begin{array}{r}
\Delta\overrightarrow{f} = \ \begin{bmatrix}
n_{2} - n_{1} \\
n_{3} - n_{1} \\
n_{4} - n_{1} \\
\end{bmatrix}\ \#(7) \\
\end{array}$$

Then by using (5) the following is found:

$$\begin{array}{r}
E\left\lbrack \left\lbrack \begin{matrix}
n_{2} - n_{1} \\
n_{3} - n_{1} \\
n_{4} - n_{1} \\
\end{matrix}\  \right\rbrack\begin{bmatrix}
n_{2} - n_{1} & n_{3} - n_{1} & n_{4} - n_{1} \\
\end{bmatrix} \right\rbrack\ \#(7.1) \\
\end{array}$$

After expanding and removing cross terms we are left with the
uncorrelated noise (7.2).

$$\begin{array}{r}
 \rightarrow E\begin{bmatrix}
\left( n_{2}^{2} + n_{1}^{2} \right) & n_{1}^{2} & n_{1}^{2} \\
n_{1}^{2} & \left( n_{3}^{2} + n_{1}^{2} \right) & n_{1}^{2} \\
n_{1}^{2} & n_{1}^{2} & \left( n_{4}^{2} + n_{1}^{2} \right) \\
\end{bmatrix}\ \#(7.2) \\
\end{array}$$

Notice that each noise term is ∝ $\sigma_{f}^{2}$ resulting in the
following:

$$\begin{array}{r}
\Sigma = \begin{pmatrix}
2\sigma_{f}^{2} & \sigma_{f}^{2} & \sigma_{f}^{2} \\
\sigma_{f}^{2} & 2\sigma_{f}^{2} & \sigma_{f}^{2} \\
\sigma_{f}^{2} & \sigma_{f}^{2} & 2\sigma_{f}^{2} \\
\end{pmatrix} = \ \sigma_{f}^{2}\begin{pmatrix}
2 & 1 & 1 \\
1 & 2 & 1 \\
1 & 1 & 2 \\
\end{pmatrix} = \ \sigma_{f}^{2}Q\ \#(7.3) \\
\end{array}\ $$

$$\begin{array}{r}
Q = \begin{pmatrix}
2 & 1 & 1 \\
1 & 2 & 1 \\
1 & 1 & 2 \\
\end{pmatrix}\ \#(7.4) \\
\end{array}$$

To solve this, we are interested in finding unknown parameters on which
our measurements depend. By using maximum likelihood estimation (MLE) we
locate the values that assign the highest probability to the
measurements. The MLE is defined as follows:

$$\begin{array}{r}
\widehat{\overrightarrow{\theta}} = \max_{\theta}{p\left( \mathcal{D} \middle| \overrightarrow{\theta} \right)}\ \#(8) \\
\end{array}\ $$

Where $p\left( \mathcal{D} \middle| \overrightarrow{\theta} \right)$ is
the probability distribution function of 𝒟 given
$\overrightarrow{\theta}$. When the measurements in
$\mathcal{D = \lbrack}\overrightarrow{y_{1}},\ \overrightarrow{y_{2}},\ldots,\ \overrightarrow{y_{n}}\rbrack$
are identically distributed and independent, the objective function
becomes (8.1). Taking the log of both sides gives:

$$\begin{array}{r}
\log{(p\left( \mathcal{D} \middle| \overrightarrow{\theta} \right)) = \ \log\left( \Sigma_{i = 1}^{m}p\left( \overrightarrow{y_{i}} \middle| \overrightarrow{\theta} \right) \right)}\ \#(8.1) \\
\end{array}$$

In the case where G is a measurement matrix and ϵ ∼ 𝒩(ϵ; 0, Σ),
measurements $\overrightarrow{y}$ can be written as:

$$\begin{array}{r}
\overrightarrow{y} = \left\lbrack y_{1},\ y_{2},\ \ldots,\ y_{m} \right\rbrack^{T} = G\overrightarrow{\theta} + \overrightarrow{\epsilon}\ \ \#(8.2) \\
\end{array}$$

The measurements thus follow a Gaussian distribution
$\mathcal{N}\left( \overrightarrow{y};G\overrightarrow{\theta},\ \Sigma \right)$
as shown in (6).

$$\begin{array}{r}
p\left( \overrightarrow{y} \middle| \overrightarrow{\theta} \right) = \mathcal{N}\left( \overrightarrow{y};G\overrightarrow{\theta},\ \Sigma \right)\#(8.3) \\
\end{array}$$

The log of the multivariate Gaussian distribution removes the
exponential from the equation giving the following equation.

$$\begin{array}{r}
\log\left( p\left( \overrightarrow{y} \middle| \overrightarrow{\theta} \right) \right) \propto - \left( \overrightarrow{y} - G\overrightarrow{\theta} \right)^{T}\Sigma^{- 1}\left( \overrightarrow{y} - G\overrightarrow{\theta} \right)\#(8.4) \\
\end{array}$$

To find the signal source the maximum needs to be found. The equation
can be further simplified by removing the negative sign and finding the
minimum. This is shown in the equation below.

$$\begin{array}{r}
\widehat{\theta} = \max_{\overrightarrow{\theta}}{( - \log{\left( p\left( \overrightarrow{y} \middle| \overrightarrow{\theta} \right) \right))}} = \min_{\overrightarrow{\theta}}{\left( \overrightarrow{y} - \overrightarrow{g}\left( \overrightarrow{\theta} \right) \right)^{T}\Sigma^{- 1}\left( \overrightarrow{y} - \overrightarrow{g}\left( \overrightarrow{\theta} \right) \right)}\ \#(8.5) \\
\end{array}$$

The measurement vector in this case is $\overrightarrow{f}$ , which
gives.

$$\begin{array}{r}
 \rightarrow \widehat{\overrightarrow{\theta}} = \min_{\overrightarrow{\theta}}{\left( \overrightarrow{f} - \overrightarrow{g}\left( \overrightarrow{\theta} \right) \right)^{T}\Sigma^{- 1}\left( \overrightarrow{f} - \overrightarrow{g}\left( \overrightarrow{\theta} \right) \right)}\ \#(9) \\
\end{array}\begin{array}{r}
 \rightarrow \widehat{\overrightarrow{\theta}} = \min_{\overrightarrow{\theta}}{\left( \overrightarrow{f} - \overrightarrow{g}\left( \overrightarrow{\theta} \right) \right)^{T}\Sigma^{- 1}\left( \overrightarrow{f} - \overrightarrow{g}\left( \overrightarrow{\theta} \right) \right)}\ \#(8.6) \\
\end{array}$$

The objective function in this form is in fact a non-linear least
squares (NLS) problem.

# Grid Search 

To develop a map of the objective function (the nonlinear least squares
problem) the Doppler shift of the signal received from each satellite
needed to be found. This was done using Equation (3). As the frequency
of the source is not known the difference in Doppler shift between
satellites 1 and 2, 3, and 4 uses Equation 3.1, this also allows for the
generation of the source vector at (5, 10, 0) as shown in Vector 3.2. To
get the measured frequency vector for each grid point the coordinates of
the 40 by 40 grid were iterated through and solved as shown in Vector 4.
At each coordinate the expected frequency vector and the measured
frequency vector were substituted into Equation (6) to give. This would
then be recorded into a matrix that was then plotted as a contour map
(using Python) to get a visual representation of the objective function.
The code that calculates for each plot is shown in the Python code in
Appendix A. The contour plot is shown in Figure 1 below.

![](./image1.png){width="4.173332239720035in" height="3.13in"}

Figure 1: Contour plot of the signal with no noise.

Figure 1 shows that there is a minimum at $g(5,10)$ which is the
location of the source. Although this contour graph does not have noise,
the following contour map in Figure 2 does. The noise added to the
simulation for Figure 2 was 1Hz which resulted in the minimum being
found 600 meters away from the true minimum.

![A chart with a gradient of green and blue Description automatically
generated](./image2.png){width="4.177419072615923in"
height="3.1330643044619424in"}

Figure 2: Contour plot of the signal with 1Hz of noise.

To find the minimum within the data a grid search algorithm was
implemented. This is done by iterating over each grid point and checking
if the new grid point is smaller. Once every grid point is checked, the
coordinates at that minimum are returned. This algorithm was implemented
in the function for calculating as it already steps through each
coordinate with a step size of 0.01. This allows a 0.01-degree
accuracy/resolution. The algorithm is shown in the Python code in
Appendix B. The grid search algorithm was tested with different noise
levels. The results of this test are shown in Table 1.

  ----------------------------------------------------------------------------------------------------------------------------------------------
  Nosie deviation (Hz)                Minima
                                      ($\mathbf{B}^{\mathbf{0}}\mathbf{,\ }\mathbf{L}^{\mathbf{0}}\mathbf{,}\mathbf{H}^{\mathbf{0}}\mathbf{)}$
  ----------------------------------- ----------------------------------------------------------------------------------------------------------
  0                                   (5, 10, 0)

  1                                   (5, 10, 0)

  2                                   (5, 10, 0)

  4                                   (5, 10.1, 0)

  8                                   (4.9, 9.9, 0)

  16                                  (5.1, 10.1, 0)
  ----------------------------------------------------------------------------------------------------------------------------------------------

  : Table 1: Minima found with different noise levels.

# Find FDOA With Unknown Frequency

In a real-world situation, the original frequency of the transmitter is
unknown which means that the wavelength of the signal cannot be found.
To get around this an estimation of the original frequency needs to be
made to get an estimate of the wavelength. This can be obtained by
evaluating the frequency at satellite j using the following equation.

$$\begin{array}{r}
f_{j} = \frac{1}{\lambda}\frac{\left( u^{0} - s_{i} \right)^{T}\dot{s_{i}}}{\left\| u^{0} - s_{i} \right\|} + f_{0},\ \ for\ j = 1,2,3,4\#(10) \\
\end{array}$$

The frequencies at each satellite can then be averaged to get an
estimate of the frequency as shown in equation 11.

$$\begin{array}{r}
f_{e} = \frac{f_{1} + f_{2} + f_{3} + f_{4}}{4}\#(11) \\
\end{array}$$

Using this estimate of the frequency the wavelength can be found as
shown below where $c$ is the speed of light and $\lambda$ is the
wavelength.

$$\begin{array}{r}
\frac{1}{\lambda} = \frac{c}{f_{e}}\#(12) \\
\end{array}$$

This can be rearranged to give lambda as shown in equation 12.1.

$$\begin{array}{r}
\lambda = \frac{f_{e}}{c}\#(12.1) \\
\end{array}$$

This can then be substituted back into equation 3.1 to give the
frequency difference between satellites as shown in equation 13.

$$\begin{array}{r}
f_{j1} = \frac{f_{e}}{c}\ \frac{\left( \overrightarrow{u^{0}} - \overrightarrow{s_{j}} \right)^{T}\dot{s_{j}}}{\left\| \overrightarrow{u^{0}} - \overrightarrow{s_{j}} \right\|} - \frac{f_{e}}{c}\ \frac{\left( \overrightarrow{u^{0}} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left\| \overrightarrow{u^{0}} - \overrightarrow{s_{1}} \right\|} + \left( n_{j} - n_{1} \right),\ \ for\ j = 2,3,4\ \#(13) \\
\end{array}$$

This now allows for the frequency difference between satellites to be
found so the objective function can be formulated. Results for this
frequency estimation can be seen in Figure 3. This was also verified
using the grid search algorithm to check the minimum was at (5,10) which
it was.

![A chart with a gradient of green and blue colors Description
automatically
generated](./image3.png){width="2.9784951881014874in"
height="2.2338713910761157in"}![](./image4.png){width="2.973333333333333in"
height="2.23in"}

Figure 3: Two plots showing the difference between plots of the know
frequency and the estimated frequency.

**Note:** Using more FDOA measurements such as $f_{23},\ f_{13},f_{43}$
does not result in any performance improvements. The additional
measurements are the result of values being shifted but the data in the
matrix retains the same information with a different reference
satellite, this can be shown mathematically:

$$\begin{array}{r}
\begin{pmatrix}
 - 1 & 1 & 0 & 0 \\
 - 1 & 0 & 1 & 0 \\
 - 1 & 0 & 0 & 1 \\
\end{pmatrix}\begin{pmatrix}
f_{1} \\
f_{2} \\
f_{3} \\
f_{4} \\
\end{pmatrix} = \begin{pmatrix}
f_{21} \\
f_{31} \\
f_{41} \\
\end{pmatrix} = \ \begin{pmatrix}
0 & 1 & - 1 & 0 \\
1 & 0 & - 1 & 0 \\
0 & 0 & - 1 & 1 \\
\end{pmatrix}\begin{pmatrix}
f_{1} \\
f_{2} \\
f_{3} \\
f_{4} \\
\end{pmatrix} = \begin{pmatrix}
f_{23} \\
f_{13} \\
f_{43} \\
\end{pmatrix}\ \#(13.1) \\
\end{array}$$

# Gradient Derivation

To apply gradient-based optimization methods, the gradient must be
computed or approximated. From Equation 3.1 we have:

$$\begin{array}{r}
f_{j1} = \frac{1}{\lambda\ }\ \frac{\left( \overrightarrow{u^{0}} - \overrightarrow{s_{j}} \right)^{T}\dot{s_{j}}}{\left\| \overrightarrow{u^{0}} - \overrightarrow{s_{j}} \right\|} - \frac{1}{\lambda\ }\ \frac{\left( \overrightarrow{u^{0}} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left\| \overrightarrow{u^{0}} - \overrightarrow{s_{1}} \right\|}\ ,\ \frac{\partial f_{j1}}{\partial{\overrightarrow{u}}^{0}} = 0\ \#(14) \\
\end{array}$$

$\ $Let

$$\begin{array}{r}
\overrightarrow{a} = {\overrightarrow{u}}^{0} - \overrightarrow{s_{j}}\ \#(14.1) \\
\end{array}$$

$$\begin{array}{r}
\overrightarrow{a_{1}} = {\overrightarrow{u}}^{0} - \overrightarrow{s_{1}}\ \#(14.2) \\
\end{array}$$

$$\begin{array}{r}
f_{j1} = \frac{1}{\lambda\ }\ \frac{{\overrightarrow{a}}^{T}\dot{s_{j}}}{\left\| \overrightarrow{a} \right\|} - \frac{1}{\lambda\ }\ \frac{{\overrightarrow{a_{1}}}^{T}\dot{s_{1}}}{\left\| \overrightarrow{a_{1}} \right\|}\ ,\ \frac{\partial f_{j1}}{\partial{\overrightarrow{u}}^{0}},\ \ \ \frac{\partial f_{j1}}{\partial{\overrightarrow{u}}^{0}} = \frac{\partial f_{j1}}{\partial\overrightarrow{a}}\frac{\partial\overrightarrow{a}}{\partial{\overrightarrow{u}}^{0}} - \frac{\partial f_{j1}}{\partial\overrightarrow{a}}\frac{\partial\overrightarrow{a}}{\partial{\overrightarrow{u}}^{0}}\ \#(14.3) \\
\end{array}$$

Let

$$\begin{array}{r}
\left\| \overrightarrow{a} \right\| = b^{\frac{1}{2}} = \ \sqrt{\sum a_{i}^{2}}\ \#(14.4) \\
\end{array}$$

$$\begin{array}{r}
f_{j1} = \frac{1}{\lambda}\frac{\left( {\overrightarrow{a}}^{T}{\dot{s}}_{j} \right)}{b^{\frac{1}{2}}} - \frac{1}{\lambda}\frac{\left( {\overrightarrow{a_{1}}}^{T}{\dot{s}}_{1} \right)}{b^{\frac{1}{2}}}\ \#(14.5) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{j1}}{\partial\overrightarrow{a}} = \ \frac{\partial}{\partial\overrightarrow{a}}\left( \frac{1{{\overrightarrow{a}}_{1}}^{T}{\dot{s}}_{j}b^{- \frac{1}{2}}}{\lambda} \right)\ \#(14.6) \\
\end{array}$$

Using product rule and chain rule gives:

$$\begin{array}{r}
\frac{\partial f_{j1}}{\partial{\overrightarrow{u}}^{0}} = \left( \frac{1}{\lambda}\left( \dot{s_{j}}b^{- \frac{1}{2}} - \frac{a^{T}\dot{s_{j}}}{2b^{\frac{3}{2}}}\frac{\partial b}{\partial\overrightarrow{a}} \right) \right)\frac{\partial\overrightarrow{a}}{\partial{\overrightarrow{u}}^{0}}\ \#(14.7) \\
\end{array}$$

$$\begin{array}{r}
\left( b^{\frac{1}{2}} \right)^{2} = \left( \sqrt{\sum a_{i}^{2}} \right)^{2} = \sum a_{i}^{2} = b,\ \ \frac{\partial b}{\partial a} = \sum 2a_{i}\ \#(14.8) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial\overrightarrow{a}}{\partial{\overrightarrow{u}}^{0}} = \frac{\partial}{\partial{\overrightarrow{u}}^{0}}\left( {\overrightarrow{u}}^{0} - {\overrightarrow{s}}_{j} \right) = \begin{bmatrix}
1 \\
1 \\
1 \\
\end{bmatrix}\ \#(14.9) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{j1}}{\partial{\overrightarrow{u}}^{0}} = \left( \frac{1}{\ \lambda}\left( \frac{\dot{s_{j}}}{b^{\frac{1}{2}}} - \frac{a^{T}\dot{s_{j\ }}}{2b^{\frac{3}{2}}}\left( 2\overrightarrow{a} \right) \right) \right)\begin{bmatrix}
1 \\
1 \\
1 \\
\end{bmatrix} - \left( \frac{1}{\ \lambda}\left( \frac{\dot{s_{1}}}{b^{\frac{1}{2}}} - \frac{a^{T}\dot{s_{1\ }}}{2b^{\frac{3}{2}}}\left( 2{\overrightarrow{a}}_{1} \right) \right) \right)\begin{bmatrix}
1 \\
1 \\
1 \\
\end{bmatrix}\ \#(14.10) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{j1}}{\partial{\overrightarrow{u}}^{0}} = \left( \frac{1}{\lambda}\left( \frac{\dot{s_{j}}}{\left\| \overrightarrow{a} \right\|} - \frac{a^{T}{\dot{s}}_{j}}{2\left\| \overrightarrow{a} \right\|^{3}}\left( 2\overrightarrow{a} \right) \right) \right) - \left( \frac{1}{\lambda}\left( \frac{\dot{s_{1}}}{\left\| \overrightarrow{a_{1}} \right\|} - \frac{a_{1}^{T}{\dot{s}}_{1}}{2\left\| \overrightarrow{a_{1}} \right\|^{3}}\left( 2\overrightarrow{a_{1}} \right) \right) \right)\ \#(14.11) \\
\end{array}$$

Substituting $\overrightarrow{a}$ back in gives the derivative:

$$\begin{array}{r}
\frac{\partial f_{j1}}{\partial{\overrightarrow{u}}^{0}} = \ \ \left( \frac{1}{\lambda}\left( \frac{\dot{s_{j}}}{\left\| {\overrightarrow{u}}^{0} - \overrightarrow{s_{j}} \right\|} - \frac{\left( {\overrightarrow{u}}^{0} - \overrightarrow{s_{j}} \right)^{T}{\dot{s}}_{j}}{\left\| {\overrightarrow{u}}^{0} - \overrightarrow{s_{j}} \right\|^{3}}{(\overrightarrow{u}}^{0} - \overrightarrow{s_{j}}) \right) \right) - \left( \frac{1}{\lambda}\left( \frac{\dot{s_{1}}}{\left\| {\overrightarrow{u}}^{0} - \overrightarrow{s_{1}} \right\|} - \frac{\left( {\overrightarrow{u}}^{0} - \overrightarrow{s_{1}} \right)^{T}{\dot{s}}_{1}}{\left\| {\overrightarrow{u}}^{0} - \overrightarrow{s_{1}} \right\|^{3}}\left( {\overrightarrow{u}}^{0} - \overrightarrow{s_{1}} \right) \right) \right)\ \#(14.12) \\
\end{array}$$

The resulting Jacobian can be written in the following form:

$$\begin{array}{r}
J_{\overrightarrow{f}}\left( \overrightarrow{u} \right) = \begin{pmatrix}
\frac{\partial f_{21}}{\partial\overrightarrow{x}} & \frac{\partial f_{21}}{\partial\overrightarrow{y}} & \frac{\partial f_{21}}{\partial\overrightarrow{z}} \\
\frac{\partial f_{31}}{\partial\overrightarrow{x}} & \frac{\partial f_{31}}{\partial\overrightarrow{y}} & \frac{\partial f_{31}}{\partial\overrightarrow{z}} \\
\frac{\partial f_{41}}{\partial\overrightarrow{x}} & \frac{\partial f_{41}}{\partial\overrightarrow{y}} & \frac{\partial f_{41}}{\partial\overrightarrow{z}} \\
\end{pmatrix} = \ \begin{pmatrix}
\nabla f_{21}\left( \overrightarrow{u} \right)^{T} \\
\nabla f_{31}\left( \overrightarrow{u} \right)^{T} \\
\nabla f_{41}\left( \overrightarrow{u} \right)^{T} \\
\end{pmatrix}\ \#(14.13) \\
\end{array}$$

The derivatives of the coordinate conversion functions (Eq. 2, 2.1, 2.2,
2.3, 2.4) are used to find the derivative of $f_{j1}$ with respect to
$\overrightarrow{p}$ (Eq. 1.1). The resulting Jacobian
$J_{\overrightarrow{f}}(\overrightarrow{p})$ is shown below, with
derived terms.

$$\begin{array}{r}
J_{f}\left( \overrightarrow{p} \right) = \ \begin{pmatrix}
\frac{\partial x}{\partial B} & \frac{\partial x}{\partial L} & \frac{\partial x}{\partial H} \\
\frac{\partial y}{\partial B} & \frac{\partial y}{\partial L} & \frac{\partial y}{\partial H} \\
\frac{\partial z}{\partial B} & \frac{\partial z}{\partial L} & \frac{\partial z}{\partial H} \\
\end{pmatrix}\ \#(15) \\
\end{array}$$

The desired Jacobian for
$\frac{\partial\overrightarrow{f}}{\partial\overrightarrow{p}}$ is:

$$\begin{array}{r}
J_{f}\left( \overrightarrow{u},\ \overrightarrow{p} \right) = J_{\overrightarrow{f}}\left( \overrightarrow{u} \right)J_{\overrightarrow{f}}\left( \overrightarrow{p} \right)\#(15.1) \\
\end{array}$$

The derivation results follow:

$$\begin{array}{r}
\frac{\partial R_{e}}{\partial B} = \frac{R_{0}e^{2}\cos(B)\sin(B)}{(1 - e^{2}{\sin^{2}{(B))}}^{\frac{2}{3}}}\ \#(15.2) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial x^{0}}{\partial B^{0}} = \cos\left( L^{0} \right)\left\lbrack \frac{\partial R_{e}}{\partial B}\cos\left( B^{0} \right) - R_{e}\left( B^{0} \right)\sin\left( B^{0} \right) \right\rbrack\ \#(15.3) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial x^{0}}{\partial L^{0}} = R_{e}\left( B^{0} \right)\cos\left( B^{0} \right)\sin\left( L^{0} \right)\ \#(15.4) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial y^{0}}{\partial B^{0}} = \sin\left( L^{0} \right)\left\lbrack \frac{\partial R_{e}}{\partial B}\cos\left( B^{0} \right) - R_{e}\left( B^{0} \right)\sin\left( B^{0} \right) \right\rbrack\ \#(15.5) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial y^{0}}{\partial L^{0}} = R_{e}\left( B^{0} \right)\cos\left( B^{0} \right)\cos\left( L^{0} \right)\ \#(15.6) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial z^{0}}{\partial B^{0}} = \left( 1 - e^{2} \right)\left( \frac{\partial R_{e}}{\partial B}\sin\left( B^{0} \right) + R_{e}\left( B^{0} \right)\cos(B) \right)\ \#(15.7) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial z^{0}}{\partial L^{0}} = 0\ \#(15.8) \\
\end{array}$$

$$\begin{array}{r}
J_{f}\left( \overrightarrow{p} \right) = \begin{bmatrix}
\frac{\partial x^{0}}{\partial B^{0}} & \frac{\partial x^{0}}{\partial L^{0}} \\
\frac{\partial y^{0}}{\partial B^{0}} & \frac{\partial y^{0}}{\partial L^{0}} \\
\frac{\partial z^{0}}{\partial B^{0}} & \frac{\partial z^{0}}{\partial L^{0}} \\
\end{bmatrix}\#(15.9) \\
\end{array}$$

Note that the H column is omitted as the altitude in degrees (H) is
zero. Using Equation 14.13 and Equation 15.9 in Equation 15.1 results in
the gradient used for the gradient-based optimization method in section
7.

# Derivative Verification

This section provides a verification of the $f_{21}(B,\ L)$ derivative
using the finite-difference and complex step methods.

The finite-difference method formula is:

$$\begin{array}{r}
\frac{\partial f}{\partial x_{j}} = \ \frac{f\left( x + h{\widehat{e}}_{j} \right) - f(x)}{h} + \mathcal{O}(h)\#(16) \\
\end{array}$$

From Equation 3.1 and Equation 16 we obtain:

$$\begin{array}{r}
\frac{\partial f_{21}}{\partial B} = \frac{\left( \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} + h{\widehat{e}}_{B} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} + h{\widehat{e}}_{B} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} + h{\widehat{e}}_{B} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} + h{\widehat{e}}_{B} - \overrightarrow{s_{1}} \right| \right|}\  \right) - \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} - \overrightarrow{s_{1}} \right| \right|}\  \right) \right)}{h}(16.1) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{21}}{\partial L} = \frac{\left( \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} + h{\widehat{e}}_{L} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} + h{\widehat{e}}_{L} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} + h{\widehat{e}}_{L} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} + h{\widehat{e}}_{L} - \overrightarrow{s_{1}} \right| \right|}\  \right) - \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} - \overrightarrow{s_{1}} \right| \right|}\  \right) \right)}{h}(16.2) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{21}}{\partial H} = \frac{\left( \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} + h{\widehat{e}}_{H} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} + h{\widehat{e}}_{H} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} + h{\widehat{e}}_{H} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} + h{\widehat{e}}_{H} - \overrightarrow{s_{1}} \right| \right|}\  \right) - \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} - \overrightarrow{s_{1}} \right| \right|}\  \right) \right)}{h}(16.3) \\
\end{array}$$

The complex step formula is:

$$\begin{array}{r}
\frac{\partial f}{\partial x_{j}} = \ \frac{\mathfrak{I}\left( f\left( x + ih{\widehat{e}}_{j} \right) \right)}{h} + \mathcal{O}(h^{2})\ \#(17) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{21}}{\partial B} = \frac{\mathfrak{I}\left( \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} + ih{\widehat{e}}_{B} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} + ih{\widehat{e}}_{B} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} + ih{\widehat{e}}_{B} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} + ih{\widehat{e}}_{B} - \overrightarrow{s_{1}} \right| \right|}\  \right) \right)}{h}\ \#(17.1) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{21}}{\partial L} = \frac{\mathfrak{I}\left( \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} + ih{\widehat{e}}_{L} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} + ih{\widehat{e}}_{L} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} + ih{\widehat{e}}_{L} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} + ih{\widehat{e}}_{L} - \overrightarrow{s_{1}} \right| \right|}\  \right) \right)}{h}\ \#(17.2) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f_{21}}{\partial H} = \frac{\mathfrak{I}\left( \frac{1}{\lambda}\left( \frac{\left( \overrightarrow{u} + ih{\widehat{e}}_{H} - \overrightarrow{s_{2}} \right)^{T}\dot{s_{2}}}{\left| \left| \overrightarrow{u} + ih{\widehat{e}}_{H} - \overrightarrow{s_{2}} \right| \right|} - \frac{\left( \overrightarrow{u} + ih{\widehat{e}}_{H} - \overrightarrow{s_{1}} \right)^{T}\dot{s_{1}}}{\left| \left| \overrightarrow{u} + ih{\widehat{e}}_{H} - \overrightarrow{s_{1}} \right| \right|}\  \right) \right)}{h}\ \#(17.3) \\
\end{array}$$

Where ${\overrightarrow{e}}_{j}$ is a unit vector in the $j_{th}$
direction (B, L, H for longitude, latitude, and altitude respectively),
$h$ is a step size, and $\overrightarrow{u}$ is the point about which to
differentiate (in B, L, H). The approximation of the derivatives (Eq.
16.1, 16.2, 16.3, 17.1, 17.2, 17.3) and the analytical derivative in
Equation 13.12 is compared, with results shown in Table 2:

  -----------------------------------------------------------------------------------------------------------------
  **Method**     **Derivative**                           **Step size**    **Point (B, L,     **Result**
                                                                           H)**               
  -------------- ---------------------------------------- ---------------- ------------------ ---------------------
  Analytical     $$\frac{\partial f_{21}}{\partial B}$$   $$10^{- 8}$$     $$(5,\ 10,\ 0)$$   -0.57761319

  Finite         $$\frac{\partial f_{21}}{\partial B}$$   $$10^{- 8}$$     $$(5,\ 10,\ 0)$$   -0.5775973477284424
  Difference                                                                                  

  Complex Step   $$\frac{\partial f_{21}}{\partial B}$$   $$10^{- 200}$$   $$(5,\ 10,\ 0)$$   -0.5834981372799554

  Analytical     $$\frac{\partial f_{21}}{\partial L}$$   $$10^{- 8}$$     $$(5,\ 10,\ 0)$$   3.90670032

  Finite         $$\frac{\partial f_{21}}{\partial L}$$   $$10^{- 8}$$     $$(5,\ 10,\ 0)$$   3.9067515444912715
  Difference                                                                                  

  Complex Step   $$\frac{\partial f_{21}}{\partial L}$$   $$10^{- 200}$$   $$(5,\ 10,\ 0)$$   3.905921531938494

  Analytical     $$\frac{\partial f_{21}}{\partial H}$$   $$10^{- 8}$$     $$(5,\ 10,\ 0)$$   -3.58224626

  Finite         $$\frac{\partial f_{21}}{\partial H}$$   $$10^{- 8}$$     $$(5,\ 10,\ 0)$$   -3.5822438348986907
  Difference                                                                                  

  Complex Step   $$\frac{\partial f_{21}}{\partial H}$$   $$10^{- 200}$$   $$(5,\ 10,\ 0)$$   -3.5822462633136434
  -----------------------------------------------------------------------------------------------------------------

  : Table 2: Comparison of derivative approximations and analytical
  solutions.

It is interesting that the complex step method is very accurate on the
altitude derivative but less accurate on the latitude and longitude than
the finite difference method. The error resulting from both
approximations also appears to be larger than the $\mathcal{O}(h)$ and
$\mathcal{O(}h^{2})$ remainder terms implied by Equations 16 & 17. This
could be due to the computer's finite arithmetic or how functions are
being handled "behind the scenes" in the Python code. Still, the
accuracy of the approximations is enough to verify the
$\frac{\partial f_{21}}{\partial\overrightarrow{u}}$ derivatives. It
would follow that the other FDOA function derivatives are correct due to
the implementation being the same.

# Gradient Based algorithm

The Optimization method used to find the source location was the
Levenberg-Marquardt algorithm (LMA). This is a surrogate-based
optimization method that estimates the hessian. LMA was developed to
solve non-linear least-square equations. The LMA interpolates between
the Gauss-Newton algorithm (GNA) and gradient-based optimization. GNA is
very fast at converging to the minimum but is not very accurate once it
gets close. This is where the gradient-based algorithm works best. LMA
will get close to the solution with GNA and then switch to the steepest
decent method to get a more precise solution. A derivation of the update
formula (Eq. 18) can be found in Appendix A.

$$\begin{array}{r}
s = \left( J^{T}J - \mu D \right)^{- 1}J^{T}\left( y - f(x) \right)\ \#(18) \\
\end{array}$$

Where $J$ is the Jacobian of the frequency derivatives for each
satellite at the current position, $y$ is the evaluation of the
frequency at the source location using equation 4, $f(x)$ is the
evaluation of the frequency at the estimated position that the
satellites are measuring, $\mu$ is the control of the direction ranging
between GNA and steepest decent. We can improve the scaling of the
function by multiplying by $D$. Where $\ D\ $is:

$$\begin{array}{r}
D = diag\left( J^{T}J \right)\#(18.1) \\
\end{array}$$

The algorithm works by initializing the code with an initial position,
$\rho$, $\mu$, and initial position $x$. The Jacobian and the residual
will then be evaluated. The first step is calculated and added to the
current position, and the error is then found using Equation 18.2. This
is the error for the next step ($e_{s}$).

$$\begin{array}{r}
es = \left\| e_{s} \right\|^{2}\ \#(18.2) \\
\end{array}$$

The difference in the current error $e$ and the previous error $e_{s}$
is then found using Equation 18.3.

$$\begin{array}{r}
\mathrm{\Delta} = es - e\ \#(18.3) \\
\end{array}$$

If the difference is negative then the new position will be accepted and
e, the Jacobian and the residual will be updated. For every fifth
iteration, the $\mu$ value will be divided by $\rho$ (the damping
factor) making it act more like GNA. If the difference in error is not
negative, $\mu$ will be multiplied by $\rho$ to make the algorithm act
more like a steepest decent. This code will then loop and check if the
difference in errors is less than a tolerance of $10^{- 3}$. If the
error is less than the tolerance it will return the position, if not the
code will loop again. This process is shown in the flow diagram in
Figure 4.

![A diagram of a process Description automatically
generated](./image5.png){width="4.851161417322834in"
height="3.0967607174103238in"}

Figure 4: Flow chart of the algorithim.

Figure 5 shows the convergence of the LMA algorithm from multiple
starting positions. This shows that the algorithm converges to the same
position from many different start
positions.![](./image6.png){width="4.225807086614173in"
height="3.1693547681539807in"}

Figure 5: Plot of the contour and the path of the algorithm.

Delta (Eq. 18.3) is used as the stopping criteria for the algorithm.
Once delta is less than $10^{- 3}$ the algorithm stops looping and
returns the position. The LMA algorithm converges in 15- 30 iterations.
It will evaluate the Jacobian and the residual $k + 1$ times where k is
the number of iterations it takes to converge. The step calculation will
be evaluated once every iteration so will have the same number of
evaluations as iterations. The relation between the number of
evaluations, the difference in the current error and next step error
(delta) can be seen in Figure 6 below.

![](./image7.png){width="3.967742782152231in"
height="2.9758070866141733in"}

Figure 6: Plot showing the convergence of the algorithm.

# Monty Carlo Comparison

To determine the accuracy of the grid-search algorithm and LMA, a Monte
Carlo simulation is run to estimate the root-mean-squared-error. The
formula for the calculation is as follows:

$$\begin{array}{r}
RMSE = \sqrt{\frac{1}{L}\sum_{l - 1}^{L}\left. \ \left\| u_{l} - u^{o} \right.\  \right\|^{2}}\ \#(19) \\
\#\  \\
\end{array}$$

The equation uses the total number of Monte Carlo simulations, L, and
the summation of the 2-norm squared difference between the position
vector found by the grid search and the actual position vector of the
source. The number of Monte Carlo runs used was $10^{4}$. The way that
the algorithm was run is shown in the Monte Carlo function in Appendix
B. While testing the grid search, noise levels of
$\sigma_{f} = \ ($`<!-- -->`{=html}1, 2, 4, 8, 16) Hz were used. This
process gave the results shown in Table 3 below.

  -----------------------------------------------------------------------
  Noise level (Hz)        LMA (km)                Grid search (km)
  ----------------------- ----------------------- -----------------------
  1                       1.0941220938193204      0.6197978901842434

  2                       2.1879750170699395      1.0849477826620166

  4                       4.365063017360187       2.2123866587402907

  8                       8.725234012378314       4.227109425169162

  16                      17.44703333621689       8.175176490117455
  -----------------------------------------------------------------------

  : Table 3: Table of the Monty Carlo results.

The results from the RMSE were plotted in Figure 7. This shows that the
results are linear and that the LMA algorithm follows the same trend as
the grid search. The relationship between the grid search and the LMA
was that the error doubled what the grid search error was.

Figure 7 7: Plot showing the RMSE evaluation of the grid search
algorithm and the LMA algorithm in km.

#  Appendix A: Levenberg-Marquardt Update Derivation {#appendix-a-levenberg-marquardt-update-derivation .unnumbered}

Let $y$ be the true function value, $f(B,\ L)$ the next function
estimate and $f\left( B_{0},\ L_{0} \right)$ the current function
estimate.

Let $\overrightarrow{x}$ denote $(B,\ L)$ and $\overrightarrow{x_{0}}$
denote $(B_{0},L_{0})$.

$$\begin{array}{r}
\left( y - f(B,\ L) \right)^{T}\left( y - f(B,L) \right)\#(A1) \\
\end{array}$$

$$\begin{array}{r}
f(B,L) = f\left( B_{0},\ L_{0} \right) + J\left( B_{0},L_{0} \right)\begin{pmatrix}
B - B_{0} \\
L - L_{0} \\
\end{pmatrix}\#(A1.1) \\
\end{array}$$

$$\begin{array}{r}
\left( y - \ f\left( B_{0},\ L_{0} \right) + J\left( B_{0},L_{0} \right)\begin{pmatrix}
B - B_{0} \\
L - L_{0} \\
\end{pmatrix} \right)^{T}\left( y - \ f\left( B_{0},\ L_{0} \right) + J\left( B_{0},L_{0} \right)\begin{pmatrix}
B - B_{0} \\
L - L_{0} \\
\end{pmatrix} \right)\ \#(A1.2) \\
\end{array}$$

$$\begin{array}{r}
g = y - f\left( \overrightarrow{x_{0}} \right) + J\left( \overrightarrow{x_{0}} \right)^{T}\overrightarrow{x_{0}}\ \#(A1.3) \\
\end{array}$$

$$\begin{array}{r}
\left( g - J\left( \overrightarrow{x_{0}} \right)^{T}\overrightarrow{x} \right)^{T}\left( g - J\left( \overrightarrow{x_{0}} \right)^{T}\overrightarrow{x} \right)\#(A1.4) \\
\end{array}$$

$$\begin{array}{r}
Let\ u_{i} = g - J\left( \overrightarrow{x_{0}} \right)^{T}\overrightarrow{x}\ \#(A1.5) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial f}{\partial x} = \frac{\partial f}{\partial u}\ \ \frac{\partial u}{\partial x}\ \#(A1.6) \\
\end{array}$$

$$\begin{array}{r}
 \rightarrow \frac{\partial}{\partial x}\Sigma u_{i}^{2} \rightarrow 2u_{i}\frac{\partial u_{i}}{\partial x}\ \#(A1.7) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial u_{i}}{\partial x} = \  - J\left( \overrightarrow{x_{0}} \right) + \frac{\partial g}{\partial x}\ \#(A1.8) \\
\end{array}$$

$$\begin{array}{r}
\frac{\partial g}{\partial x} = 0\ \#(A1.9) \\
\end{array}$$

$$\begin{array}{r}
 \rightarrow \ \frac{\partial f}{\partial x} = 2\left( g - J\left( \overrightarrow{x_{0}} \right)^{T}\overrightarrow{x} \right)\left( - J\left( \overrightarrow{x_{0}} \right) \right)\#(A1.10) \\
\end{array}$$

$$\begin{array}{r}
 \rightarrow \  - 2g^{T}J\left( \overrightarrow{x_{0}} \right) + 2J\left( \overrightarrow{x_{0}} \right)^{T}\overrightarrow{x}J\left( \overrightarrow{x_{0}} \right)\#(A1.11) \\
\end{array}$$

$$\begin{array}{r}
 \rightarrow J\left( \overrightarrow{x_{0}} \right)^{T}\overrightarrow{x}J\left( \overrightarrow{x_{0}} \right) = g^{T}J\left( \overrightarrow{x_{0}} \right)\#(A1.12) \\
\end{array}$$

$$\begin{array}{r}
\overrightarrow{x} = \left( J\left( \overrightarrow{x_{0}} \right)^{T}J\left( \overrightarrow{x_{0}} \right) \right)^{- 1}g^{T}J\left( \overrightarrow{x_{0}} \right)\#(A1.13) \\
\end{array}$$

$$\begin{array}{r}
\begin{pmatrix}
B \\
L \\
\end{pmatrix} = \left( J\left( B_{0},L_{0} \right)^{T}J\left( B_{0},L_{0} \right) \right)^{- 1}g^{T}J\left( B_{0},L_{0} \right)\#(A1.14) \\
\end{array}$$

For the Levenberg-Marquardt Algorithm a factor of $\mu D$ is added,
where $\mu$ is a scaling constant.

$$\begin{array}{r}
D = diag\left( J\left( B_{0},L_{0} \right)^{T}J\left( B_{0},L_{0} \right) \right)\#(A1.15) \\
\end{array}$$

$$\begin{array}{r}
\begin{pmatrix}
B \\
L \\
\end{pmatrix} = {- \left( J\left( B_{0},L_{0} \right)^{T}J\left( B_{0},L_{0} \right) + \mu D \right)}^{- 1}g^{T}J\left( B_{0},L_{0} \right)\ \#(A1.16) \\
\end{array}$$

# Appendix B: Algorithm and Grid Search Code. {#appendix-b-algorithm-and-grid-search-code. .unnumbered}

import numpy as np

import matplotlib.pyplot as plt

import random

fs = 1\*10\*\*9      #frequency of transmitter

wl = (3\*10\*\*5)/fs #wave length

#source location

u0 = np.array(\[\[2\], \[1\], \[0\]\])

#Data matrix

theta = np.empty(shape = (501, 501))

#Weighting function

Q = np.array(\[\[2, 1, 1\],

            \[1, 2, 1\],

            \[1, 1, 2\]\])

Qinv = np.linalg.inv(Q)

#making the noise vector

def nf(lvl):

    n1 = random.gauss(0, lvl)

    n2 = random.gauss(0, lvl)

    n3 = random.gauss(0, lvl)

    n4 = random.gauss(0, lvl)

    nf = np.array(\[\[n2-n1\], \[n3-n1\], \[n4-n1\]\])

    return nf

#satalites

s1 = np.array(\[\[7378.1, 0, 0\]\]).T

s1d = np.array(\[\[0.0001, 4.4995, 5.3623\]\]).T

s2 = np.array(\[\[7377.5, 100, 0\]\]).T

s2d = np.array(\[\[-0.0671, 4.9493, 4.9497\]\]).T

s3 = np.array(\[\[7377.5, -100, 0\]\]).T

s3d = np.array(\[\[0.0610, 4.4991, 5.3623\]\]).T

s4 = np.array(\[\[7377.5, 0, 100\]\]).T

s4d = np.array(\[\[-0.0777, 4.0150, 5.7335\]\]).T

#Should convert from ECEF to LLA

def LLA_to_ECEF(lat, long):

    R0 = 6378.137

    e = 0.081819198425

    lat = np.radians(lat)

    long = np.radians(long)

    Re = R0/(np.sqrt(1-(e)\*\*2\*(np.sin(lat))\*\*2))

    x = (Re\*np.cos(lat)\*np.cos(long))

    y = (Re\*np.cos(lat)\*np.sin(long))

    z = ((1-e\*\*2)\*Re\*np.sin(lat))

    u = np.array(\[\[x\], \[y\], \[z\]\])

    return u

u0 = LLA_to_ECEF(u0\[0,0\], u0\[1,0\])\[:2\]

#grid search doppler shift

def dopler_shift(lat, long):

    u = LLA_to_ECEF(lat, long)

    f0 = 1/wl\*np.dot((u - s1).T , s1d)/(np.linalg.norm(u-s1))

    f1 = 1/wl\*np.dot((u - s2).T , s2d)/(np.linalg.norm(u-s2)) - f0

    f2 = 1/wl\*np.dot((u - s3).T , s3d)/(np.linalg.norm(u-s3)) - f0

    f3 = 1/wl\*np.dot((u - s4).T , s4d)/(np.linalg.norm(u-s4)) - f0

    f = np.array(\[f1\[0\], f2\[0\], f3\[0\]\])

    return f

#grid search calc

def grid_search():

    long = 0

    lat  = 0

    min = float(\'inf\')

    coord = 0

    noise = nf(2)  

    f = dopler_shift(2, 1)

    for i in range(0, 501):

        for j in range(0, 501):

            g = dopler_shift(lat, long)

            g += noise

            theta\[j\]\[i\] = np.log((f - g).T @ Qinv @ (f - g))

            if theta\[j\]\[i\] \< min: #grid search stuff

                min = theta\[j\]\[i\]

                coord = np.array(\[\[j\], \[i\], \[0\]\])\*0.01

            lat += 0.01

        lat = 0

        long += 0.01

    xyz = LLA_to_ECEF(coord\[0,0\], coord\[1,0\])

    return coord, xyz\[:2\]

#monte carlo simulation to find the accuracy of the data

def monte_carlo():

    sumation = 0

    L = 500

    for l in range(0, L):

        coord, xyz = grid_search()

       

        sumation += np.linalg.norm(xyz - u0)\*\*2

        print(l)

    RMSE = np.sqrt((1/L)\*sumation)

    print(RMSE)

#plotting stuff

def graph():

    coord = grid_search()

    tmin = np.min(theta)

    tmax = np.max(theta)

    levels = np.linspace(tmin, tmax, 100)

    feature_x = np.arange(0, 5.01, 0.01)

    feature_y = np.arange(0, 5.01, 0.01)

    \[X, Y\] = np.meshgrid(feature_x, feature_y)

    contour = plt.contourf(X, Y, theta, levels = levels)

    plt.colorbar(contour, label=\'log(theta-values)\')

    plt.xlabel(\"Longitude (Degrees)\")

    plt.ylabel(\"Latitude (Degrees)\")

    plt.title(\"FDOA plot\")

    #print(\"Minmum is found at: \", coord.T)

    #plt.annotate(\'Min\',xy=(10,5),xytext=(5,10),arrowprops={})

   

    plt.show()

def main():

    monte_carlo()

    #graph()

main()

This code is also in our Git hub in project.py.

https://github.com/Mangin4/ENEL445_Optimization
