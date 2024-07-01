




ENEL445: Engineering Optimization

Source Geo-Location

Final Report

Submitted by: 

Daegan Rose-Love (33480267)

Ben Mangin (92699294)

Date: 31/05/2024


# Table of Contents
[1.	Background	3](#_toc168070731)

[1.1.	Brief	3](#_toc168070732)

[1.2.	Math Background	3](#_toc168070733)

[1.3 Assumptions	4](#_toc168070734)

[2.	FDOA-Only Objective Function Formulation	5](#_toc168070735)

[3.	Grid Search	6](#_toc168070736)

[4.	Find FDOA With Unknown Frequency	8](#_toc168070737)

[5.	Gradient Derivation	10](#_toc168070738)

[6.	Derivative Verification	12](#_toc168070739)

[7.	Gradient Based algorithm	14](#_toc168070740)

[8.	Monty Carlo Comparison	17](#_toc168070741)

[Appendix A: Levenberg-Marquardt Update Derivation	19](#_toc168070742)

[Appendix B: Algorithm and Grid Search Code.	21](#_toc168070743)




1. # <a name="_toc168070731"></a>Background
   1. ## <a name="_toc168070732"></a>Brief
Radio source geo-location involves locating the coordinates of where a particular signal is broadcast by using arbitrary measurements. This project uses four satellites to locate a ground-based signal using the frequency-difference of arrival (FDOA) measurements. Source geo-location has a range of applications, from military to civil, including surveillance, navigation, and search and rescue.
1. ## <a name="_toc168070733"></a>Math Background
When calculating the Doppler shift of the signal the coordinates need to be in the Earth-centred Earth Fixed (ECEF) coordinate system as shown in vector (1) as x, y, and z. The position and velocity of the satellites are already in the ECEF system which is shown in vectors (1.2) and (1.3) but the grid points in the search area (grid points) are geodetic positions as shown in vector (1.1) where B is latitude, L is longitude, and H is altitude. 

u0= x0y0z0T #1

p0=B0L0H0T #1.1

si= sx,isy,isz,iT #1.2

si=sx,isy,isz,iT #1.3

To convert the coordinates the geodetic coordinates, need to be converted into radians. Once the coordinates were converted to radians the equatorial radius is found as shown in equation (2). REB0= R0(1-e2sin2B0 #2REB0= R0(1-e2sin2B0 #1

This uses the Earth's equatorial radius R0=6378.137 km and the Earth's eccentricity e=0.081819198425. The resultant can then be used in the following equations which will give the converted x, y, and z coordinates.

x0=REB0+HcosB0cosL0#2.1

y0=REB0+HcosB0sinL0#2.2

z0=1+e2REB0+HsinB0#2.3

The equation shown below in (3) calculates the frequency of the signal received by the satellite. This is not the same frequency as the transmitted frequency. The transmitted signal will be subjected to Doppler shift and channel noise. fi= f0+1λu0-sisiu0-si+ni #3

The equation uses the frequency of the source and the position vector of the point that is being measured, the position vector of the satellite, the velocity vector of the satellite, and additive noise. As the frequency of the source is unknown this equation can’t be used. To get around this the frequency of the received signal can be compared to the frequency of the first satellite as shown in the following equation (3.1). This will result in a vector of frequencies as shown below in vector (3.2).

f=fj-f1 = 1λu0-sjsju0-sj-1λu0-s1s1u0-s1+(nj-n1) , for j=2,3,4 #3.1

→f= f21f31f41T #3.2

The following vector (3.3) is the noise vector used to add noise to each frequency vector at each coordinate. Each n value is a random value chosen within a Gaussian function.

→n=n2-n1n3-n1n4-n1 #3.3

When creating a contour plot of the objective function the frequency at each grid point needs to be evaluated. This will be done as shown in the vector below (4). Where each value is the Doppler shift a satellite receives minus the Doppler shift of satellite one.

gB0, L0= gB0, L0, s2,s2-gB0,L0,s1,s1gB0, L0, s3,s3-gB0,L0,s1,s1gB0, L0, s4,s4-gB0,L0,s1,s1 #4

The following equation (5) allows us to find the Q matrix for our weighting function used to find theta.

Covxi, xj=Exi-µixj-µj=EΔfΔfT#5

The function for a multivariate Gaussian probability function (PDF) is shown below. 

Nx;μ, Σ=12πΣexp-x-μTΣ-1x-μ2 #6



X is the known measurement vector, μ is the mean vector, and the Σ is the covariance matrix composed of σ2Q.
## <a name="_toc168070734"></a>1.3 Assumptions 
The assumptions that were made for this simulation are as follows:

- Altitude is zero and does not change.
- Measurements do not account for special relativity.
- The source frequency is known.
- The location of the source is known for the grid-search and Monte-Carlo simulation.
1. # <a name="_toc168070735"></a>FDOA-Only Objective Function Formulation
For convenience, set the noise vector found in (3.3) to be:

Δf= n2-n1n3-n1n4-n1 #7

Then by using (5) the following is found:

En2-n1n3-n1n4-n1 n2-n1n3-n1n4-n1 #7.1

After expanding and removing cross terms we are left with the uncorrelated noise (7.2).

→En22+n12n12n12n12n32+n12n12n12n12n42+n12 #7.2

Notice that each noise term is ∝ σf2 resulting in the following:

Σ=2σf2σf2σf2σf22σf2σf2σf2σf22σf2= σf2211121112= σf2Q #7.3 

Q=211121112 #7.4

To solve this, we are interested in finding unknown parameters on which our measurements depend. By using maximum likelihood estimation (MLE) we locate the values that assign the highest probability to the measurements. The MLE is defined as follows:

θ=maxθpDθ #(8) 

Where pDθ is the probability distribution function of 𝒟 given θ. When the measurements in D=[y1, y2,…, yn] are identically distributed and independent, the objective function becomes (8.1). Taking the log of both sides gives: 

log(pDθ)= logΣi=1mpyiθ #8.1

In the case where G is a measurement matrix and ϵ ∼ 𝒩(ϵ; 0, Σ), measurements y can be written as:

y=y1, y2, …, ymT=Gθ+ϵ  #(8.2)

The measurements thus follow a Gaussian distribution Ny;Gθ, Σ as shown in (6).

pyθ=Ny;Gθ, Σ#8.3

The log of the multivariate Gaussian distribution removes the exponential from the equation giving the following equation.

logpyθ∝-y-GθTΣ-1y-Gθ#8.4

To find the signal source the maximum needs to be found. The equation can be further simplified by removing the negative sign and finding the minimum. This is shown in the equation below.

θ=maxθ(-logpyθ)=minθy-gθTΣ-1y-gθ #8.5

The measurement vector in this case is f , which gives.

→θ=minθf-gθTΣ-1f-gθ #9→θ=minθf-gθTΣ-1f-gθ #8.6

The objective function in this form is in fact a non-linear least squares (NLS) problem.
1. # <a name="_toc168070736"></a>Grid Search 
To develop a map of the objective function (the nonlinear least squares problem) the Doppler shift of the signal received from each satellite needed to be found. This was done using Equation (3). As the frequency of the source is not known the difference in Doppler shift between satellites 1 and 2, 3, and 4 uses Equation 3.1, this also allows for the generation of the source vector at (5, 10, 0) as shown in Vector 3.2. To get the measured frequency vector for each grid point the coordinates of the 40 by 40 grid were iterated through and solved as shown in Vector 4. At each coordinate the expected frequency vector and the measured frequency vector were substituted into Equation (6) to give. This would then be recorded into a matrix that was then plotted as a contour map (using Python) to get a visual representation of the objective function. The code that calculates for each plot is shown in the Python code in Appendix A. The contour plot is shown in Figure 1 below.

![](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.001.png)

<a name="_ref165032741"></a>*Figure 1: Contour plot of the signal with no noise.*

Figure 1 shows that there is a minimum at g(5,10) which is the location of the source. Although this contour graph does not have noise, the following contour map in [Figure 2](#mergeformat) does. The noise added to the simulation for Figure 2 was 1Hz which resulted in the minimum being found 600 meters away from the true minimum.

![A chart with a gradient of green and blue

Description automatically generated](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.002.png)

<a name="_ref165033056"></a>*Figure 2: Contour plot of the signal with 1Hz of noise.*

To find the minimum within the data a grid search algorithm was implemented. This is done by iterating over each grid point and checking if the new grid point is smaller. Once every grid point is checked, the coordinates at that minimum are returned. This algorithm was implemented in the function for calculating as it already steps through each coordinate with a step size of 0.01. This allows a 0.01-degree accuracy/resolution. The algorithm is shown in the Python code in Appendix B. The grid search algorithm was tested with different noise levels. The results of this test are shown in Table 1. 

<a name="_ref165045289"></a><a name="_ref165045263"></a>*Table 1: Minima found with different noise levels.*

|**NOSIE DEVIATION (HZ)**|**MINIMA (B0, L0,H0)**|
| :- | :- |
|0|(5, 10, 0)|
|1|(5, 10, 0)|
|2|(5, 10, 0)|
|4|(5, 10.1, 0)|
|8|(4.9, 9.9, 0)|
|16|(5.1, 10.1, 0)|
1. # <a name="_toc168070737"></a>Find FDOA With Unknown Frequency
In a real-world situation, the original frequency of the transmitter is unknown which means that the wavelength of the signal cannot be found. To get around this an estimation of the original frequency needs to be made to get an estimate of the wavelength. This can be obtained by evaluating the frequency at satellite j using the following equation.

fj=1λu0-siTsiu0-si+f0,  for j=1,2,3,4#(10)

The frequencies at each satellite can then be averaged to get an estimate of the frequency as shown in equation 11. 

fe=f1+f2+f3+f44#(11)

Using this estimate of the frequency the wavelength can be found as shown below where c is the speed of light and λ is the wavelength.

1λ=cfe#(12)

This can be rearranged to give lambda as shown in equation 12.1.

λ=fec#(12.1)

This can then be substituted back into equation 3.1 to give the frequency difference between satellites as shown in equation 13.

fj1=fec u0-sjTsju0-sj-fec u0-s1Ts1u0-s1+nj-n1,  for j=2,3,4 #(13)

This now allows for the frequency difference between satellites to be found so the objective function can be formulated. Results for this frequency estimation can be seen in Figure 3. This was also verified using the grid search algorithm to check the minimum was at (5,10) which it was.

![A chart with a gradient of green and blue colors

Description automatically generated](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.003.png)![](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.004.png)

*Figure 3: Two plots showing the difference between plots of the know frequency and the estimated frequency.*

**Note:** Using more FDOA measurements such as f23, f13,f43 does not result in any performance improvements. The additional measurements are the result of values being shifted but the data in the matrix retains the same information with a different reference satellite, this can be shown mathematically:

-1100-1010-1001f1f2f3f4=f21f31f41= 01-1010-1000-11f1f2f3f4=f23f13f43 #13.1

1. # <a name="_toc168070738"></a>Gradient Derivation
To apply gradient-based optimization methods, the gradient must be computed or approximated. From Equation 3.1 we have:

fj1=1λ  u0-sjTsju0-sj-1λ  u0-s1Ts1u0-s1 , ∂fj1∂u0=0 #14

` `Let 

a=u0-sj #14.1

a1=u0-s1 #14.2

fj1=1λ  aTsja-1λ  a1Ts1a1 , ∂fj1∂u0,   ∂fj1∂u0=∂fj1∂a∂a∂u0-∂fj1∂a∂a∂u0 #14.3

Let

a=b12= ∑ai2 #14.4

fj1=1λaTsjb12-1λa1Ts1b12 #14.5

∂fj1∂a= ∂∂a1a1Tsjb-12λ #14.6

Using product rule and chain rule gives:

∂fj1∂u0=1λsjb-12-aTsj2b32∂b∂a∂a∂u0 #14.7

b122=∑ai22=∑ai2=b,  ∂b∂a=∑2ai #14.8

∂a∂u0=∂∂u0u0-sj=111 #14.9

∂fj1∂u0=1 λsjb12-aTsj 2b322a111-1 λs1b12-aTs1 2b322a1111 #14.10

∂fj1∂u0=1λsja-aTsj2a32a-1λs1a1-a1Ts12a132a1 #14.11

Substituting a back in gives the derivative:

∂fj1∂u0=  1λsju0-sj-u0-sjTsju0-sj3(u0-sj)-1λs1u0-s1-u0-s1Ts1u0-s13u0-s1 #14.12

The resulting Jacobian can be written in the following form:

Jfu=∂f21∂x∂f21∂y∂f21∂z∂f31∂x∂f31∂y∂f31∂z∂f41∂x∂f41∂y∂f41∂z= ∇f21uT∇f31uT∇f41uT #14.13

The derivatives of the coordinate conversion functions (Eq. 2, 2.1, 2.2, 2.3, 2.4) are used to find the derivative of fj1 with respect to p (Eq. 1.1). The resulting Jacobian Jf(p) is shown below, with derived terms.

Jfp= ∂x∂B∂x∂L∂x∂H∂y∂B∂y∂L∂y∂H∂z∂B∂z∂L∂z∂H #15

The desired Jacobian for ∂f∂p is:

Jfu, p=JfuJfp#15.1

The derivation results follow:

∂Re∂B=R0e2cosBsinB(1-e2sin2(B))23 #15.2

∂x0∂B0=cosL0∂Re∂BcosB0-ReB0sinB0 #15.3

∂x0∂L0=ReB0cosB0sinL0 #15.4

∂y0∂B0=sinL0∂Re∂BcosB0-ReB0sinB0 #15.5

∂y0∂L0=ReB0cosB0cosL0 #15.6

∂z0∂B0=1-e2∂Re∂BsinB0+ReB0cosB #15.7

∂z0∂L0=0 #15.8

Jfp=∂x0∂B0∂x0∂L0∂y0∂B0∂y0∂L0∂z0∂B0∂z0∂L0#15.9

Note that the H column is omitted as the altitude in degrees (H) is zero. Using Equation 14.13 and Equation 15.9 in Equation 15.1 results in the gradient used for the gradient-based optimization method in section 7.
1. # <a name="_toc168070739"></a>Derivative Verification
This section provides a verification of the f21(B, L) derivative using the finite-difference and complex step methods.

The finite-difference method formula is:

∂f∂xj= fx+hej-fxh+Oh#16

From Equation 3.1 and Equation 16 we obtain:

∂f21∂B=1λu+heB-s2Ts2u+heB-s2-u+heB-s1Ts1u+heB-s1 -1λu-s2Ts2u-s2-u-s1Ts1u-s1 h16.1

∂f21∂L=1λu+heL-s2Ts2u+heL-s2-u+heL-s1Ts1u+heL-s1 -1λu-s2Ts2u-s2-u-s1Ts1u-s1 h16.2

∂f21∂H=1λu+heH-s2Ts2u+heH-s2-u+heH-s1Ts1u+heH-s1 -1λu-s2Ts2u-s2-u-s1Ts1u-s1 h16.3

The complex step formula is:

∂f∂xj= Ifx+ihejh+O(h2) #17

∂f21∂B=I1λu+iheB-s2Ts2u+iheB-s2-u+iheB-s1Ts1u+iheB-s1 h #17.1

∂f21∂L=I1λu+iheL-s2Ts2u+iheL-s2-u+iheL-s1Ts1u+iheL-s1 h #17.2

∂f21∂H=I1λu+iheH-s2Ts2u+iheH-s2-u+iheH-s1Ts1u+iheH-s1 h #17.3

Where ej is a unit vector in the jth direction (B, L, H for longitude, latitude, and altitude respectively), h is a step size, and u is the point about which to differentiate (in B, L, H).  The approximation of the derivatives (Eq. 16.1, 16.2, 16.3, 17.1, 17.2, 17.3) and the analytical derivative in Equation 13.12 is compared, with results shown in [Table 2](#mergeformat):



<a name="_ref167998168"></a><a name="_ref167998157"></a>*Table 2: Comparison of derivative approximations and analytical solutions.*

|**Method**|**Derivative**|**Step size** |**Point (B, L, H)**|**Result**|
| :- | :- | :- | :- | :- |
|Analytical|∂f21∂B|10-8|(5, 10, 0)|-0.57761319|
|Finite Difference|∂f21∂B|10-8|(5, 10, 0)|-0.5775973477284424|
|Complex Step|∂f21∂B|10-200|(5, 10, 0)|-0.5834981372799554|
|Analytical|∂f21∂L|10-8|(5, 10, 0)|3\.90670032|
|Finite Difference|∂f21∂L|10-8|(5, 10, 0)|3\.9067515444912715|
|Complex Step|∂f21∂L|10-200|(5, 10, 0)|3\.905921531938494|
|Analytical|∂f21∂H|10-8|(5, 10, 0)|-3.58224626|
|Finite Difference|∂f21∂H|10-8|(5, 10, 0)|-3.5822438348986907|
|Complex Step|∂f21∂H|10-200|(5, 10, 0)|-3.5822462633136434|

It is interesting that the complex step method is very accurate on the altitude derivative but less accurate on the latitude and longitude than the finite difference method. The error resulting from both approximations also appears to be larger than the Oh and O(h2) remainder terms implied by Equations 16 & 17. This could be due to the computer’s finite arithmetic or how functions are being handled “behind the scenes” in the Python code. Still, the accuracy of the approximations is enough to verify the ∂f21∂u derivatives. It would follow that the other FDOA function derivatives are correct due to the implementation being the same. 
1. # <a name="_toc168070740"></a>Gradient Based algorithm
The Optimization method used to find the source location was the Levenberg-Marquardt algorithm (LMA). This is a surrogate-based optimization method that estimates the hessian. LMA was developed to solve non-linear least-square equations. The LMA interpolates between the Gauss-Newton algorithm (GNA) and gradient-based optimization. GNA is very fast at converging to the minimum but is not very accurate once it gets close. This is where the gradient-based algorithm works best. LMA will get close to the solution with GNA and then switch to the steepest decent method to get a more precise solution. A derivation of the update formula (Eq. 18) can be found in Appendix A. 

s=JTJ-μD-1JTy-fx #18

Where J is the Jacobian of the frequency derivatives for each satellite at the current position, y is the evaluation of the frequency at the source location using equation 4, f(x) is the evaluation of the frequency at the estimated position that the satellites are measuring, μ is the control of the direction ranging between GNA and steepest decent. We can improve the scaling of the function by multiplying by D. Where  D is:

D=diagJTJ#18.1

The algorithm works by initializing the code with an initial position, ρ, μ, and initial position x. The Jacobian and the residual will then be evaluated. The first step is calculated and added to the current position, and the error is then found using Equation 18.2. This is the error for the next step (es).

es=es2 #18.2

The difference in the current error e and the previous error es is then found using Equation 18.3.

∆=es-e #18.3

If the difference is negative then the new position will be accepted and e, the Jacobian and the residual will be updated. For every fifth iteration, the μ value will be divided by ρ (the damping factor) making it act more like GNA. If the difference in error is not negative, μ will be multiplied by ρ to make the algorithm act more like a steepest decent. This code will then loop and check if the difference in errors is less than a tolerance of 10-3. If the error is less than the tolerance it will return the position, if not the code will loop again. This process is shown in the flow diagram in Figure 4.

![A diagram of a process

Description automatically generated](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.005.png)

*Figure 4: Flow chart of the algorithim.*

Figure 5 shows the convergence of the LMA algorithm from multiple starting positions. This shows that the algorithm converges to the same position from many different start positions.![](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.006.png)

*Figure 5: Plot of the contour and the path of the algorithm.*

Delta (Eq. 18.3) is used as the stopping criteria for the algorithm. Once delta is less than 10-3 the algorithm stops looping and returns the position. The LMA algorithm converges in 15- 30 iterations. It will evaluate the Jacobian and the residual k+1 times where k is the number of iterations it takes to converge. The step calculation will be evaluated once every iteration so will have the same number of evaluations as iterations. The relation between the number of evaluations, the difference in the current error and next step error (delta) can be seen in Figure 6 below.

![](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.007.png)

*Figure 6: Plot showing the convergence of the algorithm.*

1. # <a name="_toc168070741"></a>Monty Carlo Comparison
To determine the accuracy of the grid-search algorithm and LMA, a Monte Carlo simulation is run to estimate the root-mean-squared-error. The formula for the calculation is as follows:

RMSE=1Ll-1Lul-uo2 #(19)# 

The equation uses the total number of Monte Carlo simulations, L, and the summation of the 2-norm squared difference between the position vector found by the grid search and the actual position vector of the source. The number of Monte Carlo runs used was 104.  The way that the algorithm was run is shown in the Monte Carlo function in Appendix B. While testing the grid search, noise levels of σf= (1, 2, 4, 8, 16) Hz were used. This process gave the results shown in Table 3 below. 

*Table 3: Table of the Monty Carlo results.*

|**NOISE LEVEL (HZ)**|**LMA (KM)**|**GRID SEARCH (KM)**|
| :- | :- | :- |
|**1**|1\.0941220938193204|0\.6197978901842434|
|**2**|2\.1879750170699395|1\.0849477826620166|
|**4**|4\.365063017360187|2\.2123866587402907|
|**8**|8\.725234012378314|4\.227109425169162|
|**16**|17\.44703333621689|8\.175176490117455|

The results from the RMSE were plotted in Figure 7. This shows that the results are linear and that the LMA algorithm follows the same trend as the grid search. The relationship between the grid search and the LMA was that the error doubled what the grid search error was.

![](Aspose.Words.bc329e17-135e-4133-9302-76a16f0a93f5.008.png)

*Figure 7 7: Plot showing the RMSE evaluation of the grid search algorithm and the LMA algorithm in km.*
# <a name="_toc168070742"></a>Appendix A: Levenberg-Marquardt Update Derivation
Let y be the true function value, fB, L the next function estimate and fB0, L0 the current function estimate. 

Let x denote (B, L) and x0 denote (B0,L0).

y-fB, LTy-fB,L#A1

fB,L=fB0, L0+JB0,L0B-B0L-L0#A1.1

y- fB0, L0+JB0,L0B-B0L-L0Ty- fB0, L0+JB0,L0B-B0L-L0 #A1.2

g=y-fx0+Jx0Tx0 #A1.3

g-Jx0TxTg-Jx0Tx#A1.4

Let ui=g-Jx0Tx #A1.5

∂f∂x=∂f∂u  ∂u∂x #A1.6

→∂∂xΣui2→2ui∂ui∂x #A1.7

∂ui∂x= -Jx0+∂g∂x #A1.8

∂g∂x=0 #A1.9

→ ∂f∂x=2g-Jx0Tx-Jx0#A1.10

→ -2gTJx0+2Jx0TxJx0#A1.11

→Jx0TxJx0=gTJx0#A1.12

x=Jx0TJx0-1gTJx0#A1.13

BL=JB0,L0TJB0,L0-1gTJB0,L0#A1.14

For the Levenberg-Marquardt Algorithm a factor of μD is added, where μ is a scaling constant.

D=diagJB0,L0TJB0,L0#A1.15

BL=-JB0,L0TJB0,L0+μD-1gTJB0,L0 #A1.16


# <a name="_toc168070743"></a>Appendix B: Algorithm and Grid Search Code.
import numpy as np

import matplotlib.pyplot as plt

import random

fs = 1\*10\*\*9      #frequency of transmitter

wl = (3\*10\*\*5)/fs #wave length 

#source location

u0 = np.array([[2], [1], [0]])

#Data matrix

theta = np.empty(shape = (501, 501)) 

#Weighting function

Q = np.array([[2, 1, 1],

`            `[1, 2, 1],

`            `[1, 1, 2]])

Qinv = np.linalg.inv(Q)

#making the noise vector

def nf(lvl):

`    `n1 = random.gauss(0, lvl)

`    `n2 = random.gauss(0, lvl)

`    `n3 = random.gauss(0, lvl)

`    `n4 = random.gauss(0, lvl)

`    `nf = np.array([[n2-n1], [n3-n1], [n4-n1]])

`    `return nf

#satalites

s1 = np.array([[7378.1, 0, 0]]).T

s1d = np.array([[0.0001, 4.4995, 5.3623]]).T

s2 = np.array([[7377.5, 100, 0]]).T

s2d = np.array([[-0.0671, 4.9493, 4.9497]]).T

s3 = np.array([[7377.5, -100, 0]]).T

s3d = np.array([[0.0610, 4.4991, 5.3623]]).T

s4 = np.array([[7377.5, 0, 100]]).T

s4d = np.array([[-0.0777, 4.0150, 5.7335]]).T

#Should convert from ECEF to LLA

def LLA\_to\_ECEF(lat, long):

`    `R0 = 6378.137

`    `e = 0.081819198425

`    `lat = np.radians(lat)

`    `long = np.radians(long)

`    `Re = R0/(np.sqrt(1-(e)\*\*2\*(np.sin(lat))\*\*2))

`    `x = (Re\*np.cos(lat)\*np.cos(long))

`    `y = (Re\*np.cos(lat)\*np.sin(long))

`    `z = ((1-e\*\*2)\*Re\*np.sin(lat))

`    `u = np.array([[x], [y], [z]])

`    `return u

u0 = LLA\_to\_ECEF(u0[0,0], u0[1,0])[:2]

#grid search doppler shift

def dopler\_shift(lat, long):

`    `u = LLA\_to\_ECEF(lat, long)

`    `f0 = 1/wl\*np.dot((u - s1).T , s1d)/(np.linalg.norm(u-s1))

`    `f1 = 1/wl\*np.dot((u - s2).T , s2d)/(np.linalg.norm(u-s2)) - f0

`    `f2 = 1/wl\*np.dot((u - s3).T , s3d)/(np.linalg.norm(u-s3)) - f0

`    `f3 = 1/wl\*np.dot((u - s4).T , s4d)/(np.linalg.norm(u-s4)) - f0 

`    `f = np.array([f1[0], f2[0], f3[0]])

`    `return f

#grid search calc

def grid\_search():

`    `long = 0

`    `lat  = 0

`    `min = float('inf')

`    `coord = 0

`    `noise = nf(2)  

`    `f = dopler\_shift(2, 1)

`    `for i in range(0, 501):

`        `for j in range(0, 501):

`            `g = dopler\_shift(lat, long)

`            `g += noise

`            `theta[j][i] = np.log((f - g).T @ Qinv @ (f - g))

`            `if theta[j][i] < min: #grid search stuff

`                `min = theta[j][i]

`                `coord = np.array([[j], [i], [0]])\*0.01

`            `lat += 0.01

`        `lat = 0

`        `long += 0.01

`    `xyz = LLA\_to\_ECEF(coord[0,0], coord[1,0])

`    `return coord, xyz[:2]

#monte carlo simulation to find the accuracy of the data

def monte\_carlo():

`    `sumation = 0

`    `L = 500

`    `for l in range(0, L):

`        `coord, xyz = grid\_search()



`        `sumation += np.linalg.norm(xyz - u0)\*\*2

`        `print(l)

`    `RMSE = np.sqrt((1/L)\*sumation)

`    `print(RMSE)

#plotting stuff

def graph():

`    `coord = grid\_search()

`    `tmin = np.min(theta)

`    `tmax = np.max(theta)

`    `levels = np.linspace(tmin, tmax, 100)

`    `feature\_x = np.arange(0, 5.01, 0.01) 

`    `feature\_y = np.arange(0, 5.01, 0.01) 

`    `[X, Y] = np.meshgrid(feature\_x, feature\_y)

`    `contour = plt.contourf(X, Y, theta, levels = levels) 

`    `plt.colorbar(contour, label='log(theta-values)')

`    `plt.xlabel("Longitude (Degrees)")

`    `plt.ylabel("Latitude (Degrees)")

`    `plt.title("FDOA plot")

`    `#print("Minmum is found at: ", coord.T)

`    `#plt.annotate('Min',xy=(10,5),xytext=(5,10),arrowprops={})



`    `plt.show()

def main():

`    `monte\_carlo()

`    `#graph()

main()

This code is also in our Git hub in project.py.

https://github.com/Mangin4/ENEL445\_Optimization


