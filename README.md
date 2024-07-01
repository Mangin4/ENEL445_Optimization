ENEL445: Engineering Optimization 

Source Geo-Location Final Report 

Submitted by:               Daegan Rose-Love (33480267) Ben Mangin (92699294) 

Date: 31/05/2024 

Table of Contents 

1. [Background 	 3 ](#_page2_x69.00_y72.00)
1. [Brief 	 3 ](#_page2_x69.00_y103.00)
1. [Math Background 	 3 ](#_page2_x69.00_y221.00)

   3. [Assumptions 	 4 ](#_page3_x69.00_y652.00)
2. [FDOA-Only Objective Function Formulation 	 5 ](#_page4_x69.00_y131.00)
3. [Grid Search 	 6 ](#_page5_x69.00_y435.00)
4. [Find FDOA With Unknown Frequency 	 8 ](#_page7_x69.00_y311.00)
5. [Gradient Derivation 	 10 ](#_page9_x69.00_y72.00)
6. [Derivative Verification 	 12 ](#_page11_x69.00_y425.00)
7. [Gradient Based algorithm 	 14 ](#_page13_x69.00_y473.00)
8. [Monty Carlo Comparison 	 17 ](#_page16_x69.00_y97.00)

[Appendix A: Levenberg-Marquardt Update Derivation 	 19 ](#_page18_x69.00_y72.00)

[Appendix B: Algorithm and Grid Search Code.	 21 ](#_page20_x69.00_y72.00)

1. Background<a name="_page2_x69.00_y72.00"></a> 
1. Brief<a name="_page2_x69.00_y103.00"></a> 

Radio source geo-location involves locating the coordinates of where a particular signal is broadcast by using arbitrary measurements. This project uses four satellites to locate a ground- based signal using the frequency-difference of arrival (FDOA) measurements. Source geo- location has a range of applications, from military to civil, including surveillance, navigation, and search and rescue. 

2. Math<a name="_page2_x69.00_y221.00"></a> Background 

When calculating the Doppler shift of the signal the coordinates need to be in the Earth-centred Earth Fixed (ECEF) coordinate system as shown in vector (1) as x, y, and z. The position and velocity of the satellites are already in the ECEF system which is shown in vectors (1.2) and (1.3) but the grid points in the search area (grid points) are geodetic positions as shown in vector (1.1) where B is latitude, L is longitude, and H is altitude.  

⃗ 0 = [ 0 0 0]  (1) 0 = [ 0 0 0]  (1.1) 

- [ , , , ]  (1.2) 
- [ , , , ]  (1.3) 

To convert the coordinates the geodetic coordinates, need to be converted into radians. Once the coordinates were converted to radians the equatorial radius is found as shown in equation (2). 

( 0) =  √(1− 2 s0in2( 0) (2)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.001.png)

This uses the Earth's equatorial radius  0 = 6378.137  and the Earth's eccentricity  = 0.081819198425. The resultant can then be used in the following equations which will give the converted x, y, and z coordinates. 

0 = ( ( 0) + )cos( 0)cos( 0) (2.1) 0 = ( ( 0) + )cos( 0)sin( 0) (2.2) 0 = [(1 + 2) ( 0) + ]sin( 0) (2.3) 

The equation shown below in (3) calculates the frequency of the signal received by the satellite. This is not the same frequency as the transmitted frequency. The transmitted signal will be subjected to Doppler shift and channel noise. 

- + 1 (⃗ 0−  )  + (3)

  0 ‖⃗ 0−  ‖

[^1]The equation uses the frequency of the source and the position vector of the point that is being measured, the position vector of the satellite, the velocity vector of the satellite, and additive noise. As the frequency of the source is unknown this equation can’t be used. To get around this the frequency of the received signal can be compared to the frequency of the first satellite as shown in the following equation (3.1). This will result in a vector of frequencies as shown below in vector (3.2). 

- − = 1 (⃗ 0 −  )  − 1 (‖⃗ ⃗00−− [^2] 1)‖1 + ( − ) ,  = 2,3,4  (3.1) 1 ‖⃗ 0 −  ‖ 1
- = [ 21 31 41]  (3.2) 

The following vector (3.3) is the noise vector used to add noise to each frequency vector at each coordinate. Each n value is a random value chosen within a Gaussian function. 

2 − 1

- ⃗  = [ 3 − 1]  (3.3) 4 − 1

When creating a contour plot of the objective function the frequency at each grid point needs to be evaluated. This will be done as shown in the vector below (4). Where each value is the Doppler shift a satellite receives minus the Doppler shift of satellite one. 

( 0, 0, 2, 2) − ( 0, 0, 1, 1)

` `( 0, 0) = [ ( 0, 0, 3, 3) − ( 0, 0, 1, 1)]  (4) 

( 0, 0, 4, 4) − ( 0, 0, 1, 1)

The following equation (5) allows us to find the Q matrix for our weighting function used to find theta. 

( , ) = [( − µ )( − µ )] = [   ] (5) The function for a multivariate Gaussian probability function (PDF) is shown below.  

1 ( − ) Σ−1( − )

( ; ,Σ) =![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.002.png)exp(− )  (6) 

√2 Σ 2

- Measurements do not account for special relativity. 
- The source frequency is known. 
- The location of the source is known for the grid-search and Monte-Carlo simulation. 
2. FDOA-Only<a name="_page4_x69.00_y131.00"></a> Objective Function Formulation 

For convenience, set the noise vector found in (3.3) to be: 

2 − 1

Δ  = [ 3 − 1]  (7) 

4 − 1

Then by using (5) the following is found: 

2 − 1

[[ 3 − 1 ] [ 2 − 1 3 − 1 4 − 1]]  (7.1) 

4 − 1

After expanding and removing cross terms we are left with the uncorrelated noise (7.2). 

( 22 + 12) 12 12

- [ 12 ( 32 + 12) 12 ]  (7.2) 

12 12 ( 42 + 12)

Notice that each noise term is ∝  2 resulting in the following: 

2 2 2 2 1 1

2

Σ = ( 2 2 2 2 ) =  2 (1 2 1) =  2 (7.3)  

2 2 2 2 1 1 2

2 1 1

- (1 2 1)  (7.4) 1 1 2

To solve this, we are interested in finding unknown parameters on which our measurements depend. By using maximum likelihood estimation (MLE) we locate the values that assign the highest probability to the measurements. The MLE is defined as follows: 

- max ( |  )  (8)  

Where  ( |  ) is the probability distribution function of  given  . When the measurements in 

- [⃗⃗⃗⃗1 ,⃗⃗⃗⃗2 ,…,⃗⃗⃗⃗ ] are identically distributed and independent, the objective function becomes 

(8.1). Taking the log of both sides gives:  

log( ( |  )) = log(Σ =1 (⃗⃗⃗ | ))  (8.1) 

In the case where G is a measurement matrix and ϵ ∼  (ϵ; 0, Σ), measurements   can be written as: 

- [ 1, 2,…, ] =  +    (8.2) 

The measurements thus follow a Gaussian distribution  (  ; ,Σ) as shown in (6). 

(  | ) = (  ; ,Σ) (8.3) 

The log of the multivariate Gaussian distribution removes the exponential from the equation giving the following equation. 

log( (  | )) ∝ −(  − ) Σ−1(  −  ) (8.4) 

To find the signal source the maximum needs to be found. The equation can be further simplified by removing the negative sign and finding the minimum. This is shown in the equation below. 

- max(−log( (  |  ))) = min(  − ( )) Σ−1 (  − ( ))  (8.5) ⃗⃗  ⃗⃗ 

The measurement vector in this case is   , which gives. 

- = min(  − ( )) Σ−1 (  − ( ))  (9) ⃗⃗ 

The objective function in this form is in fact a non-linear least squares (NLS) problem. 

3. Grid<a name="_page5_x69.00_y435.00"></a> Search  

To develop a map of the objective function (the nonlinear least squares problem) the Doppler shift of the signal received from each satellite needed to be found. This was done using Equation (3). As the frequency of the source is not known the difference in Doppler shift between satellites 1 and 2, 3, and 4 uses Equation 3.1, this also allows for the generation of the source vector at (5, 10, 0) as shown in Vector 3.2. To get the measured frequency vector for each grid point the coordinates of the 40 by 40 grid were iterated through and solved as shown in Vector 4. At each coordinate the expected frequency vector and the measured frequency vector were substituted into Equation (6) to give. This would then be recorded into a matrix that was then plotted as a contour map (using Python) to get a visual representation of the objective function. The code that calculates for each plot is shown in the Python code in Appendix A. The contour plot is shown in Figure 1 below. 

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.003.jpeg)

*Figure 1: Contour plot of the signal with no noise.* 

Figure 1 shows that there is a minimum at  (5,10) which is the location of the source. Although this contour graph does not have noise, the following contour map in[ Figure 2 ](#_page6_x69.00_y677.00)does. The noise added to the simulation for Figure 2 was 1Hz which resulted in the minimum being found 600 meters away from the true minimum. 

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.004.jpeg)

<a name="_page6_x69.00_y677.00"></a>*Figure 2: Contour plot of the signal with 1Hz of noise.* 

[^3]To find the minimum within the data a grid search algorithm was implemented. This is done by iterating over each grid point and checking if the new grid point is smaller. Once every grid point is checked, the coordinates at that minimum are returned. This algorithm was implemented in the function for calculating as it already steps through each coordinate with a step size of 0.01. This allows a 0.01-degree accuracy/resolution. The algorithm is shown in the Python code in Appendix B. The grid search algorithm was tested with different noise levels. The results of this test are shown in Table 1.  

*Table 1: Minima found with different noise levels.* 

**NOSIE DEVIATION (HZ)  MINIMA (** , , ) **![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.005.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.006.png)**0  (5, 10, 0) 

1  (5, 10, 0) 

||2 ||(5, 10, 0) ||
| :- | - | :- | - | :- |
|4 ||(5, 10.1, 0) |||
||8 ||(4.9, 9.9, 0) ||
|16 ||(5.1, 10.1, 0) |||

4. Find<a name="_page7_x69.00_y311.00"></a> FDOA With Unknown Frequency 

In a real-world situation, the original frequency of the transmitter is unknown which means that the wavelength of the signal cannot be found. To get around this an estimation of the original frequency needs to be made to get an estimate of the wavelength. This can be obtained by evaluating the frequency at satellite j using the following equation. 

1( 0 − )

- + ,  = 1,2,3,4 (10) 
- 0 − ‖ 0

The frequencies at each satellite can then be averaged to get an estimate of the frequency as shown in equation 11.  

- [^4] + 2 + 3 + 4 (11) 4
- ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.008.png) (12.1) 

This can then be substituted back into equation 3.1 to give the frequency difference between satellites as shown in equation 13. 

(⃗⃗⃗⃗0  − ⃗⃗ ) (⃗⃗⃗⃗0  − ⃗⃗⃗1 ) 1

- ` `− + ( − ),  = 2,3,4  (13) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.009.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.010.png)

1 ‖⃗⃗⃗⃗0  − ⃗⃗ ‖ ‖⃗⃗⃗⃗0  − ⃗⃗⃗1 ‖ 1

This now allows for the frequency difference between satellites to be found so the objective function can be formulated. Results for this frequency estimation can be seen in Figure 3. This was also verified using the grid search algorithm to check the minimum was at (5,10) which it was. 

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.011.jpeg)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.012.jpeg)

*Figure 3: Two plots showing the difference between plots of the know frequency and the estimated frequency.* 

**Note:** Using more FDOA measurements such as  23, 13, 43 does not result in any performance improvements. The additional measurements are the result of values being shifted but the data in the matrix retains the same information with a different reference satellite, this can be shown mathematically: 

−1 1 0 0 1 0 1 −1 0 1

(−1 0 1 0) ( 2) = ( 21) = (1 0 −1 0) ( 2) = ( 1233)  (13.1) 

−1 0 0 1 3 4311 0 0 −1 1 3 43

4 4

5. Gradient<a name="_page9_x69.00_y72.00"></a> Derivation 

To apply gradient-based optimization methods, the gradient must be computed or approximated. From Equation 3.1 we have: 

1 (⃗⃗⃗⃗0  − ⃗⃗ ) 1 (⃗⃗⃗⃗0  − ⃗⃗⃗1 ) 1 1

1 = ‖⃗⃗⃗⃗0  − ⃗⃗ ‖ − ‖⃗⃗⃗⃗0  − ⃗⃗⃗1 ‖ , ⃗ 0 = 0  (14) 

` `Let  

- ⃗ 0 − ⃗⃗   (14.1) 

⃗⃗⃗⃗1  = ⃗  0 − ⃗⃗⃗1   (14.2) 

1 1 ⃗⃗⃗⃗1  1 1 1 1 1 ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.013.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.014.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.015.png)

1 = ‖  ‖ − ‖⃗⃗⃗⃗1 ‖ , ⃗ 0 ,    ⃗ 0 = ⃗ 0 − ⃗ 0 14.3) 

(

Let 

1![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.016.png)

- ‖ = 2 = √∑ 2 (14.4) 

1 = 1 (  1 − 1 (⃗⃗⃗⃗1  1 1) (14.5) 

)

2 2

1  −12

[^5] = ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.017.png)( 1 )  (14.6) 

Using product rule and chain rule gives: 

⃗ 01 = (1 ( −1[^6] − )) (14.7) 

3 ⃗ 0![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.018.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.019.png)

2 2

( 12) = (√∑ 2)2 = ∑ 2 = , = ∑2  (14.8) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.020.png)

2![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.021.png)

1

- (⃗ 0 −  ) = [1]  (14.9) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.022.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.023.png)

⃗ 0 ⃗ 0

1

1 = (1( − (2  ))) − (1( 1 − 1 1 (2⃗⃗⃗⃗ )))  (14.11) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.024.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.025.png)

⃗ 0 ‖  ‖ 2‖  ‖3 ‖⃗⃗⃗⃗1 ‖ 2‖⃗⃗⃗⃗1 ‖3 1

Substituting   back in gives the derivative: 

⃗ 0 =  (1 (‖⃗ 0 − ⃗⃗ ‖ − (‖⃗ ⃗0 − ⃗⃗ ‖3 (⃗ 0 − ⃗⃗ ))) − (1 (‖⃗  1 − (⃗ 0 − ⃗⃗⃗1 ) 1 (⃗ 0 − ⃗⃗⃗1 ))) (14.12) 

1 − ⃗⃗ )![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.026.png)

0 0 − ⃗⃗⃗1 ‖ ‖⃗  0 − ⃗⃗⃗1 ‖3

The resulting Jacobian can be written in the following form: 

21 21 21

- (⃗ )

(⃗ ) = 31 31 31 = (∇ 2311(⃗ ) )  (14.13) 

- 41(⃗ )

41 41 41![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.027.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.028.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.029.png)

(  )

The derivatives of the coordinate conversion functions (Eq. 2, 2.1, 2.2, 2.3, 2.4) are used to find the derivative of  1 with respect to   (Eq. 1.1). The resulting Jacobian  ( ) is shown below, 

with derived terms. ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.030.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.031.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.032.png)

( ) =  ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.033.png) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.034.png) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.035.png) (15) 

( )![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.036.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.037.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.038.png)

The desired Jacobian for  is: ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.039.png)

(⃗ , ) = (⃗ ) ( ) (15.1) The derivation results follow: 

0 2 cos( )sin( )![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.040.png)

- 2 (15.2) (1 − 2 sin2( ))3

00 = cos( 0)[ cos( 0) − ( 0)sin( 0)]  (15.3) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.041.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.042.png)

[^7]0

- ( 0)cos( 0)sin( 0)  (15.4) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.043.png)

0

0

0 = sin( 0)[ cos( 0) − ( 0)sin( 0)]  (15.5) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.044.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.045.png)00 = ( 0)cos( 0)cos( 0)  (15.6) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.046.png)

0

0 = (1 − 2)( sin( 0) + ( 0)cos( ))  (15.7) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.047.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.048.png)00 = 0  (15.8) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.049.png)

0 0

0 0![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.050.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.051.png)

0 0

( ) = (15.9) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.052.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.053.png)

0 0

0 0

[ 0 0]![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.054.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.055.png)

Note that the H column is omitted as the altitude in degrees (H) is zero. Using Equation 14.13 and Equation 15.9 in Equation 15.1 results in the gradient used for the gradient-based optimization method in section 7. 

6. Derivative<a name="_page11_x69.00_y425.00"></a> Verification 

This section provides a verification of the  21( , ) derivative using the finite-difference and complex step methods. 

The finite-difference method formula is: 

( + ℎ ) − ( )

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.056.png)=  + (ℎ) (16) 

ℎ

From Equation 3.1 and Equation 16 we obtain: 

1 (⃗  + ℎ − ⃗⃗⃗~~2~~ ) ~~2~~ (⃗  + ℎ − ⃗⃗⃗~~1~~ ) ~~1~~ 1 (⃗  − ⃗⃗⃗~~2~~ ) ~~2~~ (⃗  − ⃗⃗⃗~~1~~ ) ~~1~~

( ( − ) − ( − ))

||⃗  + ℎ − ⃗⃗⃗2 || ||⃗  + ℎ − ⃗⃗⃗1 || ||⃗  − ⃗⃗⃗2 || ||⃗  − ⃗⃗⃗1 ||

[^8] = (16.1)

ℎ

(1 ((⃗  + ℎ − ⃗⃗⃗~~2~~ ) ~~2~~ − (⃗  + ℎ − ⃗⃗⃗~~1~~ ) ~~1~~ ) − 1 ((⃗  − ⃗⃗⃗~~2~~ ) ~~2~~ − (⃗  − ⃗⃗⃗~~1~~ ) ~~1~~ ))

||⃗  + ℎ − ⃗⃗⃗2 || ||⃗  + ℎ − ⃗⃗⃗1 || ||⃗  − ⃗⃗⃗2 || ||⃗  − ⃗⃗⃗1 ||

21 = (16.3)

ℎ

The complex step formula is: 

ℑ( ( + ℎ ))

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.057.png)=  + (ℎ2)  (17) 

ℎ

- ℑ (1 ((⃗| |⃗++ℎℎ −~~2~~⃗⃗⃗2 ||~~2~~ − ||⃗  + ℎ 1 (17.1)
- ⃗⃗⃗ ) (⃗  + ℎ − ⃗⃗⃗~~1~~ ) ~~1~~ ))

21 − ⃗⃗⃗ ||

ℎ

- ℑ(1((⃗| |⃗++ℎℎ −−⃗⃗⃗~~2~~ ⃗⃗)⃗ ||~~2~~ − (|⃗ |⃗++ℎℎ ⃗⃗⃗1 ||~~1~~ )) (17.2)
- ⃗⃗⃗~~1~~ )

21 2 −

ℎ

- ℑ( ( ||⃗  + ℎ − ⃗⃗⃗2 ||~~2~~ − (|⃗ |⃗++ℎℎ 1 (17.3) 1 (⃗  + ℎ − ⃗⃗⃗~~2~~ ) − ⃗⃗⃗~~1~~ ) ~~1~~ ))

21 − ⃗⃗⃗ ||

ℎ

Where   is a unit vector in the  ℎ direction (B, L, H for longitude, latitude, and altitude respectively), ℎ is a step size, and ⃗  is the point about which to differentiate (in B, L, H).  The approximation of the derivatives (Eq. 16.1, 16.2, 16.3, 17.1, 17.2, 17.3) and the analytical derivative in Equation 13.12 is compared, with results shown in[ Table 2:](#_page13_x69.00_y97.00) 

<a name="_page13_x69.00_y97.00"></a>*Table 2: Comparison of derivative approximations and analytical solutions.* 



|**Method** |**Derivative** |**Step size**  |**Point (B, L, H)** |**Result** |
| - | - | - | - | - |
|Analytical |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.058.png)|<p>−8</p><p>10</p>|(5, 10, 0) |-0.57761319 |
|Finite Difference |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.059.png)|<p>−8</p><p>10</p>|(5, 10, 0)|-0.5775973477284424 |
|Complex Step |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.060.png)|<p>−200</p><p>10</p>|(5, 10, 0)|-0.5834981372799554 |
|Analytical |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.061.png)|<p>−8</p><p>10</p>|(5, 10, 0)|3\.90670032 |
|Finite Difference |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.062.png)|<p>−8</p><p>10</p>|(5, 10, 0)|3\.9067515444912715 |
|Complex Step |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.063.png)|<p>−200</p><p>10</p>|(5, 10, 0)|3\.905921531938494 |
|Analytical |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.064.png)|<p>−8</p><p>10</p>|(5, 10, 0)|-3.58224626 |
|Finite Difference |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.065.png)|<p>−8</p><p>10</p>|(5, 10, 0)|-3.5822438348986907 |
|Complex Step |21![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.066.png)|<p>−200</p><p>10</p>|(5, 10, 0)|-3.5822462633136434 |

It is interesting that the complex step method is very accurate on the altitude derivative but less accurate on the latitude and longitude than the finite difference method. The error resulting from both approximations also appears to be larger than the  (ℎ) and  (ℎ2) remainder terms implied by Equations 16 & 17. This could be due to the computer’s finite arithmetic or how functions are being handled “behind the scenes” in the Python code. Still, the accuracy of the approximations is enough to verify the  21 derivatives. It would follow that the other FDOA function derivatives 

⃗ 

are correct due to the implementation being the same.  

7. Gradient<a name="_page13_x69.00_y473.00"></a> Based algorithm 

The Optimization method used to find the source location was the Levenberg-Marquardt algorithm (LMA). This is a surrogate-based optimization method that estimates the hessian. LMA was developed to solve non-linear least-square equations. The LMA interpolates between the Gauss-Newton algorithm (GNA) and gradient-based optimization. GNA is very fast at converging to the minimum but is not very accurate once it gets close. This is where the gradient-based algorithm works best. LMA will get close to the solution with GNA and then switch to the steepest decent method to get a more precise solution. A derivation of the update formula (Eq. 18) can be found in Appendix A.  

- ( − )−1 ( − ( ))  (18) 

Where  is the Jacobian of the frequency derivatives for each satellite at the current position,  is the evaluation of the frequency at the source location using equation 4,  ( ) is the evaluation of the frequency at the estimated position that the satellites are measuring,  is the control of the direction ranging between GNA and steepest decent. We can improve the scaling of the function by multiplying by  . Where   is: 

- ( ) (18.1) 

The algorithm works by initializing the code with an initial position,  ,  , and initial position  . The Jacobian and the residual will then be evaluated. The first step is calculated and added to the current position, and the error is then found using Equation 18.2. This is the error for the next step ( ). 

- ‖ ‖2  (18.2) 

The difference in the current error  and the previous error  is then found using Equation 18.3. 

∆= −  (18.3) 

If the difference is negative then the new position will be accepted and e, the Jacobian and the residual will be updated. For every fifth iteration, the  value will be divided by  (the damping factor) making it act more like GNA. If the difference in error is not negative,  will be multiplied by  to make the algorithm act more like a steepest decent. This code will then loop and check if the difference in errors is less than a tolerance of 10−3. If the error is less than the tolerance it will return the position, if not the code will loop again. This process is shown in the flow diagram in Figure 4. 

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.067.jpeg)

*Figure 4: Flow chart of the algorithim.* 

Figure 5 shows the convergence of the LMA algorithm from multiple starting positions. This shows that the algorithm converges to the same position from many different start positions.

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.068.jpeg)

*Figure 5: Plot of the contour and the path of the algorithm.* 

Delta (Eq. 18.3) is used as the stopping criteria for the algorithm. Once delta is less than 10−3 the algorithm stops looping and returns the position. The LMA algorithm converges in 15- 30 iterations. It will evaluate the Jacobian and the residual  + 1 times where k is the number of iterations it takes to converge. The step calculation will be evaluated once every iteration so will have the same number of evaluations as iterations. The relation between the number of evaluations, the difference in the current error and next step error (delta) can be seen in Figure 6 below. 

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.069.jpeg)

*Figure 6: Plot showing the convergence of the algorithm.* 

8. Monty<a name="_page16_x69.00_y97.00"></a> Carlo Comparison 

To determine the accuracy of the grid-search algorithm and LMA, a Monte Carlo simulation is run to estimate the root-mean-squared-error. The formula for the calculation is as follows: ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.070.png)

1

- √ ∑‖ − ‖2 (19)

−1

The equation uses the total number of Monte Carlo simulations, L, and the summation of the 2- norm squared difference between the position vector found by the grid search and the actual position vector of the source. The number of Monte Carlo runs used was 104.  The way that the algorithm was run is shown in the Monte Carlo function in Appendix B. While testing the grid search, noise levels of  = (1, 2, 4, 8, 16) Hz were used. This process gave the results shown in Table 3 below.  

*Table 3: Table of the Monty Carlo results.* 

**NOISE LEVEL (HZ)  LMA (KM)  GRID SEARCH (KM) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.071.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.072.png)1**  1.0941220938193204  0.6197978901842434 ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.073.png)

**2**  2.1879750170699395  1.0849477826620166 

||**4** ||4\.365063017360187 |||2\.2123866587402907 ||
| :- | - | :- | - | :- | :- | - | :- |
|**8** ||8\.725234012378314 |||4\.227109425169162 |||
||**16** ||17\.44703333621689 |||8\.175176490117455 ||
The results from the RMSE were plotted in Figure 7. This shows that the results are linear and that the LMA algorithm follows the same trend as the grid search. The relationship between the grid search and the LMA was that the error doubled what the grid search error was. 

Geolocational Accuracy![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.074.png)

20 15 10

|||||
| :- | :- | :- | :- |
|||||
|||||
|||||
RMSE (km) 5

0

0 5 10 15 20 Noise Level (Hz)

LMA (km) Grid search (km)

*Figure 7 7: Plot showing the RMSE evaluation of the grid search algorithm and the LMA algorithm in km.* 

<a name="_page18_x69.00_y72.00"></a>Appendix A: Levenberg-Marquardt Update Derivation 

Let  be the true function value,  ( , ) the next function estimate and  ( 0, 0) the current function estimate.  

Let   denote ( , ) and ⃗⃗⃗⃗0  denote ( 0, 0). 

( − ( , )) ( − ( , )) ( 1) ( , ) = ( 0, 0) + ( 0, 0)( −− 00) ( 1.1) 

( −  ( 0, 0) + ( 0, 0)( −− 00)) ( −  ( 0, 0) + ( 0, 0)( −− 00))  ( 1.2) 

- − (⃗⃗⃗⃗0 ) + (⃗⃗⃗⃗0 ) ⃗⃗⃗⃗0   ( 1.3) 

( − (⃗⃗⃗⃗0 )  ) ( − (⃗⃗⃗⃗0 )  ) ( 1.4) 

- − (⃗⃗⃗⃗0 )   ( 1.5) 

![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.075.png)=![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.076.png) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.077.png) ( 1.6) 

- ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.078.png)Σ 2 → 2 ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.079.png) ( 1.7) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.080.png)= − (⃗⃗⃗⃗0 ) +![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.081.png) ( 1.8) ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.082.png)= 0  ( 1.9) 
- ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.083.png)= 2( − (⃗⃗⃗⃗0 )  )(− (⃗⃗⃗⃗0 )) ( 1.10) 
  - −2 (⃗⃗⃗⃗0 ) + 2 (⃗⃗⃗⃗0 )  (⃗⃗⃗⃗0 ) ( 1.11) 
    - (⃗⃗⃗⃗0 )  (⃗⃗⃗⃗0 ) = (⃗⃗⃗⃗0 ) ( 1.12) 
- ( (⃗⃗⃗⃗0 ) (⃗⃗⃗⃗0 ))−1 (⃗⃗⃗⃗0 ) ( 1.13) 

( ) = ( ( 0, 0) ( 0, 0))−1 ( 0, 0) ( 1.14) For the Levenberg-Marquardt Algorithm a factor of  is added, where  is a scaling constant. 

- ( ( 0, 0) ( 0, 0)) ( 1.15) 

( ) = −( ( 0, 0) ( 0, 0) + )−1 ( 0, 0)  ( 1.16) 

<a name="_page20_x69.00_y72.00"></a>Appendix B: Algorithm and Grid Search Code. 

import numpy as np ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.084.png)

import matplotlib.pyplot as plt import random 

fs = 1\*10\*\*9      #frequency of transmitter wl = (3\*10\*\*5)/fs #wave length  

#source location 

u0 = np.array([[2], [1], [0]]) 

#Data matrix 

theta = np.empty(shape = (501, 501))  

#Weighting function Q = np.array([[2, 1, 1],             [1, 2, 1], 

`            `[1, 1, 2]]) 

Qinv = np.linalg.inv(Q) 

#making the noise vector def nf(lvl): 

`    `n1 = random.gauss(0, lvl)     n2 = random.gauss(0, lvl)     n3 = random.gauss(0, lvl) 

n4 = random.gauss(0, lvl) 

nf = np.array([[n2-n1], [n3-n1], [n4-n1]]) return nf 

#satalites 

s1 = np.array([[7378.1, 0, 0]]).T 

s1d = np.array([[0.0001, 4.4995, 5.3623]]).T s2 = np.array([[7377.5, 100, 0]]).T 

s2d = np.array([[-0.0671, 4.9493, 4.9497]]).T s3 = np.array([[7377.5, -100, 0]]).T 

s3d = np.array([[0.0610, 4.4991, 5.3623]]).T s4 = np.array([[7377.5, 0, 100]]).T 

s4d = np.array([[-0.0777, 4.0150, 5.7335]]).T 

#Should convert from ECEF to LLA def LLA\_to\_ECEF(lat, long): 

`    `R0 = 6378.137 

e = 0.081819198425 ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.085.png)

lat = np.radians(lat) 

long = np.radians(long) 

Re = R0/(np.sqrt(1-(e)\*\*2\*(np.sin(lat))\*\*2)) x = (Re\*np.cos(lat)\*np.cos(long)) 

y = (Re\*np.cos(lat)\*np.sin(long)) 

z = ((1-e\*\*2)\*Re\*np.sin(lat)) 

u = np.array([[x], [y], [z]]) 

return u 

u0 = LLA\_to\_ECEF(u0[0,0], u0[1,0])[:2] 

#grid search doppler shift 

def dopler\_shift(lat, long): 

`    `u = LLA\_to\_ECEF(lat, long) 

`    `f0 = 1/wl\*np.dot((u - s1).T , s1d)/(np.linalg.norm(u-s1)) 

`    `f1 = 1/wl\*np.dot((u - s2).T , s2d)/(np.linalg.norm(u-s2)) - f0     f2 = 1/wl\*np.dot((u - s3).T , s3d)/(np.linalg.norm(u-s3)) - f0     f3 = 1/wl\*np.dot((u - s4).T , s4d)/(np.linalg.norm(u-s4)) - f0      f = np.array([f1[0], f2[0], f3[0]]) 

`    `return f 

#grid search calc 

def grid\_search(): 

`    `long = 0 

`    `lat  = 0 

`    `min = float('inf') 

`    `coord = 0 

`    `noise = nf(2)   

`    `f = dopler\_shift(2, 1) 

`    `for i in range(0, 501): 

`        `for j in range(0, 501): 

`            `g = dopler\_shift(lat, long) 

`            `g += noise 

`            `theta[j][i] = np.log((f - g).T @ Qinv @ (f - g)) 

`            `if theta[j][i] < min: #grid search stuff 

`                `min = theta[j][i] 

`                `coord = np.array([[j], [i], [0]])\*0.01 

`            `lat += 0.01 

`        `lat = 0 

`        `long += 0.01 

`    `xyz = LLA\_to\_ECEF(coord[0,0], coord[1,0])     return coord, xyz[:2] 

#monte carlo simulation to find the accuracy of the data ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.086.png)def monte\_carlo(): 

`    `sumation = 0 

`    `L = 500 

`    `for l in range(0, L): 

`        `coord, xyz = grid\_search() 

`        `sumation += np.linalg.norm(xyz - u0)\*\*2         print(l) 

`    `RMSE = np.sqrt((1/L)\*sumation) 

`    `print(RMSE) 

#plotting stuff 

def graph(): 

`    `coord = grid\_search() 

`    `tmin = np.min(theta) 

`    `tmax = np.max(theta) 

`    `levels = np.linspace(tmin, tmax, 100)     feature\_x = np.arange(0, 5.01, 0.01)      feature\_y = np.arange(0, 5.01, 0.01)  

[X, Y] = np.meshgrid(feature\_x, feature\_y) contour = plt.contourf(X, Y, theta, levels = levels)  plt.colorbar(contour, label='log(theta-values)') 

plt.xlabel("Longitude (Degrees)") 

plt.ylabel("Latitude (Degrees)") 

plt.title("FDOA plot") 

#print("Minmum is found at: ", coord.T) #plt.annotate('Min',xy=(10,5),xytext=(5,10),arrowprops={}) 

plt.show() 

def main(): 

`    `monte\_carlo()     #graph() 

main() 

This code is also in our Git hub in project.py. https://github.com/Mangin4/ENEL445\_Optimization ![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.087.png)![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.088.png)

[^1]: X is the known measurement vector,  is the mean vector, and the Σ is the covariance matrix 

    composed of  2 . 
[^2]: <a name="_page3_x69.00_y652.00"></a>.3 Assumptions   The assumptions that were made for this simulation are as follows: 

    ￿  Altitude is zero and does not change. 
[^3]: Using this estimate of the frequency the wavelength can be found as shown below where  is the 

    speed of light and  is the wavelength. 
[^4]: =![](Aspose.Words.88ada844-83fa-40da-8cba-dd702e56bbff.007.png) (12) 

    This can be rearranged to give lambda as shown in equation 12.1. 
[^5]: 1 1 1

    ⃗ 01 = ( ( 1 − 3 (2  ))) [1] − ( ( 11 − 31  (2 1))) [1]  (14.10) 

    2 2 2 1 2 2 2 1
[^6]: 
[^7]: (1 ((⃗  + ℎ − ⃗⃗⃗~~2~~ ) ~~2~~ − (⃗  + ℎ − ⃗⃗⃗~~1~~ ) ~~1~~ ) − 1 ((⃗  − ⃗⃗⃗~~2~~ ) ~~2~~ − (⃗  − ⃗⃗⃗~~1~~ ) ~~1~~ ))

    ||⃗  + ℎ − ⃗⃗⃗2 || ||⃗  + ℎ − ⃗⃗⃗1 || ||⃗  − ⃗⃗⃗2 || ||⃗  − ⃗⃗⃗1 ||
[^8]: = (16.2)

    ℎ
