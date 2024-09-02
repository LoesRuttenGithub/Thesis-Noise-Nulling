# Online thesis text


Here I will update my work on my thesis as we go. It is intended to be structured but not as formal as the final text. To be discussed: Which content can be in a public folder.

## Table of contents
[Literature_List](#Literature_List)

[Meetings](#Meetings)

[Planning](#Planning)

[Quick_Summary](#Quick_Summary)


## Literature_List
<details>
<summary> Quantz et al. 2022: LIFE Paper 1: Improved exoplanet detection yield estimates </summary>

Description of LIFE space mission with 4 formation flying infrared telescopes. Paper presents an intstrument simulator with all sources of astrophysical noise (but not yet including instrumental noise) which was coupled with the P-Pop Monte Carlo tool of a synthetic exoplanet population. The gain (yield per unit time) was averaged over 500 MC realisations.

Based on 4 2m apertures, it was estimated that within 20 p of the sun, 25-45 rocky exoplanets within the habitable zones could be detected, and doubling these numbers with 3,5m apertures. In a first observation, estimates of radii and effective temperatures of exoplanets could be made, followed by estimates of thermal emission spectra in a second visit.

It is known from Kepler, Tess and Radial Velocity surveys that planets similar to Earth should be very abundant. In the future, James Webb Space Telescope may reveal if planets around red dwarfs can retain their atmospheres despite high levels of activity of their host stars. 

ESA's Ariel will provide transmission and emission spectra of exoplanets with warm hydrogen-dominated atmospheres. The ~30m ground based Extremely Large Telescopes will be able to detect thermal spectra of small planets around nearby stars using the mid-infrared METIS spectrograph. In addition, the PCS and HIRES will detect via reflected light.

The US-led Habitable Exoplanet Observatory HabEx and the Large UV/Optical/IR Surveyor LUVOIR  is intended to study the atmospheres of several dozen small exoplanets in Habitable Zones via reflected light.

Several high-precision RV instruments are planned: CARMENES, NIRPS, ESPRESSO, MAROON-X, HARPS3, EXPRES, so it is expected that a significant fraction of exoplanets within 20 pc will be uncovered. In that case, more LIFE observation time can go to characterisation.

LIFE explores the development of a space-based mid-infrared nulling interferometer within ESA's Voyage 2050 program, building upon studies of the Darwin mission and Terrestrial Planet Finder - Interferometer. The study would exist of a direct detection phase, followed by a characterisation phase. This paper focusses on the detection phase. 

Advantages of working in the mid-infrared: More direct constraints on temperature and size of the objects, and spectrum includes molecular absorption lines including biosignatures.

The yield simulation includes single main sequence stars and wide separation binaries <span style="color:green">(what happens to close binaries?)</span> and NASA's ExoPaG SAG13 underlying exoplanet occurrence rates. The signal as measured from a planetary system in a 6:1 X-array configuration of nulling-interferometry is simulated using LIFEsim. The short baselines are 'nulling baselines' responsible for dark fringes, the long baselines are 'imaging baselines' responsible for high frequency modulation of the transmission map perpendicular to the modulation map. The 6:1 baseline seems more robust against instability noise than the 2:1 ratio <span style="color:green"> (!) </span>

The simulations assumed separations between the spacecrafts of 10 to 600m, 2m aperture size and 5% optical throughput. All major astrophysical noise terms were included: Planetary photon noise assuming black-body emission, photon noise due to stellar leakage into the dark fringe, noise from exozodi disks, assumed to be optically thin, smooth and face-on, and photon noise from local zodiacal light aka dust in our Solar System, which is overcome by pointing in anti-sunward direction. 

Instrumental effects were neglected. What is unknown for now is 1) the impact of phase and/or amplitude variations as systematic noise sources, and 2) thermal background from the aperture mirrors. The S/N calculations assume the use of single-mode fibers (?). 

The detection criterion was S/N>7 where the S/N was taken as the square root of the sum of the instrumental and astrophysical noise. Here the instrumental noise was assumed lower than the astrophysical noise. (Why astrophysical signal of S/N>5?). It seems that if the S/N is sufficiently high, the radius, effective temperature and a very rough SED can be estimated by first observation to first order. 

There are two ways in which the scenario can be optimised: maximising the *total* number of detected exoplanets, or maximising the number of **rocky exoplanets within the eHZ Habitable Zone*, which leads to a decrease in the total number of detectable planets. Temperate rocky exoplanets in the eHZ are fainter than objects closer to the star and therefore require >3x more integration time. In absolute numbers, rocky eHZ exoplanets occur more around M dwarfs because M dwarfs are more numerous in the solar neighbourhood. 

Comparing the theoretical number of planets in the simulation to those that are detected, the required sensitivity is the main factor rather than the spatial resolution. The detections approximately double/half with the aperture size from 1/2/3.5m. The yield is comparable to LUVOIR. 

It has been argued by previous studies that 30 to 50 exoplanets need to be studied to get statistically robust results on the fraction of habitable rocky HZ exoplanets. It is unknown if exoplanets around M-dwarfs can sustain atmospheres, but hopefully JWST will give insights on this. Depending on this result, LIFE will be an excellent next step with a detection bias around M dwarfs, and otherwise M stars should be deprioritised. HabEx and LUVOIR have detection bias for solar-type stars. 

</details>


<details>
<summary> Laugier et al. 2022: Asgard/NOTT: L-band nulling interferometry at the VLTI </summary>

NOTT = Nulling Observations of ExoplaneTs and dusT, beam combiner for VLTI which aims to resolve young giant exoplanets around nearby stars down to 5 mas with $10^5$ contrast in the L-band. Down to mag 7, errors are dominated by correlated instrumental errors. Beyond that magnitude, the thermal background of the relay and telescope are dominant. Nulling interferometry is based on the technique to tune phases f acombined beam, creating a dark fringe. The path length difference must be matched within a fraction of the wavelength. The technique was first proposed by Bracewell (1978) and in particular the sin-chop architecture of Angel & Woolf (1997) is considered. Phase errors due to errors in the optical path difference caused by turbulence in the atmosphere can be cancelled out using *kernel-nulling*, as proposed by Martinache & Ireland (2018). 

*SCIFIsim* was developed to simulate the effects of instrumental noise on nulling-interferometry. It takes into account wavelength dependence, via chromatic combination, which means that the coupler matrix, representing the architecture, is wavelength dependent. For each wavelength bin and microsecond order subexposure the following steps are executed by SCIFYsim:
- A vector of the the real spectrum collected by the telescope is transferred into a complex phasor, encoding also the geometric position, optical abberations , transmission map of the waveguide and optical transmission
- This is combined with the complex amplitude of each source, to compute the complex amplitude of light entering the integrated beam combiners
- This is multiplied with the matrix of the beam combiner to get the complex amplitude at the outputs, namely the exact signal, the errors and the total signal.
The spectrograph modules convolves the intensity at the output with a spectroscopic point spread function via a MonteCarlo simulation to form a 2D array of pixels. Throughput and thermal background are evaluated based on transmisson-emission objects. Thermal background is computed via Planck's law. 

Main sources of errors included in the model:
- Aberrations in the beam produce variations in the focal plane complex amplitude upon injection into single-mode waveguides
- Residual optical path errors, computed via the GRAVITY fringe tracker
- Internal combiner chromatic errors, modelled using order 6 polynomial, which leads to imperfect nulling. Mitigated using variable thickness planes.
- Not included: Longitudinal atmospheric dispersion, assumed to be compensated by ZnSe corrector plates to first order.
There is a higher correlation in erros closer to the nominal scheme. 

The distribution of errors is studied, and the Shapiro-Wilk test for 'Gaussian noise' fails beyond 400 samples for 3 seconds detector integration time. But Gaussian model is still practical to use. 

The covariance matrix of the total errors is the sum of read-out noise, photon noise and instrumental noise. It is assumed that detection integration time measurements are statistically independent. The data is averaged over several chunks or 'observing blocks' whose errors are assumed uncorrelated.

Result: the inner part of the spectrum is correlated together. the outer part of the spectrum is anti-correlated. 

A series of statistical detection tests is implemented based on Ceau et al. 2019. The covariance matrix is 'block-diagonal' with size $n_{chunks}\times n_{channels} \times n_{outputs}$. The measured signal can be thought of as the sum of the theoretical signal plus errors: $z=z_{theory}+\epsilon'$. This can be multiplied with a whitening matrix such that the noise is normally distributed on the diagonal: $x=W z = \Sigma^{-1/2}$ so that $y=x+\epsilon$ with $\epsilon \sim N(0,I)$. The target amplitude must satisfy some threshold $\xi$ of a probability to measure a false alarm and the probability of a true detection. The two statistical test are the *energy detector test* and the *Neyman-Pearson test*.

Results: A simulation for a K-dwarf gave sensitivity maps for different magnitudes. In the center, where the lobes of transition maps of different wavelengths overlap, there is most sensitivity. Except for bright stars, for which in the center, the planet pattern looks similar to the instrumental noise, therefore the performance in outer regions is better.
</details>

<details>
<summary> Martinache & Ireland 2018: Kernel-nulling for a robust direct interferometric detection of extrasolar planets </summary>
Nulling interferometry for exoplanet has great potential, but the performance is limited by instrumental phase and background errors that keep the instrument off the null. A modified nuller architecture is presented that is robust against piston excursions. A concept of kernel is applied to the outputs to make them robust to second order pupil phase error, under the name VIKiNG: the VLTI Infrared Kernel Nulling.

The resolving power is limited by diffraction, which introduces featuerwhose contribution tot he data doinates that of faint structures of interest. High contrast image uses coronography, but is limited by the atmosphere, even despite extreme adaptoive optics sistems. The contrast of a speckle at wavelength labda is directly related to the amplitude of the modulaton. This can be be inverted to the requir4ement that a wavefront quality requirement better than 25 nanometer is needed to to reach a contrast of $10^{-6}$ inn the H-band (1.6 um). Recent studies used angular differential imaging as post-processing technique to disentangle the diffraction features from genuine structures. The approach in this paper is to improve the *pre-processing* phase, to make the architecture more robust against small perturbations.

In nulling interferometry, the main constraints are residual background errors, and at high-frequency, residual phase errors. The idea is to improve the rejection of the nuller, for example by combining multiple apertures like Angel & Woolf 1997. This paper presents a true self-calibration technique, which takes advantage of the coupling between atmospheric induced piston errors along baselines forming a triangle, to produce a subset of clean observable quantities robust against residual piston errors from a finite set of polluted raw phase measurements. The trick is to look for linear combinations of polluted data that reside in a space orthogonal to the source of perturbation.

A nuller design with one bright output and three dark ones is considered.. The four-beam nuller can be represented by a 4x4 matrix N, producing the outputs from the inputs, normalised such that the total flux is preserved. The raw interferometric phase per baseline is linearly related to instrumental phase, while the *output* of a nuller is a quadratic function of piston excursions. Setting the phase of one aperture as reference phase, one can construct a three-parameter correlated piston vector p. By first order the Taylor expansion of the pistons dependence on the input electric field is $E+k=e^{-j \phi_k} \approx 1 - j \phi_k$. The three nulled intensities can then be written as $x=\frac{1}{4} \times \left[ (+\phi_1-\phi_2-\phi_3)^2,  (-\phi_1+\phi_2-\phi_3)^2,  (-\phi_1-\phi_2+\phi_3)^2 \right]$. The piston-induced leak of the nuller is a function of six parameters: $\phi_k^2$ and $\phi_k \times \phi_{l \neq k}$. This problem is underconstrained, but this can be solved with help of a split-and-mix operation, that breaks down each nuller output into two non-symmetric outputs with help of a pre-defined phase offset $\theta$ that helps to discriminate variations in the two parts of the complex visibilities. The final function now is a six-component intensity vector $x$. 

In reality, the nuller also deals with a temporal variation due to fluctuating background and detector noise. For the longterm amplitude balance between inputs, modulated phase shifters are inserted between the coupling and split-and-mix operations. However, when these background fluctuations are ignored the nulling and sensing functions can be combined into a single six-by-four operator $M$ trhat transforms four telescope inputs to six nulled output complex amplitudes. The difference with a classical nuller is that a classical nuller has symmetric response curves with respect to the axis, and this nuller has anti-symmetric response curves so that positive and negative offset positions can be distinguished. This gives stronger constraints on the position of a companion from a single observation. 

The goal is to build a sub-set of observables that exhibit robustness against residual piston errors. When everything is in phase, the system sits on the null and the first order derivative terms of phase and amplitude are zero but the leaked intensity $\Delta x$ will be dominated by six second order terms, measuring the local curvature: $\Delta x = A \cdot \left[ \frac{\partial^2 x}{\partial \phi_1^2}, \frac{\partial^2 x}{\partial \phi_2^2}, \frac{\partial^2 x}{\partial \phi_3^2}, \frac{\partial^2 x}{\partial \phi_1 \partial \phi_2} ,\frac{\partial x}{\partial \phi_1 \partial \phi_3} ,\frac{\partial^2 x}{\partial \phi_2 \partial \phi_3} \right] ^T$. It may be possible to identify a sub-set of linear combinations of rows of A such that a kernel operator $K$ exists for which $K \cdot A =0$, which then also has the property that the observables $K \cdot x$ are independent of second-order phase differences in the pupil plane.  

One approach to construct kernel operator K is to compute the singular value decomposition of A, where the rank of matrix A corresponds to the number of independent closure-phases. For phase shift $\theta=\pi/2$ K can be determined by hand, which erases all second instrumental phase errors for the primary observables $y=K\dot x$, whose variation along position is shown in Figure 3. This property is independent from the presence of phase noise during or in between integration times.

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/MartinacheIreland.png" width="400" height="320">

As a next step, 2D transmission maps of the modified nullers are modelled, which are half the flux of the original 3 classical nullers, and have anti-symmetric properties. In the absence of coupling losses, the 6 modified nullers should combine back to the same transmission as the original 3 nullers. The behaviour of the modified nuller is compared to the original architecture for the layout of the VLTI for a companion with contrast 0.01 at coordinates where the sensitivity of the nuller is near optimal, in the presence of 50 nm residual piston excursions. The theoretical value of the null depth should be a Dirac delta around the dotted lines, but in the presence of the residual piston error, become three plotted skewed distributions. In real life, this distribution would still be convolved with a Gaussian due to background and residual target shot-noise. The middle figure shows six instead of three similarly skewed distributions. The figure on the right shows the kernel outputs $y$, they are symmetric with uncertainties proportional to the cube of the phase errors. This means that 10% intensity fluctuations on the input translate to errors smaller than 0.001 on the kernel outputs.

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/MartinacheIreland2.png" width="1000" height="270">

The uncertainty is then analysed using a least squares solution. Depending on the configuration, the uncertainty in the best estimate for the contrast c is $\sigma_c=\frac{1}{|m|} \sigma_k$ where the median ratio is 0.8 but goes up to 1000 near the null. 

Four key fundamental sources of uncertainty constitute kernel-uncertainty $\sigma_k$:
- *Fringe tracking phase errors:* modelled as white power spectrum up to cut-off frequency $\Delta_{\nu_{FT}}^{-1/2}$ 
- *Cross-terms between fringe-tracking phase errors and intensity fluctuations on other telescopes:* a second order term, depending on the intensity fluctuation of each telescope and capacity of the adaptive optics. Their simulation suggests that under realistic (2%) intensity fluctuations, these cross-terms do not affect the performance, as also highlighted by Lay 2004.
- *Thermal background:* scales with power 1/2 of background flux, based on the Bose-Einstein distribution
- *Residual target photon noise:* scales with power -1/2 of the target flux

Assuming non-ideal but good observing conditions, the design is robust against residual wavefront aberrations and photometric fluctuations, with errors dominated by third order input phase and intensity errors.

When applying the analysis to the VIKiNG survey, it is calculated that contrasts of $10^{-5}$ can be achieved for targets brighter than $M_L=6$, reaching SNR=5 within two hours observing time, which means that a dozen nearby planets could be discovered.


</details>

<details>
<summary> Ceau et al. 2019: Kernel-phase detection limits (JWST) </summary>

Statistical tests are presented to derive contrast-detection limits for the JWST Near Infrared Imager and Slitless Spectrograph which will do kernel-phase analysis on faint brown dwarfs. 

Several statistically independent observable quantities are constructed from a collection of aberration-robust kernel phases. 

The NIRISS aperture masking interferometer observing mode takes advantage of self-calibrating observable quantities called closure phases. A *kernel phase* generalises this idea to apertures of arbitrary shapes and is reliable when aberration are withing one radian. They exploit the phase of the Fourier transform of the image, represented in the complex visibility vector $v$, with the phase vector $\Phi$ defined as its argument: $\Phi=\angle v$. 

In the small aberration regime, a linear model relates the measured phase to the true phase $\Phi_0$ and the aberration phase present across the aperture $\phi$ via a phase transfer matrix $A$ dictating how the aberrations propagate into the Fourier phase via $\Phi =\Phi_0+A \phi$. A has the dimensions of the number of subapertures times the number of distinct baselines. The kernel matrix $K$ is the corresponding matrix that verifies $KA=0$, which can be constructed using singular value decomposition, and cancels the phase perturbations to the first order. After finding the kernel matrix, the vector of kernel phases can be extracted via $k=K\Phi$. 

The main question in the paper is, given a data image, how likely is it that a companion is present. Generally, a statistical test tests if a certain test statistic $T(y)$ which is a function of the data $y$ is above or below a certain threshold $\xi$ which represents the rejection of the null hypothesis. The performance of a detect test is given by its probability of false alarm ($p_{FA}$, as low as possible) and its probability of detection ($p_{det}$, as high as possible), and the power of a test is its $p_{det}$ for a given $p_{FA}$, often represented in a ROC curve. 

The noise on the kernels can be modelled by a correlated Gaussian distribution with covariance $\Sigma$. If the matrix is known, the kernel phases can be whitened into an uncorrelated vector with standard normally distributed noise. 

Different tests are considered. 
The likelihood ratio or Neyman-Pearson test requires the companion signature $x$ to be known. The test is then the ratio between dot product of the data $y$ with the signature $x$ over the data $y$ without signature ($x=0$) compared to a threshold. If $x$ is not known, this requirement can be omitted by testing it on all possible values of $x$ within a given parameter space. The generalised likelihood ratio is then a comparison between the maximum likelihood value for $x$, which must then exceed a certain threshold $\eta$. 

If nothing is known about the signature, instead one can use an energy detector test, which uses the measured squared norm of the signal $||y||^2$. As the test does not exploit any prior knowledge of the structure of the object, it is a lower bound for the detection performance. 

Another test statistic is the operational binary test $T_B$, which compares the chi squared value of the some model fitted on data y compared to no model fitted on data y (null hypothesis). The value $T_B$ quantifies how much the sum of squared residuals is reduced compared to the null hypothesis. 


The tests are applied to a simulation of a binary system, represented by two dirac functions with different contrasts, at the origin and at coordinates $\alpha, \beta$. The maximum likelihood is found numerically, avoiding getting stuck in a local minimum. 

Noises are modelled as 1) random noises via statistical errors corresponsing to the coraviance matrix and 2) systematic errors that remain even after the subtraction of kernal phases, due to wavefront drifts bewteen observations. 

Experience shows that the covariance matrix can be accurately estimated by modelling noise from photons, readout and dark current. This can be used to whiten the noise. The covariance is estimated from Monte Carlo simulations. To account for the impact of unknown calibration errors on the contrast detection limits, the residuals of OPS maps are added to the diagonal of the covariance. 

The results show that the parameters can be estimated but that the scatter becomes larger for lower S/N. At angular separations below the diffraction limit, estimates for the contrast and angular eparaiton are strongly correlated so they cannot be constrained simultaneously in a single visit, although this may be lifted by taking a second measurement at an epoch when the companion has moved outside the diffraction limit. 

The theoretical and simulated ROC curves are plotted for the different tests. The Neyman-Pearson is the most powerful test, giving an upper performance bound. The Energy test is the least powerful, providing a lower bound. The operational binary test lies in between, closer to the upper bound. The simulated tests perform close to the theoretical maximum, and also perform well when wavefront errors are introduced. 

For the faint dwarfs considered thus far, the dominating noise effect what that of dark current and read out noise rather than photon noise of the central object. In the case of a bright object, the contrast detection limit is affected by photon noise for separations smaller than 500 mas and then plateaus as dark current and readout noise dominate. When wavefront drift is taken into account for the bright case, 85% of the error comes from calibration errors and leads to 10x worse performance. 

Conclusion: the statistical tests are useful, providing upper and lower bounds for achievable detection limits. Their false alarm rates are not affected by fluctuating aberrations within the linear regime. 
</details>

<details>
<summary> Lay 2004: Systematic errors in nulling interferometers </summary>

This paper presents an analysis of systematic errors to second order resulting from fluctuations in the null depth of a dual Bracewell interferometer. In other analyses, random noise consisting of photon shot noise and electronic noise has been studied. Systematic noise consists of 
- Imperfections in the amplitude response due to (static, linear) errors in the reflectivity of the optics.  Higher-order wave-front errors, tilt errors, focus, astigmatism, coma and beam shear all have non-linear responses which have to be modelled in detail. 
- Imperfection in the phase response of each arms, including path-length errors, dispersion and birefringence, which can be static (part of the equipment) or dynamic (drifts and vibrations). Fringe tracking corrects for this, but has its own imperfections.
- Errors in the positions of the collectors
- Errors in the polarisations, which will degrade the null depth.
- Fluctuations in the amount of stray light
- Thermally emitted photons
- Fluctuations in the detector gain

Expressions for the detected photon rate for each collector are formulated as $N=\int \int B_{sky}(s) R(s) P(s) d\Omega$ where B is the brightness distribution on the sky, P is a function representing the aperture and R is the response of the interferometer on the sky. The electric field in the single-mode spatial filter produced by a point source in the sky is given by the sum over the contributions from each collector: $r(s)=\sum_j  A_j e^{[i(\phi_j'+\frac{2 \pi}{\lambda} x_j \cdot s)]}$. The response $R(s)=r r^*=\sum_j \sum_k A_j A_k \cos \left( \Phi_j'-\Phi_k'+\frac{2 \pi}{\lambda} (x_j-x_k)\cdot s \right) \approx \sum_j \sum_k A_j A_k \cos \left( \Phi_j-\Phi_k+\frac{2 \pi}{\lambda} [x_{jk}\theta_x+y_{jk}\theta_y] \right) $, where the position of the collectors has been translated to a phase shift and where the small-angle approximation has been applied. The expression for $R(s)$ can be rewritten using $\cos(\alpha + \beta)=\cos \alpha \cos \beta - \sin \alpha \sin \beta$. Note that the first term is an even function and the second term is an odd function. When the sky brightness distribution is split up in a symmetric and asymmetric component, the full integral for N can be expressed as an even and an odd term, where $\cos(\Phi_j-\Phi_k)$ and $\sin(\Phi_j-\Phi_k)$ can be taken out of the first and second integral as they do not depend on solid angle $d\Omega=d \theta_x d \theta_y$.

Finally the total photon rate is the sum over all possible pairs of collectors, j&k.

Because of the properties of the sine and cosine, all terms with phase differences $\Phi_j-\Phi_k$ of integer multiples of half a period, $\pi$, adds up the cosine parts, but cancels the sine parts, so that only the symmetric part of the input signal remains. For all terms with phase differences $\Phi_j-\Phi_k$ of integer multiples of a quarter period, $\pi/2$, the cosine exactly cancels while the sine parts add up, so that the asymmetric part of the input signal remains. The 16 different permutations for a dual Bracewell are shown in Figure 2, together with the total output signal, for a signal of a star, local zodiacal dust, exozodiacal dust and a planet.

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/Lay2.png" width="450" height="600">


Signals of different astronomical sources:
- Planet: The planet has both a symmetric and asymmetric component. [? planet and antiplanet thing?]
- Star: The star is $~10^7$ times brighter than the planet and dominates the signal. Assuming that the star is center symmetric, there is only a cosine term. If in addition the star is modelled as a uniform disk, the Fourier transform $B_{*jk}$ can be represented analytically by a Bessel function. 
	- The stellar photon rate can then be perturbed using a Taylor series to first order  ($f(x)|_0=\approx f(0) + f'(0) \cdot x$) to its variables $\delta A_j, \delta \Phi_j, \delta x_j, \delta y_j$,. But when the perturbation is centered on the minimum, both $f(0)=0$ and $f'(0)=0$ so these first order terms are small. 
	- The second order terms introduce additional $\frac{1}{2}\frac{\partial^2 N_*}{\partial a \partial b} \delta a \delta b$ terms. 
	- [?] Lay drops the perturbations in x and y? Or does he absorb them into A and phi?
	- One can then substitute the partial derivatives with sensitivity coefficients. As there is an analytical expression for $N^*$, there is an analytical expression for each sensitivity coefficient too.
	- No derivation is made beyond second order. If the coefficients $C$ start to become large or if the deviations from the minimum become large, the higher order perturbations may introduce large terms. Is it interesting to test this?
	- The effects of asymmetric star spots on the noise is investigated, which contributes to the signal of the first harmonic of rotation frequency, which can omitted by considering only higher order uneven harmonics.
- Exo-zodiacal dust: It extends beyond the null, but if its contribution is assumed to be center symmetric, it will introduce a constant signal which can be subtracted
	- The emission from dust is easily 350 times the flux of the planet. Asymmetric inhomogeneities in the disk can be confused with the planet signal. With long imaging baselines, corresponding to high angular resolution, one should be able to distinguish between the point-like planet signal and a more spatially extended exo-zodi variation. 
- Local zodiacal dust: It is a uniform incoherent foreground source of photons. Its incoherence means that there is no phase or baseline dependence. The amplitude sensitivity coefficient is computed to be half that of the star.

The star, exo and local zodiacal dust combine into the perturbation leakage, that can be distinguished into 

1) Geometric leakage, caused by the linear terms from to the star, LZ and EZ extending beyond the null of the interferometer. These can be separated from the planet signal in post-processing because of being linear.

2) Null floor leakage, dominated by the star, arising from non-linear second order terms, causing the minimum to be above zero. The contributions are unique for each baseline pair. The non-linearity of the errors means that they cannot be calibrated in post-processing, so them must be minimised in real time. 

[!] When one assumes uncorrelated amplitude and phase perturbations, then linear and cross terms are dropped because they average to zero, so the *mean photon rate* is $\langle N \rangle \approx N_0 + \sum_j [C_{AAjj}\langle \delta \alpha_j^2\rangle + C_{\Phi \Phi j j }\langle \delta \Phi_j^2 \rangle ]$.

**Rotation:**

As the array is rotated, the planet signal is modulated. It can be matched with a family of template functions representing a grid of possible planet positions. The signal can be studied in the time-domain as well as in the Fourier domain. Ten measurements per rotation are sufficient to sample signals up of the fifth harmonic of rotation frequency. The star contributes a direct current, while centrally symmetric exozodi contribute to even harmonics, so these can be removed from the template function creating the modified template function $\eta$.

The demodulated (original, removing the effect of rotation) output signal of the planet is given by $O_{planet}=\frac{1}{T} \int_0^T N_{planet} \eta dt = \widehat{N_{planet}}$, so the standard deviation of the planet photon rate, but the total output also contains random and systematic errors, where the random error is $\frac{\langle N \rangle}{T}$, the square root of the total photon rate, which can include polarisation errors, stray light and thermal photon rates, and the systematic errors are related to the perturbation terms derived before, which, [!] assuming independence in the perturbations [!] is the sum over the individual variances. The constants are noise-coupling terms which quantify how much each term mimics the planet signal. $\langle \widehat{\delta a_j}^2\rangle$ are linear terms, averaged before the square, $\langle \widehat{\delta a_j \delta \Phi_k}^2\rangle$ are bilinear terms, and $\langle \widehat{\delta a_j^2}^2\rangle$ the average square value, squared. The perturbations are written as Fourier series to allow their analysis in the frequency domain, where their shape is determined by the power spectral density function. The amplitude fluctuations can be filtered out.

The bilinear cases with contributions from both an amplitude perturbation and a phase perturbation is more complicated because the power spectrum is a convolution of the amplitude and phase power spectra. For quadratic terms in amplitude, phase and polarsation, noise perturbations of one harmonic contribute other harmonics, also called frequency mixing. These cannot be removed in post-processing and therefore need to be limited in real time. 

It is systematic errors, in particular the amplitude-phase cross term, that drive the system requirements. Although stray light, thermal emission and changes in detector gain are not dominant sources of noise, their slow variations can mimic the planet signal even in small quantities. It is possible to remove them via phase chopping. 


**Phase chopping**

The output is split in two beams with opposite phases. The symmetric component is the same for both while the antisymmetric is the opposite. Subtracting one from the other, only the antisymmetric term $O_{asym}$ remains, which contains the planet signal, and errors at the chop frequency $O_{chop}$. The planet signal modulates differently between the left and right chop, but the chopped signal is similar for both signals. The variance of the systematic errors is the variance of the asymmetric signal and the chop signal. 

$\langle O_{asym} ^2 \rangle$

The variance of the asymmetric signal is correlated and can be redefined as independent variables containing the part that is common in both chop states (+) and the part that is different (-), for all combinations of perturbations up to second order. An ideal chopper has only (+) terms, whirch results in sensitivities to only first-order phase effects and second order amplitude-phase cross term. 

$\langle O_{chop}^2 \rangle$
The DC component cancels while the perturbations remain. One can expand the terms to arrive to an expression with all explicit noise-couplking terms. This resutls in equations that specify how to calculate the error at the chop frequency given the power spectra of the perturbatoins. 

An important caviat of phase chopping is that the amplitude-phase cross terms remain, which were the dominant error. It can be shown mathematically that it is impossible to remove the remaining systematic noise without removing the planet signal.

**Combining spectral channels**

The different channels can be added upt to obtain an overall signal-to-noise ratio. The systematic errors of adjacent wavelength bins are highly correlated. The different errors have different wavelength dependencies: path-length errors $\propto \lambda^{-1}$ in phase, wavefront error $\propto \lambda^{-2}$ in amplitude, , amplitude-phase  error $\propto \lambda^\{-3}$, shear error $\propto \lambda^0$ in amplitude and coating disperion $\propto \lambda^0$ in phase. These would require detailed modelling. Instead, a workaround is to define a number of independent channels and base the total SNR on the square-root of the sum of the squared SNR of the independent channels. 

**Wavelength dependence**
The planet signal risres towards long wavelengths. The systematic noise increases rapidly at short wavelengths. The random noise is higher at short wavelengths due to star light and at long wavelengths due to dust. As a result, spectroscopy becomes challenging below 10 um and has a maximum sensitivity at 14 um. 
<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/Lay10.png" width="500" height="450">

Moreover, as the target is closer to us, the angular size of the star increases beyond the null so that setllar leakage increases, so that the co-phasing error becomes the dominant limitation on the SNR. 

(Idea: colour code terms through thesis)
</details>

<details>
<summary> Bracewell 1978: Detecting nonsolar planets by spinning infrared interferometer </summary>

Written before discovery of exoplanets. Bracewell emphasises potential of infrared to study exoplanets, because in the visible, a planet emits $\sim 10^{-9}$ the flux of starlight, but in the Rayleigh-Jeans regime of IR, the factor would be $\sim 10^{-4}$ which is a *relative difference of* $10^5$ *between visible and infrared*. Bracewell proposes placing an interference null on the star. The star+planet signal can be disentangled by *rotating the interferometer*, which will modulate the (asymmetric) planet signal but will keep the stellar signal constant. Effects of pointing errors can be mitigated by designing the architecture such that a planet is $k$ fringe spacings away, so that the star signal modulation from misalignment has a different period than the planet modulation. Photon noise mainly originating from zodiacal light is pointed out as the main challenge when studying planetary photons. At the time of writing, not enough was known about the infrared environment in Earth's atmosphere, the solar system, as well as galactic and extragalactic contributions, to be sure that exoplanets could be measured over other incoming environmental infrared radiation.

</details>



<details>
<summary> Huber et al. 2024: Analytical and Numerical Instrumental Noise for LIFE </summary>
For LIFE's anticipated mission featuring four formation flying detectors, which aims assess the habitability of 20-50 temperate terrestrial exoplanets within ~20 pc., understanding all noise sources is necessary. In this paper we analyse two methods:

1) Numerical solution obtained by propagating instability-induced errors in Monte Carlo simulations

2) Analytical second order approximation of leakage from instrumental instability.

Based on flexibility towards different architectures, fidelity in modelling correlations of errors and fewer limitations concerning the variety of sources, the numerical methods is preferred.

All previous studies for LIFE described the interferometer using a mathematical framework known as a beam combination transfer matrix as proposed by Guyon et al. 2013. Hansen et al. 2023 used a Monte Carlo simulation to model reflectance and phase shift errors for a 5-input kernel nulling beam combiner. Lay 2004 provided an analytical framework using a second order approximation initially developed for the Terrestrial Planet Finder Interferometer, which was applied to the setup of LIFE by Dannert et al. 2024. In parallel this framework was applied by Laugier t al. 2023 on the ground-based VLTI-NOTT instrument. Dannert et al. 2024 proposed the 'numerical method', a matrix-based framework in combination with a Monte Carlo simulation. 

The measured number of photons at the output of an interferometer at a wavelength bin is an integral including the width of that bin, the sky brightness distribution and the instrument response. If needed, information on the two polarisation states can be passed on via a two component vector. The instrument response depends on the architecture of the k-input j-output interferometer which is represented by a $j \times k$ beam combination transfer matrix. It is key for interferometry to represent the incoming signals at the inputs as a vector of complex amplitudes. After multiplication with M, one obtains the outputs for each channel j, whose square modulus represents the instrument response $R_j$. Different combinations of output channels can be summed or subtracted to build differential maps. The photon rate can be written as a function of the amplitudes, phase shifts, polarisation angles and collector positions. The effect of instrumental perturbations on the photon rate $\delta N_j$ can then by described by perturbing the arguments of function $f_j$ via $\delta A_k, \delta \Phi_jk, \delta \theta_k, \delta x_k, \delta y_k$. 

The numerical simulation of Huber et al. 2024 takes the following steps: Each instrument perturbation term is represented by a random time series generated based on a specific power spectrum. Based on these inputs, the perturbed instrument response is calculated. This is multiplied with the true signal so the photon rate per bin can be calculated. This is repeated for every timestep, so the result is a 'heatmap' showing the intensity at each wavelength for each point in time. 

Benefits of the numerical method:
- Different architectures are possible, whereas not all architectures have an analytical framework ready for performance evaluations
- There are numerical methods that can investigate correlations between different wavelength channels from the simulations. Lay mentions the existence of correlations in the systematic noise between channels but does not provide an analytical framework for them. 
- Allows to investigate Lay's assumption that the instrument perturbations have zero mean over the ensemble average: $Var(\delta A_k)=\langle\delta A_k^2 \rangle$. This assumption helps to simplify the analytical expression, but its not quantified for which level of perturbation this 'small instrument perturbation regime' is valid

The analytical method by Lay 2004 writes out a perturbation of the analytical equation up to second order. Several terms are neglected and the remaining terms are used for an expression of the photon rate perturbation. Direct current terms and even harmonics can be taken out by a Fourier decomposition. The remaining expression for the variance of the signal at an individual output $j$ contains several noise-coupling terms which express the extend to which the specific instrumental perturbation mimics a planet signal. A method called phase chopping is applied, where the induced phase shift is swapped sign, which allows for the removal of the contribution of stray light, thermal emission and detector gain. Each remaining noise coupling term can be written as a Fourier transform of the best fitting planet template signal and the power spectral distribution of that noise term, which need to be defined by the user. This then feeds into signal to noise calculations.

Benefits of analytical method:
- 'Curse of dimensionality' of numerical integration can be omitted (proper sampling for many dimensions requires a lot of computational time), so the analytical method allows for more simulated observations at fewer computational costs.
- Provides framework for detection statistics and noise distributions, so that meaningful confidence intervals and false-positive-fractions can be defined, whereas a numerical model is uninvertible. 





</details>

<details>
<summary> Dannert et al. 2022: LIFE paper 2: Signal simulation, extraction and exoplanet parameters from single epoch observations </summary>

This paper presents LIFEsim, a tool to simulate observations for the Large Interferometer For Exoplanets, with the aim to quantify the mission's science potential. This version of LIFEsim includes fundamental noise, and has the option to include instrumental noise in the future. From the fundamental noise levels, already some limits can be placed on the permitted instrumental noise terms. Based on the current version of LIFEsim applied on a Monte-Carlo simulation of a mock survey, the signal to noise ratio of the spectra is dominated by photon noise.


Many of the existing and upcoming missions targeted at characterising exoplanet atmospheres are based on transiting exoplanets (e.g. occultation spectroscopy) but this limits the characterization space to close-in and large planets. Direct detection missions in the optical and infrared regimes would allow for a much larger number of detections omitting this bias. The idea of the life mission consists of a detection phase followed by a more detailed characterisation phase for objects of interest. This paper focusses on the first observation phase. Although detection is the main goal of this phase, it is also crucial to measure the radius, effective temperature and separation from the host in order to prioritise objects for the following phase.



**Concept of Nulling Interferometry**

LIFEsim is based on nulling interferometry using a double Bracewell setup. 
First, the paper illustrates the main idea via a single Bracewell: two collector apertures separated by baseline b can be combined to produce a sinusoidal fringe pattern. If the beams are offset by a $\pi$ phase difference, all on-axis light is cancelled, which in the context of exoplanets detection allows to cancel the light of the host star while keeping the light emitted by the off-axis exoplanets. The spatial resolution of an interferometer with a large baseline b is fundamentally better than that of a single aperture telescope of size D, with $\lambda/2b$ compared to $\lambda/D$. Coverage of the uv-plane is maximised by rotating the array, so that the exoplanet signal modulates as a function of angle as it 'moves' over the fringes in the transmission map. Repeating this observation at different wavelengths gives another type of extra information. A single Bracewell has some limitations: the transmission pattern is symmetric, so there is an ambiguity of 180 &deg to the position angle of the exoplanet. Moreover, it has a steep central minimum ('null'), which makes it sensitive to leakage from star light as soon as the star becomes a resolved disk.

**X-array architecture**

The LIFE mission proposes an X-array architecture with a dual chopped Bracewell combiner, where a 6:1 baseline has been found to best suited for post-processing techniques to remove instrumental 'instability noise'. A fifth spacecraft would collect and combine the beams. 


Mathematically, a the signal of the interferometer can be described via an input amplitude vector V representing the signal collected by each aperture: $V_k=r_k e^{i 2 \pi (x_k \alpha + y_k \beta) / \lambda}$. The architecture of the combination of beams by the interferometer can be described by matrix U so that the output amplitude vector is $W=UV$. For a dual chopped Bracewell combiner, the matrix U is


$ U = \frac{1}{\sqrt{4}} \left( \begin{array}
0 & 0 & \sqrt{2} & \sqrt{2}\\
\sqrt{2} & \sqrt{2} & 0 & 0 \\
1 & -1 & -e^{i \pi /2} &e^{i \pi /2}\\
1 & -1 & e^{i \pi /2} &-e^{i \pi /2}\\
\end{array} \right)$

$
U = \frac{1}{\sqrt{4}}$

$ \left( \right) 
$ 

$\begin{pmatrix}
0 & 0 & \sqrt{2} & \sqrt{2} \\ 
\sqrt{2} & \sqrt{2} & 0 & 0 \\ 
1 & -1 & -e^{i \pi /2} & e^{i \pi /2} \\ 
1 & -1 & e^{i \pi /2} & -e^{i \pi /2} \\ 
\end{pmatrix} $

```math
\begin{bmatrix}X\\Y\end{bmatrix}
```

The top two rows of the matrix denote the combination of the apertures constructively into two single Bracewell combiners without phase delays. The bottom two rows represent the 1&2 and 3&4 undergoing an $e^{i \pi}$ phase difference (=\times -1), making their interference destructive. With the goal of phase chopping, one of the destructive outputs receives an additional phase shift of $\pi / 2$  which reduces the susceptibility to instrumental noise effects. 

One can then construct 2D transmission maps $T_m$ relating the signal from sky position ($\alpha,\beta$) collected by the k apertures to the output m. The paper includes an analytical example (some questions about this). The transmission map shows the sensitivity of the instrument to incoming photons from different locations, with zero sensitivity in the center of the transmission map. Output 3 and 4 can be subtracted from one another to create a differential transmission map $T_dif$. The symmetry of this map means that any point symmetric source will not be transmitted to this differential map, although its statistical shot noise will be. The analytical expression for the transmission map can be used to illustrate how the light from different sources modulates as a function of rotation angle of this transmission map: an off-axis source modulates as it 'moves' over positive and negative fringes, while a spherically symmetric source gives a constant contribution independent of angle, which can be subtracted from the signal. The sensitivity of the transmission map can be summarised via the modulation efficiency $\xi$, which is the root-mean-square of the differential map over a full rotation of $\phi$. The angular distance from the center $\theta$ where the transmission map is most sensitive ($\theta_{\xi_{max}}$) depends on both baseline and wavelength via $\theta_{\xi_{max}}=0.59 \frac{\lambda}{b}$. 

**Astrophysical sources**
- True signal of the star, modelled as black body radiation, neglecting limb darkening i.e. assuming constant brightness across the disk
- True signal of the planet, modelled as black body radiation
- Geometric stellar leakage: not all stellar flux is suppressed when the angular extent of the disk of the star is wider than the 'null' of the interferometer
- Local zodiacal dust: Dust inside our own solar system adds radiation background via scattering of visible sunlight (not relevant for a MIR instrument) and thermal emission. Earlier work on DarwinSIM included a model of spectral surface brightness. Local zodiacal light is *diffuse* ([?] which means it comes from different directions, with different phases and wavelengths?), so it cannot be brought to destructive interference. The most effective way of minimizing this contribution is by observing away from local dust.
- Exozodiacal dust: Generally modelled as optically thin, symmetric with a power law surface density (or homogeneous? p6). The temperature of the dust depends on distance from the central star. Most of the radiation originates from the (hottest) central region. High-luminosity stars are expected to have larger disks. Larger disks have larger surface areas, thus larger fluxes. A symmetric disk can be filtered out, but its light does contribute to shot noise.
- Instrumental noise 
	- Instrumental noise terms exists of intensity variations A and optical path difference errors (related to the phase response $\phi$), and uncertainties in polarizaion rotation $\theta$ and space craft positions $x$ and $y$, which introduce additional stellar leakages and instability noise. These terms are systematic perturbations which can resemble a signal. 
	- In addition, detector dark current $I_D$ and thermal background noise deteriorate the measurement via photon noise. The total instrumental noise is $\sigma_{inst}=\sqrt{\sigma^2_{s,inst} + \sigma^2_{p,inst}}$, and the *fundamental noise limited case* is defined the regime where astrophysical noise is larger than instrumental noise.
	- The errors on the interferometric response, taht is detected photon rate $N$, can be found from the sum of responses of te individual baselines. The sensitivity of photon rate N to perturbations in amplitude and phase can be written as a 2nd order taylor expansion containing 2nd order terms and amplitude-phase cross-terms. Assuming a symmetric noise source, the asymmetric part can be set to zero so the remaining photon rate is $N*=\sum_j \sum_k Aj Ak cos(\phi_j - \phi_k)B_{*,jk}$. The first-order amplitude and phase perturbation can be derived from that, followed by second-order terms and cross terms, which together formulate an expression for the simplified Taylor series that constitutes $\delta N$. The photon noise of the instrument is then $\sigma_p = \sqrt{N+<\delta N>}. 
	- *Something about phase chopping: Up to second order, the method of phase-chopping is able to remove all systematic noise contributions except for one first order phase* $\delta \phi_n$ *and one amplitude-phase cross term* $\delta a_n \delta \phi_m$.

Comparison of the fundamental noise sources and instrumental noise, the short wavelength section is dominated by astronomical noise from the host star, the low wavelength section is dominated by local zodiacal dust. Systematic noise is the largest at shorter wavelengths, which becomes relevant when considering that the fundamental noise limit must hold for all wavelength bins. Regarding photon noise sources, both contributions are modelled as constant, where the requirements on the dark current seem to be one order of magnitude more stringent than what has been achieved so far.

For a closer distance to the target or hotter stars, the phase error requirement becomes stricter while the amplitude error can be less strict. An analysis of these effects for different types of stars is needed to investigate systematic biases.

An further analysis based on experiments could inform a more detailed model of the shape of the perturbation spectra with a component level breakdown. Current experiments are able to achieve the necessary amplitude error constraint but the OPD error is one order of magnitude larger than required. 

**SNR calculations**
   
Although the 'signal' part of symmetrical sources can be filtered using incoherent combination, the noise part remains and contributes to the statistical noise. The ratio of the signal of the planet over the total noise at a specific wavelengths bin, at output 3 or 4, is then 

$S
NR_\lambda = 
\frac{
 \int \sqrt{ \< S_{p}^2( \lambda )} \> 
 d \lambda 
}
{
 \sqrt{
 2 \int \left( S_{sym,3}(\lambda) + 
\sqrt{\< S_{p,3}^2(\lambda) \> \right)}
 d \lambda} } 
$

The total SNR assuming uncorrelated noise between the bins, is then the square-root of the sum of squares of the SNR of each wavelength bin, which scales with the square-root of the integration time, area and detection 

efficiency. Something with numerical discretisation errors? (p.6).

The SNR is baseline-dependent and can be optimised, although there is a trade-off between the amount of received planetary signal and the amount of stellar leakage. A baseline must be chosen such that the SNR is sufficient for all wavelengths of interest. A case study show that the decrease in planet signal at wavelengths shorter than 8 um with increasing shot noise causes a drop in SNR at short wavelengths. 



**Signal extraction**
The physical characteristics of the planet are extracted from the data using a maximum-likelihood method for single planets and for multi planet systems, starting with the most likely planet position $ \hat \theta_{p} $ within the domain of positive fluxes. The robustness of this method is investigated using a Monte Carlo simulation. The noise is assumed to be spectrally independent normal noise whose variance can be estimated from the data. Priors can be introduced. The SNR at the position of the planet is the flux over its standard deviation $SNR_\lambda(\theta_{p})=$,$SNR \left[ \hat F_{p} ( \theta , \lambda ) \right]=$ , $\frac{\hat F_{p} (\theta,\lambda)}{\sigma(\hat{F}_p(\theta,\lambda)}$. It is proposed that the combined cost functions of different wavelengths p(J') follow a chi-squared distribution with as many degrees of freedom as the number of wavelength bins. This is simulated and confirmed. 

Based on a chi-squared distribution, one can define the probability of a false alarm for a 5 sigma confidence level, which leads to a detection threshold of $\eta=87$.  It is better to use $p(J_\lambda'')$ which has one degree of freedom. This can then be combined for all wavelength bins. This leads to a detection threshold of $\eta=65$.

**Signal analysis**
The signal extraction method is applied to an Earth twin at 100 mas and a predicted SNR of 9.7 is found. The model is run many times to obtain representative average and standard deviation: 99.7 $\pm$ 1.5 mas and always a correct azimuthal position. The measurements of the individual wavelength bins lead to a spectrum with correct measurements but with relatively large 
uncertainty bars, especially for low wavelengths. The experiment is repeated for multi-planet systems, where the brightest planet is found first, subtracted from the signal, and the next one is found after. From the spectral shape, the temperature and radius are derived, with more uncertainty in the radius than in the temperature, with the best precision for the innermost (brightest) planet. It is possible to develop a maximum likelihood method to estimate the multiple planet positions simultaneously, which reduces error propagation. 

When applied to a simulation of nearby exoplanets, 98% of the planets were detected with separations within 15% accuracy and angles within $10^\circ$. Better SNR means better estimation of the parameters. Temperature and radius estimates are (anti)correlated. 

Machine learning based signal extraction algorithms still need to be tested. 




</details>

<details>
<summary> Angel and Woolf 1996: An imaging nulling interferometer to study extrasolar planets </summary>
Destructive interference seemed to have a trade-off with high-resolution imaging, but this paper reconciles this via a configuration with a *broad interference null* combined with a cross-correlation technique. When pointed away from the sun, the sensitivity is claimed to be only limited by the photon flux of the planet itself for a space-based interferometer with 50m baselines, 1 m telescopes and 10 hours exposure times.

At the time of writing, four exoplanets had been found. The atmospheric composition is best studied through infrared spectroscopy with CO2 at 15um, water at 7-8 um and ozone at 9.7 um. Ozone is a proxy for the presence of photosynthesis. To resolve planets with a conventional telescope, a mirror diameter of 60 would be needed. Two problems with conventional interferometry are the need for low transmission across the whole stellar disk, and the confusion of dust with planets. The solar system emits 10 um flux 300 times stronger than Earth from about the same orbital distance.

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/AngelWoolf.png" width="800" height="250">

The proposed arrangement has two superimposed Bracewell interferometers with 180$^\circ$ achromatic phase difference and different spacings. The outer pair has an amplitude of almost half of the inner pair, such that the summed amplitude has 0th, 1st and 2nd order zero. The intensity now depends on $\alpha^6$ near the axis, which leads to an exceptionally broad as well as deeper central minimum. The setup is envisioned with a fixed geometry with path lengths adjusted using fixed trombone delay elements. Changes in element spacing are inadvisable because of the risk of mechanical failure.  

Like Bracewell, the interferometer world be rotated about the line of sight. As the physics is dependent on wavelength, it will give extra information when the infrared continuum is dispersed into separate narrow-band channels. The images can be recovered by a cross-correlation method, using a Fourier transform, where a summation is taken over rotation angle and wavelength. Off-centre point sources appear twice upon rotation of the interferometer.

Calculations are done for an architecture of four 1m telescopes with a total length of 50m with a synthetic signal of four planets. All planets could be retrieved. The stronger but smooth modulation expected from an elliptical dust cloud would be expected to be straight forward to filter out because it modulates with twice the rotation frequency. Local zodiacal dust, read out noise and dark leakage were assumed negligible. 

A spectrum is modelled as well. An atmosphere with H20, O3 and CO2 deviates from a black body spectrum via absorption lines. To significantly detect them, 20 SNR is needed. , which would require 80 days of observation. If zodiacal dust would be 10 times brighter, we would need 10 times larger telescope area or 10 times longer iterations for the same SNR. 

The double bracewell outperforms a configuration with diamond as described by Angel (1990), which has transmission varying as $\alpha^4$ near the null. Leger et. al proposed a geometry with five elements in a circle which should allow to distinguish a planet from an elliptical dust cloud close to the star, but it would not allow to distinguish multiple planets from each other. The double bracewell allows for baselines beyond 50 meters so that 1m telescopes are sufficient provided that dust signals are low.
</details>

<details>
<summary> Laugier et al. 2020: Kernel nullers for an arbitrary number of apertures </summary>

Kernel-nuller architectures for nulling interferometry can be designed to produce observables that are inherently robust to tiny optical path differences between apertures. In addition to a four-beam combiner, this paper presents properties for arbitrary numbers of apertures. The combiner matrices of all-in-one kernel nullers grow in complexity for larger numbers of outputs, but still seem to be the most efficient in comparison to multiplexing nullers. Multiplexing nullers work independently over a smaller number of sub-apertures and could conserve a higher throughput per output, but have fewer robust observables and lose some of the off-axis light leading to a lower overall efficiency. 


Nulling interferometry has a great potential for imaging exoplanets but the vulnerability to optical path difference errors is an obstacle. Some designs allow for *linear combinations* of outputs that are robust to aberrations to second order, which are called kernel nulls.

The design by Martinache & Ireland is a four-beam nuller described by a real matrix. The output intensities are a *degenerate* function of the target information and input phase perturbations. 

If the outputs are fed to a second stage which introduces $\pi/2$ phase shifts to the nulled outputs, then the combiner matrix is complex, and the components of the output electric field span the complex plane. This makes it possible to disentangle the effects that environmental perturbations have on the input electric fields. There are more outputs than inputs because there are more ways to combine the different input fields. Pairs of outputs produce the same response to environmental effects while producing different response to off-axis light. 

Electric field combinations which are mirror images of one another are called enantiomorph. The outputs are insensitive to global phase shift and therefore the *phase reference input*[?] can be arbitrarily phase shifted to coincide with the real axis. 

It can then be shown, when the phasors are small ($z_k=e^{-j \phi_k}\approx 1-j\phi_k$), that  that the conjugate pair of null rows $m_1$ and $m_2$ of $M$ always produces a kernel null that is robust to these arbitrary imaginary phasors. The same is true for purely real input electric fields that would correspond to pure photometric error generated by fluctuations of the coupling efficiencies. 

[To do: understand figure 3]


It can also be shown that two input electric field vectors $z$ and $z'$ from sources located at symmetric positions in the field of view cancel out in the kernel nulls, by substituting $z'=z*$ in the expression for the kernel null $k(z)$.

The complementary observable to the kernel null is $\tau(z)$, the sum of the two outputs. It does not have the same robustness to aberrations, but it may be used to study brighter symmetrical features and encodes some information about the input phase errors. Therefore, $\tau$ may be useful to locate the minimum OPD error.  

The maximum number of nulled outputs for a given number of apertures is $n_{max}=(n_a-1)!$, although some of them are a linear combinations of the response function of other nulls. The number of independent kernel nulls is always the same of the number of independent closure phases: $\frac{(n_a-1)(n_a-2)}{2}}$. A complete nuller provides the maximum number $n_{kn}$ of independent observables.

The complex coefficients of $M$ can generally be written as $M_{l,k}=a_l \cdot e^{j\phi_{k,l}}$ where $a_l$ is a real coefficient that normalises the the column vectors. 

Our goal may be to construct complete combiners using the minimum number of combinations from the full matrix M, with the intent of increasing its throughput. Some words about ncrops and finding these combiner matrixes...?

The most simple example is the three-input nuller with two enantiomorph nulled outputs that together produce one robust observable. Another example is for a six-input nuller whose combiner matrix has 121 rows, but which an be cropped to 21 rows, which offers 20 independent nulls and 10 independent kernel nulls. Each kernel null corresponds to a map that covers the field of view differently and brings new constraints on the properties of the target.

The robustness to phase-noise is studied using a simulation using Gaussian, non-correlated phase noise with the same rms value for all inputs. The standard deviation of the output intensities of the nulls is analysed. Larger kernel-nulling combiners provide a rejection of the phase noise that is very similar to the smaller combiners. Comparing the raw nulls to the kernel nulls, shows the expected suppression of second-order phase noise provided that the input phase errors are small enough for the small phase approximation to hold. 

While the distribution of kernel nulls is close to Gaussian, the distribution of null *intensities* is not. A full performance comparison would include the coupled effects of phase and amplitude fluctuations such as formulated in Lay 2004.
</details>

<details>
<summary> Flasseur et al. 2024: Shrinkage MMSE estimators of covariances beyond the zero-mean and stationary variance assumptions </summary>

Estimating the numbers of a covariance matrix of a given dataset is challenging for small sample sizes because the variance on the estimators is high. Shrinkage methods can force regression coefficients (such as the components of the covariance matrix) towards lower values via a biased regularization. As long as the reduction in the error is large compared to the introduced bias, the resulting estimate will be more accurate. There is a trade-off between reducing the variance and introducing the bias. 

This paper extends existing shrinkage methods for a multi-variate Gaussian case to include a regularisation matrix accounting for sample variances (non-uniform diagonal values) for both centered an non-centered (mean is non-zero) samples. 

To optimise the trade-off between reducing the variance and introducing the bias, a coefficient $\rho$ can be defined, which serves as weighting between the original unbiases-high-variance estimate $\hat S$ and a low-variance biased estimate $\hat F$, which together combine into a new estimator $\hat C$, such that it minimises the error between the a *best estimate* of the covariance matrix $\hat C$ and the true value $C$. 

If the true covariance is known, one can compute the *oracle shrinkage estimator* that optimises this weighting, but in practice the true $C$ is generally unknown and the goal is to find a best estimate together with a value for the uncertainty. However, it is possible to iteratively find a near-optimal *oracle-approximating shrinkage* through iteration of information inherent to the sample. 

The paper presents analytical expressions for the *oracle-approximating shrinkage* estimator in the case of 

1) Centered samples + non-weighted sample covariances + a diagonal structure for F

2) Centered samples + regularization matrix $\hat F$ that accounts for sample variances

3) General uncentered case (unknown $\mu \neq 0$ ) + weighted sample covariances

The performance of the method is validated using a numerical simulation.

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/MMSE1.png" width="800" height="370">

The performance of the method is then compared to the Ledoit-Wolf and Rao-Blackwell-Ledoit-Wolf shrinkage methods for the case of i.d.d. regularisation matrices $\hat F$, applied datasets with 3 levels of underlying covariances $C$, and compared to the results for the average variance solely (case 1) and sample variances (case 2). 

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/MMSE2.png" width="700" height="660">

Finally the results are compared to high-contrast imaging data of VLT-SHPERE where stellar leakage causes spatially correlated speckles. Stellar speckles are quasi-static upon rotation of the telescope, while a planet signal rotates in correspondence with Earth's relative rotation. A match-filter looks for an exoplanet. The mehtod is applied to a real dataset with N=92 esposures with $10^6$ pixels each, based on which P=113 circulare patches are chosen. In addition to the two real sources known to be around the star, 30 more are added. The model is forced to make 60 detection, of which thus 28 are false, so that the SNR for true and false detections can be compared for the different methods.
<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/MMSE4.png" width="800" height="600">

</details>


Anticipated papers:
- Photobombing

## Meetings 
- **Romain & Loes 3 May**
	- Start in July (tbc with Rose)
	- Start with literature study, first version around august
	- Then define perimeter and objective. 
	- General scope: Calculate effect of instrumental errors on interferometry measurements.
	- Collaboration of TU Delft + KU Leuven + ESA + AMOS + Liège to evaluate potential of mission similar to LIFE but without formation flying. 
	- HWO: In visible range, using coronography. As stars emit stronger in visible than IR, so stronger SNR needed: SNR 10^-10. NASA led, specific timeline.
	- LIFE: Not funded yet, but group of scientists interested and hoping for funding. In mid-IR good thing is that for Earth-like exoplanets you need SNR 10^-7. But longer wavelengths also require longer baselines.
	- Architecture 1: Double bracewell, based on Bracewell ‘detecting nonsolar planet’ paper. Allows for longer baseline. If your resolution is good enough to resolve the star as a disk rather than a point, you need to cancel the light of the whole disk. Double bracewell allows wider baseline and still cancel light of star.	If there is an error in Output 1, then the resulting errors in Output 3 and Output 4 can be correlated and taken out.
	- Architecture 2: Kernel-nuller (Martinache & Ireland 2018). More general concept of combinations of different outputs, of which double bracewell configuration is a special case. More possibilities to cancel out correlated noise. 
	- Existing work on noise: LIFESim from Philipp Huber in Zurich, including fundamental errors. Future ambition to add Lay aspect to it. Romain has paper with analytical approach for kernel nullers and one with more numerical Monte Carlo approach for chromatic errors.
	- Possible output of thesis: 
		- Library of errors. Helpful in the future (and would be helpful now) to make approach to calculate errors more standardised. Idea is to collect errors from different instrumental sources and allow a wider group of scientists to build simulations for different instrumental configurations. Idea to explore: get in touch with many different types of engineers e.g. at teams of AMOS, with support of Jérôme.
		- Second output: Impact of instrumental noise, possibly also correlated noise. 3rd order error is still there. We don’t have good model to express covariance. Maybe work out semi-analytical/semi-numerical estimate, build on Lian’s results.
	- Lian's work: Present practical benefits of Kernel architecture via more complete figure of merit than just 'yield'.


- **Romain & Loes 29 May**
	- Overview of Loes' statistics course: 
		- First half is on frequentist method, ordinary/weighted/total/non-linear least squares, robust regression, regularised regression, bootstrapping, model selection
		- Second half on Bayesian statistics: priors, likelihood, numerical methods including Monte Carlo simulations, Bayesian model comparison, Gaussian processes
		- Most relevant part of the course for the thesis will be the part about the frequentist method. The framework for covariance and uncertainties is useful, and in this field it is not common to insert priors.
	- Discussion of Lay 2004: Fundamental noise is statistical variation in light from planets, host star, background. Instrumental noise is noise due to drift, mirror (mis-)alignment and generally correlated. The errors in interferometry are generally not Gaussian, because the intensity of the fringe pattern varies with $cos^2$. Lay 2004 presents an *analytical approach* with as input error a phase and amplitude misalignment, and a 2nd order Taylor series as output. Some terms that are summed in the Taylor series are correlated. You want to obtain the covariance matrix. With this, you can make a whitening matrix to rescale and rotate covariance matrix. Challenge is that there are only few exposures. You want to know beforehand what the covariance matrix looks like. A technique similar to bootstrapping   is used. 
	- Numerical approach: Generate a sample with a distribution of errors. But it is hard to link this backwards to make statements about the necessary stability of input variables. You can make iterations,  but it is not rooted in physics. It *does* keep track of derivatives via the autodiff tool, so you can get the Jacobian of how 'Output Y' changes as a function of 'Input X'. Machine learning could also be explored. Goal is to get something 'smart-numerical'.
	- Low hanging fruit: extraction of signal with and without whitening
	- Error trifecta: Lay assumes that statistics are centred on zero. But what if you have a bias? In that case many symmetric terms of Lay's Taylor series don't cancel out, analytical mess.  

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/ErrorTrifacta.png" width="250" height="300">


- **Leuven/Delft/Liege consortium meeting 10 July**
	- Introduction by Jerôme Loic: In October, as part of task 2, project should deliver one or two preliminary designs with performance assessments, to be then explored further, while keeping existence of alternatives in mind. Mian designs: Linear (DCB-A), X-array, triangle, diamond shape (NR-array). 
	- Presentation literature review by Denis Defrère: Literature review to be presented to ESA. Best currently operating equivalent example of nulling interferometry on Earth is Keck-nuller, consisting of 2 telescopes with splitted mirrors, adjustable optical paths and phase modulators. Using integrated optics may be interesting because very compact, but may be too sophisticated for space. Integration time and length of the mission are not strongly constrained. The Large Binocular Telescope Interferometer (ground-based single-mount telescope connecting two 8-meter telescopes in Arizona) may also be interesting to look into.
	- Presentation of Ida Sigusch, PhD at TU Delft: Did literature review to see if choice of beam combiner and phase shifter influences choice of architecture. Defined metric weighing wavelength/null-depth, stability and throughput and maturity. Placed 9 achromatic phase shifter on grid of broadband vs. dispersive/narrow and fixed vs. arbitrary phase shift. Also looked at bulk optics versus integrated optics. Integrated optics as a high loss, birefringence and is not tested on cryotemperatures. Conclusion for now that architecture determines choice of beam combiner and phase shifter but not the other way around. 
	- Presentation on yield by Colin Dandumont (PhD at Université de Liège): Used p-pop from Kammerer fto simulate yield. Master student Rogier Norbruis made 6 configurations. All configurations have a similar order of yield. Discovery rate was studied, starts dropping after 4 years as most easily detectable exoplanets in the neighbourhood have then been found. Finding that *if* you have a longer baseline, then you can really benefit from larger aperture, but otherwise you improve within the same order of magnitude. Important to still take into account is the characterisation phase, which comes after the detection phase. Then is becomes relevant to be able to detect the same object in all wavelengths of interest.
	- Presentation on instrumental noise by Romain Laugier (Postdoc KU Leuven): Studied transmission maps of different architectures. Main take-away: sensitivity for longer and shorter wavelengths peaks at different angular separations. For characterisation, you need sensitivities of all wavelengths to be high enough at the separation of interest.
	- Presentation by Jasper Bouwmeester (Assistant Professor TU Delft): Trade-offs of architecture shapes. Performance criteria are detection of new exoplanets and characterisation of detected exoplanets. Main practical metric is complexity of deployment. Main shapes considered are linear, X-array and NR kernel nuller. Box schemes of criteria, some can be calculated scientifically while others are based on gut feeling. Discussion: Variable baseline would be interesting but practically challenging. Important criterion is ability to characterise. Linear design would already be fantastic. Kernel-nuller is promising but more risky. Way forward: model performances of different designs, then among best performances choose design that is easiest to build.




- **WG 3.2 meeting 11 July**
	- Jens Kammerer (topic: yield estimation), Daniel Angerman (project manager), Felix Dannert (PhD on instrumental noise), Angelica Kovacevic (extragalactic time astronomy), Philipp Huber (PhD on instrumentation) 
	- Presentation on PHIRINGE, photon count data generator for nulling interferometers. Takes config files based on astropy and scipy. Uncertainty in flux estimated from spread in (MC simulated?) data. Open question is to find the cleanest way to estimate uncertainty. 
	- SCIFYsim, initially for NOTT. Efforts made to include symbolic math, but so far not as useful as hoped. Instrumental noise is correlated, interesting to extract. Tip: old gravitational wave literature has analysis on the noise introduced by different tools. 
	- Discussion afterwards: Loes has been added to LIFE Slack. Most relevant groups are 3.2 LifeSIM, 3.3 signal procext, 4.1 nulling tech. Thesis could include some joint interphase with PHIRINGE, SCIFYsim and LIFEsim. Datastandard from NIFits group should provide 'one way of doing stuff'. Other idea: explain impact of Lay's assumption of optical alignment and correlated terms within sum for non-ideal case. With Monte Carlo you can only 'retrofit' the underlying distribution if you have strong confidence in the underlying model. Some new numerical tools could be explored, machine learning algorithms that do a symbolic search. 


- **Romain & Ida & Lian 11 July**
	- MC simulations take long. Idea to use Python functionality to dump outcomes into file and load from there.
	- Python code gives error on Windows because it doesn't have a /tmp folder as in Unix based systems.
	- Idea to have joint meeting with Delft & Leuven every other week? To be determined.
	- 3rd week of August, internal office moving in Delft should be done

- **Meeting Loes & Romain 12/7**
	- Question why linear setup is still interesting compared to 2D layout. Has to do with broader rejection, so it still works if your star is resolved as a disk. Key concept is spatial coherence, explained further in Angel & Woolf. Useful to review theory ondiffraction and point spread function of telescopes.
	- Does nulling interferometry also work for binaries, as some research concluded that the majority of stars are binary systems? Yes, as long as the other star is not in the field of view. Moreover, mostly heavy stars are binaries, and in the solar neighbourhood are many lighter stars.
	- Question about symmetric baselines: They are redundant if they are exactly the same, skew gives more results. Best layout to be studied further. Can we choose optimal baseline for specific molecular lines in the spectrum? Not needed, because planet generally (if not face-on) will move over different angular separations as projecten on the line of sight. 
	- Hint for JWST Kernel-phase paper: kernel-phase is not the same as kernel-null. It's more like a closure phase. Focus on statistics in that paper. Tutorial based on on Ceau does not include correlation and whitening yet. Main question: what to put in matrix that is representative?
- **Meeting Loes & Romain 7/8**
	- General concepts make sense but when light is described in complex form it's a bit harder to follow. Light is like a 'vibration' represented by a complex number, which cannot be measured in real time, but what can be measured is its square, the intensity, which is related to the probability to measure a photon, which in an experiment is measured from the number of photons. The diagrams with coloured arrows of Laugier et al. 2020 have black circles, which represent the complex amplitude. Circles were chosen because they scale with radius squared. 
	- You can represent multiple incoming rays as a vector z, where $z=(1,1,1,1)$ represents on-axis light, but generally different rays have different phases. Figure 3.1 represents small phase errors (within the small phase approximation) in multiple rays from the same bright on-axis source. The nulled outputs are not exactly zero, and their intensities appears in 3 pairs based on symmetries. When the pairs are subtracted from one another (kernel-nulling), the light cancels. Figure 3.2 represents small amplitudes of significantly different phases coming off-axis from a planet. The pairs of nulled outputs become decorrelated, so that if one nulled signal is subtracted from its counterpart, a signal remains, which is the goal. 
Figure 3.3 represents the null peak, the peak of the transmission map of the first kernel. The configuration there is such that output 1 constructive, output 2 is destructive, so that their difference (kernel-null) is maximal. Outputs 3&4 and 5&6 are equal and therefore cancel, so the other transmission maps are flat.
	- Discussion on transmission maps of Martinache & Ireland 2018. Transmission maps of nulled outputs go from 0 to some positive maximum. Different nulls outputs look like translations of one another. A phase shift of a nulled output is a rotated transmission map. It has some axis of symmetry. The differential maps, where one is subtracted from the other, are odd functions, with values from negative to positive.  It is the square modulus of output 1 minus output 2 etc that leads to $\sin^2 \cos^2$-like expressions as seen in LIFE paper 2. 
	- By taking the difference of two outputs you dump a part of the information. This is what's meant with $\tau(z)$ in Laugier et al. 2020. Mathematically, it has half of the information of the planet contained in a projection on the parameter space with (anti-)planets. They are like modes. Where the Zernike-Von Mittard theorem states that the Fourier transform of the intensity distribution function of a distant, incoherent source is equal to its complex visibility, in our case something analogous is true for not a FT but an eigenmode transform. $\tau(z)$ has half of the planet and anti-planet, plus half of the planet minus he anti-planet [?]. 
	- A key figure is Figure 8 from Martinache & Ireland 2018, which shows that the nullers have a non-Gaussian error in the intensity, but the but that the error in the difference between two channels (the kernel) is close to Gaussian. It is key to understand how much error there was in the first place. Possibly $\tau(z)$ can help us. In the paper of Hanot 2011 is an analysis of nulling self-calibration for a single Bracewell setup based on the distribution of the outputs. 
	- There is a notebook with analytical expressions from sympy on Romain's GitHub: kernuller. 
	- Way forward: still look at Lay 2004, SPIE paper and covariance shrinking. Look at GitHub. Goal of some written review is an overview of everything that we know now. 
- **Meeting LIFESim 8/8**
	- Introduction to Claudine, bachelor student studying error propagation of wavefront error into amplitude and phase of collected light
	- Presentation of progress on nifits by Romain: Instrumental errors data standard for nulling is making progress. Nifits contains model of instrument, IO loads and saves FITS files, backend is there for the math. Notebook was created that walks you through a simple file for a three-telescope differential nuller. You can model points in the field of view. Gives you the results of the outputs, the response maps and the intensity. Different grid shapes are possible. You can adjust the gathering power of your telescope. You can find it online: https://github.com/rlaugier/nifits. There is also the nifits documentation.
	- Update from Philipp on pipeline for signal extraction of synthetic single-epoch data. The goal was to include error bars based on the data itself, not based on prior knowledge. The steps are 1) Lay 2006-based signal extraction where the black and white lines are correlated instrumental noise. After subtracting correlated instrumental noise, it looks much more like a data that is fit for some detection test based on match filtering. T-test is used, which is based on Gaussian noise. MCMC used to estimate flux for the bins to obtain the spectrum with the error bars. Then a black body spectrum can be fitted to extract a temperature and radius. Discussion on further points of improvement: 
		- If *temporal* spectral covariance can be extracted, we could base errors on something less dependent on user input. 
		- Regarding T-test: These maps can be converted in a test statistic inherent to the data, like some log-likelihood. Idea to take inspiration from Ceau et al. 2019 and the notebook on signal detection: Energy detection test - tells you whether or not you have reason to believe that there is more than the star. The Neyman-Pearson test tells you the likelihood that you are detecting a planet of a specific position and spectrum (or T and r). For these tests, the shape of the ROC curve is known mathematically.
		- There is also the optimal binary test TB, when you know that the signature has a certain dictionary of possible shapes and you find the maximum likelihood within the family of signatures. Formulate a null hypothesis and detection hypothesis, based on e.g. less than 5% false alarm rate. The more parameters you allow to vary, the lesser the power of TB. The more you know about your target, the closer TB is to Neyman-Pearson. For TE and TNP, the shape is known mathematically, while for TB you may have to do some simulation. 
		- Concern that that it is computationally demanding to scale up the test procedure. But an a posteriori method may occult spatial correlation. Fitting a black body spectrum *within* the MC simulation should be more efficient and you wont loose the correlations. So you can pass on the covariance matrix to people that will analyse the spectrum.
		- A step with some simulations may allow to find the covariance matrix. Loes may look into analytical methods that can shrink the covariance matrix. 
- **Meeting Loes & Romain 15/8**
	- Questions on Lay 2004: What happens to dxi, dyi after equation 13? -> They are considered only to first order because they are small. The x and y change the transmission map but not the behaviour of the on-axis light. 
	- In the text after equation 29, Lay notes that some perturbations are straightforward to calibrate out. Is this the same concept as discussed later in the paper in the context of linear and non-linear noise-coupling terms, where some perturbations are a convolution of the amplitude and power spectra, causing frequency mixing, making calibration in postprocessing difficult? -> some terms can be subtracted via spectral analysis if the signal is known, Felix knows more, but the Lay paper is not the last word of the discipline so not needed to get stuck on it.
	- The assumption for equations 30 and 38 is that the amplitude and phase perturbations are uncorrelated between the collectors. Is this valid? -> It is not the first think to investigate. Correlation of $\Phi$ and $A$ could happen e.g. when a collector is titled. However, if you calculate the kernel-null $\kappa = K \cdot |M \cdot z|^{2}$ you end up with an expression where the phases only couple with amplitudes of *other* terms, so not for $\Phi_i A_i$, and you do not expect very large correlations between the cross-terms. The correlations that are much more important are 1) that of the spectral channels, realising that both amplitude and phase are functions of wavelength. 2) correlations due to combiners of multiple outputs, for which the optics did not yet exist at the time of writing. 3) temporal correlations: vibrations will be similar between one timeframe and the next, so that there is a temporal spectrum of input errors. If there is a way to translate a power spectral density into a temporal covariance matrix, that would be interesting. Helpful to look at Lay 2006 and Matsuo et al. 2023. Matsuo is based on fitting polynomials.
	- Question on covariance shrinkage: It is stated that OAS methods outperform the (RB)LW methods, but comparing the performance Figure 4 shows not really a difference? -> The important message is that this method is based on solid analytical arguments for approaching the best answer. The improvement in performance is incremental. In the field still very basic methods are used. 
	- On the y-axis of Figure 4, you could construct a test statistic, so it links with Ceau et al. 2019.
	- Next step is to start implementing stuff. Functions like whitening can be added to the NIFits package.






## Planning

Just a first rough outline as conversation starter! Anything can be adapted.

| Month | Activity | Hours | Comments |
| --- | --- | --- | --- |
| July & August | Exploration: Literature and research plan | 150h| Fulltime except last 2 weeks of July |
| September | Deep thinking: Extend theoretical concept | 150h| Fulltime |
| October-December | Implementation: Create new tool or integrate existing tools | 150h| ~2 days per week. Romain 5 weeks in Nice from November |
| January- mid-February| TU Delft-based applications | 150h | ~2 days per week, depends on exams|
| mid-February - mid-March| Writing and knowledge transfer | 150h | Fulltime |
| mid-March - end-March| Buffer | 150h| Fulltime |



Explore [150h] July & August

- Explore key questions in instrumental noise for nulling interferometry
- Explore which theoretical concepts are exciting and realistic to further work out in the scope of the thesis
- Explore synergies between Leuven/Romain's work and Delft/Jerôme/Ida's work
- Write down literature and research plan in coherent document

Deep thinking [150h] September

- Work out theory further for analytical or somewhat analytical+numerical component (tbd)
- How to know the covariance matrix, given our knowledge on the instrument architecture and a small bit of observations. 

Implement something [150h] October-December

- Try out some new innovative approach, and/or
- Integrate existing elements in SCIFYsim, LIFEsim?

Apply / test [150h] January - mid-February

- Bridge to work in Delft and partners

Bring together [150h] mid-Feb - mid-March

- Wrap-up
- Write down
- Present / defend thesis
- Sustainable handover

Buffer [150h] rest of March

Next thesis deadline: mid-June


## Quick_Summary

**1. Exoplanet science**

Since 1995, the number of detected exoplanets has increased from zero to over 5000. Based on large scale surveys, it is now expected that almost all stars host at least one planet. [reference to be added]. Small planets are more common than gas giants. Earth-sized planets are expected to exist in the habitable zone of approximately 20% of sunlike stars and in 40% of M-dwarfs. In addition to the rocky planets, gas giants and ice giants and their observed orbits in our own solar system, new types of planets were detected named "super-Earths" ($1.5-2 R_E$) and "mini-Neptunes", as well as the presence of 'hot jupiters': gas giants in very close orbits, with super-Earths even being the most common type of planet. The insights have led to improved theories  on the accretion of planets from proto-planetary disks. Moreover, the discovery of 'hot jupiters' let to theories on orbital migration of planets. More detailed and more complete detections of planetary systems are needed to consolidate such theories.

An open question is whether rocky exoplanets in the habitable zone of stars can sustain life. Besides suitable temperatures for water to exist in its liquid form, estimated via the habitable zone which depends on the temperature of the star and the planet's orbit, such planets need an atmosphere, which in theory could be traced by molecules such as CO2 and CH4. Should life as we know it already exist on a planet, then oxygen should be present in the atmosphere, which can be traced through O3. The measurement of the chemical signature and more detailed physics of a planet is referred to as *characterisation*. As planets are very faint compared to their host star, direct characterisation has been very challenging so far and is one of the main aims of upcoming missions. 
It has been argued by previous studies that 30 to 50 exoplanets need to be studied to get statistically robust results on the fraction of habitable rocky HZ exoplanets. 


**2. Main detection methods for exoplanets**

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/ExoplanetArchive.png" width="1000" height="400">


- Transit photometry: Surveys like TESS (Transiting Exoplanet Survey Satellite) and Kepler have detected dips in the intensity of a star due to light blocked by planets while transiting in the line of sight of the star and the observer. The next generation instrument will ESA's PLATO mission (PLAnetary Transits and Oscillations of stars), scheduled for launch in 2026. The next step in studying transits, is to focus on the light of the star as it passes through or blocks the atmosphere of the planet. It is anticipated that the James Webb Space Telescope (JWST) will help to characterise exoplanet atmospheres through  transmission and emission spectroscopy (as well as high-contrast imaging). Moreover, ESA's ARIEL infrared spectroscopy instrument is designed for this.
- Radial velocity: Planets and the host star orbit around a common centre of mass. The spectrum observed from the host star is therefore slightly Doppler-shifted, forwards and backwards, in a periodic manner. The main past missions were HARPS (High Accuracy Radial Velocity Planet Searcher) HIRES  (High Resolution Echelle Spectrometer), and CORALIE. Missions currently in development with higher precision radial velocity instruments include CARMENES, NIRPS, ESPRESSO, MAROON-X, HARPS3 and EXPRES.
- Gravitational lensing: This technique benefits from the magnification of light from a far away location due to the gravitational field of a massive object in its line of sight. The main discoveries were by OGLE (Optical Gravitational Lensing Experiment).
- Direct imaging: Detection of light from the planet itself, either by reflection of the (visible) starlight by the planet's albedo, or by recording the approximately black-body thermal emission from the planet itself. This requires very advanced instruments, such as the Hubble Space Telescope ACS and NICMOS instruments, the VLT's Spectro-Polarimetric High-contrast Exoplanet REsearch (SPHERE) and Keck's Near-Infrared Camera (NIRC2). The upcoming ELT will contain next level direct imaging possibilities through its instruments METIS, PCS and HIRES. NASA has proposed the Habitable Exoplanet Observatory and Large UV/Optical/IR Surveyor LUVOIR missions as *space-based* direct imaging instruments focussing mainly on coronography based techniques to detect visible starlight of sun-like stars reflected from the surfaces of exoplanets. On the other hand, efforts are made among a predominantly Europe-based informal LIFE consortium to include a formation-flying four-aperture infrared nulling-interferometer in the next ESA program. ESA is also exploring the feasibility of more technologically ready alternative space-based nulling interferometer with four apertures mounted at fixed positions on a single spacecraft. A general framework for instrumental noise for a family of architectures of space-based nulling interferometers is the focus of this thesis.


**3. Assets of infrared nulling interferometry**

Despite the ambitious planned missions using a range of exoplanet observation methods, there are several additional aspects that a space-based mid-infrared nulling interferometry mission would bring to the table.
- Direct imaging puts more direct constraints on temperature and size of the objects than indirect methods, which are based on several assumptions.
- Higher spatial resolution than single telescopes thanks to the dependence on baseline rather than aperture size. Compared to interferometry in the visible spectrum, longer baselines are needed, but  the necessary baselines to resolve planets around nearby stars are still feasible. 
- The infrared spectrum includes molecular absorption lines revealing information about planetary atmospheres, including biosignatures.
- More favourable contrast ratio than interferometry in the visible light: The thermal emission of exoplanets lies in infrared wavelengths, while the light of stars peaks in the visible spectrum (depending on the type). Therefore the required contrast of an instrument to disentangle the signature of the planet from that of the host star is more favourable in the infrared $c\approx 10^6$ than in the visible $c\approx 10^{10}$. 



**4. Main design and predicted yield**

The main design of the LIFE space mission consists of 4 formation flying infrared telescopes, with apertures to be determined between 2 and 3.5m in a 6:1 X-array configuration, with baselines of 10m to maximum 600m. A mission would consist of a detection phase followed by a characterisation phase which features a second visit to a detected exoplanet. However, the optimal design is an active field of study, with the relevant metrics an active topic of debate. In addition to formation flying, a single-mount configuration is an option. Besides an X-array configuration, a relevant candidate is a linear array representing a double bracewell as described by Angel and Woolf (1997), which has a broader central null to exclude an extended image of a central star. This thesis aims to describe a general framework applicable to multiple configurations. In a study for the yield, a detection criterion of S/N>7 was applied, where the S/N was taken as the square root of the sum of the instrumental and astrophysical noise. For 2m apertures and 5% throughput, in a radius of 20 pc, $\approx 550$ exoplanets with radii between 0.5 and 6 $R_E$ could be detected, of which 25-45 rocky exoplanets within the empirical habitable zone. A large number of these planets would be around nearby M-dwarfs, to which IR interferometry is more sensitive.

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/Artists-impression-LIFE.png" width="700" height="300">

**5. Physics of nulling interferometry**
	- To be included!
	- Optics, diffraction, PSF, wavefronts etc.

Direct observation of a star and a planetary companion has several challenges. 

**1) Diffraction limit:** The separation between the star and the exoplanet is very small (in the order of AUs) compared to the distance between us and the star. Therefore the *angular separation* is very small, smaller than existing telescopes can resolve. Even if a source is infinitely small, its image on a telescope is not infinitely small, but given by the diffraction limit:
$$\theta \propto \lambda/D$$
A point source of wavelength $\lambda$ measured by a telescope with aperture size D results in an image with an characteristic angular size of $\lambda/D$. This effect is called diffraction and is relevant when imaging very small objects. The focal properties of the telescope can be described as an aperture function $A(x,y)$ in the spatial domain. The Point Spread Function in the image plane then is the squared magnitude of the Fourier transform of the aperture function. 

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/PSFdummy.png" width="700" height="300">

Single mirrors of telescopes cannot practically become much larger than 8 meters. This can be partly overcome by using segmented mirrors, but the all kinds of technical and financial issues grow with telescope size. 

**2) Contrast problem:** The light of the star is 6 to 9 orders of magnitude brighter than the planet, so it is difficult to have an instrument with sensitivity to the exoplanet while not being saturated by the light of the host star.

**3) Atmospheric window:** In case of ground-based observations, a part of the spectrum is absorbed by Earth's atmosphere. It is exactly the spectral lines of some of the molecules of interest to study exoplanet atmospheres, that are blocked for ground-based telescopes. However, as this thesis considers a space-based telescope, this will not be a problem.  

**4) Optical aberrations:**
The light does not arrive as perfect plane waves due to interaction with the interstellar medium between the star an the telescope. This problem is largely mitigated for a space-based telescope compared to a ground-based telescope, for which turbulence in the atmosphere is the main source of wavefront distortion. The result is optical aberrations, where each short exposure is slightly different so that adding it up to a long exposure doesn't improve the image but just adds up to a blurry image. For ground-based telescopes, it is common to include adaptive optics. This instrument measures the wavefront errors and adapts the shape of a mirror to correct for them.

The combination of working on a high contrast problem within the diffraction limit is as follows: When studying a hypothetical similar object to the host star with the same PSF, that PSF cannot simply be subtracted from the system of interest to end up with the planet signal, because the arrival of photons is a fundamentally random Poisson distributed process. The standard deviation of the photon count scales with the square of the photon count, which means that the error between to images of the same object taken at different moment is many orders of magnitude larger than the planet signal of interest.

Therefore it is key to remove starlight before recording it. One method is coronography, a method to block on-axis light before the light is recorded by the telescope. This method is not discussed in detail here. A fundamental limitation to coronography is still the diffraction limit of the aperture, which means the spatial resolution is limited. This is a limiting factor for the study of terrestrial exoplanets in habitable zones, because they have smaller angular separations from their host stars than gas giants and ice giants. Angular resolution can be enhanced one step further via interferometry. 


**Interferometry**

In interferometry, the baseline between the two telescopes determines the diffraction limit. An interferometer probes spatial frequency coordinates $u, v$ of an object.
As stated by the Van Cittert-Zernike theorem, the Fourier transform of the intensity distribution function of a distant, incoherent source is equal to its complex visibility $V(u,c)$:
$$V(u,v)=|V|e^{i\Phi}=\frac{\int \int I(\alpha, \beta)e^{-2 i \pi (\alpha u + \beta v)} d(\alpha, \beta)}{\int \int I(\alpha, \beta) d(\alpha, \beta)}$$
Different spatial shapes and symmetries translate to different spatial frequencies. With knowledge on the expected source, the best fit with some template function can be found, based on which the original image can be reconstructed.

Light of different places is combined. Under the condition that light is coherent, it creates a fringe pattern of constructive and destructive interference bands. The coherence condition means that the phase difference between the light measured at the different subapertures is constant over time. This means several things: The wavefront is intrinsically spatially coherent, so the source light is in phase with itself along the distance of the baseline. This also implies that the light is confined to a narrow wavelength bandwidth. In practice this can be worked around by channelling chromatic light into separate channels. This is generally true for a distant source that passes through space in the absence of turbulence. It also requires the optical path from both subapertures to the combiner to be exactly the same length, to a fraction of the wavelength of interest. 

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/OPDdummy.png" width="600" height="300">

From the fringe pattern created by multiple subapertures, one an extract an amplitude and a phase of the input signal. Besides these, there are indirect parameters specific to the combination or differences of signals from the different subapertures, that contain extra information on the source, such as the closure phase in classic interferometery and the kernel phase in nulling interferometry. It is common to rotate the interferometer, and sometimes one can also study different wavelengths of baselines. As mentioned before, different spatial shapes and symmetries translate to different spatial frequencies, and these also vary characteristically with rotation angle, wavelength or baseline. Based on the pattern in the $(u,v)$-plane, the original image can be reconstructed. For example, a point source leaves a different image than an extended central disk, a thin ring leaves a different image than a thick ring, and binary system with one on-axis and one off-axis point source leaves a different image than an on-axis point source. 



**Nulling interferometry**

A symmetric source will have light that perfectly interferes, which means that the minima of the fringe pattern cancel to zero flux. However, if the object of interest also contains off-axis light (for example, a planet), then the flux from the asymmetric source does not exactly equal zero at the bottom fringe. Nulling interferometry studies this 'almost null' in detail. However, in practice, the non-zero signal at the 'null' is also plagued by several sources of astrophysical and instrumental noise. This is predominantly due to leakage of random photon noise from the host star, which scales with the square-root of the host star's intensity. Therefore, further solutions are needed in the optical design of interferometers to separate the noise contributions from the off-axis signal of the planet, where possible before detection, and otherwise in post-processing, by applying statistical methods that are grounded in a physical understanding of the noise. 

**Single Bracewell**

The concept of nulling interferometry for exoplanet detection was first proposed by Bracewell in 1978, well before the existence of exoplanets had been proven. He proposed centering an interference null on the star in infrared wavelengths, and to subsequently disentangle the remaining star versus planet signal by rotating the interferometer, which would modulate the (asymmetric) planet signal on top of a constant stellar flux. 

**Double Bracewell**

Since then, different architectures have been proposed, including a double Bracewell by Angel & Woolf in 1996. The architecture of the single Bracewell was such that an improved angular resolution due to longer baselines also resulted in a larger image of the stellar disk, beyond the angular distance that was effectively nulled, thereby compensating the better resolution with stellar leakage. The double Bracewell consists of two pairs of nulls, that are offset by half a phase. When combining the two signals, this results in one differential null whose amplitude is zero, also in both the first and second order sensitivity to errors. This leads to an interference pattern with a central null that is broader in its angular extent, as shown in the Figure below. Upon rotation of the interferometer, contributions from an elliptical dust cloud would modulate at twice the frequency of a planet due to its symmetry axis, and therefore could be distinguished from a planet signal. Angel and Woolf showed that their architecture allowed for baselines extended beyond 50 meters, and found through simulations that this would allow for the detection of exoplanets at distances of 10 parsecs provided sufficient exposure time and exozodiacal dust levels similar to our own solar system.
<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/AngelWoolf.png" width="800" height="250">

**Kernel nulling using a four-beam combiner**
In addition to stellar photon noise due to leakage of central light into the null, on of the main practical limitations of nulling interferometry are piston errors. These errors originate from imperfect relative alignment of the apertures, leading to optical path difference errors that result in phase differences, causing the interference to be imperfectly destructive. Martinache & Ireland (2018) proposed a kernel-nulling architecture that is independent of path length errors. The term 'kernel' implies a quantity that is robust to certain errors. Although their architecture of four inputs, one bright output and three dark outputs, was proposed for the ground-based Very Large Telescope, the same concept works for space. They present a method for self-calibration from linear combinations of polluted data that reside in a space orthogonal to the source of perturbation. As shown in the figure below, the four inputs T 1 to 4 are transferred into one bright and three dark channels via coupling matrix N. Assuming small phase errors, a first order Taylor expansion allows one to write the electric field as a linear function of phase: $E_k=e^{(- j \phi_k)} \approx 1 - \phi_k$. The phases of the three nulled intensities can then be defined with respect to the phase of aperture $T_0$ as $x=\frac{1}{4} \times \left[ (+\phi_1-\phi_2-\phi_3)^2,  (-\phi_1+\phi_2-\phi_3)^2,  (-\phi_1-\phi_2+\phi_3)^2 \right]$. The piston-induced leak leak of the nuller is a function of the six parameters $\phi_k^2, \phi_k \times \phi_{l \neq k}$. To overcome this underconstrained problem, each nuller output is split into two non-symmetric outputs with help of a pre-defined phase offset $\theta$ as indicated by the 'phase shifters' in the figure below, so that the detector records an intensity vector of six components.

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/MartinacheIreland1.png" width="800" height="250">

The sensitivity of the instrument to second-order phase differences can then be collected in a 6x6 matrix $A$. The next step is to find a kernel operator $K$ for which $K \dot A = 0$, so that the kernel applied to the six channels of the intensity vectors is independent of second-order phase differences: $K \cdot x$. As an example shown in the figure below, a phase shift of $\theta/2$ applied to the three nulled outputs results in six outputs with sensitivities as shown in the middle panel. The combination of signals 1-2, 3-4 and 5-6, leads the three kernel-nulls in the right panel, which are independent of phase shift errors to second order. 
<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/MartinacheIreland23.png" width="800" height="320">
Given the spatial set-up of the detectors, the sensitivity of the 6 channels can be transferred into a transmission map, mapping the sensitivity of each channel to a 2D position in the sky, which now contain several anti-symmetries that constrain more information on location than the classically nulled channels. 

Another key result of kernel-nulling is that the error in the null depth of the kernel-null outputs is symmetric, with uncertainties proportional to the cube of the phase errors: The figure below shows the results of a simulation of 50nm RMS residual piston excursions. In the presence of piston errors, the theoretical values become skewed distributions for both the 3 classically nulled outputs and the 6 modified outputs. This poses problems in post-processing [elaborate], while the right panel shows the symmetric kernel outputs whose mean value coincides with the true null depth and whose variance reveals the phase error.

**Kernel nulling using an arbitrary number of apertures**
Laugier et al. (2020) present a method to find combiner matrices for an arbitrary number of apertures, allowing to compare many different architectures.

They first argue that, even though combiner matrices of all-in-one kernel nullers grow in complexity for larger numbers of outputs, they still seem to be the most efficient in comparison to multiplexing nullers, which work independently over a smaller number of sub-apertures, conserving a higher throughput per output at the cost of fewer robust observables and lower overall efficiency. The number of independent kernel nulls is always the same of the number of independent closure phases: $\frac{(n_a-1)(n_a-2)}{2}}$. 

Using a Mont Carlo simulation with Gaussian input piston errors, the standard deviation was evaluated for a number of architectures: The simple Bracewell nuller, a 3-telescope kernel nuller (one bright channel, two dark channels, one independent kernel null), a 4-telescope kernel nuller (one bright channel, three dark channels, six independent nulls, and three independent kernel nulls) and a 6-telescope kernel nuller (one bright channel, 20 independent nulls, and 10 independent kernel nulls). The figure below shows the suppression of second-order phase noise by the kernel nulls (solid lines) compared to the classical nulls (dashed lines) provided that the input piston errors remain within the small phase approximation (left part of the plot). 

<img src="https://github.com/LoesRuttenGithub/Thesis-Noise-Nulling/blob/main/Figures/Laugier20208.png" width="500" height="320">


**Nulling interferometry in practice**

Nulling interferometry has been shown to work in practice, namely through the Multiple Mirror Telescope in Arizona (Hinz et al. 1998) and the Keck Interferometer Nuller using a double Bracewell setup, the Large Binocular Telescope Nulling Optimised Mid-Infrared Camera, GLINT integrated optics nulling combiner at Subaru/SCExAO, VLT/NACO from which a spectrum was obtained, ExoGRAVITY Lacour et al. 2019 extremely good astrometric position of the planet. 

Example figures from Lagadec et al. 2018

**Sources of noise.**
Following the notation of Huber et al. 2024, the photon rate at a certain bandwidth at a certain output channel $j$ of an interferometer is the sum of the photons collected within a wavelength span due to the instrument response $R$ to the image on the sky $I_{sky}$. In particular, the instrument response is related to the square modulus of the output electric field of output channel $j$, which is a function of the amplitudes, phase shifts, polarisation angles and positions of the $k$ collectors.
$$ N_j(t, \lambda, \alpha, \beta) = w(\lambda)\int_\Omega R_j(t,\lambda,\alpha,\beta) I_{sky}(t,\lambda,\alpha, \beta) d\Omega \\
\propto f_j (A,\phi_j,\theta, x(t), y(t), \lambda, \alpha, \beta $$

The effect of instrumental perturbations on the photon rate can then be written as $$ (N_j+\delta N_j)(t, \lambda, \alpha, \beta) \propto f_j (A \circ \delta A(t,\lambda),\phi_j + \delta \phi_j(t,\lambda),\theta + \delta \theta (t,\lambda) , (x+\delta x)(t), (y+\delta y)(t), t, \lambda, \alpha, \beta $$



**6. Sources of fundamental noise & LIFEsim**
Fundamental noise, also referred to as astrophysical noise, is an umbrella term for all true light signals that are contained in the measurement, that are not from the exoplanet of interest. In particular, these fundamental sources include, as described in Dannert et al. (2022):
- True signal of the star, modelled as black body radiation, neglecting limb darkening i.e. assuming constant brightness across the disk
- True signal of the planet, modelled as black body radiation
- Geometric stellar leakage: not all stellar flux is suppressed when the angular extent of the disk of the star is wider than the 'null' of the interferometer
- Local zodiacal dust: Dust inside our own solar system adds radiation background via scattering of visible sunlight (not relevant for a MIR instrument) and thermal emission. Earlier work on DarwinSIM included a model of spectral surface brightness. Local zodiacal light is *diffuse*, which means it comes from different directions, with different phases, rather than a coherent wavefront, so it cannot be brought to destructive interference. The most effective way of minimizing this contribution is by observing away from local dust.
- Exozodiacal dust: Generally modelled as optically thin, symmetric with a power law surface density. The temperature of the dust depends on distance from the central star. Most of the radiation originates from the (hottest) central region. High-luminosity stars are expected to have larger disks. Larger disks have larger surface areas, thus larger fluxes. A symmetric disk can be filtered out, but its light does contribute to shot noise.
A software called LIFESim has been developed in the context of the LIFE mission to simulate its yield. The code currently includes fundamental noise with the option to include instrumental noise in the future. 

P-pop MC tool: To simulate the yield of realistic parameters for the LIFE mission, it needs to be tested on a synthetic population of exoplanets. This synthetic population is built using P-Pop [add reference and some sentences].  

**7. Sources of instrumental noise**
- Sources of instrumental noise [description by Lay]


In a comparison of the fundamental noise sources and instrumental noise, Dannert et al. (2022) report that the short wavelength section is dominated by astronomical noise from the host star, and that the low wavelength section is dominated by local zodiacal dust. Systematic noise is the largest at shorter wavelengths, which becomes relevant when considering that the fundamental noise limit must hold for all wavelength bins. Regarding photon noise sources, both contributions are modelled as constant, where the requirements on the dark current seem to be one order of magnitude more stringent than what has been achieved so far.

hallohiero

- Explanation how architecture leads to correlated noise
- Study by Lay
- Study for Asgard/NOTT


**8. State of the art**
**8c. Open questions for numerical simulations**

**8b. Open questions for instrumental noise**
- Consequence of sub-optimal alignment, error trifecta
- Semi-analytical framework for underlying physics of combiner matrix
- Standardised language and approach for instrumental errors
- To be elaborated...


**8c. Further development of detection statistics**
- Statistical framework by Ceau for Kernel-phase detection limits
- Covariance matrix
- $\tau(z)$: While subtracting electric field vectors from two sources at symmetric positions in the field of view cancel out kernel nulls $k(z)$, the complementary observable to the kernel null, Laugier et al. (2020) suggest that $\tau(z)$, the sum of the two outputs, may encode information about the input phase errors. It will be a subgoal of this thesis to explore the information contained in $\tau(z)$.


## Parking lot for ideas
- Compare with literature for error propagation of LISA (Laser Interferometer Space Antenna for Gravitational Waves)
- To do: Study Principles of Stellar Interferometry pdf


## Syntax for images, links, code snippets etcetera

An **image**:

![alt text](http://picsum.photos/200/200)

What is your favourite `variable`?

A *fantastic* script: 
```javascript
let num = Math.random();
``` 

> This is a block quote

[Useful link?](https://youtube.com/watch?v=dQw4w9WgXcQ)

~~got you!~~
