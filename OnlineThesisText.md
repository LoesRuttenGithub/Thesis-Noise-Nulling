# Online thesis text


Here I will update my work on my thesis as we go. It is intended to be structured but not as formal as the final text. To be discussed: Which content can be in a public folder.

## Table of contents
[Literature list](#Literature-list)

[Meetings](#Meetings)

[Planning](#Planning)



## <a name="Literature-list"></a> Literature list 
<details>
<summary> Quantz et al. 2022: LIFE Paper 1: Improved exoplanet detection yield estimates </summary>

Description of LIFE space mission with 4 formation flying infrared telescopes. Paper presents an intstrument simulator with all sources of astrophysical noise (but not yet including instrumental noise) which was coupled with the P-Pop Monte Carlo tool of a synthetic exoplanet population. The gain (yield per unit time) was averaged over 500 MC realisations.

Based on 4 2m apertures, it was estimated that within 20 p of the sun, 25-45 rocky exoplanets within the habitable zones could be detected, and doubling these numbers with 3,5m apertures. In a first observation, estimates of radii and effective temperatures of exoplanets could be made, followed by estimates of thermal emission spectra in a second visit.

It is known from Kepler, Tess and Radial Velocity surveys that planets similar to Earth should be very abundant. In the future, James Webb Space Telescope may reveal if planets around red dwarfs can retain their atmospheres despite high levels of activity of their host stars. 

ESA's Ariel will provide transmission and emission spectra of exoplanets with warm hydrogen-dominated atmospheres. The ~30m ground based Extremely Large Telescopes will be able to detect thermal spectra of small planets around nearby stars using the mid-infrared METIS spectrograph. In addition, the PCS and HIRES will detect via reflected light.

[Add some words about PLATO]

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
- Abberations in the beam produce variations in the focal plane complex amplitude upon injection into single-mode waveguides
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

<summary> Martinache & Ireland 2018: Kernel-nulling for a robust direct interferometric detection of extrasolar planets </summary>
</details>

<details>
<summary> Ceau et al. 2019: Kernel-phase detection limits (JWST) </summary>
</details>

<details>
<summary> Lay 2004: Systematic errors in nulling interferometers </summary>
</details>

<details>
<summary> Bracewell 1978: Detecting nonsolar planets by spinning infrared interferometer </summary>
</details>

<details>
<summary> Dannert et al. 2022: LIFE paper 2: Signal simulation, extraction and exoplanet parameters from single epoch observations </summary>
</details>

<details>
<summary> Angel and Woolf 1996: An imaging nulling interferometer to study extrasolar planets </summary>
</details>

<details>
<summary> Laugier et al. 2020: Kernel nullers for an arbitrary number of apertures </summary>
</details>

<details>
<summary> Flasseur et al. 2024: Shrinkage MMSE estimators of covariances beyond the zero-mean and stationary variance assumptions </summary>
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




- WG 3.2 meeting 11 July
	- Jens Kammerer (topic: yield estimation), Daniel Angerman (project manager), Felix Dannert (PhD on instrumental noise), Angelica Kovacevic (extragalactic time astronomy), Philipp Huber (PhD on instrumentation) 
	- Presentation on PHIRINGE, photon count data generator for nulling interferometers. Takes config files based on astropy and scipy. Uncertainty in flux estimated from spread in (MC simulated?) data. Open question is to find the cleanest way to estimate uncertainty. 
	- SCIFYsim, initially for NOTT. Efforts made to include symbolic math, but so far not as useful as hoped. Instrumental noise is correlated, interesting to extract. Tip: old gravitational wave literature has analysis on the noise introduced by different tools. 
	- Discussion afterwards: Loes has been added to LIFE Slack. Most relevant groups are 3.2 LifeSIM, 3.3 signal procext, 4.1 nulling tech. Thesis could include some joint interphase with PHIRINGE, SCIFYsim and LIFEsim. Datastandard from NIFits group should provide 'one way of doing stuff'. Other idea: explain impact of Lay's assumption of optical alignment and correlated terms within sum for non-ideal case. With Monte Carlo you can only 'retrofit' the underlying distribution if you have strong confidence in the underlying model. Some new numerical tools could be explored, machine learning algorithms that do a symbolic search. 


- **Romain & Ida & Lian 11 July*
	- MC simulations take long. Idea to use Python functionality to dump outcomes into file and load from there.
	- Python code gives error on Windows because it doesn't have a /tmp folder as in Unix based systems.
	- Idea to have joint meeting with Delft & Leuven every other week? To be determined.
	- 3rd week of August, internal office moving in Delft should be done

- **Meeting Loes & Romain 12/7**
	- Question why linear setup is still interesting compared to 2D layout. Has to do with broader rejection, so it still works if your star is resolved as a disk. Key concept is spatial coherence, explained further in Angel & Woolf. Useful to review theory ondiffraction and point spread function of telescopes.
	- Does nulling interferometry also work for binaries, as some research concluded that the majority of stars are binary systems? Yes, as long as the other star is not in the field of view. Moreover, mostly heavy stars are binaries, and in the solar neighbourhood are many lighter stars.
	- Question about symmetric baselines: They are redundant if they are exactly the same, skew gives more results. Best layout to be studied further. Can we choose optimal baseline for specific molecular lines in the spectrum? Not needed, because planet generally (if not face-on) will move over different angular separations as projecten on the line of sight. 
	- Hint for JWST Kernel-phase paper: kernel-phase is not the same as kernel-null. It's more like a closure phase. Focus on statistics in that paper. Tutorial based on on Ceau does not include correlation and whitening yet. Main question: what to put in matrix that is representative?


## Planning

Just a first rough outline as conversation starter! Anything can be adapted.

| Month | Activity | Hours | Comments |
| --- | --- | --- | --- |
| July & August | Exploration: Literature and research plan | 150h| Fulltime except last 2 weeks of July |
| September | Deep thinking: Extend theoretical concept | 150h| Fulltime |
| October-December | Implementation: Create new tool or integrate existing tools | 150h| ~2 days per week. Romain 5 weeks in Nice from Nov |
| January- mid-February| TU Delft-based applications | 150h | ~2 days per week|
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

## First quick literature review
- Occurrence of exoplanets as we know so far and main methods 
- Parallel initiatives to determine exoplanets (LUVOIR, etc)
- Assets of infrared nulling interferometry
- Idea behind nulling interferometry
- Sources of fundamental noise
- Sources of instrumental noise

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
