# Online thesis text


Here I will update my work on my thesis as we go. It is intended to be structured but not as formal as the final text. To be discussed: Which content can be in a public folder.

## Table of contents
[Literature list](#Literature-list)

[Meetings](#Meetings)

[Planning](#Planning)



## <a name="Literature-list"></a> Literature list 
<details>
<summary> Quantz et al. 2022: LIFE Paper 1 </summary>
Description of LIFE space mission with 4 formation flying infrared telescopes. Paper presents an intstrument simulator with all sources of astrophysical noise (but not yet including instrumental noise) which was coupled with a Monte Carlo simulation of a synthetic exoplanet population. Based on 4 2m apertures, it was estimated that within 20 p of the sun, 25-45 rocky exoplanets within the habitable zones could be detected, and doubling these numbers with 3,5m apertures. In a first observation, estimates of radii and effective temperatures of exoplanets could be made, followed by estimates of of thermal emission spectra in a second visit.

It is known from Kepler, Tess and Radial Velocity surveys that planets similar to Earth should be very abundant. In the future, James Webb Space Telescope may reveal if planets around red dwarfs can retain their atmospheres despite high levels of activity of their host stars. ESA's Ariel will provide transmission and emission spectra of exoplanets with warm hydrogen-dominated atmospheres. The ~30m ground based Extremely Large Telescopes will be able to detect thermal spectra of small planets around nearby stars using the mid-infrared METIS spectrograph. In addition, the PCS and HIRES will detect via reflected light.

The US-led Habitable Exoplanet Observatory HabEx and the Large UV/Optical/IR Surveyor LUVOIR  is intended to study the atmospheres of several dozen small exoplanets in Habitable Zones via reflected light.

LIFE explores the development of a space-based mid-infrared nulling interferometer within ESA's Voyage 2050 program, building upon studies of the Darwin mission and Terrestrial Planet Finder - Interferometer
</details>


- Lay 
- LIFE paper 1
- LIFE paper 2
- Etc.

Anticipated papers:
- Photobombing

## Meetings 
- Romain & Loes 3 May
	- Topic 1
	- Topic 2
		1. Subtopic 1
		2. Subtopic 2
- Romain & Loes 29 May
- Leuven/Delft/Liege consortium meeting 10 July
- WG 3.2 meeting 11 July
- Romain & Loes 11 July
	- Added to LIFE Slack. Most relevant groups are 3.2, 3.3, 4.1
- Romain & Ida & Lian 11 July
	- MC simulations take long. Idea to use Python functionality to dump outcomes into file and load from there.
	- Python code gives error on Windows because it doesn't have a /tmp folder as in Unix based systems.
	- Idea to have joint meeting with Delft & Leuven every other week? To be determined.
	- 3rd week of August, internal office moving in Delft should be done


## Planning

Just a first rough outline as conversation starter! Anything can be adapted.

| Month | Activity | Hours|
| --- | --- | --- |
| July & August | Exploration: Literature and research plan | 150h|
| September | Deep thinking: Extend theoretical concept | 150h|
| October-December | Implementation: Create new tool or integrate existing tools | 150h|
| January- mid-February| TU Delft-based applications | 150h|
| mid-February - mid-March| Writing and knowledge transfer | 150h|
| mid-March - end-March| Buffer | 150h|



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
