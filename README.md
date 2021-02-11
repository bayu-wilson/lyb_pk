# A measurement of the Lyβ forest power spectrum and its cross with the Lyα forest in X-Shooter 100

This pipeline was used to create the results of Wilson et al. 2021. Read the following abstract for more physical context.

### Abstract
The Lyα forest is the large-scale structure probe for which we appear to have modeling control to the highest wavenumbers, which makes it of great interest for constraining the warmness/fuzziness of the dark matter and the timing of reionization processes. However, when using the standard statistic, the Lyα forest power spectrum, there are still large parameter degeneracies that limit inferences at these highest wavenumbers, such as between the gas temperature and particle mass in warm/fuzzy dark matter models. With the aim of breaking these degeneracies, we measure the power spectrum of the Lyβ forest and its cross correlation with the coeveal Lyα forest using the one hundred spectra of z= 3.5−4.5 quasars in the VLT/X-Shooter XQ-100 Legacy Survey, motivated by this transition’s lower absorption cross-section that makes it sensitive to somewhat higher densities relative to the Lyα transition. Our measurements of the temperature-density relation are consistent with the recent Lyα forest measurements of Gaikwad et al. (2020), providing a consistency check on IGM models that explain the Lyα forest using instead the Lyβ forest. The z=3.5−4.5 trends we find in the Lyβ forest show a similar consistent flattening of the slope of the temperature-density relation with decreasing redshift from Heii reionization. The limiting factor in our analysis is the significant uncertainty in the effective spectral resolution of X-Shooter spectrograph.  Our competitive constraints marginalize over this uncertainty. This plus the significant improvement in our constraint over the Lyα power spectrum alone from our data set suggest that a similar measurement of the Lyβ forest in another spectroscopic data set could result in a significant improvement.

### Data
We use 100 z~3.5-4.5 quasar spectra from the XQ-100 legacy survey (Lopez et al. 2016). You may download the data [here](https://www.dropbox.com/sh/eijuc5jhg4olo0x/AAAGRMf110uiSAe49L3_RJ7Ga?dl=0).
The analysis pipeline was tested on synthetic mock quasar Lyα and Lyβ forests from the [Sherwood Simulation suite](https://www.nottingham.ac.uk/astronomy/sherwood/). [Here](https://www.dropbox.com/sh/c4zr9pbd5zg8i8d/AADp6uloIl6nBWEOmGcwir2Oa?dl=0) is an example of a mock XQ-100 dataset.

### How to run the code
1) Go to desired local directory and clone the repository: `git clone git@github.com:bayu-wilson/lyb_pk.git`
2) Download observational data [here](https://www.dropbox.com/sh/eijuc5jhg4olo0x/AAAGRMf110uiSAe49L3_RJ7Ga?dl=0) and synthetic data [here](https://www.dropbox.com/sh/c4zr9pbd5zg8i8d/AADp6uloIl6nBWEOmGcwir2Oa?dl=0). 
3)
