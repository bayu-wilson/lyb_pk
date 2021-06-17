# A measurement of the Lyβ forest power spectrum and its cross with the Lyα forest in X-Shooter 100

This pipeline was used to create the results of [Wilson et al. 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv210604837W/abstract). Read corresponding abstract.

### Abstract
<!-- The Lyα forest is the large-scale structure probe for which we appear to have modeling control to the highest wavenumbers, which makes it of great interest for constraining the warmness/fuzziness of the dark matter and the timing of reionization processes. However, when using the standard statistic, the Lyα forest power spectrum, there are still large parameter degeneracies that limit inferences at these highest wavenumbers, such as between the gas temperature and particle mass in warm/fuzzy dark matter models. With the aim of breaking these degeneracies, we measure the power spectrum of the Lyβ forest and its cross correlation with the coeveal Lyα forest using the one hundred spectra of z= 3.5−4.5 quasars in the VLT/X-Shooter XQ-100 Legacy Survey, motivated by this transition’s lower absorption cross-section that makes it sensitive to somewhat higher densities relative to the Lyα transition. Our measurements of the temperature-density relation are consistent with the recent Lyα forest measurements of Gaikwad et al. (2020), providing a consistency check on IGM models that explain the Lyα forest using instead the Lyβ forest. The z=3.5−4.5 trends we find in the Lyβ forest show a similar consistent flattening of the slope of the temperature-density relation with decreasing redshift from HeII reionization. The limiting factor in our analysis is the significant uncertainty in the effective spectral resolution of X-Shooter spectrograph. Our competitive constraints marginalize over this uncertainty. This plus the significant improvement in our constraint over the Lyα power spectrum alone from our data set suggest that a similar measurement of the Lyβ forest in another spectroscopic data set could result in a significant improvement. -->
<!-- The Lyα forest is the large-scale structure probe for which we appear to have modeling control to the highest wavenumbers, which makes it of great interest for constraining the warmness/fuzziness of the dark matter and the timing of reionization processes. However, the standard statistic, the Lyα forest power spectrum, is unable to strongly constrain the IGM temperature-density relation, and this inability further limits how well other high wavenumber-sensitive parameters can be constrained. With the aim of breaking these degeneracies, we measure the power spectrum of the Lyβ forest and its cross correlation with the coeveal Lyα forest using the one hundred spectra of z=3.5-4.5 quasars in the VLT/X-Shooter XQ-100 Legacy Survey, motivated by the Lyβ transition's smaller absorption cross section that makes it sensitive to somewhat higher densities relative to the Lyα transition. Our inferences from this measurement for the IGM temperature-density relation appear to latch consistently onto the recent tight lower-redshift Lyα forest constraints of [Gaikwad et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020arXiv200900016G). The z=3.4-4.7 trends we find using the Lyα-Lyβ cross correlation show a flattening of the slope of the temperature-density relation with decreasing redshift. This is the trend anticipated from ongoing HeII reionization and there being sufficient time to reach the asymptotic temperature-density slope after hydrogen reionization completes. Furthermore, our measurements provide a consistency check on IGM models that explain the Lyα forest, with the cross correlation being immune to systematics that are uncorrelated between the two forests, such as metal line contamination (and do suggest a mild tension). -->
The Ly-α forest is the large-scale structure probe for which we appear to have modeling control to the highest wavenumbers, which makes it of great interest for constraining the warmness/fuzziness of the dark matter and the timing of reionization processes. However, the standard statistic, the Ly-α forest power spectrum, is unable to strongly constrain the IGM temperature-density relation, and this inability further limits how well other high wavenumbersensitive parameters can be constrained. With the aim of breaking these degeneracies, we measure the power spectrum of the Ly-β forest and its cross correlation with the coeveal Ly-α forest using the one hundred spectra of z=3.5−4.5 quasars in the VLT/X-Shooter XQ-100 Legacy Survey, motivated by the Ly-β transition’s smaller absorption cross-section that makes it sensitive to somewhat higher densities relative to the Ly-α transition. Our inferences from this measurement for the IGM temperature-density relation appear to latch consistently onto the recent tight lower redshift Ly-α forest constraints of [Gaikwad et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020arXiv200900016G). The z=3.4−4.7 trends we find using the Ly-α–Ly-β cross correlation show a flattening of the slope of the temperature-density relation with decreasing redshift. This is the trend anticipated from ongoing HeII reionization and there being sufficient time to reach the asymptotic temperature density slope after hydrogen reionization completes. Furthermore, our measurements provide a consistency check on IGM models that explain the Ly-α forest, with the cross correlation being immune to systematics that are uncorrelated between the two forests, such as metal line contamination.



### Data
We use 100 z~3.5-4.5 quasar spectra from the XQ-100 legacy survey [(Lopez et al. 2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...594A..91L/abstract). You may download the data [here](https://www.dropbox.com/sh/eijuc5jhg4olo0x/AAAGRMf110uiSAe49L3_RJ7Ga?dl=0).
The analysis pipeline was tested on synthetic mock quasar Lyα and Lyβ forests from the [Sherwood Simulation suite](https://www.nottingham.ac.uk/astronomy/sherwood/). [Here](https://www.dropbox.com/sh/c4zr9pbd5zg8i8d/AADp6uloIl6nBWEOmGcwir2Oa?dl=0) is an example of a mock XQ-100 dataset.

### How to run the code
1) Go to desired local directory and clone the repository: `git clone git@github.com:bayu-wilson/lyb_pk.git`
2) Download observational data [here](https://www.dropbox.com/sh/eijuc5jhg4olo0x/AAAGRMf110uiSAe49L3_RJ7Ga?dl=0) and synthetic data [here](https://www.dropbox.com/sh/c4zr9pbd5zg8i8d/AADp6uloIl6nBWEOmGcwir2Oa?dl=0). 
3) In control center for the pipeline, `lyb_pk/pipeline/inis.py`, you can decide to use observational data (set `use_obs=1`) or use synthetic data (set `use_obs=0`). This file contains other flags that may be changed. For example, you may choose whether to remove DLAs, use logarithmic k-binning, or make a continuum correction based off of Faucher-Giguére+2008.
4) In the `lyb_pk/pipeline/` directory, run `python main.py`. This calculates the mean flux in each redshift bin and power spectrum in each redshift \& wavenumber bin. This outputs the results into `lyb_pk/output/`.
5) The uncertainties of the measurements are quantified using the bootstrap method. Run `python boot_indo.py` to calculate the error bars.
6) Finally, the `lyb_pk/plot/` directory contains plotting routines to make the mean flux (`plot_meanflux.py`) and power spectra (`plot_pk.py`) figures similar to those in the paper.
