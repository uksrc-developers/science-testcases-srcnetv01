# SWF-008-T1: Image-based EoR power spectrum Test Case - SRCNet v0.1

**Test name:** SWF-008-T1: Image-based EoR power spectrum Test Case

**Documentation on confluence:** https://confluence.skatelescope.org/x/uKASEw

**Summary:** This test assesses the ability of the SRCNet v0.1 to produce a cylindrical power spectrum given an image field, in FITS format, based on data from the SKA Data Challenge 3b. The main application is in the context of the Epoch of Reionisation science case. There is a second leg to this test, SWF-008-T2, which computes the cylindrical power spectrum from a set of visibilities.

## Dependencies

This notebook is written in python and has the following dependencies:

- [ps_eor](https://gitlab.com/flomertens/ps_eor/)
- scipy
- astropy
- h5py
- healpy
- time
- matplotlib
- tables

## Running the code

The code can be run as either a jupyter notebook (`SWF-008-T1.ipynb`) or a python script (`SWF-008-T1.py`) via

```
python SWF-008-T1.py my_folder
```
where `my_folder` is the folder that includes the `teal/` folder that hosts the data necessary to run these tests (the data must have been staged beforehand, if running on a specific node).

If running via the python script, the output figures will be saved to the current working directory as PNGs in a `testcases-results/` folder, that is created if non-existant.

## Description

### Data
We use data produced in the context of the **SKA data challenge 3b**: a simulated image (IM1), the corresponding simulated instrument PSF (All_PSF), and the cylindrical power spectrum of the image, provided as reference (PS3). All data are publicly available. 

### Content
In this notebook, we
- Load the image data and PSF, provided as FITS files
- Save reduced datasets to FITS files, in order to limit computing needs
- Compute the power spectrum of the reduced dataset with the `ps_eor` package (https://gitlab.com/flomertens/ps_eor/)
- Compare the output image to PS3.

### Expected outputs
- Figure of the cylindrical power spectrum.
  
### Outputs for validation
- Output is validated by comparing images between output of test on SRCNet and test ran on known infrastructure. See documentation for expected outputs.

## SRCNet Links

### Confluence

- [Test description](https://confluence.skatelescope.org/x/uKASEw)

### Jira

Features
- [SP-5206](https://jira.skatelescope.org/browse/SP-5206)
- [SP-5593](https://jira.skatelescope.org/browse/SP-5593)
Tickets
- [TEAL-923](https://jira.skatelescope.org/browse/TEAL-923)
- [TEAL-930](https://jira.skatelescope.org/browse/TEAL-930)
- [TEAL-1045](https://jira.skatelescope.org/browse/TEAL-1045)

## Contributors

- Adélie Gorce
- Florent Mertens
- TEAL team