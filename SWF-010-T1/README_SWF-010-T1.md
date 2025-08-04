# SWF-010-T1: Source Finding Test Case – SRCNET v0.1

**Test name:** SWF-010-T1: Radio continuum source finding with PyBDSF

**Documentation on confluence:** https://confluence.skatelescope.org/x/5pUSEw

**Summary:** This notebook serves as a test case for evaluating the source-finding using **PyBDSF** on a real-world radio continuum dataset from the **LoTSS** survey. The goal is to extract and visualise radio sources from a single-frequency 2D mosaic using PyBDSF’s Gaussian decomposition approach.

## Dependencies

This notebook is written in python and has the following dependencies:

- astropy
- bdsf
- matplotlib
- numpy
- pandas
- pillow
- scipy	

## Running the code

The code can be run as either a jupyter notebook (`SWF-010-T1.ipynb`) or a python script (`SWF-010-T1_workflow.py`) via

```
python SWF-010-T1_workflow.py my_folder
```
where `my_folder` is the folder that includes the `teal/` folder that hosts the data necessary to run these tests (the data must have been staged beforehand, if running on a specific node).

If running via the python script, the output figures will be saved to the current working directory as PNGs in a `testcases-results/` folder, that is created if non-existant.

## Description

### Data
We use a small-area FITS mosaic from **LoTSS-DR2**, publicly available via the [LOFAR Surveys website](https://lofar-surveys.org/releases.html). The image has a resolution of ~6 arcsec and includes typical survey noise and structure.

### Content
In this test we perform:

- FITS header inspection
- Radio source detection:
  - Background estimation
  - RMS noise map creation
  - Island detection and Gaussian fitting
- Extraction of source and Gaussian component parameters (position, flux, shape, ect)
- Overlay of source positions on the original radio image
- Basic output validation through visual inspection by confirming source positions match radio emission regions that must the same as showed on the plot in the confluence page

### Expected outputs
- Catalogues that contain information about the sources and the Gaussians (not exported to file currently)
- Each source and Gaussian contain basic information such as
  - RA, Dec
  - Peak and total flux
  - Source size and shape (major, minor axis, position angle)
  
### Outputs for validation
- Plot showing detected sources overlapping the radio maps

## SRCNet Links

### Confluence

- [Test description](https://confluence.skatelescope.org/display/SRCSC/SWF-010-T1+Science+test+case)

### Jira

Features
- [SP-5206](https://jira.skatelescope.org/browse/SP-5206)
- [SP-5593](https://jira.skatelescope.org/browse/SP-5593)
Tickets
- [TEAL-924](https://jira.skatelescope.org/browse/TEAL-924)
- [TEAL-931](https://jira.skatelescope.org/browse/TEAL-931)
- [TEAL-1045](https://jira.skatelescope.org/browse/TEAL-1045)

## Contributors

- Lara Alegre
- Adélie Gorce
- TEAL team