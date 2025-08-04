# SWF-008-T2

This test loads simulated visibilities stored as uvfits files containing a mock EoR signal and foregrounds and estimates the EoR power spectrum from these visibilities using hydra-pspec. The test notebook walks through all of the steps required to run a hydra-pspec analysis and produces plots of the cylindrically and spherically averaged power spectrum of the EoR.

For more details, please see this [confluence page](https://confluence.skatelescope.org/x/o2FoEw).

## Dependencies

This notebook is written in python and has the following dependencies:

- numpy
- python
- pyuvdata
- matplotlib
- [hydra-pspec](https://github.com/HydraRadio/hydra-pspec.git)

A `mamba/conda` yaml has been provided for convenience.  All of these dependencies can be installed via

```
mamba env create -f environment.yaml
```

To install with `conda`, replace `mamba` with `conda` in the above command.

## Running the code

The code can be run as either a jupyter notebook (`SWF-008-T2.ipynb`) or a python script (`SWF-008-T2.py`) via

```
python SWF-008-T2.py
```

If running via the python script, the output figures will be saved to the current working directory as PDFs following a `figure-#.pdf` syntax.  For example, after running the script, Figure 1 in the notebook can be found in `figure-1.pdf`.

## SRCNet Links

### Confluence

- [Test description](https://confluence.skatelescope.org/x/o2FoEw)

### Jira

- [TEAL-954](https://jira.skatelescope.org/browse/TEAL-954)
- [TEAL-967](https://jira.skatelescope.org/browse/TEAL-967)
- [TEAL-1045](https://jira.skatelescope.org/browse/TEAL-1045)

## Contributors

- [Burba, Jacob](https://github.com/jburba)
