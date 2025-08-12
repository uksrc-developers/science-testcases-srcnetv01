# SWF-005-T1

This notebook uses `astroquery` to search for an image from the [Faint Images of the Radio Sky at Twenty cm (FIRST)](https://sundog.stsci.edu/) VLA radio survey and plots it using `matplotlib`.  The source was chosen to match the [documentation](https://astroquery.readthedocs.io/en/latest/image_cutouts/first/first.html) for `astroquery.image_cutouts.first`.

For more details, please see the [SWF-005-T1 confluence page](https://confluence.skatelescope.org/x/AqkSEw).

## Dependencies

The notebook and script are written in python and has the following dependencies:

- astropy
- astroquery
- python
- matplotlib

A `mamba/conda` yaml has been provided for convenience.  All of these dependencies can be installed via

```
mamba env create -f environment.yaml
```

To install with `conda`, replace `mamba` with `conda` in the above command.

## Running the code

The code can be run as either a jupyter notebook (`SWF-005-T1.ipynb`) or a python script (`SWF-005-T1.py`) via

```
python SWF-005-T1.py
```

If running via the python script, the output figure is generated as a PDF in the current working directory named `image.pdf`.

## SRCNet Links

### Confluence

- [Test description](https://confluence.skatelescope.org/x/AqkSEw)

### Jira

- [TEAL-921](https://jira.skatelescope.org/browse/TEAL-921)
- [TEAL-929](https://jira.skatelescope.org/browse/TEAL-929)
- [TEAL-1045](https://jira.skatelescope.org/browse/TEAL-1045)

## Contributors

- [Burba, Jacob](https://github.com/jburba)
