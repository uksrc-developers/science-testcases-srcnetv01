# Profiling workflows extracted from test notebooks with psrecord

**Documentation on Confluence:** https://confluence.skatelescope.org/x/dWf7Ew

**Summary:** The scripts in this sub-repository allow one to run `psrecord` to profile python scripts adapted from the SRCNet v0.1 test campaign and described [here](https://confluence.skatelescope.org/display/SRCSC/v0.1+Test+campaign+-+Test+descriptions) and [here](https://confluence.skatelescope.org/display/SRCSC/Framework+Reporting). The instructions to perform these profiling tasks can be found, with more details [here](https://confluence.skatelescope.org/x/dWf7Ew).

## Dependencies

Running the profiling script requires one to have installed `psrecord` (see [profiling instructions](https://confluence.skatelescope.org/x/dWf7Ew).
Additional dependencies might be required, depending on which test you want to profile (see [test descriptions](https://confluence.skatelescope.org/x/LG61Ew)).


## Instructions

1. If needed, install `psrecord` (see [profiling instructions](https://confluence.skatelescope.org/x/dWf7Ew).
2. If you have not done it yet, clone the Git repository which contains the tests and download the required data (see [testing instructions](https://confluence.skatelescope.org/x/NW61Ew) if needed).
3. Update the paths in the configuration file `config/config.yaml` inside the test repo you have cloned.
        `data_path` is the absolute path to the data you have just downloaded or staged.
        `result_path` is the absolute path to the directory where you want to store the outputs of both the tests and their profiling
        `test_path` is the absolute path to the test-case repo you have (just) cloned.
        `psrecord_path` is the path to your `.local/bin/` directory to add to your path and use `psrecord` from the command line.
4. Run the profiling script by typing in your command line: `./profiling_script.sh test_name` where `test_name` is the code describing the test you want to profile (e.g., SWF-002-T1). Note that you might need to run `chmod u+x profiling_script.sh` to make the script executable and you must have `psrecord` installed, with python>=3.7.
5. Read the outputs from `result_path/test_name/profiling/`: the logfile summarises the profiling results in a format that is easy to input in the [reporting form](https://forms.gle/ZQwsdQNxTnsC1fRG6) and image files are ready to be uploaded in the form.

Note that `psrecord` cannot profile I/O on Mac OSX -- you need to remove the `--include-io` in the profiling script and all I/O-related outputs from the `monitor_activity` script to run the profiling workflow on Mac OSX.

## SRCNet Links

### Confluence

- [Guidelines for running psrecord](https://confluence.skatelescope.org/x/dWf7Ew)
- [Framework reporting](https://confluence.skatelescope.org/display/SRCSC/Framework+Reporting)

### Jira

Features
- [SP-5622](https://jira.skatelescope.org/browse/SP-5158)
Tickets
- [TEAL-1058](https://jira.skatelescope.org/browse/TEAL-1058)

## Contributors

- Adélie Gorce
- TEAL team