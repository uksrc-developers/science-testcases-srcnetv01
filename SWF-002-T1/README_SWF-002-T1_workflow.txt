*SWF-002-T1: Positional Cross-Match Test*

**Test Name:** SWF-002-T1: Positional cross-matching on multi-wave data

**Authors:** Bonny Barkus, Adélie Gorce, and Teal team

**Confluence description:** https://confluence.skatelescope.org/pages/viewpage.action?pageId=319985104

**Usage**

To run the script, type `python SWF-002-T1_workflow.py my_folder"`, where `my_folder` is the folder that includes the `teal/` folder that hosts the data necessary to run these tests (the data must have been staged beforehand, if running on a specific node).

**Summary** 
This is a short script to run a simple positional cross-match between two catalogues, and links to the following confluence page ([Positional cross-matching on multi-wave data test](https://confluence.skatelescope.org/pages/viewpage.action?pageId=319985104)). It is based on the notebook with the same name. 
This example uses the LOFAR Virgo Cluster Survey data set for the radio catalogue ([Virgo Cluster Survey](https://lofar-surveys.org/virgo_data.html)), and PanSTARRS ([PanSTARRS home page](https://outerspace.stsci.edu/display/PANSTARRS/)) for the optical data. The user is welcome to use any catalogues they have to hand.

**Description:** 
This notebook uses `astropy SkyCoord`'s `match_to_catalog_sky` to find the closest on-sky sources by RA and DEC and filter this down to those matches which are less than *1.5 arcsecs* in distance. The notebook is broken down into the following sections: 
- Imports
    - This cell contains all the imports required.
- Load in files
    - This section loads in the files that are stored in the data system. For this example, the LOFAR Virgo Cluster survey catalogue is available, and a small filter PanSTARRs catalogue has been downloaded. If the PanSTARRS catalogue is not available, this example requires the use of `astroquery` to run a search through PanSTARRS around the given area. Be aware **this will max at the first 50 rows** and therefore may not produce may matches inside the given area. The area for the radio catalogue has been chosen as a 2 degree circle at the centre of the survey. The PanSTARRS data, is selected to be slightly larger to allow for boundary effects. If the user choses their own catalogues the appropriate paths ways or astroquery must be used at this stage.
- Catalogue information check and setup
    - This is a quick check on what the columns and units of the catalogues are, and follows by setting the column names for the rest of the notebook. If the user chooses to use different catalogues the column names need to be changed at this location.
- Sky area plots to check coverage
    - These plots will check the coverage of the catalogues. This is important to make sure that the loaded data overlaps, and that the filtered radio catalogue covers a slightly smaller area than optical catalogue. In the example these plots should be as on the confluence page. In the case the user has their own data, it is a check that the catalogues cover the same sky area. These can be saved.
- Positional cross match
    - This section does the cross-matching using the RA and DEC positions in each catalogue. It find the closest on-sky optical source for the filtered radio catalogue. These matches are then filtered to only those counterparts which are within 1.5 arcsecs of each other. A new table of the matches is made, comtaining all the columns from both the radio and optical catalogues, and a column containing the on-sky separation between the two sources in arcseconds.
- Plot the results
    - The results are visualised in two plots. The first plot is a bar chart showing the number of matching pairs which have an on-sky separation within in each 0.1 arcsec bin. The second shows the position of the matched radio source, and these are shaded according to their on-sky separation from their matched optical source. In this plot, the darker the radio source position the further the optical match is.
    