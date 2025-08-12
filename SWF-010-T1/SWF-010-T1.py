#!/usr/bin/env python
# coding: utf-8
 
# SWF-010-T1 – Radio continuum source finding with PyBDSF
# Authors: Lara Alegre, Adélie Gorce, and the Teal team 
# See README for details.

# Imports 
import os
import yaml
import warnings
import bdsf
import sys
import numpy as np 
import pandas as pd
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, ZScaleInterval, AsinhStretch, SqrtStretch

logger = logging.getLogger(__name__)

def main():

    # set up paths        
    config_file = './config/config.yml'
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)
    else:
        warnings.warn(f"! Configuration file '{config_file}' not found. Using default settings.", UserWarning)
        config = {}
        
    default_base_path = "../datasets"
    base_path = config.get("data_path")
    if base_path is None:
        warnings.warn(f"! 'data_path' not found in '{config_file}'. Using default path '{default_base_path}'.", UserWarning)
        base_path = default_base_path
    datafolder = base_path +"/" # where to find the downloaded data - here, on Azimuth
    
    default_result_path = "../results"
    result_path = config.get("result_path")
    if result_path is None:
        warnings.warn(f"! 'result_path' not found in '{config_file}'. Using default: '{default_result_path}'.", UserWarning)
        result_path = default_result_path
    result_path += "/SWF-010-T1/"
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    # Add path to save any plots and tables
    save_path = result_path 
    
    # Path to the folders and FITS file (inside the container in this case)
    fits_path = base_path + '/P020_39-mosaic-blanked.fits'
    # Path to the output directory
    output_dir = result_path + '/interim'
    # Path to results
    results_dir = result_path
    
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set up log file
    if not os.path.exists(f'{result_path}profiling/'):
        os.makedirs(f'{result_path}profiling/')
    logging.basicConfig(
        filename=f'{result_path}profiling/SWF-010-T1.log',
        filemode='w',
        level=logging.INFO,
        format="%(message)s",
    )
    logger.info('\nStarted running the workflow.')

    # Run PyBDSF
    
    # Process the patched file with PyBDSF (this takes a about 5 minutes) 
    # Got parameters from https://github.com/mhardcastle/ddf-pipeline/blob/5d11d5ec212da29826433fed44bab329c037855e/scripts/sourcefind.py
    
    # define frequency since the FITS file does not have frequency information and it is needed to run PyBDSF
    restfrq=143650000.0 #144000000.0
    
    img = bdsf.process_image(fits_path, thresh_isl=3.0, thresh_pix=5, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)
    
    # Confirm success without matplotlib
    logger.info("Image shape:", img.ch0_arr.shape)
    logger.info("Beam:", img.beam)
    # Number of detected sources
    logger.info("Sources (regions of radio emission):", len(img.sources))
    logger.info("Fitted Gaussians:", len(img.gaussians))

    # Save source to FITS format 
    outfile_srl=os.path.join(output_dir, 'SWF-010-T1_source_catalogue.fits')
    img.write_catalog(outfile=outfile_srl, format='fits',  catalog_type='srl', clobber = True)
    # Save Gaussian catalogue to FITS format 
    outfile_gaul=os.path.join(output_dir, 'SWF-010-T1_gaussian_catalogue.fits')
    img.write_catalog(outfile=outfile_gaul, format='fits',  catalog_type='gaul', clobber = True)
     
    # # Plotting the results 

    # We will plot the results over the original radio map
    with fits.open(fits_path) as hdul:
        data = hdul[0].data.squeeze()
        header = hdul[0].header
        wcs = WCS(header)
    
    # Set up WCS plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=wcs)
    # Normalise the image
    norm = ImageNormalize(data*1000, interval=ZScaleInterval(), stretch=SqrtStretch())
    im = ax.imshow(data*1000, origin='lower', norm=norm)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.set_title("Original LoTSS mosaic")
    plt.colorbar(im, ax=ax, pad=0.05, shrink=0.7, aspect=20, label='Intensity [mJy/beam]')
    output_mosaic = os.path.join(results_dir, 'SWF-010-T1_LoTSS_mosaic_radec.png')
    plt.savefig(output_mosaic, dpi=300)
    logger.info(f'{output_mosaic} figure saved.')

    # Read catalogues that contain PyBDSF information
    
    # Source catalogue (SRL)
    srl_table = Table.read(outfile_srl)
    srl_df = srl_table.to_pandas()
    logger.info("Columns in the source catalogue:", srl_df.columns.values)    
    # Gaussian catalogue (GAUS)
    gaus_table = Table.read(outfile_gaul)
    gaus_df = gaus_table.to_pandas()
    
    # Zoom in into a region and plot PyBDSF source and Gaussian positions and ellipses 
    
    # --- CUTOUT SETUP ---
    cutout_size = 512 #can define a wider area 
    
    ny, nx = data.shape
    
    # Can shift values (in pixels)
    shift_x = -400   # positive = right, negative = left
    shift_y = -1100  # positive = up, negative = down
    
    # Calculate new central pixel with shift
    yc = ny // 2 + shift_y
    xc = nx // 2 + shift_x
    
    half = cutout_size // 2
    y1, y2 = yc - half, yc + half
    x1, x2 = xc - half, xc + half
    
    data_cut = data[y1:y2, x1:x2]
    wcs_cut = wcs.slice((slice(y1, y2), slice(x1, x2)))
    
    # --- RA/DEC BOUNDING BOX ---
    cutout_corners = [(x1, y1), (x2, y2)]
    cutout_world_coords = wcs.wcs_pix2world(cutout_corners, 0)
    ra_min, dec_min = cutout_world_coords[0]
    ra_max, dec_max = cutout_world_coords[1]
    
    ra_low = min(ra_min, ra_max)
    ra_high = max(ra_min, ra_max)
    dec_low = min(dec_min, dec_max)
    dec_high = max(dec_min, dec_max)
    
    # --- FILTER CATALOGUESS ---
    srl_cut = srl_df[
        (srl_df['RA'] >= ra_low) & (srl_df['RA'] <= ra_high) &
        (srl_df['DEC'] >= dec_low) & (srl_df['DEC'] <= dec_high)
    ]
    
    gaus_cut = gaus_df[
        (gaus_df['RA'] >= ra_low) & (gaus_df['RA'] <= ra_high) &
        (gaus_df['DEC'] >= dec_low) & (gaus_df['DEC'] <= dec_high)
    ]
    
    # --- PLOT ---
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': wcs_cut})
    
    norm = ImageNormalize(data_cut*1000, interval=ZScaleInterval(), stretch=SqrtStretch())
    im = ax.imshow(data_cut*1000, origin='lower', norm=norm)
    
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.set_title("Cutout of LoTSS mosaic with PyBDSF source and Gaussian positions and ellipses", pad=10)
    
    # Add colorbar 
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.05, shrink=0.7, aspect=20)
    cbar.set_label('Intensity [mJy/beam]')
    
    
    # --- OVERLAY SOURCE POSITIONS & ELLIPSES ---
    for _, row in srl_cut.iterrows():
        ra, dec = row['RA'], row['DEC']
        maj, min_, pa = row['Maj']*3600, row['Min']*3600, row['PA']
    
        # Convert RA/DEC to pixel coordinates
        x, y = wcs_cut.wcs_world2pix(ra, dec, 0)
    
        # Convert FWHM from arcsec → pixels
        arcsec_per_pixel = np.abs(header['CDELT1']) * 3600  # degrees → arcsec
        
        width = maj / arcsec_per_pixel
        height = min_ / arcsec_per_pixel
        #logger.info(f"Maj={maj}, RA={ra}, Dec={dec}, width={width:.2f}, height={height:.2f}, PA={pa}")
        # Draw the ellipse
        ell = Ellipse((x, y), width=width, height=height, angle=pa,
                      edgecolor='red', facecolor='none', lw=1)
        ax.add_patch(ell)
        # mark the center
        ax.plot(x, y, 'rx', markersize=4, label='Sources')
    
    # --- OVERLAY GAUSSIANS ---
    
    for _, row in gaus_cut.iterrows():
        ra, dec = row['RA'], row['DEC']
        maj, min_, pa = row['Maj']*3600, row['Min']*3600, row['PA']
        x, y = wcs_cut.wcs_world2pix(ra, dec, 0)
    
        width = maj / arcsec_per_pixel
        height = min_ / arcsec_per_pixel
    
        ell = Ellipse((x, y), width=width, height=height, angle=pa,
                     edgecolor='blue', facecolor='none', lw=1)
        ax.plot(x, y, 'x', color ='blue',  markersize=4, label='Gaussians')
        ax.add_patch(ell)
    
    
    legend_elements = [
        Line2D([0], [0], color='red', lw=2, label='Sources'),
        Line2D([0], [0], color='blue', lw=2, label='Gaussians')
    ]
    
    ax.legend(handles=legend_elements, loc='upper right')
    output_results = os.path.join(results_dir, 'SWF-010-T1_LoTSS_zoomed_in_validation_file.png')
    plt.savefig(output_results, dpi=300)
    logger.info(f'{output_results} figure saved.')

    logger.info('Finished.')

if __name__ == '__main__':
    main()