# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 13:56:08 2024

@author: mwegorski
"""

from photutils.aperture import CircularAperture, ApertureStats

import numpy as np

def Circle_Astrometry(position, radius, data):
    """
    Creates a set number of apertures given a set of positions and radii for 
    said apertures

    Parameters
    ----------
    position : LIST
        list of coordinates for apertures.
    radius : LIST
        List of radii for the .
    data : FITS
        fits data of image the apertures and fluxes to be calculated of.

    Returns
    -------
    apertures : LIST
        List of aperture objects.
    phot_table : TABLE
        Table of photometry data.
    fluxes : LIST
        List of fluxes for the corresponding apertures.

    """
    
    #create appeatures for each coordinate and radius in arrays
    apertures = [CircularAperture(a, r) for a, r in zip(position, radius)]
    
    #list to store all fluxes from apertures
    fluxes = []
    
    #calculate fluxs for all apertures
    for i in apertures:
        phot_table = ApertureStats(data, i)
        fluxes.append(phot_table.sum)
    
    return apertures, phot_table, fluxes
    

def flux_to_intens(flux1, flux2, intens_standard):
    """
    Converts flux to intensity given the flux of the object which the intensity
    is to be calculated of and the flux of a known standard star and its 
    magnitude

    Parameters
    ----------
    flux1 : Float
        Flux of known standard object.
    flux2 : Float
        Flux of unknown magnitude object.
    intens_standard : Float
        magnitude of standard object.

    Returns
    -------
    intens : Float
        magnitude of unknown object.

    """
    d_intens = -2.5* np.log10(flux1/flux2)
    intens = intens_standard + d_intens
    
    return intens