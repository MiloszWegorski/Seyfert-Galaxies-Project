import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import CCDData
import ccdproc


def saveFitsFile(file, Name):
    """
    Saves FITS data as a FITS file.
    If a file with a certain name already exists this will save the file as the
    file name with _rev(revision number) at the end of it

    Parameters
    ----------
    file : FITS Data
        Data of the FITS file to be saved.
    Name : String
        Destination and name of the file to be saved.

    Raises
    ------
    e
        In case another exception occurs of unknown origin ¯\_(ツ)_/¯.

    Returns
    -------
    none.

    """
    NewName = Name
    while True:
        break
        with open(NewName + '.FITS', 'a'):    
            file.write(NewName + '.FITS', overwrite='True')
        break 

def MasterBias():
    """
    gets master Bias from files

    Returns
    -------
    FITS
        Master Bias.

    """
    bias_dir = 'DBF_night1/Bias'
    biases = sorted(glob.glob(bias_dir+'/*.FITS'))
    
    totalbias = []
    
    for i in biases:
        totalbias.append(CCDData(fits.getdata(i), unit='adu'))
    return ccdproc.combine(totalbias, bias_dir+'/Masterbias.FITS', 'median')

def MasterFlat():
    """
    Combines all flats into one master flat

    Returns
    -------
    Saves master Flat into file.

    """
    flat_dir = ['DBF_night1/Flat/B','DBF_night1/Flat/Ha',
    'DBF_night1/Flat/O3','DBF_night1/Flat/S2','DBF_night1/Flat/U',
    'DBF_night1/Flat/V']
    
    flats_lib = []
    
    for i in range(len(flat_dir)):
        
        flats_lib.append(sorted(glob.glob(flat_dir[i]+'/*.FITS')))
    
    totalFlats = [[],[],[],[],[],[]]
    
    MasterFlats = [6]
    
    for i in range(len(flats_lib)):
        for j in range(len(flats_lib[i])):
            #total flats is array of arrays of results indexes being for
            #different filters
            
            #i is array of file names of flats for each filter
            #j is used to
            totalFlats[i].append(CCDData(fits.getdata(flats_lib[i][j]),
                                         unit='adu'))
    
    MasterFlats[0] = ccdproc.combine(totalFlats[0],'DBF_night1/Flat/\
FlatMasterFlat_B.FITS', 'median')
    MasterFlats[1] = ccdproc.combine(totalFlats[1],'DBF_night1/Flat/\
FlatMasterFlat_Ha.FITS', 'median')
    MasterFlats[2] = ccdproc.combine(totalFlats[2],'DBF_night1/Flat/\
FlatMasterFlat_O3.FITS', 'median')
    MasterFlats[3] = ccdproc.combine(totalFlats[3],'DBF_night1/Flat/\
FlatMasterFlat_S2.FITS', 'median')
    MasterFlats[4] = ccdproc.combine(totalFlats[4],'DBF_night1/Flat/\
FlatMasterFlat_U.FITS', 'median')
    MasterFlats[5] = ccdproc.combine(totalFlats[5],'DBF_night1/Flat/\
FlatMasterFlat_V.FITS', 'median')

def getMasterFlat():
    """
    Returns the master flat file for all filters
    order of files = B, Ha, O3, S2, U, V

    Parameters
    ----------
    fil : TYPE
        DESCRIPTION.

    Returns
    -------
    FITS files of flat field images.

    """
    raw_fil = glob.glob(
       'C:/Users/mwegorski/Desktop/Astro_Project/Astro_Project/DBF_night1/Flat' 
                       + '**/*.FITS')
    
    fin_flats= []
    
    for i in raw_fil:
        fin_flats.append(fits.getdata(i))
    
    return fin_flats

def getMasterBias():
    """
    Gets master bias from files

    Returns
    -------
    FITS
        DESCRIPTIOmaster bias fits image.

    """
    raw_file = glob.glob('../DBF_night1/Bias/Masterbias.FITS')
    
    return fits.getdata(raw_file[0], unit='adu')

def getRawData():
    """
    Returns the raw images of the galaxy with the filters V, B*2, S2, Ha in 
    that order

    Returns
    -------
    fin_Raw : TYPE
        Raw fits images.

    """
    raw_fil = glob.glob(
        'C:/Users/mwegorski/Desktop/\
        Astro_Project/Astro_Project/NGC1275'+ '**/*.FITS')
    
    fin_Raw= []
    
    for i in raw_fil:
        fin_Raw.append(fits.getdata(i))
    
    return fin_Raw

def processFinalIMG(MasterBias, MasterFlat, Images):
    """
    Corrects a set of raw FITS images using a MasterBias and Master Flat

    Parameters
    ----------
    MasterBias : Fits Image
        DESCRIPTION.
    MasterFlat : Fits Image
        DESCRIPTION.
    Images :  lisi
        list of FITS images.

    Returns
    -------
    list
        List of Processed FITS images.

    """
    #order to process = V, B , b, S2, Ha
    
    for i in range(len(Images)):
        Images[i] = CCDData(Images[i], unit='adu')
        
    #getting master flats
    MasterFlat_V = CCDData(MasterFlat[5], unit='adu', 
                           mask=np.zeros(MasterFlat[5].shape))
    MasterFlat_B = CCDData(MasterFlat[0], unit='adu', 
                           mask=np.zeros(MasterFlat[0].shape))
    MasterFlat_S2 = CCDData(MasterFlat[3], unit='adu', 
                           mask=np.zeros(MasterFlat[3].shape))
    MasterFlat_Ha = CCDData(MasterFlat[1], unit='adu', 
                           mask=np.zeros(MasterFlat[1].shape))
    
    masterBias = CCDData(MasterBias, unit='adu')
    
    #bias substraction
    MasterFlat_V = MasterFlat_V.subtract(masterBias)
    MasterFlat_B = MasterFlat_B.subtract(masterBias)
    MasterFlat_S2 = MasterFlat_S2.subtract(masterBias)
    MasterFlat_Ha = MasterFlat_Ha.subtract(masterBias)
    
    for i in range(len(Images)):
        Images[i] = Images[i].subtract(masterBias)
    
    #Flat field correction
    final_image1 = ccdproc.flat_correct(Images[0], MasterFlat_V)
    final_image2 = ccdproc.flat_correct(Images[1], MasterFlat_B)
    final_image3 = ccdproc.flat_correct(Images[2], MasterFlat_B)
    final_image4 = ccdproc.flat_correct(Images[3], MasterFlat_S2)
    #final_image5 = ccdproc.flat_correct(Images[4], MasterFlat_Ha)
    
    final_comb_b = ccdproc.Combiner([final_image2, final_image3])
    
    saveFitsFile(final_image1, 'Final_filter_V')
    saveFitsFile(final_image2, 'Final_filter_B')
    saveFitsFile(final_image3, 'Final_filter_B_Image2')
    saveFitsFile(final_image4, 'Final_filter_S2')
    saveFitsFile(final_comb_b, 'Final_Combined_B')
    
    return [final_image1, final_image2, final_image3, final_image4]

#process image
fin_img = processFinalIMG(getMasterBias(), getMasterFlat(), getRawData())

#plot image to check if processing is done correctly
fig, axes = plt.subplots(2, 2)

axes[0,0].imshow(np.log(fin_img[0])
           , cmap='gray',vmin=5.3, vmax=np.log(900))
#axes[0,0].colorbar()

axes[0,1].imshow(np.log(fin_img[1]), cmap='gray',
           vmin=np.log(50), vmax=np.log(900), label='')
#axes[0,1].colorbar()

axes[1,0].imshow(np.log(fin_img[2]), cmap='gray',
           vmin=np.log(110), vmax=np.log(900))
#axes[1,0].colorbar()
axes[1,1].imshow(np.log(fin_img[3]), cmap='gray',
           vmin=np.log(50), vmax=np.log(900))
#axes[1,1].colorbar()