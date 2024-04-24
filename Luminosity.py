# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:32:27 2024

@author: mwegorski
"""

import glob
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from photutils.isophote import EllipseGeometry, build_ellipse_model, Ellipse
from photutils.aperture import EllipticalAperture, CircularAperture
from astropy.stats import sigma_clipped_stats


from calc_ellipse import Circle_Astrometry,  flux_to_intens

#Standard objects of known intensity
#TYC 2856-334-1
Int_TYC = 12.32
#Pul -3 270377
Int_PUL = 13.3

#fontsize for plots
plt.rcParams['font.size'] = 17

#set image name
image = 'Final_filter_V.FITS'

    

#get image and convert into FITS format
data = glob.glob(image)
data = fits.getdata(data[0], unit='adu')

#calculating the mean value
_, mean, _ = sigma_clipped_stats(data)
print('Mean of the image = ', mean)

#ensuring no negative values in Fits file from processing
for i in range(len(data)):
    for j in range(len(data[i])):
        if data[i][j] < 0:
            data[i][j] = 0
#--------------------------------Standard stars-------------------------------#

TYC_2856_Name = "TYC 2856-334-1"
TYC_2856_V_Int = 12.32
TYC_2856_B_Int = 12.84

GPM_5041_name = "GPM 50.092290+41.539668"
GPM_5041_V_Int = 11.9
GPM_5041_B_Int = 12.2

GPM_4941_name = "GPM 49.949861+41.429504"
GPM_4941_V_Int = 13.1
GPM_4941_B_Int = 13.8

GPM_4941_2_name = "GPM 49.963274+41.523785"
GPM_4941_2_V_Int = 13.76

#-----------------------------------------------------------------------------#

#define mask aperture
mask = CircularAperture((576, 421), 4.0)
mask2 = CircularAperture((542, 422), 4.0)

#convert aperture object to mask object
masks = mask.to_mask(method='exact')
masks2 = mask2.to_mask(method='exact')

#multiply mask by data to create image of object
data_weighted = masks.multiply(data)
data_weighted2 = masks2.multiply(data)

#coordinates of respective corners of mask on image
top_left_x = 570
top_left_y = 419

top_left_x_2 = 536
top_left_y_2 = 419
#substract masks from data
for i in range(9):
    for j in range(9):
        data[top_left_y+j][top_left_x+i] -= (data_weighted[i][j])/2
        data[top_left_y_2+j][top_left_x_2+i] -= (data_weighted2[i][j])/2


#-----------------------------------------------------------------------------#

#positions of standard stars
positions = [(548.5, 259.5), #12.32 TYC 2856-334-1
             (416, 357.5), #Pul -3 270377 13.3 B BAND
             (452, 521), #bg count
             (53, 292.5), #11.9 V GPM 50.092290+41.539668
             (556, 810), #13.1 V band GPM 49.949861+41.429504
             (505, 368), #13.76 V GPM 49.963274+41.523785   
             (551, 424)
             ]

#radii which best fit the standard satrs
radii = [8, 7, 12, 7.5, 6, 6, 35]

#get all fluxes and create aperture objects of standard stars
apers, isolist, fluxes = Circle_Astrometry(positions, radii, data)

#-----------------------------------------------------------------------------#
#get ellipse to fit galaxy
geometry = EllipseGeometry(551, 424, 30, 0.3, 130 * np.pi / 180.0)

#create aperture from said ellipse
aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                          geometry.sma * (1 - geometry.eps),
                          geometry.pa)

#obtaining data from ellipse
ellipse = Ellipse(data, geometry)
isolist_NGC1275 = ellipse.fit_image()

#create ellipse model
model_image = build_ellipse_model(data.shape, isolist_NGC1275)
residual = data - model_image

#-----------------------------------------------------------------------------#

#create circle over star to calculate HWHM
geometry_see = EllipseGeometry(548.5, 259.5, 11, 0, 130 * np.pi / 180.0)

#create aperture object
aper_see = EllipticalAperture((geometry_see.x0, geometry_see.y0), 
                              geometry_see.sma, 
                              geometry_see.sma * (1 - geometry_see.eps),
                              geometry_see.pa)

#get data from ellipse
ellipse_see = Ellipse(data, geometry_see)
seeing_star = ellipse_see.fit_image()

model_image = build_ellipse_model(data.shape, seeing_star)

#plot the brightness curve with respect to radius (here radius is equivalent to
#the SMA because its a circle)
plt.figure()
plt.plot(seeing_star.sma[0:40], seeing_star.intens[0:40])
plt.axhline(y=(seeing_star.intens[0]/2), color= 'red',xmin=0, xmax = 2.09/6.5, 
            linestyle='--')
plt.axvline(x=2.09,ymin=0, ymax=1/2, color='red', linestyle='--')
plt.xlabel('Radius (Pixels)')
plt.ylabel('Flux')
plt.xlim(0, 6.5)
plt.grid()

#-----------------------------------------------------------------------------#

#plot for apertures
plt.figure()
#names and separate colors for each aperture
names = ['12.32 TYC 2856-334-1', 'Pul -3 270377', 'Background', 
         'GPM 50.092290+41.539668', 'GPM 49.949861+41.429504', 
         'GPM 49.963274+41.523785 ']
colors = ['r', 'g', 'b', 'y', 'purple']

#plot image
ax = plt.subplot()
ax.imshow(data, cmap='gray', vmin=20, vmax= 1000 )

#plot the seeing calculation star aperture
#aper_see.plot(color='green')

#plot all apertures
for i, j, col in zip(apers, names, colors):
    i.plot(label=j, color=col)
    aper.plot(color='green')

#plot masks
mask.plot(color= 'red')
mask2.plot(color= 'red')
results = []

#plot titles and legend
ax.set_ylabel('Declination (pix)')
ax.set_xlabel('Right Ascention (pix)')
#ax.legend()

#plot for magnitude vs SMA of NGC1275
plt.figure()

#calculating data to be plotted
for i in isolist_NGC1275.intens:
    results.append(flux_to_intens( i/(0.76**2), fluxes[0], TYC_2856_B_Int))

#Saving data to text file
#Uncomment if intens vs SMA data needs to be saved to text file

Save_intens = open("Intens_vs_SMA.txt", 'w')
Save_intens.write(str(list(isolist_NGC1275.sma*(0.76))))
Save_intens.write('\n')
Save_intens.write(str(results))
Save_intens.close()

#plot the magnitude vs SMA of NGC 1275
plt.scatter((isolist_NGC1275.sma*(0.76**2))[0:-1], results[0:-1])
plt.grid()
plt.legend()

#plot labels and legend
plt.ylabel("Apparent Magnitude")
plt.xlabel("semi-major axis (arcsec)")
plt.gca().invert_yaxis()
plt.show()

def plots():
    
    #plot of ellipticity vs SMA
    plt.figure()
    plt.errorbar(isolist_NGC1275.sma*0.76, isolist_NGC1275.eps, yerr=isolist_NGC1275.ellip_err,
                 fmt='o', markersize=4)
    plt.xlabel('Semimajor Axis Length (arcsec)')
    plt.ylabel('Ellipticity')
    plt.grid()
    plt.axvline(x=3, color='red', linestyle='--')
    plt.axvline(x=9.6, color='red', linestyle='--')
    
    #Plot ofposition angle VS SMA
    plt.figure()
    plt.errorbar(isolist_NGC1275.sma*0.76, isolist_NGC1275.pa / np.pi * 180.0,
                 yerr=isolist_NGC1275.pa_err / np.pi * 80.0, fmt='o', markersize=4)
    plt.xlabel('Semimajor Axis Length (arcsec)')
    plt.ylabel('PA (deg)')
    plt.grid()
    plt.axvline(x=6.8, color='red', linestyle='--')
    plt.axvline(x=19, color='red', linestyle='--')
    
    #x0 position of galaxy vs sma
    plt.figure()
    plt.errorbar(isolist_NGC1275.sma*0.76, isolist_NGC1275.x0, yerr=isolist_NGC1275.x0_err, fmt='o',
                 markersize=4)
    plt.xlabel('Semimajor Axis Length (arcsec)')
    plt.ylabel('x0')
    plt.grid()
    plt.axvline(x=6, color='red', linestyle='--')
    plt.axvline(x=18, color='red', linestyle='--')
    
    #y0 position of galaxy vs sma
    plt.figure()
    plt.errorbar(isolist_NGC1275.sma*0.76, isolist_NGC1275.y0, yerr=isolist_NGC1275.y0_err, fmt='o',
                 markersize=4)
    plt.xlabel('Semimajor Axis Length (arcsec)')
    plt.ylabel('y0')
    plt.grid()
    plt.axvline(x=4, color='red', linestyle='--')
    plt.axvline(x=10, color='red', linestyle='--')
    
    #plots to plot Ellipse model, data and residual
    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(14, 5), nrows=1, ncols=3)
    fig.subplots_adjust(left=0.04, right=0.98, bottom=0.02, top=0.98)
    
    #plotting galaxy data
    ax1.imshow(data, origin='lower', cmap='grey', vmin=10, vmax= 2000)
    ax1.set_title('Data')
    ax1.set_ylabel('Declination (pix)')

    smas = np.linspace(10, 50, 5)
    for sma in smas:
        iso = isolist_NGC1275.get_closest(sma)
        x, y, = iso.sampled_coordinates()
        ax1.plot(x, y, color='white')
    
    #plotting image and ellipse model on top of it
    ax2.imshow(model_image, origin='lower', cmap='grey', vmin=20, vmax= 1000 )
    ax2.set_title('Ellipse Model')
    ax2.set_xlabel('Right Ascention (pix)')
    
    #plotting residual of image
    ax3.imshow(residual, origin='lower', cmap='grey', vmin=20, vmax= 1000 )
    ax3.set_title('Residual')

plots()