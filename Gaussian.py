import numpy as np
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#------------------------------Global-Variables-------------------------------#
gauss_bar = []
gauss_disk = []
gauss_bulge = []

sums = []

#---------------------------Data for processing-------------------------------#
#get data from file which contains 3 arrays in separate lines in order of 
#x vals, yvals, yvals 2
data = open('Intens_vs_SMA.txt', 'r')
arrs = []
temp = data.readlines()

for i in temp:
    arrs.append(i.split(', '))
    
    
for i in arrs:
    i[0] = i[0][1:-1]
    i[-1] = i[-1][0:-2]
    for j in range(len(i)):
        i[j] = float(i[j])

xvals = np.linspace(arrs[0][0], arrs[0][-1], 1000)


plt.rcParams['font.size'] = 21
def func2(r, I_e1, r_e1, n1, I_e2, r_e2, n2, I_e3, r_e3):
    """
    Function to calculate Magnitude of seircic parameter for a three parameter
    fit

    Parameters
    ----------
    r : Float
        SMA values of galaxy (x axis values).
    I_e(1,2,3) : Float
        Seirsic parameter I_e.
    r_e(1,2,3) : Float
        Seirsic paramete R_e.
    n(1,2) : Float
        Seirsic parameter n.

    Returns
    -------
    Float
        Seirsic magnitude for 3 component fit.

    """
    B_n1 = 1.997*n1-0.327
    B_n2 = 1.997*n2-0.327
    B_n3 = 1.997-0.327
    
    comp1 = -2.5* np.log10(I_e1*np.exp(-B_n1*(((r/r_e1)**(1/n1))-1)))
    comp2 = -2.5* np.log10(I_e2*np.exp(-B_n2*(((r/r_e2)**(1/n2))-1)))
    comp3 = -2.5* np.log10(I_e3*np.exp(-B_n3*(((r/r_e3)**(1))-1)))
    
    tot = 10**(-comp1) + 10**(-comp2) + 10**(-comp3)
    
    return -np.log10(tot)

def func3(r, I_e, r_e, n):
    """
    Function to calculate Magnitude of seircic parameter for a single parameter
    fit

    Parameters
    ----------
    r : Float
        SMA values of galaxy (x axis values).
    I_e : Float
        Seirsic parameter I_e.
    r_e : Float
        Seirsic paramete R_e.
    n : Float
        Seirsic parameter n.

    Returns
    -------
    Float
        Seirsic magnitude for one component fit.

    """
    B_n = 1.997*n-0.327
    return -2.5* np.log10(I_e*np.exp(-B_n*(((r/r_e)**(1/n))-1)))

def func(I_e, B_n, r_e, n, r):
    """
    calculates the resulting intensity of one seirsic parameter

    Parameters
    ----------
    I_e : Float
        Seirsic parameter I_e.
    B_n : Float
        Seirsic parameter B_n.
    r_e : Float
        Seirsic paramete R_e.
    n : Float
        Seirsic parameter n.
    r : Float
        SMA values of galaxy (x axis values).

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return I_e*np.exp(-B_n*(((r/r_e)**(1/n))-1))

def update_bar(val):
    """
    Updates the Y values of the bar parameter

    Parameters
    ----------
    val : method
        to obtain the values bar the sliders.

    Returns
    -------
    Updates graph of fit and bar parameter.

    """
    # Convert back from logarithmic form
    I_e = Ie_slider_bar.val  
    B_val = (1.997*N_slider_bar.val)-(0.327)
   
    gaussian_bar.set_ydata([-2.5*np.log10(func(I_e, B_val, Re_slider_bar.val,
                                               N_slider_bar.val, i))
                            for i in xvals])
    
    
    sum_profiles(gaussian_bar, gaussian_disk, gaussian_bulge)
    fig.canvas.draw_idle()

def update_bulge(val):
    """
    Updates the Y values of the Bulge parameter

    Parameters
    ----------
    val : method
        to obtain the values from the sliders.

    Returns
    -------
    Updates graph of fit and Bulge parameter.

    """
    # Convert back from logarithmic form
    I_e = Ie_slider_bulge.val  
    B_val = (1.997*N_slider_bulge.val)-(0.327)
    
    gaussian_bulge.set_ydata([-2.5*np.log10(func(I_e, B_val, 
                                                 Re_slider_bulge.val, 
                                                 N_slider_bulge.val, i)) 
                                                               for i in xvals])
    
    sum_profiles(gaussian_bar, gaussian_disk, gaussian_bulge)
    fig.canvas.draw_idle()

def update_disk(val):
    """
    Updates the Y values of the disk parameter

    Parameters
    ----------
    val : method
        to obtain the values from the sliders.

    Returns
    -------
    Updates graph of fit and disk parameter.

    """
    # Convert back from logarithmic form
    I_e = Ie_slider_disk.val  
    B_val = (1.997*N_slider_disk.val)-(0.327)
    
    gaussian_disk.set_ydata([-2.5*np.log10(func(I_e, B_val, Re_slider_disk.val
                                                , N_slider_disk.val, i))
                                                 for i in xvals])
    
    sum_profiles(gaussian_bar, gaussian_disk, gaussian_bulge)
    fig.canvas.draw_idle()

def sum_profiles(Bar, Disk, Bulge):
    """
    Sums the 3 profiles as magnitudes by converting out of logarythmic form and
    them converting back to it

    Parameters
    ----------
    Bar : plt.plot line
        bar parameter for Sersic fit.
    Disk : plt.plot line
        Disk  parameter for Sersic fit.
    Bulge : plt.plot line
        Bulge  parameter for Sersic fit.

    Returns
    -------
    Sets Sum_mags on the plot to be the sum of the 3 parameters.

    """
    result = []
    
    for Br, D, Bu in zip(Bar._y, Disk._y, Bulge._y):
        result.append(10**(-Br) + 10**(-D) + 10**(-Bu))
        
        sum_mags.set_ydata((-np.log10(result)))

def app_to_abs_mag(app_mag, dist):
    """
    Converts from apparent to absolute magnitude

    Parameters
    ----------
    app_mag : Float
        Apparent magnitude of object.
    dist : float
        Distance to object.

    Returns
    -------
    Float
        Absolute magnitude.

    """
    return app_mag - 5*np.log10(dist/10)

#conversion from apparent to absolute magnitudes
arrs[1] = app_to_abs_mag(arrs[1], 68.2e6)
arrs[2] = app_to_abs_mag(arrs[2], 68.2e6)

fig = plt.figure()
ax = fig.add_subplot(111)

# Adjust the subplots region to leave some space for the sliders and buttons
#comment out if curvefit is being preformed
fig.subplots_adjust(left=0.25, bottom=0.35)

#commenting out the data which is not being used is recommended as the plot 
#gets crowded with data really quickly
ax.scatter(arrs[0], arrs[2], label='V magnitude', zorder= 3)
ax.scatter(arrs[0], arrs[1], label='B Magnitude', zorder= 3)

#---------------------Multi-Component-Fit-Parameters--------------------------#
#comment out the parameters which are not being used

#V data
I_E_bulge = 0.326e7
r_e_bulge = 2.4
n_val_bulge = 0.637
B_val_bulge = 1.997*n_val_bulge-0.327

I_E_disk = 0.1383e6
r_e_disk = 21.8
n_val_disk = 1
B_val_disk = 1.997*n_val_disk-0.327

I_E_bar = 0.537e6
r_e_bar = 7.23
n_val_bar = 0.659
B_val_bar = 1.997*n_val_bar-0.327

#B data
I_E_bulge = 0.206e7
r_e_bulge = 2.25
n_val_bulge = 0.637
B_val_bulge = 1.997*n_val_bulge-0.327

I_E_disk = 0.84e5
r_e_disk = 18.6
n_val_disk = 1
B_val_disk = 1.997*n_val_disk-0.327

I_E_bar = 0.4967e6
r_e_bar = 5.01
n_val_bar = 0.659
B_val_bar = 1.997*n_val_bar-0.327

#--------------------Single-Component-Fit-Parameters--------------------------#
#single comp fit B
I_E_bulge = 1e5
r_e_bulge = 15.36
n_val_bulge = 4.153
B_val_bulge = 1.997*n_val_bulge-0.327
#single comp fit V
I_E_bulge = 1.4e5
r_e_bulge = 19.03
n_val_bulge = 3.73
B_val_bulge = 1.997*n_val_bulge-0.327

#-------------------Calculated fits using curvefit----------------------------#
#if manual fit is attempted to find starting values comment out entire section

#for single fit to be preformed comment out multi component fit and vice versa

#multiple component fit
popt, pcov = curve_fit(func2, arrs[0][0:-1], arrs[2][0:-1], 
                       p0=(I_E_bar, r_e_bar, n_val_bar,
                           I_E_bulge, r_e_bulge, n_val_bulge, 
                           I_E_disk, r_e_disk
                           ), maxfev=int(1e9))

I_E_bar = popt[0]
r_e_bar = popt[1]
n_val_bar = popt[2]
B_val_bar = 1.997*n_val_bar-0.327

I_E_bulge = popt[3]
r_e_bulge = popt[4]
n_val_bulge = popt[5]
B_val_bulge = 1.997*n_val_bulge-0.327   

I_E_disk = popt[6]
r_e_disk = popt[7]
B_val_disk = 1.997*n_val_disk-0.327


#Single component fit
popt, pcov = curve_fit(func3, arrs[0][1:-1], arrs[1][1:-1], 
                       p0=(I_E_bulge, r_e_bulge, n_val_bulge,
                           ), maxfev=int(1e9))

I_E_bar = popt[0]
r_e_bar = popt[1]
n_val_bar = popt[2]
B_val_bar = 1.997*n_val_bar-0.327

#arrays to store each fit 1 and 3 parameter fits comment out whichever is not 
#being used
sers_fit_1_param = [func3(i, popt[0], popt[1], popt[2]) for i in arrs[0]]


sers_fit_3_param = [func2(i, popt[0], popt[1], popt[2], 
                  popt[3], popt[4], popt[5], popt[6], popt[7]) 
            for i in arrs[0]]

#----------------Calculating each component and plotting it-------------------#

#calculating the parameters
for i in xvals:
    gauss_bar.append(-2.5*np.log10(func(I_E_bar, B_val_bar, r_e_bar, 
                                        n_val_bar, i)))
    gauss_disk.append(-2.5*np.log10(func(I_E_disk, B_val_disk, r_e_disk, 
                                         n_val_disk, i)))
    gauss_bulge.append(-2.5*np.log10(func(I_E_bulge, B_val_bulge, r_e_bulge, 
                                          n_val_bulge, i)))

#for Bu in gauss_bulge:
    #sums.append(-1*np.log10(10**(-Bu)))        


for Br, D, Bu in zip(gauss_bar, gauss_disk, gauss_bulge):
    sums.append(-1*np.log10(10**(-Br) + 10**(-D) + 10**(-Bu)))
    
#plotting the parameters
gaussian_bar, = ax.plot(xvals, gauss_bar, color='g', label='bar',
                        linestyle='--')
gaussian_disk, = ax.plot(xvals, gauss_disk, color='r', label='disk',
                         linestyle='--')
gaussian_bulge, = ax.plot(xvals, gauss_bulge, color='y', label='bulge',
                          linestyle='--')

#comment out if an automatic fit with curvefit is being done
sum_mags, = ax.plot(xvals, sums, color='black', label='Sercic fit',zorder=2)

#---------------------------Root-Mean-Square-Deviation------------------------#
rms_delta = 0

#inside for loop which calculates the square sum of the deviations switch out 
#the sers_fit array to the one which is approperiate for the fit being done
for i, j in zip(sers_fit_1_param[1:-1], arrs[2][1:-1]):
    rms_delta += (i-j)**2

rms_delta = (1/(len(sums)-1))*rms_delta
rms_delta = np.sqrt(rms_delta)

#plotting error bars
ax.errorbar(arrs[0], arrs[2], yerr=rms_delta, color='red', capsize = 3, 
            fmt='.', markersize=0, zorder=1)

#comment out whichever is not being used aka if 3 param fit is done comment out
# the fuction which contains sers_fit_1_param
ax.plot(arrs[0], sers_fit_1_param ,color='black', label='Sercic fit', zorder=2)

ax.plot(arrs[0], sers_fit_3_param, color='black', label='Sercic fit', zorder=2)

print('rms error = ', rms_delta)
print('Intensity bar = ', popt[0])
print('R value bar = ', popt[1])
print('N bar = ', popt[2])
#comment out all print statements below here if single param fit is being done
print('_________________________________________')
print('Intensity Bulge = ', popt[3])
print('R value Bulge = ', popt[4])
print('N Bulge = ', popt[5])
print('_________________________________________')
print('Intensity disk = ', popt[6])
print('R value disk = ', popt[7])


#-----------------------------Sliders-----------------------------------------#
#comment out if curvefit fit is being done
#------------------------------bulge------------------------------------------#
axIe_bulge = fig.add_axes([0.4, 0, 0.45, 0.03])
axN_bulge = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
axRe_bulge = fig.add_axes([0.4, 0.05, 0.45, 0.03])

Ie_slider_bulge = Slider(
    ax=axIe_bulge,
    label='I_E_bulge',
    valmin=0,
    valmax=10000000,#change these limits if different range is necessary
    valinit=I_E_bulge,
    )

Re_slider_bulge = Slider(
    ax=axRe_bulge,
    label='R_e_bulge',
    valmin=0,
    valmax=20,#change these limits if different range is necessary
    valinit=r_e_bulge,
    )

N_slider_bulge = Slider(
    ax=axN_bulge,
    label="n_bulge",
    valmin=0,
    valmax=5,#change these limits if different range is necessary
    valinit=n_val_bulge,
    orientation="vertical"
)

#checking for change in sliders and updating the line if any
Ie_slider_bulge.on_changed(update_bulge)
Re_slider_bulge.on_changed(update_bulge)
N_slider_bulge.on_changed(update_bulge)
#------------------------------disk-------------------------------------------#
axIe_disk = fig.add_axes([0.4, 0.1, 0.45, 0.03])
axN_disk = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
axRe_disk = fig.add_axes([0.4, 0.15, 0.45, 0.03])

Ie_slider_disk = Slider(
    ax=axIe_disk,
    label='I_E_disk',
    valmin=2000,
    valmax=1000000,#change these limits if different range is necessary
    valinit=I_E_disk,
    )

Re_slider_disk = Slider(
    ax=axRe_disk,
    label='R_e_disk',
    valmin=0,
    valmax=100,#change these limits if different range is necessary
    valinit=r_e_disk,
    )

N_slider_disk = Slider(
    ax=axN_disk,
    label="n_disk",
    valmin=0,
    valmax=5,#change these limits if different range is necessary
    valinit=n_val_disk,
    orientation="vertical"
)

#checking for change in sliders and updating the line if any
Ie_slider_disk.on_changed(update_disk)
Re_slider_disk.on_changed(update_disk)
N_slider_disk.on_changed(update_disk)


#------------------------------Bar--------------------------------------------#
axIe_bar = fig.add_axes([0.4, 0.2, 0.45, 0.03])
axRe_bar = fig.add_axes([0.4, 0.25, 0.45, 0.03])
axN_bar = fig.add_axes([0.15, 0.25, 0.0225, 0.63])

Ie_slider_bar = Slider(
    ax=axIe_bar,
    label='I_E_bar',
    valmin=5000,
    valmax=1000000,#change these limits if different range is necessary
    valinit=I_E_bar,
    ) 

Re_slider_bar = Slider(
    ax=axRe_bar,
    label='R_e_bar',
    valmin=0,
    valmax=10,#change these limits if different range is necessary
    valinit=r_e_bar,
    )

N_slider_bar = Slider(
    ax=axN_bar,
    label="n_bar",
    valmin=0,
    valmax=2,#change these limits if different range is necessary
    valinit=n_val_bar,
    orientation="vertical"
)
#checking for change in sliders and updating the line if any
Ie_slider_bar.on_changed(update_bar)
Re_slider_bar.on_changed(update_bar)
N_slider_bar.on_changed(update_bar)

ax.set_ylim(-18, -12)
ax.invert_yaxis()
ax.set_ylabel('Surface Brightness (mag/arcsec²)')
ax.set_xlabel('Semi Major axis (arcsec)')
ax.grid()
ax.legend()