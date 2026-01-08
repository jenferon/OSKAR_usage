import matplotlib.pyplot as plt
import numpy as np
from astropy import units as un
import py21cmsense as p21sense
from py21cmsense import GaussianBeam, Observation, Observatory, PowerSpectrum, beam
from astropy.cosmology import Planck13 as cosmoP
from astropy.cosmology import FlatLambdaCDM, LambdaCDM
from py21cmsense.theory import EOS2021
import astropy.units as u
import math
cosmo = FlatLambdaCDM(H0=71 * u.km / u.s / u.Mpc, Om0=0.27)
import astropy.io.fits as fits

base = "/gpfs01/home/ppxjf3/OSKAR/"

def zFromNu(nu):
    """
    Convert frequency of 21cm line to redshift
    
    Input: nu [MHz]
    """
    nu21 = 1.420405e3  #MHz
    return nu21/nu - 1.0

def NuFromz(z):
    """
    Convert frequency of 21cm line to redshift
    
    Input: nu [MHz]
    """
    nu21 = 1.420405e3  #MHz
    return nu21/(z + 1.0)


nu = 155.0
c = 2.998e8
eta_a = 1.0 #antenna efficiency
pi = 3.14159
k_b = 1.38e-23
eta_s=0.9 #system efficiency
T_sys = 40.0 + 1.1*60.0*(nu*1e6/c)**(-2.55) # from Table 3 https://astronomers.skatelescope.org/wp-content/uploads/2016/12/SKA-TEL-SKO-DD-001-1_BaselineDesign1.pdf

zs = zFromNu(nu)

baseline_data = np.loadtxt(base+"SKA_all_stations_Rev3.tm/layout_ecef.txt")
aa4 =  Observatory(antpos=baseline_data* u.m, beam=beam.GaussianBeam(
           frequency=nu*u.MHz, dish_size=8*u.m), latitude = -27.0*u.deg, Trcv=T_sys*u.K) #check beam and dish size

"""layout = obs.antpos.value[:, :2]  # exclude z position
print("obs has", layout.shape[0], "antennas")
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.scatter(layout[:, 0], layout[:, 1])
plt.legend()
plt.ylabel("Y position [m]")
plt.xlabel("X position [m]")
plt.show()
plt.savefig('/home/ppxjf3/RSD_LC/layout.pdf')
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
im = ax.hexbin(
    obs.baselines_metres[:, :, 0].ravel(),
    obs.baselines_metres[:, :, 1].ravel(),
    gridsize=20,
    bins="log",
    vmin=100,
)
plt.colorbar(im, ax=ax, label="Number of baselines")
ax.set_xlabel("X baseline [m]")
ax.set_ylabel("Y baseline [m]")
plt.savefig('/home/ppxjf3/RSD_LC/baselines.pdf')
plt.close()
"""
#observation
observation_params = {}
observation_params["ndays"] = 166.7
observation_params["cosmo"] = cosmo
observation_params["h"] = cosmo.H0.value / 100.0
observation_params["freq_bands"] = nu
observation_params["redshifts"] = zs
observation_params["time_per_day_hrs"] = 6.0
observation_params["bandwidth"] = 0.1e6

observation = p21sense.Observation(
    observatory=aa4,
    lst_bin_size=observation_params["time_per_day_hrs"] * un.hour,
    time_per_day=observation_params["time_per_day_hrs"] * un.hour,
    n_days=observation_params["ndays"],
    bandwidth=observation_params["bandwidth"]*un.Hz,
    coherent=True, #add baselines coherantly or not
    cosmo=observation_params["cosmo"],
    baseline_filters=p21sense.BaselineRange(bl_max=3400* un.m),
#    n_channels=1,
)




ska_aa4_senses1 =PowerSpectrum(horizon_buffer=0.0/u.Mpc,foreground_model="foreground_free", 
            observation=observation,
        ).at_frequency(nu * un.MHz)
sense1d_both = ska_aa4_senses1.calculate_sensitivity_1d(thermal=True, sample=False)
print(ska_aa4_senses1.k1d)
try:
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    plt.plot(
        ska_aa4_senses1.k1d,
        sense1d_both
        
    )
    plt.xscale("log")
    plt.yscale("log")
    ax.set_ylabel("sensitivity [mK^2]")
    ax.set_xlabel("k [Mpc^-1]")
    plt.savefig(base+"test_senses.pdf")
    plt.close()
except:
    print('could not plot')

import tools21cm as t2c

def get_lengths(nu_low, nu_hi, nu_mid, theta_FOV):
    z_lo = zFromNu(nu_low)
    z_mid = zFromNu(nu_mid)
    z_hi = zFromNu(nu_hi)
    
    L_para = cosmo.comoving_distance(z_lo) - cosmo.comoving_distance(z_hi)
    L_perp = cosmo.comoving_distance(z_mid) * (np.pi * theta_FOV / 180.0)
    return L_para/u.Mpc, L_perp/u.Mpc
"""
with fits.open('Noise_N512_FOV1.000_155.0MHz_SKA_SKA_central_area_Rev3_EOR0_1.0_0512_natural_I_I.fits', memmap=True) as hdu:
    NOISE =np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()
    
with fits.open('Noise_N512_FOV1.000_155.0MHz_SKA_SKA_central_area_Rev3_EOR0_1.0_0512_natural_I_I_isotropic_beam.fits', memmap=True) as hdu:
    NOISE_ISO =np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()

with fits.open('Noise_N512_FOV1.000_155.0MHz_SKA_SKA_central_area_Rev3_EOR0_1.0_0512_natural_I_I_Aperture_beam.fits', memmap=True) as hdu:
    NOISE_AP =np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()
    
    
NOISE = NOISE *1e3
NOISE_ISO = NOISE_ISO*1e3
NOISE_AP = NOISE_AP *1e3
NOISE = np.append(NOISE,NOISE)
NOISE_ISO = np.append(NOISE_ISO,NOISE_ISO)
NOISE_AP = np.append(NOISE_AP,NOISE_AP)
L_para, L_perp = get_lengths(150, 160, 155, 1.0)
box_dims = [L_perp, L_perp, L_para]
V_NOISE = L_perp*L_perp*L_para

p_NOISE, k_NOISE, n_NOISE = t2c.power_spectrum_1d(NOISE, box_dims=box_dims, kbins=15, binning = 'log', return_n_modes=True)
d_NOISE = (p_NOISE*V_NOISE*k_NOISE**3)/(2*np.pi**2)

p_NOISE_ISO, k_NOISE_ISO, n_NOISE_ISO = t2c.power_spectrum_1d(NOISE_ISO, box_dims=box_dims, kbins=15, binning = 'log', return_n_modes=True)
d_NOISE_ISO = (p_NOISE_ISO*V_NOISE*k_NOISE_ISO**3)/(2*np.pi**2)

p_NOISE_AP, k_NOISE_AP, n_NOISE_AP = t2c.power_spectrum_1d(NOISE_AP, box_dims=box_dims, kbins=15, binning = 'log', return_n_modes=True)
d_NOISE_AP = (p_NOISE_AP*V_NOISE*k_NOISE_AP**3)/(2*np.pi**2)

plt.plot(k_NOISE, d_NOISE, label='gauss beam')
plt.plot(k_NOISE_ISO, d_NOISE_ISO, label='iso beam')
plt.plot(k_NOISE_AP, d_NOISE_AP, label='aperture beam')
plt.yscale('log')
plt.xlabel('k')
plt.xscale('log')
plt.ylabel(r'$\Delta(k)$')
plt.legend()
plt.savefig('oskar_noise.pdf')
plt.close()

"""
with fits.open('/gpfs01/home/ppxjf3/peculiar_vel/data/NOISE/Noise_N512_FOV1.000_155.0MHz_SKA_SKA_all_stations_Rev3_EOR0_1.0_0512_natural_I_I.fits', memmap=True) as hdu:
    NOISE_all =np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()
NOISE_all=NOISE_all*3
#NOISE_all= np.append(NOISE_all,NOISE_all)
print(NOISE_all.shape)
p_NOISE_all, k_NOISE_all, n_NOISE_all = t2c.power_spectrum_1d(NOISE_all, box_dims=box_dims, kbins=15, binning = 'log', return_n_modes=True)
d_NOISE_all = (p_NOISE_all*V_NOISE*k_NOISE_all**3)/(2*np.pi**2)
"""with fits.open('Noise_N512_FOV1.000_155.0MHz_SKA_SKA_core_area_Rev3_EOR0_1.0_0512_natural_I_I.fits', memmap=True) as hdu:
    NOISE_core =np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()
NOISE_core=NOISE_core*3
NOISE_core= np.append(NOISE_core,NOISE_core)

p_NOISE_core, k_NOISE_core, n_NOISE_core = t2c.power_spectrum_1d(NOISE_core, box_dims=box_dims, kbins=15, binning = 'log', return_n_modes=True)
d_NOISE_core = (p_NOISE_core*V_NOISE*k_NOISE_core**3)/(2*np.pi**2)
plt.plot(k_NOISE, d_NOISE, label='core')
plt.plot(k_NOISE_all, d_NOISE_all, label='All')
plt.plot(k_NOISE_AP, d_NOISE_AP, label='central')
plt.yscale('log')
plt.xlabel('k')
plt.xscale('log')
plt.ylabel(r'$\Delta(k)$')
plt.legend()
plt.savefig('stations_noise.pdf')
plt.close()
"""
t_int_h = 1000 # number of hours to normalise noise to
freq = 155.0 # first frequency map / MHz

n_tiles = 1 # per station
n_dipoles = 256 # per tile
Di = 8 # from baseline design.
const_PSF = 0

n_stations = 512
D =  3400



t_int = t_int_h*3600. # integration time / ss
FWHM_arcmin = 4.0

W = 1.3 #weighting factor (1.3-2) #WE SET THIS AS 1 FOR SUSSEX but also didnt have factor of 2. with equation as they are now, set to 2 to get 67 mK at 150 MHz for 48 stattions. 1.3 89 mK for 24 stations . set to 1.3 to match website.
lamb = c/(freq*1e6)


# Calculate effective area at each frequency
eta_rad = (0.056 * freq +82.2)/100.0 # worked out since 85 at 50 MHz and 99 at 300. asssumed linear.
A_di=lamb**2/(4*pi)*eta_rad*Di #from Table 3 https://astronomers.skatelescope.org/wp-content/uploads/2016/12/SKA-TEL-SKO-DD-001-1_BaselineDesign1.pdf foot note
if (A_di > 3.2):
    A_di=3.2 # Lower frequencies limited by antenna area
A=A_di*n_dipoles*n_tiles # Area per station.
K=(A)/(2.0*k_b)
SEFD=T_sys/K # in J m^-2=W m^-2 Hz^-1
SEFD=SEFD*1.0e26 # in Jy
 # Calculate noise sensitivities in Jy 
noise_Jy=(W/eta_s)*(SEFD/math.sqrt(2*n_stations*(n_stations-1)*0.1e6*t_int))
FWHM = 1.22 * lamb / D  # radians

beamarea = pi * FWHM**2 / (4.0*math.log(2))
rms =0.
print(NOISE_all.shape)
for ii in range(0,512):
    for jj in range(0,512):
        rms = rms + NOISE_all[0,ii,jj]**2
rms = math.sqrt(rms/(512*512)) 
NOISE_all = NOISE_all * (noise_Jy/rms) # the slice is now in Janskys
NOISE_all = NOISE_all *1e-26*(c/(freq*1.0e6))**2*(1.0/(2.0*k_b*beamarea)) # converting from Jy to K
  
  
p_NOISE_all, k_NOISE_all, n_NOISE_all = t2c.power_spectrum_1d(NOISE_all, box_dims=box_dims, kbins=15, binning = 'log', return_n_modes=True)
d_NOISE_all = (p_NOISE_all*V_NOISE*k_NOISE_all**3)/(2*np.pi**2)
  
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
plt.plot(
        ska_aa4_senses1.k1d,
        sense1d_both, label='21cmSENSE'
        
)
plt.plot(k_NOISE_all, d_NOISE_all, label='OSKAR')
plt.xscale("log")
plt.yscale("log")
ax.set_ylabel("sensitivity [mK^2]")
ax.set_xlabel("k [Mpc^-1]")
plt.legend()
plt.savefig(base+"tried_to_normalise.pdf")
plt.close()
