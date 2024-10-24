import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord,FK5,get_sun
from astropy.time import Time
import argparse
import os
import psrqpy
import ephem


parser = argparse.ArgumentParser(description='N_earth vs Time plot')


parser.add_argument('-psrname', type=str, required=False,help='Name of the pulsar')
parser.add_argument('-raj', type=str, required=False,help='RA of the source')
parser.add_argument('-decj', type=str, required=False,help='DEC of the source')
parser.add_argument('-parfile', type=str, required=False,help='par file')
parser.add_argument('-dm_sw_file', type=str, required=True,help='DM SW file. This file should have three columns MJD, DM and DM error. The columns should be separated by space and the comments should start with #')

args = parser.parse_args()

psrname = args.psrname
parfile = args.parfile
raj = None
decj = None

db = psrqpy.QueryATNF()


def solar_angle(ra,dec,mjd):
    coo_psr = SkyCoord(ra,dec, unit = 'hourangle,degree',frame=FK5)
    t = Time(mjd, format='mjd')
    sunpos = get_sun(t)
    sep = sunpos.separation(coo_psr)
    return sep.value#.degree[0]

sun = ephem.Sun()

def calculate_distance_from_sun(mjd):
    # Specify the Modified Julian Date (MJD)
    date = Time(mjd, format='mjd').iso
    sun.compute(date)
    distance = sun.earth_distance

    return distance


if args.parfile:
    raj = [line.split()[1] for line in open(parfile) if line.startswith("RAJ ")][0]
    decj = [line.split()[1] for line in open(parfile) if line.startswith("DECJ ")][0]
    if raj==None or decj==None:
        print("RAJ and DECJ not found in the par file. So trying to retrieve from the ATNF catalogue")
        raj = db[psrname]['RAJ'].value[0]
        decj = db[psrname]['DECJ'].value[0]
        
else:
    if not args.raj or not args.decj:
        print("RAJ and DECJ not provided. So trying to retrieve from the ATNF catalogue")
        raj = db[psrname]['RAJ'].value[0]
        decj = db[psrname]['DECJ'].value[0]
    else:
        raj = args.raj
        decj = args.decj

if raj==None or decj==None:
    print("RAJ and DECJ not found in the par file or it is not provided from the pulsar name as well")
    exit()  

print("RAJ: ", raj)
print("DECJ: ", decj) 

distance_from_sun=[]


# Read the DM SW file
dm_sw_file = args.dm_sw_file
dm_sw_data = np.loadtxt(dm_sw_file, dtype='str', delimiter=' ',comments='#')

mjd_rec = dm_sw_data[:,0].astype(float)
dm_rec = dm_sw_data[:,1].astype(float)
dmerr_rec = dm_sw_data[:,2].astype(float)

# Calculate the distance from the Sun
for i in range(len(mjd_rec)):
    distance_from_sun.append(calculate_distance_from_sun(mjd_rec[i]))


# Convert the distance to Astronomical Units (AU)
distance_au = distance_from_sun * u.au
sangle = solar_angle(raj,decj,mjd_rec)*np.pi/180.0
frac = (np.pi-sangle)/np.sin(np.pi-sangle)
recovered_nesw=(dm_rec)*distance_au/(4.84814e-6*u.au*frac)
recovered_nesw_err = dmerr_rec*distance_au/(4.84814e-6*u.au*frac)


mask = np.where(sangle<np.pi/4.0)[0]

plt.errorbar(mjd_rec,recovered_nesw,recovered_nesw_err,fmt = '.',color='k',label='n$_e$ values greater than 45$^\circ$')
plt.errorbar(mjd_rec[mask],recovered_nesw[mask],recovered_nesw_err[mask],fmt = '.',color = 'r', label='n$_e$ values less than 45$^\circ$')
plt.axhline(y=np.median(recovered_nesw),color='k',linestyle='--',label='Median n$_e$')
plt.xlabel('MJD')
plt.ylabel('n$_e$(cm$^{-3}$)')
plt.title('n$_e$ vs Time')
plt.legend()
plt.show()
