#----------------------------------------------------------------------#
# filename: pyfits+psycopg2_ex.py 
# author: Peter Nugent
# date: 10/30/2014
# ---------------------------------------------------------------------#
# Function: Read in a Arjun fits binary table from standard in and load it
# into the desi calib pg table database with psycopg2.
# ---------------------------------------------------------------------#

# First, we'll load up the dependencies:

# psycopg2 is an open-source postgres client for python. 
# We may want db access at some point and, of course, pyfits & sys

import psycopg2 

import astropy
from astropy.io import fits

import sys, os, re, glob, distutils
from distutils.util import strtobool
import numpy as np


# Read in the fits binary table 

fitsbin = str(sys.argv[1])

# First, open the table using pyfits:

table = fits.open( fitsbin )


# Access the data part of the table.

tbdata = table[1].data
tbdata = np.asarray(tbdata)

# determine the number of elements 

nrows = tbdata.shape[0] 
newdata = []

# Fire up the db

con = psycopg2.connect(host='scidb2.nersc.gov', user='desi_admin', password='L00cy-1959', database='desi')
cursor = con.cursor()

# Re-cast as strings so the load is easy 

for i in range(0, nrows):

    if tbdata['has_g'][i] == 70:
        g = 0
    else:
	g = 1

    if tbdata['has_r'][i] == 70:
        r = 0
    else:
	r = 1

    if tbdata['has_z'][i] == 70:
        z = 0
    else:
	z = 1

    line = [ tbdata['brickname'][i], tbdata['brickid'][i], tbdata['brickrow'][i], tbdata['brickcol'][i], tbdata['brickq'][i], tbdata['ra'][i], tbdata['dec'][i], tbdata['ra1'][i], tbdata['dec1'][i], tbdata['ra2'][i], tbdata['dec2'][i], bool(g), bool(r), bool(z) ]

    newdata.append(line) 

## Re-cast as strings so the load is easy 
#
for i, f in enumerate(newdata):
##
   query = 'INSERT INTO bricks ( brickname, brickid, brickrow, brickcol, brickq, ra, dec, ra1, dec1, ra2, dec2, has_g, has_r, has_z) VALUES ( %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s )' 
##
   cursor.execute( query, tuple( [str(elem) for elem in newdata[i]] ) ) 
#
#

con.commit()

# That's it!

