#----------------------------------------------------------------------#
# filename: astropy+psycopg2_ex.py 
# author: Peter Nugent
# date: 2/18/2015
# ---------------------------------------------------------------------#
# Function: Read in a Dustin's fits tractor binary table from standard 
# in and load it into the desi candidate pg table database with psycopg2.
# ---------------------------------------------------------------------#

# First, we'll load up the dependencies:

# psycopg2 is an open-source postgres client for python. 
# We may want db access at some point and, of course, fits & sys

import psycopg2 

import astropy
from astropy.io import fits

import sys, os, re, glob
import numpy as np

# Read in the image name (ooi) and create the mask (ood) and weight (oow) names

fitsbin = str(sys.argv[1])

fimage = os.path.abspath(sys.argv[1])

# First, open the table using pyfits:

table = fits.open( fitsbin )

#
# Fire up the db
#
con = psycopg2.connect(host='scidb2.nersc.gov', user='desi_admin', password='L00cy-1959', database='desi')

cursor = con.cursor()

prihdr=table[0].header

tractver = prihdr['TRACTORV']
hdrbrickid = prihdr['brickid']

##
cursor.execute( "SELECT loaded from bricks where brickid = %s", (hdrbrickid,)  )
loaded = cursor.fetchone()[0]

if loaded:
  print fimage, 'is already loaded'
  sys.exit(0)

##
cursor.execute( "UPDATE bricks set filename = %s, loaded = 'true', tractorvr = %s where brickid = %s", (fimage, tractver, hdrbrickid,)  )

# Access the data part of the table.

tbdata = table[1].data
tbdata = np.asarray(tbdata)

# determine the number of elements 

nrows = tbdata.shape[0] 

newdata_cand = []
newdata_decam = []
newdata_wise = []
newdata_decap = []



for i in range(0, nrows):

   if  tbdata['brick_primary'][i] == 84:

      if tbdata['left_blob'][i] == 70:
         lb = 0
      else:
         lb = 1

      if tbdata['out_of_bounds'][i] == 70:
         oob = 0
      else:
         oob = 1
 
   
      line_cand = [ tbdata['brickid'][i], tbdata['objid'][i], tbdata['blob'][i], tbdata['type'][i], tbdata['ra'][i], tbdata['ra_ivar'][i], tbdata['dec'][i], tbdata['dec_ivar'][i], tbdata['bx'][i], tbdata['by'][i], tbdata['bx0'][i], tbdata['by0'][i], bool(lb), bool(oob), tbdata['ebv'][i], tbdata['dchisq'][i][0], tbdata['dchisq'][i][1], tbdata['dchisq'][i][2], tbdata['dchisq'][i][3], tbdata['fracDev'][i], tbdata['fracDev_ivar'][i], tbdata['shapeExp_r'][i], tbdata['shapeExp_r_ivar'][i], tbdata['shapeExp_e1'][i], tbdata['shapeExp_e1_ivar'][i], tbdata['shapeExp_e2'][i], tbdata['shapeExp_e2_ivar'][i], tbdata['shapeDev_r'][i], tbdata['shapeDev_r_ivar'][i], tbdata['shapeDev_e1'][i], tbdata['shapeDev_e1_ivar'][i], tbdata['shapeDev_e2'][i], tbdata['shapeDev_e2_ivar'][i] ]

      line_decam = [ tbdata['decam_flux'][i][0], tbdata['decam_flux_ivar'][i][0], tbdata['decam_fracflux'][i][0], tbdata['decam_fracmasked'][i][0], tbdata['decam_fracin'][i][0], tbdata['decam_rchi2'][i][0], tbdata['decam_nobs'][i][0], tbdata['decam_anymask'][i][0], tbdata['decam_allmask'][i][0], tbdata['decam_mw_transmission'][i][0],  tbdata['decam_flux'][i][1], tbdata['decam_flux_ivar'][i][1], tbdata['decam_fracflux'][i][1], tbdata['decam_fracmasked'][i][1], tbdata['decam_fracin'][i][1], tbdata['decam_rchi2'][i][1], tbdata['decam_nobs'][i][1], tbdata['decam_anymask'][i][1], tbdata['decam_allmask'][i][1], tbdata['decam_mw_transmission'][i][1],  tbdata['decam_flux'][i][2], tbdata['decam_flux_ivar'][i][2], tbdata['decam_fracflux'][i][2], tbdata['decam_fracmasked'][i][2], tbdata['decam_fracin'][i][2], tbdata['decam_rchi2'][i][2], tbdata['decam_nobs'][i][2], tbdata['decam_anymask'][i][2], tbdata['decam_allmask'][i][2], tbdata['decam_mw_transmission'][i][2],  tbdata['decam_flux'][i][3], tbdata['decam_flux_ivar'][i][3], tbdata['decam_fracflux'][i][3], tbdata['decam_fracmasked'][i][3], tbdata['decam_fracin'][i][3], tbdata['decam_rchi2'][i][3], tbdata['decam_nobs'][i][3], tbdata['decam_anymask'][i][3], tbdata['decam_allmask'][i][3], tbdata['decam_mw_transmission'][i][3],  tbdata['decam_flux'][i][4], tbdata['decam_flux_ivar'][i][4], tbdata['decam_fracflux'][i][4], tbdata['decam_fracmasked'][i][4], tbdata['decam_fracin'][i][4], tbdata['decam_rchi2'][i][4], tbdata['decam_nobs'][i][4], tbdata['decam_anymask'][i][4], tbdata['decam_allmask'][i][4], tbdata['decam_mw_transmission'][i][4],  tbdata['decam_flux'][i][5], tbdata['decam_flux_ivar'][i][5], tbdata['decam_fracflux'][i][5], tbdata['decam_fracmasked'][i][5], tbdata['decam_fracin'][i][5], tbdata['decam_rchi2'][i][5], tbdata['decam_nobs'][i][5], tbdata['decam_anymask'][i][5], tbdata['decam_allmask'][i][5], tbdata['decam_mw_transmission'][i][5] ]

      line_wise = [ tbdata['wise_flux'][i][0], tbdata['wise_flux_ivar'][i][0], tbdata['wise_fracflux'][i][0], tbdata['wise_rchi2'][i][0], tbdata['wise_nobs'][i][0], tbdata['wise_mw_transmission'][i][0], tbdata['wise_flux'][i][1], tbdata['wise_flux_ivar'][i][1], tbdata['wise_fracflux'][i][1], tbdata['wise_rchi2'][i][1], tbdata['wise_nobs'][i][1], tbdata['wise_mw_transmission'][i][1], tbdata['wise_flux'][i][2], tbdata['wise_flux_ivar'][i][2], tbdata['wise_fracflux'][i][2], tbdata['wise_rchi2'][i][2], tbdata['wise_nobs'][i][2], tbdata['wise_mw_transmission'][i][2], tbdata['wise_flux'][i][3], tbdata['wise_flux_ivar'][i][3], tbdata['wise_fracflux'][i][3], tbdata['wise_rchi2'][i][3], tbdata['wise_nobs'][i][3], tbdata['wise_mw_transmission'][i][3] ]

      line_decap = [ tbdata['decam_apflux'][i][0], tbdata['decam_apflux'][i][1], tbdata['decam_apflux'][i][2], tbdata['decam_apflux'][i][3], tbdata['decam_apflux'][i][4], tbdata['decam_apflux'][i][5], tbdata['decam_apflux'][i][6], tbdata['decam_apflux'][i][7], tbdata['decam_apflux_resid'][i][0], tbdata['decam_apflux_resid'][i][1], tbdata['decam_apflux_resid'][i][2], tbdata['decam_apflux_resid'][i][3], tbdata['decam_apflux_resid'][i][4], tbdata['decam_apflux_resid'][i][5], tbdata['decam_apflux_resid'][i][6], tbdata['decam_apflux_resid'][i][7], tbdata['decam_apflux_ivar'][i][0], tbdata['decam_apflux_ivar'][i][1], tbdata['decam_apflux_ivar'][i][2], tbdata['decam_apflux_ivar'][i][3], tbdata['decam_apflux_ivar'][i][4], tbdata['decam_apflux_ivar'][i][5], tbdata['decam_apflux_ivar'][i][6], tbdata['decam_apflux_ivar'][i][7], tbdata['decam_apflux'][i][8], tbdata['decam_apflux'][i][9], tbdata['decam_apflux'][i][10], tbdata['decam_apflux'][i][11], tbdata['decam_apflux'][i][12], tbdata['decam_apflux'][i][13], tbdata['decam_apflux'][i][14], tbdata['decam_apflux'][i][15], tbdata['decam_apflux_resid'][i][8], tbdata['decam_apflux_resid'][i][9], tbdata['decam_apflux_resid'][i][10], tbdata['decam_apflux_resid'][i][11], tbdata['decam_apflux_resid'][i][12], tbdata['decam_apflux_resid'][i][13], tbdata['decam_apflux_resid'][i][14], tbdata['decam_apflux_resid'][i][15], tbdata['decam_apflux_ivar'][i][8], tbdata['decam_apflux_ivar'][i][9], tbdata['decam_apflux_ivar'][i][10], tbdata['decam_apflux_ivar'][i][11], tbdata['decam_apflux_ivar'][i][12], tbdata['decam_apflux_ivar'][i][13], tbdata['decam_apflux_ivar'][i][14], tbdata['decam_apflux_ivar'][i][15], tbdata['decam_apflux'][i][16], tbdata['decam_apflux'][i][17], tbdata['decam_apflux'][i][18], tbdata['decam_apflux'][i][19], tbdata['decam_apflux'][i][20], tbdata['decam_apflux'][i][21], tbdata['decam_apflux'][i][22], tbdata['decam_apflux'][i][23], tbdata['decam_apflux_resid'][i][16], tbdata['decam_apflux_resid'][i][17], tbdata['decam_apflux_resid'][i][18], tbdata['decam_apflux_resid'][i][19], tbdata['decam_apflux_resid'][i][20], tbdata['decam_apflux_resid'][i][21], tbdata['decam_apflux_resid'][i][22], tbdata['decam_apflux_resid'][i][23], tbdata['decam_apflux_ivar'][i][16], tbdata['decam_apflux_ivar'][i][17], tbdata['decam_apflux_ivar'][i][18], tbdata['decam_apflux_ivar'][i][19], tbdata['decam_apflux_ivar'][i][20], tbdata['decam_apflux_ivar'][i][21], tbdata['decam_apflux_ivar'][i][22], tbdata['decam_apflux_ivar'][i][23], tbdata['decam_apflux'][i][24], tbdata['decam_apflux'][i][25], tbdata['decam_apflux'][i][26], tbdata['decam_apflux'][i][27], tbdata['decam_apflux'][i][28], tbdata['decam_apflux'][i][29], tbdata['decam_apflux'][i][30], tbdata['decam_apflux'][i][31], tbdata['decam_apflux_resid'][i][24], tbdata['decam_apflux_resid'][i][25], tbdata['decam_apflux_resid'][i][26], tbdata['decam_apflux_resid'][i][27], tbdata['decam_apflux_resid'][i][28], tbdata['decam_apflux_resid'][i][29], tbdata['decam_apflux_resid'][i][30], tbdata['decam_apflux_resid'][i][31], tbdata['decam_apflux_ivar'][i][24], tbdata['decam_apflux_ivar'][i][25], tbdata['decam_apflux_ivar'][i][26], tbdata['decam_apflux_ivar'][i][27], tbdata['decam_apflux_ivar'][i][28], tbdata['decam_apflux_ivar'][i][29], tbdata['decam_apflux_ivar'][i][30], tbdata['decam_apflux_ivar'][i][31], tbdata['decam_apflux'][i][32], tbdata['decam_apflux'][i][33], tbdata['decam_apflux'][i][34], tbdata['decam_apflux'][i][35], tbdata['decam_apflux'][i][36], tbdata['decam_apflux'][i][37], tbdata['decam_apflux'][i][38], tbdata['decam_apflux'][i][39], tbdata['decam_apflux_resid'][i][32], tbdata['decam_apflux_resid'][i][33], tbdata['decam_apflux_resid'][i][34], tbdata['decam_apflux_resid'][i][35], tbdata['decam_apflux_resid'][i][36], tbdata['decam_apflux_resid'][i][37], tbdata['decam_apflux_resid'][i][38], tbdata['decam_apflux_resid'][i][39], tbdata['decam_apflux_ivar'][i][32], tbdata['decam_apflux_ivar'][i][33], tbdata['decam_apflux_ivar'][i][34], tbdata['decam_apflux_ivar'][i][35], tbdata['decam_apflux_ivar'][i][36], tbdata['decam_apflux_ivar'][i][37], tbdata['decam_apflux_ivar'][i][38], tbdata['decam_apflux_ivar'][i][39], tbdata['decam_apflux'][i][40], tbdata['decam_apflux'][i][41], tbdata['decam_apflux'][i][42], tbdata['decam_apflux'][i][43], tbdata['decam_apflux'][i][44], tbdata['decam_apflux'][i][45], tbdata['decam_apflux'][i][46], tbdata['decam_apflux'][i][47], tbdata['decam_apflux_resid'][i][40], tbdata['decam_apflux_resid'][i][41], tbdata['decam_apflux_resid'][i][42], tbdata['decam_apflux_resid'][i][43], tbdata['decam_apflux_resid'][i][44], tbdata['decam_apflux_resid'][i][45], tbdata['decam_apflux_resid'][i][46], tbdata['decam_apflux_resid'][i][47], tbdata['decam_apflux_ivar'][i][40], tbdata['decam_apflux_ivar'][i][41], tbdata['decam_apflux_ivar'][i][42], tbdata['decam_apflux_ivar'][i][43], tbdata['decam_apflux_ivar'][i][44], tbdata['decam_apflux_ivar'][i][45], tbdata['decam_apflux_ivar'][i][46], tbdata['decam_apflux_ivar'][i][47] ]

      newdata_cand.append(line_cand)
      newdata_decam.append(line_decam)
      newdata_wise.append(line_wise)
      newdata_decap.append(line_decap)

#
#
## Re-cast as strings so the load is easy 
#
for i, f in enumerate(newdata_cand):
##
   query = 'INSERT INTO candidate ( brickid, objid, blob, type, ra, ra_ivar, dec, dec_ivar, bx, by, bx0, by0, left_blob, out_of_bounds, ebv, dchisq1, dchisq2, dchisq3, dchisq4, fracdev, fracdev_ivar, shapeexp_r, shapeexp_r_ivar, shapeexp_e1, shapeexp_e1_ivar, shapeexp_e2, shapeexp_e2_ivar, shapedev_r, shapedev_r_ivar, shapedev_e1, shapedev_e1_ivar, shapedev_e2, shapedev_e2_ivar ) VALUES ( %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s ) RETURNING id' 
##
   cursor.execute( query, tuple( [str(elem) for elem in newdata_cand[i]] ) ) 
   id = cursor.fetchone()[0]
##
   query = 'INSERT INTO decam ( cand_id, uflux, uflux_ivar, ufracflux, ufracmasked, ufracin, u_rchi2, unobs, u_anymask, u_allmask, u_ext, gflux, gflux_ivar, gfracflux, gfracmasked, gfracin, g_rchi2, gnobs, g_anymask, g_allmask, g_ext, rflux, rflux_ivar, rfracflux, rfracmasked, rfracin, r_rchi2, rnobs, r_anymask, r_allmask, r_ext, iflux, iflux_ivar, ifracflux, ifracmasked, ifracin, i_rchi2, inobs, i_anymask, i_allmask, i_ext, zflux, zflux_ivar, zfracflux, zfracmasked, zfracin, z_rchi2, znobs, z_anymask, z_allmask, z_ext, yflux, yflux_ivar, yfracflux, yfracmasked, yfracin, y_rchi2, ynobs, y_anymask, y_allmask, y_ext) VALUES ( %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s )' 
##
   cursor.execute( query, tuple( [str(elem) for elem in [ id ] + newdata_decam[i]] ) )
##
##
   query = 'INSERT INTO wise (cand_id, w1flux, w1flux_ivar, w1fracflux, w1_rchi2, w1nobs, w1_ext, w2flux, w2flux_ivar, w2fracflux, w2_rchi2, w2nobs, w2_ext, w3flux, w3flux_ivar, w3fracflux, w3_rchi2, w3nobs, w3_ext, w4flux, w4flux_ivar, w4fracflux, w4_rchi2, w4nobs, w4_ext ) VALUES ( %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s )' 
##
   cursor.execute( query, tuple( [str(elem) for elem in [ id ] + newdata_wise[i]] ) )
#
##
   query = 'INSERT INTO decam_aper ( cand_id, uapflux_1, uapflux_2, uapflux_3, uapflux_4, uapflux_5, uapflux_6, uapflux_7, uapflux_8, uapflux_resid_1, uapflux_resid_2, uapflux_resid_3, uapflux_resid_4, uapflux_resid_5, uapflux_resid_6, uapflux_resid_7, uapflux_resid_8, uapflux_ivar_1, uapflux_ivar_2, uapflux_ivar_3, uapflux_ivar_4, uapflux_ivar_5, uapflux_ivar_6, uapflux_ivar_7, uapflux_ivar_8, gapflux_1, gapflux_2, gapflux_3, gapflux_4, gapflux_5, gapflux_6, gapflux_7, gapflux_8, gapflux_resid_1, gapflux_resid_2, gapflux_resid_3, gapflux_resid_4, gapflux_resid_5, gapflux_resid_6, gapflux_resid_7, gapflux_resid_8, gapflux_ivar_1, gapflux_ivar_2, gapflux_ivar_3, gapflux_ivar_4, gapflux_ivar_5, gapflux_ivar_6, gapflux_ivar_7, gapflux_ivar_8, rapflux_1, rapflux_2, rapflux_3, rapflux_4, rapflux_5, rapflux_6, rapflux_7, rapflux_8, rapflux_resid_1, rapflux_resid_2, rapflux_resid_3, rapflux_resid_4, rapflux_resid_5, rapflux_resid_6, rapflux_resid_7, rapflux_resid_8, rapflux_ivar_1, rapflux_ivar_2, rapflux_ivar_3, rapflux_ivar_4, rapflux_ivar_5, rapflux_ivar_6, rapflux_ivar_7, rapflux_ivar_8, iapflux_1, iapflux_2, iapflux_3, iapflux_4, iapflux_5, iapflux_6, iapflux_7, iapflux_8, iapflux_resid_1, iapflux_resid_2, iapflux_resid_3, iapflux_resid_4, iapflux_resid_5, iapflux_resid_6, iapflux_resid_7, iapflux_resid_8, iapflux_ivar_1, iapflux_ivar_2, iapflux_ivar_3, iapflux_ivar_4, iapflux_ivar_5, iapflux_ivar_6, iapflux_ivar_7, iapflux_ivar_8, zapflux_1, zapflux_2, zapflux_3, zapflux_4, zapflux_5, zapflux_6, zapflux_7, zapflux_8, zapflux_resid_1, zapflux_resid_2, zapflux_resid_3, zapflux_resid_4, zapflux_resid_5, zapflux_resid_6, zapflux_resid_7, zapflux_resid_8, zapflux_ivar_1, zapflux_ivar_2, zapflux_ivar_3, zapflux_ivar_4, zapflux_ivar_5, zapflux_ivar_6, zapflux_ivar_7, zapflux_ivar_8, yapflux_1, yapflux_2, yapflux_3, yapflux_4, yapflux_5, yapflux_6, yapflux_7, yapflux_8, yapflux_resid_1, yapflux_resid_2, yapflux_resid_3, yapflux_resid_4, yapflux_resid_5, yapflux_resid_6, yapflux_resid_7, yapflux_resid_8, yapflux_ivar_1, yapflux_ivar_2, yapflux_ivar_3, yapflux_ivar_4, yapflux_ivar_5, yapflux_ivar_6, yapflux_ivar_7, yapflux_ivar_8 ) VALUES ( %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s )'  
##
   cursor.execute( query, tuple( [str(elem) for elem in [ id ] + newdata_decap[i]] ) )
##
    
#
#
con.commit()

print fimage, 'has been loaded'

# That's it!

