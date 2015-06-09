"""
Tools for DESI target selection from DECALS Database.
"""

import psycopg2
import sys
import numpy as np
from astropy.io import fits
import os


import sys; sys.path.insert(0, '/global/project/projectdirs/m779/yfeng1/source/imaginglss')

from imaginglss.analysis import cuts

#dictionaries with general target info

target_types = ['ELG', 'LRG', 'QSO', 'STDSTAR', 'GAL', 'OTHER']
target_priority = {'ELG': 4, 'LRG': 3, 'QSO': 1}
target_nobs = {'ELG':1, 'LRG':2, 'QSO': 2}

def targetsintile(ra=330, dec=0.0, radius=0.0):
    """
    Selects the targets in a circular tile by making a DB query.
    Uses the selection cuts defined in imaginglss.analysis.

    Parameters
    ----------
    ra : RA for the tile's center.
    dec : DEC for the tile's center.
    radius : radius of the circular tile, in degrees.

    Returns 
    -------

    target_ra : numpy array
                RA for the targets.
    target_dec : numpy array
                DEC for the targets
    target_type: numpy float array
                target type
                
    """

    #open DB session
    con = psycopg2.connect(host='scidb2.nersc.gov', user='desi_user', database='desi')
    cur = con.cursor()
    
    #execute
    cur.execute("select candidate.id, candidate.ra, candidate.dec, decam.gflux, decam.rflux, decam.zflux, wise.w1flux, wise.w2flux, decam.g_ext, decam.r_ext, decam.z_ext, wise.w1_ext, wise.w2_ext from  candidate, decam, wise where q3c_radial_query(candidate.ra, candidate.dec, %f, %f, %f) and decam.cand_id = candidate.id and wise.cand_id = candidate.id and decam.g_anymask =0  and decam.r_anymask =0 and decam.z_anymask=0;"%(ra,dec,radius))
    
    m=cur.fetchall()
    data = np.array(m)

    data_dic = dict([('ID', data[:,0]), ('RA', data[:,1]), ('DEC', data[:,2]), 
                     ('GFLUX', data[:,3]), ('RFLUX', data[:,4]), ('ZFLUX', data[:,5]), 
                     ('W1FLUX', data[:,6]), ('W2FLUX', data[:,7]),
                     ('GFRAC', data[:,8]), ('RFRAC',data[:,9]), ('ZFRAC', data[:,10]), 
                     ('W1FRAC', data[:,11]),  ('W2FRAC', data[:,12])])

    not_zero = np.where((data_dic['GFLUX']!=0) & (data_dic['RFLUX']!=0) & (data_dic['ZFLUX']!=0) 
                    & (data_dic['W1FLUX']!=0) & (data_dic['W2FLUX']!=0))

    not_zero = not_zero[0]


    print np.size(not_zero)
    if con:
        con.close()

