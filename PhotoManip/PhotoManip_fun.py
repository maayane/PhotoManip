"""*******************************************************
Photomanip:
-conversion of data from the Marshal into the following format: 'jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter','instr'
-calculation of the epoch t0 at which the flux in a chosen band crosses zero (often used as explosion date)
******************************************************"""
print(__doc__)

__author__ = 'maayanesoumagnac'

import numpy as np
import pdb
import csv
import os
from scipy.integrate import trapz
import pylab
from . import fitter_general
import pandas as pd

__all__=['magAB_in_filter_to_flux_in_filter','tref_from_P48','read_data_Marshall_simple']

def magAB_in_filter_to_flux_in_filter(mag_in_filter,Filter_vector=None,Filter_object=None,filters_directory=None,verbose=False):
    """Description: given an array of apparent magnitude_AB in a given filter, calculates the array of fluxes, ith the equation
        m_AB(P)=-2.5logf(P)-2.5log(lambda_P^2/c)-48.6
        Input  :- an apparent mag or an array of apparent mag
                - either
                    Filter vector: N-long array of of arrays [[filter family, filtername],[filter family, filtername],[filter family, filtername],...etc]
            the family and names are the ones provided by the pyphot library
                    or Filter object: a filter defined as in pyphot, e.g. P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon',
                  unit='Angstrom')
                  The allowed combinations in Filter_vectors are: ['ztf_p48','r_p48'],['ztf_p48','g_p48'],['ztf_p48','i_p48'],['sdss','z_sdss'],['sdss','r_sdss'],['sdss','g_sdss'],['sdss','i_sdss'],
            ['sdss','u_sdss'],['swift','uvw1'],['swift','uvw2'],['swift','uvm2'],['swift','u_swift'],['swift','v_swift'],['swift','b_swift'],['galex','galex_nuv'],['2mass','j_2mass'],['2mass','h_2mass'],
            ['2mass','k_2mass'],['cousin','r_cousin'],['cousin','i_cousin'],['johnson','u_johnson'],['johnson','b_johnson'],['johnson','v_johnson']
                - the path to the filters directory
        Output : flux value or array of fluxes in erg/s/AA/cm^2
        Tested : yes
             By : Maayane T. Soumagnac Nov 2016
        Example:magAB_in_filter_to_flux_in_filter(mag_value,Filter_vector=np.array([['sdss','SDSS_g']]))"""
    #print('filters directory,',filters_directory)
    fluxes=np.zeros(np.shape(mag_in_filter))
    tens = 10.0 * np.ones(np.shape(mag_in_filter))
    if Filter_vector is not None:
        for i,j in enumerate(Filter_vector):
            [P,wav]=make_filter_object([j],filters_directory=filters_directory)
            if verbose==True:
                print("The filter AB_zero_flux is",P['filter_object'][0].AB_zero_flux())
                print("The filter AB_zero_mag is", P['filter_object'][0].AB_zero_mag)
                print("The filter lambda_pivot is", P['filter_object'][0]._lpivot)
            fluxes=P['filter_object'][0].AB_zero_flux()*np.power(tens,-mag_in_filter/2.5)
    elif Filter_object is not None:
        if isinstance(Filter_object, (list,)) and len(Filter_object)==1:
            if verbose == True:
                print('The filter AB_zero_flux is',Filter_object[0].AB_zero_flux())
            fluxes = Filter_object[0].AB_zero_flux() * np.power(tens, -mag_in_filter / 2.5)
        else:
            if verbose == True:
                print('Filter_object.AB_zero_flux is', Filter_object.AB_zero_flux())
            fluxes = Filter_object.AB_zero_flux() * np.power(tens, -mag_in_filter / 2.5)
    else:
        raise Exception('ERROR: you need to define a filter, either a Filter_vector or a filter object')
        #pdb.set_trace()
    return fluxes

def error_on_mag_into_error_on_flux(magerr,flux):
    if isinstance(magerr,np.ndarray):
        fluxerr=np.abs(-2.303/2.5*magerr*flux) ##1.0855405992184108=2.5/2.303 car ln=2.303log
    else:
        fluxerr = abs(-2.303 / 2.5 * magerr * flux)
    return fluxerr

def dico_transmissions(filters_directory):
    #print('I am loading the dic')
    #pdb.set_trace()
    dic_transmissions={}
    dic_transmissions['ptf_g_p48']=np.genfromtxt(filters_directory + '/PTF_G.rtf', delimiter=None)
    dic_transmissions['ptf_r_p48']=np.genfromtxt(filters_directory + '/P48_R_T.rtf', delimiter=None)
    dic_transmissions['ztf_g_p48']=np.genfromtxt(filters_directory + '/ZTF_g_fromgithub_AA.txt', delimiter=None)
    dic_transmissions['ztf_r_p48']=np.genfromtxt(filters_directory + '/ZTF_r_fromgithub_AA.txt', delimiter=None)
    dic_transmissions['ztf_i_p48']=np.genfromtxt(filters_directory + '/ZTF_i_fromgithub_AA_reorder.txt',
                                             delimiter=None)
    dic_transmissions['uvw1_swift']=np.genfromtxt(filters_directory + '/Swift_UVW1.rtf', delimiter=None)
    dic_transmissions['uvw2_swift']=np.genfromtxt(
                    filters_directory + '/Swift_UVW2.rtf', delimiter=None)
    dic_transmissions['uvm2_swift'] =np.genfromtxt(
        filters_directory + '/Swift_UVM2.rtf', delimiter=None)
    dic_transmissions['u_swift']=np.genfromtxt(
        filters_directory + '/Swift_u.rtf', delimiter=None)
    dic_transmissions['v_swift'] =np.genfromtxt(
        filters_directory + '/Swift_V.rtf', delimiter=None)
    dic_transmissions['b_swift'] =np.genfromtxt(
        filters_directory + '/Swift_B.rtf', delimiter=None)
    dic_transmissions['nuv_galex']=np.genfromtxt(
                    filters_directory + '/GALEX_NUV.dat', delimiter=None)
    dic_transmissions['r_sdss']=np.genfromtxt(
        filters_directory + '/SDSS_r.txt', delimiter=',')
    dic_transmissions['z_sdss']=np.genfromtxt(filters_directory + '/SDSS_z.txt', delimiter=',')
    dic_transmissions['i_sdss'] = np.genfromtxt(
        filters_directory + '/SDSS_i.txt', delimiter=',')
    dic_transmissions['g_sdss'] = np.genfromtxt(
        filters_directory + '/SDSS_g.txt', delimiter=',')
    dic_transmissions['u_sdss'] = np.genfromtxt(
        filters_directory + '/SDSS_u.txt', delimiter=',')

    dic_transmissions['j_2mass'] =np.genfromtxt(filters_directory + '/2MASS_J.txt',
                  delimiter=',')
    dic_transmissions['h_2mass']=np.genfromtxt(filters_directory + '/2MASS_H.txt',
                  delimiter=',')
    dic_transmissions['k_2mass']=np.genfromtxt(
        filters_directory + '/2MASS_K.txt',
        delimiter=',')
    dic_transmissions['i_cousin']=np.genfromtxt(filters_directory + '/cousin_i.txt',
                  delimiter=',')
    dic_transmissions['r_cousin'] =np.genfromtxt(filters_directory + '/cousin_r.txt',
                  delimiter=',')
    dic_transmissions['u_johnson']=np.genfromtxt(filters_directory + '/johnson_u.txt',
                  delimiter=',')
    dic_transmissions['b_johnson']=np.genfromtxt(filters_directory + '/johnson_b.txt',
                  delimiter=',')
    dic_transmissions['v_johnson'] =np.genfromtxt(
        filters_directory + '/johnson_v.txt',
        delimiter=',')
    return dic_transmissions

def make_filter_object(Filter_vector, dic_transmissions=None, central=True, verbose=False,filters_directory=None):
    """Description: from a filter vector where each element is a couple [filter family,filter name], create a filter object P as in pyphoy
        Input  :- a filter vector: can be given in two shapes:

        OPTION 1:
        Filter_vector = np.empty([2, 2], dtype=object)
        Filter_vector[0] = [str('GALEX'), str('GALEX_NUV')]
        Filter_vector[1]=[str('ptf_p48'),str('r_p48')]

        OPTION 2:
        Filter_vector_2=[['swift','UVW1'],['ztf_p48','p48_r']]
                - dic tranmission: a dictionnary where all the transmission curves have been loaded as nparrays
                - central. If true, gives pyphot .cl wavelength, which corresponds to Eran's AstFilter.get('family','band').eff_wl
                            else, gives phyphot .eff wavelength, which I am not sure what it is..
        Output :- a dictionnary P where
            P['filter_object'] is a list  with all the filters,
            P['filtername'] is a numpy array with all the filters names
            P['filter_family'] is a numpy array with all the families
        with the corresponding data
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example:Filter_vector = np.empty([2, 2], dtype=object)
                Filter_vector[0] = [str('GALEX'), str('GALEX_NUV')]
                Filter_vector[1]=[str('ptf_p48'),str('r')]
                [P, wav]=make_filter_object(Filter_vector)
        Reliable:  """
    wavelength_filter_effective = dict()  # np.empty(np.shape(Filter_vector)[0])
    wavelength_filter_central = dict()
    wavelength_filter_pivot = dict()
    P_vector = dict()
    #print('Filter_vector is', Filter_vector)
    #print(isinstance(Filter_vector, (list,)))
    # pdb.set_trace()
    if isinstance(Filter_vector, (list,)):
        Filter_vectorx = np.empty([len(Filter_vector), 2], dtype=object)
        for i, j in enumerate(Filter_vector):
            Filter_vectorx[i, 0] = Filter_vector[i][0]
            Filter_vectorx[i, 1] = Filter_vector[i][1]
        P_vector['filter_family'] = Filter_vectorx[:, 0]
        P_vector['filter_name'] = Filter_vectorx[:, 1]
        P_vector['filter_object'] = []

    else:
        Filter_vectorx = Filter_vector
        P_vector['filter_family'] = Filter_vector[:, 0]
        P_vector['filter_name'] = Filter_vector[:, 1]
        P_vector['filter_object'] = []
    # print("P_vector['filter_object'] is",P_vector['filter_object'])
    # pdb.set_trace()
    if dic_transmissions is None:
        if filters_directory is None:
            raise Exception('ERROR you need to give either a dic_transmissions or a filters_directory variable')
            #pdb.set_trace()
        else:
            #print('filters_directory',filters_directory)
            dic_transmissions=dico_transmissions(filters_directory)
    for i, j in enumerate(Filter_vectorx):
        #print('j is', j)
        if j[0].lower() == 'ptf_p48':
            # print(j[1])
            if j[1].lower() == 'g_p48':
                # if verbose == True:
                print('You gave the G filter of the PTF_P48 family')
                # pdb.set_trace()
                Transmission = dic_transmissions['ptf_g_p48']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_G', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'r_p48':
                # if verbose==True:
                print('You gave the R filter of the PTF_P48 family')
                Transmission = dic_transmissions['ptf_r_p48']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'ztf_p48':
            # print(j[1])
            if j[1].lower() == 'g_p48':
                if verbose == True:
                    print('You gave the G filter of the ztf_P48 family')
                # pdb.set_trace()
                Transmission = dic_transmissions['ztf_g_p48']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='ZTF_P48_G', dtype='photon',
                                  unit='Angstrom')
                # print('P is',P)
                # pdb.set_trace()
            elif j[1].lower() == 'r_p48':
                if verbose == True:
                    print('You gave the R filter of the ZTF_P48 family')
                Transmission = dic_transmissions['ztf_r_p48']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='ZTF_P48_R', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'i_p48':
                if verbose == True:
                    print('You gave the I filter of the ZTF_P48 family')
                Transmission = dic_transmissions['ztf_i_p48']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='ZTF_P48_I', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'swift':
            if j[1].lower() == 'uvw1':
                if verbose == True:
                    print('You gave the uvw1 filter of the swift family')
                # pdb.set_trace()
                Transmission = dic_transmissions['uvw1_swift']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_UVW1', dtype='photon',
                                  unit='Angstrom')
                # elif j[0].lower() == 'swift':
            elif j[1].lower() == 'uvw2':
                if verbose == True:
                    print('You gave the uvw2 filter of the swift family')
                # pdb.set_trace()
                Transmission = dic_transmissions['uvw2_swift']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_UVW2', dtype='photon',
                                  unit='Angstrom')
                # elif j[0].lower() == 'swift':
            elif j[1].lower() == 'uvm2':
                if verbose == True:
                    print('You gave the uvm2 filter of the swift family')
                # pdb.set_trace()
                Transmission = dic_transmissions['uvm2_swift']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_UVM2', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'u_swift':
                if verbose == True:
                    print('You gave the u filter of the swift family')
                # pdb.set_trace()
                Transmission = dic_transmissions['u_swift']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_u', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'v_swift':
                if verbose == True:
                    print('You gave the v filter of the swift family')
                # pdb.set_trace()
                Transmission = dic_transmissions['v_swift']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_V', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'b_swift':
                if verbose == True:
                    print('You gave the b filter of the swift family')
                # pdb.set_trace()
                Transmission = dic_transmissions['b_swift']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_B', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'galex':
            if j[1].lower() == 'galex_nuv':
                if verbose == True:
                    print('You gave the nuv filter of the galex family')
                # pdb.set_trace()
                Transmission = dic_transmissions['nuv_galex']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='galex_nuv', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'sdss':
            if verbose == True:
                print('You gave the sdss family')
            # print(j[1].lower())
            if j[1].lower() == 'r_sdss':
                if verbose == True:
                    print('You gave the r filter of the sdss family')
                # pdb.set_trace()
                Transmission = dic_transmissions['r_sdss']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='r_sdss', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'g_sdss':
                if verbose == True:
                    print('You gave the g filter of the sdss family')
                # pdb.set_trace()
                Transmission = dic_transmissions['g_sdss']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='g_sdss', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'i_sdss':
                if verbose == True:
                    print('You gave the i filter of the sdss family')
                # pdb.set_trace()
                Transmission = dic_transmissions['i_sdss']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='i_sdss', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'u_sdss':
                if verbose == True:
                    print('You gave the u filter of the sdss family')
                # pdb.set_trace()
                Transmission = dic_transmissions['u_sdss']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='u_sdss', dtype='photon',
                                  unit='Angstrom')

            elif j[1].lower() == 'z_sdss':
                if verbose == True:
                    print('You gave the z filter of the sdss family')
                # pdb.set_trace()
                Transmission = dic_transmissions['z_sdss']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='z_sdss', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == '2mass':
            if verbose == True:
                print('You gave the 2MASS family')
                print(j[1].lower())
            if j[1].lower() == 'j_2mass':
                if verbose == True:
                    print('You gave the J filter of the 2MASS family')
                    # pdb.set_trace()
                Transmission = dic_transmissions['j_2mass']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='2MASS_J', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'h_2mass':
                if verbose == True:
                    print('You gave the h filter of the 2MASS family')
                    # pdb.set_trace()
                Transmission = dic_transmissions['h_2mass']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='h_2mass', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'k_2mass':
                if verbose == True:
                    print('You gave the k filter of the 2MASS family')
                    # pdb.set_trace()
                Transmission = dic_transmissions['k_2mass']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='k_2mass', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'cousin':
            if verbose == True:
                print('You gave the cousin family')
                # print(j[1].lower())
            if j[1].lower() == 'i_cousin':
                if verbose == True:
                    print('You gave the i filter of the cousin family')
                    # pdb.set_trace()
                Transmission = dic_transmissions['i_cousin']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='i_counsin', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'r_cousin':
                if verbose == True:
                    print('You gave the r filter of the cousin family')
                Transmission = dic_transmissions['r_cousin']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='r_cousin', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'johnson':
            if verbose == True:
                print('You gave the johnson family')
                print(j[1].lower())
            if j[1].lower() == 'u_johnson':
                if verbose == True:
                    print('You gave the u filter of the johnson family')
                    # pdb.set_trace()
                Transmission = dic_transmissions['u_johnson']
                # print('the shape of transission is',Transmission)
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='u_johnson', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'b_johnson':
                if verbose == True:
                    print('You gave the b filter of the johnson family')
                    # pdb.set_trace()
                Transmission = dic_transmissions['b_johnson']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='b_johnson', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'v_johnson':
                if verbose == True:
                    print('You gave the v filter of the johnson family')
                    # pdb.set_trace()
                Transmission = dic_transmissions['v_johnson']
                P = Filter(Transmission[:, 0], Transmission[:, 1], name='v_johnson', dtype='photon',
                                  unit='Angstrom')
        else:
            print('I HAVE NOT RECOGNIZE THE FILTER')
            pdb.set_trace()
            #lib = pyphot.get_library()
            #f = lib.find(j[0].lower())  # filter family
            #if verbose == True:
            #    for name in f:
            #        lib[name].info(show_zeropoints=True)
            #P = lib[j[1]]  # filter name
            # min_w=np.min(Transmission[:,0])
            # max_w=np.max(Transmission[:,0])
            # print(P_vector)
            # print('Filter_vector is',Filter_vector)
            # print('P is',P)
            # print("P_vector['filter_object'] is",P_vector['filter_object'])
            # if isinstance(P_vector['filter_object'],(list,)) is False:
            #   P_vector['filter_object']=[P_vector['filter_object']].append(P)
            # print("P_vector['filter_object'] is", P_vector['filter_object'])
            # print("P_vector['filter_name'] is", P_vector['filter_name'])
        if isinstance(P_vector['filter_object'], (list,)) is True:
            #print('oui, liste de longueur', len(P_vector['filter_object']))
            if len(P_vector['filter_object']) > 0:
                ##print('oui')
                # print('I am trying to append {0} to {1}'.format(P,P_vector['filter_object']))
                P_vector['filter_object'] = P_vector['filter_object'] + [P]
                P_vector['filter_name'] = P_vector['filter_name'] + [j[1].lower()]
                P_vector['filter_family'] = P_vector['filter_family'] + [j[0].lower()]
                # print('ca a marche?',P_vector['filter_object'])
            else:
                # print('faison une liste')
                P_vector['filter_object'] = [P]
                # print(P_vector['filter_object'])
                P_vector['filter_name'] = [j[1].lower()]
                P_vector['filter_family'] = [j[0].lower()]
                # pdb.set_trace()

                # print("P_vector['filter_object'] is",P_vector['filter_object'])
                # print("P_vector['filter_name'] is",P_vector['filter_name'])
                # print('i is',i)
        #wavelength_filter_effective[P_vector['filter_name'][i]] = P.leff.item()
        wavelength_filter_central[P_vector['filter_name'][i]] = P.cl()
        # wavelength_filter_min[P_vector['filter_name'][i]]=np.min(Transmission[:,0])
        # wavelength_filter_max[P_vector['filter_name'][i]] = np.max(Transmission[:, 0])

        # if isinstance(Filter_vector, (list,)):
        # P_vector['filter_object']=P_vector['filter_object'][0]

        # print(P_vector)
        # print(P)
        # P_vector['filter_object'].append(P)
        # wavelength_filter_effective[P_vector['filter_name'][i]]=P.leff.item()
        # wavelength_filter_central[P_vector['filter_name'][i]]=P.cl.item()
        # if isinstance(Filter_vector, (list,)):
        #    P_vector['filter_object']=P_vector['filter_object'][0]

    # print('Filter_vector was', Filter_vector)
    # print(' and P_vector is', P_vector)
    # pdb.set_trace()

    if central == True:
        return P_vector, wavelength_filter_central  # same as Eran eff_wl
    else:
        return P_vector, wavelength_filter_effective

class Filter:
    def __init__(self, wavelength, transmit, name='', dtype="photon",unit='Angstrom'):
        """a class inspired by the pyphot class with same name by Morgan Fouesneau"""
        self.name       = name
        self.set_dtype(dtype)
        self._wavelength = wavelength
        self.transmit   = np.clip(transmit, 0., np.nanmax(transmit))
        self.norm       = trapz(self.transmit, self._wavelength)#int(T dlambda)
        self._lT        = trapz(self._wavelength * self.transmit, self._wavelength)#int(lambda*T dlambda)
        self._lpivot    = self.calculate_lpivot()
        self._cl        = self._lT / self.norm
        self.unit=unit
    def set_dtype(self, dtype):
        """ Set the detector type (photon or energy)"""
        _d = dtype.lower()
        if "phot" in _d:
            self.dtype = "photon"
        elif "ener" in _d:
            self.dtype = "energy"
        else:
            raise ValueError('Unknown detector type {0}'.format(dtype))
    def calculate_lpivot(self):
        if 'photon' in self.dtype:
            lpivot2 = self._lT / trapz(self.transmit / self._wavelength, self._wavelength)
        else:
            lpivot2 = self.norm / trapz(self.transmit / self._wavelength ** 2, self._wavelength)
        return np.sqrt(lpivot2)
    def cl(self):
        return self._cl
    def info(self, show_zeropoints=True):
        """ display information about the current filter"""
        print("""Zeropoints
          AB: {0} magAB,
              {1} erg/s/cm^2/AA
        """.format(self.AB_zero_mag(),self.AB_zero_flux()))
    def get_flux(self, slamb, sflux, axis=-1,verbose=True):
        """getFlux
        Integrate the flux within the filter and return the integrated energy
        Parameters
        ----------
        slamb: ndarray(dtype=float, ndim=1)
            spectrum wavelength definition domain
        sflux: ndarray(dtype=float, ndim=1)
            associated flux
        Returns
        -------
        flux: float
            Energy of the spectrum within the filter
        """
        if verbose==True:
            print('in Filter.get_flus, WARNING: T and the flux wavelengths should be given in AA')
        _slamb=slamb
        # clean data for inf values by interpolation
        if True in np.isinf(sflux):
            indinf = np.where(np.isinf(sflux))
            indfin = np.where(np.isfinite(sflux))
            sflux[indinf] = np.interp(_slamb[indinf], _slamb[indfin], sflux[indfin])

        # reinterpolate transmission onto the same wavelength def as the data
        ifT = np.interp(_slamb, self._wavelength, self.transmit, left=0., right=0.)
        # if the filter is null on that wavelength range flux is then 0
        # ind = ifT > 0.
        nonzero = np.where(ifT > 0)[0]
        if nonzero.size <= 0:
            return 0.
        nonzero_start = max(0, min(nonzero) - 5)
        nonzero_end = min(len(ifT), max(nonzero) + 5)
        ind = np.zeros(len(ifT), dtype=bool)
        ind[nonzero_start:nonzero_end] = True
        if True in ind:
            try:
                _sflux = sflux[:, ind]
            except:
                _sflux = sflux[ind]
            # limit integrals to where necessary
            # ind = ifT > 0.
            if 'photon' in self.dtype:
                a = np.trapz(_slamb[ind] * ifT[ind] * _sflux, _slamb[ind], axis=axis)
                b = np.trapz(_slamb[ind] * ifT[ind], _slamb[ind])
            elif 'energy' in self.dtype:
                a = np.trapz(ifT[ind] * _sflux, _slamb[ind], axis=axis)
                b = np.trapz(ifT[ind], _slamb[ind])
            if (np.isinf(a).any() | np.isinf(b).any()):
                print(self.name, "Warn for inf value")
            return a / b
        else:
            return 0.
    def AB_zero_mag(self):
        """ AB magnitude zero point
        ABmag = -2.5 * log10(f_nu) - 48.60
              = -2.5 * log10(f_lamb) - 2.5 * log10(lpivot ** 2 / c) - 48.60
              = -2.5 * log10(f_lamb) - zpts
        """
        #C1 = unit[self.wavelength_unit].to('AA').magnitude ** 2 / unit['c'].to('AA/s').magnitude
        c1 = self._lpivot ** 2 /(2.99792458e18)
        m = 2.5 * np.log10(c1) + 48.6
        return m
    def AB_zero_flux(self):
        """ AB flux zero point in erg/s/cm2/AA """
        return 10 ** (-0.4 * self.AB_zero_mag())
    def get_mag(self, slamb, sflux):
        print('in Filter.get_mag, WARNING: T and the flux wavelengths should be given in AA')
        return -2.5 * np.log10(self.get_flux(slamb, sflux) / self.AB_zero_flux())

def read_data_Marshall_simple(path,no99=False,filters_directory=None,output_path=None):
    """Description: Reads data from a file downloaeded from the Marshall and outputs as many dictionnaries as filters,
    with the relevant data in them
        Input  :- a path to an ascii file downloaded from the Marshall
                - plot_all: if True, plots all the light curves for each filter
        Output :-data_dict: a dictionnary with keys ['jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter','instr'] and values a numpy array
        with the corresponding data
        Tested : yes
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example:
        Reliable:  """
    data_dic=dict()
    #from the Marshal "export light curve"
    #print('path is',path)
    print('**** I am converting the Marshall data into the right format ****')
    with open(path) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',',skipinitialspace=True)
            headers = next(readCSV)#skip the header
            dates_mjd = []
            filt = []
            absmag = []
            mag = []
            magerr = []
            instr = []
            for row in readCSV:
                date_mjd = float(row[1])
                fi = row[2]
                absm = float(row[3])
                ma = float(row[4])
                maerr = float(row[5])
                inst = row[7]
                dates_mjd.append(date_mjd)
                filt.append(fi)
                absmag.append(absm)
                mag.append(ma)
                magerr.append(maerr)
                instr.append(inst)
    dates=np.array(dates_mjd,dtype=float)
    filt=np.array(filt,dtype=object)
    absmag=np.array(absmag,dtype=float)
    mag=np.array(mag,dtype=float)
    magerr=np.array(magerr,dtype=float)
    instr=np.array(instr,dtype=str)

    ############### Swift ###############
    data_full=dict()
    if no99==True:
        print('no99 is true: I am removing 99 mags. Set this parameter to False if you want to leave 99 mags.')
        data_full['mag']=mag[mag!=99.0]
        data_full['jd']=dates[mag!=99.0]
        data_full['magerr'] = magerr[mag!=99.0]
        data_full['instr']= instr[mag!=99.0]
        data_full['filter']=filt[mag!=99.0]
        data_full['absmag']=absmag[mag!=99.0]
        data_full['absmagerr']=magerr[mag!=99.0]
        #data_full['absmagerr']=absmagerr[mag!=99.0]
    else:
        data_full['mag']=mag
        data_full['jd']=dates
        data_full['magerr'] = magerr
        data_full['instr']= instr
        data_full['filter']=filt
        data_full['absmag']=absmag
        data_full['absmagerr']=magerr
    #print('data_full[filter] is',data_full['filter'])
    #print(type(data_full['filter']))
    #pdb.set_trace()
    for i,j in enumerate(data_full['filter']):
        if (data_full['filter'][i]=='r') & (data_full['instr'][i]=='P48+ZTF'):
            data_full['filter'][i]='r_p48'
        if (data_full['filter'][i]=='g') & (data_full['instr'][i]=='P48+ZTF'):
            data_full['filter'][i]='g_p48'
        if (data_full['filter'][i]=='i') & (data_full['instr'][i]=='P48+ZTF'):
            data_full['filter'][i]='i_p48'
        if (data_full['filter'][i]=='r') & (data_full['instr'][i]=='P60+SEDM'):
            data_full['filter'][i]='r_sdss'
        if (data_full['filter'][i]=='g') & (data_full['instr'][i]=='P60+SEDM'):
            data_full['filter'][i]='g_sdss'
        if (data_full['filter'][i]=='i') & (data_full['instr'][i]=='P60+SEDM'):
            data_full['filter'][i]='i_sdss'
        if (data_full['filter'][i]=='u') & (data_full['instr'][i]=='P60+SEDM'):
            data_full['filter'][i]='u_sdss'
        if (data_full['filter'][i].lower()=='u') & (data_full['instr'][i]in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='u_swift'
        if (data_full['filter'][i].lower()=='uuu') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='u_swift'
        if (data_full['filter'][i].lower()=='v') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='v_swift'
        if (data_full['filter'][i].lower()=='uvv') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='v_swift'
        if (data_full['filter'][i].lower()=='b') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='b_swift'
        if (data_full['filter'][i].lower()=='ubb') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='b_swift'
        if (data_full['filter'][i].lower()=='uw1') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='uvw1'
        if (data_full['filter'][i].lower()=='uw2') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='uvw2'
        if (data_full['filter'][i].lower()=='um2') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            data_full['filter'][i]='uvm2'
        if (data_full['filter'][i].lower()=='r') & (data_full['instr'][i]=='LT+IOO'):
            data_full['filter'][i]='r_sdss'
        if (data_full['filter'][i].lower()=='g') & (data_full['instr'][i]=='LT+IOO'):
            data_full['filter'][i]='g_sdss'
        if (data_full['filter'][i].lower()=='i') & (data_full['instr'][i]=='LT+IOO'):
            data_full['filter'][i]='i_sdss'
        if (data_full['filter'][i].lower()=='u') & (data_full['instr'][i]=='LT+IOO'):
            data_full['filter'][i]='u_sdss'
        if (data_full['filter'][i].lower()=='z') & (data_full['instr'][i]=='LT+IOO'):
            data_full['filter'][i]='z_sdss'



    data_full['flux']=np.zeros(np.shape(data_full['mag'])[0])
    for i,j in enumerate(data_full['filter']):
        if (data_full['filter'][i]=='r_p48') & (data_full['instr'][i]=='P48+ZTF'):
            [P, wav] = make_filter_object([[str('ztf_p48'),'r_p48']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='g_p48') & (data_full['instr'][i]=='P48+ZTF'):
            [P, wav] = make_filter_object([[str('ztf_p48'), 'g_p48']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='i_p48') & (data_full['instr'][i]=='P48+ZTF'):
            [P, wav] = make_filter_object([[str('ztf_p48'), 'i_p48']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='r_sdss') & (data_full['instr'][i]=='P60+SEDM'):
            [P, wav] = make_filter_object([[str('sdss'), 'r_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='g_sdss') & (data_full['instr'][i]=='P60+SEDM'):
            [P, wav] = make_filter_object([[str('sdss'), 'g_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='i_sdss') & (data_full['instr'][i]=='P60+SEDM'):
            [P, wav] = make_filter_object([[str('sdss'), 'i_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='u_sdss') & (data_full['instr'][i]=='P60+SEDM'):
            [P, wav] = make_filter_object([[str('sdss'), 'u_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i].lower()=='u_swift') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            [P, wav] = make_filter_object([[str('swift'), 'u_swift']],filters_directory=filters_directory)
        if (data_full['filter'][i].lower()=='v_swift') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            [P, wav] = make_filter_object([[str('swift'), 'v_swift']],filters_directory=filters_directory)
        if (data_full['filter'][i].lower()=='b_swift') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            [P, wav] = make_filter_object([[str('swift'), 'b_swift']],filters_directory=filters_directory)
        if (data_full['filter'][i].lower()=='uvw1') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            [P, wav] = make_filter_object([[str('swift'), 'uvw1']],filters_directory=filters_directory)
        if (data_full['filter'][i].lower()=='uvw2') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            [P, wav] = make_filter_object([[str('swift'), 'uvw2']],filters_directory=filters_directory)
        if (data_full['filter'][i].lower()=='uvm2') & (data_full['instr'][i] in ['Swift+UVOT','SWIFT+UVOT']):
            [P, wav] = make_filter_object([[str('swift'), 'uvm2']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='u_sdss') & (data_full['instr'][i]=='LT+IOO'):
            [P, wav] = make_filter_object([[str('sdss'), 'u_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='r_sdss') & (data_full['instr'][i]=='LT+IOO'):
            [P, wav] = make_filter_object([[str('sdss'), 'r_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='g_sdss') & (data_full['instr'][i]=='LT+IOO'):
            [P, wav] = make_filter_object([[str('sdss'), 'g_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='i_sdss') & (data_full['instr'][i]=='LT+IOO'):
            [P, wav] = make_filter_object([[str('sdss'), 'i_sdss']],filters_directory=filters_directory)
        if (data_full['filter'][i]=='z_sdss') & (data_full['instr'][i]=='LT+IOO'):
            [P, wav] = make_filter_object([[str('sdss'), 'z_sdss']],filters_directory=filters_directory)

        data_full['flux'][i] = magAB_in_filter_to_flux_in_filter(data_full['mag'][i], Filter_vector=None,
                                                                 Filter_object=P['filter_object'])#.magnitude
    data_full['fluxerr']= error_on_mag_into_error_on_flux(data_full['magerr'],data_full['flux'])
    #data_full['absmag']=absmag
    #data_full['absmagerr']=magerr
    #array_to_save=np.array(list(zip(data_full['jd'],data_full['flux'],data_full['fluxerr'],data_full['mag'],data_full['magerr'],data_full['absmag'],data_full['absmagerr'],data_full['filter'])))
    #np.savetxt('data_photofit_format.txt',array_to_save,header='jd,flux,fluxerr,mag,magerr,absmag,absmagerr,filter')
    header = ['jd', 'mag', 'magerr', 'flux', 'fluxerr', 'absmag', 'absmagerr', 'filter', 'instr']
    if output_path is None:
        output_path='./data_photofit_format.txt'
    with open(output_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(i for i in header)
        for row in zip(data_full['jd'], data_full['mag'], data_full['magerr'], data_full['flux'],
                       data_full['fluxerr'], data_full['absmag'], data_full['absmagerr'],
                       data_full['filter'], data_full['instr']):
            writer.writerow(row)
    return data_full

def read_data_into_numpy_array(path,delimiter=',',header=False,skiprows=None,no_repeat_rows=True):
    """Description: numpy genfrom txt but with mixed types
        Input  :- a path to an ascii file
                - delimiter
                - header: if false, then the output is simply a numpy array with mized types
                         if true, then there are 3 outputs: the numpy array, a numpy array with the fields names, and a dictionnary of which the keys
                         are the fields names and the values are the columns
        Output :- if header=True: an array with the data, an array with the header,
        and a dictionnary where the keys are the header values and the values are the corresponding columns
                - if header=False: the data array only
        a numpy array with mixed types
        with the corresponding data
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example: Swift=read_data_from_file.read_data_into_numpy_array('Swift.txt',delimiter=' ',header=True)[0]
                header_array=read_data_from_file.read_data_into_numpy_array('Swift.txt',delimiter=' ',header=True)[1]

        Reliable:  """
    if header==True:
        pd_array = pd.read_csv(path, header=0,
                                     index_col=False, delimiter=delimiter,skiprows=skiprows)
        pd_numpy_array=pd_array.values
        #print(pd_numpy_array)
        #print(pd_numpy_array[:, :-2])
        #print('the shape of pd_numpy_array is',np.shape(pd_numpy_array))
        #unique, index = np.unique(pd_numpy_array.astype("<U22"), axis=0, return_index=True)
        #print(index)
        #print('the shape of index is',np.shape(index))
        #print('the indexes sorted are:',np.argsort(pd_numpy_array.astype("<U22"), axis=0)[:,0])
        #pdb.set_trace()
        #pdb.set_trace()
        if no_repeat_rows==True:
            unique, index = np.unique(pd_numpy_array.astype("<U22"),axis=0,return_index=True)
            pd_numpy_array_norepeat=pd_numpy_array[index]
            #print('the shape after unique is', np.shape(pd_numpy_array_norepeat))
            #pdb.set_trace()


        pd_keys=np.array(pd_array.keys())
        pd_dict=dict()
        for i,k in enumerate(pd_keys):
            if no_repeat_rows==True:
                pd_dict[k]=pd_numpy_array_norepeat[:,i]
            else:
                pd_dict[k]=pd_numpy_array[:,i]
        return pd_numpy_array,pd_keys,pd_dict



    else:
        pd_array = pd.read_csv(path, index_col=False, delimiter=delimiter,skiprows=skiprows)
        pd_numpy_array = pd_array.values
        return pd_numpy_array

class model_exponent_concav_pos_withref(object): #given eo,to,t,tref an array of eo*[1-exp(-[t-tref]/to)] (=0 for t<=tref)
    def __init__(self, eo,to,tref,t):
        self.eo = eo
        self.to = to
        self.tref=tref
        self.t=t
    def model_array(self):
        g = np.zeros((np.shape(self.t)[0], 2))
        g[:, 0] = self.t
        #print((-(self.t[self.t[:] - self.tref >= 0]-self.tref)/self.to).dtype)
        #print(np.exp((-(self.t[self.t[:]- self.tref >= 0]-self.tref)/self.to).astype(float)))
        g[self.t[:] - self.tref >= 0, 1] = self.eo*(1-np.exp(np.array(-(self.t[self.t[:] - self.tref >= 0]-self.tref)/self.to,dtype=np.float32)))
        g[self.t[:] - self.tref < 0, 1] = 0.
        return g

class model_powerlaw_withref(object):  # given a,n, an array of a(x-xref)^n
    def __init__(self, a, n,xref, x):
        self.a = a
        self.n = n
        self.x = x
        self.xref = xref
    def model_array(self):
        g = np.zeros((len(self.x), 2))
        g[:, 0] = self.x[:]
        g[self.x[:]-self.xref>=0, 1] = self.a * (self.x[self.x[:]-self.xref>=0]-self.xref)**self.n
        g[self.x[:] - self.xref < 0, 1] = 0.
        return g

def tref_from_P48(path_to_data, tref_assumed=None, tref_priors=None,band='r',days_rising=10, already_run_mcmc_exp=False,
                  already_run_calc_all_chis_exp=False, already_plot_exp=False, already_run_mcmc_power=False,
                  already_run_calc_all_chis_power=False, already_plot_power=False,nwalkers=200,
                  num_steps=1000,eo_prior=[0, 1e-14],to_prior=[0,50],a_prior=[1e-17, 1e-15],n_prior=[0,1]):
    """Description: fit the P48_R light curve with an exponent or a power law to deduce the explosion date
    Input  :- a path to a file where the data is in the following format, 'jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter', e.g. after running
    Read_data.read_data_Marshall_simple() on some Marshall extracted data.
            -the assumed explosion date. If None, it will be taken as the minimum of the data
            -the number of days to fit (i.e. the number of days during which the LC is rising)
            -the band to use. the default is 'r': P48_r. If set to g, will be 'p48_g'
            -already_run.... these are False by default, but you can set them to True if you ran the analysis once already and only want to see the output
            -number of walkers in mcmc
            -number of steps in mcmc
            -prior on eo (exponential law)
            -prior on to (exponential law)
            -prior on a (power law)
            -prior on n (power law)
    Output :-
    Tested : ?
         By : Maayane T. Soumagnac Nov 2019
        URL :
    Example: samples=fitter_general.emcee_n_param(ndim=2,model_nparam=model_lo_to,prior_param=[prior_alpha,prior_tref],data=my_data,uncertainties=errors,initial_conditions=[alpha_true,tref_true],flatchain_path='output_test_fitter_general/flatchain_test.txt',already_run=False)
    Reliable:  """
    if band=='r':
        if os.path.exists('results_tref_calculator_from_P48R'):
            print('the output directory,results_tref_calculator_from_P48R exists already')
        else:
            os.mkdir('results_tref_calculator_from_P48R')

        output_mcmc_exp = 'results_tref_calculator_from_P48R/exp'
        output_mcmc_power = 'results_tref_calculator_from_P48R/power'
    elif band=='g':
        if os.path.exists('results_tref_calculator_from_P48G'):
            print('the output directory,results_tref_calculator_from_P48G exists already')
        else:
            os.mkdir('results_tref_calculator_from_P48G')

        output_mcmc_exp = 'results_tref_calculator_from_P48G/exp'
        output_mcmc_power = 'results_tref_calculator_from_P48G/power'

    if os.path.exists(output_mcmc_exp) == False:
        os.mkdir(output_mcmc_exp)

    if os.path.exists(output_mcmc_power) == False:
        os.mkdir(output_mcmc_power)

    # print(dict_ZTF_Marshal)


    ########################################## plots the bolometric luminosity and R flux ###########################################

    dict_all = read_data_from_file.read_data_into_numpy_array(path_to_data, header=True)[
        2]  # ['jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter']

    print(dict_all)
    condition = (dict_all['mag'] != 99.0) & (dict_all['absmagerr'] != 99.0)

    # remove 99
    for j in dict_all.keys():
        if j not in ['mag', 'absmagerr']:
            print('the key is', j)
            print('the dimensions of data_dicts[j] is', np.shape(dict_all[j]))
            print("the dimension of data_dicts['mag'] is", np.shape(dict_all['mag']))
            print('the dimension of the boolean is', np.shape([dict_all['mag'] != 99.0]))
            dict_all[j] = dict_all[j][condition]
    dict_all['mag'] = dict_all['mag'][condition]
    dict_all['absmagerr'] = dict_all['absmagerr'][condition]

    data_P48 = dict()

    if band=='r':
        data_P48['days'] = dict_all['jd'][np.asarray(dict_all['filter']) == 'r_p48']
        data_P48['flux'] = dict_all['flux'][np.asarray(dict_all['filter']) == 'r_p48']
        data_P48['fluxerr'] = dict_all['fluxerr'][np.asarray(dict_all['filter']) == 'r_p48']
    elif band=='g':
        data_P48['days'] = dict_all['jd'][np.asarray(dict_all['filter']) == 'g_p48']
        data_P48['flux'] = dict_all['flux'][np.asarray(dict_all['filter']) == 'g_p48']
        data_P48['fluxerr'] = dict_all['fluxerr'][np.asarray(dict_all['filter']) == 'g_p48']

    pylab.figure()
    pylab.errorbar(data_P48['days'], data_P48['flux'], yerr=data_P48['fluxerr'])
    pylab.grid()
    pylab.show()
    # pdb.set_trace()
    # data_P48=np.zeros((np.shape(data_P48_R_days)[0],3))


    condition_rise = data_P48['days'] - np.min(data_P48['days']) <= days_rising
    data_P48_rise = dict()
    data_P48_rise['days'] = data_P48['days'][condition_rise]
    data_P48_rise['flux'] = data_P48['flux'][condition_rise]
    data_P48_rise['fluxerr'] = data_P48['fluxerr'][condition_rise]

    pylab.figure()
    pylab.plot(data_P48['days'], data_P48['flux'], 'bo', label='all the data')
    pylab.plot(data_P48_rise['days'], data_P48_rise['flux'], 'ro', label='data to be fitted')
    pylab.axvline(np.max(data_P48_rise['days']))
    pylab.grid()
    pylab.legend()
    pylab.show()

    ########################################## fit the rising chunk with an exponent (0 for x<xref) ###########################################


    my_data = np.zeros((np.shape(data_P48_rise['days'])[0], 2))
    my_data[:, 0] = data_P48_rise['days']
    my_data[:, 1] = data_P48_rise['flux']

    ########################################## fit the rising chunk with an exponent ###########################################

    eo_assumed = 1e-15
    to_assumed = 15.

    if tref_assumed == None:
        tref_assumed = np.min(data_P48_rise['days'])

    if tref_priors == None:
        prior_tref= np.array([tref_assumed - 20, tref_assumed + 1])
    else:
        prior_tref = np.array(tref_priors)


    prior_eo = np.array(eo_prior)
    prior_to = np.array(to_prior)



    pylab.figure()
    pylab.errorbar(data_P48_rise['days'], data_P48_rise['flux'], yerr=0.5 * data_P48_rise['fluxerr'], fmt='r.')
    pylab.plot(models.model_exponent_concav_pos_withref(eo=eo_assumed, to=to_assumed, tref=tref_assumed,
                                                        t=data_P48_rise['days']).model_array()[:, 0],
               models.model_exponent_concav_pos_withref(eo=eo_assumed, to=to_assumed, tref=tref_assumed,
                                                        t=data_P48_rise['days']).model_array()[:, 1])
    pylab.xlabel(r"time (JD)")
    pylab.ylabel(r"flux $[erg/s/cm^2/\AA]$")
    # pylab.axvline(JD_peak,color='blue',linestyle='--')
    pylab.title('An attempt to check the prior range')
    pylab.title('P48 R photometry and exponnential fit')
    pylab.axvline(tref_assumed)

    print("********** FIT WITH A POSITIVE CONCAVE EXPONENT ************")

    output_mcmc = output_mcmc_exp

    samples = fitter_general.emcee_n_param(3, models.model_exponent_concav_pos_withref,
                                           prior_param=[prior_eo, prior_to, prior_tref],
                                           data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                           initial_conditions=[eo_assumed, to_assumed, tref_assumed],
                                           nwalkers=nwalkers, num_steps=num_steps, flatchain_path=output_mcmc + '/flatchain.txt',
                                           already_run=already_run_mcmc_exp)

    best_exp = fitter_general.calc_best_fit_n_param(ndim=3, model_nparam=models.model_exponent_concav_pos_withref,
                                                    flatchain_path=output_mcmc + '/flatchain.txt',
                                                    data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                                    winners=100,
                                                    output_file_path=output_mcmc,
                                                    bounds=[prior_eo, prior_to, prior_tref],
                                                    already_run_calc_all_chis=already_run_calc_all_chis_exp,
                                                    show_plots=False)
    # print('reduced chi2 is',best_exp[-1]/(np.shape(my_data)[0]-3))
    # pdb.set_trace()
    if np.shape(my_data)[0] - 3 !=0:
        chi2_exp = best_exp[-1] / (np.shape(my_data)[0] - 3)
    bests = best_exp[:-1]
    if already_plot_exp != True:
        fitter_general.plot_opt_fit_n_param(3, models.model_exponent_concav_pos_withref, bests, my_data,
                                            flatchain_path=output_mcmc + '/flatchain.txt',
                                            uncertainties=data_P48_rise['fluxerr'], output_file_path=output_mcmc,
                                            xlabel=None, ylabel=None)

        triangle = fitter_general.plot_2D_distributions(
            flatchain_path=output_mcmc + '/flatchain.txt', bests=bests, title='triangle',
            output_file_path=output_mcmc, parameters_labels=['eo', 'to', 'tref'])
        #histos = fitter_general.plot_1D_marginalized_distribution(
        #    flatchain_path=output_mcmc + '/flatchain.txt', bests=bests,
        #    output_pdf_file_path=output_mcmc, output_txt_file_path=output_mcmc, parameters_labels=['eo', 'to', 'tref'],
        #    number_bins=1000)
        histos = fitter_general.plot_1D_marginalized_distribution_mini_interval(flatchain_path=output_mcmc + '/flatchain.txt', bests=bests, output_pdf_file_path=output_mcmc,
                                                        output_txt_file_path=output_mcmc,
                                                        parameters_labels=['eo', 'to', 'tref'], number_bins=1000)
    print('********************************************')
    print('when fitting whith the exponent, the best fit is {0}, with sigma is {1}'.format(bests[:], np.genfromtxt(
        output_mcmc + '/1sigma_mini.txt')[:, 2]))
    # print('when fitting whith the power law, the best fit is {0}, with sigma is {1}'.format(bests[:-1],np.genfromtxt('tref_calculator_results_powerlaw/1sigma.txt')[:,2])
    best_tref = bests[2]
    print('tmin-tref={0}'.format(-bests[2] + np.min(data_P48_rise['days'])))
    print('tbo=tpeak-tref={0}'.format(np.max(data_P48_rise['days']) - bests[2]))
    print('the characteristic timescale to of the exponent is {0}'.format(bests[1]))
    print('********************************************')
    best_fit = models.model_exponent_concav_pos_withref(bests[0], bests[1], bests[2],
                                             np.linspace(np.min(prior_tref), np.max(my_data[:, 0]),1000)).model_array()

    pylab.figure()
    pylab.plot(best_fit[:,0],best_fit[:,1],'b-')
    pylab.plot(my_data[:, 0], my_data[:, 1],'ro')
    ax=pylab.gca()
    ax.axvspan(np.min(prior_tref),np.max(prior_tref), alpha=0.5, color='yellow')
    #pylab.show()
    #pdb.set_trace()

    #pylab.show()
    #pdb.set_trace()

    ########################################## fit the rising chunk with a power law ###########################################

    a_assumed = 1e-16
    n_assumed = 0.5
    if tref_assumed==None:
        tref_assumed = np.min(data_P48_rise['days'])

    prior_a = np.array(a_prior)
    prior_n = np.array(n_prior)
    # prior_tref=np.array([tref_assumed-20,tref_assumed+1])

    pylab.figure()
    pylab.errorbar(data_P48_rise['days'], data_P48_rise['flux'], yerr=0.5 * data_P48_rise['fluxerr'], fmt='r.')
    pylab.plot(models.model_powerlaw_withref(a=a_assumed, n=n_assumed, xref=tref_assumed,
                                             x=data_P48_rise['days']).model_array()[:, 0],
               models.model_powerlaw_withref(a=a_assumed, n=n_assumed, xref=tref_assumed,
                                             x=data_P48_rise['days']).model_array()[:, 1])
    pylab.xlabel(r"time (JD)")
    pylab.ylabel(r"flux $[erg/s/cm^2/\AA]$")
    # pylab.axvline(JD_peak,color='blue',linestyle='--')
    pylab.title('An attempt to check the prior range')
    pylab.title('P48 R photometry and power law fit')
    pylab.axvline(tref_assumed)
    #pylab.show()
    #pdb.set_trace()

    print("********** FIT WITH A POWER LAW ************")

    output_mcmc = output_mcmc_power

    samples = fitter_general.emcee_n_param(3, models.model_powerlaw_withref, prior_param=[prior_a, prior_n, prior_tref],
                                           data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                           initial_conditions=[a_assumed, n_assumed, tref_assumed],
                                           nwalkers=nwalkers, num_steps=num_steps, flatchain_path=output_mcmc + '/flatchain.txt',
                                           already_run=already_run_mcmc_power)

    best_power = fitter_general.calc_best_fit_n_param(ndim=3, model_nparam=models.model_powerlaw_withref,
                                                      flatchain_path=output_mcmc + '/flatchain.txt',
                                                      data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                                      winners=100,
                                                      output_file_path=output_mcmc,
                                                      bounds=[prior_a, prior_n, prior_tref],
                                                      already_run_calc_all_chis=already_run_calc_all_chis_power,
                                                      show_plots=False)
    # print('reduced chi2 is',best[-1]/(np.shape(my_data)[0]-3))
    # pdb.set_trace()
    if np.shape(my_data)[0] - 3 != 0:
        chi2_power = best_power[-1] / (np.shape(my_data)[0] - 3)
    bests = best_power[:-1]
    if already_plot_power == False:
        fitter_general.plot_opt_fit_n_param(3, models.model_powerlaw_withref, bests, my_data,
                                            flatchain_path=output_mcmc + '/flatchain.txt',
                                            uncertainties=data_P48_rise['fluxerr'], output_file_path=output_mcmc,
                                            xlabel=None, ylabel=None)

        triangle = fitter_general.plot_2D_distributions(
            flatchain_path=output_mcmc + '/flatchain.txt', bests=bests, title='triangle',
            output_file_path=output_mcmc, parameters_labels=['a', 'n', 'tref'])

        histos = fitter_general.plot_1D_marginalized_distribution_mini_interval(
            flatchain_path=output_mcmc + '/flatchain.txt', bests=bests, output_pdf_file_path=output_mcmc,
            output_txt_file_path=output_mcmc,
            parameters_labels=['a', 'n', 'tref'], number_bins=1000)

        #histos = fitter_general.plot_1D_marginalized_distribution(
        #    flatchain_path=output_mcmc + '/flatchain.txt', bests=bests,
        #    output_pdf_file_path=output_mcmc, output_txt_file_path=output_mcmc, parameters_labels=['a', 'n', 'tref'],
        #    number_bins=1000)
    print('********************************************')
    print('when fitting whith the power law, the best fit is {0}, with sigma is {1}'.format(bests[:], np.genfromtxt(
        output_mcmc + '/1sigma_mini.txt')[:, 2]))
    # print('when fitting whith the power law, the best fit is {0}, with sigma is {1}'.format(bests[:-1],np.genfromtxt('tref_calculator_results_powerlaw/1sigma.txt')[:,2])
    best_tref = bests[2]
    print('tmin-tref={0}'.format(-bests[2] + np.min(data_P48_rise['days'])))
    print('tbo=tpeak-tref={0}'.format(np.max(data_P48_rise['days']) - bests[2]))
    print('the index of the power law is {0}'.format(bests[1]))
    print('********************************************')

    best_fit = models.model_powerlaw_withref(bests[0], bests[1], bests[2],
                                             np.linspace(np.min(prior_tref), np.max(my_data[:, 0]),
                                                         1000)).model_array()

    pylab.figure()
    pylab.plot(best_fit[:,0],best_fit[:,1],'b-')
    pylab.plot(my_data[:, 0], my_data[:, 1],'ro')
    ax=pylab.gca()
    ax.axvspan(np.min(prior_tref),np.max(prior_tref), alpha=0.5, color='yellow')
    #pylab.show()
    #pdb.set_trace()

    print('************* Summary of tref calculation ************')
    print('exponential:')
    print('t_ref={0}, chi2_reduced={1}'.format(best_exp[-2], chi2_exp))
    print('power law:')
    print('t_ref={0}, chi2_reduced={1}'.format(best_power[-2], chi2_power))
    if np.shape(my_data)[0] - 3 != 0:
        if chi2_exp < chi2_power:
            print('the best fit is obtained with the exponential, i.e. tref={0}'.format(best_exp[-2]))
            chi2_winner = chi2_exp
            tref_winner = best_exp[-2]
        else:
            print('the best fit is obtained with the power law, i.e. tref={0}'.format(best_power[-2]))
            chi2_winner = chi2_power
            tref_winner = best_power[-2]

        pylab.show()
        return tref_winner, chi2_winner
    else:
        print('there is not enough data to calculate chi2/dof and return a winner')
        pylab.show()
    pylab.show()


'''
WORK IN PROGRESS
def tref_from_P48(path_to_data, tref_assumed=None, band='r',days_rising=10):
    """Description: fit the P48_R light curve with an exponent or a power law to deduce the explosion date
    Input  :- a path to a file where the data is in the following format, 'jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter', e.g. after running
    Read_data.read_data_Marshall_simple() on some Marshall extracted data.
            -the assumed explosion date. If None, it will be taken as the minimum of the data
            -the number of days to fit (i.e. the number of days during which the LC is rising)
            -the band to use. the default is 'r': P48_r. If set to g, will be 'p48_g'
    Output :-
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: samples=fitter_general.emcee_n_param(ndim=2,model_nparam=model_lo_to,prior_param=[prior_alpha,prior_tref],data=my_data,uncertainties=errors,initial_conditions=[alpha_true,tref_true],flatchain_path='output_test_fitter_general/flatchain_test.txt',already_run=False)
    Reliable:  """
    if band=='r':
        if os.path.exists('results_tref_calculator_from_P48R'):
            print('the output directory,results_tref_calculator_from_P48R exists already')
        else:
            os.mkdir('results_tref_calculator_from_P48R')

        output_mcmc_exp = 'results_tref_calculator_from_P48R/exp'
        output_mcmc_power = 'results_tref_calculator_from_P48R/power'
    elif band=='g':
        if os.path.exists('results_tref_calculator_from_P48G'):
            print('the output directory,results_tref_calculator_from_P48G exists already')
        else:
            os.mkdir('results_tref_calculator_from_P48G')

        output_mcmc_exp = 'results_tref_calculator_from_P48G/exp'
        output_mcmc_power = 'results_tref_calculator_from_P48G/power'

    if os.path.exists(output_mcmc_exp) == False:
        os.mkdir(output_mcmc_exp)

    if os.path.exists(output_mcmc_power) == False:
        os.mkdir(output_mcmc_power)

    print('****************************************')
    print('******* Plotting the data to fit *******')
    print('****************************************')

    ########################################## plots the bolometric luminosity and R flux ###########################################

    dict_all = read_data_into_numpy_array(path_to_data, header=True)[
        2]  # ['jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter']

    #print(dict_all)
    condition = (dict_all['mag'] != 99.0) & (dict_all['absmagerr'] != 99.0)

    # remove 99
    for j in dict_all.keys():
        if j not in ['mag', 'absmagerr']:
            #print('the key is', j)
            #print('the dimensions of data_dicts[j] is', np.shape(dict_all[j]))
            #print("the dimension of data_dicts['mag'] is", np.shape(dict_all['mag']))
            #print('the dimension of the boolean is', np.shape([dict_all['mag'] != 99.0]))
            dict_all[j] = dict_all[j][condition]
    dict_all['mag'] = dict_all['mag'][condition]
    dict_all['absmagerr'] = dict_all['absmagerr'][condition]

    data_P48 = dict()

    if band=='r':
        data_P48['days'] = dict_all['jd'][np.asarray(dict_all['filter']) == 'r_p48']
        data_P48['flux'] = dict_all['flux'][np.asarray(dict_all['filter']) == 'r_p48']
        data_P48['fluxerr'] = dict_all['fluxerr'][np.asarray(dict_all['filter']) == 'r_p48']
    elif band=='g':
        data_P48['days'] = dict_all['jd'][np.asarray(dict_all['filter']) == 'g_p48']
        data_P48['flux'] = dict_all['flux'][np.asarray(dict_all['filter']) == 'g_p48']
        data_P48['fluxerr'] = dict_all['fluxerr'][np.asarray(dict_all['filter']) == 'g_p48']

    condition_rise = data_P48['days'] - np.min(data_P48['days']) <= days_rising
    data_P48_rise = dict()
    data_P48_rise['days'] = data_P48['days'][condition_rise]
    data_P48_rise['flux'] = data_P48['flux'][condition_rise]
    data_P48_rise['fluxerr'] = data_P48['fluxerr'][condition_rise]

    pylab.figure()
    pylab.plot(data_P48['days'], data_P48['flux'], 'bo', label='{0}-band data'.format(band),alpha=0.5)
    pylab.plot(data_P48_rise['days'], data_P48_rise['flux'], 'ro', label='data used for the fit',alpha=0.5)
    pylab.axvline(np.max(data_P48_rise['days']))
    pylab.grid()
    pylab.title('{0}-band light curve and points used for fit'.format(band))
    pylab.xlabel(r"time (JD)")
    pylab.ylabel(r"flux $[erg/s/cm^2/\AA]$")
    pylab.tight_layout()
    pylab.legend()
    if band=='r':
        pylab.savefig('results_tref_calculator_from_P48R/lc_and_rising_piece.png')
    elif band=='g':
        pylab.savefig('results_tref_calculator_from_P48R/lc_and_rising_piece.png')



    print('***************************************************')
    print('******* IMPORTANT CHECK BEFORE YOU PROCEED ********')
    print('The red points shown on the plot are the points which will be used for the fit.'
          'Adjust the "days_rising" parameter of Photomanip_fun.tref_from_P48, untill the red points cover all (and do not go beyond) the rising part of the light curve!' )
    print('***************************************************')

    pylab.show()

    print('*************************************************************')
    print('****** CALCULATION OF tref (FIT WITH CONCAVE EXPONENT) ******')
    print('****************  eo*[1-exp(-[t-tref]/to)]   ****************')
    print('*************************************************************')

    ########################################## fit the rising chunk with an exponent (0 for x<xref) ###########################################

    my_data = np.zeros((np.shape(data_P48_rise['days'])[0], 2))
    my_data[:, 0] = data_P48_rise['days']
    my_data[:, 1] = data_P48_rise['flux']

    ########################################## fit the rising chunk with an exponent ###########################################

    eo_assumed = 1e-15
    to_assumed = 15.
    if tref_assumed==None:
        tref_assumed = np.min(data_P48_rise['days'])

    prior_eo = np.array([1e-19, 1e-14])
    prior_to = np.array([0, 50])
    prior_tref = np.array([tref_assumed - 2.5, tref_assumed + 1])


    #print("********** FIT WITH A POSITIVE CONCAVE EXPONENT ************")

    output_mcmc = output_mcmc_exp

    samples = fitter_general.emcee_n_param(3, model_exponent_concav_pos_withref,
                                           prior_param=[prior_eo, prior_to, prior_tref],
                                           data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                           initial_conditions=[eo_assumed, to_assumed, tref_assumed],
                                           nwalkers=100, num_steps=500, flatchain_path=output_mcmc + '/flatchain.txt',
                                           already_run=False)

    best_exp = fitter_general.calc_best_fit_n_param(ndim=3, model_nparam=model_exponent_concav_pos_withref,
                                                    flatchain_path=output_mcmc + '/flatchain.txt',
                                                    data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                                    winners=10,
                                                    output_file_path=output_mcmc,
                                                    bounds=[prior_eo, prior_to, prior_tref],
                                                    already_run_calc_all_chis=False,
                                                    show_plots=False)
    if np.shape(my_data)[0] - 3 !=0:
        chi2_exp = best_exp[-1] / (np.shape(my_data)[0] - 3)
    bests = best_exp[:-1]
    already_plot_exp=False
    if already_plot_exp != True:
        fitter_general.plot_opt_fit_n_param(3, model_exponent_concav_pos_withref, bests, my_data,
                                            flatchain_path=output_mcmc + '/flatchain.txt',
                                            uncertainties=data_P48_rise['fluxerr'], output_file_path=output_mcmc,
                                            xlabel=None, ylabel=None)

        triangle = fitter_general.plot_2D_distributions(
            flatchain_path=output_mcmc + '/flatchain.txt', bests=bests, title='triangle',
            output_file_path=output_mcmc, parameters_labels=['eo', 'to', 'tref'])
        histos = fitter_general.plot_1D_marginalized_distribution(
            flatchain_path=output_mcmc + '/flatchain.txt', bests=bests,
            output_pdf_file_path=output_mcmc, output_txt_file_path=output_mcmc, parameters_labels=['eo', 'to', 'tref'],
            number_bins=1000)
    print('********************************************')
    print('when fitting whith the exponent, the best fit is {0}, with sigma is {1}'.format(bests[:], np.genfromtxt(
        output_mcmc + '/1sigma.txt')[:, 2]))
    # print('when fitting whith the power law, the best fit is {0}, with sigma is {1}'.format(bests[:-1],np.genfromtxt('tref_calculator_results_powerlaw/1sigma.txt')[:,2])
    best_tref = bests[2]
    #print('tmin-tref={0}'.format(-bests[2] + np.min(data_P48_rise['days'])))
    #print('tbo=tpeak-tref={0}'.format(np.max(data_P48_rise['days']) - bests[2]))
    #print('the characteristic timescale to of the exponent is {0}'.format(bests[1]))
    #print('********************************************')
    # pylab.show()
    # pdb.set_trace()

    ########################################## fit the rising chunk with a power law ###########################################

    a_assumed = 1e-16
    n_assumed = 0.5
    if tref_assumed==None:
        tref_assumed = np.min(data_P48_rise['days'])
    prior_a = np.array([1e-18, 1e-15])
    prior_n = np.array([0, 1])
    # prior_tref=np.array([tref_assumed-20,tref_assumed+1])


    print('******************************************************')
    print('****** CALCULATION OF tref (FIT WITH POWER LAW) ******')
    print('*************        a(t-tref)^n        **************')
    print('******************************************************')

    output_mcmc = output_mcmc_power

    samples = fitter_general.emcee_n_param(3, model_powerlaw_withref, prior_param=[prior_a, prior_n, prior_tref],
                                           data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                           initial_conditions=[a_assumed, n_assumed, tref_assumed],
                                           nwalkers=100, num_steps=500, flatchain_path=output_mcmc + '/flatchain.txt',
                                           already_run=False)

    best_power = fitter_general.calc_best_fit_n_param(ndim=3, model_nparam=model_powerlaw_withref,
                                                      flatchain_path=output_mcmc + '/flatchain.txt',
                                                      data=my_data, uncertainties=data_P48_rise['fluxerr'],
                                                      winners=10,
                                                      output_file_path=output_mcmc,
                                                      bounds=[prior_a, prior_n, prior_tref],
                                                      already_run_calc_all_chis=False,
                                                      show_plots=False)
    # print('reduced chi2 is',best[-1]/(np.shape(my_data)[0]-3))
    # pdb.set_trace()
    if np.shape(my_data)[0] - 3 != 0:
        chi2_power = best_power[-1] / (np.shape(my_data)[0] - 3)
    bests = best_power[:-1]
    already_plot_power=False
    if already_plot_power == False:
        fitter_general.plot_opt_fit_n_param(3, model_powerlaw_withref, bests, my_data,
                                            flatchain_path=output_mcmc + '/flatchain.txt',
                                            uncertainties=data_P48_rise['fluxerr'], output_file_path=output_mcmc,
                                            xlabel=None, ylabel=None)

        triangle = fitter_general.plot_2D_distributions(
            flatchain_path=output_mcmc + '/flatchain.txt', bests=bests, title='triangle',
            output_file_path=output_mcmc, parameters_labels=['a', 'n', 'tref'])
        histos = fitter_general.plot_1D_marginalized_distribution(
            flatchain_path=output_mcmc + '/flatchain.txt', bests=bests,
            output_pdf_file_path=output_mcmc, output_txt_file_path=output_mcmc, parameters_labels=['a', 'n', 'tref'],
            number_bins=1000)
    print('********************************************')
    print('when fitting whith the power law, the best fit is {0}, with sigma is {1}'.format(bests[:], np.genfromtxt(
        output_mcmc + '/1sigma.txt')[:, 2]))
    # print('when fitting whith the power law, the best fit is {0}, with sigma is {1}'.format(bests[:-1],np.genfromtxt('tref_calculator_results_powerlaw/1sigma.txt')[:,2])
    best_tref = bests[2]
    #print('tmin-tref={0}'.format(-bests[2] + np.min(data_P48_rise['days'])))
    #print('tbo=tpeak-tref={0}'.format(np.max(data_P48_rise['days']) - bests[2]))
    #print('the index of the power law is {0}'.format(bests[1]))
    #print('********************************************')
    #pylab.show()
    # pdb.set_trace()

    print('**********************************************************')
    print('****** SUMMARY OF tref CALCULATION WITH BOTH MODELS ******')
    print('**********************************************************')

    print('exponential:')
    print('t_ref={0}, chi2_reduced={1}'.format(best_exp[-2], chi2_exp))
    print('power law:')
    print('t_ref={0}, chi2_reduced={1}'.format(best_power[-2], chi2_power))
    if np.shape(my_data)[0] - 3 != 0:
        if chi2_exp < chi2_power:
            print('the best fit is obtained with the exponential, i.e. tref={0}'.format(best_exp[-2]))
            chi2_winner = chi2_exp
            tref_winner = best_exp[-2]
        else:
            print('the best fit is obtained with the power law, i.e. tref={0}'.format(best_power[-2]))
            chi2_winner = chi2_power
            tref_winner = best_power[-2]
        pylab.show()
        return tref_winner, chi2_winner
    else:
        print('there is not enough data to calculate chi2/dof and return a winner')
        pylab.show()
'''

