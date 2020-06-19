# ==========
#
# Library of functions that aim to aid access and use to elegant, genesis1.3, and ICS code. 
# Includes particle distribution manipulation functions. 
# Able to work with SDDS IO. 
# Deals with the difference in particle distribution definitions between elegant and genesis.
# 
# ==========

import sdds # only python 3.6 is supported for SDDS (2020.04.15)
import os
import subprocess
import csv
import re
import numpy as np
import scipy.constants as pc
import pandas as pd

# path of the module 
file_path = os.path.dirname(os.path.abspath(__file__))

#
# === === === helper functions === === ===
#

def gaussian(x, mu, sig):
    """ Gaussian
    x - ndarray
    mu - float. mean
    sig - float. std dev
    
    """
    y = (1 / (np.sqrt(2*pc.pi)  * sig)) * np.exp( - (x - mu)**2 / (2 * sig**2))
    
    return y

def gaussianDx(x, mu, sig):
    """ First derivative of a Gaussian. Normalized so that its maximum is 1.
    x - ndarray
    mu - float. mean
    sig - float. std dev
    
    """
    
    y = (-(x - mu)/sig**2) * gaussian(x, mu, sig)
    # NOTE: normalization to maximum value = 1
    y = y/np.abs(y).max()
    
    return y

def gaussianDx2(x, mu, sig):
    """ Second derivative of a Gaussian (inverted sombrero). Normalized so that its maximum is 1.
    x - ndarray
    mu - float. mean
    sig - float. std dev
    
    """
    
    y = ( -1/sig**2 + ((x - mu)/sig**2)**2 ) * gaussian(x, mu, sig)
    # NOTE: normalization to maximum value = 1
    y = y/np.abs(y).max()
    return y



#
# === === === genesis === === ===
#

class GenControl:
    """ This is a class that controls and interacts with genesis.

    genControl is intended to provide python features with which to interact with genesis. The intent is to eliminate the need to manually run genesis and gather its output. 
 
    """
        
    def __init__(self, index):

        self.FILENAME_IN = ''
        self.FILENAME_OUT = ''
        self.FILENAME_DF = ''
        self.FILENAME_PAR = ''

        self.PARAMS_OK = []
        try:
            with open(file_path + r'\gen_accepted_param_names.txt', newline='') as fileok:
                okreader = csv.reader(fileok, delimiter=' ', skipinitialspace=True)
                for row in okreader:
                    self.PARAMS_OK = self.PARAMS_OK + [row[0].upper().strip(', ')]
        except FileNotFoundError:
            print('file `gen_accepted_param_names.txt` is not present in the `../src/features` directory. This will effect error checking parameter names.')

    def make_infile(self, inputfn):
        """ Make a genesis input file by running geneis with default settings, extracting the input parameters from the out file and putting them into a text file with the provided name.
        
        Specific parametrs can be edited with paramSet.
        """
        # print(inputfn)

        # === Run `genesis` with its default parameters. ===
        print('Running `genesis` with its default parameters.')
        
        # list of commands executed outside python script
        cmdlist = [
            'echo $newrun > temp.in',
            "echo outputfile = 'temp.out' >> temp.in",
            'echo $end >> temp.in',
            'echo temp.in > runme',
            'echo n >> runme',
            'genesis < runme'
            ]

        for cmdstr in cmdlist:
            subprocess.run(cmdstr, shell=True)

        ### extract the input file parameters from the temp.out file
        print('Create an input file from the default `genesis` parameters.')

        with open('temp.out', 'r') as fout:
            outtext = fout.readlines()
        # find the line numbers of the first and last parameter lines
        for i in range(len(outtext)):
            if '$newrun' in outtext[i]:
                rstart = i
            elif '$end' in outtext[i]:
                rend = i+1
                break
            elif ('=' not in outtext[i]) and ('$' not in outtext[i]):
                # lout line in the .out file runs to the next line for some reason. 
                # why genesis? 
                outtext[i-1] = outtext[i-1][:-1] + outtext.pop(i)
        # filter out the 'itram##' lines for readability
        inlist = list(filter(lambda x: ('itram' not in x), outtext[rstart:rend]))
        # write to file given by inputfn
        fh = open(inputfn, 'w') 
        for row in inlist:
            fh.write(row)
        fh.close()

        # cleaning up workspace
        os.remove('runme')
        os.remove('temp.in')
        os.remove('temp.out')

        print('New input file was created\nfile name = '+inputfn)
        print('this was the `makeinfile` function. \nbye.')
        return 0

    def param_set(self, inputfn, pnv, verbose=False):
        """ Set the parameter with name 'pname' to the value 'pvalue', inside the input filename 'inputfn'. Makes the changes inplace.
        `pnv` is a dictionary of names and values for parameters:
        pnv = dict([('pname0',pvalue0), ('pname1',pvalue1), ...])
        pnv = {'pname0':pvalue0, 'pname0':pvalue0, ...}
        'pname' is a str (should not be case sensitive)
        type of 'pvalue' depends on the type that is expected by the parameter as per the genesis input file. 

        This feature is best suited for setting new values of parameters on the fly and executing genesis with them. For example, you would like to run a script that scans multiple genesis parameters.
        If you are setting up a new genesis simulation after creating a template input file with `makeinfile()`, it could be more straight forward to open the .in file in a text editor and make your changes there. 
        """

        # Check validity of parameter names in pnv variable
        for pn in pnv.keys():
            if pn.upper() not in self.PARAMS_OK:
                print('WARNING:\nThe parameter `' + pn + '` is not found in the list of accepted genesis input parameters. ')
        
        # load the current version of the in file
        with open(inputfn, 'r') as fout:
            oldin_text = fout.readlines()
        
        # regex a pattern for the lines of the in file
        pattern = re.compile(r'(\w*)\s*=\s*(.*)',re.IGNORECASE)
        # make a dict from in file
        in_dict = {}
        for ii in oldin_text:
            s = pattern.search(str(ii))
            try:
                # regex group(2) will have a different interpretation based on the parameter it defines. int, float, string, list of ints/floats.
                in_dict[str(s.group(1)).upper()] = s.group(2)
            except AttributeError:
                # catch exceptions to the regex pattern search+grouping
                # this should happen only for the $newrun and $end lines of a properly written in file.
                pass
        
        for pn in pnv.keys():
            # delete parameter if 'delete' and key is in dict 
            # (i.e. don't throw an error if parameter is already missing)
            if (str(pnv[pn]).lower() == 'delete'):
                if (pn.upper() in in_dict):
                    del in_dict[pn.upper()]
                    if verbose:
                        print(pn.upper() + ' has been DELETED from input file.')
            else:
                # update the parameter value inplace
                in_dict[pn.upper()] = pnv[pn]

        # write to file given by inputfn
        try:
            fh = open(inputfn, 'w')
            fh.write('$newrun\n')
            for pn in in_dict.keys():
                fh.write(' ' + pn.upper() + ' = ' + str(in_dict[pn]) + '\n' )
            fh.write('$end\n')
            fh.close()
            if verbose:
                print('The input file : ' + inputfn + ' has been updated by the `.param_set()` method.')
            
        except:
            print('`.param_set()` method has failed to write to `' + inputfn + '`.')

        return 0


    
    def param_get(self, inputfn, pnl='all'):
        """ Get and return the value of the requested parameter inside the file given by filename.

        inputfn - (str) path to the input file from which to get parameter value
        pnl - (str) parameter names e.g. 'GAMMA0', e.g. ['gamma0', 'AW0', 'lout']
        
        This function loads the whole input file into memory and makes a dict out of the parameter names and parameter values. This could possibly be inefficient if a bunch of calls are made in succession. However, instead of that, the user should use a list of parameter names instead of calling the function multiple times.
        """

        # catch single string object instead of list of strings
        if type(pnl) is not list: pnl = [ pnl ]

        # Check validity of parameter names in pnl variable
        for pn in pnl:
            if (pn.upper() not in self.PARAMS_OK) and (pn != 'all'):
                print('WARNING:\nThe parameter `' + pn + '` is not found in the list of accepted genesis input parameters. ')
        
        # load the current version of the in file
        with open(inputfn, 'r') as fout:
            oldin_text = fout.readlines()
        
        # regex a pattern for the lines of the .in file
        pattern = re.compile(r'(\w*)\s*=\s*(.*)',re.IGNORECASE)
        # make a dict from in file
        in_dict = {}
        for ii in oldin_text:
            s = pattern.search(str(ii))
            try:
                # regex group(2) will have a different interpretation based on the parameter it defines. int, float, string, list of ints/floats.
                in_dict[str(s.group(1)).upper()] = s.group(2)
            except AttributeError:
                # catch exceptions to the regex pattern search+grouping
                # this should happen only for the $newrun and $end lines of a properly written in file.
                pass
        
        # if any of the parameter names are 'all' then return all parameters
        if 'all' in pnl:
            return list(in_dict.values())
        else:
            pval = []
            keyerrmsg = in_dict[pn.upper()] + ' was not found as a parameter in the input file: ' + inputfn
            for pn in pnl:
                try:
                    pval = pval + [ in_dict[pn.upper()] ]
                except KeyError as keyerrmsg:
                    print(keyerrmsg)
                    break
                
            return pval

    def run_genesis(self, inputfn):
        """ Runs genesis with the specified input file.
        inputfn - str. name of the input file to genesis. No path, just the file name.
        It is assummed that .run_genesis() is executed inside a directory that contains the file.
        This function effectively executes the following shell commands:
            $ echo inputfn > runme
            $ echo y >> runme
            $ genesis<runme
        """
        # === Run `genesis` with its default parameters. ===
        print('Run `genesis`.')
        
        # list of commands executed outside python script
        cmdlist = [
            'echo ' + inputfn + ' > runme',
            'echo y >> runme',
            'genesis < runme'
            ]

        for cmdstr in cmdlist:
            subprocess.run(cmdstr, shell=True)
        
        # Note: It would be nice to indicate when genesis is done running. Currently, the console output is the only indication of running.


class GenOut():
    """Class handles the import of genesis output files.
    Extracts and makes available relevant info like, number of z records, parameter names and values, etc.
    Requires sdds.
    """
    
    def __init__(self):
        self.filename = 'no filename'
        self.ncol = 0
        self.nslice = 0
        self.nzrec = 0
        self.data = pd.DataFrame([])
        self.colnames = np.array([])
        # self.parnames = np.array([])
        # self.parvalues = np.array([])
        self.paramsin = {}
    
    def load(self, gen_out_filename, verbose=False):
        """ Load the genout sdds file, make data array, get parameters."""
        self.filename = gen_out_filename

        genout = sdds.SDDS(0)
        genout.load(self.filename)
        
        # use these column names to slice the np array
        self.colnames = [cn.lower() for cn in np.array(genout.columnName)]
                
        # get size
        [self.ncol, self.nslice, self.nzrec] = np.shape(genout.columnData)

        # convert to pd DataFrame
        if self.nslice == 1:
            self.data = pd.DataFrame(np.array(genout.columnData)[:, 0, :].T, columns=self.colnames)
            # add a column for the slice number to handle the case of multiple slices
            self.data['slice'] = np.zeros(self.nzrec)

        elif self.nslice > 1:
            for slicenumber in range(self.nslice):
                # make a dataframe for each slice
                datatemp = pd.DataFrame(np.array(genout.columnData)[:, slicenumber, :].T, columns=self.colnames)
                # add a column for the slice number to handle the case of multiple slices
                datatemp['slice'] = slicenumber * np.ones(self.nzrec)
                # append slice dataframe to self.data
                self.data = self.data.append(datatemp)
        
        self.colnames = self.colnames + ['slice']
         
        parnames = np.array(genout.parameterName)
        parvalues = np.array(genout.parameterData)
        self.paramsin = {parnames[ii]:parvalues[ii] for ii in range(len(parnames))}
        
        if verbose:
            print('==========')
            print('output filename = '+self.filename)
            print('The genesis output file contains:'+
                '\n{} columns\n{} pages\n{} z-records'.format(self.ncol, self.nslice, self.nzrec))
            print('==========')
        
        return 1

    def calc_gain(self, ss=0, sig=1.0, nsig=6):
        
        """ Attempts to estimate the gain-length for the specified slice.

        Convolve with d2(gaussian)/dx2 to simultaneously smooth and differentiate log(power). The output response will have a maximum at the beginning of the exponential gain and a minimum when saturation takes over and the gain diminishes. 
        Make sure that the first and last points of response are dropped from the max/min search to avoid boundary effects from the filtering. This can be achieved by using mode='valid'
        Mode ‘valid’ returns output of length max(M, N) - min(M, N) + 1.
        
        ss - int. the indexing number of the genesis slice. if there is only one slice, nslice=0 
        sig=1.0, float. The width of the gaussianDx2 filter.
        nsig=6, int. The width of the interval over which to compute the filter array.

        """

        # the log(P) signal
        xvec = self.data[self.data.slice == ss].z.to_numpy()
        yvec = np.log(self.data[self.data.slice == ss].power.to_numpy())

        delz = np.abs(xvec[1] - xvec[0])
        mu = 0
        
        xvecfilter = np.arange(-nsig*sig, nsig*sig, delz)
        yvecfilter = gaussianDx2(xvecfilter, mu, sig)
        nptsfilter = yvecfilter.shape[0]

        # convolve
        dy2vec = np.convolve(yvecfilter, yvec, mode='valid')
        
        # find min and max indecies and account for mode='valid'
        
        imax = dy2vec.argmax()
        imin = dy2vec.argmin()
        if imax.shape != ():
            imax = imax[0]
        if imin.shape != ():
            imin = imin[0]
        
        imax = imax + int(nptsfilter/2)
        imin = imin + int(nptsfilter/2)

        if imax > imin:
            print('Gain-length calculation failed. The max and min of the second derivative of log(power) are out of order. Power gain seems to be ill defined. ')
            print('=====')
            return 0
        
        # slope of log(P) = 1/Lgain
        dyvec = (yvec[1:] - yvec[:-1]) / delz
        gain = dyvec[imax:imin].mean()
        lgain = 1 / gain

        return nptsfilter, imax, imin, lgain


class GenPar():
    """Class handles the import of genesis .par files.
    Extracts and makes available relevant info like, number of particles, etc.
    Requires sdds.
    The FEL parameters can be found in the genout() class above. 
    """
    
    def __init__(self):
        self.filename = 'no filename'
        self.ncol = 0
        self.nzrec = 0
        self.nparticles = 0
        self.data = pd.DataFrame([])
        self.colnames = []
    
    def load(self, gen_par_filename):
        """ load the .par sdds file, make data array, get parameters."""
        self.filename = gen_par_filename

        genpar = sdds.SDDS(0)
        genpar.load(self.filename)

        # use these column names to slice the np array
        self.colnames = [cn.lower() for cn in np.array(genpar.columnName)]
        # get size
        [self.ncol, self.nzrec, self.nparticles] = np.shape(genpar.columnData)

        # make into DataFrame and add a column for the zrecord
        for zz in range(self.nzrec):
            zzpd = pd.DataFrame(np.array(genpar.columnData)[:,zz,:].T, columns=self.colnames)
            zzpd['zrec'] = np.zeros(self.nparticles) + zz
            # append to big DF object
            self.data = self.data.append(zzpd)
        # add the zrec column 
        self.colnames = self.colnames + ['zrec']
        
        print('==========')
        print('Loaded .par filename = '+self.filename)
        print('The genesis .par file contains:'+
              '\n{} columns\n{} z-records\n{} particles'.format(self.ncol, self.nzrec, self.nparticles))
        print('Column names:\n'+str(self.colnames))
        print('==========')

#
# === === === elegant === === ===
#

class EleControl():

    def __init__(self):
        self.FILENAME_ELE = 'no filename'

    def run_elegant(self, elepath, inputfn, verbose=False):
        """ Runs elegant with the specified input file.
        inputfn - str. 
        This function effectively executes the following shell commands:
            $ elegant `inputfn`
        """
        # remember where the call is coming from
        origincw = os.getcwd()
        # step into the folder that contains the inital particle distributions
        os.chdir(elepath)

        # === Run `elegant`. ===
        if verbose:
            print('Running `elegant`.')

        # list of commands executed outside python script
        cmdlist = [
            'elegant ' + inputfn
            ]

        for cmdstr in cmdlist:
            if verbose:
                print(cmdstr)
                print('===')
            subprocess.run(cmdstr, shell=True)

        # get back to origin CW
        os.chdir(origincw)

# 
# === === === Particle Distributions === === ===
#

class ParticleDist():

    """ Class with functionality to create a custom 6D phase-space for particle distributions compatible with elegant's SDDS format and/or genesis' text based distfile.

    """

    def __init__(self):
        self.ebeamsdds = sdds.SDDS(0)
        self.ebeamnp = np.array([])
        self.colnames = ['x', 'xp', 'y', 'yp', 't', 'p', 'particleID']
        self.collookup = {self.colnames[ii]:ii for ii in range(len(self.colnames))}
        self.data = pd.DataFrame([])

    def load_sdds(self, sddsfilename, makedf=True):
        """ loads a SDDS file class instance for custom 6d phas-space manipulation.
        sddsfilename - str. path of the sdds file being loaded. Can be an absolute path or just the filename if the current working directory contains the file.

        The default column names output by elegant to a SDDS file are:
        sdds.SDDS(0).columnName = ['x', 'xp', 'y', 'yp', 't', 'p', 'particleID']

        elegant: (x,xp,y,yp,t,p), where x and y are in meters, xp = x′ and xp = y′ are dimensionless, t is in seconds, and p = gamma*beta is the dimensionless momentum. 
        """
        try:
            # load the beam SDDS file
            self.ebeamsdds.load(sddsfilename)
            # convert to numpy array
            self.ebeamnp = np.array(self.ebeamsdds.columnData)
            if makedf:
                self.data = pd.DataFrame(self.ebeamnp, columns=self.colnames)
        except FileNotFoundError:
            print(sddsfilename + ' could not be found. Please check that the file is present in the specified directory.')

    def xE_parabola(self, ascale=0.12, momres=138.5, momoffset=0.0, absx=False,zeropx=True):
        """ Creates a parabolic correlation between the x position and relative energy of the ebeam macro-particles. 
        
        ascale - float. parabola scale
        momres - float. the resonant momentum inside the FEL
        momoffset - float. offset of the parabola. 
        note: the output 'p' array is always made such that the minimum energy particle has a resonant momentum given by 'momres'. This is true even if the input 'p' array has a minimum momentum different from 'momres'. The variable 'momoffset' can used to shift this output minimum. 
        
        zeropx - bool. zero out the transverse momentum in the x direction for the particles. 
        

        """
        
        if ascale >= 0:
            # ensure that the min energy particle is on resonance
            # this shift ensures that there is no negative eta under the sqrt
            colp = self.collookup['p']
            self.ebeamnp[colp, 0, :] = (self.ebeamnp[colp, 0, :] - self.ebeamnp[colp, 0, :].min() ) + momres

            # define the position and relative energy vectors
            pos = self.ebeamnp[self.collookup['x'], 0, :]
            eta = (self.ebeamnp[self.collookup['p'], 0, :] + momoffset - momres) / momres

            # make correlation only if ascale is not zero
            if ascale > 0:
                pos = 1e-6 * ascale * np.sqrt( eta )

            # update the position column of the file
            if absx:
                # put all particles on positive x values even if not conditioned
                self.ebeamnp[self.collookup['x'], 0, :] = np.abs(pos)
            else:
                self.ebeamnp[self.collookup['x'], 0, :] = pos
            
            # change the energy to match the x-position in the BC parabola
            self.ebeamnp[self.collookup['p'], 0, :] = eta * momres + momres
        
        else:
            pass
                
        # zero px. make x momentum zero for matching into lattice?
        if zeropx:
            self.ebeamnp[self.collookup['xp'], 0, :] = np.zeros(self.ebeamnp[0,0,:].shape)

        # update the sdds holder for the 6D phase-space
        for col in self.colnames:
            self.ebeamsdds.setColumnValueLists(col, [self.ebeamnp[self.collookup[col], 0, :].tolist()] )

    def xE_condition(self, ascale=0.12, momres=138.5, momoffset=0.0,absx=False,zeropx=True):
        """ Creates a correlation between the x position and relative energy of the ebeam macro-particles. The formula is derived from equating the de-phasing (slippage) due to energy spread to the added path-length from the betatron oscillations in the strong FODO lattice. The equation is:

        ascale - float. parabola scale
        momres - float. the resonant momentum inside the FEL
        momoffset - float. offset of the parabola. 
        note: the output 'p' array is always made such that the minimum energy particle has a resonant momentum given by 'momres'. This is true even if the input 'p' array has a minimum momentum different from 'momres'. The variable 'momoffset' can used to shift this output minimum. 
        
        zeropx - bool. zero out the transverse momentum in the x direction for the particles. This does not perserve emittance.
        

        """
        if ascale >= 0:
            # ensure that the min energy particle is on resonance
            # this shift ensures that there is no negative eta under the sqrt
            colp = self.collookup['p']
            self.ebeamnp[colp, 0, :] = (self.ebeamnp[colp, 0, :] - self.ebeamnp[colp, 0, :].min() ) + momres

            # define the position and relative energy vectors
            pos = self.ebeamnp[self.collookup['x'], 0, :]
            mom1 = self.ebeamnp[self.collookup['p'], 0, :] + momoffset

            # make correlation only if ascale is not zero
            if ascale > 0:
                
                pos = 1e-6 * ascale * np.sqrt( (mom1**2 - momres**2) / (momres**2 * mom1) )

            # update the position column of the file
            if absx:
                # put all particles on positive x values even if not conditioned
                self.ebeamnp[self.collookup['x'], 0, :] = np.abs(pos)
            else:
                self.ebeamnp[self.collookup['x'], 0, :] = pos
            
            # change the energy to match the x-position in the BC parabola
            self.ebeamnp[self.collookup['p'], 0, :] = mom1
        
        else:
            if absx:
                # put all particles on positive x values even if not conditioned
                self.ebeamnp[self.collookup['x'], 0, :] = np.abs(pos)
                
        # zero px. make x momentum zero for matching into lattice?
        if zeropx:
            self.ebeamnp[self.collookup['xp'], 0, :] = np.zeros(self.ebeamnp[0,0,:].shape)

        # update the sdds holder for the 6D phase-space
        for col in self.colnames:
            self.ebeamsdds.setColumnValueLists(col, [self.ebeamnp[self.collookup[col], 0, :].tolist()] )


    def save_sdds(self, outputfilename, genpscale=True):
        """ Saves the current self.ebeamsdds to SDDS file.
        Apllies energy scaling if the distribution will be used for genesis.
        The difference is in the transverse momentum definitions. Elegant uses x' and y' which are just angles (i.e. x' = beta_x / beta). Genesis uses the transverse momentum in units of mc (i.e. px = gamma*beta_x = p*x', where p is the elegant momentum).
        """
        if genpscale:
            # print('inside genpscale -----')
            self.ebeamnp[self.collookup['xp'],0,:] = (
                self.ebeamnp[self.collookup['p'],0,:] 
                * self.ebeamnp[self.collookup['xp'],0,:] 
                )
            self.ebeamnp[self.collookup['yp'],0,:] = (
                self.ebeamnp[self.collookup['p'],0,:] 
                * self.ebeamnp[self.collookup['yp'],0,:]
                )
        # update the sdds holder for the 6D phase-space
        for col in self.colnames:
            self.ebeamsdds.setColumnValueLists(col, [self.ebeamnp[self.collookup[col], 0, :].tolist()] )

        # save to outputfilename
        self.ebeamsdds.save(outputfilename)

    
    def create_genesis_distfile(self, sddsfilename, gendistfilename,charge=1200e-12, verbose=False):
        """ Create a distfile for genesis input from the SDDS file given in sddsfilename.
        """
        # helper sddsprintout filename
        sddsprintoutFilename = 'sdds_prinout.txt'
        if os.path.isfile(sddsprintoutFilename):
            pass
        else:
            print('=== === ===')
            print('~~~!~!~!~~~')
            print('Please make a temporary empty file named `sdds_printout.txt` inside the genesis output directory.\n Sorry for inconvenience. :[ ')
            print('~~~!~!~!~~~')
            print('=== === ===')
        # these column names apply to a SDDS file that comes out of elegant
        colnamesprintout = '(p,t,x,y,xp,yp)'
        
        # sddsprintout command call outside py script
        cmdsddsprintout = 'sddsprintout ' + sddsfilename + ' ' + sddsprintoutFilename + ' -columns=' + colnamesprintout
        print(cmdsddsprintout)
        subprocess.run(cmdsddsprintout, shell=True)

        # open the file from which to read data
        fh_sddsprintout = open(sddsprintoutFilename,'r')

        line_count = 0
        with open(gendistfilename, mode='w') as genDist:
            for row in fh_sddsprintout:
                if line_count == 0:
                    genDist.write(' ? VERSION = 1.3\n')
                    line_count += 1
                elif line_count == 1:
                    # change this for a different beam-charge
                    genDist.write(' ? CHARGE = ' + str(charge) + '\n')
                    line_count += 1
                elif line_count == 2:
                    # make these column names match the genesis expectations
                    genDist.write(' ? COLUMNS P T X Y PX PY\n')
                    line_count += 1
                elif line_count >= 5:
                    genDist.write(row) #don't need a '\n' newline here
                    line_count += 1
                else:
                    line_count += 1
        # close the read file
        fh_sddsprintout.close()

        if verbose:
            print('=====')
            print('genesis 6D filename\n>>> ' + gendistfilename)
            print('\n\n==========')


    def create_ics_distfile(self, sddsfilename, icsdistfilename, verbose=False):
        """ Create a distfile for genesis input from the SDDS file given in sddsfilename.

        The LBNL ICS code takes in a particle file that is has the rows: x (cm), dx/dz (mrad), y (cm), dy/dz (mrad), phase (deg), energy (MeV)

        phase is the RF phase of the accelerator and it is used to specify the time of the particles as a position inside the RF wave. 
        For a S-band accelerator we have:
        f_RF = 2.8559 GHz
        wavelength_RF = 10.47 cm

        Transformation from a particle's z coordinate to phi:
        phi = Mod( z / wavelength_RF - pi, 2pi ) - pi

        """
        # helper sddsprintout filename
        sddsprintoutFilename = 'sdds_prinout.txt'
        if os.path.isfile(sddsprintoutFilename):
            pass
        else:
            print('=== === ===')
            print('~~~!~!~!~~~')
            print('Please make a temporary empty file named `sdds_printout.txt` inside the genesis output directory.\n Sorry for inconvenience. :[ ')
            print('~~~!~!~!~~~')
            print('=== === ===')
        # these column names apply to a SDDS file that comes out of elegant
        colnamesprintout = '(p,t,x,y,xp,yp)'
        
        # sddsprintout command call outside py script
        cmdsddsprintout = 'sddsprintout ' + sddsfilename + ' ' + sddsprintoutFilename + ' -columns=' + colnamesprintout
        print(cmdsddsprintout)
        subprocess.run(cmdsddsprintout, shell=True)

        # open the file from which to read data
        fh_sddsprintout = open(sddsprintoutFilename,'r')


#
# === === === ICS LLNL code === === ===
#

class ICSControl():

    def __init__(self):
        self.FILENAME_ICS = ''
        self.PARAMS_OK = []
        try:
            with open(file_path + r'\ics_accepted_param_names.txt', newline='') as fileok:
                okreader = csv.reader(fileok, delimiter=' ', skipinitialspace=True)
                for row in okreader:
                    self.PARAMS_OK = self.PARAMS_OK + [row[0].upper().strip(', ')]
        except FileNotFoundError:
            print('file `ics_accepted_param_names.txt` is not present in the `../src/features` directory. This will effect error checking parameter names.')

    def run_ics(self, icspath, verbose=False):
        """ Runs an instance of the ICS code from LLNL. 
        There is an executable inside the folder `Compton_Code` named `Compton.exe`. The input file to that .exe must be in the same directory and named `Compton.ini`. 

        icspath - str. directory path to the folder in which the relevant .ini file is located. The directory must also include the .exe file that runs the code.

        This function effectively executes the following shell commands:
            $ Compton.exe
        """
        # remember where the call is coming from
        origincw = os.getcwd()
        # step into the folder that contains the inital particle distributions
        os.chdir(icspath)

        # === Run `Compton.exe`. ===
        if verbose:
            print('Running `Compton.exe`.')

        # list of commands executed outside python script
        cmdlist = [
            'Compton.exe'
                    ]

        for cmdstr in cmdlist:
            if verbose:
                print(cmdstr)
                print('===')
            subprocess.run(cmdstr, shell=True)

        # get back to original CW
        os.chdir(origincw)
        
        # NOTE: It would be nice to indicate when genesis is done running. Currently, the console output is the only indication of running.

        return 1

    def param_set(self, pnv, inputfn='Compton.ini', verbose=False):
        """ Set the parameter with name 'pname' to the value 'pvalue', inside the input filename 'inputfn'. Makes the changes inplace.
        
        inputfn - str. 
            The LLNL ICS code expects an input file that is named 'Compton.ini', so that is the default inputfn. 
        `pnv` - dictionary of names and values for parameters:
            pnv = dict([('pname0',pvalue0), ('pname1',pvalue1), ...])
            pnv = {'pname0':pvalue0, 'pname0':pvalue0, ...}
            'pname' is a str (should not be case sensitive)
            'pvalue' type depends on the type that is expected by the parameter as per the LLNL input file.
                'pvalue' can be set to 'delete' in order to remove the parameter

        This feature is best suited for setting new values of parameters on the fly and executing genesis with them. For example, you would like to run a script that scans multiple genesis parameters.
        If you are setting up a new genesis simulation after creating a template input file with `makeinfile()`, it could be more straight forward to open the .in file in a text editor and make your changes there. 
        """

        # Check validity of parameter names in pnv variable
        for pn in pnv.keys():
            if pn.upper() not in self.PARAMS_OK:
                print('WARNING:\nThe parameter `' + pn + '` is not found in the list of accepted LLNL input parameters. ')
        
        # load the current version of the .ini file
        with open(inputfn, 'r') as fout:
            oldin_text = fout.readlines()
        
        # regex a pattern for the lines of the in file
        pattern = re.compile(r'(\w*)\s*(.*)',re.IGNORECASE)
        # make a dict from .ini file
        in_dict = {}
        for ii in oldin_text:
            s = pattern.search(str(ii))
            try:
                # do not add the exiting string to the parameter dictionary
                if str(s.group(1)) == 'EX':
                    pass
                else:
                    # regex group(2) will have a different interpretation based on the parameter it defines. int, float, string, list of ints/floats.
                    in_dict[str(s.group(1)).upper()] = s.group(2)

            except AttributeError:
                # catch exceptions to the regex pattern search+grouping
                # should happen for commented out lines
                print(s)
                pass
        
        # update/delete the parameter value inside the dictionary
        for pn in pnv.keys():
            # delete parameter if 'delete' and key is in dict 
            # (i.e. don't throw an error if parameter is already missing)
            if (str(pnv[pn]).lower() == 'delete'):
                if (pn.upper() in in_dict):
                    del in_dict[pn.upper()]
                    if verbose:
                        print(pn.upper() + ' has been DELETED from input file.')
            else:
                # update the parameter value inplace based on type
                if (type(pnv[pn]) is list) or (type(pnv[pn]) is np.ndarray):
                    # format the list into a string with spaces
                    liststring = str(pnv[pn]).strip('[ ]').replace(',',' ')
                    in_dict[pn.upper()] = liststring
                
                else:
                    # just write the value
                    in_dict[pn.upper()] = pnv[pn]
                

        # write to file given by inputfn
        try:
            fh = open(inputfn, 'w')
            
            for pn in in_dict.keys():
                fh.write(pn.upper() + ' ' + str(in_dict[pn]) + '\n' )
            fh.write('EX\n')
            fh.close()
            if verbose:
                print('The input file : ' + inputfn + ' has been updated by the `.param_set()` method.')
            
        except:
            print('`.param_set()` method has failed to write to `' + inputfn + '`.')

        return 1
        
    def param_get(self, inputfn='Compton.ini', pnl='all'):
        """ Get and return the value of the requested parameter inside the file given by filename.

        inputfn - (str) path to the input file from which to get parameter value
        pnl - (str) parameter names e.g. 'WAVELENGTH', e.g. ['wavelength', 'gen_part_basic', 'LBASE']
        
        This function loads the whole input file into memory and makes a dict out of the parameter names and parameter values. This could possibly be inefficient if a bunch of calls are made in succession. However, instead of that, the user should use a list of parameter names instead of calling the function multiple times.
        """

        # catch single string object instead of list of strings
        if type(pnl) is not list: pnl = [ pnl ]

        # Check validity of parameter names in pnl variable
        for pn in pnl:
            if (pn.upper() not in self.PARAMS_OK) and (pn != 'all'):
                print('WARNING:\nThe parameter `' + pn + '` is not found in the list of accepted genesis input parameters. ')
        
        # load the current version of the .ini file
        with open(inputfn, 'r') as fout:
            oldin_text = fout.readlines()
        
        # regex a pattern for the lines of the .in file
        pattern = re.compile(r'(\w*)\s*(.*)',re.IGNORECASE)
        # make a dict from in file
        in_dict = {}
        for ii in oldin_text:
            s = pattern.search(str(ii))
            try:
                # regex group(2) will have a different interpretation based on the parameter it defines. int, float, string, list of ints/floats.
                in_dict[str(s.group(1)).upper()] = s.group(2)
            except AttributeError:
                # catch exceptions to the regex pattern search+grouping
                # this should happen only for the $newrun and $end lines of a properly written in file.
                pass
        
        # if any of the parameter names are 'all' then return all parameters
        if 'all' in pnl:
            return list(in_dict.values())
        else:
            pval = []
            keyerrmsg = in_dict[pn.upper()] + ' was not found as a parameter in the input file: ' + inputfn
            for pn in pnl:
                try:
                    pval = pval + [ in_dict[pn.upper()] ]
                except KeyError as keyerrmsg:
                    print(keyerrmsg)
                    break
                
            return pval

    
    # TODO: add make ini file? This may not be that useful...

    
class ICSOut(object):

    def __init__(self, filename=None):
        """ init class """

        if filename:
            self.filename = filename
        else:
            self.filename = ''

        self.simdata = pd.DataFrame([])
        self.ntheta = 0
        self.nenergy = 0
        self.ecentral = 0 # xray energy with maximum flux for angles close to axis 0.1*max(theta)
        self.photons_per_ev = np.array([]) # xray spectrum (inegrated flux over angles)
        self.photons_per_mrad = np.array([]) # xray angular distribution (integrated flux over energy)

    def load_data(self, filename=None, mode='photon_ang_spect', verbose=False):
        """ Loads the specified file name as a pandas DF and updats;
        .ntheta, .nenergy 
        
        TODO:
         - make loading feature for different code calculations:
            - angular flux
            - total flux 
            - energy spectrum
            - ...

        """
        
        if mode == 'photon_ang_spect':
            
            if (filename == None) & (self.filename == ''):
                outputfn = 'Graph_Angle_Spectrum2.txt'
            elif (filename == None) & (self.filename != ''):
                outputfn = self.filename
            else:
                outputfn = filename
                self.filename = outputfn

            colnames = ['theta', 'energy', 'flux']
            self.simdata = pd.read_csv(outputfn, 
                    delimiter='\s+', names=colnames, skiprows=1, dtype=np.float32)

            # get the header line and extract the number of points in the theta and energy mesh
            with open(outputfn,'r') as f:
                headerline = f.readline().split()
            
            if verbose:
                print(headerline)

            self.ntheta = np.float(headerline[-2])
            self.nenergy = np.float(headerline[-1])


        elif mode == 'energy_ang_spect':

            if (filename == None) & (self.filename == ''):
                outputfn = 'Graph_Angle_Spectrum.txt'
            elif (filename == None) & (self.filename != ''):
                outputfn = self.filename
            else:
                outputfn = filename
                self.filename = outputfn


            colnames = ['theta', 'energy', 'flux']
            self.simdata = pd.read_csv(outputfn, 
                    delimiter='\s+', names=colnames, skiprows=1, dtype=np.float32)

            # get the header line and extract the number of points in the theta and energy mesh
            with open(outputfn,'r') as f:
                headerline = f.readline().split()
            
            if verbose:
                print(headerline)
            
            self.ntheta = np.float(headerline[-2])
            self.nenergy = np.float(headerline[-1])
        
        else:
            print('ERROR: No simulation data loaded.')
            print('Please specify the mode of the simulation output: {"photon_ang_spect", "energy_ang_spect"} ')
            return 0

        return 1

    def integrate_flux(self, anglelims, energylims, jacobian='sphere'):
        """ Integrates the photon/energy flux over angle and energy defined by the limits.
        multiply by the differential crosssection 
            in solid angle (2pi factor for the full phi range) 
            and energy bandwidth (1e3 factor for keV to eV conversion)
        
        self.simdata - df with the angle energy and photon/mrad^2/eV or keV/mrad^2/eV
        anglelims - [angle_min, angle_max] integration limits for the off-axis angle
        energylims - [energy_min, energy_max] integration limits for the photon energy
        jacobian = 'cart', 'sphere' specifies the integration jacobian
            if 'sphere' -> sin(theta) dtheta dphi
            if 'cart' -> dtheta dphi (!!!not correct for most cases!!!)
        
        """
                
        # integrate over angles
        spect_en = self.calc_spectrum(anglelims=anglelims,energylims=energylims, jacobian=jacobian, updateself=False)
        
        envec = spect_en[:,0]
        spectvec = spect_en[:,1]

        try:
            # differential element vector
            delenvec = np.abs(envec[1:] - envec[:-1])
            delenvec = np.append(delenvec, delenvec[-1]) # make same length
        
            photonnumberinrange = np.sum(spectvec * delenvec * 1e3)
            
            return photonnumberinrange

        except IndexError:
            print('ERROR: Not enough resolution in energy for integration.')
            print('integrated flux is zero')
            return 0.0

    def calc_spectrum(self, anglelims, energylims=None, jacobian='sphere',updateself=True):
        """ Integrates the photon/energy flux over angles defined by the limits.
        
        multiply by the differential crosssection 
            in solid angle (2pi factor for the full phi range) 
            and energy bandwidth (1e3 factor for keV to eV conversion)
        
        self.simdata - df with the angle energy and photon/mrad^2/eV or keV/mrad^2/eV
        anglelims - [angle_min, angle_max] integration limits for the off-axis angle
        energylims=None, [energy_min, energy_max] cut-off limits for the photon energy. If `None`, no cut is made on the energy
        
        jacobian = 'cart', 'sphere' specifies the integration jacobian
            if 'sphere' -> sin(theta) dtheta dphi
            if 'cart' -> dtheta dphi (!!!not correct for most cases!!!)
        updateself - True, flag to update the self attribute photons_per_ev
        """
        amin, amax = anglelims
        # apply angular cut
        simdf_angcut = (self.simdata[
                            (self.simdata['theta'] >= amin) 
                            & (self.simdata['theta'] <= amax)]
                            )

        # apply energy cut if supplied and define envec
        if energylims:
            emin, emax = energylims
            simdf_cut = (self.simdata
                            [(self.simdata['energy'] >= emin) 
                            & (self.simdata['energy'] <= emax)]
                            )
            envec = simdf_cut['energy'].unique()

        else:
            simdf_cut = simdf_angcut
            envec = simdf_cut['energy'].unique()

        # init output array
        photons_sumoverangles = np.zeros( (len(envec), 2) )
        photons_sumoverangles[:,0] = envec
        
        for i, en in enumerate(envec):

            # theta vector
            thvec = (simdf_cut
                        ['theta']
                        [simdf_cut['energy'] == en]
                        .to_numpy()
                        )
            # differential element
            delthvec = np.abs(thvec[1:] - thvec[:-1])
            delthvec = np.append(delthvec, delthvec[-1]) # make it the same length as thvec

            if jacobian == 'sphere':
                jacobian_factor = 2 * pc.pi * np.sin(thvec) * delthvec
            elif jacobian == 'cart':
                jacobian_factor = 2 * pc.pi * delthvec
            else:
                print('WARNING! Invalid Jacobian setting provided by user. Integration will be performed with the default settin: jacobian="sphere". This setting can be changed to "cart" to remove the factor of sin(theta)')
                jacobian_factor = 2 * pc.pi * np.sin(thvec) * delthvec
            
            # temp flux vector for the energy, en
            fluxvec = simdf_cut['flux'][simdf_cut['energy'] == en]
            
            # sum over flux vec with the jacobian factor
            photons_sumoverangles[i, 1] = np.sum(jacobian_factor * fluxvec)

        if updateself:
            # update .intspect
            self.photons_per_ev = photons_sumoverangles

        return photons_sumoverangles

    def calc_angdist(self, energylims, anglelims=None):
        """ Integrates the photon/energy flux over angles defined by the limits.
        
        multiply by the differential crosssection 
            in solid angle (2pi factor for the full phi range) 
            and energy bandwidth (1e3 factor for keV to eV conversion)
        
        self.simdata - df with the angle energy and photon/mrad^2/eV or keV/mrad^2/eV
        anglelims - [angle_min, angle_max] integration limits for the off-axis angle
        energylims=None, [energy_min, energy_max] cut-off limits for the photon energy. If `None`, no cut is made on the energy
        
        jacobian = 'cart', 'sphere' specifies the integration jacobian
            if 'sphere' -> sin(theta) dtheta dphi
            if 'cart' -> dtheta dphi (!!!not correct for most cases!!!)
            
        """


        emin, emax = energylims
        # apply cut on energy
        simdf_encut = (self.simdata
                        [(self.simdata['energy'] >= emin) 
                        & (self.simdata['energy'] <= emax)]
                        )

        # apply angular cut if supplied and define thvec
        if anglelims:
            amin, amax = anglelims
            simdf_cut = (simdf_encut
                        ['theta']
                        [(self.simdata['theta'] >= amin) 
                        & (self.simdata['theta'] <= amax)]
                        )
            # vector of theta angles
            thvec = simdf_cut.unique()

        else:
            simdf_cut = simdf_encut
            thvec = simdf_cut['theta'].unique()

        # init output array
        photons_sumoverenergy = np.zeros( (len(thvec), 2) )
        photons_sumoverenergy[:,0] = thvec

        for i, th in enumerate(thvec):

            # energy vector
            envec = (simdf_cut
                        ['energy']
                        [simdf_cut['theta'] == th]
                        .to_numpy()
                        )
            delenvec = np.abs(envec[1:] - envec[:-1])
            delenvec = np.append(delenvec, delenvec[-1]) # make the same len as envec

            # temp flux vector for theta
            fluxvec = simdf_cut['flux'][simdf_cut['theta'] == th]

            # sum over fluxvec with differential element
            photons_sumoverenergy[i,1] = np.sum(delenvec * fluxvec)

        # update .photons_per_mrad
        self.photons_per_mrad = photons_sumoverenergy

        return photons_sumoverenergy

