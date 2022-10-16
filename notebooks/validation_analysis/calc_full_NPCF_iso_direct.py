#########################################################################
#                 Jiamin Hou Jan 26 2021

#   The code is based on a set of isotropic function basis.
#   It is served as a small scale test of the NPCF functions.
#   N=2 ..5 are given in arXiv:2010.14418 A.1-A.4


#########################################################################

import sys, os, time, pickle
import collections, numpy
from scipy import spatial, stats
from itertools import combinations 

class calc_NPCF(object):
    
    def __init__(self, npcf=None, data=None, nbins=None, lbox=None, rmax=None, lls=None, verbose=False):
        
        self.npcf  = npcf
        self.nbins = nbins
        self.lbox  = lbox
        self.rmax  = rmax
        self.verbose = verbose
        
        if data is not None:
            self.ngals = data.shape[0]
            self.data = data
        else:
            raise AssertionError("Give me Data!")

        if lls is not None:
            self.lls = lls
        else:
            if self.npcf == 3:
                self.lls = ['00', '11', '22', '33', '44']
            elif self.npcf == 4:
                self.lls = ['110', '111', '112', '222']
            elif self.npcf == 5:
                self.lls = ['11(0)11', '10(1)01', '21(1)12', '21(2)12',
                            '11(2)11', '11(1)11', '21(1)10', '21(3)12',
                            '21(2)11', '21(1)11']
            elif self.npcf == 6:
                self.lls = ['00(0)0(0)0', '11(0)0(0)0', '00(0)1(1)10','10(1)1(1)01', '10(1)1(2)02']

    
    def P_lambda(self, gals_np, order='11(0)11'):
        '''
            define basis function:
            
            Slepian&Cahn 2020 Appendix A
        '''
        
        if self.npcf == 3:
            r1_hat = gals_np[0]/numpy.linalg.norm(gals_np[0], ord=2)
            r2_hat = gals_np[1]/numpy.linalg.norm(gals_np[1], ord=2)
            
            x = numpy.dot(r1_hat, r2_hat)
            
            if order == '00':
                return 1./4/numpy.pi
            
            elif order == '11':
                return -numpy.sqrt(3)/4./numpy.pi * x
            
            elif order == '22':
                return 3./2*numpy.sqrt(5./(4*numpy.pi)**2) * (x**2 - 1./3)
            
            elif order == '33':
                factor = - (numpy.sqrt(7) / (4 * numpy.pi)) * (3/2)
                return factor * ( (5/3) * x**3 - x )
            
            elif order == '44':
                factor = (3 / (4 * numpy.pi)) * (1/8)
                return factor * (35 * x**4 - 30 * x**2 + 3)
                        
        elif self.npcf == 4:
            r1_hat = gals_np[0]/numpy.linalg.norm(gals_np[0], ord=2)
            r2_hat = gals_np[1]/numpy.linalg.norm(gals_np[1], ord=2)
            r3_hat = gals_np[2]/numpy.linalg.norm(gals_np[2], ord=2)
         
            if order == '110':
                return -1./4/numpy.pi*numpy.sqrt(3./4/numpy.pi) * numpy.dot(r1_hat, r2_hat)

            elif order == '111':
                return -3j/numpy.sqrt(2)/(4 * numpy.pi)**(3./2) * numpy.dot(r1_hat, numpy.cross(r2_hat, r3_hat))

            elif order == '112':
                return numpy.sqrt(27/2/(4*numpy.pi)**3) * (numpy.dot(r1_hat, r3_hat)*numpy.dot(r2_hat, r3_hat)-
                                                          1./3 * numpy.dot(r1_hat, r2_hat))
            elif order == '222':
                return -45/numpy.sqrt(14*(4*numpy.pi)**3) * (numpy.dot(r1_hat, r2_hat)*numpy.dot(r1_hat, r3_hat)*numpy.dot(r2_hat, r3_hat) 
                                                              -1./3*numpy.dot(r1_hat, r2_hat)**2 
                                                              -1./3*numpy.dot(r1_hat, r3_hat)**2
                                                              -1./3*numpy.dot(r2_hat, r3_hat)**2 + 2./9)
        
        else:
            raise AssertionError("Unacceptable npcf value (3/4 only).")
                       


    def init_coeff(self):
        
        '''initialize coefficients'''

        self.zeta = {}
        for il in self.lls:
            self.zeta[il] = numpy.zeros(([self.nbins]*int(self.npcf-1)), dtype='c16')
            
    def save(self, sname, infile=None):
        
        '''save option:
        
           if no input: save coefficients as pickle file
        
           if input given: save as .txt file
        
        '''
        
        if infile is None:       
            with open(sname + '.pkl', 'wb') as f:
                pickle.dump(self.zeta, f, -1)
            print(">> save file to", sname + '.pkl')
                
        else:
            numpy.savetxt(sname + '.txt', infile)
            print(">> save file to", sname + '.txt')
        
            
    def make_catalog(self):
        
        '''
           make catalog based on the galcoords
           columns:
                   pos1, pos2, pos3, weight
        '''
        
        self.catalog = numpy.zeros([len(self.galcoords), 4])
        self.catalog[:,:3] = self.galcoords
        self.catalog[:,-1] = 1.


    def run(self):
        
        numpy.random.seed(10)
        eps = 1e-8
        
        # self.galcoords = self.lbox * numpy.random.rand(self.ngals, 3)
        self.galcoords = self.data[:,0:3] # x,y,z coords
        
        def shiftbox(gal_coords, boxsize,xshift,yshift,zshift):
            galxsh = gal_coords[:,0] + xshift*boxsize
            galysh = gal_coords[:,1] + yshift*boxsize
            galzsh = gal_coords[:,2] + zshift*boxsize
            gal_coords_new = numpy.array([galxsh, galysh, galzsh]).transpose()
            return gal_coords_new

        galcoords = self.galcoords
        #now loop over shiftbox to append virtual boxes to coordinate list.
        for xshift in (-1,1,0):
            for yshift in (-1,1,0):
                for zshift in (-1,1,0):
                    # print(xshift,yshift,zshift)
                    if xshift == 0 and yshift == 0 and zshift == 0:
                        print("no shift executed because in genuine box")
                    else:
                        galcoords = numpy.append(galcoords, shiftbox(self.galcoords,self.lbox,xshift,yshift,zshift),axis=0) 

        self.galcoords = galcoords
        
        tree = spatial.cKDTree(self.galcoords,leafsize=self.lbox)
        
        self.centralgals = self.galcoords[:self.ngals,:] 
        
        self.sides = []
        for ii in range(self.ngals):
            if self.verbose: print (f">> sit on {ii}th galaxy")
            prime_coord_3d = self.centralgals[ii]
            # ballct: list of indices of the secondary galaxies of
	    # the current primary galaxy
            ballct = tree.query_ball_point(prime_coord_3d, self.rmax+eps)
            ballct.remove(ii)
            # coordinate of the secondary galaxies 
            second_coord_3d = self.galcoords[ballct] - prime_coord_3d
            second_coord_radial = numpy.sqrt(numpy.sum(second_coord_3d**2, axis=1))
            # sort the secondary galaxies according to their distance
            second_coord_3d = second_coord_3d[second_coord_radial.argsort()]
            second_coord_3d = second_coord_3d[::-1]
            second_coord_radial = second_coord_radial[second_coord_radial.argsort()]
            second_coord_radial = second_coord_radial[::-1]
            # find all possible combinations for N-1 neighbouring galaxies 
            # out of total secondary galaxies. comb_ballct is an object that 
            # contains a list of tuples of the indices of secondary galaxies
            comb_ballct = combinations(numpy.arange(len(second_coord_3d)), int(self.npcf-1))
            for icomb_ballct in list(comb_ballct):
                if self.verbose: print("  local secondary indices:", icomb_ballct)
                # get the current N-1 secondary galaxies
                second_coord_3d_np = second_coord_3d[list(icomb_ballct)]
                second_coord_radial_np = second_coord_radial[list(icomb_ballct)]
                
                self.second_coord_3d_np = second_coord_3d_np
                if self.verbose: print(f"  position of {self.npcf-1:1d} secondary galaxies:\n", second_coord_3d_np)
                # get the radial bin index of the N-1 secondary galaxies
                # proceed to add the result to zeta coefficients only if
                # these N-1 secondaries are not in the same bin
                sides = [int(x) for x in list(second_coord_radial_np/self.rmax*self.nbins)]
                if len(sides) == len(set(sides)):
                    for il in self.lls:
                        self.sides.append(sides)
                        if self.verbose: print("    order", il, "coefficients", 
                                               self.P_lambda(second_coord_3d_np, order=il))
                        self.zeta[il][tuple(sides)] += self.P_lambda(second_coord_3d_np, order=il) 
                else:
                    if self.verbose: print(f"    {len(sides)-len(set(sides))+1:1d} galaxies fall in the same radial bin")
