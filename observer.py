from numpy import *

class Human_Observer:
    def __init__(self, zeta, viewing_distance = 24, display_luminance = 1):
        # viewing distance input is in inches and converted to mm
        self.viewing_dist = viewing_distance*25.4
        self.disp_lum = display_luminance
        self.zeta = zeta
        self.t_eye = 0.1

# teye is probably about 10e-3 seconds
    def calc_B(self, zeta):
        return exp(-2.2*(log10(abs(zeta)))**2)
        

    def calc_H_eye(self,luminance=1):
        #print("guessing luminance =1 for calculating eye MTF")
        L = luminance
        zeta = absolute(self.zeta)
        pup_diam = self.__pupil_diameter__(L)
        f_0 = self.__f_0__(pup_diam)
        i_0 = self.__i_0__(pup_diam)
        MTF_optics = exp(-(zeta/f_0)**(i_0))
        MTF_retina = exp(-0.375*zeta**(1.21))
        MTF_tremor = exp(-0.4441*zeta**(2.0))
        MTF_eye = MTF_optics*MTF_retina*MTF_tremor
        return MTF_eye

    def calc_H_eye_scaled(self,zeta,luminance=1):
        L = luminance
        zeta = absolute(zeta)
        pup_diam = self.__pupil_diameter__(L)
        f_0 = self.__f_0__(pup_diam)
        i_0 = self.__i_0__(pup_diam)
        MTF_optics = exp(-(zeta/f_0)**(i_0))
        MTF_retina = exp(-0.375*zeta**(1.21))
        MTF_tremor = exp(-0.4441*zeta**(2.0))
        MTF_eye = MTF_optics*MTF_retina*MTF_tremor
        return MTF_eye

    
    def __pupil_diameter__(self, light_level):
        diam = [7.0,6.2,5.6,4.9,4.2,3.6,3.0,2.5]
        log_fL = [-4,-3,-2,-1,0,1,2,3]
        res_diam = interp(log10(light_level),log_fL,diam)
        return res_diam

    def __f_0__(self,pup_diam):
        diam = [1.5,2.0,2.4,3.0,3.8,4.9,5.8,6.6]
        f_0 = [36,39,35,32,25,15,11,8]
        out_f_0 = interp(pup_diam,diam,f_0)
        return out_f_0

    def __i_0__(self,pup_diam):
        diam = [1.5,2.0,2.4,3.0,3.8,4.9,5.8,6.6]
        i_0 = [0.9,0.8,0.8,0.77,0.75,0.72,0.69,0.66]
        out_i_0 = interp(pup_diam,diam,i_0)
        return out_i_0
    
    def calc_ctf(self,zeta,Luminance = 1):
        L = Luminance
        w = 15
        zeta = absolute(zeta)
        a = 540*(1+(0.2/L))**(-0.2)/(1+(12.0/((w**2)*(1+5.8*zeta)**2)))
        b = 5.24*(1+(29.2/L))**(0.15)
        return (a*zeta*exp(-b*zeta)*sqrt(1+0.06*exp(b*zeta)))**(-1)

'''zeta = linspace(0,50,num=100)
SMAG = (15.0/24.0)/0.017453293
zeta1 = zeta/SMAG

ctf_response = ctf(zeta1)

import matplotlib
from matplotlib import pyplot
pyplot.plot(zeta1,ctf_response)'''
