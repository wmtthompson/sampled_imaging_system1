from numpy import *

class Target_Therm:
    def __init__(self,target_size=(sqrt(9.7),sqrt(9.7)),Teff=10e3):
        self.target_size = target_size
        self.Teff = Teff
        self.radiance_mean = 0
        self.radiance_variance = 0

    def calc_Wb(self,temp,wl):
        C1= 11934
        C2 = 14388
        # units W cm^-2 micrometer^-1
        return C1/(wl**5*(exp((C2/(wl*temp)))-1))

    def __Wb_Tdiff__(self,temp,wl):
        T = 300
        diff = (self.calc_Wb(wl,T+0.1)-self.calc_Wb(wl,T))/0.1
        return diff



class Background_Therm:
    def __init__(self,background_area=(sqrt(9.7),sqrt(9.7)),Teff=0.08):
        self.background_area = background_area
        self.Teff = Teff
        self.radiance_mean = 0
        self.radiance_variance = 0


class Scene_Therm:
    def __init__(self,target,background,Rng=16):
        self.target = target
        self.background = background
        # Range in meters from target to detector
        self.Rng = Rng
        # Ctgt is target signature radiance in watts/m^2
        # Ctgt target signature is now in watts/cm^2
        #self.Ctgt = 0.004*(1/(1e-2)**2)
        self.Ctgt = 10
        self.Ltgt = target.target_size[0]*target.target_size[1]

    def calc_delta_scene(self,detector,wl1=3,wl2=5):
       T =  300
       C1= 11934
       C2 = 14388
       wls = linspace(wl1,wl2)
       vals = zeros(shape(wls))
       for i, w in enumerate(wls):
           #vals[i] = self.target.__Wb_Tdiff__(T,w)*detector.calc_S(w)
           #vals[i] = self.target.__Wb_Tdiff__(T,w)*1
           vals[i] = (C1*C2*exp((C2/(300*w))))/((300**2)*(w**6)*(exp(C2/(300*w))-1)**2)
       delta_scene = trapz(vals,wls)
       return delta_scene
       
        


