#Reflective Optics class

 

class refl_choicelist:
    good = 'good'
    typical = 'typical'
    WFOV = 'WFOV'

from numpy import *
from matplotlib import pyplot

class Reflective_Optics:

    def __init__(self,zeta,choice=refl_choicelist.good):
        self.choice = choice
        # fpa pixel pitch in meters
        self.fpa_h_pitch = 25e-6
        # fpa pixel pitch in meters
        self.fpa_v_pitch = 25e-6
        # fpa horizontal number of pixels
        self.fpa_hpixn = 640
        # fpa vertical number of pixels
        self.fpa_vpixn = 480
        # fpa vertical size in meters
        self.fpa_v_size = self.fpa_vpixn*self.fpa_v_pitch
        # fpa horizontal size in meters
        self.fpa_h_size = self.fpa_hpixn*self.fpa_h_pitch
        # vertical field of view in degrees
        self.vfov = 2.3
        # horizontal field of view in degrees
        self.hfov = 2.5
        # wavelength of choice
        self.lamd = 0.75e-6
        self.zeta = zeta
        self.fnum = 3
        self.hmtf = zeros(shape(zeta))
        self.vmtf = zeros(shape(zeta))
        self.hmtf = self.calc_reflective_Hxoptic(self.choice,self.hfov,self.fnum,shift=0)
        self.vmtf = self.calc_reflective_Vxoptic(self.choice,self.vfov,self.fnum,shift=0)
        self.transmission = 0.8

 

# def __init__(self,detector_props,choice):

#       self.choice = choice

#       self.fpa_v_size = detector_props.fpa_v_size

#       self.fpa_h_size = detector_props.fpa_h_size


    def plot_mtf(self):
        pyplot.plot(zeta,self.hmtf)
        hold(True)
        pyplot.plot(zeta,self.vmtf)
        hold(False)

    def calc_shifted_Hxoptic(self,shift=0):
        return self.calc_reflective_Hxoptic(self.choice,self.hfov,self.fnum,shift)

    def calc_shifted_hmtf(self,shift=0):
        return self.calc_reflective_Hxoptic(self.choice,self.hfov,self.fnum,shift)

    def calc_focal_length(self):
        w = self.fpa_h_size
        hp = w/2
        f = hp/tan(deg2rad(self.hfov/2))
        self.f_eff = f
        return f

    def calc_reflective_Hxoptic(self,choice,hfov=2.5,fnum=3.0,shift=0):
        zeta = self.zeta
        # shift is the frequency offset
        # the reciporical is the angle offset
        w = self.fpa_h_size
        hp = w/2
        #fov = self.hfov
        f = hp/tan(deg2rad(hfov/2))
        D = f/fnum
        lamd = self.lamd
        outa = zeros(shape(zeta))
        if choice == refl_choicelist.good:
            omega = 0.7
        elif choice == refl_choicelist.typical:
            omega = 0.35     
        elif choice == refl_choicelist.WFOV:
            omega = 0.15     
        else:
            omega = 0.7
            print 'not a valid choice, horizontal mtf default to good!'
        for i, z in enumerate(zeta):
            z = abs(z-shift)
           # z = float(z)
            #z = abs(z)
            if (1000*z*lamd)/(D*omega)>=1:
                outa[i] = 0
            else:
                outa[i] = (2/pi)*(arccos((1000*z*lamd)/(D*omega))-((1000*z*lamd)/(D*omega))*(sqrt(1-((1000*z*lamd)/(D*omega))**2)))
        return outa

    def calc_shifted_Vxoptic(self,shift=0):
        return self.calc_reflective_Vxoptic(self.choice,self.vfov,self.fnum,shift)

    def calc_shifted_vmtf(self,shift=0):
        return self.calc_reflective_Vxoptic(self.choice,self.vfov,self.fnum,shift)

    def calc_reflective_Vxoptic(self,choice,vfov=2.3,fnum=3,shift=0):
        zeta = self.zeta
        w = self.fpa_v_size
        hp = w/2
        #fov = self.vfov
        f = hp/tan(deg2rad(vfov/2))
        D = f/fnum
        lamd = self.lamd
        outa = zeros(shape(zeta))
        if choice == refl_choicelist.good:
            omega = 0.7
        elif choice == refl_choicelist.typical:
            omega = 0.35     
        elif choice == refl_choicelist.WFOV:
            omega = 0.15     
        else:
            omega = 0.7
            print 'not a valid choice, vertical mtf default to good!'
        for i, z in enumerate(zeta):
            z = abs(z-shift)
            if (1000*z*lamd)/(D*omega)>=1:
                outa[i] = 0
            else:
               outa[i] = (2/pi)*(arccos((1000*z*lamd)/(D*omega))-((1000*z*lamd)/(D*omega))*(sqrt(1-((1000*z*lamd)/(D*omega))**2)))
        return outa
