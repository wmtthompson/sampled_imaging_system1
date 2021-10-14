#Thermal Optics Class

 

from numpy import *
from matplotlib import pyplot

class therm_choicelist:
    good = 'good'
    typical = 'typical'
    ideal = 'ideal'

class thermal_optics:

    def __init__(self,zeta,choice=therm_choicelist.good):
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
        self.fpa_v_size = fpa_vpixn*fpa_v_pitch
        # fpa horizontal size in meters
        self.fpa_h_size = fpa_hpixn*fpa_h_pitch
        # vertical field of view in degrees
        self.vfov = 2.3
        # horizontal field of view in degrees
        self.hfov = 2.5
        # wavelength of choice
        self.lamd = 5e-6
        self.zeta = zeta
         #frequency units are in cycles/milliradian we chose to go to 40
        # choicelist = 'good','typical','ideal'
        self.fnum = 3
        self.hmtf = self.calc_thermal_Hxdiff(choice,0,self.hfov,self.fnum)
        self.vmtf = self.calc_thermal_Vxdiff(choice,0,self.hfov,self.fnum)
        self.transmision = 0.8

    def plot_mtf(self):
        pyplot.plot(zeta,hmtf)
        hold(True)
        pyplot.plot(zeta,vmtf)
        hold(False)

    def calc_shifted_Hxdiff(self,shift=0):
        return self.calc_thermal_Hxdiff(self.choice,shift,self.vfov,self.fnum)

    def calc_shifted_hmtf(self,shift=0):
        return self.calc_thermal_Hxdiff(self.choice,shift,self.vfov,self.fnum)

    def calc_thermal_Hxdiff(self,choice,shift=0,hfov=2.5,fnum=3.0):
        # the Horizontal MTF due to diffraction
        # diameter D in meters
        # f/# number is 3
        # fov is 2.5 degree
        # hp is the h prime from tan(theta1/2) = hp/f
        # assuming chip is about 0.016 wide
        # light sensitive area of chip is determined from pixel pitch
        # pixel pitch multiplied by number of pixels
        #w = 640*25e-6
        w = self.fpa_h_size
        hp = w/2
        #fov = self.hfov
        f = hp/tan(deg2rad(hfov/2))
        D = f/fnum
        lamd = self.lamd
        outa = zeros(shape(zeta))
        for i, z in enumerate(zeta):
            z=abs(z-shift)
            if (1000*z*lamd)/D>=1:
                outa[i] = 0
            else:
               outa[i] = (2/pi)*(arccos((1000*z*lamd)/D)-((1000*z*lamd)/D)*(sqrt(1-((1000*z*lamd)/D)**2)))
        if choice == therm_choicelist.good:
            out_hmtf = outa**1.2
        elif choice == therm_choicelist.typical:
            out_hmtf = outa**1.75      
        elif choice == therm_choicelist.ideal:
            out_hmtf = outa       
        else:
            print 'not a valid choice, horiz mtf was not calculated!'
        return out_hmtf

    def calc_shifted_Vxdiff(self,shift=0):
        return self.calc_thermal_Vxdiff(self.choice,shift,self.vfov,self.fnum)

    def calc_shifted_vmtf(self,shift=0):
        return self.calc_thermal_Vxdiff(self.choice,shift,self.vfov,self.fnum)
    

    def calc_thermal_Vxdiff(self,choice,shift=0,vfov=2.3,fnum=3.0):
        # the Horizontal MTF due to diffraction
        # diameter D in meters
        # f/# number is 3
        # fov is 2.5 degree
        # hp is the h prime from tan(theta1/2) = hp/f
        # assuming chip is about 0.016 wide
        # light sensitive area of chip is determined from pixel pitch
        # pixel pitch multiplied by number of pixels
        #w = 640*25e-6
        w = self.fpa_v_size
        hp = w/2
        #fov = self.vfov
        f = hp/tan(deg2rad(vfov/2))
        D = f/fnum
        lamd = self.lamd
        outa = zeros(shape(zeta))
        for i, z in enumerate(zeta):
            z= abs(z-shift)
            if (1000*z*lamd)/D>=1:
                outa[i] = 0
            else:
               outa[i] = (2/pi)*(arccos((1000*z*lamd)/D)-((1000*z*lamd)/D)*(sqrt(1-((1000*z*lamd)/D)**2)))
          
        if self.choice == therm_choicelist.good:
            out_vmtf = outa**1.2
        elif self.choice == therm_choicelist.typical:
            out_vmtf = outa**1.75      
        elif self.choice == therm_choicelist.ideal:
            out_vmtf = outa      
        else:
            print 'not a valid choice, vertical mtf was not calculated!'
        return out_vmtf
