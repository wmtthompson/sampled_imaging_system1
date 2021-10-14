# Display class
from numpy import *
#zeta = linspace(0,40,num=100)

class Display:
    def __init__(self,zeta,hfov=2.5,vfov=2.3,disp_diag=6,asprat=(4.0/3.0),hres=640,vres=480):
        #zeta = linspace(0,40,num=100)
        #test display with AMLCD resolution 1024x768
        #display diagonal is in mm by conversion from inches
        # Ldisp = display average luminance in fL
        # Lmin = display minimum luminance in fL
        self.Ldisp = 10.0
        self.Lmin = 0.0
        self.Mdsp = (self.Ldisp-self.Lmin)/self.Ldisp
        # deltaL is the display luminance difference peak to peak that the target results in
        # look at Barten to look at the relationship between display modulation and vision
        # there is a contrast modulation that you want to target in order to get the best visual acuity
        # this should be related to the 
        self.deltaL = 8
        self.disp_diag = disp_diag*25.4
        self.disp_asrat = asprat
        self.disp_horz_dim = self.disp_diag/(asprat**2+1)**(1/2)
        self.disp_vert_dim = (self.disp_diag*asprat)/(asprat**2+1)**(1/2)
        self.zeta = zeta
        self.hmag_disp = hfov/self.disp_horz_dim
        self.vmag_disp = vfov/self.disp_vert_dim
        self.hres = hres
        self.vres = vres
        # display horizontal pixel pitch
        self.Hpit = 0.29
        # display vertical pixel pitch
        self.Vpit = 0.31
        # for an LCD display
        self.hmtf = zeros(shape(zeta))
        self.vmtf = zeros(shape(zeta))
        self.calc_horiz_mtf(self.hmag_disp,self.Hpit)
        self.calc_vert_mtf(self.vmag_disp,self.Vpit)

    def calc_hmag_disp():
        pass

    def calc_vmag_disp():
        pass

    def calc_pixel_size():
        pass

    def calc_horiz_mtf(self,hmag_disp,Hpit):
        zeta = self.zeta
        zeta1 = zeta*hmag_disp*17.45
        #Hpit = self.Hpit
        a = sin(0.29*pi*zeta1*Hpit)
        b = 0.29*pi*zeta1*Hpit
        c = (0.68*0.32*cos(0.64*pi*zeta1*Hpit))
        ab = a/b
        abc = ab*c
        self.hmtf = abc

    def calc_scaled_hmtf(self,zeta):
        hmag_disp = self.hmag_disp
        zeta1 = zeta*hmag_disp*17.45
        Hpit = self.Hpit
        a = sin(0.29*pi*zeta1*Hpit)
        b = 0.29*pi*zeta1*Hpit
        c = (0.68*0.32*cos(0.64*pi*zeta1*Hpit))
        ab = a/b
        abc = ab*c
        hmtf = abc
        return hmtf

    def calc_vert_mtf(self,vmag_disp,Vpit):
        #if you want to change the class properties change them before or after
        zeta = self.zeta
        zeta1 = zeta*vmag_disp*17.45
        #Vpit = self.Vpit
        self.vmtf = sin(0.85*pi*zeta1*Vpit)/(1e-20+0.85*pi*zeta1*Vpit)

    def calc_scaled_vmtf(self,zeta):
        vmag_disp = self.vmag_disp
        zeta1 = zeta*17.45
        Vpit = self.Vpit
        a = sin(0.29*pi*zeta1*Vpit)
        b = 0.29*pi*zeta1*Vpit
        c = (0.68*0.32*cos(0.64*pi*zeta1*Vpit))
        ab = a/b
        abc = ab*c
        vmtf = abc
        return vmtf
