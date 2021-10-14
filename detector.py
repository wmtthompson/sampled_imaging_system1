#Detector class
from numpy import *
from matplotlib import pyplot

class Detector_Properties:
    fpa_h_pitch = 30e-5
    # fpa pixel pitch in meters
    fpa_v_pitch = 30e-5
    # fpa horizontal number of pixels
    fpa_v_pixel_size = 28e-5
    fpa_h_pixel_size = 28e-5
    fpa_hpixn = 640
    # fpa vertical number of pixels
    fpa_vpixn = 480
    # fpa vertical size in meters
    fpa_v_size = fpa_vpixn*fpa_v_pitch
    # fpa horizontal size in meters
    fpa_h_size = fpa_hpixn*fpa_h_pitch
    # vertical field of view in degrees
    # dff is the distance that photon penetrates into field-free material
    fpa_dff = 0.0001
    # detector diffusion length in centimeters
    fpa_Ln = 0.01
    # detector noise equivalent power
    # wavelength is in microns
    NEP = (1e-5*ones((200)),linspace(3,5,num=200))
    # number of fields or frames per second
    T_CCD = 60
    # tint is the integration time
    t_int = 0.01
    

class Detector:

    def __init__(self,zeta,detector_props=Detector_Properties(),hfov=2.5,vfov=2.3):
        self.fpa_h_pitch = detector_props.fpa_h_pitch
        self.fpa_v_pitch = detector_props.fpa_v_pitch
        self.fpa_hpixn = detector_props.fpa_hpixn
        self.fpa_vpixn = detector_props.fpa_vpixn
        self.fpa_v_size = detector_props.fpa_v_size
        self.fpa_h_size = detector_props.fpa_h_size
        self.fpa_v_pixel_size = detector_props.fpa_v_pixel_size
        self.fpa_h_pixel_size = detector_props.fpa_h_pixel_size
        self.fpa_dff = detector_props.fpa_dff
        self.fpa_Ln = detector_props.fpa_Ln
        self.fpa_hmag = hfov/(self.fpa_h_size/1e-3)
        self.fpa_vmag = vfov/(self.fpa_v_size/1e-3)
        self.zeta = zeta
        self.fpa_size_units = "meters"
        self.hmtf = self.calc_det_spatial_hmtf(self.fpa_hmag,self.fpa_h_pixel_size)
        self.vmtf = self.calc_det_spatial_vmtf(self.fpa_vmag,self.fpa_v_pixel_size)
        #self.calc_det_diffusion_hmtf(self.fpa_hmag,self.fpa_Ln,self.fpa_dff)
        #self.calc_det_diffusion_vmtf(self.fpa_vmag,self.fpa_Ln,self.fpa_dff)
        #self.calc_comb_det_hmtf()
        #self.calc_comb_det_vmtf()
        self.NEP = detector_props.NEP
        self.T_CCD = detector_props.T_CCD
        self.t_int = detector_props.t_int

    def calc_eta_stare(self):
        H_det = self.fpa_h_pixel_size
        V_det = self.fpa_v_pixel_size
        H_pit = self.fpa_h_pitch
        V_pit = self.fpa_v_pitch
        eta_stare = (self.t_int*self.T_CCD*H_det*V_det)/(H_pit*V_pit)
        return eta_stare

    def plot_mtf(self):
        pyplot.plot(zeta,self.hmtf)
        hold(True)
        pyplot.plot(zeta,self.vmtf)
        hold(False)

    def calc_fpa_mag(fpa_dim,fov_dim):
        return fov_dim/(fpa_dim/1e-3)

    def calc_pixel_size():
        pass

    def calc_D_star(self,wavelength):
        #wl = self.NEP[1]
        #D_lambda = 1/self.NEP[0]
        # NEP units are in watts/hz^(1/2)
        D_lambda = 1/2.71e-11
        # A_det is in m^2?
        v_pix_size_cm = self.fpa_v_pixel_size*1e-2
        h_pix_size_cm = self.fpa_h_pixel_size*1e-2
        # A_det size is now in cm^2
        A_det = v_pix_size_cm * h_pix_size_cm
        # delta_f is temporal bandwidth based on interlaced system
        delta_f = 30
        D_star = D_lambda*(A_det*delta_f)**(.5)
        #return interp(wavelength,wl,D_star)
        return D_star

    def calc_S(self,wavelength):
        wl = self.NEP[1]
        vals = self.calc_D_star(wl)
        peak = vals.max()
        return self.calc_D_star(wavelength)/peak



    def calc_horiz_sample_freq(self,hfov):
        # calculates sample frequency in samples per milliradian
        # given field of view in degrees
        # 17.45 is milliradians per degree
        return (self.fpa_hpixn/hfov)*(1/17.45)

    def calc_vert_sample_freq(self,vfov):
        # calculates sample frequency in samples per milliradian
        # given field of view in degrees
        # 17.45 is milliradians per degree
        return (self.fpa_vpixn/vfov)*(1/17.45)

    def calc_det_spatial_hmtf(self,fpa_hmag,fpa_h_pixel_size,shift=0):
        #need the horizontal display size
        #cy/mm = cycles/milrad * deg/mm * 17.45 milrads/deg
        #deg/mm is fpa mag factor
        #mult by 0.1 to get from mm to cm
        zeta = self.zeta
        zeta = (zeta-shift)
        zeta1 = zeta*fpa_hmag*17.45*10
        # convert horz detector size to centimeters
        h_det = fpa_h_pixel_size/1e-6
        spatial_hmtf = sin(0.0001*pi*zeta1*h_det)/(0.0001*pi*zeta1*h_det)
        return spatial_hmtf

    def calc_det_spatial_vmtf(self,fpa_vmag,fpa_v_pixel_size,shift=0):
        #need the horizontal display size
        #cy/mm = cycles/milrad * deg/mm * 17.45 milrads/deg
        #deg/mm is displ mag factor
        #mult by 0.1 to get from mm to cm
        zeta = self.zeta
        zeta = (zeta-shift)
        zeta1 = zeta*fpa_vmag*17.45*10
        # convert horz detector size to centimeters
        h_det = fpa_v_pixel_size/1e-6
        spatial_vmtf = sin(0.0001*pi*zeta1*h_det)/(0.0001*pi*zeta1*h_det)
        return spatial_vmtf

    def calc_det_diffusion_hmtf(self,fpa_hmag,fpa_Ln,fpa_dff,shift=0):
        zeta = self.zeta
        zeta = (zeta-shift)
        zeta1 = zeta*fpa_hmag*17.45*0.1
        L = (((1/fpa_Ln**2)+4*pi**2*zeta1**2)**-1)**(1/2)
        diffusion_hmtf = (exp((fpa_dff/L))+exp((-fpa_dff/L)))/(exp((fpa_dff/L))+exp((-fpa_dff/L)))
        return diffusion_hmtf

    def calc_det_diffusion_vmtf(self,fpa_vmag,fpa_Ln,fpa_dff,shift=0):
        zeta = self.zeta
        zeta = (zeta-shift)
        zeta1 = zeta*fpa_vmag*17.45*0.1
        L = (((1/fpa_Ln**2)+4*pi**2*zeta1**2)**-1)**(1/2)
        diffusion_vmtf = (exp((fpa_dff/L))+exp((-fpa_dff/L)))/(exp((fpa_dff/L))+exp((-fpa_dff/L)))
        return diffusion_vmtf

    def calc_comb_det_hmtf(self):
        # this is the last method to be called for horizontal mtf
        zeta = self.zeta
        zz = ones(shape(zeta))
        hmtflist = list()
        hmtflist.append(self.spatial_hmtf)
        #hmtflist.append(self.diffusion_hmtf)
        self.hmtflist = hmtflist
        for el in self.hmtflist:
            for i, z in enumerate(el):
                zz[i] = zz[i]*z
        self.hmtf = zz

    def calc_comb_det_vmtf(self):
        # this is the last method to be called for vertical mtf
        zeta = self.zeta
        zz = zeros(shape(zeta))
        vmtflist = list()
        vmtflist.append(self.spatial_vmtf)
        #vmtflist.append(self.diffusion_vmtf)
        self.vmtflist = vmtflist
        for el in self.vmtflist:
            for i, z in enumerate(el):
                zz[i] = zz[i]*z
        self.vmtf = zz

    def calc_shifted_hmtf(self,shift=0):
        return self.calc_shifted_comb_det_hmtf(shift)

    def calc_shifted_vmtf(self,shift=0):
        return self.calc_shifted_comb_det_vmtf(shift)

    def calc_shifted_comb_det_hmtf(self,shift=0):
        shifted_hmtf_list = list()
        spatial_shifted_hmtf = self.calc_det_spatial_hmtf(self.fpa_hmag,self.fpa_h_pixel_size,shift)
        #diffusion_shifted_hmtf = self.calc_det_diffusion_hmtf(self,self.fpa_hmag,self.fpa_Ln,self.fpa_dff,shift)
        shifted_hmtf_list.append(spatial_shifted_hmtf)
        #shifted_hmtf_list.append(diffusion_shifted_hmtf)
        zeta = self.zeta
        zz = ones(shape(zeta))
        #hmtflist.append(self.diffusion_hmtf)
        for el in shifted_hmtf_list:
            for i, z in enumerate(el):
                zz[i] = zz[i]*z
        shifted_hmtf = zz
        return shifted_hmtf

    def calc_shifted_comb_det_vmtf(self,shift=0):
        shifted_vmtf_list = list()
        spatial_shifted_vmtf = self.calc_det_spatial_vmtf(self.fpa_vmag,self.fpa_v_pixel_size,shift)
        #diffusion_shifted_vmtf = self.calc_det_diffusion_vmtf(self.fpa_vmag,self.fpa_Ln,self.fpa_dff,shift)
        shifted_vmtf_list.append(spatial_shifted_vmtf)
        #shifted_hmtf_list.append(diffusion_shifted_hmtf)
        zeta = self.zeta
        zz = ones(shape(zeta))
        #hmtflist.append(self.diffusion_hmtf)
        for el in shifted_vmtf_list:
            for i, z in enumerate(el):
                zz[i] = zz[i]*z
        shifted_vmtf = zz
        return shifted_vmtf

   
