# This contains the class for calculating the targeting task performance metric

from numpy import *
import numpy as np


class TTP:
    def __init__(self,sampled_system,target_scene,observer):
        # SMAG is system magnification
        #disp1.disp_vert_dim
        disp_dim_mm = sampled_system.display.disp_vert_dim
        viewing_distance_mm = observer.viewing_dist
        vfov = sampled_system.optics.vfov
        self.SMAG = sampled_system.calc_SMAG(disp_dim_mm,viewing_distance_mm,vfov)
        # deltaL is display luminance of target from peak to peak
        deltaL = float(sampled_system.display.deltaL)
        # Ldisp is the average display luminance
        Ldisp = float(sampled_system.display.Ldisp)
        # Mdsp = (L-Lmin)/L where L is average display luminance in fL
        # L is minimum luminance (Lmin can be affected by ambient light)
        self.Mdsp = sampled_system.display.Mdsp
        Ctgt = target_scene.Ctgt
        self.Stmp = (Ctgt/deltaL)*Ldisp
        #self.Stmp = (deltaL/Ctgt)*(1/Ldisp)
        self.zeta = sampled_system.zeta
        self.sampled_system = sampled_system
        self.target_scene = target_scene
        self.observer = observer
        # alpha_thermal is the proportionality constant for detectivity models
        self.alpha_thermal = 862
        self.Tcon = float(Ctgt)/(2.0*self.Stmp)
        # in this model LACE (Local Area Contrast Enhancement) is modeled by decreasing the Stmp
        # digital boost is there to compensate for diffraction blur
        # the effect of digital boost can be modeled by modifying diffraction blur such that higher
        # spatial frequencies are incorporated.

    def calc_QH_hor(self):
        # noise is integrated over one sqare milliradian
        # one milliradian corresponds to the sample frequency
        hfov = self.sampled_system.optics.hfov
        horiz_samp_freq = self.sampled_system.detector.calc_horiz_sample_freq(hfov)
        zz = np.absolute(self.zeta)
        SMAG = self.SMAG
        qh_hor = np.zeros(shape(zz))
        zint = np.linspace(0.001,horiz_samp_freq)
        H_dsp = self.sampled_system.display.calc_scaled_hmtf(zint)
        H_eye = self.observer.calc_H_eye_scaled(zint/SMAG)
        dsp_eye = H_dsp*H_eye
        for i, z in enumerate(zz):
            val1 = self.__B__((zint/z))*dsp_eye
            val2 = (abs(val1))**2.0
            integ = np.trapz(val2,zint)
            qh_hor[i] = integ
        return qh_hor

    def calc_QV_ver(self):
        zz = np.absolute(self.zeta)
        vfov = self.sampled_system.optics.vfov
        vert_samp_freq = self.sampled_system.detector.calc_vert_sample_freq(vfov)
        SMAG = self.SMAG
        qv_ver = np.zeros(shape(zz))
        zint = np.linspace(0.001,vert_samp_freq)
        H_eye = self.observer.calc_H_eye_scaled(zint/SMAG)
        V_dsp = self.sampled_system.display.calc_scaled_vmtf(zint)
        dsp_eye = H_eye*V_dsp
        for i, z in enumerate(zz):
            val1 = self.__B__((zint/z))*dsp_eye
            val2 = (abs(val1))**2.0
            integ = np.trapz(val2,zint)
            qv_ver[i] = integ
        return qv_ver

    def calc_QV_hor(self):
        SMAG = self.SMAG
        vfov = self.sampled_system.optics.vfov
        vert_samp_freq = self.sampled_system.detector.calc_vert_sample_freq(vfov)
        zint = np.linspace(0.001,vert_samp_freq)
        H_eye = self.observer.calc_H_eye_scaled(zint/SMAG)
        V_dsp = self.sampled_system.display.calc_scaled_vmtf(zint)
        val1 = V_dsp*H_eye
        val2 = (np.absolute(val1))**2.0
        qv_hor = np.trapz(val2,zint)
        return qv_hor

    def calc_QH_ver(self):
        #zz = np.absolute(self.zeta)
        SMAG = self.SMAG
        hfov = self.sampled_system.optics.hfov
        horiz_samp_freq = self.sampled_system.detector.calc_horiz_sample_freq(hfov)
        zint = np.linspace(0.001,horiz_samp_freq)
        #qh_ver =np.zeros(shape(zz))
        H_eye = self.observer.calc_H_eye_scaled(zint/SMAG)
        H_dsp = self.sampled_system.display.calc_scaled_hmtf(zint)
        val1 = H_dsp*H_eye
        val2 = (np.absolute(val1))**2.0
        qh_ver = np.trapz(val2,zint)
        return qh_ver
   
    def __B__(self,z):
        return self.observer.calc_B(z)

    def __calc_gamma_det__(self):
       fnum = self.sampled_system.optics.fnum
       #scn.calc_delta_scene(self.sampled_system.detector)
       delta_scene = self.target_scene.calc_delta_scene(self.sampled_system.detector)
       #detector1.calc_eta_stare()
       eff_focal_len = self.sampled_system.optics.calc_focal_length()
       # effective focal length comes from the optics calculation in meters
       # but the units used for gamma det are centimeters
       eff_focal_len_cm = eff_focal_len*1e2
       eta_stare = self.sampled_system.detector.calc_eta_stare()
       t_eye = self.observer.t_eye
       trans = self.sampled_system.optics.transmission
       #detector1.calc_D_star(4e-6)
       #d_star_peak = self.sampled_system.detector.calc_D_star(4e-6)
       #d_star_peak is in cm Hz^(0.5) W^-1
       d_star_peak = 1e11
       eff_term = sqrt(eta_stare*t_eye)
       gamma_det = (4*fnum**2)/(1e-3*delta_scene*eff_focal_len_cm*eff_term*d_star_peak*pi*trans)
       return gamma_det
         
    # This is the horizontal aliasing bandwidth
    def calc_QH_alias(self):
        zz = self.zeta
        SMAG = self.SMAG
        qh_alias = zeros(shape(zz))
        H_eye = self.observer.calc_H_eye_scaled(zz/SMAG)
        hfov = self.sampled_system.optics.hfov 
        samp_freq = self.sampled_system.detector.calc_horiz_sample_freq(hfov)
        ah = self.sampled_system.calc_horiz_aliasing_noise(samp_freq)
        for i, z in enumerate(zz):
            val1 = ah[i]*self.__B__(zz/z)*H_eye
            val2 = (abs(val1))**2
            integ = trapz(val2,zz)
            qh_alias[i] = integ
        return qh_alias

    # This is the vertical aliasing bandwidth
    def calc_QV_alias(self):
        zz = self.zeta
        SMAG = self.SMAG
        qv_alias = zeros(shape(zz))
        H_eye = self.observer.calc_H_eye_scaled(zz/SMAG)
        vfov = self.sampled_system.optics.vfov 
        samp_freq = self.sampled_system.detector.calc_vert_sample_freq(vfov)
        av = self.sampled_system.calc_vert_aliasing_noise(samp_freq)
        for i, z in enumerate(zz):
            val1 = av[i]*self.__B__(zz/z)*H_eye
            val2 = (abs(val1))**2
            integ = trapz(val2,zz)
            qv_alias[i] = integ
        return qv_alias

    def calc_CTFHsys(self,SMAG,Rng=16):
        zeta = self.zeta
        ctf_magd = self.observer.calc_ctf(zeta/SMAG)
        hsys = self.sampled_system.calculate_system_hmtf()
        adj_ctf = ctf_magd/(self.Mdsp*hsys)
        qh_hor = self.calc_QH_hor()
        qh_alias = self.calc_QH_alias()
        qv_hor = self.calc_QV_hor()
        gamma_det = self.__calc_gamma_det__()
        Ltgt = self.target_scene.Ltgt
        det_noise_num = ((self.alpha_thermal**2.0)*(gamma_det**2.0)*qh_hor*qv_hor)
        det_noise_dem = self.Stmp**2
        detector_noise = det_noise_num/det_noise_dem
        alpha = self.alpha_thermal*sqrt(self.observer.t_eye)
        alias_noise_num = ((alpha**2.0)*(self.Tcon**2.0)*(Rng**2.0)*(qh_alias))
        alias_noise_dem = (Ltgt*(self.Stmp**2.0))
        aliasing_noise = alias_noise_num/alias_noise_dem
        noise_factor = (1+detector_noise+aliasing_noise)**(0.5)
        ctfh_sys = adj_ctf*noise_factor
        return ctfh_sys

    def calc_CTFVsys(self,SMAG,Rng=16):
        zeta = self.zeta
        ctf_magd = self.observer.calc_ctf(zeta/SMAG)
        hsys = self.sampled_system.calculate_system_vmtf()
        adj_ctf = ctf_magd/(self.Mdsp*hsys)
        qh_ver = self.calc_QH_ver()
        qv_alias = self.calc_QV_alias()
        qv_ver = self.calc_QV_ver()
        # need to calculate detector gamma
        gamma_det = self.__calc_gamma_det__()
        Ltgt = self.target_scene.Ltgt
        alpha = self.alpha_thermal*sqrt(self.observer.t_eye)
        detector_noise = ((self.alpha_thermal)*(gamma_det**2)*qh_ver*qv_ver)/(self.Stmp**2)
        aliasing_noise = ((alpha)*(self.Tcon**2)*(Rng**2)*(qv_alias))/(Ltgt*(self.Stmp**2))
        noise_factor = (1+detector_noise+aliasing_noise)**(0.5)
        ctfv_sys = adj_ctf*noise_factor
        return ctfv_sys

    def calc_alpha(self,ctgt,ctf):
        k = 3
        beta = (k/0.87)
        self.beta = beta
        zeta = self.zeta
        #ctfres = interp(z,zeta,ctf)
        alpha = ctf/(log(2)**(1/beta))
        return alpha

    def calc_p(self,ctgt,ctf):
        alpha = self.calc_alpha(ctgt,ctf)
        prob = (1-exp(-1*(ctgt/alpha)**self.beta))
        prob[np.isnan(prob)] = 0
        return prob

    # TTP output is resolution in cycles per meter
    def calc_TTP(self,CTFH,CTFV,Ctgt,Rng):
        zeta = self.zeta
        #Rng_arry = np.linspace(0,30)
        #self.Rng_arry = Rng_arry
        hprobs = self.calc_p(Ctgt,CTFH)
        hsq = hprobs*(Ctgt/CTFH)*(1/Rng**(2.0))
        hsqrt = np.sqrt(hsq)
        hint = np.trapz(hsqrt,zeta)
        vprobs = self.calc_p(Ctgt,CTFV)
        vsq = vprobs*(Ctgt/CTFV)*(1/Rng**(2.0))
        vsqrt = np.sqrt(vsq)
        vint = np.trapz(vsqrt,zeta)
        vhmult = hint*vint
        ttp = np.sqrt(vhmult)
        return ttp

    def calc_own_TTP(self,SMAG,Rng):
        ctfh =  self.calc_CTFHsys(SMAG,Rng)
        ctfv = self.calc_CTFVsys(SMAG,Rng)
        Ctgt = self.target_scene.Ctgt
        ttp = self.calc_TTP(ctfh,ctfv,Ctgt,Rng)
        return ttp

    def calc_PID(self,ttp):
        # 12 vehicle set phi84 is 14.8'
        ratio = ttp/14.8
        t = np.linspace(0,ratio,num=1000)
        integrd = (2/sqrt(np.pi))*np.exp(-1*t**2.0)
        reslt= np.trapz(integrd,t)
        return reslt

    def calc_PID_wRange(self,Rng_Array):
        PID_arry = np.zeros(shape(Rng_Array))
        for i, r in enumerate(Rng_Array):
            ttp = self.calc_own_TTP(self.SMAG,r)
            PID_arry[i] = self.calc_PID(ttp)
        return PID_arry
        

    
        
    
        
