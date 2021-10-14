#Sampled imaging system

# The following imports all the packages 
from time import clock, time
start = clock()
from numpy import *
import numpy as np
import matplotlib as mp
import pylab 
from matplotlib import pyplot
# The remaining import statement relate to 'modules' that go with the sampled system
import thermal_optics as to 
import reflective_optics as ro 
import detector as detr 
import display as displ
import observer as obsvr 
import scene_thermal_model as scene 
import TTP_metric 
import samp_interpolation as sintp
import copy
# Reload the packages to make sure that we have to most recent
reload(obsvr)
reload(displ)
reload(detr)
reload(scene)
reload(TTP_metric)
reload(to)
reload(ro)

# zeta is the spatial frequency space in cy/milliradian
zeta = np.linspace(-50,50,num=400)
# the element that is zero is replaced to avoid divide by zero issues later on
zeta[zeta==0] = 1e-3

# The following is the class Sampled System which is the class that defines the whole sampled system
# This does not include the observer or the scene. The constructor __init__ shows th at the sampled
# system is composed of opics, detector and display as a minimum.
class Sampled_System:
    def __init__(self,zeta,optics,detector,display):
        self.optics = optics
        self.detector = detector
        self.display= display
        # the elementlist is a list of elements that are contained in the system
        # in future implementations it would be nice to add elements to a system once
        # instead of adding it separately by equating it to the property say display
        # and also adding that display to the element list
        # this could be done by having one 'append' method and then the method checks what kind of
        # class is being added
        self.elementlist = list()
        self.zeta = zeta
        # presample lists are divided by horz and vert for pre and post
        # presample items would be optics and detector and and any other filtering or conversions
        # post sample items would be items like interpolation and display, as well as any other electronic filters
        self.presamplehmtflist = list()
        self.presamplevmtflist = list()
        self.postfilter_hmtflist = list()
        self.postfilter_vmtflist = list()
        # there are optional spaces for display mtfs but may not be used
        self.display_hmtf = np.zeros(shape(zeta))
        self.display_vmtf = np.zeros(shape(zeta))

    # This method calculates the system magnification, this method is called in the TTPmetric class
    # vert_disp is the vertical display size
    # viewing_dist is the viewing distance from the display
    # vfov_deg is the vertical field of view in degrees
    def calc_SMAG(self, vert_disp,viewing_dist,vfov_deg=2.5):
        vfov_rad = np.deg2rad(vfov_deg)
        return (vert_disp/viewing_dist)/vfov_rad
        
    # This calculates the post-sample horiz MTF based on internal properties and sample frequency
    # Postsample means after sampling by the detector
    # allows exploration of postsample hmtf by changing sample frequency (in samples/milliradian)
    def calc_post_sample_hmtf(self,sample_freq):
        hmtflist = list()
        temparray = np.ones(shape(self.zeta))
        nums = range(-2,3)
        for n in nums:
            for el in self.presamplehmtflist:
                temparray= temparray*el.calc_shifted_hmtf(n*sample_freq)
            hmtflist.append(temparray)
        shifted_mtx = asmatrix(hmtflist)
        postsample_hmtf = np.sum(shifted_mtx, axis=0)
        return postsample_hmtf
    
    # calculates the spurious response depending on the sample frequency, which allows calculation of system
    # spurious response but also allows for exploring of spurious response for the the same elements in the system
    # but with a changed sample frequency (in samples /milliradian)
    # the difference between this and post-sample mtf is that the center or unshifted mtfs are not factored in
    # this is made to simply show the intrusion of shifted spectrum into the transfer response of the postfilter mtfs
    # postfilter mtfs would be those from interpolation and/or display
    def calc_spr_hresp(self,sample_freq):
        hmtflist = list()
        temparray = np.ones(shape(self.zeta))
        nums = range(-2,3)
        nums.remove(0)
        for n in nums:
            for el in self.presamplehmtflist:
                temparray= temparray*el.calc_shifted_hmtf(n*sample_freq)
            hmtflist.append(temparray)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(hmtflist)
        hsum = np.sum(shifted_mtx, axis=0)
        hsum = np.array(hsum[0,:], dtype=float)
        hsum = hsum[0,:]
        hspur = hsum
        for el in self.postfilter_hmtflist:
            hspur = el.hmtf*hspur
        return hspur


    # spurious response for vertical response 
    def calc_spr_vresp(self,sample_freq):
        vmtflist = list()
        temparray = np.ones(shape(self.zeta))
        nums = range(-2,3)
        nums.remove(0)
        for n in nums:
            for el in self.presamplevmtflist:
                temparray= temparray*el.calc_shifted_vmtf(n*sample_freq)
            vmtflist.append(temparray)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(vmtflist)
        vsum = np.sum(shifted_mtx, axis=0)
        vsum = np.array(vsum[0,:], dtype=float)
        vsum = vsum[0,:]
        vspur = vsum
        for el in self.postfilter_vmtflist:
            vspur = vspur*el.vmtf
        return vspur

    # calculates the total response of the system without removing the center or transfer response
    def calc_tot_hresp(self,sample_freq):
        hmtflist = list()
        temparray = np.ones(shape(self.zeta))
        nums = range(-2,3)
        #nums.remove(0)
        for n in nums:
            for el in self.presamplehmtflist:
                temparray= temparray*el.calc_shifted_hmtf(n*sample_freq)
            hmtflist.append(temparray)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(hmtflist)
        hsum = np.sum(shifted_mtx, axis=0)
        hsum = np.array(hsum[0,:], dtype=float)
        hsum = hsum[0,:]
        hresp = hsum
        for el in self.postfilter_hmtflist:
            hresp = el.hmtf*hresp
        return hresp

    # calculates the total vertical response without removing the center or transfer response
    def calc_tot_vresp(self,sample_freq):
        vmtflist = list()
        temparray = np.ones(shape(self.zeta))
        nums = range(-2,3)
        #nums.remove(0)
        for n in nums:
            for el in self.presamplevmtflist:
                temparray= temparray*el.calc_shifted_vmtf(n*sample_freq)
            vmtflist.append(temparray)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(vmtflist)
        vsum = np.sum(shifted_mtx, axis=0)
        vsum = np.array(vsum[0,:], dtype=float)
        vsum = vsum[0,:]
        vresp = vsum
        for el in self.postfilter_vmtflist:
            vresp = vresp*el.vmtf
        return vresp

    # calculates horizontal aliasing noise based on the sampling frequency and other factors internally
    def calc_horiz_aliasing_noise(self,sample_freq):
        hmtflist = list()
        temparray = np.ones(shape(self.zeta))
        nums = range(-2,3)
        nums.remove(0)
        #hold(True)
        for n in nums:
            for el in self.presamplehmtflist:
                temparray= temparray*el.calc_shifted_hmtf(n*sample_freq)
            hmtflist.append(temparray**2)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(hmtflist)
        postsample_hmtf = sum(shifted_mtx, axis=0)
        # power(postsample_hmtf,0.5) would also work
        # power(postsample_hmtf,(1/2)) doesnt work because
        # (1/2) is evaluated to 0
        postsample_hmtf = sqrt(postsample_hmtf)
        disp_hmtf_mtx = asmatrix(self.display.hmtf)
        aliasing_noise = multiply(postsample_hmtf,disp_hmtf_mtx)
        aliasing_noise_arry = np.array(aliasing_noise[0,:], dtype=float)
        aliasing_noise_arry = aliasing_noise_arry[0,:]
        #aliasing_noise = postsample_hmtf*asmatrix(self.display_hmtf)
        #hold(False)
        #aliasing_noise_array
        return aliasing_noise_arry

    # calculates the vertical aliasing noise with sample frequency taken into account
    def calc_vert_aliasing_noise(self,sample_freq):
        vmtflist = list()
        temparray = ones(shape(self.zeta))
        nums = range(-2,3)
        nums.remove(0)
        for n in nums:
            for el in self.presamplevmtflist:
                temparray= temparray*el.calc_shifted_vmtf(n*sample_freq)
            vmtflist.append(temparray**2)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(vmtflist)
        postsample_vmtf = sum(shifted_mtx, axis=0)
        postsample_vmtf = sqrt(postsample_vmtf)
        disp_vmtf_mtx = asmatrix(self.display.vmtf)
        aliasing_noise = multiply(postsample_vmtf,disp_vmtf_mtx)
        aliasing_noise_arry = np.array(aliasing_noise[0,:], dtype=float)
        aliasing_noise_arry = aliasing_noise_arry[0,:]
       #aliasing_noise = postsample_vmtf*asmatrix(self.display_vmtf)
        return aliasing_noise_arry

    
    # post sample mtf calculates mtf after sampling, note that this result is not multiplied by the display mtf
    def calc_post_sample_vmtf(self,sample_freq):
        vmtflist = list()
        temparray = ones(shape(self.zeta))
        nums = range(-2,3)
        for n in nums:
            for el in self.presamplevmtflist:
                temparray= temparray*el.calc_shifted_vmtf(n*sample_freq)
            vmtflist.append(temparray)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(vmtflist)
        postsample_vmtf = sum(shifted_mtx, axis=0)
        return postsample_vmtf
    
    # this method was written for future expansion where there may be an element in the system that converts
    # to analog, and then an element that then samples that analog signal and converts it back into a digital signal
    # for display.
    def calc_post_sample_hmtf2(self,first_sample_freq,second_sample_freq):
        #this function is for inserting an extra sampling step that results from conversion to RS170
        hmtflist = list()
        temparray = ones(shape(self.zeta))
        nums = range(-2,3)
        mums = range(-2,3)
        # mums is the second sample operation
        for m in mums:
            # nums is the first sample operation
            for n in nums:
                for el in self.presamplehmtflist:
                    temparray = temparray*el.calc_shifted_hmtf(n*first_sample_freq+m*second_sample_freq)
            for el in self.postfilter_hmtflist:
                temparray = temparray*el.calc_shifted_hmtf(m*second_sample_freq)
            hmtflist.append(temparray)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(hmtflist)
        postsample_hmtf = sum(shifted_mtx,axes=0)
        return postsample_hmtf

     # this method was written for future expansion where there may be an element in the system that converts
    # to analog, and then an element that then samples that analog signal and converts it back into a digital signal
    # for display.
    def calc_post_sample_vmtf2(self,first_sample_freq,second_sample_freq):
        #this function is for inserting an extra sampling step that results from conversion to RS170
        vmtflist = list()
        temparray = ones(shape(self.zeta))
        nums = range(-2,3)
        mums = range(-2,3)
        # mums is the second sample operation
        for m in mums:
            # nums is the first sample operation
            for n in nums:
                for el in self.presamplevmtflist:
                    temparray = temparray*el.calc_shifted_vmtf(n*first_sample_freq+m*second_sample_freq)
            for el in self.postfilter_vmtflist:
                temparray = temparray*el.calc_shifted_vmtf(m*second_sample_freq)
            vmtflist.append(temparray)
            temparray = ones(shape(self.zeta))
        shifted_mtx = asmatrix(vmtflist)
        postsample_vmtf = sum(shifted_mtx,axes=0)
        return postsample_vmtf
                
    # calculates the entire system vmtf using every element in the elementlist, this does not include any sampling
    # this is simply the system transfer response without any aliasing. If the detector mtf is in the element list,
    # the detector's mtf will be used as part of the system vmtf. It could be what you want because you can use the
    # detector's mtf without looking at the result of aliasing from sampling.
    def calculate_system_vmtf(self):
        zeta = self.zeta
        zz = ones(shape(zeta))
        for el in self.elementlist:
            for i, z in enumerate(el.vmtf):
                zz[i] = zz[i]*z
        self.vmtf = zz
        return zz

    # calculates the entire system hmtf using every element in the elementlist, this does not include any sampling
    # this is simply the system transfer response without any aliasing. If the detector mtf is in the element list,
    # the detector's mtf will be used as part of the system hmtf. It could be what you want because you can use the
    # detector's mtf without looking at the result of aliasing from sampling.
    def calculate_system_hmtf(self):
        zeta = self.zeta
        zz = ones(shape(zeta))
        for el in self.elementlist:
            for i, z in enumerate(el.hmtf):
                zz[i] = zz[i]*z
        self.hmtf = zz
        return zz

# Now we define the elements to go into an instance of our system
# setting up the optics with default characteristics of a good MTF with field of view of 2.5 degrees horiz 2.3 vert
refloptics1 = ro.Reflective_Optics(zeta)
# using the default detector properties with 640x480 resolution
Det_Props = detr.Detector_Properties()
detector1 = detr.Detector(zeta)
hfov = refloptics1.hfov
vfov = refloptics1.vfov
# getting the sample frequencies from the detector properties and the optics field of view
hsamp_freq = detector1.calc_horiz_sample_freq(hfov)
vsamp_freq = detector1.calc_vert_sample_freq(vfov)
# The first display has diagonal of 6 aspect ratio of 4/3 and hres=640 vres=480
disp1 = displ.Display(zeta,hfov,vfov,6,(4.0/3.0),640,480)
# Now we put these elements into a system to instantiate the system
system1 = Sampled_System(zeta,refloptics1,detector1,disp1)
system1.elementlist.append(refloptics1)
system1.elementlist.append(detector1)
system1.presamplehmtflist.append(refloptics1)
system1.presamplehmtflist.append(detector1)
# we have made system1 up to the point we want to so that we can copy it.

import copy
# system2 is a deep copy of system1 up to the interpolation and the display
# deep copy is used to separate the associations of one system from affecting the other system
system2 = copy.deepcopy(system1)
# a different display is now defined with larger size and with higher resolution than disp1
disp2 = displ.Display(zeta,hfov,vfov,15,(4.0/3.0),1024,768)
# display1 is added to the element list of system1
system1.elementlist.append(disp1)
# display 1 is added to the post-mtf lists for system1
system1.postfilter_hmtflist.append(disp1)
system1.postfilter_vmtflist.append(disp1)
# interpolation is going to be used in system 2 because system2 display has more pixels than the detector resolution
interp = sintp.Interpolation(zeta,hsamp_freq,vsamp_freq)
# display and interpolation are added
system2.elementlist.append(interp)
system2.elementlist.append(disp2)
system2.display = disp2
system2.postfilter_hmtflist.append(interp)
system2.postfilter_vmtflist.append(interp)
system2.postfilter_hmtflist.append(disp2)
system2.postfilter_vmtflist.append(disp2)

system1.calculate_system_hmtf()
system1.calculate_system_vmtf()
system2.calculate_system_hmtf()
system2.calculate_system_vmtf()

# spurious response for system1
sr1 = system1.calc_spr_hresp(hsamp_freq)

# spurious response for system2
sr2 = system2.calc_spr_hresp(hsamp_freq)

# create elements for calculating the Targeting Task Performance metric
tgt = scene.Target_Therm()
observer1 = obsvr.Human_Observer(zeta)
bckgrd = scene.Background_Therm()
scn = scene.Scene_Therm(tgt,bckgrd)
sys1TTP = TTP_metric.TTP(system1,scn,observer1)
sys2TTP = TTP_metric.TTP(system2,scn,observer1)

Rng_array = linspace(1,37,num=100)

PID1 = sys1TTP.calc_PID_wRange(Rng_array)
PID2 = sys2TTP.calc_PID_wRange(Rng_array)

pyplot.plot(Rng_array,PID1)
pyplot.hold(True)
pyplot.plot(Rng_array,PID2)
pyplot.hold(False)

# What is the diffference in PID at a distance of 37 km?
diffs = PID2[99]-PID1[99]
print "difference in PID is "+str(diffs)

# What is the difference in resolution?
ttp1_37 = sys1TTP.calc_own_TTP(sys1TTP.SMAG,37)
ttp2_37 = sys2TTP.calc_own_TTP(sys2TTP.SMAG,37)
print "system1 resolution at 20 Nautical Miles= "+str(ttp1_37)+"cycles/meter"
print "system2 resolution at 20 Nautical Miles= "+str(ttp2_37)+"cycles/meter"
diffr = ttp2_37-ttp1_37
print "difference btwn resolutions = "+str(diffr)

elapsed = (clock() -start)
print "Time elapsed = "+str(elapsed)
