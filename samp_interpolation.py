import numpy as np
#import matplotlib as mp
#from matplotlib import pyplot as pp

class Interpolation:
    def __init__(self,zeta,hsamp_freq,vsamp_freq):
        self.kernel = {1:0.6213,2:-0.1704,3:0.0693,4:-0.0276,5:0.01,6:-0.003}
        self.zeta = zeta
        self.hsamp_freq = hsamp_freq
        self.vsamp_freq = vsamp_freq
        self.hmtf = self.calc_hmtf()
        self.vmtf = self.calc_vmtf()

    def calc_hmtf(self):
        kernel = self.kernel
        zeta = self.zeta
        sf = self.hsamp_freq*1000.0
        X = 1.0/sf
        res = 0.5
        for j in range(1,7):
            res = res+kernel[j]*np.cos(((2*j-1)*X*(1000*zeta))*3)
        return res

    def calc_vmtf(self):
        kernel = self.kernel
        zeta = self.zeta
        sf = self.vsamp_freq*1000.0
        X = 1.0/sf
        res = 0.5
        for j in range(1,7):
            res = res+kernel[j]*np.cos(((2*j-1)*X*(1000*zeta))*3)
        return res


