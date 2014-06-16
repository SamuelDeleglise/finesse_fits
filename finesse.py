
from pyinstruments.pyhardwaredb import instrument
from pyinstruments.curvestore import models
from curve.fitting import Fit

import pylab
import json
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy.interpolate import splrep, splev
import numpy as np

c = 3e8

def full_analysis(curve_id, ramp_id):
    fsr = FSRScan(curve_id, ramp_id)
    fsr.make_portions()
    fsr.fit_peaks()
    fsr.plot_fits()
    fsr.fit_dfdt_of_t()
    fsr.fit_length()
    fsr.make_summary()
    return fsr

def fit_rebonds_sb(curve, mod_freq_mhz=250):
    """
    From a curve containing 2 sidebands and a rebond,
    guess the scan speed from the sidebands and fit tau
    on the rebonds with fixed df_dt.
    A child curve containing only the central rebonds is created
    """
    
    fit = Fit(curve.data, 'lorentzSB')
    guess = fit._guesslorentzSB()
    
    child = models.CurveDB.create(curve.data.index, fit.lorentzSB(**guess), name='guess_lorentz_sb')
    curve.add_child(child)
    
    lower = guess['x0'] - guess['SBwidth']/2
    upper = guess['x0'] + guess['SBwidth']/2
    new_dat = curve.data[lower:upper]
    child = models.CurveDB.create(new_dat, name='rebonds_only')
    curve.add_child(child)
    dfdt_mhz_per_us = mod_freq_mhz/(guess['SBwidth']*1e6)
    
    fitcurve = child.fit('rebonds_reflexion', 
                         autosave=True,
                         manualguess_params={'tcav_us':1./(2*np.pi*dfdt_mhz_per_us*2*guess['bandwidth']*1e6)},  fixed_params={'dfdt_mhz_per_us':dfdt_mhz_per_us}, graphicalfit=True)
   
    
    print mod_freq_mhz/(guess['SBwidth']*1e6)
    curve.params['dfdt_mhz_per_us'] = fitcurve.params['dfdt_mhz_per_us']
    curve.params['tcav_us'] = fitcurve.params['tcav_us']
    curve.params['x0'] = fitcurve.params['x0']
    

class FSRScan(object):
    """
    Extract the finesse, and length of a cavity from an FSR scan
    """
    def __init__(self, curve_id, curve_id_ramp, mod_freq=250e6):
        self.curve = models.CurveDB.objects.get(id=curve_id)
        self.ramp  = models.CurveDB.objects.get(id=curve_id_ramp)
        curve = self.curve
        self.mod_freq = mod_freq
        if 'summary_json' in curve.params:
            dic = json.loads(curve.params["summary_json"])
            self.x0 = dic['x0']
            self.dfdt = dic['dfdt']
            self.kappa = dic['kappa']
            if 'amplitude' in dic:
                #self.phi = dic['phi']
                self.amplitude = dic['amplitude']
                self.pzt_freq = dic['pzt_freq']
                
    def acquire(self):
        scope = instrument('IfoScope')
        curve = scope.get_curve()
        self.curve = curve.clever_downsample()
        self.curve.save()
        
    
    def make_portions(self, threshold=0.5, time=35e-6):
        """
        All peaks that are stronger than threshold are saved in a separated child curve
        """
        
        if self.curve.has_childs:
            yn = raw_input("delete childs ?")
            if yn=="y" or yn=="Y":
                for child in self.curve.childs.all():
                    child.delete()
        portions = []
        data = self.curve.data
        mag = (data - data.mean()).abs()
        threshold = threshold*mag.max()
        next = mag.index[0]
        for x, val in zip(mag.index, mag):
            if val>threshold:
                if x>next:
                    portion = data[x - time:x + time]
                    portions.append(portion)
                    next = x + time
                    child = models.CurveDB.create(portion.index, portion, name="peak@"+str(x*1e6)+"us")
                    self.curve.add_child(child)
        
    def fit_peaks(self):
        peaks = self.curve.childs.filter(_name__contains='peak@')
        dfdt = []
        x0 = []
        kappa = []
        times = []
        for child in peaks:
            #fit_rebonds_sb(child)
            child = child.fit("lorentzSB", autosave=True, maxiter=10000)
            temp = self.mod_freq/child.params['SBwidth']
            dfdt.append(temp)
            x0.append(child.params['x0'])
            kappa.append(2*abs(child.params['bandwidth'])*temp)

        self.dfdt = dfdt
        self.x0 = x0
        self.kappa = kappa
        self.curve.params["summary_json"] = json.dumps({'x0':x0,'dfdt':dfdt,'kappa':kappa})
        self.curve.save()
        return dfdt, x0, kappa
 
 
    def plot_fits(self):
        pylab.subplot(412)
        pylab.title('dfdt')
        pylab.plot(self.x0, self.dfdt,'*')
        #pylab.subplot(413)
        #pylab.title('FSR')
        #pylab.plot(self.x0, range(len(self.x0)), '*')
        pylab.subplot(414)
        pylab.title('kappa')
        pylab.plot(self.x0, self.kappa,'*')


    def direction_pzt(self):
        times = np.array(self.times)
        inter_diff = interp1d(np.array(self.ramp.data.index[1:], dtype=float), np.diff(self.res_freq_of_t_filtered))
        return np.sign(inter_diff(times))

    def dfdt_of_t(self, args, times):
        from scipy import signal
        b, a = signal.butter(2, 0.005)
        [amplitude,] = args
        self.res_freq_of_t_filtered = signal.filtfilt(b, a, amplitude*self.ramp.data, padlen=150)
        dt = self.ramp.data.index[1] - self.ramp.data.index[0]
        inter_diff = interp1d(np.array(self.ramp.data.index[1:], dtype=float), np.diff(self.res_freq_of_t_filtered)/dt)
        return np.abs(inter_diff(times))
        #[amplitude, phi] = args
        #return np.abs(pzt_freq*amplitude*np.sin(phi + 2*np.pi*pzt_freq*times))
    
    def res_freq_of_t(self, times):
        #inter = interp1d(np.array(self.ramp.data.index,dtype=float), self.res_freq_of_t_filtered)
        inter = interp1d(np.array(self.ramp.data.index,dtype=float), self.res_freq_of_t_filtered)
        return inter(times)
    
    def error_function_dfdt(self, args, times):
        return self.dfdt_of_t(args, times) - np.array(self.dfdt)


    def fit_dfdt_of_t(self, pzt_freq=25.):
        times = np.array(self.x0)
        dfdt = np.array(self.dfdt)
        
        res = leastsq(self.error_function_dfdt, [1e11],self.x0)
        
        self.pzt_freq = pzt_freq
        self.amplitude, = res[0]
        
        pylab.subplot(411)
        times = np.linspace(min(times), max(times),1000)
        pylab.plot(times, self.res_freq_of_t(times))
        
        pylab.subplot(412)
        
        
        pylab.plot(times, self.dfdt_of_t(res[0],times))
        self.curve.params["summary_json"] = json.dumps({'x0':self.x0,'dfdt':self.dfdt,'kappa':self.kappa, 'amplitude':self.amplitude, 'pzt_freq':self.pzt_freq})
        self.curve.save()
        
        
    def mode_number(self):
        current_mode = 0
        list_mode = []
        self.times = self.x0
        for signe in self.direction_pzt():
            #np.sign(np.sin(self.phi + 2*np.pi*self.pzt_freq*x))
            if signe>0:
                current_mode+=signe
            list_mode.append(current_mode)
            if signe<0:
                current_mode+=signe
        return abs(np.array(list_mode))

    def linear_freq(self, args, mode_number):
        L, offset = args
        mode_number = np.array(mode_number)
        return c/(2*L)*mode_number + offset 

    def error_freq_of_t(self, args, times):
        times = np.array(times)
        return self.linear_freq(args, self.mode_number()) - np.array(self.res_freq_of_t(times))
        
    def fit_length(self):
        """
        def accumulated_freq(L):
            #times = np.array(self.x0)
            #return self.amplitude*np.cos(self.phi + 2*np.pi*self.pzt_freq*times) - self.amplitude*np.cos(self.phi + 2*np.pi*self.pzt_freq*times[0])
            return self.res_freq_of_t()/1064e-9*np.pi/(2*L)
        """
                    
        pylab.subplot(413)
        pylab.plot(self.mode_number(), self.res_freq_of_t(self.x0), '*') #

        res = leastsq(self.error_freq_of_t, [200e-6, -1e11], self.x0)
        self.L = res[0][0]
        self.offset = res[0][1]
        pylab.subplot(413)
        x = np.linspace(0, max(self.mode_number()))
        
        pylab.plot(x, self.linear_freq([self.L, self.offset], x))
    
    def make_summary(self):
        """
        
        """
        
        self.mean_kappa = np.mean(self.kappa)
        self.mean_finesse = np.pi*3e8/(self.L*self.mean_kappa*2*np.pi)
        
        print """=====================SUMMARY========================"""
        print "Length: ", self.L
        print "Kappa/2pi: ", np.mean(self.kappa)
        print "Mean finesse: ", self.mean_finesse
        self.curve.params["length"] = self.L
        self.curve.params["kappa_over_2pi"] = np.mean(self.kappa)
        self.curve.params["finesse"] = self.mean_finesse
        self.curve.save()
        print """====================================================="""
