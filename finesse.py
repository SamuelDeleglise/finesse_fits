from pyinstruments.curvestore.curve_create_widget import CurveCreateWidget
from pyinstruments.pyhardwaredb import instrument
from pyinstruments.curvestore import models
from curve.fitting import Fit

from PyQt4 import QtCore, QtGui
import pylab
import json
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy.interpolate import splrep, splev
import numpy as np
from pylab import subplot, plot

c = 2.99e8

class FinesseMeasurementWidget(CurveCreateWidget):
    def __init__(self, parent=None):
        super(FinesseMeasurementWidget, self).__init__(**self._get_defaults())
        self.save_button.setText('acquire')
        self.save_button.clicked.connect(self.acquire)
        self.curve_modified.connect(self._save_defaults)
        self.show()
    
    
    def _get_initial_defaults(self):
        return {"default_name":'multi_FSR',
                 "default_window":'default_win',
                 "tags":[],
                 "comment":""}
        
    def _get_defaults(self):
        settings = QtCore.QSettings("pyinstruments", "pyinstruments")
        kwds_str = str(settings.value("default_session").toString())
        if kwds_str != "":
            kwds = json.loads(kwds_str)
        else:
            kwds = self._get_initial_defaults()
        return kwds

    def _save_defaults(self):
        kwds = {"default_name":self.name,
                 "default_window":self.window,
                 "tags":self.tags,
                 "comment":self.comment}
        settings = QtCore.QSettings("pyinstruments", "pyinstruments")
        settings.setValue("default_session", json.dumps(kwds))

    def acquire(self):
        scope = instrument('IfoScope')
        print 'acquiring signal (channel 1)'
        scope.channel_idx = 1
        curve = scope.get_curve()
        curve.clever_downsample(100)
        self.apply_values_to_curve(curve)
        curve.params["name"] = curve.params["name"]
        curve.save()

        print 'acquiring ramp (channel 2)'
        scope.channel_idx = 2
        curve2 = scope.get_curve()
        curve2.downsample()
        self.apply_values_to_curve(curve2)
        curve2.params["name"] = "ramp"
        curve.add_child(curve2)
        #curve2.save()
        
FINESSE_WIDGET = FinesseMeasurementWidget()

class DataPeaks(object):
    def __init__(self, parent):
        self.parent = parent
        self.curve = parent.curve
        
        self._times = None
        self._fit_curves = None
        self._dfdt_abs = None
        self._kappa_over_2pi = None
        self._slopes = None
        self._bandwidth = None

    def make_portions(self, threshold=0.5, time=35e-6):
        """
        All peaks that are stronger than threshold are saved in a separated child curve
        """
        
        peaks = self.child_peaks
        if peaks.count()!=0:
            yn = raw_input("delete childs ?")
            if yn=="y" or yn=="Y":
                for child in peaks:
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

    @property
    def child_peaks(self):
        return self.curve.childs.filter(_name__contains='peak@')
    
    def fit_peaks(self):
        peaks = self.curve.childs.filter(_name__contains='peak@')
        dfdt = []
        x0 = []
        kappa = []
        times = []
        for child in peaks:
            child = child.fit("lorentzSB", autosave=True)


    
    @property
    def peaks(self):
        res = self.curve.childs.filter(_name__contains="peak@")
        if res.count()==0:
            print "Partitioning curve into peaks..."
            self.make_portions()
            print "Fitting individual peaks..."
            self.fit_peaks()
            return self.curve.childs.filter(_name__contains="peak@")
        return res

    @property
    def fit_curves(self):
        if self._fit_curves is None:
            fits = []
            for curve in self.peaks:
                if curve.childs.all().count()!=1:
                    raise ValueError("Curve " + curve + " should have exactly 1 child!")
                else:
                    fits.append(curve.childs.all()[0])
            self._fits = fits
        return self._fits

    @property
    def times(self):
        if self._times is None:
            self._times = np.array([f.params['x0'] for f in self.fit_curves])
        return self._times
    
    @property
    def dfdt_abs(self):
        if self._dfdt_abs is None:
            self._dfdt_abs = np.array(
            [self.parent.mod_freq/f.params["SBwidth"] for f in self.fit_curves])
        return self._dfdt_abs
    
    @property
    def bandwidth_in_s(self):
        if self._bandwidth is None:
            self._bandwidth = np.array([np.abs(f.params['bandwidth']) for f in self.fit_curves])
        return self._bandwidth
    
    @property
    def kappa_over_2pi(self):
        return 2*self.bandwidth_in_s*self.dfdt_abs
    
    def sign_of_slope(self):
        from scipy import signal
        b, a = signal.butter(2, 0.005)
        self.voltage_of_t_filtered = signal.filtfilt(b, a, self.parent.ramp.data, padlen=150)
        dt = self.parent.ramp.data.index[1] - self.parent.ramp.data.index[0]
        inter_diff = interp1d(np.array(self.parent.ramp.data.index[1:], dtype=float), np.diff(self.voltage_of_t_filtered)/dt)
        return np.array(np.sign(inter_diff(self.times))).astype(int)
    
    @property
    def slopes(self):
        if self._slopes is None:
            if 'slope' in self.peaks[0].params:
                self._slopes = [p.params['slope'] for p in self.peaks]
            else:
                self._slopes = self.sign_of_slope()
                for curve, slope in zip(self.peaks, self._slopes):
                    curve.params['slope'] = slope
                    curve.save()
                    print 'slope saved', slope
        return np.array(self._slopes)
    
    @property
    def dfdt(self):
        return self.dfdt_abs*self.slopes
    
    @property
    def mode_number(self):
        current_mode = 0
        list_mode = []
        for signe in self.slopes:
            #np.sign(np.sin(self.phi + 2*np.pi*self.pzt_freq*x))
            if signe>0:
                current_mode+=signe
            list_mode.append(current_mode)
            if signe<0:
                current_mode+=signe
        return abs(np.array(list_mode))
    
    @property
    def finesse(self):
        return 2.99e8/(2*self.parent.length*self.kappa_over_2pi)
    
        
class ModelPeaks(object):
    def __init__(self, parent):
        self.parent = parent
        self.pzt_freq = parent.pzt_freq
        self._fourier_coeffs = None
    
    
    @property
    def nice_x(self):
        start = np.min(np.array(self.parent.curve.data.index))
        stop = np.max(np.array(self.parent.curve.data.index))
        return np.linspace(start, stop, 500)
    
    @property
    def peak_times(self):
        return self.parent.data_peaks.times
    
    def _ampl_phase(self, coeffs):
        """
        Separate the fourier coeffs into amplitudes and phases coefficients
        """
        return  coeffs[:len(coeffs)/2], coeffs[len(coeffs)/2:]
    
    def dfdt(self, coeffs, times):
        y = np.zeros(len(times))
        amplitudes, phases = self._ampl_phase(coeffs)
        for n, (ampl, phase) in enumerate(zip(amplitudes, phases)):
            omega_n = 2*np.pi*self.pzt_freq*(n+1)
            y-=omega_n*ampl*np.sin(omega_n*times + phase)
        return y
    
    def fit_fourier_coeffs(self):
        def error_dfdt(coeffs):
            return self.dfdt(coeffs, self.peak_times) - self.parent.data_peaks.dfdt
        res = leastsq(error_dfdt, [1e14,1e14,1.,1.])
        return res[0]
    
    @property
    def fourier_coeffs(self):
        """
        These are the coefficients of the fourier series for pzt_position(t)
        """
        
        if self._fourier_coeffs is None:
            self._fourier_coeffs = self.fit_fourier_coeffs()
        return self._fourier_coeffs
    
    def freq_of_t(self, times):
        """
        returns the infered resonance frequency as a function of time
        """
        y = np.zeros(len(times))
        amplitudes, phases = self._ampl_phase(self.fourier_coeffs)
        for n, (ampl, phase) in enumerate(zip(amplitudes, phases)):
            omega_n = 2*np.pi*self.pzt_freq*(n+1)
            y+=ampl*np.cos(omega_n*times + phase)
        return y
        
        
        
class FSRScan(object):
    """
    Extract the finesse, and length of a cavity from an FSR scan
    """
    
    def __init__(self, curve_id, pzt_freq=25, mod_freq=250e6):
        self.curve = models.CurveDB.objects.get(id=curve_id)
        self.ramp  = self.curve.childs.get(_name='ramp')
        self.pzt_freq = pzt_freq
        self.mod_freq = mod_freq
        self._length = None
        self._offset = None
        
        self.data_peaks = DataPeaks(self)
        self.model_peaks = ModelPeaks(self)
    
    def clear_all(self):
        for child in self.data_peaks.child_peaks:
            child.delete()
            self.data_peaks = DataPeaks(self)
            self.model_peaks = ModelPeaks(self)
    
    def linear_freq(self, args, mode_number):
        length, offset = args
        mode_number = np.array(mode_number)
        return c/(2*length)*mode_number + offset 

        
    def fit_length(self):
        """
        def accumulated_freq(L):
            #times = np.array(self.x0)
            #return self.amplitude*np.cos(self.phi + 2*np.pi*self.pzt_freq*times) - self.amplitude*np.cos(self.phi + 2*np.pi*self.pzt_freq*times[0])
            return self.res_freq_of_t()/1064e-9*np.pi/(2*L)
        """
                 
        def error_freq_of_t(args):
            return self.linear_freq(args, self.data_peaks.mode_number) \
                - np.array(self.model_peaks.freq_of_t(self.data_peaks.times))
    

        res = leastsq(error_freq_of_t, [200e-6, -1e11])
        self._length, self._offset = res[0]
        
    @property
    def length(self):
        if self._length is None:
            self.fit_length()
        return self._length
    @property
    def offset(self):
        if self._offset is None:
            self.fit_length()
        return self._offset
    
    def plot_summary(self):
        pylab.suptitle('Finesse/Length analysis of curve #' \
                    + str(self.curve.id) + ' with ramp #' + str(self.ramp.id))
        subplot(311)
        pylab.xlabel('time[s]')
        pylab.ylabel('df/dt[Hz/s]')
        plot(self.data_peaks.times, self.data_peaks.dfdt_abs, 'o', label='measured')
        plot(self.data_peaks.times, self.data_peaks.dfdt, '.', label='infered')
        plot(self.model_peaks.nice_x, 
             self.model_peaks.dfdt(self.model_peaks.fourier_coeffs, 
                                   self.model_peaks.nice_x),
            label='model')
        pylab.legend(loc='best')
    
        subplot(312)
        pylab.xlabel('mode number')
        pylab.ylabel('freq[Hz]')
        plot(self.data_peaks.mode_number,
             self.model_peaks.freq_of_t(self.data_peaks.times),
             '+', 
             label='mode freq. (model)')
        
        abcisse = np.linspace(min(self.data_peaks.mode_number)-0.5,
                                       max(self.data_peaks.mode_number)+0.5)
        plot(abcisse, 
             self.linear_freq([self.length, self.offset], abcisse),
             label='mode freq. for L=' + "{0:.1f}".format(self.length*1e6) + 'um)')
        pylab.legend(loc='best')
    
        subplot(313)
        pylab.xlabel('time[s]')
        pylab.ylabel('kappa/2pi[Hz]')
        pylab.plot(self.data_peaks.times,
                   self.data_peaks.kappa_over_2pi,
                   'o',
                   label='measured linewidth')
        mean_kappa = np.mean(self.data_peaks.kappa_over_2pi)
        pylab.plot([min(self.data_peaks.times), max(self.data_peaks.times)],
                   [mean_kappa, mean_kappa],
                   label="%.1f"%(mean_kappa*1e-6) + ' MHz -> F='+ \
                   "{0:.0f}".format(self.mean_finesse))
        pylab.legend(loc='best')
        pylab.show()
        
    @property
    def mean_finesse(self):
        return np.mean(self.data_peaks.finesse)
    

        
