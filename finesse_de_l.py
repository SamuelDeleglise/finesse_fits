import finesse

from pyinstruments.curvestore import models
from pylab import errorbar, show

c = models.CurveDB.objects.filter_tag("process/recuit")

finesses = []
finesses_error_low = []
finesses_error_up = []
lengths = []
for curve in c:
    if not curve.params['name']=='ramp':
        print curve.params['name']
        fsr = finesse.FSRScan(curve.params['id'])
        m1 = fsr.data_peaks.finesse.min()
        m2 = fsr.data_peaks.finesse.max()
        m3 = fsr.data_peaks.finesse.mean()
        print m1,m2,m3
        length = fsr.length
        finesses.append(m3)
        finesses_error_low.append(m3 - m1)
        finesses_error_up.append(m2 - m3)
        lengths.append(length)

errorbar(lengths, finesses, yerr=[finesses_error_low, finesses_error_up], fmt='.')
show()
