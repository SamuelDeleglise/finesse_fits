from pyinstruments.curvestore import models

from pylab import *



curve = models.CurveDB.objects.get(id=11963)
#figure(' large scan L18(uroc) L20(plan) ')
subplot(211)
curve.plot()
subplots_adjust(hspace=0.45, wspace=0.3)
xlabel('time (s)')
ylabel('Voltage (V)')


curve = models.CurveDB.objects.get(id=11955)
#figure(' FSR with large scan')
subplot(234)
curve.plot()
title(' TEM00')
xlabel('time (s)')
ylabel('Voltage (V)')

curve = models.CurveDB.objects.get(id=11958)
subplot(235)
curve.plot()
title('TEM01')
xlabel('time (s)')
ylabel('Voltage (V)')

curve = models.CurveDB.objects.get(id=11961)
subplot(236)
curve.plot()
title('TEM10')
xlabel('time (s)')
ylabel('Voltage (V)')





""""
plot(curve.data.index, curve.data)
"""

show()
