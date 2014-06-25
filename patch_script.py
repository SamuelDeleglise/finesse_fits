


for ch in c.childs.filter_param('name', value__contains='largescan'):
    fsr = finesse.FSRScan(ch.id)
    ch.params['length'] = fsr.length
    ch.params['finesse'] = fsr.mean_finesse
    ch.save()
for ch in c.childs.filter_param('name', value__contains='largescan'):
    lowscan = c.childs.filter_param('name', value__contains=ch.name.split('_')[1] + '_lowscan')[0]
    lowscan.params['length'] = ch.params['length']
    lowscan.save()
    print lowscan

for child in c.childs.filter_param('name', value__contains='lowscan'):
    query = child.childs.filter_param('name', value__contains='manualfit')
    if query.count()>=1:
        fit = query[0]
        print fit
    else:
        query = child.childs.filter_param('name', value__contains='portion')
        fit = query[0].childs.all()[0]
        print fit
    finesse_ = abs(get_finesse(fit.id, child.params['length']))
    print finesse_
    child.params['finesse'] = finesse_
    child.save()



for child in c.childs.filter_param('name', value__contains='lowscan'):
    query = child.childs.filter_param('name', value__contains='manualfit')
    if query.count()>=1:
        fit = query[0]
        print fit
    else:
        query = child.childs.filter_param('name', value__contains='portion')
        fit = query[0].childs.all()[0]
        print fit
    scale = abs(fit.params['scale'])
    print scale
    child.params['scale'] = scale
    child.save()

    
finesses = [curve.params['finesse'] for curve in c.childs.filter_param('name', value__contains='lowscan')]

lengths = [curve.params['length'] for curve in c.childs.filter_param('name', value__contains='lowscan')]

scales = [curve.params['scale'] for curve in c.childs.filter_param('name', value__contains='lowscan')]


from pylab import *
import numpy as np
figure()
subplot(211)
xlabel('length($\mu$m)')
ylabel('finesse')
plot(np.array(lengths)*1e6, finesses, 'o')
subplot(212)
ylabel('peak dip (arb)')
xlabel('length($\mu$m)')
plot(np.array(lengths)*1e6, scales, 'o')

