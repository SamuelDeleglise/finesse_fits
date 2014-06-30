import finesse
from pyinstruments.curvestore import models
from modele import get_finesse

def plot_finesse(id):
    curve = models.CurveDB.objects.get(id=id)
    test = curve.childs.filter_param('name', value__contains='lowscan')[0]
    try:
        test.params['finesse']
        test.params['length']
    except KeyError:
        for ch in curve.childs.filter_param('name', value__contains='largescan'):
            fsr = finesse.FSRScan(ch.id)
            fsr.save_summary()
            #ch.params['length'] = fsr.length
            #ch.params['finesse'] = fsr.mean_finesse
            #ch.save()

        for ch in curve.childs.filter_param('name', value__contains='largescan'):
            lowscans = curve.childs.filter_param('name', value__contains=ch.name.split('_')[1] + '_lowscan')
            for lowscan in lowscans.all():
                lowscan.params['length'] = ch.params['length']
                lowscan.save()
                print lowscan

        
        for child in curve.childs.filter_param('name', value__contains='lowscan'):
            query = child.childs.filter_param('name', value__contains='manualfit')
            if query.count()>=1:
                fit = query[0]
                print fit
            else:
                query = child.childs.filter_param('name', value__contains='portion')
                fit = query[0].childs.all()[0]
                print fit
            finesse_ = abs(get_finesse(fit.id))
            print finesse_
            child.params['finesse'] = finesse_
            child.save()

    finesses = [c.params['finesse'] for c in curve.childs.filter_param('name', value__contains='lowscan')]

    lengths = [c.params['length'] for c in curve.childs.filter_param('name', value__contains='lowscan')]


    
    from pylab import *
    import numpy as np
    figure('length scan')
    xlabel('length($\mu$m)')
    ylabel('finesse')
    plot(np.array(lengths)*1e6, finesses, 'o')
   
