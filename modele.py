import math as m










def get_finesse(id):
    from pyinstruments.curvestore import models
    curve = models.CurveDB.objects.get(id=id)
    finesse = 3e8/(2*curve.params['length']*2*curve.params["bandwidth"]*250e6/curve.params["SBwidth"])
    curve.params['finesse'] = finesse
    curve.save()    
    return finesse




def radius(id):
    from pyinstruments.curvestore import models
    curve = models.CurveDB.objects.get(id=id)
    dnu_mean = (curve.params['dnu_low'] + curve.params['dnu_hi'])/2
    Rc = curve.params['length']/(1. - m.cos(m.pi*dnu_mean/curve.params['dnu_isl'])**2)
    return Rc
