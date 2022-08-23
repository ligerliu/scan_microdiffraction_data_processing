import numpy as np
from lmfit import Model

def xy_pos_cal(deg,xc,yc,th,r):
    #x = xc + r*np.cos(np.radians(deg+th))
    y = yc + r*np.sin(np.radians(deg+th))
    return y#np.append(x,y)

fmod = Model(xy_pos_cal)

fmod.set_param_hint('xc',value=-10.5)#,min=-12,max=-10)
fmod.set_param_hint('yc',value=-66)#,min=-68,max=-64)
fmod.set_param_hint('th',value=0)#,min=0,max=360)
fmod.set_param_hint('r',value=0.5)#,min=0.001,max=1)

params = fmod.make_params()

def test():
    x = np.array([-60.53,-60.48,-60.43,-60.40,-60.34])
    y = np.array([-11.6345,-11.5880,-11.5345,-11.458,-11.3945])
    deg = np.array([90,80,70,60,50])
    #res = fmod.fit(np.append(x,y),params,deg=deg)
    res = fmod.fit(y,params,deg=deg)
    print(res.params)
    print(f'\n{res.best_fit}')
