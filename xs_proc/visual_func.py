import numpy as np
from skimage.draw import polygon
import matplotlib.pyplot as plt
import h5py

def poly_vertex_draw():
    # choose points from plot, maximum number of point is 100
    # choose need to finished within 120s
    coord = plt.ginput(100,timeout=120)
    # change all the coordinate of point into integer
    coord = np.array(coord,dtype=int)
    return coord

def draw_polygon(args,shape=None):
    if shape is None:
        print("shape info needed")
        return
    l = list(args)
    if len(l) < 3:
        print("polygon has at least three vertex")
    r = []
    c = []
    for _ in l:
        #here r correlates to the x coordinate, c correlates to y
        r.append(_[1])
        c.append(_[0])
    return polygon(r,c,shape)

def polygon_with_int_thrshd(int_map,
                            vertex,
                            low_int_thrshd = 0,
                            high_int_thrshd=np.inf,
                            ):
    #if input contains multiple vertex array, vertex = [vertex1,vertex2,...]
    shape = int_map.shape
    #m1 = (int_map>=low_int_thrshd)
    #m1[m1>high_int_thrshd] = False
    m1 = (int_map<=low_int_thrshd)
    m1[int_map>=high_int_thrshd] = True
    if not isinstance(vertex,type(None)):
        if not isinstance(vertex[0],list):
            vertex = [vertex]
        for _ in range(len(vertex)):
            p = draw_polygon(vertex[_],shape=shape)
            if _ == 0:
                m2 = np.zeros(shape).astype(bool)
            m2[p] = True
    #m = (m1* m2).astype(bool)
    m = m1+m2
    idxx,idxy = np.meshgrid(np.arange(shape[1]),np.arange(shape[0]))
    idx_total = np.arange(shape[0]*shape[1]).reshape(shape)
    return ((idxx[m],idxy[m]),idx_total[m],m)
    
def mask_making(int_map,
                vertex=None,
                low_int_thrshd = 0,
                high_int_thrshd=np.inf,
                            ):
    #if input contains multiple vertex array, vertex = [vertex1,vertex2,...]
    shape = int_map.shape
    #m1 = (int_map>=low_int_thrshd)
    #m1[m1>high_int_thrshd] = False
    m1 = (int_map<=low_int_thrshd)
    m1[int_map>=high_int_thrshd] = True
    if not isinstance(vertex,type(None)):
        # if not isinstance(vertex[0],list):
        #    vertex = [vertex]
        # the definition of vertex[0] is confusing the program to work correctly
        #if len(vertex[0]) < 3:
        #    vertex = [vertex] 
        for _ in range(len(vertex)):
            p = draw_polygon(vertex[_],shape=shape)
            if _ == 0:
                m2 = np.zeros(shape).astype(bool)
            m2[p] = True
    #m = (m1* m2).astype(bool)
    m = m1+m2
    return m
