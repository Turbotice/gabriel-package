import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def animation_pcolormesh_gld(data):
    #%matplotlib qt
    fig = plt.figure()
    im = plt.imshow(data[0,...], animated=True)
    def updatefig(*args):
        global i
        if (i<len(data[:,1,1])):
            i += 1
        else:
            i=0
        im.set_array(data[i,...])
        return im,
    ani = animation.FuncAnimation(fig, updatefig,  blit=True)
    plt.show()
