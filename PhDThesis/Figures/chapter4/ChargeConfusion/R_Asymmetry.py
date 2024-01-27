import numpy as np
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show

def function(x1,x2):
    return (x1-x2)/(x1+x2)

#x1=np.arange(0,150,0.5)
#x2=np.arange(0,150,0.5)

x1=np.arange(-150,150,0.5)
x2=np.arange(-150,150,0.5)


X1,X2 = meshgrid(x1, x2) # grid of point
Z=function(X1,X2)

im = imshow(Z,cmap=cm.RdBu) # drawing the function
colorbar(im) # adding the colobar on the right
show()



