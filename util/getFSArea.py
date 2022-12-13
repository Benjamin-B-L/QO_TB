import matplotlib.pyplot as plt
import numpy as np 

########LINEBUILDER CLASS
class LineBuilder:
    """
    Object that allows to build lines point per point by clicking with a mouse on the plot. At each click ON the plot,
    a data point is added to the current line. A new line is created when clicking outside of the plot. 
    Once a new line has been created, the previous one can not be modified again. 
    """
    def __init__(self, line,params):
        """
        Simple initialization, create with an empty line. 
        """
        self.nlines = 1
        self.line = line
        self.xs = {self.nlines: list(line[self.nlines].get_xdata())}
        self.ys = {self.nlines: list(line[self.nlines].get_ydata())}
        self.cid = line[self.nlines].figure.canvas.mpl_connect('button_press_event', self)
        self.params = params
        print("\nLine builder initialized, we start with first line")

    def __call__(self, event):
        """
        Event control function: when mouse clicking this function is called. If the click is outside
        of the plot, create a new line, otherwise append the current line and replot. 
        """

        ##if outside the plot
        if event.inaxes!=self.line[self.nlines].axes:
            self.add_line() 
            return

        ##if inside, read position
        xmouse, ymouse = event.xdata, event.ydata
        ##append the current line
        self.xs[self.nlines].append(xmouse)
        self.ys[self.nlines].append(ymouse)
        self.line[self.nlines].set_data(self.xs[self.nlines], self.ys[self.nlines])
        ##replot
        self.line[self.nlines].figure.canvas.draw()
        print("new point in line at ("+str(xmouse)+","+str(ymouse)+")")

    def add_line(self):
        """
        Add a new line to the dictionary: we increse nlines by 1 which defines the index of the line. 
        The new created line starts empy. 
        """
        print("\nAdding a new line!  Total nbr of lines=",self.nlines+1)
        self.nlines +=1
        self.line[self.nlines], = ax.plot([], [], ls=self.params.ls,lw=self.params.lw,marker=self.params.marker,color=self.params.color)
        self.xs[self.nlines] = list(self.line[self.nlines].get_xdata())
        self.ys[self.nlines] = list(self.line[self.nlines].get_ydata())
        self.cid = self.line[self.nlines].figure.canvas.mpl_connect('button_press_event', self)

###########LINEPARAMS CLASS
class LineParams:
    """
    A simple class that allows to control the display of the created lines. 
    """

    def __init__(self, ls='-', lw = 2, marker='+', color='magenta'):
        """
        Simple initialization, create with an empty line. 
        """
        self.ls=ls
        self.lw = lw
        self.marker=marker
        self.color=color
        print("Line parameters are defined as:\n\
Linestyle = "+ls+", Linewidth = "+str(lw)+", marker = "+marker+" and color = "+color+".")

###################################################################################################################
####  FUNCTIONS
###################################################################################################################

#### AREA CALCULATOR ####
def getAreas(line_collection, x, y):
    """
    Compute the area defined by each line, and the total area of the BZ. 
    WARNING: assumes RECTANGULAR shape of the BZ. For other shapes, just use a standard
    line to get the full BZ area. 
    """
    #First total area
    total_area = (x[len(x)-1]-x[0])*(y[len(y)-1]-y[0])
    
    #Then areas for each line using shoelace algorithm
    areas = np.zeros(line_collection.nlines,dtype=float)
    for iline in range(line_collection.nlines):
        xpoly = line_collection.xs[iline+1]
        ypoly = line_collection.ys[iline+1]
        areas[iline] = 0.5*np.abs(np.dot(xpoly,np.roll(ypoly,1))-np.dot(ypoly,np.roll(xpoly,1)))
    
    return areas,total_area

#### READ FS from calculation ####
def read_FSfromfile(filename, Nkx, Nky):
    """
    Read FS file from QO_TB calculation, and put it in format suitable for this small code. 
    """

    #open file
    f = open(filename,'r')
    dat = f.readlines()
    ncols = len(dat[1].strip().split())-2
    #define output arrays
    kx = np.zeros(Nkx,dtype=float)
    ky = np.zeros(Nky,dtype=float)
    FS = np.zeros([Nkx,Nky,ncols],dtype=float)

    ##Read data now
    for ikx in range(Nkx):
        #read ky
        if ikx==0:
            for iky in range(Nky):
                ky[iky] = float(dat[1+iky].strip().split()[1])
        #read kx
        kx[ikx] = float(dat[1+ikx*Nkx].strip().split()[0])
        #read FS
        for iky in range(Nky):
            for icol in range(ncols):
                FS[ikx,iky,icol] = float(dat[1+ikx*Nkx+iky].strip().split()[2+icol])

    f.close()
    return kx, ky, FS



######### PARAMETERS ##########
filename='myCALC.FS'
Nkx = 200
Nky = 200
kx, ky, FS = read_FSfromfile(filename,Nkx,Nky)
FStoplot = FS[:,:,0]                              ## here restrict to the column you want

##Initialize line parameters
lineparams = LineParams()


######### Line building ########
fig, ax = plt.subplots()
tolerance = 0.01                    ###How large is the area considered around the click
im = ax.pcolormesh(kx, ky, np.transpose(FStoplot), shading='auto', rasterized='True',vmin=0,vmax=3,picker=tolerance)
line_collection = LineBuilder({1: ax.plot([], [],ls=lineparams.ls,lw=lineparams.lw,marker=lineparams.marker,color=lineparams.color)[0]},lineparams)
plt.show()

######### Calculate all areas and replot ########
areas, total_area = getAreas(line_collection,kx,ky)               ##get all areas, including that of the BZ

##print the relative areas
print("\nRelative areas for each line drawn (in %):")
cnt = 1
for area in areas:
    print("  - line "+str(cnt)+": "+str(100*(area/total_area)))

##replot all the lines on top of FS
fig, ax = plt.subplots()
tolerance = 1                      ###How large is the area considered around the click
im = ax.pcolormesh(kx, ky, np.transpose(FStoplot), shading='auto', rasterized='True',vmin=0,vmax=3,picker=tolerance)
for iline in range(line_collection.nlines):
    x = np.concatenate([line_collection.xs[iline+1],[line_collection.xs[iline+1][0]]])    ##do this to close the loops
    y = np.concatenate([line_collection.ys[iline+1],[line_collection.ys[iline+1][0]]])
    ax.plot(x,y,ls=lineparams.ls,lw=lineparams.lw,marker=lineparams.marker,color=lineparams.color)
plt.show()