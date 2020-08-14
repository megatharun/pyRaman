
''' 
Copyright (c) 2018, Megat Harun Al Rashid bin Megat Ahmad and Wilfred@Sylvester Paulus
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of conditions
    and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
    and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
    or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
'''

# Import all the necessary libraries
import wdfReader as wdf
import numpy as np
from scipy import misc
import matplotlib.pyplot as plt
import imageio
from scipy.optimize import curve_fit
from bokeh.plotting import figure, output_file, output_notebook, show

colorlist = ['blue','green','red', 'aqua', 'magenta', 'yellow', 'black',\
             'chocolate','violet', 'purple', 'olive']

output_notebook()

#-------------------------------------general functions--------------------------------------

def iplot():
        
    global plot1
    global counter
    counter = 0
    labelStr = '\u0394\u03C9, cm\u207b\u00B9'
    
    plot1 = figure(plot_width = 600, plot_height = 400, \
                   x_axis_label = labelStr, y_axis_label = 'Counts')
    
def showplot(displayhtml = 'no', htmlfile = "output.html"):
   
    global plot1
    
    if displayhtml == 'yes':
        output_file(htmlfile)
        
    plot1.legend.location = "top_left"
    plot1.legend.background_fill_alpha = 0.5

  
    show(plot1)
    
def fitLorentz(xVar,IVar, gammaVar, x0Var, mVar,cVar):
    
    return IVar*(gammaVar**2/((xVar-x0Var)**2+gammaVar**2)) + (mVar*xVar + cVar)

#-------------------------------------------------------------------------------------------

class scandata:
    
    def __init__(self, filename, size = [20,20], spectrumlabel = 'spectrum 1'):
    
        # Reading .wdf files using wdfReader
        fileContent = wdf.wdfReader(filename)
        # the x array is flipped because it starts from end
        self.arrayX = np.flip(fileContent.get_xdata(),0); self.arrayY = fileContent.get_spectra()
        
        # Pixel dimension
        self.size = size
        
        # Array shape of data extracted
        self.arrayXshape = self.arrayX.shape
        self.arrayYshape = self.arrayY.shape
        
        self.spectrumlabel = spectrumlabel
        
        # Flipping (because the array start from end
        # and reshaping the array to select spectrum at specific pixel coordinate
        self.arrayYP = np.flip(self.arrayY.reshape(size[0],size[1],self.arrayXshape[0]),2)
        self.totalNumSpectrum = int(self.arrayYshape[0]/self.arrayXshape[0])
        
        # Number of pixels
        print ('Total number of spectrums :', self.totalNumSpectrum)
        print ('Pixel area :', self.size[0],'x',self.size[1])
        print ('Total number of x-axis data :', self.arrayXshape[0])
        
    def spectrumplot(self, spot = [0,0], ylabel = 'Counts'):
        
        global counter
        global colorlist
        global plot1
        
        # By default the spectrum at coordinate (0,0) is shown
        self.yData = self.arrayYP[spot[0],spot[1],:]
        
        plot1.line(x=self.arrayX, y=self.yData, color = colorlist[counter], \
                   line_width=2, alpha=0.50, legend = self.spectrumlabel+str(spot))
            
        if counter == 10:
            counter = 0
        else:
            counter = counter + 1
            
    def checkfitting(self, spot = [0,0], xrange = [750,800], \
                     IntLor = 800.0, mU = 777.0, gammaL = 2.5, num=250):
                
        global counter
        global colorlist
        global plot1
        
        labelUnit = 'cm\u207b\u00B9'
        
        # Setting the range based on raman wavelength
        limMin = ((self.arrayX <= xrange[0])*1).sum()
        limMax = ((self.arrayX <= xrange[1])*1).sum()
        
        xData = self.arrayX[limMin:limMax]
        yData = self.arrayYP[spot[0],spot[1],limMin:limMax]
        
        # Creating starting values for background slope and constant (based on linear equation)
        mVal = (yData[-1]-yData[0])/(xData[-1]-xData[0])
        cVal = yData[0] - mVal*xData[0]
        
        dataX = np.linspace(xrange[0], xrange[1], num) # Creating array for plotting fitting
        Lor = IntLor*(gammaL/((dataX-mU)**2 + gammaL)) # Lorentz equation
        Lin = mVal*dataX + cVal # Background linera equation
        dataY =  Lor + Lin # Total fitting
        
        # -------------------------Ploting------------------------------------
        plot1.circle(x = xData, y = yData, \
                     fill_color = colorlist[counter], fill_alpha=0.5, alpha=0.5, size = 8,
                    legend = self.spectrumlabel+str(spot))
        plot1.patch(x = dataX, y = dataY, \
                   color = colorlist[counter], line_width=2, alpha=0.25, \
                    legend = 'Lorentz '+self.spectrumlabel+str(spot))
        
        if counter == 10:
            counter = 0
        else:
            counter = counter + 1

    def functionfitting(self,spot = [0,0], xrange = [750,800], \
                        IntLor = 800.0, mU = 777.0, gammaL = 2.5, num=250):
        
        global counter
        global colorlist
        global plot1
        
        labelUnit = 'cm\u207b\u00B9'
        
        # Setting the range based on raman wavelength
        limMin = ((self.arrayX <= xrange[0])*1).sum()
        limMax = ((self.arrayX <= xrange[1])*1).sum()
        
        xData = self.arrayX[limMin:limMax]
        yData = self.arrayYP[spot[0],spot[1],limMin:limMax]
        
        # Creating starting values for background slope and constant (based on linear equation)
        mVal = (yData[-1]-yData[0])/(xData[-1]-xData[0])
        cVal = yData[0] - mVal*xData[0]
        
        dataX = np.linspace(xrange[0], xrange[1], num) # Creating array for plotting fitting
        
        # It returns best fit values for parameters
        # Pass the function name, x data, y data, constant of parameters (in self.limMin:self.limMax)
        popt, pcov = curve_fit(fitLorentz, self.arrayX[limMin:limMax],\
                               self.arrayYP[spot[0],spot[1],limMin:limMax],\
                               [IntLor, gammaL, mU, mVal,cVal],\
                               sigma = self.arrayYP[spot[0],spot[1],limMin:limMax])
        
        # Standard deviation for all parameters
        stddvn = np.sqrt(np.diag(pcov))
        
        # Array values based on fitting parameters
        dataY = popt[0]*(popt[1]**2/((dataX-popt[2])**2+popt[1]**2)) +\
        (popt[3]*dataX + popt[4])
        
        # Error for IntLor and mU
        EIntLor = popt[0]+popt[3]*popt[2]+popt[4]
        
        EIntLor1 = EIntLor - 3*stddvn[0]
        EIntLor2 = EIntLor + 3*stddvn[0]
        EmU1 = popt[2] - 3*stddvn[2]
        EmU2 = popt[2] + 3*stddvn[2]
        
        # -------------------------Ploting------------------------------------
        plot1.circle(x = xData, y = yData, \
                     fill_color = colorlist[counter], fill_alpha=0.5, alpha=0.5, size = 8,
                    legend = self.spectrumlabel+str(spot))
        plot1.patch(x = dataX, y = dataY, \
                   color = colorlist[counter], line_width=2, alpha=0.25, \
                    legend = 'Lorentz '+self.spectrumlabel+str(spot))
        
        plot1.line(x = [popt[2],popt[2]], y = [EIntLor1,EIntLor2], \
                   color = colorlist[counter], line_width=3, alpha=0.5)
        
        plot1.line(x = [EmU1,EmU2], y = [EIntLor,EIntLor], \
                   color = colorlist[counter], line_width=3, alpha=0.5)
        plot1.line(x = [popt[2],popt[2]], y = [0,EIntLor+10], \
                   color = colorlist[counter], line_width=1, line_dash="4 4", alpha=0.5)
        
        if counter == 10:
            counter = 0
        else:
            counter = counter + 1
            
    def fitallspot(self, xrange = [750,800], IntLor = 800.0, mU = 777.0, gammaL = 2.5):
        
        # Setting the range based on raman wavelength
        limMin = ((self.arrayX <= xrange[0])*1).sum()
        limMax = ((self.arrayX <= xrange[1])*1).sum()
        
        xData = self.arrayX[limMin:limMax]
        yData = self.arrayYP[:,:,limMin:limMax]
        
        newyData = yData.reshape(-1, yData.shape[-1]) # flatten 3D into 2D
        
        # Creating starting values for background slope and constant (based on linear equation)
        mVal = (newyData[:,-1]-newyData[:,0])/(xData[-1]-xData[0])
        cVal = newyData[:,0] - mVal*xData[0]
        
        # Arrays for fitted parameter values and errors
        LDataAll = np.array([])
        LErrorAll = np.array([])
        
        # Fitting process
        for i,j in enumerate(newyData):
            popt, pcov = curve_fit(fitLorentz, xData,j,\
                                   [IntLor, gammaL, mU, mVal[i],cVal[i]], sigma = j)
    
            LDataAll = np.append(LDataAll,popt)
            LErrorAll = np.append(LErrorAll,np.sqrt(np.diag(pcov)))

        # Reshaping the data into (400,5) - pixel vs parameters
        LData = LDataAll.reshape(self.totalNumSpectrum,5)
        LError = LErrorAll.reshape(self.totalNumSpectrum,5)
        
        surfInt = LData[:,0].reshape(self.size[0],self.size[1]) # Intensity
        surfStress = LData[:,2].reshape(self.size[0],self.size[1]) # Stress
        
        # -------------------------Ploting------------------------------------
        plt.figure(figsize=(12,4))

        plt.subplot(1,2,1)
        plt.imshow(surfInt,cmap=plt.cm.gnuplot2, interpolation='spline36')
        plt.colorbar()
        title1 = ('Surface intensity at '+str(mU)+' cm$^{-1}$')
        plt.title(title1)
        
        plt.subplot(1,2,2)
        plt.imshow(surfStress,cmap=plt.cm.gnuplot2, interpolation='spline36')
        plt.colorbar()
        title2 = ('Surface stress at '+str(mU)+' cm$^{-1}$')
        plt.title(title2)
        
        plt.show()
        
        
                   