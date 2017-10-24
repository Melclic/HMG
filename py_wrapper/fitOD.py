##
# @file melPywafo.py 
#
# Contains all the functions required for the analysis of OD growth curves. Includes visual inspection of the gowth curve with first and second derivative to identify the exponential section of the growth curves. This contains user defined parameters for the smoothing parameter and the minimal length required for the growth curve. Using smoothing spline function, the volumetric changes of the growth curve, its instantaneous growth rate post exponential and exponential growth rate are returned. Volumetric changes are calculated based on the Volkmer et al. mean cell volume based on growth rate function. The smoothing spline (also called piecewise polynolmial) has been extracted from the pywafo packges
#
# @version 1.0
# @author Melchior du Lac
#

#TODO: use the expoentnial section of the growth curve detection as an indicator of section of the growth curve that is exponential, and have user defined input of start and end of the growht curve.

from __future__ import division
import numpy as np
import scipy.signal
import scipy.sparse as sparse
from numpy import ones, zeros, prod, sin, diff, pi, inf, vstack, linspace
from scipy.interpolate import interp1d
import math
import csv

import polynomial as pl

import matplotlib.pyplot as plt

#########################################################################################
################################### pywafo objects ######################################
#########################################################################################

class PPform(object):

    """The ppform of the piecewise polynomials
                    is given in terms of coefficients and breaks.
    The polynomial in the ith interval is
        x_{i} <= x < x_{i+1}

    S_i = sum(coefs[m,i]*(x-breaks[i])^(k-m), m=0..k)
    where k is the degree of the polynomial.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> coef = np.array([[1,1]]) # unit step function
    >>> coef = np.array([[1,1],[0,1]]) # linear from 0 to 2
    >>> coef = np.array([[1,1],[1,1],[0,2]]) # linear from 0 to 2
    >>> breaks = [0,1,2]
    >>> self = PPform(coef, breaks)
    >>> x = linspace(-1,3)
    >>> h=plt.plot(x,self(x))
    """

    def __init__(self, coeffs, breaks, fill=0.0, sort=False, a=None, b=None):
        if sort:
            self.breaks = np.sort(breaks)
        else:
            self.breaks = np.asarray(breaks)
        if a is None:
            a = self.breaks[0]
        if b is None:
            b = self.breaks[-1]
        self.coeffs = np.asarray(coeffs)
        self.order = self.coeffs.shape[0]
        self.fill = fill
        self.a = a
        self.b = b

    def __call__(self, xnew):
        saveshape = np.shape(xnew)
        xnew = np.ravel(xnew)
        res = np.empty_like(xnew)
        mask = (self.a <= xnew) & (xnew <= self.b)
        res[~mask] = self.fill
        xx = xnew.compress(mask)
        indxs = np.searchsorted(self.breaks[:-1], xx) - 1
        indxs = indxs.clip(0, len(self.breaks))
        pp = self.coeffs
        dx = xx - self.breaks.take(indxs)

        v = pp[0, indxs]
        for i in range(1, self.order):
            v = dx * v + pp[i, indxs]
        values = v

        res[mask] = values
        res.shape = saveshape
        return res

    def linear_extrapolate(self, output=True):
        '''
        Return 1D PPform which extrapolate linearly outside its basic interval
        '''

        max_order = 2

        if self.order <= max_order:
            if output:
                return self
            else:
                return
        breaks = self.breaks.copy()
        coefs = self.coeffs.copy()
        # pieces = len(breaks) - 1

        # Add new breaks beyond each end
        breaks2add = breaks[[0, -1]] + np.array([-1, 1])
        newbreaks = np.hstack([breaks2add[0], breaks, breaks2add[1]])

        dx = newbreaks[[0, -2]] - breaks[[0, -2]]

        dx = dx.ravel()

        # Get coefficients for the new last polynomial piece (a_n)
        # by just relocate the previous last polynomial and
        # then set all terms of order > maxOrder to zero

        a_nn = coefs[:, -1]
        dxN = dx[-1]

        a_n = pl.polyreloc(a_nn, -dxN)  # Relocate last polynomial
        # set to zero all terms of order > maxOrder
        a_n[0:self.order - max_order] = 0

        # Get the coefficients for the new first piece (a_1)
        # by first setting all terms of order > maxOrder to zero and then
        # relocate the polynomial.

        # Set to zero all terms of order > maxOrder, i.e., not using them
        a_11 = coefs[self.order - max_order::, 0]
        dx1 = dx[0]

        a_1 = pl.polyreloc(a_11, -dx1)  # Relocate first polynomial
        a_1 = np.hstack([zeros(self.order - max_order), a_1])

        newcoefs = np.hstack([a_1.reshape(-1, 1), coefs, a_n.reshape(-1, 1)])
        if output:
            return PPform(newcoefs, newbreaks, a=-inf, b=inf)
        else:
            self.coeffs = newcoefs
            self.breaks = newbreaks
            self.a = -inf
            self.b = inf

    def derivative(self):
        """
        Return first derivative of the piecewise polynomial
        """

        cof = pl.polyder(self.coeffs)
        brks = self.breaks.copy()
        return PPform(cof, brks, fill=self.fill)

    def integrate(self):
        """
        Return the indefinite integral of the piecewise polynomial
        """
        cof = pl.polyint(self.coeffs)

        pieces = len(self.breaks) - 1
        if 1 < pieces:
            # evaluate each integrated polynomial at the right endpoint of its
            # interval
            xs = diff(self.breaks[:-1, ...], axis=0)
            index = np.arange(pieces - 1)

            vv = xs * cof[0, index]
            k = self.order
            for i in range(1, k):
                vv = xs * (vv + cof[i, index])

            cof[-1] = np.hstack((0, vv)).cumsum()

        return PPform(cof, self.breaks, fill=self.fill)


class SmoothSpline(PPform):

    """
    Cubic Smoothing Spline.

    Parameters
    ----------
    x : array-like
        x-coordinates of data. (vector)
    y : array-like
        y-coordinates of data. (vector or matrix)
    p : real scalar
        smoothing parameter between 0 and 1:
        0 -> LS-straight line
        1 -> cubic spline interpolant
    lin_extrap : bool
        if False regular smoothing spline
        if True a smoothing spline with a constraint on the ends to
        ensure linear extrapolation outside the range of the data (default)
    var : array-like
        variance of each y(i) (default  1)

    Returns
    -------
    pp : ppform
        If xx is not given, return self-form of the spline.

    Given the approximate values

        y(i) = g(x(i))+e(i)

    of some smooth function, g, where e(i) is the error. SMOOTH tries to
    recover g from y by constructing a function, f, which  minimizes

      p * sum (Y(i) - f(X(i)))^2/d2(i)  +  (1-p) * int (f'')^2


    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(0,1)
    >>> y = np.exp(x)+1e-1*np.random.randn(x.size)
    >>> pp9 = SmoothSpline(x, y, p=.9)
    >>> pp99 = SmoothSpline(x, y, p=.99, var=0.01)
    >>> h=plt.plot(x,y, x,pp99(x),'g', x,pp9(x),'k', x,np.exp(x),'r')

    See also
    --------
    lc2tr, dat2tr


    References
    ----------
    Carl de Boor (1978)
    'Practical Guide to Splines'
    Springer Verlag
    Uses EqXIV.6--9, self 239
    """

    def __init__(self, xx, yy, p=None, lin_extrap=True, var=1):
        coefs, brks = self._compute_coefs(xx, yy, p, var)
        super(SmoothSpline, self).__init__(coefs, brks)
        if lin_extrap:
            self.linear_extrapolate(output=False)

    def _compute_coefs(self, xx, yy, p=None, var=1):
        x, y = np.atleast_1d(xx, yy)
        x = x.ravel()
        dx = np.diff(x)
        must_sort = (dx < 0).any()
        if must_sort:
            ind = x.argsort()
            x = x[ind]
            y = y[..., ind]
            dx = np.diff(x)

        n = len(x)

        # ndy = y.ndim
        szy = y.shape

        nd = prod(szy[:-1])
        ny = szy[-1]

        if n < 2:
            raise ValueError('There must be >=2 data points.')
        elif (dx <= 0).any():
            raise ValueError('Two consecutive values in x can not be equal.')
        elif n != ny:
            raise ValueError('x and y must have the same length.')

        dydx = np.diff(y) / dx

        if (n == 2):  # % straight line
            coefs = np.vstack([dydx.ravel(), y[0, :]])
        else:

            dx1 = 1. / dx
            D = sparse.spdiags(var * ones(n), 0, n, n)  # The variance

            u, p = self._compute_u(p, D, dydx, dx, dx1, n)
            dx1.shape = (n - 1, -1)
            dx.shape = (n - 1, -1)
            zrs = zeros(nd)
            if p < 1:
                # faster than yi-6*(1-p)*Q*u
                Qu = D * diff(vstack([zrs, diff(vstack([zrs, u, zrs]),
                                                axis=0) * dx1, zrs]), axis=0)
                ai = (y - (6 * (1 - p) * Qu).T).T
            else:
                ai = y.reshape(n, -1)

            # The piecewise polynominals are written as
            # fi=ai+bi*(x-xi)+ci*(x-xi)^2+di*(x-xi)^3
            # where the derivatives in the knots according to Carl de Boor are:
            #    ddfi  = 6*p*[0;u] = 2*ci;
            #    dddfi = 2*diff([ci;0])./dx = 6*di;
            #    dfi   = diff(ai)./dx-(ci+di.*dx).*dx = bi;

            ci = np.vstack([zrs, 3 * p * u])
            di = (diff(vstack([ci, zrs]), axis=0) * dx1 / 3)
            bi = (diff(ai, axis=0) * dx1 - (ci + di * dx) * dx)
            ai = ai[:n - 1, ...]
            if nd > 1:
                di = di.T
                ci = ci.T
                ai = ai.T
            if not any(di):
                if not any(ci):
                    coefs = vstack([bi.ravel(), ai.ravel()])
                else:
                    coefs = vstack([ci.ravel(), bi.ravel(), ai.ravel()])
            else:
                coefs = vstack(
                    [di.ravel(), ci.ravel(), bi.ravel(), ai.ravel()])

        return coefs, x

    def _compute_u(self, p, D, dydx, dx, dx1, n):
        if p is None or p != 0:
            data = [dx[1:n - 1], 2 * (dx[:n - 2] + dx[1:n - 1]), dx[:n - 2]]
            R = sparse.spdiags(data, [-1, 0, 1], n - 2, n - 2)

        if p is None or p < 1:
            Q = sparse.spdiags(
                [dx1[:n - 2], -(dx1[:n - 2] + dx1[1:n - 1]), dx1[1:n - 1]],
                [0, -1, -2], n, n - 2)
            QDQ = (Q.T * D * Q)
            if p is None or p < 0:
                # Estimate p
                p = 1. / \
                    (1. + QDQ.diagonal().sum() /
                     (100. * R.diagonal().sum() ** 2))

            if p == 0:
                QQ = 6 * QDQ
            else:
                QQ = (6 * (1 - p)) * (QDQ) + p * R
        else:
            QQ = R

        # Make sure it uses symmetric matrix solver
        ddydx = diff(dydx, axis=0)
        # sp.linalg.use_solver(useUmfpack=True)
        u = 2 * sparse.linalg.spsolve((QQ + QQ.T), ddydx)  # @UndefinedVariable
        return u.reshape(n - 2, -1), p

#########################################################################################
#########################################################################################

##
# @brief Fill the gaps in the OD input data through linear interpolation
#
# Because the smoothing spline algorithm does not deal well with gaps in data, if that is the case we fill the gaps assuming that the difference between rwo points is lienar. Two points is determined two have a gap if distance between the two is larger than the mean of the sampling intervals. The gap is then filled assuming linear interpolation with equal sampling frequency.
#
# @param time 1D array of x-axis of the OD growth curve
# @param od 1D array of the y-axis of the OD growth curve
#
# @return time Filled and parsed (minutes) x-axis of the growth curve
# @return od Filled and parsed y-axis of the growth curve
# @return parcentMissing Calculated percentage of the data that hase been filled
#

#TODO: fix WARNING: Buggy, if there is nothing to fill, then this method attempts to divide by 0
def gapFill(time, od):
 oldODLen = len(od)
 oldTimeLen = len(time)
 #calculate the intervals between each time points
 gapArray = []
 for i in range(len(time)-1):
  gapArray.append(abs(time[i]-time[i+1]))
 meanGap = sum(gapArray)/float(len(gapArray))
 #gap values for time and OD that are bigger than twice the mean where we will fill them
 gapsTime = []
 gapsOD = []
 #gap values that value is less than the mean and we will not fill
 nonGaps = []
 count = 0
 for i in gapArray:
  if i>(2.0*meanGap):
   gapsTime.append([time[count],time[count+1]])
   gapsOD.append([od[count],od[count+1]])
  elif i<meanGap:
   nonGaps.append(i)
  count += 1
 nonGapMean = 0.0
 if len(nonGaps)==0:
  return time, od, 0.0
 nonGapMean = sum(nonGaps)/float(len(nonGaps))
 #Fill the gaps using linear interpolation
 for i in range(len(gapsTime)):
  linearFit = interp1d(gapsTime[i], gapsOD[i])
  filledTime = np.arange(gapsTime[i][0]+nonGapMean, gapsTime[i][1], nonGapMean) #+1 to include the last value
  filledOD = linearFit(filledTime)
  for i in range(len(time)-1):
   if time[i]<filledTime[0] and time[i+1]>filledTime[-1]:
    posAdd = i
    for y in range(len(filledTime)):
     posAdd += 1
     time.insert(posAdd,filledTime[y])
     od.insert(posAdd,filledOD[y])
    break
 percentMissing = 0.0
 if(oldTimeLen<len(time)):
  percentMissing = math.ceil(100.0*(len(od)-oldODLen)/len(od))
 return time, od, percentMissing

##
# @brief Fill the end of the growth curve that is descending with flat OD changes
#
# Because our modelling strategy cannot handle diminishing growth curves, we replace the death phase with linear unchanging growth. Future plans involve the implementation of cell death so as to be able to model decreasing OD values.
#
# @param od 1D array with the OD values
#
# @return od 1D array with parsed OD values
#
def linearEnds(od):
 startPin = np.median(od[0:5])
 endPin = np.median(od[len(od)-5:-1])
 if od[0]!=min(od):
  for i in range(5):
   od.insert(i, startPin)
 if od[-1]!=max(od):
  for i in range(len(od), len(od)-5, -1):
   od.insert(i, endPin)
 return od

#################################################################
############################# Plot ##############################
#################################################################

##
# @brief Visual inspection of the growth curve with currently calculated exponential phase, with first and second order polynomials
#
# To check that the detected exponential section of the growth curve indeed is valid, we plot the OD growth curve (log scale), its first and second orde differential. Theoretically, the section of the growth curve that is linear should be the section of the growth curve where the second order differential is linear. 
#
# @param fitTime Spline interpolation fit time output with dt as time step
# @param fitOD Spline interpolation fit OD output with dt as time step
# @param time Original 1D array of the time of the growth curve
# @param od Original 1D array of the OD of the growth curve
# @param diffOD First order differential of fitOD (not used)
# @param isNegOD Boolean determining if there are any parts of the fit growth curve that is negative and thus invalid
# @param isNegDeriv Boolean determining if there are any parts of the first order differential of the growth curve that is negative and thus inavalid
# @param expoTau Calculated doubling rate of the exponential section of the growth curve 
# @param startLinear Time in minutes of the start of the linear phase
# @param stopLinear Time in minutes of the stop of the linear phase
#
def plotIt(fitTime, fitOD, time, od, diffOD, isNegOD, isNegDeriv, expoTau, startLinear, stopLinear):
 plt.subplot(3,1,1)
 #plt.semilogy(time, od, 'o', label='Original Data')
 #plt.semilogy(fitTime, fitOD, label='Spline Fit')
 plt.plot(time, np.log(od), 'o', label='Original Data')
 plt.plot(fitTime, np.log(fitOD), label='Spline Fit')
 #plt.axvline(maxOD, color='red', label='Max OD ('+str(round(expoTau,2))+'min-1)')
 if not startLinear==-1.0 and not stopLinear==-1.0:
  plt.axvline(startLinear, color='red', label='Start Linear (tau: '+str(round(expoTau,2))+')')
  plt.axvline(stopLinear, color='yellow', label='Stop Linear')
 plt.xlabel('Time (min)')
 plt.ylabel('ln(OD)')
 plt.title('Valid OD: '+str(not isNegOD)+' | Valid diff(OD): '+str(not isNegDeriv))
 plt.legend(loc=4)
 plt.subplot(3,1,2)
 #plt.plot(fitTime, diffOD)
 #startLoc = -1.0
 #for i in range(len(fitTime)):
 # if fitTime[i]>=startLinear:
 #  startLoc = i
 #  break
 #ft = fitTime[:]
 #ft = ft[0:-1]
 plt.plot(fitTime[0:-1], np.diff(np.log(fitOD)))
 #plt.plot(ft[startLoc:len(ft)], np.gradient(np.log(fitOD)[startLoc:-1]))
 plt.xlabel('Time (min)')
 plt.ylabel('F1(ln(OD))')
 plt.subplot(3,1,3)
 plt.plot(fitTime[0:-2], np.diff(np.log(fitOD), 2))
 plt.xlabel('Time (min)')
 plt.ylabel('F2(ln(OD))')
 plt.show()

 """
 plt.semilogy(time, od, 'o', label='Original Data')
 plt.semilogy(fitTime, fitOD, label='Spline Fit')
 #plt.plot(time, od, 'o', label='Original Data')
 #plt.plot(fitTime, fitOD, label='Spline Fit')
 plt.axvline(maxOD, color='red', label='Max OD')
 plt.xlabel('Time (min)')
 plt.ylabel('OD')
 plt.title('Valid OD: '+str(not isNegOD)+' | Valid diff(OD): '+str(not isNegDeriv))
 plt.legend(loc=4)
 plt.show()
 """

################################################################
################################################################
################################################################

##
# @brief Try a new smoothing parameter 
#
# Because this is user defined smoothing and exponential window determined function, we must recalculate the fitOD
#
# @param fitTime 1D array of the fitted time using dt as tim step
# @param time Original time array
# @param od Original od array
# @param p Smoothing parameter
# @param dt Time step
# @param confidenceBound confidence of the linear fit to the OD growth curve
# @param minFlatSize Minimal size of the array that is the exponential section of the growth curve 
#
# @return fitOD Spline fit OD values
# @return fitDerivOD Spline fit first order derivative OD values
# @return maxOD Maximal OD values of the spline fit
# @return isNegOD Boolean determining if the spline fit contains negative OD values
# @return isNegDeriv Boolean determining the the first order derivative of the spline fit contains decreasing values
# @return expoTau Growth rate of the exponential section of the growth curve
# @return ss SmoothingSpline object 
# @return startLinear Start of the linear section of the growth curve
# @return stopLinear End of the linear section of the growth curve
#
def tryNewP(fitTime, time, od, p, dt, confidenceBound, minFlatSize):
 ss = SmoothSpline(time, od, p=p)
 fitOD = ss(fitTime)
 ssDeriv = ss.derivative()
 #fitDerivOD = ssDeriv(fitTime).tolist()
 fitDerivOD = np.gradient(np.log(fitOD))
 #lnDeriv = np.diff(np.log(fitOD)).tolist()
 lnDeriv = np.gradient(np.log(fitOD)).tolist()
 #cannot use max() function because it may be at time 0 --> data we cannot trust
 #need to identify and ignore the lag phase
 #remove anythin with an OD of X
 maxOD = fitTime[lnDeriv.index(max(lnDeriv))]
 #the linear growth rate
 numSide = int(10.0/dt)
 fitTime = fitTime.tolist()
 expoTime = fitTime[fitTime.index(maxOD)-numSide:fitTime.index(maxOD)+numSide]
 expoOD = fitOD[fitTime.index(maxOD)-numSide:fitTime.index(maxOD)+numSide]
 #test: identify the linear part of growth: second derivative of the natural lof of the od curve. If it is 0 then it is linear. However because of the data must have a confidence interval
 secondDeriv = np.diff(np.log(fitOD), 2)
 startLinear = -1.0
 stopLinear = -1.0
 linearTime = []
 linearOD = []
 count = 0
 for i in range(len(secondDeriv)):
  if confidenceBound>=secondDeriv[i]>=-confidenceBound:
   if startLinear==-1.0:
    startLinear = fitTime[i]
   else:
    if count>=minFlatSize:
     stopLinear = fitTime[i]
    count += 1
   linearTime.append(fitTime[i])
   linearOD.append(fitOD[i])
  else:
   if not startLinear==-1.0 and not stopLinear==-1.0:
    count = 0
    break
   if not startLinear==-1.0 and stopLinear==-1.0:
    linearTime = []
    linearOD = []
    startLinear = -1.0
    count = 0    
   if startLinear==-1.0 and stopLinear==-1.0:
    linearTime = []
    linearOD = [] 

 startLoc = -1.0
 for i in range(len(fitTime)):
  if fitTime[i]>=stopLinear:
   startLoc = i
   break
 #[ unicode(x.strip()) if x is not None else '' for x in row ]
 fitDerivOD = [i if not 0.0>i>-0.0000005 else 0.0 for i in fitDerivOD]
 fitDerivOD = [i*100.0 for i in fitDerivOD]
   
 isNegOD = any(i<0.0 for i in fitOD[startLoc:len(fitOD)])
 isNegDeriv = any(y<0.0 for y in fitDerivOD[startLoc:len(fitDerivOD)])
 #expoTau = math.log(2)/np.polyfit(linearTime, np.log(linearOD), 1)[0]
 expoTau = math.log(2)/fitDerivOD[startLoc]
 return fitOD, fitDerivOD, maxOD, isNegOD, isNegDeriv, expoTau, ss, startLinear, stopLinear

#############################################################
###################### USER SMOOTHIG ########################
#############################################################

##
# @brief Matplotlib plotting of the OD spline fit with user input of smoothing 
#
# Given that the smoothing parameter, minimal length of the linear section of growth, and the confidence bounds, this functions plots the results and requires the user inout of all of these three parameters. This also checks that there are no decreasing OD values and that the fit is valid overall 
#
# @param time Orginial measured time values
# @param od Original measured OD values
# @param dt Time step
#
# @return fitTime Spline fit time
# @return fitOD Spline fit OD values
# @return fitDerivOD Spline fit first order derivative OD values
# @return p Smoothing parameter of the 
# @return maxOD Maximal OD values of the spline fit
# @return expoTau Growth rate of the exponential section of the growth curve
# @return fitObj SmoothingSpline object 
# @return startLinear Start of the linear section of the growth curve
# @return stopLinear End of the linear section of the growth curve
# @return minFlatSize Minimal size of the array that is the exponential section of the growth curve 
# @return confidenceBound confidence of the linear fit to the OD growth curve
#
#TODO: Make this in something else (javascript), not matplotlib
def userSmoothing(time, od, dt):
 spacing = np.mean([math.fabs(time[i+1]-time[i]) for i in range(len(time)-2)])
 p = 1.0/(1.0+math.pow(np.mean(spacing),3.0)/6.0)
 fitTime = np.arange(time[0], time[-1], dt)
 confidenceBound = 0.000000001 
 minFlatSize = 3000 

 #identify maxOD, and if not OD[-1] fill it with maxOD until next Flow Measurement
 if max(od)!=od[-1]:
  maxODIndex = od.index(max(od))
  od = od[0:maxODIndex]
  time = time[0:maxODIndex]
 
 time, od, percentMissing = gapFill(time[:], od[:])
 fitOD, fitDerivOD, maxOD, isNegOD, isNegDeriv, expoTau, fitObj, startLinear, stopLinear = tryNewP(fitTime, time, od, p/(1+1), dt, confidenceBound, minFlatSize)

 #isValid = False
 #for i in range(5):
 # fitOD, fitDerivOD, maxOD, isNegOD, isNegDeriv, expoTau, fitObj, startLinear, stopLinear = tryNewP(fitTime, time, od, p/(i+1), dt, confidenceBound)
 # if not isNegOD and not isNegDeriv:
 #  isValid = True
 #  break
 #if not isValid:
 # time, od, percentMissing = gapFill(time[:], od[:])

 fitTime = np.arange(time[0], time[-1], dt)
 
 plt.ion()
 fitOD, fitDerivOD, maxOD, isNegOD, isNegDeriv, expoTau, fitObj, startLinear, stopLinear = tryNewP(fitTime, time, od, p, dt, confidenceBound, minFlatSize)
 print(startLinear)
 print(stopLinear)
 plotIt(fitTime, fitOD, time, od, fitDerivOD, isNegOD, isNegDeriv, expoTau, startLinear, stopLinear)   
 
 while True:
  print('p: '+str(p))
  inP = input('Input different p value?: ') 
  print('confidenceBound: '+str(confidenceBound))
  inConfB = input('Input different confidenceBound value?: ') 
  print('minFlatSize: '+str(minFlatSize))
  mfs = input('Input different minFlatSize?: ')
  try:
   p = float(inP)
   confidenceBound = float(inConfB)
   minFlatSize = float(mfs)
   plt.clf()
   fitOD, fitDerivOD, maxOD, isNegOD, isNegDeriv, expoTau, fitObj, startLinear, stopLinear = tryNewP(fitTime, time, od, p, dt, confidenceBound, minFlatSize)
   plotIt(fitTime, fitOD, time, od, fitDerivOD, isNegOD, isNegDeriv, expoTau, startLinear, stopLinear)
   print(startLinear)
   print(stopLinear)
  except ValueError:
   if not isNegOD and not isNegDeriv:
    plt.close()
    break
   else:
    print('The smoothing parameters are invalid... Try again')
    print('isNegOD: '+str(isNegOD)) 
    print('isNegDeriv: '+str(isNegDeriv)) 
 return fitTime, fitOD, fitDerivOD, p, maxOD, expoTau, fitObj, startLinear, stopLinear, minFlatSize, confidenceBound

##############################################################
######################## MAIN   ##############################
##############################################################

##
# @brief Main function that calculates the different 
#
# @param inputFile Path to the input file with the histograms
# @param inputTime Time points from the OD in minutes
# @param inputOD Growth rate input OD
# @param polyChannel parameters of the first order polynomial that converts the channel DNA content
# @param dt Time step
#
# @return flowTimes Ouput flow cytometry times in minutes
# @return flowData Output flow cytometry histograms
# @return flowScale Output flow cytometry DNA x-axis scale
# @return injectionOD Spline fit OD data based on dt time step
# @return injectionGR Spline fit instateneous growth rate based on dt
# @return injectionGrownMass Calculated spline fit total volumetric changes of the growth curve based on dt
# @return injectionTime Spline fit time based on dt
# @return expoTau Linear section of the growth curve growth rate
#
def genData(inputTime, inputOD, dt):
 # p Smoothing parameter
 # confidenceBound Confidence bounds of the linear fit of the growth curve
 # minFlatSize Minimal size of the linear section of the growth curve
 fitTime, fitOD, fitDerivOD, p, maxOD, expoTau, fitObj, startLinear, endLinear, minFlatSize, confidenceBound = userSmoothing(inputTime, inputOD, dt)
 
 injectionTime = [np.arange(flowTimes[i], flowTimes[i+1]+dt, dt) for i in range(len(flowTimes)-1)]
 injectionOD = [list(fitObj(i)) for i in injectionTime]
 injectionGR = [list(np.gradient(np.log(i))*100.0) for i in injectionOD]
 injectionGrownMass = [[(3.6*y)*math.pow(10.0, 9.0) for y in i] for i in injectionOD]

 return flowTimes, flowData, flowScale, injectionOD, injectionGR, injectionGrownMass, injectionTime, expoTau
