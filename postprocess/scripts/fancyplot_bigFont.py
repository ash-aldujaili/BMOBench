#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import palettable as brewer2mpl  # for colors, see: Documentation @ https://github.com/jiffyclub/brewer2mpl/wiki
import numpy as np


class PlotWithFancyLegend:
    """ This class implements the plotting of some statistic of the data with confidence interval region.
    """
    def __init__(self, figsize=(7, 5), position=[0.1, 0.11, 0.87, 0.81]):
        #======================================#
        # Setting up sizes, positioning params #
        #======================================#
        self.figsize = figsize
        self.position = position
        
        self.title = "Fancy plot with legend on the right"
        self.xlabel = "xlabel"
        self.ylabel = "ylabel"
        self.legends = []
        #self.xlimPerc = [0/100., 100/100.]
        self.lineWidth = 1.5
        self.markerSize = 7.5
        self.nMarkers = 5
        self.alpha = 0.9
        self.xExtRatio    = 23/100.
        self.xSegLenRatio = 17/100.
        self.yShrinkRatio = 30/100.
        self.labelBottomTop = [0.5, 1.5] # [1, 1] # # a factor determines the the spread of the right legends
        self.rightLegend = True
        self.keepBox = True
        self.axesLineWidth = 1
        # Newly added properties
        self.topFontSize        = 16 # fontsize for plot title
        self.bottomFontSize     = 15 # fontsize for x-label
        self.leftFontSize       = 15 # fontsize for y-label
        self.rightFontSize      = 12 # fontsize for legends on the right
        self.annotationFontSize = 15 # fontsize for additional texts inside the plot
        self.xTickLabelSize  = 12
        self.yTickLabelSize  = 12
        # For grid lines
        self.gridLineStyle  = ':'
        self.gridLineWeight = 0.35
    
    
    def plotting(self, filename, data, xstart, logx=False, yfull=True, rainbow=None, markers=None, pngdpi=300):
        """ data is expected to be a dictionary of n profiles,
        each profile is again a dictionary {'x': [], 'y': []}.
        Each profile will be plotted in a curve.
        """
        if not self.legends:
            for i in xrange(len(data)):
                self.legends.append("legends[%d]" % (i))
        if rainbow is None:
            rainbow = ['b','g','r','c','m','orange','y','k','silver','coral','lime','brown','violet','navy','greenyellow']
            # rainbow = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors  # see: Documentation @ https://github.com/jiffyclub/brewer2mpl/wiki
        if markers is None:
            markers = ['o','s','*','v','p','D','^','+','<','>','d','H','x']
        zorder0 = 5 # the lowest zorder for major curves
        
        # Sort the keys list as dictionary is non-ordered
        # and taking care of the ordering in legends as well
        data_keys = data.keys()
        data_keys.sort()
        order_legends = []
        for key in data_keys:
            idx = data.keys().index(key)
            order_legends.append(self.legends[idx])
        
        #====================================#
        #    Plotting commands start here    #
        #====================================#
        #----------------------------------------------
        import matplotlib
        # Use a non-interactive backend such as Agg (for PNGs), PDF, SVG, or PS.
        matplotlib.use('Agg') # make sure to call this before pyplot
        #----------------------------------------------
        import matplotlib.pyplot as plt
        # plt.rc('text', usetex = True)
        fig = plt.figure(1, figsize=self.figsize)
        
        # ax = fig.add_subplot(111)
        ax = fig.add_axes(self.position) # self.position = [left, bottom, width, height]
        
        # Plot the main data:
        #=====================
        for ialg, alg in enumerate(data_keys):
            if not logx:
                ax.plot(data[alg]['x'], data[alg]['y'],
                        drawstyle='steps-post', clip_on=False,
                        color=rainbow[ialg % len(rainbow)],
                        lw=self.lineWidth, alpha=self.alpha, zorder=zorder0+ialg)
            else:
                ax.semilogx(data[alg]['x'], data[alg]['y'],
                            drawstyle='steps-post', clip_on=False,
                            color=rainbow[ialg % len(rainbow)],
                            lw=self.lineWidth, alpha=self.alpha, zorder=zorder0+ialg)
        
        xLim = [xstart, max([max(data[alg]['x']) for alg in data_keys])]
        if logx and xLim[0]==0: xLim[0] = 1
        
        # Plot the horizontal 'extended' line from the end of each curve to xLim[1]
        # and also put a cross to mark the end of the curve
        #==========================================================================
        for ialg, alg in enumerate(data_keys):
            if data[alg]['x'][-1] < xLim[1]:
                if not logx:
                    # Extended horizontal line
                    ax.plot([data[alg]['x'][-1], xLim[1]], [data[alg]['y'][-2], data[alg]['y'][-2]],
                            clip_on=False,
                            color=rainbow[ialg % len(rainbow)],
                            lw=0.9*self.lineWidth, alpha=self.alpha, zorder=zorder0+ialg)
                    # Add a cross to mark the end of the curve
                    ax.plot([data[alg]['x'][-1]], [data[alg]['y'][-2]],
                            marker='x', markeredgecolor=rainbow[ialg % len(rainbow)],
                            markersize=2.7*self.markerSize,
                            markeredgewidth=0.8*self.lineWidth,
                            clip_on=False,
                            color=rainbow[ialg % len(rainbow)],
                            lw=self.lineWidth, alpha=1, zorder=zorder0+ialg)
                else:
                    # Extended horizontal line
                    ax.semilogx([data[alg]['x'][-1], xLim[1]], [data[alg]['y'][-2], data[alg]['y'][-2]],
                                clip_on=False,
                                color=rainbow[ialg % len(rainbow)],
                                lw=0.9*self.lineWidth, alpha=self.alpha, zorder=zorder0+ialg)
                    # Add a cross to mark the end of the curve
                    ax.semilogx([data[alg]['x'][-1]], [data[alg]['y'][-2]],
                                marker='x', markeredgecolor=rainbow[ialg % len(rainbow)],
                                markersize=2.7*self.markerSize,
                                markeredgewidth=0.8*self.lineWidth,
                                clip_on=False,
                                color=rainbow[ialg % len(rainbow)],
                                lw=self.lineWidth, alpha=1, zorder=zorder0+ialg)
        
        # Plot markers for the curves:
        #==============================
        #xarr = np.arange(xstart, data[data_keys[0]]['x'][-1] + 1) # generate data for the x-axis
        #print xarr[-1]
        #print data[data_keys[0]]['x'][-1] # TODO: debug ecdf.py: why 40000.2 ???
        nMarkers = self.nMarkers  # 5 markers on each line
        for ialg, alg in enumerate(data_keys):
            if not logx:
                lenx   = xLim[1] - xLim[0]
                alt    = int(lenx / (nMarkers*len(data))) # alternate among the first markers over lines (in idx unit)
                offset = int(lenx / nMarkers) # offset between 2 consecutive markers of a line (in idx unit)
                # Generate estimated x's for markers
                estxMarkers = int(alt/2) + np.arange(start=ialg*alt, stop=lenx-0.5*offset, step=offset, dtype=int)
            else:
                lenx   = np.log10(xLim[1]) - np.log10(xLim[0])
                alt    = lenx / (nMarkers*len(data)) # alternate among the first markers over lines (in idx unit)
                offset = lenx / nMarkers # offset between 2 consecutive markers of a line (in idx unit)
                # Generate estimated x's for markers
                estxMarkers = alt/2 + np.log10(xLim[0]) + np.arange(start=ialg*alt, stop=lenx-0.5*offset, step=offset)
                estxMarkers = 10 ** estxMarkers
                if ialg==0:
                    estxMarkers = np.delete(estxMarkers, 0) # skip the very first (not so visible) marker on the semilogx scale
            # Sample the real x and y of the markers from the curve
            markerCoord = {'x': [], 'y': []}
            icur = 0
            for estx in estxMarkers:
                for idx, x in enumerate(data[alg]['x']):
                    if idx < icur: continue
                    
                    # TODO: this is added to avoid a list index error with DIFF maxfevals
                    try: data[alg]['x'][idx+1]
                    except: continue
                    
                    if estx == x or (estx > x and estx < data[alg]['x'][idx+1]):
                        markerCoord['x'].append(x)
                        markerCoord['y'].append(data[alg]['y'][idx])
                        icur = idx + 1
                        break
            # Plot the sampled markers. TODO: alpha doesn't work for marker!
            for x, y in zip(markerCoord['x'], markerCoord['y']):
                if not logx:
                    ax.plot(x, y, ls='', clip_on=False, markerfacecolor='none',
                            marker=markers[ialg % len(markers)], markeredgecolor=rainbow[ialg % len(rainbow)],
                            markersize=self.markerSize*(1 if markers[ialg % len(markers)]!='*' else 1.45),
                            markeredgewidth=0.7*self.lineWidth, alpha=self.alpha, zorder=zorder0+ialg)
                else:
                    ax.semilogx(x, y, ls='', clip_on=False, markerfacecolor='none',
                                marker=markers[ialg % len(markers)], markeredgecolor=rainbow[ialg % len(rainbow)],
                                markersize=self.markerSize*(1 if markers[ialg % len(markers)]!='*' else 1.45),
                                markeredgewidth=0.7*self.lineWidth, alpha=self.alpha, zorder=zorder0+ialg)
        
        
        ax.set_xlim(xLim)
        if yfull:
            ax.set_ylim([0, 1])
        
        ax.grid(True, which='both', color="gray", alpha=0.6, ls=self.gridLineStyle, lw=self.gridLineWeight) # print visibility?
        #ax.grid(True, which='both', color="gray", alpha=0.1, ls='-', lw=0.3)
        #ax.grid(True, which='both', axis='x', color="gray", alpha=0.15, ls='-', lw=0.2)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(self.axesLineWidth)
        
        
        #=============================#
        #  Make legends on the right  #
        #=============================#
        if self.rightLegend:
            xExtRatio      = self.xExtRatio
            xSegLenRatio   = self.xSegLenRatio
            yShrinkRatio   = self.yShrinkRatio
            labelBottomTop = self.labelBottomTop
            
            yLim = ax.get_ylim()
            xTicks = ax.get_xticks()  # backup xticks for use after making legends on the right
            #if xTicks[0] < xLim[0]:   # the 1st element usually doesn't appear on the axis
            #    xTicks = np.delete(xTicks, 0)  # thus remove it of the backed up xTicks
            while xTicks[-1] > xLim[1]:
                xTicks = np.delete(xTicks, -1)
            if logx:
                while xTicks[0] < xLim[0]:
                    xTicks = np.delete(xTicks, 0)
            
            endData = [data[key]['y'][-1] for key in data_keys]
            for idx, key in enumerate(data_keys):
                j = -1
                while np.isnan(endData[idx]):
                    endData[idx] = data[key]['y'][j]
                    j = j - 1
            idx = np.argsort(endData)
            xLength = xLim[1] - xLim[0]
            yLength = yLim[1] - yLim[0]
            if not logx:
                normFactor = (xLim[1] - xLim[0]) / xLim[1]
                xExt = xLim[1] + xExtRatio*xSegLenRatio * normFactor*xLength
            else:
                normFactor = (np.log10(xLim[1]) - np.log10(xLim[0])) / np.log10(xLim[1]) # helps adjust the extension parts when xstart is large
                xExt = 10 ** (np.log10(xLim[1]) + xExtRatio*xSegLenRatio * normFactor*np.log10(xLength))
            yExt = np.linspace(yLim[0] + labelBottomTop[0]*yShrinkRatio*yLength,
                               yLim[1] - labelBottomTop[1]*yShrinkRatio*yLength,
                               num=len(data), endpoint=True)
            
            # Plot all extension segments:
            #==============================
            for k, alg in enumerate(data_keys):
                if not logx:
                    ax.plot(np.array([xLim[1], xExt]),  np.array([endData[idx[k]], yExt[k]]), clip_on=False,
                            ls='-', lw=self.lineWidth, solid_capstyle="round",
                            color=rainbow[idx[k] % len(rainbow)], alpha=self.alpha,
                            marker=markers[idx[k] % len(markers)], markeredgecolor=rainbow[idx[k] % len(rainbow)],
                            markersize=self.markerSize*(1 if markers[idx[k] % len(markers)]!='*' else 1.45),
                            markerfacecolor='none', markeredgewidth=0.7*self.lineWidth, zorder=zorder0+idx[k])
                    ax.text(xExt*1.015*((self.markerSize/7.5)**1)*normFactor, yExt[k],
                            r'%s' % (order_legends[idx[k]]), verticalalignment='bottom', fontsize=self.rightFontSize)  # verticalalignment='center'
                else:
                    ax.semilogx(np.array([xLim[1], xExt]),  np.array([endData[idx[k]], yExt[k]]), clip_on=False,
                                ls='-', lw=self.lineWidth, solid_capstyle="round",
                                color=rainbow[idx[k] % len(rainbow)], alpha=self.alpha,
                                marker=markers[idx[k] % len(markers)], markeredgecolor=rainbow[idx[k] % len(rainbow)],
                                markersize=self.markerSize*(1 if markers[idx[k] % len(markers)]!='*' else 1.45),
                                markerfacecolor='none', markeredgewidth=0.7*self.lineWidth, zorder=zorder0+idx[k])
                    ax.text(xExt * 10**(0.015*((self.markerSize/7.5)**1)*normFactor*np.log10(xExt)), yExt[k],
                            r'%s' % (order_legends[idx[k]]), verticalalignment='bottom', fontsize=self.rightFontSize)  # verticalalignment='center'
            
            # Plot the vertical separation line:
            #===================================
            if self.keepBox:
                #if not logx:
                #    ax.plot(np.array([xLim[1], xLim[1]]), np.array([yLim[0],yLim[1]]), 'k-', lw=0.7, clip_on=False, zorder=1)
                #else:
                #    ax.semilogx(np.array([xLim[1], xLim[1]]), np.array([yLim[0],yLim[1]]), 'k-', lw=0.7, clip_on=False, zorder=1)
                #ax.spines['right'].set_visible(False)
                ax.spines['right'].set_linewidth(0.7)
                
                # Plot extension part of bottom and top bars:
                if not logx:
                    ax.plot(np.array([xLim[1], xLim[1] + xExtRatio * normFactor*xLength]), np.array([yLim[0], yLim[0]]), 'k-', lw=self.axesLineWidth, clip_on=False)
                    ax.plot(np.array([xLim[1], xLim[1] + xExtRatio * normFactor*xLength]), np.array([yLim[1], yLim[1]]), 'k-', lw=self.axesLineWidth, clip_on=False)
                else:
                    ax.semilogx(np.array([xLim[1], 10 ** (np.log10(xLim[1]) + xExtRatio * normFactor*np.log10(xLength))]), np.array([yLim[0], yLim[0]]), 'k-', lw=self.axesLineWidth, clip_on=False)
                    ax.semilogx(np.array([xLim[1], 10 ** (np.log10(xLim[1]) + xExtRatio * normFactor*np.log10(xLength))]), np.array([yLim[1], yLim[1]]), 'k-', lw=self.axesLineWidth, clip_on=False)
            else:
                if not logx:
                    ax.plot(np.array([xLim[1], xLim[1]]), np.array([yLim[0], yExt[-1]+(yExt[0]-yLim[0])]), clip_on=False, ls='-', c='k', lw=0.6)
                else:
                    ax.semilogx(np.array([xLim[1], xLim[1]]), np.array([yLim[0], yExt[-1]+(yExt[0]-yLim[0])]), clip_on=False, ls='-', c='k', lw=0.6)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
            
            ax.yaxis.set_ticks_position('left') # turn of tick on right side
            # The following is not necessary anymore, thanks to: clip_on=False, self.position, and extension of bottom and top bars
            # ax.set_xlim(xLim[0], xLim[1] + xLength*xExtRatio)
            ax.set_xticks(xTicks)
            ax.set_ylim(yLim)
        else:
            print "TODO: Make traditional legends"
            pass
        
        
        # Write text annotations:
        try:
            if type(self.note) != type(list()):
                annotations = [self.note]
            else:
                annotations = self.note
            yrange = yLim[1] - yLim[0]
            for i, string in enumerate(annotations):
                if not logx:
                    ax.text(xLim[0] + 0.025*normFactor*xExt, yLim[1] - 0.035*yrange - i*(0.06*((self.annotationFontSize/15)**1)*yrange),
                            r"%s" % string, verticalalignment='top', fontsize=self.annotationFontSize,
                            bbox=dict(facecolor='white', edgecolor='none')) # verticalalignment='bottom'
                else:
                    ax.text(xLim[0] * 10**(0.025*normFactor*np.log10(xExt)), yLim[1] - 0.035*yrange - i*(0.06*((self.annotationFontSize/15)**1)*yrange),
                            r"%s" % string, verticalalignment='top', fontsize=self.annotationFontSize,
                            bbox=dict(facecolor='white', edgecolor='none')) # verticalalignment='bottom'
        except:
            pass
        
        
        #ax.xaxis.tick_bottom()
        #ax.yaxis.tick_left()
        #ax.tick_params(direction='inout')  # for both axes
        # ax.xaxis.set_tick_params(direction='inout')
        
        #[line.set_zorder(3) for line in ax.lines]
        import matplotlib as mpl
        mpl.rcParams['axes.unicode_minus'] = False
        #from matplotlib.ticker import ScalarFormatter
        #majorFormatter = ScalarFormatter(useMathText=True, useOffset=False)
        #majorFormatter.set_scientific(True)
        #majorFormatter.set_powerlimits((-4,4))
        #ax.xaxis.set_major_formatter(majorFormatter)
        #ax.yaxis.set_major_formatter(majorFormatter)
        
        # Change fontsize for x and y ticks
        ax.tick_params(axis='x', labelsize=self.xTickLabelSize)
        ax.tick_params(axis='y', labelsize=self.yTickLabelSize)
        
        ax.set_xlabel(self.xlabel, fontsize=self.bottomFontSize)
        ax.set_ylabel(self.ylabel, fontsize=self.leftFontSize)
        ax.set_title(self.title, fontsize=self.topFontSize)
        
        # Save the plot to file:
        ext = filename[-4:]
        if not ext in ['.pdf', '.eps', '.png']:
            filename = filename + '.pdf'
        fig.savefig(filename)
        #fig.savefig(filename[:-4] + '.png', dpi=(pngdpi))
        
        # plt.show()
        plt.close(fig)
    


if __name__ == "__main__":
    print "I am a module and should only be called from another program!" 
    # raw_input()
