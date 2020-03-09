from matplotlib import pyplot as plt
import pymultinest
import numpy as np

from config import parameters
from util import show

class CornerPlot:
    def __init__(self,
                 analyzer = None,
                 name_in = None,
                 input_directory = "",
                 output_directory = "",
                 ndim = 1):
        self.a = analyzer
        self.name_in = name_in
        self.input_directory = input_directory
        self.output_directory = output_directory
        self.ndim = ndim
        
    def basic_plot(self,output_name = ""):
        p = pymultinest.PlotMarginalModes(self.a)
        plt.figure(figsize=(5*self.ndim, 5*self.ndim))
        #plt.subplots_adjust(wspace=0, hspace=0)
        for i in range(self.ndim):
            plt.subplot(self.ndim, self.ndim, self.ndim * i + i + 1)
            p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=50)
            plt.ylabel("Probability")
            plt.xlabel(parameters[i])
            
            for j in range(i):
                plt.subplot(self.ndim, self.ndim, self.ndim * j + i + 1)
                #plt.subplots_adjust(left=0, bottom=0, right=0, top=0, wspace=0, hspace=0)
                p.plot_conditional(i, j, with_ellipses = False, with_points = True, grid_points=30)
                plt.xlabel(parameters[i])
                plt.ylabel(parameters[j])
                
        plt.savefig(self.output_directory + output_name + ".pdf") #, bbox_inches='tight')
        #show(self.output_directory + output_name +".pdf")
