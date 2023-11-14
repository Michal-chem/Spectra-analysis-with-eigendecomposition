
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import stats
from tkinter.filedialog import askopenfilename 
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from cycler import cycler
###########Function definitions#####################################

def read_as_array(filename):
    data = pd.read_csv(filename, delimiter = ';',decimal=',')
    data_head = list(data.columns.values)
    data_head.pop(0)
    arr = np.array(data, dtype=float)
    return arr, data_head

def clear_plot(axis, canvas):
    axis.cla()
    canvas.draw()

class SpectralAnalyze:
    
    def __init__(self, filename):
        self.filename = filename
        self.data = self.filename
        self.data_head = self.filename
        self.num_columns = self.data.shape[1]
        
    def get_plot_init(self, axis, canvas):
        clear_plot(axis, canvas)
        
        plt.tight_layout(pad=2)
        [[plt.plot(self.data[:,0], self.data[:,i])] for i in range(1, self.num_columns)]
        plt.subplots_adjust(left=0.1, right=0.85)
        plt.margins(x=0)
        plt.title('Initial spectra')
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Molar extinction coefficient [-]')
        plt.legend(data_head, loc="upper right", bbox_to_anchor=(0., 0.88, 1.17, .135),fontsize="15")
        plt.savefig("Initial spectra.pdf", format="pdf", bbox_inches="tight")
        canvas.draw()
        
    def get_plot_centr(self, axis, canvas):
        clear_plot(axis, canvas)
        centr, _ = self.get_centered(self.data)
        plt.tight_layout(pad=2)
        [[plt.plot(self.data[:,0], centr[:,i])] for i in range(0, self.num_columns-1)]
        plt.subplots_adjust(left=0.1, right=0.85)
        plt.margins(x=0)
        plt.title('Centered spectra')
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Centered molar extinction coefficient [-]')
        plt.legend(data_head, loc="upper right", bbox_to_anchor=(0., 0.88, 1.17, .135),fontsize="15")
        plt.savefig("Centered spectra.pdf", format="pdf", bbox_inches="tight")
        canvas.draw()
        
    def get_plot_autos(self, axis, canvas):
        clear_plot(axis, canvas)
        autos, _ = self.get_autoscaled(self.data)
        plt.tight_layout(pad=2)
        [[plt.plot(self.data[:,0], autos[:,i])] for i in range(1, self.num_columns)]
        plt.subplots_adjust(left=0.1, right=0.85)
        plt.margins(x=0)
        plt.title('Autoscaled spectra')
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Autoscaled molar extinction coefficient [-]')
        plt.legend(data_head, loc="upper right", bbox_to_anchor=(0., 0.88, 1.17, .135),fontsize="15")
        plt.savefig("Autoscaled spectra.pdf", format="pdf", bbox_inches="tight")
        canvas.draw()
        
    def get_plot_mean(self, axis, canvas):
        clear_plot(axis, canvas)
        plt.plot(data_head , mean_val, 'r')
        plt.subplots_adjust(left=0.15)
        plt.margins(x=0)
        plt.xlabel('pH [-]')
        plt.ylabel('Mean [-]')
        plt.savefig("mean.pdf", format="pdf", bbox_inches="tight")
        canvas.draw()
      
    def get_plot_std(self, axis, canvas):
        clear_plot(axis, canvas)
        plt.plot(data_head, std[1:], 'r')
        plt.subplots_adjust(left=0.15)
        plt.margins(x=0)
        plt.xlabel('pH [-]')
        plt.ylabel('STD [-]')
        plt.savefig("std.pdf", format="pdf", bbox_inches="tight")
        canvas.draw()

    def get_centered(self, init_spectra):
        mean_val = np.mean(init_spectra[:,1:], axis=0, dtype=np.float64)
        centr = init_spectra[:,1:]-mean_val
        return centr, mean_val
    
    def get_autoscaled(self, init_spectra):
        std = np.std(init_spectra, axis=0, dtype=np.float64)
        autos = init_spectra/std
        return autos, std

    def get_plot_corr(self, axis, canvas):
        clear_plot(axis, canvas)
        correlation_list_1 = []
        for i in range (1, self.num_columns):
            correltaion_1 = (round(stats.pearsonr(self.data[:,1], self.data[:,i])[0],4))
            correlation_list_1.append(correltaion_1)
        
        correlation_list_2 = []
        for i in range (1, self.num_columns):
            correltaion_2 = (round(stats.pearsonr(self.data[:,num_columns-1], self.data[:,i])[0],4))
            correlation_list_2.append(correltaion_2)
        
        plt.plot(data_head,correlation_list_1, 'ro-')
        plt.plot(data_head,correlation_list_2, 'bo-')
        plt.subplots_adjust(left=0.15)
        plt.margins(x=0)
        plt.xlabel('pH [-]')
        plt.ylabel(' Corr [-]')
        plt.savefig("corr.pdf", format="pdf", bbox_inches="tight")
        canvas.draw()
        
    def get_plot_comp(self):
        virek={}
        for i in range(self.num_columns-2):
            virek[i] = components[:,i]*np.transpose(evecs[:,(i)])
        plt.rcParams['axes.linewidth'] = 2   
        figure.tight_layout(pad=2.5)
        for j in range(2):
            for k in range(4):
                    axis[j,k].cla()
                    axis[j,k].margins(x=0)
                    axis[j,k].plot(spec_analysis.data[:,0],(components[:,(j*4+k)]*np.transpose(evecs[:,(j*4+k)])))
                    axis[j,k].title.set_text('Centered components ' +str((j*4+k)+1))
                    axis[j,k].set_xlabel('Wavelength [nm]')
                    axis[j,k].set_ylabel('Centered molar extinction coefficient [-]')
        plt.savefig("Components.pdf", format="pdf", dpi=600)
        canvas.draw()
        
    def get_plot_rec(self):  
        recon0 = virek[0]
        recon1 = recon0 + virek[1]
        recon2 = recon1 + virek[2]
        recon3 = recon2 + virek[3]
        recon4 = recon3 + virek[4]
        recon5 = recon4 + virek[5]
        recon6 = recon5 + virek[6]
        recon7 = recon6 + virek[7]
        recon8 = recon7 + virek[8]
        recon9 = recon8 + virek[9]
        recon = {'recon0':recon0, 'recon1':recon1, 'recon2':recon2, 'recon3':recon3, 'recon4':recon4, 'recon5':recon5,
                  'recon6':recon6, 'recon7':recon7, 'recon8':recon8, 'recon9':recon9} 
        figure.tight_layout(pad=2.5)
        for j in range(2):
            for k in range(4):
                    axis[j,k].cla()
                    axis[j,k].margins(x=0)
                    axis[j,k].plot(spec_analysis.data[:,0],(recon['recon'+str((j*4+k))]))
                    axis[j,k].title.set_text('Centered components ' +str((j*4+k)+1))
                    axis[j,k].set_xlabel('Wavelength [nm]')
                    axis[j,k].set_ylabel('Centered molar extinction coefficient [-]')
        plt.savefig("Reconstructed spectra.pdf", format="pdf", dpi=600)
        canvas.draw()
        
    def get_plot_resid(self): 
        residual={}
        for i in range(9):
            residual[i] = centered - recon['recon'+str(i)]
        figure.tight_layout(pad=2.5)
        for j in range(2):
            for k in range(4):
                    axis[j,k].cla()
                    axis[j,k].margins(x=0)
                    axis[j,k].plot(spec_analysis.data[:,0],(residual[(j*4+k)]))
                    axis[j,k].title.set_text('Centered residual spectra ' +str((j*4+k)+1))
                    axis[j,k].set_xlabel('Wavelength [nm]')
                    axis[j,k].set_ylabel('Centered molar extinction coefficient [-]')
        plt.savefig("Residuals.pdf", format="pdf", dpi=600)
        canvas.draw()
        
#########Loading data#####################################
        
filename, data_head = read_as_array(askopenfilename())
spec_analysis = SpectralAnalyze(filename)
#########Colors and initial plot settings#################################################

num_rows = np.shape(spec_analysis.data)[0]
num_columns = np.shape(spec_analysis.data)[1]

rainbow = cm.get_cmap('rainbow', num_columns-1)
newcolors = list(rainbow(np.linspace(1, 0, num_columns-1)))
plt.rc('axes', prop_cycle=(cycler('color', newcolors)))
plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8) 
plt.rc('xtick.major', size=6, width=2)
plt.rc('ytick.major', size=6, width=2)
plt.rc('lines',linewidth=1.5)
#########Initial spectra#######################################

centered, mean_val = spec_analysis.get_centered(spec_analysis.data)
autoscaled, std = spec_analysis.get_autoscaled(spec_analysis.data)

df = pd.DataFrame(data_head)

#################PCA###########################

cov_mat =np.cov(np.transpose(centered))
u, s, vh = np.linalg.svd(cov_mat)

############Components spectra##################

evecs = np.matrix(u)
components = np.matrix(centered)*evecs

virek={}
for i in range(num_columns-2):
    virek[i] = components[:,i]*np.transpose(evecs[:,(i)])

############Reconstructed spectra##################

recon0 = virek[0]
recon1 = recon0 + virek[1]
recon2 = recon1 + virek[2]
recon3 = recon2 + virek[3]
recon4 = recon3 + virek[4]
recon5 = recon4 + virek[5]
recon6 = recon5 + virek[6]
recon7 = recon6 + virek[7]
recon8 = recon7 + virek[8]
recon9 = recon8 + virek[9]
recon = {'recon0':recon0, 'recon1':recon1, 'recon2':recon2, 'recon3':recon3, 'recon4':recon4, 'recon5':recon5,
          'recon6':recon6, 'recon7':recon7, 'recon8':recon8, 'recon9':recon9}
####################GUI##################################
#Initialize    
#tk functions
def new_win():
    win2 = tk.Toplevel()
    figure2, axis2 = plt.subplots(figsize=(10,8))
    canvas2 = FigureCanvasTkAgg(figure2, master = win2)
    canvas2.get_tk_widget().pack()
    toolbar = NavigationToolbar2Tk(canvas2, win2, pack_toolbar = False)
    toolbar.update()
    toolbar.pack()
    button_frame2 = tk.Frame(win2)
    button_frame2.pack()
    tk.Button(button_frame2, text="Initial spectra", 
    command=lambda: spec_analysis.get_plot_init(axis2, canvas2)).pack(side='left', padx=10)
    tk.Button(button_frame2, text="Centered spectra", 
    command=lambda: spec_analysis.get_plot_centr(axis2, canvas2)).pack(side='left', padx=10)
    tk.Button(button_frame2, text="Autoscaled spectra", 
    command=lambda: spec_analysis.get_plot_autos(axis2, canvas2)).pack(side='left', padx=10)
    tk.Button(button_frame2, text="Mean values", 
    command=lambda: spec_analysis.get_plot_mean(axis2, canvas2)).pack(side='left', padx=10)
    tk.Button(button_frame2, text="Standard deviations", 
    command=lambda: spec_analysis.get_plot_std(axis2, canvas2)).pack(side='left', padx=10)
    tk.Button(button_frame2, text="Correlation coefficients", 
    command=lambda: spec_analysis.get_plot_corr(axis2, canvas2)).pack(side='left', padx=10)

def quit():
    window.destroy()
    window.quit()
####################################################################    
window = tk.Tk()
figure, axis = plt.subplots(2,4,figsize=(18,8))

window.title("Let's do some analysis!")

canvas = FigureCanvasTkAgg(figure, master = window)
canvas.get_tk_widget().pack()
toolbar = NavigationToolbar2Tk(canvas, window, pack_toolbar = False)
toolbar.update()
toolbar.pack()

button_frame = tk.Frame(window)
button_frame.pack()

tk.Button(button_frame, text="Components", command=spec_analysis.get_plot_comp).pack(side='left', padx=10) 
tk.Button(button_frame, text="Reconstructed spectra", command=spec_analysis.get_plot_rec).pack(side='left', padx=10)
tk.Button(button_frame, text="Residual spectra", command=spec_analysis.get_plot_resid).pack(side='left', padx=10) 

button = tk.Button(button_frame, text="One window mode", command=new_win).pack(side='left', padx=10)
button = tk.Button(button_frame, text="Quit", command=quit).pack(side='left', padx=10)

window.mainloop()
plt.close()
