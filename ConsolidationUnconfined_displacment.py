import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import brewer2mpl as brew
import math
import sympy as sy
import scipy as sp
#from scipy import integrate

##############################################
data = pd.read_csv('num_consolidation.csv')
ana_x = np.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0])
ana_u = np.array([-1.1347e-18,-3.5531e-7,-9.2513e-7,-2.0187e-6,-4.1593e-6,-8.2537e-6,-1.5837e-5,-2.9411e-5,-6.1052e-5,-9.2168e-5,-0.00015568,-0.00025504,-0.00040547,-0.00062602,-0.00093928,-0.0013707,-0.0019415,-0.002695,-0.0036382,-0.0047959,-0.0061804])

###############################################
def setupMatplotlib(height=8., width=6.):
    set1 = brew.get_map('Set1', 'Qualitative', 8).mpl_colors
    dark = brew.get_map('Dark2', 'Qualitative', 8).mpl_colors
    paired = brew.get_map('Paired', 'Qualitative', 12).mpl_colors
    reds = brew.get_map('Reds', 'Sequential', 8).mpl_colors
    blues = brew.get_map('Blues', 'Sequential', 9, reverse='True').mpl_colors
    spec = brew.get_map('Spectral', 'Diverging', 8).mpl_colors
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.color'] = 'black'
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.formatter.limits'] = -3, 3
    plt.rcParams['axes.color_cycle'] = set1
    # For ipython notebook display set default values.
    # plt.rcParams['lines.markersize'] = 12
    plt.rcParams['figure.figsize'] = (height, width)
    plt.rcParams['grid.linewidth'] = 1

    # General settings used by display and print contexts.
    plt.rcParams['axes.axisbelow'] = True
    grid_line_color = '0.5'
    plt.rcParams['grid.color'] = grid_line_color
    plt.rcParams['grid.linestyle'] = '-'


###############################################
def commonFormat(ax_el: object, centerx: object = None, centery: object = None) -> object:
    # ax_el.set_xlim(0,0.08)
    # ax_el.grid(True)
    # nur einfache Achsen, kein Rahmen
    ax_el.spines['top'].set_visible(0)
    ax_el.spines['right'].set_visible(0)
    ax_el.spines['bottom'].set_linewidth(0.5)
    ax_el.spines['left'].set_linewidth(0.5)
    ax_el.xaxis.set_ticks_position('bottom')
    ax_el.yaxis.set_ticks_position('left')
    if ((centerx is not None) and (centery is not None)):
        ax_el.spines['left'].set_position(('data', centerx))
        ax_el.spines['bottom'].set_position(('data', centery))
        ax_el.spines['right'].set_position(('data', centerx - 1))
        ax_el.spines['top'].set_position(('data', centery - 1))
    # Shink current axis's height by 10% on the bottom
    ax_el.legend(loc='upper right')
    # box = ax_el.get_position()
    # ax_el.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.9])
    # Put a legend below current axis
    # ax_el.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)





############################################
def plot_temp_ana(name):
    # parameters
    eta = 0.001            # viscosity
    E = 30000              # Young's modulus  
    v = 0.2                # poisson's ratio
    mu_s = E/(2+2*v)   #lame coefficient
    lambda_s = 2*mu_s*v/(1-2*v)  # lame coefficient
    HA = lambda_s + 2*mu_s   # aggregate modulus
    
    

    #np.savetxt('test',R)
    setupMatplotlib()
    plt.close('all')
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    # ax.plot(r, s_tt, label='analytical')
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    lns1 = ax.plot(data.arc_length, data.displacement1, "o", mec="red", mfc="none", label='ogs6')
    lns2 = ax.plot(ana_x, ana_u, "b--", label='analytical')
    ax.set_xlabel('length_z (m)')
    ax.set_ylabel('Displacment_z (K)', color = 'b')
    legend = ax.legend(loc='upper right', shadow=True)
    plt.title('Confined Consolidation')
    ax.grid()
  
    '''fig = plt.figure(figsize=(6, 4))
    # fig, axes = plt.subplots(nrows=2, ncols=2)
    ax = fig.add_subplot(111)
    ay = ax.twinx()
    ax1 = fig.add_subplot(221)
    # fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_title('2 months')
    #lns1 = ax1.plot(data1_2m.arc_length, data1_2m.TEMPERATURE1, "b", label='ogs5')
    ax1.plot(data2_2m.arc_length, T2months_num, "o", mec="red", mfc="none", label='ogs6')
    ax1.plot(x, T2months, "b--", label='analytical')
    ax2.plot(x, R1, "black", label='error')


    ax3 = fig.add_subplot(222)
    ax4 = ax3.twinx()
    ax3.set_title('1 year')
    #ax3.plot(data1_1y.arc_length, data1_1y.TEMPERATURE1, "b", label='ogs5')
    ax3.plot(data2_1y.arc_length, T1year_num, "o", mec="red", mfc="none", label='ogs6')
    ax3.plot(x, T1year, "g--", label='analytical')
    ax4.plot(x, R2, "black", label='error')
    ax5 = fig.add_subplot(223)
    ax6 = ax5.twinx()
    ax5.set_title('2 years')
    #ax5.plot(data1_2y.arc_length, data1_2y.TEMPERATURE1, "b", label='ogs5')
    ax5.plot(data2_2y.arc_length, T2years_num, "o", mec="red", mfc="none", label='ogs6')
    ax5.plot(x, T2years, "g--", label='analytical')
    ax6.plot(x, R3, "black", label='error')
    ax7 = fig.add_subplot(224)
    ax8 = ax7.twinx()
    ax7.set_title('4 years')
    #ax7.plot(data1_4y.arc_length, data1_4y.TEMPERATURE1, "b", label='ogs5')
    ax7.plot(data2_4y.arc_length, T4years_num, "o", mec="red", mfc="none", label='ogs6')
    ax7.plot(x, T4years, "g--", label='analytical')
    ax8.plot(x, R4, "black", label='error')

    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ay.spines['top'].set_color('none')
    ay.spines['bottom'].set_color('none')
    ay.spines['left'].set_color('none')
    ay.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ay.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Temperature (K)')
    ay.set_ylabel('Error')
    ax1.legend()
    ax3.legend()
    ax5.legend()
    ax7.legend()
    # ax1.set_title('2 months')
    # commonFormat(ax)
    # ax.set_xlabel('Distance', fontsize=16)
    # ax.xaxis.set_label_coords(0.5,-0.05)
    # ax.set_ylabel('Temperature', fontsize=18)
    # ax.set_xlim(right=10)
    # ax.set_ylim(bottom=3.e-2)
    # ax.set_ylim(top = 7.e-2)
    # ax.set_yscale('log')
    # ax.legend(loc='lower right')
    # ax.yaxis.set_label_coords(-0.05, 0.5) '''
    fig.tight_layout()
    fig.savefig(name)
    return None


print(
    "The file result_OGS.csv is created by exporting the time series data of one of the upper nodes in the simple heat transport from paraview.")

plot_temp_ana('validation.pdf')
