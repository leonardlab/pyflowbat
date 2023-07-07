from matplotlib import cycler
import numpy as np

std_pfb_style = {
    'font.family': ['Nunito', 'Lato', 'Arial'],

    'lines.linewidth': 4,
    'lines.solid_capstyle': 'butt',

    'legend.fancybox': True,

    'axes.prop_cycle': cycler('color', ['0267c1', 'efa00b', '00af54', 'e2cfea',  'd65108', '6915E0', '3DA5D9', 'FEC601', '097109', 'EA7317', 'AF0BA5']),
    'axes.facecolor': 'white',
    'axes.labelsize': 'large',
    'axes.axisbelow': True,
    'axes.grid': True,
    'axes.edgecolor': 'B0B0B0',
    'axes.linewidth': 3.0,
    'axes.titlesize': 'x-large',
    'axes.spines.left':   False,
    'axes.spines.bottom': False,
    'axes.spines.top':    False,
    'axes.spines.right':  False,

    'patch.edgecolor': 'white',
    'patch.linewidth': 0.5,

    'svg.fonttype': 'path',

    'grid.linestyle': ':',
    'grid.linewidth': 1.0,
    'grid.color': '#B0B0B0',

    'xtick.major.size': 0,
    'xtick.minor.size': 0,
    'ytick.major.size': 0,
    'ytick.minor.size': 0,

    'font.size':14.0,

    'savefig.edgecolor': 'white',
    'savefig.facecolor': 'white',

    'figure.subplot.left': 0.08,
    'figure.subplot.right': 0.95,
    'figure.subplot.bottom': 0.07,
    'figure.facecolor': 'white'
}

std_lims = {
    'FSC-A': [0, 250000],
    'SSC-A': [0, 250000],
    'FSC-H': [0, 150000]
}

# standard bead conversions from Spherotech RCP-30-5 Rainbow Calibration Beads
std_beads_conversions = {'MEFLs': np.array([2.96174193e+02, 4.16526369e+03, 1.38208041e+04, 4.15681797e+04,
       1.32189078e+05, 4.01564084e+05, 1.09439899e+06, 3.37475176e+06,
       8.42886081e+06]), 'MEPEs': np.array([2.55346084e+02, 6.18159397e+03, 2.08165858e+04, 6.08673939e+04,
       1.92590235e+05, 5.77860955e+05, 1.48562109e+06, 4.51508525e+06,
       1.20203770e+07]), 'MEPTRs': np.array([6.15257213e+02, 1.41420304e+04, 5.13220397e+04, 1.30837663e+05,
       4.02586642e+05, 1.14348797e+06, 2.75212427e+06, 7.98238992e+06,
       2.03980478e+07]), 'MEPCYs': np.array([4.36859676e+01, 1.44240949e+03, 6.78612600e+03, 2.76897057e+04,
       1.27998622e+05, 5.67522412e+05, 2.10264642e+06, 1.02166094e+07,
       5.71324419e+07]), 'MEPCY7s': np.array([2.07469132e+03, 4.42095781e+03, 1.00367893e+04, 2.78049020e+04,
       8.35203355e+04, 2.85550403e+05, 9.24370467e+05, 4.41855836e+06,
       2.56252821e+07]), 'MEAPs': np.array([1.98698849e+02, 9.59618067e+02, 1.70798410e+03, 6.16525488e+03,
       1.49365050e+04, 5.63634150e+04, 1.58226231e+05, 3.95707136e+05,
       9.23106927e+05]), 'MEA750s': np.array([5.60949820e+01, 1.26011070e+03, 2.55370504e+03, 1.09025476e+04,
       2.87236563e+04, 1.22695407e+05, 3.69795325e+05, 1.00850308e+06,
       2.79725808e+06]), 'MECSBs': np.array([6.57121454e+01, 9.08748429e+02, 2.81780205e+03, 7.60638456e+03,
       2.20734537e+04, 5.97992299e+04, 2.35662066e+05, 6.65551507e+05,
       1.56299483e+06]), 'MEBFPs': np.array([4.24492186e+02, 2.26445309e+03, 6.73356158e+03, 1.88097941e+04,
       6.04525847e+04, 1.86755351e+05, 2.91316552e+05, 9.22856088e+05,
       2.79197334e+06])}
