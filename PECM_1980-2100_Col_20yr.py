import numpy as np
from sympy.solvers import solve
from scipy.optimize import fsolve
from sympy import Symbol
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import newton
from scipy.optimize import bisect
import six 
from six.moves import zip
import csv
from numpy import genfromtxt
import xlsxwriter as xlwt
from matplotlib.ticker import FuncFormatter

# ========================================== parameters ========================================== #
gamma = 0.5                                # Share of children's welbeing in Utility function of parents
alpha = 0.55/2                              # Agricultural share in Consumption function
eps = 0.5                                   # Elasticity of substitution in Consumption function
T = 7                                       # Time horizon
lf1980 = 65.5
lf2010 = 73.3
# ========================================== RCP Scenarios =========================================== # 
RCP = ([1000, 1000, 1000, 1000, 1000, 1000, 1000], [1000, 1000, 1000, 1063, 1060, 1034, 1007], [1000, 1000, 1000, 1139, 1266, 1318, 1335], [1000, 1000, 1000, 1113, 1286, 1518, 2476], [1000, 1000, 1000, 1271, 1665, 2183, 2785])
# ========================================== Temperature =========================================== #
# Pole temperature 0 and Equator temperature 28
# Temp = - 28/(pi/2) * (lat - pi/2) 
nu1 = 0.21                                  #in Desmet it is 0.0003 but they report it after dividing by 1000/pi/2
nu2 = 0.5
nu3 = 0.0238

# ========================================== Damages =========================================== #
# D = g0 + g1 * T + g2 * T^2

# Agricultural parameters
g0a = -2.24
g1a = 0.308
g2a = -0.0073

# Manufacturing parameters
g0m = 0.3
g1m = 0.08
g2m = -0.0023

# ========================================== Variables =========================================== #                    

# == temperature == #
Temp = np.zeros((T, 5))                     # Temperature

# == Age matrix == #
nu = np.zeros((T, 5))                       # number of unskilled children
ns = np.zeros((T, 5))                       # number of skilled children
L = np.zeros((T, 5))                        # Number of unskilled parents
H = np.zeros((T, 5))                        # Number of skilled parents
h = np.zeros((T, 5))                        # Ratio of skilled labor to unskilled labor h=H/L

R = np.zeros((T, 5))                        # Retiree population
Ah = [0] * T
Age = np.zeros((T, 5))
N = np.zeros((T, 5))                        # Adult population
Pop = np.zeros((T, 5))                      # Total population

# == Prices == #
pa = np.zeros((T, 5))                       # Pice of AgricuLtural good
pm = np.zeros((T, 5))                       # Pice of Manufacturing good
pr = np.zeros((T, 5))

# == Wages == #
wu = np.zeros((T, 5))                       # Wage of unskilled labor
ws = np.zeros((T, 5))                       # Wage of skilled labor

# == Technology == #
Aa = np.zeros((T, 5))                       # Technological growth function for Agriculture
Am = np.zeros((T, 5))                       # Technological growth function for Manufacurng
Aag = 0                                     # Technological growth rate for Agriculture
Amg = 0                                     # Technological growth rate for Manufacurng                     
Ar = np.zeros((T, 5))                       # ratio of Technology in Manufacurng to Agriculture

# == Output == #
Y = np.zeros((T, 5))                        # Total output
Ya = np.zeros((T, 5))                       # AgricuLtural output
Ym = np.zeros((T, 5))                       # Manufacturing output
Yr = np.zeros((T, 5))                       # Ratio of Manufacturing output to Agricultural output
YP = np.zeros((T, 5))                       # Output per capita
YA = np.zeros((T, 5))                       # Output per adult

# == Output == #
Da = np.zeros((T, 5))                       # AgricuLtural damage
Dm = np.zeros((T, 5))                       # Manufacturing damage
Dr = np.zeros((T, 5))                       # Ratio of Manufacturing damages to Agricultural damages

# == Consumption == #
cau = np.zeros((T, 5))                      # consumption of agricultural good unskilled
cas = np.zeros((T, 5))                      # consumption of agricultural good skilled
cmu = np.zeros((T, 5))                      # consumption of manufacturing good unskilled
cms = np.zeros((T, 5))                      # consumption of manufacturing good skilled
cu = np.zeros((T, 5))                       # consumption of all goods unskilled
cs = np.zeros((T, 5))                       # consumption of all goods skilled
   
# ============================================== Country Calibration ============================================== #
# hx: Ratio of skilled to unskilled labor in 2100 hx = Hx/Lx
# popgx: population growth rate in 2100
# nux: number of unskilled children per parent in 2100
# nsx: number of skilled children per parent in 2100
# Arx: Ratio of technology in manufacturing to technology in agriculture in 2100
# H0: skilled labor in 1980 (milion people)
# L0: unskilled labor in 1980 (milion people)
# h0: Ratio of skilled to unskilled labor in 1980 h0 = H0/L0
# nu0: number of unskilled children per parent in 1980
# ns0: number of skilled children per parent in 1980
# Am0: technology in manufacturing in 1980
# Aa0: technology in agriculture in 1980
# Ar0: Ratio of technology in manufacturing to technology in agriculture in 1980
# N0: Total labor in 1980 (milion people)
# Y0: GDP in 1980 (bilion 1990 GK$)
# C0: population of children in 1980
# r0: ratio of ts/tu
# lat, latid: latitude
# =========== COLOMBIA ============

Ndata = [7.8172, 13.0324, 16.31322, 16.9904, 15.89158, 14.21424, 11.94903]
Popdata = [26.87429, 39.76427, 52.29907, 61.46505, 65.70817, 65.35975, 61.78512]
hdata = [0.322684895, 0.754331175, 2.019781123, 5.241238961, 14.36205013, 40.39748369, 104.3056314]
#N00 = 5.20
Y0 = 104.112
latid = 45

### ============================================== Model Calibration ============================================== #
hx = hdata[T-1]
popgx = (Ndata[T-1] - Ndata[T-2])/Ndata[T-2]
N1 = Ndata[1]
N0 = Ndata[0]

nux = (1 + popgx) / (1 + hx)
nsx = (1 + popgx) * hx / (1 + hx)

lat = latid * np.pi/180
Temp0 = - 28/(np.pi/2) * (lat - np.pi/2)
Da0 = max(0.001, g0a + g1a * Temp0 + g2a * Temp0**2)
Dm0 = max(0.001, g0m + g1m * Temp0 + g2m * Temp0**2)
Dr0 = Dm0/Da0

h0 = hdata[0]
h1 = hdata[1]
popg0 = (N1 - N0)/N0

nu0 = (1 + popg0) / (1 + h1)
ns0 = (1 + popg0) * h1 / (1 + h1)

tu = (nsx - ns0)/(nsx * nu0 - ns0 * nux) * gamma
ts = (nux - nu0)/(nux * ns0 - nu0 * nsx) * gamma
r0 = ts/tu

L0 = nu0 * N0
H0 = ns0 * N0

Ar0 = np.exp((eps * np.log((1 - alpha)/alpha) - eps * np.log(r0) - np.log(h1)) /(1 - eps) - np.log(Dr0))
Am0 = Y0/((alpha * (L0 / Ar0)**((eps - 1)/eps) + (1 - alpha) * H0**((eps - 1)/eps))**(eps/(eps - 1)))
Aa0 = Am0/Ar0

Arx = np.exp((eps * np.log((1 - alpha)/alpha) - eps * np.log(r0) - np.log(hx)) /(1 - eps) - np.log(Dr0))
Arg = (Arx/Ar0)**(1/((2100 - 2000)/20)) - 1

Amg = (1 + 0.02)**20 - 1
Aag = (1 + Amg)/(1 + Arg) - 1

Ya0 = Aa0 * L0 * Da0
Ym0 = Am0 * H0 * Dm0
Yr0 = Ym0 / Ya0
Y0 = (alpha * Ya0**((eps -1)/eps) + (1 - alpha) * Ym0**((eps -1)/eps))**(eps/(eps - 1))


################################# Life Expectancy #############################

Ahr = (lf2010/lf1980) - 1

# RCP 2.6 
Con1 = RCP[1][2] - RCP[0][0]
Temp1 = Temp0 + nu1 * Con1**nu2 * (1 - nu3 * Temp0)
Age1 = lf2010 + 3 * 3

# RCP 4.5 
Con2 = RCP[2][2] - RCP[0][0]
Temp2 = Temp0 + nu1 * Con2**nu2 * (1 - nu3 * Temp0)
Age2 = lf2010 + 2 * 3

theta = -np.log(Age1/Age2)/np.log(Temp1/Temp2)
Ah0 = np.log(lf1980/Temp0**(-theta))
Ah[0] = np.exp(Ah0)
Age[0, :] = lf1980

#######################################################################

R0 = N0 * (Age[0, 0] - 65)/25

pr0 = (Yr0)**(-1/eps) * (1 - alpha) / alpha

ca = Ya0 / N0
cmu0 = Ym0 / (H0 * r0 + L0)
cms0 = cmu0 * r0
cau0 = Ya0 / (H0 * r0 + L0)
cas0 = cau0 * r0    
cu0 = (alpha * cau0**((eps - 1)/eps) + (1 - alpha) * cmu0**((eps - 1)/eps))**(eps/(eps - 1))
cs0 = (alpha * cas0**((eps - 1)/eps) + (1 - alpha) * cms0**((eps - 1)/eps))**(eps/(eps - 1))
wu0 = cu0 / (1 - gamma)
ws0 = cs0 / (1 - gamma)
pa0 = wu0 / (Da0 * Aa0)
pm0 = ws0 / (Dm0 * Am0)

# ============================================== Transition Function ============================================== #

for j in range(5):
    N[0, j] = N0
    R[0, j] = R0
    Temp[0, j] = Temp0
    h[0, j] = hdata[0]
    Da[0, j] = Da0
    Dm[0, j] = Dm0
    Dr[0, j] = Dr0
    H[0, j] = H0
    L[0, j] = L0
    Y[0, j] = Y0
    YA[0, j] = Y0/N0
    Ya[0, j] = Ya0
    Ym[0, j] = Ym0
    Yr[0, j] = Yr0
    Aa[0, j] = Aa0
    Am[0, j] = Am0
    cmu[0, j] = cmu0
    cms[0, j] = cms0
    cau[0, j] = cau0
    cas[0, j] = cas0
    cu[0, j] = cu0
    cs[0, j] = cs0
    wu[0, j] = wu0
    ws[0, j] = ws0
    pa[0, j] = pa0
    pm[0, j] = pm0
    
    for i in range(T - 1):
        Con = RCP[j][i + 1] - RCP[j][0]
        Temp[i + 1, j] = Temp0 + nu1 * Con**nu2 * (1 - nu3 * Temp0)
        Da[i + 1, j] = max(0.001, g0a + g1a * Temp[i + 1, j] + g2a * Temp[i + 1, j]**2)
        Dm[i + 1, j] = max(0.001, g0m + g1m * Temp[i + 1, j] + g2m * Temp[i + 1, j]**2)
        Dr[i + 1, j] = Dm[i + 1, j]/Da[i + 1, j]
        
        Aa[i + 1, j] = Aa0 * (1 + Aag)**i
        Am[i + 1, j] = Am0 * (1 + Amg)**i
        Ar[i + 1, j] = Am[i + 1, j]/Aa[i + 1, j]
        
        h[i + 1, j] = np.exp(eps * (np.log((1 - alpha)/alpha) - np.log(r0) - (1 - eps)/eps * np.log(Dr[i + 1, j] * Ar[i + 1, j])))
        
        nu[i, j] = gamma / (h[i + 1, j] * ts + tu)
        ns[i, j] = nu[i, j] * h[i + 1, j]
      
        H[i + 1, j] = ns[i, j] * N[i, j]
        L[i + 1, j] = nu[i, j] * N[i, j]
        N[i + 1, j] = L[i + 1, j] + H[i + 1, j]
        Ya[i + 1, j] = Aa[i + 1, j] * L[i + 1, j] * Da[i + 1, j]
        Ym[i + 1, j] = Am[i + 1, j] * H[i + 1, j] * Dm[i + 1, j]
        Yr[i + 1, j] = Ym[i + 1, j] / Ya[i + 1, j]
        pr[i + 1, j] = (Yr[i + 1, j])**(-1/eps) * (1 - alpha) / alpha
        ca = Ya[i + 1, j] / N[i + 1, j]
        cmu[i + 1, j] = Ym[i + 1, j] / (H[i + 1, j] * r0 + L[i + 1, j])
        cms[i + 1, j] = cmu[i + 1, j] * r0
        cau[i + 1, j] = Ya[i + 1, j] / (H[i + 1, j] * r0 + L[i + 1, j])
        cas[i + 1, j] = cau[i + 1, j] * r0    
        cu[i + 1, j] = (alpha * cau[i + 1, j]**((eps - 1)/eps) + (1 - alpha) * cmu[i + 1, j]**((eps - 1)/eps))**(eps/(eps - 1))
        cs[i + 1, j] = (alpha * cas[i + 1, j]**((eps - 1)/eps) + (1 - alpha) * cms[i + 1, j]**((eps - 1)/eps))**(eps/(eps - 1))
        wu[i + 1, j] = cu[i + 1, j] / (1 - gamma)
        ws[i + 1, j] = cs[i + 1, j] / (1 - gamma)
        pa[i + 1, j] = wu[i + 1, j] / (Da[i + 1, j] * Aa[i + 1, j])
        pm[i + 1] = ws[i + 1, j] / (Dm[i + 1, j] * Am[i + 1, j])    
        Y[i + 1, j] = (alpha * Ya[i + 1, j]**((eps -1)/eps) + (1 - alpha) * Ym[i + 1, j]**((eps -1)/eps))**(eps/(eps - 1))
        YA[i + 1, j] = Y[i + 1, j] / N[i + 1, j]

        Ah[i + 1] = Ah[i] * (1 + Ahr)
        Age[i + 1, j] = max(65, min(90, Ah[i + 1] * Temp[i + 1, j]**(-theta)))
        R[i + 1, j] = N[i, j] * (Age[i + 1, j] - 65)/25
        
        Pop[i, j] = N[i + 1, j] + N[i, j] + R[i, j]
        YP[i, j] = Y[i, j] / Pop[i, j]
        
# ===================================================== Output ===================================================== #    
x = [i for i in range(1980, 2120, 20)]

plt.plot(x, Ndata, 'r:', label = "Data")
plt.plot(x, N[:, 0], 'blue', label = "Baseline")
plt.plot(x, N[:, 1], 'green', label = "RCP 2.6")
plt.plot(x, N[:, 2], 'cyan', label = "RCP 4.5")
plt.plot(x, N[:, 3], 'orange', label = "RCP 6.0")
plt.plot(x, N[:, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('millions')
plt.title('Adult PopuLation')
axes = plt.gca()
plt.xticks(np.arange(min(x), max(x)+1, 20))
plt.legend(loc=2, prop={'size':8})
plt.show()

plt.plot(x[0:T - 1], Popdata[0:T - 1], 'r:', label = "Data")
plt.plot(x[0:T - 1], Pop[0:T - 1, 0], 'blue', label = "Baseline")
plt.plot(x[0:T - 1], Pop[0:T - 1, 1], 'green', label = "RCP 2.6")
plt.plot(x[0:T - 1], Pop[0:T - 1, 2], 'cyan', label = "RCP 4.5")
plt.plot(x[0:T - 1], Pop[0:T - 1, 3], 'orange', label = "RCP 6.0")
plt.plot(x[0:T - 1], Pop[0:T - 1, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('millions')
plt.title('Total PopuLation')
axes = plt.gca()
plt.xticks(np.arange(min(x), max(x) -20 +1, 20))
plt.legend(loc=2, prop={'size':8})
plt.show()

plt.plot(x, hdata, 'r:', label = "Data")
plt.plot(x, h[:, 0], 'blue', label = "Baseline")
plt.plot(x, h[:, 1], 'green', label = "RCP 2.6")
plt.plot(x, h[:, 2], 'cyan', label = "RCP 4.5")
plt.plot(x, h[:, 3], 'orange', label = "RCP 6.0")
plt.plot(x, h[:, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('Ratio')
plt.title('Ratio of skilled to unskilled adults')
axes = plt.gca()
plt.xticks(np.arange(min(x), max(x)+1,20))
plt.legend(loc=2, prop={'size':8})
plt.show()

plt.plot(x, Age[:, 0], 'blue', label = "Baseline")
plt.plot(x, Age[:, 1], 'green', label = "RCP 2.6")
plt.plot(x, Age[:, 2], 'cyan', label = "RCP 4.5")
plt.plot(x, Age[:, 3], 'orange', label = "RCP 6.0")
plt.plot(x, Age[:, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('years')
plt.title('Life expectancy')
axes = plt.gca()
plt.xticks(np.arange(min(x), max(x)+1,20))
plt.legend(loc=2, prop={'size':8})
plt.show()

plt.plot(x, Dr[:, 0], 'b--', label = "Baseline")
plt.plot(x, Dr[:, 1], 'green', label = "RCP 2.6")
plt.plot(x, Dr[:, 2], 'cyan', label = "RCP 4.5")
plt.plot(x, Dr[:, 3], 'orange', label = "RCP 6.0")
plt.plot(x, Dr[:, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('Ratio')
plt.title('Ratio of Manufacturing to Agricultural Damages')
axes = plt.gca()
plt.xticks(np.arange(min(x), max(x)+1, 20))
plt.legend(loc=2, prop={'size':8})
plt.show()
#
#plt.plot(x[0:T - 1], YP[0:T - 1, 0], 'b--', label = "Baseline")
#plt.plot(x[0:T - 1], YP[0:T - 1, 1], 'green', label = "RCP 2.6")
#plt.plot(x[0:T - 1], YP[0:T - 1, 2], 'cyan', label = "RCP 4.5")
#plt.plot(x[0:T - 1], YP[0:T - 1, 3], 'orange', label = "RCP 6.0")
#plt.plot(x[0:T - 1], YP[0:T - 1, 4], 'brown', label = "RCP 8.5")
#plt.xlabel('Time')
#plt.ylabel('1990 GK$/person')
#plt.title('Output per capita')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), xend - 30, 30))
#plt.legend(loc=2, prop={'size':8})
#plt.show()

#plt.plot(x[0:T - 1], YP[0:T - 1], 'darkgreen')
#plt.xlabel('Time')
#plt.ylabel('1000 GK$ per person')
#plt.title('GDP per capita')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), max(x), 20))
#plt.show()

#plt.plot(x, Dr, 'brown')
#plt.xlabel('Time')
#plt.ylabel('Ratio')
#plt.title('Damage ratio (manufacturing to agricultural)')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), xend, 20))
#plt.show()
#
#plt.plot(x[0:T - 1], nu[0:T - 1], 'c')
#plt.xlabel('Time')
#plt.ylabel('Population')
#plt.title('Number of unskilled children per parent')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), max(x), 20))
#plt.show()
#
#plt.plot(x[0:T - 1], ns[0:T - 1], 'g')
#plt.xlabel('Time')
#plt.ylabel('Population')
#plt.title('Number of skilled children per parent')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), max(x), 20))
#plt.show()
#
#plt.plot(x, L, 'r')
#plt.xlabel('Time')
#plt.ylabel('Population (millions)')
#plt.title('Unskilled labor')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), xend, 20))
#plt.show()
#
#plt.plot(x, H, 'y')
#plt.xlabel('Time')
#plt.ylabel('Population (millions)')
#plt.title('Skilled labor')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), xend, 20))
#plt.show()
#
#plt.plot(x, pr, 'brown')
#plt.xlabel('Time')
#plt.ylabel('Ratio')
#plt.title('Price ratio (manufacturing to agricultural)')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), xend, 20))
#plt.show()
#
#plt.plot(x, Y, 'darkcyan')
#plt.xlabel('Time')
#plt.ylabel('GK$')
#plt.title('GDP (1990 billion GK$)')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), xend, 20))
#plt.show()
#
#plt.plot(x, Ar, 'chocolate')
#plt.xlabel('Time')
#plt.ylabel('level')
#plt.title('Technology ratio (Manufacturing to Agriculture)')
#axes = plt.gca()
#plt.xticks(np.arange(min(x), xend, 20))
#plt.show()