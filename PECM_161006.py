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
gamma = 0.9                                # Share of children's welbeing in Utility function of parents
alpha = 0.55/2                              # Agricultural share in Consumption function
eps = 0.5                                   # Elasticity of substitution in Consumption function
T = 6                                       # Time horizon
# ========================================== RCP Scenarios =========================================== # 
RCP = ([882, 882, 882, 882, 882, 882], [882, 998, 1063, 1060, 1034, 1007], [882, 1001, 1139, 1266, 1318, 1335], [882, 995, 1113, 1286, 1518, 2476], [882, 1049, 1271, 1665, 2183, 2785])
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
# =========== COLOMBIA ============ #4.5709,25,45,65
Ndata = [26.73, 38.17, 44.33, 44.88, 43.71, 39.43]
hdata = [0.49, 0.95, 1.89, 3.88, 6.91, 16.40]

H00 = 2.693
L00 = 13.198
Y00 = 114.667
latid = 4.5709

### ============================================== Model Calibration ============================================== #
hx = hdata[5]
popgx = (Ndata[5] - Ndata[4])/Ndata[4]
N0 = Ndata[0]

nux = (1 + popgx) / (1 + hx)
nsx = (1 + popgx) * hx / (1 + hx)

N00 = H00 + L00
h00 = H00/L00

lat = latid * np.pi/180
Temp0 = - 28/(np.pi/2) * (lat - np.pi/2)
Da0 = max(0.001, g0a + g1a * Temp0 + g2a * Temp0**2)
Dm0 = max(0.001, g0m + g1m * Temp0 + g2m * Temp0**2)
Dr0 = Dm0/Da0

def hsolve(rx):
    Arx = np.exp((eps * np.log((1 - alpha)/alpha) - eps * np.log(rx) - np.log(hx)) /(1 - eps) - np.log(Dr0))
    Ar0 = np.exp((eps * np.log((1 - alpha)/alpha) - eps * np.log(rx) - np.log(h00)) /(1 - eps) - np.log(Dr0))
    Argx = np.exp((np.log(Arx/Ar0))/((2100 - 1980)/20)) - 1
    h1 = np.exp(eps * (np.log((1 - alpha)/alpha) - np.log(rx) - (1 - eps)/eps * np.log((Dr0) * Ar0 * (1 + Argx))))
    tu1x = gamma * (1 + h1) / (h1 * rx + 1) * N00/N0
    tu2x = gamma / (rx * nsx + nux)
    return tu1x - tu2x

r0 = bisect(hsolve, 100, 0.01)
tu = gamma / (r0 * nsx + nux)
ts = tu * r0

Ar00 = np.exp((eps * np.log((1 - alpha)/alpha) - eps * np.log(r0) - np.log(h00)) /(1 - eps) - np.log(Dr0))
Am00 = Y00/((alpha * (L00 / Ar00)**((eps - 1)/eps) + (1 - alpha) * H00**((eps - 1)/eps))**(eps/(eps - 1)))
Aa00 = Am00/Ar00

Arx = np.exp((eps * np.log((1 - alpha)/alpha) - eps * np.log(r0) - np.log(hx)) /(1 - eps) - np.log(Dr0))
Arg = np.exp((np.log(Arx/Ar00))/((2100 - 1980)/20)) - 1

Amg = (1 + 0.02)**20 - 1
Aag = (1 + Amg)/(1 + Arg) - 1

Ar0 = Ar00 * (1 + Arg)
Aa0 = Aa00 * (1 + Aag)
Am0 = Am00 * (1 + Amg)

h0 = np.exp(eps * (np.log((1 - alpha)/alpha) - np.log(r0) - (1 - eps)/eps * np.log(Dr0 * Ar0)))

nu0 = gamma / (h0 * ts + tu)
ns0 = nu0 * h0

L0 = N0 / (1 + h0)
H0 = L0 * h0

Ya0 = Aa0 * L0 * Da0
Ym0 = Am0 * H0 * Dm0
Yr0 = Ym0 / Ya0
Y0 = (alpha * Ya0**((eps -1)/eps) + (1 - alpha) * Ym0**((eps -1)/eps))**(eps/(eps - 1))

pr0 = (Yr0)**(-1/eps) * alpha / (1 - alpha)

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
    Temp[0, j] = Temp0
    h[0, j] = h0
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
        
        Aa[i + 1, j] = Aa[i, j] * (1 + Aag)
        Am[i + 1, j] = Am[i, j] * (1 + Amg)
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
        pr[i + 1, j] = (Yr[i + 1, j])**(-1/eps) * alpha / (1 - alpha)
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
        
        Pop[i, j] = N[i + 1, j] + N[i, j]
        YP[i, j] = Y[i, j] / Pop[i, j]
        
# ===================================================== Output ===================================================== #    
x = [2000, 2020, 2040, 2060, 2080, 2100]
xend = 2120

plt.plot(x, Ndata, 'r:', label = "Data")
plt.plot(x, N[:, 0], 'blue', label = "Baseline")
plt.xlabel('Time')
plt.ylabel('millions')
plt.title('Adult PopuLation')
axes = plt.gca()
plt.xticks(np.arange(min(x), max(x), 20))
plt.legend(loc=2, prop={'size':8})
plt.show()

plt.plot(x, N[:, 0], 'b--', label = "Baseline")
plt.plot(x, N[:, 1], 'green', label = "RCP 2.6")
plt.plot(x, N[:, 2], 'cyan', label = "RCP 4.5")
plt.plot(x, N[:, 3], 'orange', label = "RCP 6.0")
plt.plot(x, N[:, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('Percentage Change')
plt.title('PopuLation')
axes = plt.gca()
axes.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y))) 
plt.xticks(np.arange(min(x), max(x), 20))
plt.legend(loc=2, prop={'size':8})
plt.show()

plt.plot(x, hdata, 'k:', label = "Data")
plt.plot(x, h[:, 0], 'blue', label = "Baseline")
plt.plot(x, h[:, 1], 'green', label = "RCP 2.6")
plt.plot(x, h[:, 2], 'cyan', label = "RCP 4.5")
plt.plot(x, h[:, 3], 'orange', label = "RCP 6.0")
plt.plot(x, h[:, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('Ratio')
plt.title('Ratio of skilled to unskilled adults')
axes = plt.gca()
plt.xticks(np.arange(min(x), xend, 20))
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
plt.xticks(np.arange(min(x), xend, 20))
plt.legend(loc=2, prop={'size':8})
plt.show()

plt.plot(x[0:T - 1], YP[0:T - 1, 0], 'b--', label = "Baseline")
plt.plot(x[0:T - 1], YP[0:T - 1, 1], 'green', label = "RCP 2.6")
plt.plot(x[0:T - 1], YP[0:T - 1, 2], 'cyan', label = "RCP 4.5")
plt.plot(x[0:T - 1], YP[0:T - 1, 3], 'orange', label = "RCP 6.0")
plt.plot(x[0:T - 1], YP[0:T - 1, 4], 'brown', label = "RCP 8.5")
plt.xlabel('Time')
plt.ylabel('1990 GK$/person')
plt.title('Output per capita')
axes = plt.gca()
plt.xticks(np.arange(min(x), xend - 20, 20))
plt.legend(loc=2, prop={'size':8})
plt.show()

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