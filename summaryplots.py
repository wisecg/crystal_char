import numpy as np
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys

rows_to_skip= [1,55,78] #only if there is incomplete entry
df = pd.read_csv("PATH",skiprows= rows_to_skip,usecols=(4,5,6,7,8,9,10,18))
#csv file downloaded from ELOG using 'Find'

'''
Plot the histograms
'''

'''Histogram #1: plot histogram of voltages when log(gain)=4, interpolated from gain curve'''
#data are all strings, need to convert column "Gain Offset" to a float array without uncertainty values
def offset_cleanup():
    offset=df['Gain Offset']
    offsetPar = np.zeros(len(offset)) #create an empty set to store data in for loop
    for i in range(len(offset)):
        line = offset[i]
        if isinstance(line,str):  #numbers are stored as string, "NaN"'s are float.
            count = line.index('+') #find index of '+'
            offsetPar[i] = line[:count]  #only keep digits before the found index, to get rid of uncertainty value
        else:
            offsetPar[i] = offset[i]
        i+=1
    return offsetPar
    #now offsetPar list contains no uncertainty value

def slope_cleanup():
    slope=df['Gain Slope']
    slopePar = np.zeros(len(slope))
    i = 0
    for i in range(len(slope)):
        line = slope[i]
        if isinstance(line,str):
            if len(line)>=8:
                count = line.index('+')
                slopePar[i] = line[:count]
            else:   #some don't have "+", so we can directly use the data
                slopePar[i] = line[:8]
        else:
            slopePar[i] = slope[i]
        i+=1
    return slopePar

    #now slopePar list contains no  uncertainty value

def curvature_cleanup():
    sat=df['Gain Saturation']
    curvaturePar = np.zeros(len(sat))
    for i in range(len(sat)):
        line = sat[i]
        if isinstance(line,str):
            count = line.index('e')+4
            curvaturePar[i] = line[:count]
        else:
            curvaturePar[i] = sat[i]
        i+=1
    return curvaturePar
    #now curvaturePar list contains no uncertainty value

#to find voltage when log(gain)=4
def get_voltages(gain):   #function output depends on gain value
    voltList = []
    for i in range(df.shape[0]):
        curvature = curvature_cleanup()
        slope = slope_cleanup()
        offset = offset_cleanup()
        a= curvature[i]
        b=slope[i]
        c=offset[i] - gain#gain is log(gain)
        coeff = [a, b, c]
        roots = np.roots(coeff)
        voltList.append(roots[1])
    return voltList  #a list of interpolated voltage values for all crystals is created

#plot the voltage value, at log(gain=4), as a histogram.
#Crytals with similar voltage values will be grouped together in experiment later on.
def make_volt_hist():
    volt_data = get_voltages(4)
    plt.hist(volt_data, bins=range(620,950,10), stacked=True, facecolor='blue',ec='black')
    plt.title("Voltages (Log(G)=4)")
    plt.xlabel('Voltage')
    plt.xticks(range(620,950,50))    #Range can be changed
    plt.ylabel('Counts')
    plt.savefig('PATH')
    plt.show()



'''Histogram #2:137Cs Position 3 Peak Resolution'''
def res_cleanup():
    res=df['137Cs Position 3 Peak Resolution']
    resPar= np.zeros (len(res))
    for i in range(len(res)):     #getting rid of uncertainty values
        line = res[i]
        if isinstance(line,str):
            count = line.index('+')
            resPar[i] = line[:count]
        else:
            resPar[i] = res[i]
        i+=1
    resList = resPar.astype(np.float)  #convert string array to float array
    return resList

#this function makes plot#2, a histogram of energy resolution
def make_res_hist():
    resList=res_cleanup()
    plt.hist(resList, bins=range(20,46,2), stacked=True, facecolor='blue',ec='black')
    plt.title("137Cs Energy Resolution")
    plt.xlabel('Calibrated Energy (keV)')
    plt.xticks(range(20,46,2))
    plt.ylabel('Counts')
    plt.savefig('PATH')
    plt.show()



'''Histogram #3:137Cs Total Energy Variation'''
def var_cleanup():
    var=df['137Cs Total Energy Variation']
    varPar= np.zeros(len(var))
    for i in range(len(var)):   #getting ride of uncertainty values
        line = var[i]
        if line[-1:]=='V':
            count = line.index('V')
            varPar[i] = line[:count-2]
        else:
            varPar[i] = line
        i+=1
    varList = varPar.astype(np.float)
    return varList

def make_var_hist():
    varList=var_cleanup()
    plt.hist(varList, bins=range(0,56,2),stacked=True, facecolor='blue',ec='black')
    plt.title("137Cs Peak Energy Variation with Source Position")
    plt.xlabel('Calibrated Energy (keV)')
    plt.xticks(range(0,56,5))
    plt.ylabel('Counts')
    plt.yticks(np.arange(0, 25,step=5))
    plt.savefig('PATH')
    plt.show()





'''
make 2D plots
'''

'''2D plot#1: Resolution vs Variation '''
def res_var_2d():
    resList= res_cleanup()
    varList= var_cleanup()
    plt.hist2d(resList,varList, bins=(11,12),range=[[20,42], [0,60]],cmap='binary')
#     print(max(resList))
#     print(varList)
#     print(max(varList))
    plt.colorbar(label='number of xtals',cmap='binary')
    plt.xlabel("137Cs Position 3 Peak Resolution(keV)")
    plt.ylabel("137Cs Total Energy Variation(keV)")
    plt.xticks(range(20,44,2))
    plt.yticks(range(0,65,5))
    plt.title("Energy Variation vs Peak Resolution")
    plt.show()
    plt.savefig('PATH',dpi=100)

'''2D plot #2: Resolution vs voltage '''
def res_volt_2d():
    resList= res_cleanup()
    plt.hist2d(resList,get_voltages(4), bins=(11,10),range=[[20,42], [660,960]],cmap='binary')
    plt.colorbar(label='number of xtals',cmap='binary')
    plt.xlabel("137Cs Position 3 Peak Resolution(keV)")
    plt.ylabel("Voltages (Log(G)=4)")
    plt.xticks(range(20,44,2))
    plt.yticks(range(660,990,30))
    plt.title("Gain vs Peak Resolution")
    plt.show()
    plt.savefig('PATH',dpi=100)


'''2D plot #3: variation vs voltage '''
def var_volt_2d():
    varList= var_cleanup()
    plt.hist2d(varList,get_voltages(4), bins=(12,10),range=[[0,60], [660,960]],cmap='binary')
    #print(max(varList))
    plt.colorbar(label='number of xtals',cmap='binary')
    plt.xlabel("137Cs Total Energy Variation(keV)")
    plt.ylabel("Voltages (Log(G)=4)")
    plt.xticks(range(0,65,5))
    plt.yticks(range(660,990,30))
    plt.title("Gain vs Peak Variation")
    plt.show()
    plt.savefig('PATH',dpi=100)



'''
Identify the xtals with high resolution
'''
def poor_res_list(res_value):
    bad_xtal= []
    resList=res_cleanup()
    for i in range(len(resList)):
        if resList[i] >= res_value:
            #print(i)
            bad_xtal.append(df['Crystal SN'][i])
    return bad_xtal




make_res_hist()
make_var_hist()
make_volt_hist()

res_var_2d()
res_volt_2d()
var_volt_2d()

print('Crystals with poor resolution are', poor_res_list(30))
