# Defect formation energy as function of temperature and other defects
# User must provide parameters for the calculations

import math
import csv
import scipy
from scipy import optimize
from scipy.stats import linregress
from scipy.interpolate import interp1d
from operator import itemgetter
import argparse

# Options
parser = argparse.ArgumentParser()
parser.add_argument('--s', type=str, required=True, help='Entropy is ON or OFF')
parser.add_argument('--tmin', type=float, required=True, help='starting temperature')
parser.add_argument('--tmax', type=float, required=True, help='maximum temperature')
parser.add_argument('--tinc', type=float, required=True, help='temperature increment')
parser.add_argument('--oh',   type=float, required=True, help='oxygen enthalpy')
parser.add_argument('--opp',  type=float, required=True, help='oxygen partial pressure')
args = parser.parse_args()

Ent_flag = args.s # entropy on or off
# universal constants
kB = 8.617E-5

#set the conditions for T and oxygen partial pressure
def main():
    filename = "energy_data.csv"
    data = data_dict(filename)
    tmin = args.tmin
    tmax = args.tmax
    tstep = args.tinc
    pp_H = args.oh
    pp_O = args.opp
    T_loop(tmin,tmax,tstep,pp_H,pp_O,data,'conc')


def Opp_vs_T(Opp_H,Opp_0,T):
    return Opp_0 * math.exp(Opp_H/kB/T)

def T_loop(T_min,T_max,T_step,Opp_H,Opp_0,data,output):
    T = T_min
    header = "T,1/T,Opp,stoich"
    for i in neutral(data,T,Opp_vs_T(Opp_H,Opp_0,T)):
        if "perf" not in i:
            header += ',' + i

    filename = ('output')
    ch_r = open('%s.csv' % filename, "w")
    ch_r.write(header + '\n')
    while T <= T_max:
        results = neutral(data,T,Opp_vs_T(Opp_H,Opp_0,T))
        OM = stoich(results)
        string = str(T)+","+str(1/T)+","+str(Opp_vs_T(Opp_H,Opp_0,T))
        string += ',' + str(OM)
        for i in results:
            if "perf" not in i and i != "ref_X_phase":
                string += ',' + str(neutral(data,T,Opp_vs_T(Opp_H,Opp_0,T))[i][output])
            if i == "ref_X_phase":
                if data[i] == '':
                    string += ',NULL'
                else:
                    string += ','+ data[i]

        ch_r.write(string + '\n')
        T+=T_step
    return

def stoich(results):
    O_excess = 2.0
    U_excess = 1.0
    for i in results:
        if "perf" not in i and i != "ref_X_phase":
            conc = results[i]['conc']
            if i == 'vO':
                O_excess -= conc
            elif i == 'iO':
                O_excess += conc
            elif i == 'vU':
                U_excess -= conc
            else:
                pass
    return (O_excess/U_excess)

def neutral(data,T,Opp):
    Ef = charge_neutral(data,T,Opp)[0]
    result = defect_concs(FormG(data,T,Opp,Ef),T)
    return result

def charge_neutral(data,T,Opp):
    def charge(Ef):
        return total_charge(defect_concs(FormG(data,T,Opp,Ef),T))
    sol = scipy.optimize.root(charge, 1.0, method='hybr')
    if total_charge(defect_concs(FormG(data,T,Opp,sol.x[0]),T)) > 1e-10:
        print("Charge not converged")
        return 0
    else:
        return sol.x

def total_charge(concs):
    Q = float(0.0)
    for i in concs:
        if "perf" not in i and i != "ref_X_phase":
            q = concs[i]['q'] * concs[i]['conc']
            Q += q
    return Q

def defect_concs(formations,T):
    for i in formations:
        if "perf" not in i and i != "ref_X_phase":
            conc = formations[i]['mult'] * math.exp(-formations[i]['G']/kB/T)
            formations[i]['conc'] = conc
    return formations

# calculates formation energy for each defect using the defect data read in by data_dict() function
def FormG(data,T,Opp,Ef):
    for i in data:
        if "perf" not in i and i != "ref_X_phase":
            G = data[i]['E'] - data[i]['nM'] * mu_M(T,Opp,data) - data[i]['nO'] * mu_O(T,Opp) + data[i]['q'] * mu_e(Ef,data) - data[i]['nX'] * mu_X(T,Opp,data)["ref_X"]
            if 'S' in data[i].keys():
                entropy = T * data[i]['S'] *kB
                G -= entropy
            data[i]['G'] = G
    data['ref_X_phase'] = mu_X(T,Opp,data)["ref_phase"]
    return data


# electron potential based on calcuation of positions in the band gap from e and h energies from DFT
def mu_e(E_F,data):
    mid_band_gap = (data['e']['E'] - data['h']['E']) / 2.0
    Eg = data['e']['E'] + data['h']['E']
    VBM = mid_band_gap - Eg / 2.0
    return VBM + E_F

# determine the oxygen potential for a set of conditions
def mu_O(T,Opp):
    muO20 = -10.55 # O2 molecule at reference conditions
    Opp0 = 1.0 # reference partial pressure in atm
    SO2 = 0.0021 # O2 entropy in eV/K
    T0 = 297.0 # reference temperature
    dT = -0.5 * (SO2) * (T - T0)
    dOpp = 0.5 * kB * T * math.log(Opp/Opp0) # change due to partial pressure
    muO = muO20 / 2.0 + dT + dOpp # O atom potential
    return muO

# determine the cation potential for a set of conditions - this is speci
def mu_M(T,Opp,data):
    mu_M = data['MO2_perf']['E'] / data['MO2_perf']['nM'] - data['MO2_perf']['nO'] / data['MO2_perf']['nM'] * mu_O(T,Opp)
    if Ent_flag == 'ON':
        mu_M -= data['MO2_perf']['S'] / data['MO2_perf']['nM'] * T * kB
    return  mu_M

# determine dopant chemical potential
def mu_X(T,Opp,data):
    mu_X_lowest = 1e10
    i_lowest = str()
    for i in data:
        if "X" in i and "perf" in i:
            mu_X = data[i]['E'] / data[i]['nX'] - data[i]['nO'] / data[i]['nX'] * mu_O(T,Opp)
            if Ent_flag == 'ON':
                mu_X -= data[i]['S'] / data[i]['nX'] * T * kB
            if mu_X < mu_X_lowest:
                mu_X_lowest = mu_X
                i_lowest = i
    return {"ref_phase":i_lowest,"ref_X":mu_X_lowest}

# read in information about defects defect needs: species, E, nM, nO, Q, mult
def data_dict(filename):
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        data = []
        for i in csv_reader:
            datum = i
            data.append(datum)
    labels = data[0]
    data.pop(0)
    labels.pop(0)
    dict = {}
    for i in range(0,len(data)):
        dict[data[i][0]] = {}
        for j in range(1,len(data[i])):
            dict[data[i][0]][labels[j-1]] = float(data[i][j])
    for i in dict:
        if "perf" not in i:
            for j in dict[i]:
                if j == 'E':
                    dict[i][j] -= dict['MO2_perf'][j]
                if j == 'S':
                    if Ent_flag == 'OFF':
                        dict[i][j] = 0.0
                    elif Ent_flag == 'ON':
                        dict[i][j] -= dict['MO2_perf'][j]
                    else:
                        print("Ent_flag not recognized")
    return dict


if __name__== "__main__":
  main()
