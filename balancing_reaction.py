# Balancing reation between NH3 and UF4. Calculates Gibbs reaction energy and its convex hull
# as function of temperature and gaseus species partial pressure

from chempy import balance_stoichiometry
from pprint import pprint
import numpy as np
import csv
from scipy import spatial

def main():
    reactants=['NH3', 'UF4']
    data_file = open('phases.dat','r')
    phase = data_file.read().split(',')
    nphase=len(phase)
    tmin = 100  # starting temperature
    tstep = 100  # temperature step
    tmax = 3000+tstep  # maximum temperature
    ppmin = 0   # minimum gas partial pressure
    ppmax = 0  # maximum gas partial pressure
    ppstep = 0.1  # partial pressure step
    nppstep = int((ppmax-ppmin)/ppstep)+1  # number of partial pressure steps

# reading reactants temperature and Gibbs energy
    r1 = data_dict('gibbs_' '%s' '.csv' % (reactants[1]))
    r2 = data_dict('gibbs_' '%s' '.csv' % (reactants[0]))
    react = r1
    for t in range(tmin,tmax,tstep):
        react[str(t)][reactants[0]] = r2[str(t)][reactants[0]]

# reading products temperature and Gibbs energy
    prod = data_dict('gibbs_' '%s' '.csv' % (phase[0]))
    for i in range(1,nphase-1):
        p = data_dict('gibbs_' '%s' '.csv' % (phase[i]))
        for t in range(tmin,tmax,tstep):
            prod[str(t)][phase[i]] = p[str(t)][phase[i]]

    all_d_G_t = {}
    product_detail = {}
    no_r = 0
# reactions having 1 product
    elP1 = prod1( reactants, phase, nphase )[0]
    P1 = prod1( reactants, phase, nphase )[1]
    R1 = prod1( reactants, phase, nphase )[2]
    for i in range(0,len(P1)):
        no_r = no_r + 1
        product_detail[no_r] = {}
        R = []
        R.append(R1[i][0])
        R.append(R1[i][1])
        product_detail[no_r]['R'] = R
        P = []
        P.append(P1[i][0])
        product_detail[no_r]['P'] = P
        elP = []
        elP.append(elP1[i][0])
        elP.append(0)
        elP.append(0)
        product_detail[no_r]['elP'] = elP

# reactions having 2 products
    elP2 = prod2( reactants, phase, nphase )[0]
    P2 = prod2( reactants, phase, nphase )[1]
    R2 = prod2( reactants, phase, nphase )[2]
    for i in range(0,len(P2)):
        no_r = no_r + 1
        product_detail[no_r] = {}
        R = []
        R.append(R2[i][0])
        R.append(R2[i][1])
        product_detail[no_r]['R'] = R
        P = []
        P.append(P2[i][0])
        P.append(P2[i][1])
        product_detail[no_r]['P'] = P
        elP = []
        elP.append(elP2[i][0])
        elP.append(elP2[i][1])
        elP.append(0)
        product_detail[no_r]['elP'] = elP

# reactions having 3 products
    elP3 = prod3( reactants, phase, nphase )[0]
    P3 = prod3( reactants, phase, nphase )[1]
    R3 = prod3( reactants, phase, nphase )[2]
    for i in range(0,len(P3)):
        no_r = no_r + 1
        product_detail[no_r] = {}
        R = []
        R.append(R3[i][0])
        R.append(R3[i][1])
        product_detail[no_r]['R'] = R
        P = []
        P.append(P3[i][0])
        P.append(P3[i][1])
        P.append(P3[i][2])
        product_detail[no_r]['P'] = P
        elP = []
        elP.append(elP3[i][0])
        elP.append(elP3[i][1])
        elP.append(elP3[i][2])
        product_detail[no_r]['elP'] = elP

# calculating gibbs energies and convex hull as function of partial pressure and temperature
    for p in range(0,nppstep):
        pp = p*ppstep
        all_d_G_t = T_loop( tmin, tmax, tstep, reactants, react, prod, product_detail, no_r, pp )

# calculating reaction convex hull
        convex_hull(tmin, tmax, tstep, all_d_G_t, no_r, reactants, pp )

# writing reaction Gibbs energy as function of temperature
#        write_d_G(tmin, tmax, tstep, all_d_G_t, no_r)

# obtaining reaction Gibbs energy as function of temperature
def T_loop( tmin, tmax, tstep, reactants, react, prod, product_detail, no_r, pp ):
    a_d_G_t = {}
    for i in range(1,no_r):
        a_d_G_t[i] = {}
        for temp in range(tmin,tmax,tstep):
            a_d_G_t[i][str(temp)] = {}
            d_G_t = delta_G_r( temp, reactants, react, prod, product_detail, pp, i )
            a_d_G_t[i][str(temp)] = d_G_t

    return a_d_G_t

# calculating reaction Gibbs energy
def delta_G_r( temp, reactants, react, prod, product_detail, pp, n_r ):
    d_G_r = 100  # defining initial reaction Gibbs energy

# Gibbs energy of UFx reactant
    ufx = float(react[str(temp)][reactants[1]])*float(product_detail[n_r]['R'][1])
# Gibbs energy of NH3 + partial pressure term
    nh3 = (float(react[str(temp)][reactants[0]]) + 8.617333262145*temp*pp/10**5)*float(product_detail[n_r]['R'][0])
# Gibbs energy of products
    prod_g = 0  # initialize Gibbs energy
    for i in range(0,len(product_detail[n_r]['P'])):
        prod_g = float(product_detail[n_r]['P'][i])*float(prod[str(temp)][str(product_detail[n_r]['elP'][i])]) + prod_g
# ratio between NH3 concentration and total reatants concentration , i.e., NH3/(UFx+NH3)
    ratio = float(product_detail[n_r]['R'][0])/(float(product_detail[n_r]['R'][1])+float(product_detail[n_r]['R'][0]))
# reaction gibbs energy
    d_G_r = (prod_g - (ufx + nh3))*96.4853/(float(product_detail[n_r]['R'][1])+float(product_detail[n_r]['R'][0]))

    n_d_G_r = []
    n_d_G_r.append(ratio)
    n_d_G_r.append(d_G_r)
    n_d_G_r.append(reactants[1])
    n_d_G_r.append(reactants[0])
    n_d_G_r.append(product_detail[n_r]['elP'][0])
    n_d_G_r.append(product_detail[n_r]['elP'][1])
    n_d_G_r.append(product_detail[n_r]['elP'][2])

    return n_d_G_r

# calculating reaction convex hull
def convex_hull(tmin, tmax, tstep, all_d_G_t, no_r, reactants, pp ):
    filename = ('convex_hull_' '%.2f' % (pp))
    ch_r = open('%s.dat' % filename, "w")
    ch_r.write(str(pp) + '\n')

    for temp in range(tmin,tmax,tstep):
        points = [ [0,0], [1,0] ]
        reaction_detail = [ reactants[1], reactants[0] ]
        for i in range(1,no_r):
            # disregarding reactions with positive Gibbs energies (not possible)
            if ( all_d_G_t[i][str(temp)][1] < 0 ):
                add_point = [ all_d_G_t[i][str(temp)][0], all_d_G_t[i][str(temp)][1] ]
                points.append(add_point)
                add_detail = [ str(all_d_G_t[i][str(temp)][4]), str(all_d_G_t[i][str(temp)][5]), str(all_d_G_t[i][str(temp)][6]) ]
                reaction_detail.append(", ".join(add_detail))

        ch_r.write(str(temp) + '\n')
        if (len(points) > 2):
            hull = spatial.ConvexHull(points, qhull_options="Qt")
            for i in range(len(hull.vertices)):
                p = []
                p.append(str(points[hull.vertices[i]][0]))
                p.append(str(points[hull.vertices[i]][1]))
                p.append(str(reaction_detail[hull.vertices[i]]))
                ch_r.write("\t".join(p) + '\n')
        else:
            p1 = []
            p1.append(str(points[0][0]))
            p1.append(str(points[0][1]))
            p1.append(str(reaction_detail[0]))
            p2 = []
            p2.append(str(points[1][0]))
            p2.append(str(points[1][1]))
            p2.append(str(reaction_detail[1]))
            ch_r.write("\t".join(p1) + '\n')
            ch_r.write("\t".join(p2) + '\n')

        ch_r.write(' ' + '\n')

# writing reaction Gibbs energies as function of temperature for all reactions
def write_d_G( tmin, tmax, tstep, all_d_G_t, no_r, pp ):
    filename = ('reaction_Gibbs_' '%.2f' % (pp))
    d_G_f = open('%s.dat' % filename, "w")
    d_G_f.write('#T(K)  d_G(eV/UF4)  ratio[NH3/(UF4+NH3)]' + '\n')
    for i in range(1,no_r):
        reaction = []
        for j in range(2,7):
            reaction.append(str(all_d_G_t[i]['100'][j]))

        d_G_f.write("\t".join(reaction) + '\n')

        for temp in range(tmin,tmax,tstep):
            temp_details = []
            temp_details.append(str(temp))
            temp_details.append(str(all_d_G_t[i][str(temp)][1]))
            temp_details.append(str(all_d_G_t[i][str(temp)][0]))
            d_G_f.write("\t".join(temp_details) + '\n')

        d_G_f.write(' ' + '\n')

    d_G_f.close()


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
    comp = ''.join(str(x) for x in labels)
    dict = {}
    for i in range(0,len(data)):
        dict[data[i][0]] = {}
        for j in range(1,len(data[i])):
            dict[data[i][0]][comp] = float(data[i][j])

    return dict

def prod1( reactants, phase, nphase ):
    p_el = []
    p_n = []
    r_n = []
    for i in range(0,nphase-1):
        A = phase[i]
        try:
            reac, prod = balance_stoichiometry({reactants[1], reactants[0]}, {A})
        except:
            pass
        else:
            R = list(reac.values())
            P = list(prod.values())
            elR = list(reac.keys())
            elP = list(prod.keys())
            if all(i > 0 for i in R) and all(i > 0 for i in P):
                p_el.append(elP)
                p_n.append(P)
                r_n.append(R)

    return p_el, p_n, r_n

def prod2( reactants, phase, nphase ):
    p_el = []
    p_n = []
    r_n = []
    for i in range(0,nphase-2):
        for j in range(i+1,nphase-1):
            A = phase[i]
            B = phase[j]
            try:
                reac, prod = balance_stoichiometry({reactants[1], reactants[0]}, {A, B})
            except:
                pass
            else:
                R = list(reac.values())
                P = list(prod.values())
                elR = list(reac.keys())
                elP = list(prod.keys())
                if all(i > 0 for i in R) and all(i > 0 for i in P):
                    p_el.append(elP)
                    p_n.append(P)
                    r_n.append(R)

    return p_el, p_n, r_n

def prod3( reactants, phase, nphase ):
    p_el = []
    p_n = []
    r_n = []
    for i in range(0,nphase-3):
        for j in range(i+1,nphase-2):
            for k in range(j+1,nphase-1):
                A = phase[i]
                B = phase[j]
                C = phase[k]
                try:
                    reac, prod = balance_stoichiometry({reactants[1], reactants[0]}, {A, B, C})
                except:
                    pass
                else:
                    R = list(reac.values())
                    P = list(prod.values())
                    elR = list(reac.keys())
                    elP = list(prod.keys())
                    try:
                        all(i > 0 for i in R)
                        all(i > 0 for i in P)
                    except:
                        pass
                    else:
                        if all(i > 0 for i in R) and all(i > 0 for i in P):
                            try:
                                float(P[0])
                                float(P[1])
                                float(P[2])
                            except:
                                pass
                            else:
                                p_el.append(elP)
                                p_n.append(P)
                                r_n.append(R)

    return p_el, p_n, r_n

if __name__== "__main__":
  main()