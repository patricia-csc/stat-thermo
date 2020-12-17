###################################################
# main.py: Program for Statistical Thermodynamics #
# Author: Patricia Popa-Mihai                     #
# Date: 12.12.2020                                #
###################################################

import math
import sys

# CONSTANTS
h = 6.62607015 * (10 ** (-34))
k = 1.38064852 * (10 ** (-23))
c = 299792458 
NA = 6.022 * (10 ** (23))
R = 8.31446261815324
P_std = 1
qt_ct = 0.02594678522
H_ionization = 13.6

# function used for printing
def print_table(d: dict, h: [str], t: str):
    print("{:<30}".format(t),end='')
    for i in range(0, len(h)):
        print("{:<20}".format(h[i]), end='')
    print()
    for v in d.items():
        label, num = v
        print("{:<30}".format(label), end='')
        for j in num:
            print("{:{width}.{prec}f}".format(j, width="<20", prec=3), end='')
        print()
    print()

# class containing all the data for a molecule
class Molecule:
    def __init__(self, name: str, coefficient: int, linear: str, sigma: float, 
                 g_el: [float], B: [float], v: [float], d: [float], m: float, 
                 el_lvl: [float]):
        self.name = name # name
        self.m = m # mass
        self.coefficient = coefficient # coefficient in reaction
        self.linear = linear # string, "linear" if linear, "non-linear" otherwise
        self.sigma = sigma # symmetry number sigma
        self.B = B # list of rotational constants, in cm^-1
        self.v = v # list of vibrational wavenumber, in cm^-1
        self.d = d # list of vibrational degeneracies
        self.el_lvl = el_lvl # list of electronic wavenumbers, in cm^-1
        self.g_el = g_el # list of electronic degeracies
        self.partition_functions = self.compute_partition_functions()
        self.internal_energies = self.compute_internal_energies()
        self.entropies = self.compute_entropies()

    # computes parititon functions: qt, qr, qv and qe
    def compute_partition_functions(self) -> dict:
        qt = qt_ct * math.pow(self.m, 3/2) * math.pow(T, 5/2) * math.pow(P / P_std, -1)

        try:
            if self.linear.startswith('linear'):
                qr = (1 / self.sigma) * k * T / (h * c * self.B[0] * 100)
            else: 
                qr = (1 / self.sigma) * math.sqrt(3.14/(self.B[0] * self.B[1] * self.B[2] * 1000000)) * math.pow(k * T / h / c, 3/2)
        except ZeroDivisionError:
            qr = 1

        try:
            qv = 1
            for s in range(0, len(self.v)):
                qv *= math.pow(1 / (1 - math.exp(- h * self.v[s] * c * 100 / k / T)), self.d[s])
        except ZeroDivisionError:
            qv = 1

        qe = 0      
        if len(g_el) > 1:
            for i in range(0, len(g_el)):
                qe += g_el[i] * math.exp(-el_lvl[i] * 100 * h * c / k / T)
        else:
            qe = 1

        partition_functions = {
            "qt": [qt],
            "qr": [qr],
            "qv": [qv],
            "qe": [qe]
        }
        return partition_functions

    # computes internat energies, in kJ/mol: Ut, Ur, Uv and Ue
    def compute_internal_energies(self) -> dict:
        ut = 3/2 * k * T * NA / 1000 

        try:
            if self.B[0] != 0:
                if self.linear.startswith('linear'):
                    ur = k * T * NA / 1000
                else:
                    ur = 3/2 * k * T * NA / 1000
            else:
                ur = 0
        except ZeroDivisionError:
            ur = 0
        
        try:
            uv = 0
            for s in range(0, len(self.v)):
                uv += self.d[s] * ((h*self.v[s]*c*100)/(math.exp((h*self.v[s]*c*100)/(k*T)) - 1))
            uv = uv * NA / 1000
        except ZeroDivisionError:
            uv = 0

        ue = 0      
        if len(g_el) > 1:
            upper_sum = 0
            lower_sum = 0
            for i in range(0, len(g_el)):
                upper_sum += (h*self.el_lvl[i]*c*100)*g_el[i]/k/T/T*math.exp(-h*self.el_lvl[i]*c*100/k/T)
                lower_sum += g_el[i]*math.exp(-h*self.el_lvl[i]*c*100/k/T)
            ue = upper_sum / lower_sum * T * T * k * NA / 1000


        internal_energies = {
            "Ut": [ut],
            "Ur": [ur],
            "Uv": [uv],
            "Ue": [ue],
        }
        return internal_energies
    
    # computes entropies, in J/mol*K: St, Sr, Sv and Se
    def compute_entropies(self) -> dict:
        st = NA * k * (math.log(self.partition_functions["qt"][0]) + 5/2)
        
        if self.B[0] != 0:
            sr = (k * math.log(self.partition_functions["qr"][0]) + self.internal_energies["Ur"][0] / NA * 1000 /T) * NA
        else:
            sr = 0
        
        sv = (k * math.log(self.partition_functions["qv"][0]) + self.internal_energies["Uv"][0] / NA * 1000 /T) * NA 
        
        if len(g_el) > 1:
            se = (k * math.log(self.partition_functions["qe"][0]) + self.internal_energies["Ue"][0] / NA * 1000 /T) * NA
        else:
            se = 0
        
        entropies = {
            "St": [st],
            "Sr": [sr],
            "Sv": [sv],
            "Se": [se]
        }
        return entropies

# class containing all the data for the reaction
class Reaction:
    def __init__(self, dU: float, molecules: [Molecule]):
        self.molecules = molecules # list of molecules in the reaction
        self.dU = dU # zero kelvin reaction energy, kJ/mol
        self.internal_energies = self.internal_energies_method() # interal energies, kJ/mol
        self.entropies = self.entropies_method() # entropies, kJ/mol*K
        self.td_eq = self.td_eq_method() # equilibrium constant K
        self.gibbs = self.gibbs_method() # Gibbs energy, kJ/mol
        self.int_en = self.int_en_method() # internal energy, kJ/mol
        self.enthalpy = self.enthalpy_method() # reaction enthalpy, kJ/mol
        self.entropy = self.entropy_method() # reaction entropy, kJ/mol, from formula G = H - TS
        self.ctrl_entropy = self.ctrl_entropy_method() # reaction entropy, kJ/mol, by adding all entropies

    def internal_energies_method(self) -> dict:
        internal_energies = {"Ut": [0], "Ur": [0], "Uv": [0], "Ue": [0]}
        for key in internal_energies:
            for m in molecules:
                internal_energies[key][0] += m.coefficient * m.internal_energies[key][0]
            internal_energies[key][0] = round(internal_energies[key][0], 2)
        return internal_energies

    def entropies_method(self) -> dict:
        entropies = {"St": [0], "Sr": [0], "Sv": [0], "Se": [0]}
        for key in entropies:
            for m in molecules:
                entropies[key][0] += m.coefficient * m.entropies[key][0]
            entropies[key][0] = round(entropies[key][0], 2)
        return entropies
    
    def td_eq_method(self) -> float:
        td_eq = 1
        for m in molecules:
            tmp = 1
            for key in m.partition_functions:
                tmp *= m.partition_functions[key][0]
            tmp = math.pow(tmp, m.coefficient)
            td_eq *= tmp
        td_eq *= math.exp(-self.dU/(R*T/1000))
        return td_eq

    def gibbs_method(self) -> float:
        gibbs = - R * T * math.log(self.td_eq) / 1000
        return gibbs
    
    def int_en_method(self) -> float:
        int_en = self.dU
        for U in self.internal_energies:
            int_en += self.internal_energies[U][0] 
        return int_en

    def enthalpy_method(self) -> float:
        coef_sum = 0
        for m in molecules:
            coef_sum += m.coefficient
        enth = self.int_en + R * T * coef_sum / 1000
        return enth

    def entropy_method(self) -> float:
        return ((self.enthalpy*1000 - self.gibbs*1000) / T)
    
    def ctrl_entropy_method(self) -> float:
        ctrl_entr = 0
        for e in self.entropies:
            ctrl_entr += self.entropies[e][0]
        return ctrl_entr

    # prints everything
    def printing(self) -> float:
        partition_functions = {"qt": [],"qr": [],"qv": [],"qe": []}
        internal_energies = {"Ut": [],"Ur": [],"Uv": [],"Ue": []}
        entropies = {"St": [],"Sr": [],"Sv": [],"Se": []}
        molecule_names = []
        
        for molecule in molecules:
            molecule_names.append(molecule.name)
            partition_functions = {
                key: partition_functions.get(key) + molecule.partition_functions.get(key)
                for key in sorted(set(partition_functions).union(molecule.partition_functions))
            }
            internal_energies = {
                key: internal_energies.get(key) + molecule.internal_energies.get(key)
                for key in sorted(set(internal_energies).union(molecule.internal_energies))
            }
            entropies = {
                key: entropies.get(key) + molecule.entropies.get(key)
                for key in sorted(set(entropies).union(molecule.entropies))
            }

        internal_energies = {
            key: internal_energies.get(key) + self.internal_energies.get(key)
            for key in sorted(set(internal_energies).union(self.internal_energies))
        }
        entropies = {
            key: entropies.get(key) + self.entropies.get(key)
            for key in sorted(set(entropies).union(self.entropies))
        }
        print()
        print_table(partition_functions, molecule_names, 'Partition Functions')
        molecule_names.append('\u0394U')
        print_table(internal_energies, molecule_names, 'Internal Energies (kJ/mol)')
        molecule_names.remove('\u0394U')
        molecule_names.append('\u0394S')
        print_table(entropies, molecule_names, 'Entropies (J/mol*K)')
        print("Equilibrium constant K: {:.3e}\n\u0394rG (kJ/mol): {:.3f}\n\u0394rU (kJ/mol): {:.3f}"
              "\n\u0394rH (kJ/mol): {:.3f}\nFormula \u0394rS (J/mol*K): {:.3f}"
              "\nControl \u0394rS (J/mol*K): {:.3f}".format(self.td_eq, self.gibbs, 
              self.int_en, self.enthalpy, self.entropy, self.ctrl_entropy))

# main function, reads all data from files
if __name__ == '__main__':
    T = int(input("Temperature, in K: "))
    P = int(input("Pressure, in bar: "))
    dU = float(input("\u0394rU(0), in kJ/mol: "))
    molecules = []

    for input_file in sys.argv[1:]:
        input_file = open(input_file)
        name = input_file.readline().rstrip()
        coefficient = int(input_file.readline().rstrip())
        linear = input_file.readline().rstrip()
        sigma = float(input_file.readline().rstrip())
        m = float(input_file.readline().rstrip())
        B = list(map(float, input_file.readline().rstrip().split(" ")))
        v = list(map(float, input_file.readline().rstrip().split(" ")))
        d = list(map(int, input_file.readline().rstrip().split(" ")))

        try:
            el_lvl = list(map(float, input_file.readline().rstrip().split(" ")))
            g_el = list(map(float, input_file.readline().rstrip().split(" ")))
        except ValueError:
            g_el = []
            el_lvl = []

        molecule = Molecule(name=name, coefficient=coefficient, linear=linear, 
                            sigma=sigma, g_el=g_el, B=B, v=v, m=m, el_lvl=el_lvl, d=d)
        molecules.append(molecule)

    reaction = Reaction(dU=dU, molecules=molecules)
    reaction.printing()