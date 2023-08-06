import numpy as np
import re
from scipy.optimize import curve_fit

# Define the path to the input files directory
path = "Input_files/"

# Specify the file names for the input data files
RFU_EXP = path + "Peptide.txt"       # The peptide file
RFU_all_model = path + "Modelled_RFU.txt"  # RFU data modelled by the modeller.
Kint = path + "K_int.txt"           # The Kint file

bound = 2   # The bound applied on the constrained Lnp values.

def read_run_count():
    """
    Read and count the number of runs in the model RFU data.
    """
    with open(RFU_all_model, 'r') as rf2:
        fil = rf2.readlines()
        global RUN
        RUN = [ii.startswith("Rmse") for ii in fil].count(True)
        print("RUN:", RUN)

def read_kint_data():
    """
    Read and store the information in the Kint file.
    """
    K_int = []
    with open(Kint, 'r') as kin:
        for i in kin:
            if i.strip():
                kkk = i.split()
                K_int.append(float(kkk[2]))
    return K_int

def read_input_data():
    """
    Read and store the input data submitted to the modeller.
    """
    with open(RFU_EXP, 'r') as read:
        for i, line in enumerate(read):
            if i == 0 and line.strip():
                global T
                T = [float(s) for s in re.findall(r'-?\d+\.?\d*', line)]

def process_data():
    """
    Process the data from RFU_exp and store them in appropriate lists.
    """
    Deterum = [[] for dtt in range(RUN)]
    RFU_exp = [[] for dtt in range(RUN)]

    res_range = []
    S_E_residue = []
    peptide_No = []

    bs = -1
    with open(RFU_all_model, 'r') as rf2:
        for i, line in enumerate(rf2):
            sp = line.split()
            if len(sp) == 1:
                bs = bs + 1
            if len(sp) != 1:
                RFU_exp[bs].append([float(jli) for jli in sp])

    with open(RFU_EXP, 'r') as read:
        for i, line in enumerate(read):
            if i != 0 and line.strip():
                sp = line.split()
                exch = int(sp[-1]) - int(sp[-2])
                S_E_residue.append(sp[-2::])
                count = str(sp[0][1:]).count("P")
                if count != 0:
                    exch = exch - count
                peptide_No.append(int(sp[1]))
                res_range.append([int(Aoo) for Aoo in sp[-2:]])
                for jonas in range(RUN):
                    RFU = RFU_exp[jonas][i - 1]
                    data = [float(b) * exch for b in RFU]
                    Deterum[jonas].append([g for g in data])

    Amino_acid = []
    Pep1 = []
    Pep2 = []
    for ii in range(len(res_range)):
        for jj in range(len(res_range) - ii):
            AB = np.absolute(np.subtract(res_range[ii], res_range[jj + ii]))
            if AB[0] <= 1 and AB[1] <= 1 and AB[0] != AB[1]:
                Pep1.append(ii)
                Pep2.append(jj + ii)
                if res_range[ii][0] != res_range[jj + ii][0]:
                    Amino_acid.append(max(res_range[ii][0], res_range[jj + ii][0]))
                else:
                    Amino_acid.append(max(res_range[ii][1], res_range[jj + ii][1]))

    Kobs = [[] for lsiss in range(RUN)]

    pep_sele = []
    for ijl in range(len(Pep1)):
        pep_sele.append(str([peptide_No[Pep1[ijl]], peptide_No[Pep2[ijl]]]))

    for onls in range(RUN):
        for lson in range(len(Amino_acid)):
            x = np.array(T)
            y = np.absolute(np.subtract(Deterum[onls][Pep1[lson]], Deterum[onls][Pep2[lson]]))

            def func(x, a):
                return 1 - np.exp(-a * x)

            popt, pcov = curve_fit(func, x, y)
            Kobs[onls].append(popt[0])
            print(pep_sele[lson], Amino_acid[lson], "---", popt[0])  # screen

    return Amino_acid, pep_sele, Kobs

def generate_const_file(Amino_acid, pep_sele, Kobs):
    """
    Generate a const file ("Const.txt"), with a tolerance level of ±2.
    """
    res_rep = []
    with open("Const.txt", 'w') as out:
        out.write("Residue	Kobs2	Kobs1\n")
        for i, (Ami, pepp) in enumerate(zip(Amino_acid, pep_sele)):
            kobs_all = []
            if Ami not in res_rep:
                if i != 0:
                    out.write("\n")
                res_rep.append(Ami)
                out.write("%d  " % (Ami))
                for ll, kob in enumerate(Kobs):
                    kobs_all.append(kob[i])
                fi = np.log(np.divide(K_int[Ami], kob[i]))
                fi1 = np.divide(K_int[Ami], np.exp(fi + bound))
                fi__1 = np.divide(K_int[Ami], np.exp(fi - bound))
                out.write("  %.15f %.15f" % (fi1, fi__1))

    print(res_rep)

if __name__ == "__main__":
    read_run_count()
    K_int = read_kint_data()
    read_input_data()
    Amino_acid, pep_sele, Kobs = process_data()
    generate_const_file(Amino_acid, pep_sele, Kobs)
