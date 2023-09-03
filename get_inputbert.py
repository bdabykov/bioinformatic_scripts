import pandas as pd
from pathlib import Path
import itertools

def split(str, num):
    return [str[start:start+num] for start in range(0, len(str), num)]

with open("aa_example.txt", "r") as f1:
    list_aa = []
    for aa in f1:
        prot_seq = list(aa)
        lp = len(prot_seq)
        position_aa = prot_seq[30:31]
        list_aa.append(position_aa)

with open("nn_example.txt", "r") as f2:  
    codon_list = []
    check=[]  # for checking the codons that do not exist
    for nuc_seq in f2:
        x=3 
        spl=[nuc_seq[y-x:y] for y in range(x, len(nuc_seq)+x,x)]  # split the nn sequences into codons
        new_list = [item.strip() for item in spl]  #removing \n from list
        noempty= [x for x in new_list if x]  # removing empty items from list
        check.append(noempty)
        
        codon_list.append(spl[30:31])  # appending the 31 codon to the list

data_f = pd.DataFrame(
{'position':
{0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7, 7: 8, 8: 9, 9: 10, 10: 11,
11: 12, 12: 13, 13: 14, 14: 15, 15: 16, 16: 17, 17: 18, 18: 19, 19: 20, 20: 21, 21: 22, 22: 23, 
23 : 24, 24: 25, 25: 26, 26: 27, 27: 28, 28: 29, 29: 30, 30: 31, 31: 32, 32: 33, 33: 34, 34: 35,
35: 36, 36: 37, 37: 38, 38: 39, 39: 40, 40: 41, 41: 42, 42: 43, 43: 44, 44: 45, 45: 46, 46: 47, 47: 48,
48: 49, 49: 50, 50: 51, 51: 52, 52: 53, 53: 54, 54: 55, 55: 56, 56: 57, 57: 58, 58: 59, 59: 60, 60: 61}, 
'codon': 
{0:"GCT",1:"GCC",2:"GCA",3:"GCG",4:"CGT",5:"CGC",6:"CGA",7:"CGG",8:"AGA",9:"AGG",10:"AAT",
11:"AAC",12:"GAT",13:"GAC",14:"TGT",15:"TGC",16:"CAA",17:"CAG",18:"GAA",19:"GAG",20:"GGT",
21:"GGC",22:"GGA",23:"GGG",24:"CAT",25:"CAC",26:"ATT",27:"ATC",28:"ATA",29:"TTA",30:"TTG",
31:"CTT",32:"CTC",33:"CTA",34:"CTG",35:"AAA",36:"AAG",37: "ATG",38:"TTT",39:"TTC",40:"CCT",
41:"CCC",42:"CCA",43:"CCG", 44:"TCT",45:"TCC",46:"TCA",47:"TCG",48:"AGT",49:"AGC",50:"ACT",
51:"ACC",52:"ACA",53:"ACG",54:"TGG",55:"TAT",56:"TAC",57:"GTT",58:"GTC",59:"GTA",60:"GTG"}, 
'aminoacid': 
{0:'A', 1:'A', 2:'A', 3:'A', 4:'R', 5:"R", 6:"R",7:"R",8:"R",9:"R",10:"N",11:"N", 12:"D",13:"D",
14:"C", 15:"C",16:"Q",17:"Q",18:"E",19:"E",20:"G",21:"G",22:"G",23:"G",24:"H",25:"H",26:"I",
27:"I",28:"I",29:"L",30:"L",31:"L",32:"L",33:"L",34:"L",35:"K",36:"K",37:"M",38:"F",39:"F",40:"P",
41:"P",42:"P",43:"P",44:"S",45:"S",46:"S",47:"S",48:"S",49:"S",50:"T",51:"T",52:"T",53:"T",54:"W",
55:"Y",56:"Y",57:"V",58:"V",59:"V",60:"V"}})


new_aalist = list(itertools.chain(*list_aa))  # removing a list from the codon_list
new_codonlist = list(itertools.chain(*codon_list)) # removing a list from the codon_list
vector_list = []    # saving the vectors into a list

codon_list = []
for c in new_codonlist:
    where = (data_f['codon'] == c)
    if where.any():
        pos = data_f.at[where.idxmax(), 'position']
        l_codon = f"codon {c}: position {pos}"
        codon_list.append(l_codon)
    else:
        print(f"codon {c} not found")

l_aa = []
for aa in new_aalist:
    where = (data_f['aminoacid'] == aa)
    if where.any():
            pos = data_f.at[where.idxmax(), 'position']
            
            loc_aa = f"aminoacid {aa}"
            vector = '0'*(pos-1) + '1' + '0'*(61-pos)  # getting the vector, str
            l_aa.append(vector)
    else:
            print(f"codon {c} not found")

#saving the vector in txt file

with open('codon_pos_onehot2.txt', 'w') as f:
    for x, y in zip(codon_list, l_aa):
        f.write(f"{x,y}\n")
print("End")