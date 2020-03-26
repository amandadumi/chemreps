'''
Functions fo reading common molecule files and creating the f2b/f3b representation.

Literature References:

Disclaimers:
    - This is an attempt at the recreation from literature and may not be
      implemented exactly as it is in the source listed above.
'''
from .utils.molecule import Molecule
from .utils.calcs import length
from itertools import combinations
import numpy as np
import pandas as pd

def resize(df):
    ## find the largest size representation.
    max_length_2b = df.Rep2b.map(len).max()
    max_length_3b = df.Rep3b.map(len).max()
    df2 = pd.DataFrame(columns = ['Name','Reps'])
    for i in df.index:
        # find the difference of the current rep and the longest rep
        difference_2b = max_length_2b-len(df['Rep2b'][i])
        difference_3b = max_length_3b-len(df['Rep3b'][i])
        # append zeros to the reps shorter than the longest
        if difference_2b !=0:
            correct_len_array_2b = np.append(df['Rep2b'][i],np.zeros(difference_2b))
        if difference_3b != 0:
            correct_len_array_3b = np.append(df['Rep3b'][i],np.zeros(difference_3b))

        final_vector = np.append(correct_len_array_2b,correct_len_array_3b)
        df2[df['Name'][i]] = final_vector
    ## resize all three-body vectors

    return df2


def f2b_f3b(mol_file):
    '''
    Parameters
    ---------
    mol : string

    Returns: ndarray
    -------
    '''
    current_molecule = Molecule(mol_file)
    f2b_part = f2b(current_molecule)
    f3b_part = f3b(current_molecule)
    return np.concatenate((f2b_part, f3b_part),axis=None)


def f2b(mol):
    '''
    Parameters
    ---------
    mol : molecule object or string
        the molecule object is passed when calculating both f2b and f3b at once.
         A molecule file can be passed directly if only two f2b descriptors are
         desired

    Returns: ndarray
        molecule file for reading in coordinates

    -------
    '''
    # This can be passed a molecule or a molecule file. as you may want # to
    # extract two bdoy contribution
    # this just turns the file into a molecule oject.
    if type(mol) == str:
        mol = Molecule(mol)
    M = 15
    feat = np.zeros((max(mol.at_num),max(mol.at_num), M))
    for i in range(mol.n_atom):
        for j in range(mol.n_atom):
            r_ij = length(mol, i, j)
            zi = mol.at_num[i]
            zj = mol.at_num[j]
            Z = np.array(sorted([zi,zj]))-1
            for m in range(1,M+1):
                if r_ij == 0:
                    feat[Z[0],Z[1],m-1] = 0
                else:
                    feat[Z[0],Z[1],m-1] += r_ij**-m
    return feat.ravel()

def f3b(mol):
    '''
    Parameters
    ---------
    mol : molecule object

    Returns:
    ??
    -------
    '''
    ## generate sequences without repetition of 3 elements
    g = combinations(range(1,7),3)
    g_len = len(list(g))
    F = np.zeros((max(mol.at_num),max(mol.at_num),max(mol.at_num),g_len,g_len,g_len))
    for i in range(mol.n_atom):
        for j in range(mol.n_atom):
            for k in range(mol.n_atom):
                r_ij = length(mol, i, j)
                r_ik = length(mol, i, k)
                r_jk = length(mol, j, k)
                Zi = mol.at_num[i]
                Zj = mol.at_num[j]
                Zk = mol.at_num[k]
                has_angle = False
                if r_ij < B(Zi, Zj) and r_ik < B(Zi, Zk):
                    has_angle = True
                elif r_ik < B(Zi, Zk) and r_jk < B(Zj, Zk):
                    has_angle = True
                elif r_ij < B(Zi, Zj) and r_jk < B(Zj, Zk):
                    has_angle = True
                if not has_angle:
                    continue
                Z = np.array(sorted([Zi,Zj, Zk]))-1
                for set in g:
                    F[Z[0],Z[1],Z[2],set[0]-1,set[1]-1,set[3]-1] += np.dot((np.dot(r_ij**-set[0],r_ik**-set[1])),r_ij**-set[3])
    return F.ravel()

def B(i,j):
    return 1.1*return_L(i,j)

def return_L(i,j):
    # Bond lengths as defined by Table six in the Literature
    Ldict={(1,1): 0.74,
    (1,6):1.08,
    (1,8):0.96,
    (1,7):1.01,
    (6,6):1.51,
    (6,8):1.43,
    (6,7):1.47,
    (7,7):1.45,
    (7,8):1.40,
    (8,8):1.48}
    ij_sort = sorted((i,j))
    return Ldict[(ij_sort[0],ij_sort[1])]
