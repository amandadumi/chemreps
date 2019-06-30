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


def f2b_f3b(mol_file):
    '''
    Parameters
    ---------
    mol : string

    Returns: ndarray
    -------
    '''
    current_molecule = Molecule(mol_file)
    f2b_part = f2b(mol)
    f3b_part = f3b(mol)
    return f2b_part, f3b_part


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
    feat = np.zeros(mol.natoms, 15)
    for i in range(mol.n_atom):
        for j in range(mol.n_atom):
            r_ij  = length(mol, i, j)
            zi = mol.at_num[i]
            zj = mol.at_num[j]
            Z = sorted([zi,zj])
            for m in range(1,M+1):
                feat[Z,m] += calc_distance**-m

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
    g = itertools.combinations(g,3)
    atom = itertools.combination(np.arange(mol.n_atom),3)
    for i in range(mol.n_atom):
        for j in range(mol.n_atom):
            for k in range(mol.n_atom):
                r_ij = length(mol, i, j)
                r_ik = length(mol, i, k)
                r_jk = length(mol, j , k)
                Zi = mol.at_num(i)
                Zj = mol.at_num(j)
                Zk = mol.at_num(k)
                has_angle = False
                if r_ij < B(Zi, Zj) and rik < B(Zi, Zk):
                    has_angle = True
                elif r_ik < B(Zi, Zk) and rjk < B(Zj, Zk):
                    has_angle = True
                elif r_ij < B(Zi, Zj) and rjk < B(Zj, Zk):
                    has_angle = True
                if not has_angle:
                    continue
                Z = sorted([Zi,Zj, Zk])
                for set in g:
                    F[Z,set[0],set[1],set[3]] += np.dot((np.dot(r_ij**-set[0],r_ik**-set[1])),r_jj**-set[3])


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
    (8,8):1.48,
    (8,7):1.40,
    (7,7):1.45}
    return Ldict[(i,j)]
