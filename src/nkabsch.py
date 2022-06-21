'''
Created on 09/11/2010

@author: carles
'''

from clases import *
from numpy import * 

def calcR (c1,c2):
    '''Calc covariance Matrix from coord of matched atoms'''
        
    R = zeros((3,3))
    L = len(c1)
        
    for k in range(L):
       for i in range(3):
           for j in range(3):
               R[i][j] += c2[k][i] * c1[k][j]

   # return R the 3x3 PSD Matrix.
    return R
    

def calcRtR( R ):
    """
    Calculate RtR (R-transpose * R).
    (CORRECT)

    @param R: Matrix
    @type  R: 3x3 Matrix

    @precondition: R is a the well formed matrix as per Kabsch.
    @postcondition: RtR is now defined

    @return: M{R^tR}
    @rtype : 3x3 matrix

    """
    
    RtR = dot(transpose(R), R)

    return RtR   

def calcEigenPairs( RtR ):
    """
    Calculate the corresponding eigenpairs for RtR, and sort them accordingly:
    M{m1 >= m2 >= m3}; set M{v3 = v1 x v2} to ensure a RHS
    (CORRECT)

    @postcondition: The eigen pairs are calculated, sorted such that M{m1 >= m2 >= m3} and
    M{v3 = v1 x v2}.

    @param RtR: 3x3 Matrix of M{R^t * R}.
    @type  RtR: 3x3 Matrix
    @return : Eigenpairs for the RtR matrix.
    @rtype  : List of stuff

    """

    eVal, eVec = linalg.eig(RtR)

    # This is cool.  We sort it using Numeric.sort(eVal)
    # then we reverse it using nifty-crazy ass notation [::-1].
    eVal2 = sort(eVal)[::-1]
    eVec2 = [[],[],[]] #Numeric.zeros((3,3), Numeric.Float64)

    # Map the vectors to their appropriate owners        
    if eVal2[0] == eVal[0]:
        eVec2[0] = eVec[0]
        if eVal2[1] == eVal[1]:
            eVec2[1] = eVec[1]
            eVec2[2] = eVec[2]
        else:
            eVec2[1] = eVec[2]
            eVec2[2] = eVec[1]
    elif eVal2[0] == eVal[1]:
        eVec2[0] = eVec[1]
        if eVal2[1] == eVal[0]:
            eVec2[1] = eVec[0]
            eVec2[2] = eVec[2]
        else:
            eVec2[1] = eVec[2]
            eVec2[2] = eVec[0]
    elif eVal2[0] == eVal[2]:
        eVec2[0] = eVec[2]
        if eVal2[1] == eVal[1]:
            eVec2[1] = eVec[1]
            eVec2[2] = eVec[0]
        else:
            eVec2[1] = eVec[0]
            eVec2[2] = eVec[1]

    eVec2[2][0] = eVec2[0][1]*eVec2[1][2] - eVec2[0][2]*eVec2[1][1]
    eVec2[2][1] = eVec2[0][2]*eVec2[1][0] - eVec2[0][0]*eVec2[1][2]
    eVec2[2][2] = eVec2[0][0]*eVec2[1][1] - eVec2[0][1]*eVec2[1][0]

    return [eVal2, eVec2]

def calcBVectors( R, eVectors ):
    """
    Calculate M{R*a_k} and normalize the first two vectors to obtain M{b_1, b_2} and set
    M{b_3 = b_1 x b_2}.  This also takes care of {m2 > m3 = 0}. 
    (CORRECT)

    @postcondition: b-Vectors are defined and returned

    @return: The three B-vectors
    @rtype: List of 3D vectors (Python LOL).
    """

    bVectors = zeros((3,3))

    for i in range(3):
        bVectors[i] = dot(R, eVectors[i])

    bVectors[0] = bVectors[0] / sqrt(add.reduce(bVectors[0]**2))
    bVectors[1] = bVectors[1] / sqrt(add.reduce(bVectors[1]**2))
    bVectors[2] = bVectors[2] / sqrt(add.reduce(bVectors[2]**2))

    bVectors[2][0] = bVectors[0][1]*bVectors[1][2] - bVectors[0][2]*bVectors[1][1]
    bVectors[2][1] = bVectors[0][2]*bVectors[1][0] - bVectors[0][0]*bVectors[1][2]
    bVectors[2][2] = bVectors[0][0]*bVectors[1][1] - bVectors[0][1]*bVectors[1][0]

    return bVectors


def calcU( eVectors, bVectors ):
    """
    Calculate M{U=(u_{ij})=(sum n b_{ki} * a_{kj})} to obtain the best rotation.  Set
    M{sigma_3 = -1 if b3(Ra3) < 0 else sigma_3 = +1}.
    (CORRECT)

    @postcondition: U is defined

    @param eVectors: Eigenvectors for the system.
    @type  eVectors: Eigenvectors

    @param bVectors: BVectors as described by Kabsch.
    @type  bVectors: BVectors

    @return: U the rotation matrix.
    @rtype  :3x3 matrix.

    """

    U = zeros((3,3))

    for k in range(3):
        for i in range(3):
            for j in range(3):
                U[i][j] += dot(bVectors[k][i], eVectors[k][j])

    return U


def get_atoms_coords(Protein, atoms_list):
    
    coord_list = []
     
    for atom in atoms_list:
        coord_list.append(Protein.get_coord(atom.get('atom_num')))
         
    return coord_list
 
def pselected_atoms(atomslist):
    for atom in atomslist:
        print atom.get('atom_num'), atom.get('atom_type'), atom.get('res_num'), atom.get('x_coord')
    
    print 'Total numer of atoms is  ', len(atomslist)
    
    return

def superimpose(Protein_ref, chain_domain_ref, Protein_target, chain_domain_target, res_ref, res_target):
    
    
    #Set the axis origin of both proteins in to the centroid, of the Domain.
    centroid_ref = Protein_ref.centroid()
    centroid_target = Protein_target.centroid()
    Protein_ref.reset_axis_origin_to(centroid_ref)
    Protein_target.reset_axis_origin_to(centroid_target)
    
    #Protein_target.save2pdb('zero_target.pdb')
    #Protein_ref.save2pdb('zero_ref.pdb')
    
    atoms_ref = Protein_ref.get_residues(res_ref, chain_domain_ref, res_type='bb')
    atoms_target = Protein_target.get_residues(res_target, chain_domain_target, res_type='bb') 

    #pselected_atoms(atoms_ref)
    #pselected_atoms(atoms_target)
    
    c1 = get_atoms_coords(Protein_ref, atoms_ref)    
    c2 = get_atoms_coords(Protein_target, atoms_target)
    #print c1
    #
    #print c2
    #print len(c1), len(c2)
    if len(c1) != len(c2):
        print 'Error, the selection of atoms is different'
        raise 
    
    #R = calcR(c1,c2)
    #print R

    #RtR = calcRtR(R)
    #print RtR
    #eval, evec =  calcEigenPairs(RtR)
    #bvec = calcBVectors(R, evec) 
    #rot_matrix = calcU(evec, bvec)
    
    #Protein_target.rotation(rot_matrix)
    #Protein_target.save2pdb('target.pdb')
    #Protein_ref.save2pdb('ref.pdb')
    
    '''OTRO METHODO PARA CONSEGUIR LA ROTATION MATRIX'''
    V, S, Wt = linalg.svd( dot( transpose(c2), c1))
 
    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(linalg.det(V) * linalg.det(Wt))))
 
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]
 
    #RMSD = E0 - (2.0 * sum(S))
    #RMSD = numpy.sqrt(abs(RMSD / L))
 
    #U is simply V*Wt
    U = dot(V, Wt)
    
    Protein_target.rotation(U)
    #Protein_target.save2pdb('target.pdb')
    #Protein_ref.save2pdb('ref.pdb')
    
    print 'Superimposition DONE'
    
    return Protein_target

#atoms_ref = [45, 46, 47, 48, 406, 407, 408, 409, 385, 386, 387, 388, 279, 280, 281, 282]
#atoms_target = [69, 70, 71, 72, 437, 438, 439, 440, 416, 417, 418, 419, 291, 292, 293, 294]
#Protein_ref = Protein('1BBZ.clean.pdb')
#Protein_target = Protein('1FYN.clean.pdb')

#align(Protein_ref, Protein_target, atoms_ref, atoms_target)    