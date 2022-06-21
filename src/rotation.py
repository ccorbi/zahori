#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python2.7-32

from clas import *
import nkabsch

def list_AA_ref(info):

	aa = info.split('-')

	return [a for a in range(int(aa[0]),int(aa[1]))]


protein_ref = Protein('dry.amber.3OBQ.pdb')
protein_ref_domain = ''


info = '1-150'
residues_ref = list_AA_ref(info) 
cero = [.0,.0,.0]
pdb_file = 'frame.pdb.9990'
prot_frame = Protein(pdb_file)

algn_complex = nkabsch.superimpose(protein_ref, protein_ref_domain, prot_frame,'', residues_ref, residues_ref )

algn_complex.save2pdb('superimpose.pdb')
