

        self.params['close_contacts_dist1_cutoff'] = 2.5
        self.params['close_contacts_dist2_cutoff'] = 4.0
        self.params['electrostatic_dist_cutoff'] = 4.0
        self.params['active_site_flexibility_dist_cutoff'] = 4.0
        self.params['hydrophobic_dist_cutoff'] = 4.0
        self.params['hydrogen_bond_dist_cutoff'] = 4.0
        self.params['hydrogen_bond_angle_cutoff'] = 40.0
        self.params['pi_padding_dist'] = 0.75
        self.params['pi_pi_interacting_dist_cutoff'] = 7.5
        self.params['pi_stacking_angle_tolerance'] = 30.0
        self.params['T_stacking_angle_tolerance'] = 30.0
        self.params['T_stacking_closest_dist_cutoff'] = 5.0
        self.params['cation_pi_dist_cutoff'] = 6.0
        self.params['salt_bridge_dist_cutoff'] = 5.5
        self.params['receptor'] = ''
        self.params['ligand'] = ''
        self.params['output_dir'] = ''
        self.params['output_file'] = ''

    def angle_between_three_points(self, point1, point2, point3): # As in three connected atoms
        vector1 = self.vector_subtraction(point1, point2)
        vector2 = self.vector_subtraction(point3, point2)
        return self.angle_between_points(vector1, vector2)



for hydrogen in hydrogens:
                            if math.fabs(180 - functions.angle_between_three_points(ligand_atom.coordinates, hydrogen.coordinates, receptor_atom.coordinates) * 180.0 / math.pi) <= parameters.params['hydrogen_bond_angle_cutoff']:



if (ligand_atom.element == "O" or ligand_atom.element == "N") and (receptor_atom.element == "O" or receptor_atom.element == "N"):
                        
                        # now build a list of all the hydrogens close to these atoms
                        hydrogens = []
                        
                        for atm_index in ligand.AllAtoms:
                            if ligand.AllAtoms[atm_index].element == "H": # so it's a hydrogen
                                if ligand.AllAtoms[atm_index].coordinates.dist_to(ligand_atom.coordinates) < 1.3: # O-H distance is 0.96 A, N-H is 1.01 A. See http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
                                    ligand.AllAtoms[atm_index].comment = "LIGAND"
                                    hydrogens.append(ligand.AllAtoms[atm_index])
                            
                        for atm_index in receptor.AllAtoms:
                            if receptor.AllAtoms[atm_index].element == "H": # so it's a hydrogen
                                if receptor.AllAtoms[atm_index].coordinates.dist_to(receptor_atom.coordinates) < 1.3: # O-H distance is 0.96 A, N-H is 1.01 A. See http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
                                    receptor.AllAtoms[atm_index].comment = "RECEPTOR"
                                    hydrogens.append(receptor.AllAtoms[atm_index])
                        
                        # now we need to check the angles
                        for hydrogen in hydrogens:
                            if math.fabs(180 - functions.angle_between_three_points(ligand_atom.coordinates, hydrogen.coordinates, receptor_atom.coordinates) * 180.0 / math.pi) <= parameters.params['hydrogen_bond_angle_cutoff']:
                                hbonds_key = "HDONOR_" + hydrogen.comment + "_" + receptor_atom.SideChainOrBackBone() + "_" + receptor_atom.structure
                                pdb_hbonds.AddNewAtom(ligand_atom.copy_of())
                                pdb_hbonds.AddNewAtom(hydrogen.copy_of())
                                pdb_hbonds.AddNewAtom(receptor_atom.copy_of())
                                self.hashtable_entry_add_one(hbonds, hbonds_key)