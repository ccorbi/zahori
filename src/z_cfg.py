'''
Created on 18/03/2010

@author: ccorbi

Configuration File. 
Define your Wetspots, that are defined by interaction spots
Define the rest of the variables

'''

'''
to do;
 Se tendra que pasar a el main , el parser de las opciones, y todo lo demas se tendra que modificar.. usar un 
 diccionario para las opciones,,,,
 hbfiles
 pdbfiles
'''


import sys, os


def read_cfg_file():
    
    #buscar lo del cwd, puede servir para leer el fichero en el folder de trabjao
    #se podria crear una funcion que pille exactemente y solo aquello que quieras
    
            #SETUP DEAFAULTS
        
    zahori_cfg = {'TRAJ_FILE':'REIMAGE_1BBZ_CHAINA',              # -T
                  'HB_FILE':'frame.pdb.hb2',         # -H
                  'FRAME_PDB':'frame.pdb',           # -N
                  'TIMEXFRAME':2,                  # -X
                  'PDBCHAIN':' ',                  # -J
                  'TOTALFRAMES':5000,              # -F
                  'DB_PATH':'v36b_5000_.ddbb',                 # -D
                  'CENTROID_DISTANCE_CUTOFF':3.5,
                  }
    
    path = os.getcwd()
    cfgfile = path+'/'+'zahori_default.cfg'
    
    try:

        
        rf = open(cfgfile)
        for line in rf.readlines():
            
            if not line[0] == '#':
                linesplited = line.split('\t')
                
                if zahori_cfg.has_key(linesplited[0]):
                    zahori_cfg[linesplited[0]]= linesplited[1].rstrip('\n')
                

    except:
        pass
	#print 'ERROR with configuration file ', cfgfile
        #print 'Use defaults'
    
    
    
    return zahori_cfg


def show_set_argv():
    ''' Print out all setup variables '''
    
    z_options = read_cfg_file()
    
    print '### Files values ###'
    print ' name of HB2 files ',  z_options['HB_FILE']
    print ' name of PDB frames files ', z_options['FRAME_PDB']
    print 'Trajectory file name ', z_options['TRAJ_FILE']
    print ' DDBB file name ', z_options['DB_PATH']
    print '\n'
    print '''### ###'''
    print '''Number of Frames ''', z_options['TOTALFRAMES']
    print '''PicoSeconds per frames ''', z_options['TIMEXFRAME']
    print '''Centroid Distance Cutoff (A) ''', z_options['CENTROID_DISTANCE_CUTOFF']
    
    return

def usage():
    ''' Print the command line help'''
    
    print '''
        ############################    WARNINGS:      ############################
        +Trajectory name without extension
        +Trajectory name must be the same as the parmtop file
        +file extension is case senstive, .prmtop and .mdcrd by default
    
    
        ############################    ACTIONS:       ############################
    -G Generate raw files follow this option by:
        [all] for pdb and hb2 from trajectory file (input mdcrd and number of frames)
        [pdb] for just Extract pdbs from trajectory file  (input mdcrd and number of frames)
        [hb2] for just hb files (input number of frames) optional name of pdb
          
    -C Generate ddbb  
        if no name is setup ... by default z_cfg name
        First check if files exist
        Create ddbb
        Fill ddbb with data
        
    -A Do Analisis of ddbb values
    '''

    sys.exit()
    



if __name__ == "__main__":
    
    opt_parser(sys.argv[1:])  




                    
