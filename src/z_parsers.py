
'''
Created on 22/03/2010

@author: ccorbi
'''

import os
import sys
import z_cfg
import re
import traceback

#import os.path
def formatExceptionInfo(maxTBlevel=5):
    
    cla, exc, trbk = sys.exc_info()
    excName = cla.__name__
    try:
        excArgs = exc.__dict__["args"]
    except KeyError:
        excArgs = "<no args>"
    excTb = traceback.format_tb(trbk, maxTBlevel)
    
    return (excName, excArgs, excTb)

def writefile(filename, content):
    
    filename ='%s' % filename
    f = open(filename, 'w')
    f.write( content )
    f.close()
    
    return
        
def add2file(filename, content):
    
    filename ='%s' % filename
    f = open(filename, 'a')
    f.write( content )
    f.close()
    
    return

def look_at_pdb(frame, centroid_coord,wetspot_data):

    water_def = ['WAT','HOH']
    #pdb_file = z_options['FRAME_PDB']+'.'+str(frame)    
    pdb_file = cfg.frame_pdbname+'.'+str(frame)
    structure = read_pdb(pdb_file)
    water_id = None
    min_distance = 10000000000
    
    try:    
        for line in structure:
            if line[18:20] is in water_def:
                        coord[0] = float(line[31:39].strip())
                        coord[1] = float(line[39:47].strip())
                        coord[2] = float(line[47:55].strip())                
                if coord[0] - centroid_coord[0] < wetspot_data.cutoff or coord[1] - centroid_coord[1] < wetspot_data.cutoff or coord[2] - centroid_coord[2] < wetspot_data.cutoff:
                        distance = euclidian_distance(centroid_coord[0],centroid_coord[1],centroid_coord[2],coord[0],coord[1],coord[2])
                        if min_distance > distance:
                            water_id = int(line[21:30].strip())

        if  min_distance < wetspot_data.cutoff:
            return water_id
        else:
            return None


  

def get_atom_coord(res_num, atom_type, frame):
    ''' get the coord of and atom and return in a list
        need the res_num and  atom_type, pdb file name
    '''
    
    z_options = z_cfg.read_cfg_file()
    #get the prefix of the pdb file to the frame
    pdb_file = z_options['FRAME_PDB']+'.'+str(frame)
    #Convert the res_num in a int , just in case
    res_num = int(res_num)
    #Read the pdb to a list of list
    structure = read_pdb(pdb_file)
    #Init the list that will keep the coordinates
    coord = [0.0,0.0,0.0]
    
    try:    
        for line in structure:
            if line[0:4] == 'ATOM':
                #the res_num of the water in big system could be a problem
                if res_num == int(line[21:30].strip()):
                    if atom_type == line[12:16].strip():
                        
                        coord[0] = float(line[31:39].strip())
                        coord[1] = float(line[39:47].strip())
                        coord[2] = float(line[47:55].strip())
                        break
                        
        return coord
                
    except:
        print 'ERROR getting the coordinates for the atom of the residue ', res_num, atom_type
        print 'In the frame  ', frame
        print formatExceptionInfo()
        
def read_pdb(pdb_file):
    ''' Read a pdb file to a list of list'''
    try:
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()
    
        return pdb

    except:
        print 'ERROR reading pdb file  ', pdb_file
        print formatExceptionInfo()

def read_HBplus(hbplusfile):
    '''Parse HBplus file to dicctionary Varialbe and return it'''
       
    try:
            
        hbplusfileinfo = hbplusfile.split('.')
        
        hb_dictionary = {}
        rf = open(hbplusfile)
        for line in rf.readlines():
            #The lines with information start with the name of the chain and then the number of residue 
            if re.match(".\d\d\d",line): 
                hb_num = int(line[71:75])
                hb_dictionary[hb_num] = { \
                    'donor_res': line[6:9],
                    'donor_res_id': int(line[1:5]),
                    'donor_res_atomtype':line[10:13].rstrip(' '),
                    'accep_res': line[20:23],
                    'accep_res_id': int(line[15:19]),
                    'accep_res_atomtype': line[24:27].rstrip(' '),
                    'DAdistance': float(line[28:32]),
                    'DAangle':float(line[46:51]),
                    'frame': int(hbplusfileinfo[3]),
                    'donor_chain': line[0:1], 
                    'accep_chain': line[14:15]}
     
    except:
         
        print "ERROR PARSING HBPLUS FILE"
        print formatExceptionInfo()
                    
    
    return hb_dictionary

def call_HBplus(pdbfile, PATH='hbplus', DAdistance='3.5',DAangle='90'):
    ''' Generate a HB2 file from a pdb. hb2 file will have the same name that pdb'''
    
    try:
        os.system('hbplus -d %s -a %s  %s -h 2.5 >> hb_gen.log' % (DAdistance,DAangle,pdbfile))
        pdbname = pdbfile.split('.')
        os.system('mv %s.%s.hb2 %s.%s.hb2.%s' % (pdbname[0],pdbname[1],pdbname[0],pdbname[1],pdbname[2]))
    
    except:
        print "ERROR CALLING HBPLUS"
        print formatExceptionInfo()
                  
    return


def getpdb(trajfile='trajfile',stepnum=0):
    
    '''Generate a ptraj file and call it to get a odb from frame'''
    ''' stepnum set the final frame from Totalframes variable'''
    
    z_options = z_cfg.read_cfg_file()
    
    trajfile = z_options['TRAJ_FILE']
    stepnum = z_options['TOTALFRAMES']
    
    print 'Generating pdb.........'

    try:
        ptrajrmsd='trajin %s.mdcrd 0 %s\n  trajout  frame.pdb PDB\n' % (trajfile, stepnum)
        #print ptrajrmsd
        writefile('getPDB.ptraj',ptrajrmsd)
        os.system('ptraj %s.prmtop getPDB.ptraj' % (trajfile))
    except:
        
        print 'ERROR generating pdb....'
        print formatExceptionInfo()
        
    return

def gethb2(totalframes=1000,pdb_name='frame'):
    ''' bucle to call hbplus for each frame'''
    
    z_options = z_cfg.read_cfg_file()
    pdb_name = z_options['FRAME_PDB']
    totalframes = z_options['TOTALFRAMES']
    totalframes = int(totalframes) + 1
    #Sum 1 to TOTALFRAMES about range method does not include the final number
    print 'Generating HB2 files......'

    try:
        
        for i in range(1,totalframes):
            file='%s.%s' % (pdb_name,i)
            call_HBplus(file)
    
    except:
        print 'ERROR generating HB2 files......'
        print formatExceptionInfo()


def main():
    
    pass
         
if __name__ == "__main__":

    main()
