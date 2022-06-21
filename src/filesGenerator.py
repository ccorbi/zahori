import os, sys
#import os.path


def call_HBplus(pdbfile, PATH='hbplus', DAdistance='3.5',DAangle='90'):
    ''' Generate a HB2 file from a pdb. hb2 file will have the same name that pdb'''
    
    os.system('hbplus -d %s -a %s  %s -h 2.5' % (DAdistance,DAangle,pdbfile))
    pdbname = pdbfile.split('.')
    os.system('mv %s.%s.hb2 %s.%s.hb2.%s' % (pdbname[0],pdbname[1],pdbname[0],pdbname[1],pdbname[2]))
              
    return

def writeinfile(filename, content):
    
            filename ='%s' % filename
            f = open(filename, 'w')
            f.write( content )
            f.close()
            return

def getpdb(trajfile,stepnum):
    '''Generate a ptraj file and call it to get a odb from frame'''
    ''' stepnum set the final frame'''
    
    ptrajrmsd='trajin %s.mdcrd 0 %s\n  trajout  frame.pdb PDB\n' % (trajfile, stepnum, stepnum)
    writeinfile('getPDB.ptraj',ptrajrmsd)
    os.system('ptraj %s.prmtop getPDB.ptraj' % (trajfile))
    
    return

def gethb2(TOTALFRAMES,pdb_name='frame.pdb'):
    ''' bucle to call hbplus for each frame'''
    
    for i in range(1,TOTALFRAMES):
        file='%s.%s' % (pdb_name,i)
        call_HBplus(file)


def usage():
    print '##################################'
    print '-a: to generate pdbs and hb2, -p just pdb and -h just hb2'
    print 'first trajfile name (prmtop and mdcrd must be the same'
    print 'second the Number of Frames do you want'
    sys.exit()

################MAIN###############


        
if len(sys.argv[1:]) != 3:
    usage()


trajfile = sys.argv[2]
TOTALFRAMES = sys.argv[3]
TOTALFRAMES = int(TOTALFRAMES) +1 

if sys.argv[1] == '-a':
    getpdb(trajfile,TOTALFRAMES)
    gethb2(TOTALFRAMES,pdb_name='frame.pdb')

elif sys.argv[1] == '-p':
    getpdb(trajfile,TOTALFRAMES)
    
else:
    gethb2(TOTALFRAMES,pdb_name='frame.pdb')


