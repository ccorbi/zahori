#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python2.7-32

'''
#!/usr/bin/python
Created on 18/03/2010

@author: ccorbi
'''

import sys, os
import argparse
import time
import math
import bisect
import re
import z_ddbb
#from clases import *
#import nkabsch
from ConfigParser import SafeConfigParser






'''
##################################
############# CLASSES ############
##################################
'''

class Wetspot(object):

        def __init__(self, id = ''):
                self.id = id
                self.interactions = []
                self.contacts = []
                self.cutoff = 2.6
                self.probs = [2,10,30,50,80,100,200,300,500,1000]
                self.coord_x = 0.0
                self.coord_y = 0.0
                self.coord_z = 0.0

        def add_interaction(self,inter_id,res_id,res_name,atom_type):
                interaction = {'interaction_id':inter_id,
                                'res_id':res_id,
                                'res_name':res_name,
                                'atom_type':atom_type}

                self.interactions.append(interaction)
                self.contacts.append(interaction)
                return

        def get_interaction(self,interaction_id):
                """Returns a interaction of the wetspot which the same interaction id"""
                for interaction in slef.interactions:
                        if interaction_id == interaction[interaction_id]:
                                return interaction
                return None

        def howmany_interactions(self):

                return len(self.interactions)

        def add_contact(nter_id,res_id,res_name,atom_type):
                contact = {'interaction_id':inter_id,
                                'res_id':res_id,
                                'res_name':res_name,
                                'atom_type':atom_type}
                self.contacts.append(contact)
                return


class Edfreq(object):

        def __init__(self, histogram, name=''):
                """Makes a edf from an unsorted sequence of (value, frequency) pairs Dictionary type."""

                self.name = name
                self.items = histogram.iteritems()
                runsum = 0
                xs = []
                cs = []

                for value, count in sorted(self.items):
                        runsum += count
                        xs.append(value)
                        cs.append(runsum)

                total = float(runsum)
                ps = [c/total for c in cs]

    
                self.xs = [] if xs is None else xs
                self.ps = [] if ps is None else ps
        

        def prob(self, x):
                """Returns edf(x), the probability that corresponds to value x.

                Args:
            x: number

                Returns:
            float probability
                """
                if x < self.xs[0]: return 0.0
                index = bisect.bisect(self.xs, x)
                p = self.ps[index-1]
                return p

        def value(self, p):
                """Returns Inverseedf(p), the value that corresponds to probability p.

        Args:
            p: number in the range [0, 1]

        Returns:
            number value
                """
                if p < 0 or p > 1:
                        raise ValueError('Probability p must be in range [0, 1]')

                if p == 0: return self.xs[0]
                if p == 1: return self.xs[-1]
                index = bisect.bisect(self.ps, p)
                if p == self.ps[index-1]:
                        return self.xs[index-1]
                else:
                        return self.xs[index]

        def render(self,filename):
            """Generates a sequence of points suitable for plotting.

        An empirical edf is a step function; linear interpolation
        can be misleading."""

            content = '0'+'\t'+'0'+'\n'
            add2file(filename,content)
            
            top = self.value(1)
            
            for a in range(0,top,2):
                content = str(a)+'\t'+str(self.prob(a))+'\n'
                add2file(filename,content)


            content = str(top)+'\t'+str(self.prob(top))+'\n'
            add2file(filename,content)

            
            return

'''
##################################
####### SETUP & CONFIGURATION ########
##################################
'''


class Cfg(object):

        def __init__(self,filename=''):
                self.filename = filename
                self.trj_filename = 'REIMAGE_1BBZ_CHAINA'
                self.hb_filename = 'frame.pdb.hb2'
                self.frame_pdbname = 'frame.pdb'
                self.ddbb = '3obx.ddbb'

                self.path_ptraj ='ptraj'
                self.path_hbplus = 'hbplus'
                self.path_naccess = 'naccess'

                ###########
                self.wetspots_file = 'wetspots.inp'
                self.timeXframe = 2
                self.first_frame = 1
                self.total_frames = 12000
                self.cutoff_1 = 3.5
                self.cutoff_2 = 4
                self.intesive = True #Problematic, Need improvme
                self.working_path = os.getcwd()
                self.home_path = os.path.dirname(sys.argv[0])
                self.wetspots_path = self.working_path+'/'+self.wetspots_file
                if self.filename != '':
                        self.load()
                ############

        def load(self):
                parser = SafeConfigParser()
                if not parser.read(self.filename):
                        print "WARNING!"
                        print "\t Configuration file does not exist, loading default values"

                else:
                        parser.read(self.filename)

                        return


def parse_wetspot_cfg_error(name,value):

        print "\t ERROR with the interaction definition  ", name, value
        print "\t One or More entries have a not valid format. remember to avoid the use of \":\" "
        print "\t \tThe valid format is:"
        print "\t \t[WETSPOT_ID]"
        print "\t \tinteraction_name = numeric_res_id - residue_name - Atom_type"
        return


def parse_wetspot_cfg(filename):

        parser = SafeConfigParser()
        if not parser.read(filename):
                print 'CRITICAL ERROR!!'
                print 'Wetspots Definition file do not found  ', filename
                sys.exit()
        else:
                parser.read(filename)

        wetspots = {}

        for wetspot_def in parser.sections():
                wetspot = Wetspot(wetspot_def)
                #print '  Options:', parser.options(wetspot_def)
                try:
                        for name, value in parser.items(wetspot_def):

                                if name == 'cutoff':
                                        wetpot.cutoff = float(value)
                                        continue

                                if name == 'coord':
                                        atom_def = value.split(',')
                                        wetspot.coord_x = float(atom_def[0])
                                        wetspot.coord_y = float(atom_def[1])
                                        wetspot.coord_z = float(atom_def[2])
                                        continue

                                else:
                                        atom_def = value.split('-')

                                if len(atom_def) != 3:
                                        print 'WARNING!  ', atom_def
                                        parse_wetspot_cfg_error(name, value)
                                        continue
                                else:
                                        interaction_id = name
                                        res_id = int(atom_def[0].strip())
                                        res_name = atom_def[1].strip()
                                        atom_type = atom_def[2].strip()
                                        if name[0] == 'C':
                                                wetpots.add_contact(interaction_id, res_id, res_name, atom_type)
                                        else:
                                                wetspot.add_interaction(interaction_id, res_id, res_name, atom_type)
                except Exception, e:
                        #print out infomration from python
                        print e
                        parse_wetspot_cfg_error(name,value)
                        continue

        #print ' wetspot_id: %s interaction name %s defintion %i %s %s' % (wetspot_def, interaction_id, res_id, res_name, res_atomtype)

                if not wetspot.howmany_interactions() == 0:
                        print 'Wetspot loaded Id:', wetspot.id,' Defined by ',wetspot.howmany_interactions(), ' points'
                        wetspots[wetspot_def] = wetspot

        return wetspots


'''
##################################
############# PROCESS ############
##################################
'''


def trj_to_ddbb(cfg,wetspots):
    ''' For each frame find waters with hbonds in the wetspot and generate the centroid of the wetspot'''
        
    totalframes = cfg.total_frames + 1
    start = cfg.first_frame
    load_centroids(wetspots,cfg)    
    for frame in range(start,totalframes):
        #Find waters with hbonds
        parse_hbplus_to_ddbb(frame,wetspots,cfg)
        #Calculate the centroid of each wetspot

        #calc_centroids(frame,wetspots,cfg)

        #Calculate the distance of each water to each centroid
        calc_waters_parms(frame,wetspots,cfg)
        
    return

def parse_hbplus_to_ddbb(frame,wetspots,cfg):
    '''Look for interaction waters from hbplus and write to ddbb'''   
     
    hbfiles = cfg.hb_filename
    

    hb_file = cfg.hb_filename+'.%s' % (frame)
    
    hbplus_dic = read_HBplus(hb_file)
    
    for wetspot_key,wetspot_data in wetspots.iteritems():
        for interaction in wetspot_data.interactions:
            
            for hb in hbplus_dic.itervalues():
                
                water_info = None
                water_info = get_interaction_waters(hb,interaction,wetspot_key) 
                
                if not water_info == None:
                    
                    water_info['coord_x'],water_info['coord_y'],water_info['coord_z'] = get_atom_coord(water_info['res_id'],'O',frame,cfg)
                    
                    z_ddbb.write_interaction_water(water_info,cfg)



def get_hb_selected_waters(wetspot_key,cfg):

    (c,db) = z_ddbb.link_db(cfg)
    totalframes = cfg.total_frames + 1
    filename = 'HB-waters_'+wetspot_key

    for frame in range(1,totalframes):

        add2file(filename,frame)

        query = c.execute('SELECT water_id FROM water_intervals WHERE start_interval<=? AND end_interval>=? AND wetspot_id=?', (frame,frame,wetspot_key,)
        water_id = query.fetchone()
        if water_id != None:
            interaction = {'res_id':int(water_id),'atom_type':'O'}
            for i in parse_hbplus_to_txt(frame,wetspot_key,cfg,interaction)
                if i != None:

                    content = '\t'+' '+str(i['res_id'])+'  '+str(i['res_name'])+' '+str(i['atom_type'])

                    add2file(filename,content)
                    
            
            add2file(filename,'\n')
    return 



def parse_hbplus_to_txt(frame,wetspot_key,cfg,interaction):

    hbfiles = cfg.hb_filename
    

    hb_file = cfg.hb_filename+'.%s' % (frame)
    
    hbplus_dic = read_HBplus(hb_file)
    
    for hb in hbplus_dic.itervalues():
                
        water_info = None
        water_info = get_interaction_res(hb,interaction,wetspot_key) 

        if not water_info == None:
            
            yield water_info 
                   

def get_interaction_res(hb,interaction,wetspot_key):
    '''Once the HBplus output it has been parsed to a dict, this routine found which waters have hydrogen bounds 
    with the interactions atoms of the wetspot and return a Dictionary {interaction_id, wetspot_id, water_id, water_name, 
    frame, bound distance, bound angle }'''

    step_info = {}
    
    if hb['donor_res'] ==  hb['accep_res']:
        
        step_info = None
        return step_info
        
    elif hb['donor_res_id'] == interaction['res_id']:
        
        if hb['donor_res_atomtype'] == interaction['atom_type']:
            
            if hb['accep_res'] != 'WAT':
                        
                        step_info = {\
                                                          'interaction_id':str(interaction['interaction_id']),
                                                          'wetspot_id':str(wetspot_key),
                                                          'res_id':hb['accep_res_id'],
                                                          'res_name':hb['accep_res'],
                                                          'frame':hb['frame'],
                                                          'DAdistance':hb['DAdistance'],
                                                          'DAangle':hb['DAangle'],
                                                          'atom_type':hb['accep_res_atomtype']} 
         
                               
                        
                        
    elif hb['accep_res_id'] == interaction['res_id']:
        
        if hb['accep_res_atomtype'] == interaction['atom_type']:
            
            if hb['donor_res'] == 'WAT':
                        
                        step_info = {\
                                                          'interaction_id':str(interaction['interaction_id']),
                                                          'wetspot_id':str(wetspot_key),
                                                          'res_id':hb['donor_res_id'],
                                                          'res_name':hb['donor_res'],
                                                          'frame':hb['frame'],
                                                          'DAdistance':hb['DAdistance'],
                                                          'DAangle':hb['DAangle'],
                                                          'atom_type':hb['donor_res_atomtype']}  
                        
                        
                        
    
    
    if step_info.has_key('interaction_id'):
        return step_info
    else:
        
        step_info = None
        return step_info

def get_interaction_waters(hb,interaction,wetspot_key):
    '''Once the HBplus output it has been parsed to a dict, this routine found which waters have hydrogen bounds 
    with the interactions atoms of the wetspot and return a Dictionary {interaction_id, wetspot_id, water_id, water_name, 
    frame, bound distance, bound angle }'''

    step_info = {}
    
    if hb['donor_res'] ==  hb['accep_res']:
        
        step_info = None
        return step_info
        
    elif hb['donor_res_id'] == interaction['res_id']:
        
        if hb['donor_res_atomtype'] == interaction['atom_type']:
            
            if hb['accep_res'] == 'WAT':
                        
                        step_info = {\
                                                          'interaction_id':str(interaction['interaction_id']),
                                                          'wetspot_id':str(wetspot_key),
                                                          'res_id':hb['accep_res_id'],
                                                          'res_name':hb['accep_res'],
                                                          'frame':hb['frame'],
                                                          'DAdistance':hb['DAdistance'],
                                                          'DAangle':hb['DAangle']} 
         
                               
                        
                        
    elif hb['accep_res_id'] == interaction['res_id']:
        
        if hb['accep_res_atomtype'] == interaction['atom_type']:
            
            if hb['donor_res'] == 'WAT':
                        
                        step_info = {\
                                                          'interaction_id':str(interaction['interaction_id']),
                                                          'wetspot_id':str(wetspot_key),
                                                          'res_id':hb['donor_res_id'],
                                                          'res_name':hb['donor_res'],
                                                          'frame':hb['frame'],
                                                          'DAdistance':hb['DAdistance'],
                                                          'DAangle':hb['DAangle']}  
                        
                        
                        
    
    
    if step_info.has_key('interaction_id'):
        return step_info
    else:
        
        step_info = None
        return step_info

def calc_centroids(frame,wetspots,cfg):
    '''Find centroids for each wetspot for this frame and write to the ddbb''' 
    
    (c,db) = z_ddbb.link_db(cfg)
    
    try:
        
        for wetspot_key,wetspot_data in wetspots.iteritems():
           
    
            coord_x = 0.0
            coord_y = 0.0
            coord_z = 0.0
            num_atoms = 0

            try:
                for interaction_spot in wetspot_data.contacts:
                
                        coord = get_atom_coord(interaction_spot['res_id'],interaction_spot['atom_type'],frame,cfg)
                        coord_x += coord[0]
                        coord_y += coord[1]
                        coord_z += coord[2]
               
                        num_atoms += 1

           
                coord_x = coord_x/num_atoms
                coord_y = coord_y/num_atoms
                coord_z = coord_z/num_atoms
           
                c.execute('INSERT INTO centroids VALUES (null, ?,?,?,?,?)',(wetspot_key,frame,coord_x,coord_y,coord_z))
                db.commit()
           
            except Exception, e:

                print "ERROR!! ", e
                print "Getting the centroids in frame ", frame,'for wetspot ', wetspot_key
                print coord_x,coord_y,coord_z, ' coordinations'
                print num_atoms, ' atoms loaded'
                pass

    except Exception, e:
        
        print "General ERROR!! ", e
        print "Getting the centroids in frame ", frame,'for wetspot ', wetspot_key

    
    db.commit()

    db.close()
def load_centroids(wetspots,cfg):
    '''load water X-ray molecules coord into the database'''
    (c,db) = z_ddbb.link_db(cfg)
    
    try:
        
        for wetspot_key,wetspot_data in wetspots.iteritems():
            c.execute('INSERT INTO centroids VALUES (null, ?,?,?,?)',(wetspot_key,wetspot_data.coord_x,wetspot_data.coord_y,wetspot_data.coord_z))
            db.commit()


    except Exception, e:
        
        print "General ERROR!! ", e
        print "Getting the centroids in frame for wetspot ", wetspot_key

    return

def calc_waters_parms(frame,wetspots,cfg):
    '''Calc parameter to next filter 
    Look for waters in ddbb that are within LIMIT to Centroid of wetspot
    there is the point to add new parameters'''
    
    
    (c,db) = z_ddbb.link_db(cfg)

    for wetspot_key,wetspot_data in wetspots.iteritems():
           

        try:
            #Get coordinations of the centroid for each frame    
            #for each water is interacting in the wetspot get the id and the number of interactions that it has
            waters = {}

            query = c.execute('SELECT centroids.coord_x,centroids.coord_y,centroids.coord_z,interaction_waters.water_id,interaction_waters.coord_x,\
                                interaction_waters.coord_y,interaction_waters.coord_z FROM  centroids LEFT JOIN interaction_waters \
                                WHERE centroids.wetspot_id=? AND interaction_waters.wetspot_id=? AND \
                                interaction_waters.frame=?', (wetspot_key,wetspot_key,frame,))
                        
            for field in  query:
                        #field = [CENTROID coord_x],[CENTROID coord_y],[CENTROID coord_y], [water_id], [WATER coord_x],[WATER coord_y],[WATER coord_y]
                d = euclidian_distance(field[0],field[1],field[2],field[4],field[5],field[6])
                        
                if field[3] in waters:
                        waters[field[3]][1] += 1
                else:
                        waters[field[3]] = [d,1]

            for water,data in waters.iteritems():
                
                # distance of the water to the centroid / Number of interactions
                
                c.execute('INSERT INTO water_filter_param VALUES (null, ?,?,?,?,?,?)',(wetspot_key,water,frame,data[0],data[1],float(data[0])/float(data[1]),))
                
                db.commit()
        
        except Exception, e:
            
            print "ERROR!! ",e
            print "\t Exception when I try to found waters with hydrogen bonds for the wetspot ", wetspot_key, "in the frame number  ", frame
            print 'or query', query
            print 'may be water?,',field[3],' distance ', d, 'k-value ? ', k_value
            print 'dictionary  ', waters
            
    db.commit()
    db.close()



'''
##################################
############# ANALYSIS ############
##################################
'''

def get_water_intervals(wetspot_key,wetspot_data,cfg):
    '''setup the interval of a water in a wetspot'''
    '''to find a new water -->best ratio distance/number of interactions'''
    '''to say the water is still there --> cutoff distance'''
  
    #reSet variables of the tracker algorithm 
    (c,db) = z_ddbb.link_db(cfg)
    
    totalframes = cfg.total_frames + 1 
    interaction = wetspot_data.interactions
    water_id = None
    start_interval = cfg.first_frame

    for frame in range(start_interval,totalframes):
        
        '''Ask if we have a water, if not, search a new one'''
        '''If we have water from late frame, check if are in the cutoff distance from the centroid'''   
        ''' if the water is not , save the interval and  search a new one to next frame'''    
        '''if the average distance to the centroid it's bigger than 4, the wetspot is broke'''
        # CHECK GEOMETRIC of the wetspot... is it unfold? NEED IMPROVMENT MORE WORK (USE RMSD?)
        #if wetspot_geo(frame, wetspot_key) > 4.5:
                
        #        start_interval, water_id = save_interval(start_interval,frame,water_id,wetspot_data,cfg)
        #        continue                
        
        ###'''When a wetspot by just two ACC/DON it is important track the distance between them to be sure that a water is coordinating a HB'''
        #Get the interaction of the wetspot
        if len(interaction) == 2:
                interaction_spot = interaction[0]
                interaction_spot_2 = interaction[1]
                coor_1 = get_atom_coord(interaction_spot['res_id'], interaction_spot['atom_type'], frame,cfg)
                coor_2 = get_atom_coord(interaction_spot_2['res_id'], interaction_spot_2['atom_type'], frame,cfg)
                distance = euclidian_distance(coor_1[0],coor_1[1],coor_1[2],coor_2[0],coor_2[1],coor_2[2])
                #If the two interaction spots are too close, it is not space for a water molecule 
                
                if distance < 3.5:
                        start_interval, water_id = save_interval(start_interval,frame,water_id,wetspot_data,cfg)
                        continue


        if water_id == None:
                water_id = search_new_water(frame,wetspot_data,cfg)
                start_interval = frame
                continue
            
        else: 
                if not check_water(frame,wetspot_data,water_id,cfg):
                #if water_id != search_new_water(frame,wetspot_data,cfg)
                        save_interval(start_interval,frame,water_id,wetspot_data,cfg)
                        water_id = search_new_water(frame,wetspot_data,cfg)
                        start_interval = frame
                        continue

    ''' Save the last water '''
     
    save_interval(start_interval,frame,water_id,wetspot_data,cfg)
    
    '''delete  dry intervals'''
    c.execute('DELETE FROM water_intervals WHERE water_id is null')
    c.execute('DELETE FROM water_intervals WHERE occ_time = 0')

    
    db.commit()
    
    db.close()
    
    return 

def save_interval(start_interval,now_frame,water_id,wetspot_data,cfg):

        (c,db) = z_ddbb.link_db(cfg)
        end_interval = now_frame -1
        interval_steps = int(now_frame) - int(start_interval)
        
        time = interval_steps * cfg.timeXframe

        c.execute('INSERT INTO water_intervals VALUES (null, ?,?,?,?,?,?)',(wetspot_data.id,wetspot_data.cutoff,water_id,start_interval, end_interval,time))

        db.commit()
        #Reboot the flags
        water_id = None
        start_interval = now_frame

        return start_interval, water_id 


def check_water(frame,wetspot_data, water_id,cfg):
    '''check in the pdb if the water is beyond the cutoff centroid and if is just in this wetspot or could be better in other
    and check if the water molecule is only in this wetspot'''
    
    try:
        #centroid_coord = get_centroid_frame(frame,wetspot_data.id,cfg)
        centroid_coord = [wetspot_data.coord_x,wetspot_data.coord_y,wetspot_data.coord_z]
        water_coord = get_atom_coord(water_id,'O', frame, cfg)
    
        distance_to_centroid =  euclidian_distance(centroid_coord[0],centroid_coord[1],centroid_coord[2],water_coord[0],water_coord[1],water_coord[2])
    
        if distance_to_centroid > wetspot_data.cutoff:
        
                return False

        if shared_water(frame,water_id,wetspot_data.id,cfg):

                return False

        else:
                return True
     
    except Exception, e:
        print 'EROOR! ',e
        print 'coor centroid ', centroid_coord, 'Frame', frame
        print 'ERROR checking water  ', water_id, 'coor', water_coord, 'distance ', distance_to_centroid

def shared_water(frame,water_id,wetspot_id,cfg):
    ''''If a water is more close to other wetspot, no accept the water molecule '''

    (c,db) = z_ddbb.link_db(cfg)   
    c.execute('SELECT water_id,distance,wetspot_id from water_filter_param where  frame=? and water_id=? GROUP BY distance', (frame,water_id,))
    frame_waters = c.fetchall()
    
    try: 
        if wetspot_id != frame_waters[0][2]:
            return True
        else:
            return False
    except:
            return False


def search_new_water(frame,wetspot_data,cfg):
    
    '''Get from the ddbb for a frame the best water, best ratio between distance and number of interaction'''
    ''' if are not any water below the cutoff return 0 as water_id'''
    
    #look in the pdb to find water molecules.....TAG to active or desactive
    if  cfg.intesive:
        water_id = look_at_pdb(frame,wetspot_data,cfg)
    
        return water_id


    (c,db) = z_ddbb.link_db(cfg)
    c.execute('SELECT water_id,distance from water_filter_param where wetspot_id=? and frame=? GROUP BY distance', (wetspot_data.id,frame,)) #k-value option
    waters = c.fetchall()
    
    db.close()
    
    water_id = None
    
    for water in waters:
        
        if water[1]> wetspot_data.cutoff:
            water_id = None
        
        elif shared_water(frame,water[0],wetspot_data.id,cfg):
            water_id = None
        
        else:
            return water[0] 
            break
    



'''
##################################
############# DDBB ANALYSIS ############
##################################
'''

def get_wetspots_stats(wetspots,cfg):
            
    totalframes = cfg.total_frames
    timeXframe = cfg.timeXframe

    (c,con) = z_ddbb.link_db(cfg)
    
    
    #jobtime('Wetspots stats start\n')
  
    for wetspot_key,wetspot_data in wetspots.iteritems():
        
        
        try:
            
            stats_data = []
            #[0] Total time where the wetspot is wet
            c.execute('select SUM(occ_time) from water_intervals where wetspot_id=?',(wetspot_key,))
            stats_data.extend(c.fetchone())

            #[1] Max continuous Time that one single water molecule was in the wetspot
            c.execute('select MAX(occ_time) from water_intervals where wetspot_id=?',(wetspot_key,))
            stats_data.extend(c.fetchone())

            #[2] Avg Time of the water molecules #[3] Standard desviation    
            c.execute('select occ_time from water_intervals where wetspot_id=?',(wetspot_key,))
            mean,std = meanstdv([e[0] for e in c])
            stats_data.append(mean)
            stats_data.append(std)

            #[4] Occupancy (Rate of trajectory time where the wetspot is occupied)
            stats_data.append(float(stats_data[0])/float(totalframes*timeXframe)*100)

            #[5] total number of water molecules that where in the wetspot during the trajectory  
            c.execute('select water_id from water_intervals where wetspot_id=?',(wetspot_key,))
            stats_data.append(len(c.fetchall()))
                
            #SAVE DATA in the DDBB 
            c.execute('INSERT INTO stats_wetspots VALUES (null, ?,?,?,?,?,?,?,?)',(wetspot_key,wetspot_data.cutoff,stats_data[0],
                                                                                stats_data[1],stats_data[2],stats_data[3],stats_data[4],
                                                                                stats_data[5]))

            #Now stats, where the high speed water molecules are not take it in the calc

            stats_data = []

            time_threshold = cfg.timeXframe * 2

            #[0] Total time where the wetspot is wet
            c.execute('select SUM(occ_time) from water_intervals where wetspot_id=? and occ_time>?',(wetspot_key,time_threshold,))
            stats_data.extend(c.fetchone())

            #[1] Max continuous Time that one single water molecule was in the wetspot
            c.execute('select MAX(occ_time) from water_intervals where wetspot_id=? and occ_time>?',(wetspot_key,time_threshold,))
            stats_data.extend(c.fetchone())

            #[2] Avg Time of the water molecules #[3] Standard desviation    
            c.execute('select occ_time from water_intervals where wetspot_id=? and occ_time>?',(wetspot_key,time_threshold,))
            mean,std = meanstdv([e[0] for e in c])
            stats_data.append(mean)
            stats_data.append(std)

            #[4] Occupancy (Rate of trajectory time where the wetspot is occupied)
            stats_data.append(float(stats_data[0])/float(totalframes*timeXframe)*100)

            #[5] total number of water molecules that where in the wetspot during the trajectory  
            c.execute('select water_id from water_intervals where wetspot_id=? and occ_time>?',(wetspot_key,time_threshold,))
            stats_data.append(len(c.fetchall()))
                
            #SAVE DATA in the DDBB 
            c.execute('INSERT INTO stats_wetspots_slow VALUES (null, ?,?,?,?,?,?,?,?)',(wetspot_key,wetspot_data.cutoff,stats_data[0],
                                                                                stats_data[1],stats_data[2],stats_data[3],stats_data[4],
                                                                                stats_data[5]))
        

        except Exception, e:
            print e
            print 'Error while I am getting the stadistics from the wetspot  ', wetspot_key, ' - it could be dry - '
            print stats_data
            c.execute('INSERT INTO stats_wetspots VALUES (null, ?,?,?,?,?,?,?,?)',(wetspot_key,wetspot_data.cutoff,0,0,0,0,0,0))
        
        con.commit()


        try:
            
            #DATA for a Histogram
            query = c.execute('SELECT occ_time FROM water_intervals WHERE wetspot_id=?', (wetspot_key,)) 

            histogram = {}
            for e in query:
                histogram[e[0]] = histogram.get(e[0],0) + 1
                edf = Edfreq(histogram,wetspot_key)
            #Save Histogram data to DDBB
            

            c.execute('INSERT INTO wetspots_histogram VALUES (null, ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(edf.name,
                                                                                                wetspot_data.probs[0],edf.prob(wetspot_data.probs[0]),
                                                                                                wetspot_data.probs[1],edf.prob(wetspot_data.probs[1]),
                                                                                                wetspot_data.probs[2],edf.prob(wetspot_data.probs[2]),
                                                                                                wetspot_data.probs[3],edf.prob(wetspot_data.probs[3]),
                                                                                                wetspot_data.probs[4],edf.prob(wetspot_data.probs[4]),
                                                                                                wetspot_data.probs[5],edf.prob(wetspot_data.probs[5]),
                                                                                                wetspot_data.probs[6],edf.prob(wetspot_data.probs[6]),
                                                                                                wetspot_data.probs[7],edf.prob(wetspot_data.probs[7]),
                                                                                                wetspot_data.probs[8],edf.prob(wetspot_data.probs[8]),
                                                                                                wetspot_data.probs[9],edf.prob(wetspot_data.probs[9]),
                                                                                                edf.value(0.25),edf.value(0.5),edf.value(0.75)))
            #Save a file with all the data of the Emperical distribution of the time.
            filename = str(wetspot_key+'.edf.dat')
            edf.render(filename)


        except Exception, e:

            print "Error %s:" % e
            print 'ERROR! while  Calc to histogram  of the wetspot  ', wetspot_key 
            print histogram
            #print edf
            #c.execute('INSERT INTO stats_wetspots VALUES (null, ?,?,?,?,?,?,?,?)',(wetspot_key,cutoff_mod,0,0,0,0,0,0))
        
        
        con.commit()
        
        
       
    con.close() 



    #jobtime('Wetspots stats end\n')

def get_interactions_stats(wetspots,cfg):

        
    (c,con) = z_ddbb.link_db(cfg)
    
    total_frames = cfg.total_frames
    timeXframe = cfg.timeXframe
    start = cfg.first_frame

    #jobtime('Interaction stats start\n')

    for wetspot_key,wetspot_data in wetspots.iteritems():
    
        try:
                hb_data = {}
                c.execute('SELECT SUM(occ_time) from water_intervals where wetspot_id=?',(wetspot_key,))
                wet_time = c.fetchone()

                for frame in range(start,total_frames):

                        query = c.execute('SELECT interaction_waters.interaction_id,interaction_waters.distance,\
                                interaction_waters.angle FROM  interaction_waters LEFT JOIN water_intervals WHERE\
                                water_intervals.start_interval<=? AND water_intervals.end_interval>=? AND interaction_waters.frame=?\
                                AND interaction_waters.water_id = water_intervals.water_id AND\
                                water_intervals.wetspot_id = ? AND interaction_waters.wetspot_id = water_intervals.wetspot_id',\
                                (frame,frame,frame,wetspot_key,))

                        for interaction in query:
                                #hb_data (type Dictionary) where index is interaction_id, value (type list) [ [distances], [angles] ]
                                if interaction[0] in hb_data:
                                        hb_data[interaction[0]][0].append(interaction[1])
                                        hb_data[interaction[0]][1].append(interaction[2])
                                else:
                                        hb_data[interaction[0]] = [[interaction[1]], [interaction[2]]]                              

                for hb,data in hb_data.iteritems():
                        #Calc Stats
                        distance_mean, distance_std = meanstdv(data[0])
                        angle_mean, angle_std = meanstdv(data[1])
                        #Percentage of wet time (time where) this hb is established

                        time = len(data[0])*cfg.timeXframe/float(wet_time[0])*100

                        #SAVE interaction
                        c.execute('INSERT INTO stats_interactions VALUES (null, ?,?,?,?,?,?,?)',(wetspot_key,hb,time,distance_mean, distance_std,angle_mean, angle_std))

                        con.commit()
        
        except Exception, e:
    
            print "Error %s:" % e
            print 'Error Interaction STATS at wetspot ', wetspot_key
            pass



    con.close()

    return


def get_centroid_frame(frame,wetspot_key,cfg):
    
    '''Search Coord for a Centroid in one Frame'''
    
    centroid_coord=[]
    
    (c,db) = z_ddbb.link_db(cfg)
    
    c.execute('SELECT coord_x,coord_y,coord_z FROM centroids WHERE frame=? AND wetspot_id=?' ,(frame,wetspot_key,))
    
    for row in c:
        for coord in row:
            
            centroid_coord.append(coord)
    
    db.close()
    
    return centroid_coord

def cluster_edf(wetspot_key,cfg):
    
    '''Experimental CLuser frequency'''
    
    cluster = []
    time = cfg.timeXframe * cfg.total_frames
    centroids = [10,20,30,40,50,60,70,80,90,100,200,500,1000]

    for i in range(13):
        cluster.append([])

    (c,db) = z_ddbb.link_db(cfg)
    
    c.execute('SELECT occ_time FROM water_intervals WHERE wetspot_id=?' ,(wetspot_key,))
    
    for e in c:

        #print int(e[0])
        
        #if int(e[0])<5:

            #cluster[0].append(int(e[0]))
            #continue
        if int(e[0])<11:
            cluster[0].append(int(e[0]))
            continue
        if int(e[0])<21:
            cluster[1].append(int(e[0]))
            continue
        if int(e[0])<31:
            cluster[2].append(int(e[0]))
            continue
        if int(e[0])<41:
            cluster[3].append(int(e[0]))
            continue
        if int(e[0])<51:
            cluster[4].append(int(e[0]))
            continue
        if int(e[0])<61:
            cluster[5].append(int(e[0]))
            continue
        if int(e[0])<71:
            cluster[6].append(int(e[0]))
            continue
        if int(e[0])<81:
            cluster[7].append(int(e[0]))
            continue
        if int(e[0])<91:
            cluster[8].append(int(e[0]))
            continue
        if int(e[0])<101:
            cluster[9].append(int(e[0]))
            continue
        if int(e[0])<201:
            cluster[10].append(int(e[0]))
            continue

        if int(e[0])<501:
            cluster[11].append(int(e[0]))
            continue

        else:
            cluster[12].append(int(e[0]))
            
            continue
            


    c.execute('SELECT SUM(occ_time) FROM water_intervals WHERE wetspot_id=?' ,(wetspot_key,))
    wet_time = c.fetchone()

    #print wet_time
    filename = cfg.ddbb+'_'+wetspot_key
    content= '#######    '+wetspot_key+'\n'
    add2file(filename,content)
    for i in range(13):
        try:
            
            N = float(wet_time[0])/float(centroids[i])
            content=str(centroids[i])+'  '+str(float(sum(cluster[i]))/float(wet_time[0]))+'\n'
            add2file(filename,content)
        except:
            content= '  0.0'+'\n'
            add2file(filename,content)
         
    
    db.close()
    
    return


'''
##################################
############# DUMPING PDB's ############
##################################
'''



def dump_pdbs(cfg,wetspots,start_frame=1,end_frame='',wetspot_id=''):
    ''' generate pdb files just with the waters molecues seleted by ZAHORI
    By Default, all the frames (filename dry.pdb.frame_id) 
                all the wetspots
    '''
    #IMPRO wetpots selection, remember to change filename of the pdb file acording the wetpot.
    #IMPRO default name for the wetpots
    try:

        if end_frame == '':
            totalframes = cfg.total_frames + 1

        
        
        (c,db) = z_ddbb.link_db(cfg)
    
    
        for frame in range(start_frame,totalframes):
            
            query = c.execute('SELECT DISTINCT water_id,wetspot_id from  water_intervals where start_interval<=? and end_interval>=?', (frame,frame,))
        
            dic_of_waters = {}
            for e in query:
                dic_of_waters[e[1]] = e[0]
        
            write_new_pdb(cfg,frame, dic_of_waters,wetspots)
    
        db.close()
    
    except Exception, e:
        print 'Error ', e
        print 'Error dumping files  ', frame #IMPRO handle the error, print more information
    
    
    return
    


def write_new_pdb(cfg,frame,dic_of_waters,wetspots,wetspot_id=''):
    '''Simple routine to get pdb for each frame with the selectec waters;
    copy the pdb files from the trajectory and remove the waters not selected by ZAHORI'''
    
    #IMPRO .. add a flag or a way to dump water from just one wetspot.

    try:

        filename =cfg.frame_pdbname+'.'+str(frame) #IMPRO this is a default name, 
        f = open(filename,'r')
        pdb_file = f.readlines()
        f.close()
        water_list = [id for id in dic_of_waters.itervalues()]
        
        new_pdb_file='dump'+wetspot_id+'.pdb.'+str(frame)
        last_line = ''
        num_of_waters = 0
        for atom in pdb_file:
            if not atom[17:20].strip() == 'WAT' or int(atom[22:27].strip()) in water_list:
                if not last_line == atom:
                    if atom[17:20].strip() == 'WAT':
                        atom = change_resid(atom,dic_of_waters)
                    add2file(new_pdb_file,atom)
                    last_line = atom


        #to load in VMD as a traj the number of atoms must be the same all the time, this method new IMPROVMENT
        #if len(water_list) < len(wetspots):
                #Add dummy atoms, that is need it to load the pdbs in to vmd like a trj
        #        generate_dummy_atom(new_pdb_file,len(wetspots) - len(water_list))


     
    except Exception, e:
        print e
        print 'Error dumping the pdb file ', new_pdb_file #add or improve this error segment IMPRO
                
    return

def change_resid(atom,dic_of_waters):
    ''' get the pdb line as string and change the resid for the wetspot name'''
    for wet_id, water_id in dic_of_waters.iteritems():
        if int(atom[22:27].strip()) == int(water_id):
            wetspot_id = wet_id
            wid = water_id
            break

    
    ATOM_STR = "%-6s%5i %c%-3s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"       
    atm = ATOM_STR % ('ATOM',int(atom[6:11].strip()),' ',atom[12:16].strip(),' ',wetspot_id[-3:],' ',int(wid),' ',float(atom[31:39].strip()),float(atom[39:47].strip()),float(atom[47:55].strip()),float('0.00'),float('0.00'),'O','O')

    return atm

def generate_dummy_atom(new_pdb_file,n):

        ATOM_STR = "%-6s%5i %4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"


        while n != 0:
                #only valid to 9 waters. after may be problem in parsing the res_id IMPROV
                
                dummy_line_1 = 'ATOM   999%i  O   DUM     %00i       0.000   0.000   0.000  0.00  0.00\n'  % (n,n)
                add2file(new_pdb_file,dummy_line_1)
                dummy_line_2 = 'ATOM   999%i  H1  DUM     %00i       0.000   0.000   0.000  0.00  0.00\n'  % (n,n) 
                add2file(new_pdb_file,dummy_line_2)
                dummy_line_3 = 'ATOM   999%i  H2  DUM     %00i       0.000   0.000   0.000  0.00  0.00\n'  % (n,n)
                add2file(new_pdb_file,dummy_line_3)
                end_dummy =    'TER       0              0\n'
                add2file(new_pdb_file,end_dummy)
                n -= 1

        return 


def found_wetspots(filename,cfg):
    '''Load pdb file, and define a the wetspot for each water'''

        
    f = open(filename,'r')
    pdb_file = f.readlines()
    f.close()
    waters = {}
    #look for waters in the pdb
    for atom in pdb_file:
        if  atom[17:20].strip() == 'WAT' or atom[17:20].strip() == 'HOH' and atom[12:16].strip() == 'O':
            coord = []
            coord.append(float(atom[31:39].strip()))
            coord.append(float(atom[39:47].strip()))
            coord.append(float(atom[47:55].strip()))
            #Save id water like index, and list with his coor.
            waters[int(atom[23:30].strip())] = coord

    
    #find the neightboors of these waters
    for water_id,coord in waters.iteritems():
        atoms = found_within(filename,cfg.cutoff_1,cfg.cutoff_2,coord)
        
        write_defintion(cfg.wetspots_file,water_id,waters[water_id],atoms)

    return

def found_within(filename,cutoff1,cutoff2,centroid_coord):

    structure = read_pdb(filename)
    coord = [0.0,0.0,0.0]
    water_names = ['HOH','WAT']
    don_acc = ['N','NE','NH1','NH2','ND1','ND2','NE1','NE2','NZ','SG','OH','OH2','OG','OG1','O','OD1','OD2','OE1','OE2']
    heavy_atom = ['CB','C','CA','CE','CE1','CE2','CD','CD1','CD2','CG','CG1','CG2','CZ','CZ1','CZ2','SD']
    heavy_atom.extend(don_acc)
    n = 0

    atoms = {}
    try:    
        for line in structure:
            if line[0:4] == 'ATOM':
                if not line[17:20].strip() in water_names:
                    if str(line[12:16].strip()) in heavy_atom:
                

                        coord[0] = float(line[31:39].strip())
                        coord[1] = float(line[39:47].strip())
                        coord[2] = float(line[47:55].strip())


                        if ((coord[0] - centroid_coord[0])**2 < cutoff2**2) and ((coord[1] - centroid_coord[1])**2 < cutoff2**2) and ((coord[2] - centroid_coord[2])**2 < cutoff2**2):
                                
                                distance = euclidian_distance(centroid_coord[0],centroid_coord[1],centroid_coord[2],coord[0],coord[1],coord[2])
                                
                                if cutoff1 > distance:
                                    res_id = str(line[23:30].strip())
                                    res_name = str(line[17:20].strip())
                                    atom_type = str(line[12:16].strip())
                                    n +=1

                                    if line[12:16].strip() in don_acc:

                                        atoms['W'+str(n)+res_id] = {'res_id':res_id,'res_name':res_name,'atom_type':atom_type}
                                        continue
                                    else:
                                        atoms['C'+str(n)+res_id] = {'res_id':res_id,'res_name':res_name,'atom_type':atom_type}
                                        continue

                                if cutoff2 > distance:
                                        res_id = line[23:30].strip()
                                        res_name = line[17:20].strip()
                                        atom_type = line[12:16].strip()
                                        n +=1

                                        atoms['C'+str(n)+res_id] = {'res_id':res_id,'res_name':res_name,'atom_type':atom_type}






    except Exception, e:
        print e
        print "ERROR!! Something happing while try to find atoms in the pdb  ", filename,


    return atoms

def write_defintion(wetfile,water_id,water_coords,atoms):

    name = '[WAT'+str(water_id)+']\n'
    add2file(wetfile,name)
    coords = 'coord='+str(water_coords[0])+','+str(water_coords[1])+','+str(water_coords[2])+'\n'
    add2file(wetfile,coords)

    for atom,atom_info in atoms.iteritems():
        point = str(atom)+'='+str(atom_info['res_id'])+'-'+str(atom_info['res_name'])+'-'+str(atom_info['atom_type'])+'\n'
        add2file(wetfile,point)
    
    point='\n'
    add2file(wetfile,point)

    return




    return
'''
##################################
############# HELP ROUTINEs ###############
##################################
'''

def euclidian_distance(x,y,z,a,b,c):

        '''Calc euclidian distance in 3D space, and return float'''

        return math.sqrt((x-a)**2+(y-b)**2+(z-c)**2)

def meanstdv(x):
    """ Calculate mean and standard deviation of data x[]:
        mean = {\sum_i x_i \over n}
        std = sqrt(\sum_i (x_i - mean)^2 \over n-1) """

    if len(x) == 1:
        return x[0],0
    
    n, mean, std = len(x), 0, 0
    
    for a in x:
        mean = mean + a
    
    mean = mean / float(n)
    
    for a in x:
        std = std + (a - mean)**2
    
    std = math.sqrt(std / float(n-1))
    
    return mean, std

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

def look_at_pdb(frame,wetspot_data,cfg):

    water_def = ['WAT','HOH']
    centroid_coord = [wetspot_data.coord_x,wetspot_data.coord_y,wetspot_data.coord_z]  
    pdb_file = cfg.frame_pdbname+'.'+str(frame)
    structure = read_pdb(pdb_file)
    water_id = None
    min_distance = wetspot_data.cutoff 
    coord = [0.0,0.0,0.0]
    try:    
        for line in structure:
            if line[17:20]  in water_def:
                if line[12:16].strip() == 'O':
                        coord[0] = float(line[31:39].strip())
                        coord[1] = float(line[39:47].strip())
                        coord[2] = float(line[47:55].strip())

                        if ((coord[0] - centroid_coord[0])**2 < wetspot_data.cutoff**2) and ((coord[1] - centroid_coord[1])**2 < wetspot_data.cutoff**2) and ((coord[2] - centroid_coord[2])**2 < wetspot_data.cutoff**2):
                                
                                distance = euclidian_distance(centroid_coord[0],centroid_coord[1],centroid_coord[2],coord[0],coord[1],coord[2])
                                
                                if min_distance > distance:
                                        water_id = int(line[21:30].strip())
                                        min_distance = distance
        
        if  min_distance < wetspot_data.cutoff:
            return water_id
        else:
            return None
    except Exception, e:
        print e
        print "ERROR!! Something happing while try to find new waters in the pdb  ", pdb_file, " for the wetspot  ", wetspot_data.id



def check_wet_definitions(cfg,wetspots):
    '''check if all the atoms existin the pdb file '''

    for wetspot_key,wetspot_data in wetspots.iteritems():
        for atom in wetspot_data.contacts:
            if not exist(atom,cfg):
                print 'Warning!!! a atom in the wetspots defintion ', wetspot_key, 'did not found in the pdbfile ', cfg.frame_pdbname+'.'+str(cfg.first_frame)
                print 'ATOM  ',atom, ' the numeration between the pdb and the trajectory could change. Also pay special atention to Histidines name.'
                sys.exit()

def exist(atom,cfg):


    pdb_file = cfg.frame_pdbname+'.'+str(cfg.first_frame)
    structure = read_pdb(pdb_file)
    for line in structure:

        if line[21:30].strip() == str(atom['res_id']):
                if line[12:16].strip() == str(atom['atom_type']):
                    if line[17:20].strip() == str(atom['res_name']):
                        return True

    return False



def get_atom_coord(res_num, atom_type, frame, cfg):
    ''' get the coord of and atom and return in a list
        need the res_num and  atom_type, pdb file name
    '''
    
    
    #get the prefix of the pdb file to the frame
    pdb_file = cfg.frame_pdbname+'.'+str(frame)
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
                
    except Exception, e:
        print e
        print 'ERROR getting the coordinates for the atom of the residue ', res_num, atom_type
        print 'In the frame  ', frame

        
def read_pdb(pdb_file):
    ''' Read a pdb file to a list of list'''
    try:
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()
    
        return pdb

    except Exception, e:
        print e
        print 'ERROR reading pdb file  ', pdb_file




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
     
    except Exception, e:
        print e 
        print "ERROR!! Parsing HBplus file   ", hbplusfile
        return hbplusfile

                    
    
    return hb_dictionary

def call_HBplus(pdbfile, PATH='hbplus', DAdistance='3.5',DAangle='90'):
    ''' Generate a HB2 file from a pdb. hb2 file will have the same name that pdb'''
    
    try:
        os.system('hbplus -d %s -a %s  %s -h 2.5 >> hb_gen.log' % (DAdistance,DAangle,pdbfile))
        pdbname = pdbfile.split('.')
        os.system('mv %s.%s.hb2 %s.%s.hb2.%s' % (pdbname[0],pdbname[1],pdbname[0],pdbname[1],pdbname[2]))
    
    except Exception, e:
        print e
        print "ERROR! Calling HBplus to parse the pdb file  ", pdbfile
                  
    return


def getpdb(trajfile='trajfile',stepnum=0):
    
    '''Generate a ptraj file and call it to get a odb from frame'''
    ''' stepnum set the final frame from Totalframes variable'''
    
    
    trajfile = cfg.trj_filename
    stepnum =  cfg.total_frames
    path_ptraj = cfg.path_ptraj

    print 'Generating pdb.........'

    try:
        ptrajrmsd='trajin %s.mdcrd 0 %s\n  trajout  frame.pdb PDB\n' % (trajfile, stepnum)
        #print ptrajrmsd
        writefile('getPDB.ptraj',ptrajrmsd)
        os.system('ptraj %s.prmtop getPDB.ptraj' % (trajfile))
    except Exception, e:
        print e
        print 'ERROR generating pdb....'
      
    return

def gethb2(totalframes=1000,pdb_name='frame'):
    ''' bucle to call hbplus for each frame'''
    
    
    pdb_name = cfg.pdb_filename
    totalframes = cfg.total_frames + 1
    
    #Sum 1 to TOTALFRAMES about range method does not include the final number
    print 'Generating HB2 files......'

    try:
        
        for i in range(1,totalframes):
            file='%s.%s' % (pdb_name,i)
            call_HBplus(file)
    
    except Exception, e:
        print e
        print 'ERROR generating HB2 files......'

def time_lapse(t0, task):

        print task,' --> Done, Elapsed time from start: ',  time.time() - t0

        return

'''
##################################
############# MAIN ###############
##################################
'''

def main():

        
    
    parser = argparse.ArgumentParser(description="ZAHORI: XXXXXXXX")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--found_wetspots", action="store", dest="pdb_filename", help="Automatic creation of the wetspots description \
                                                                from a X-ray pdb file")
    group.add_argument("-p","--process", action="store", choices=['pdbs', 'hb2', 'all', 'ddbb'], help="Process the trajectory file and \
                                                                            creation of the waters database. CHOICES: pdb (Extract pdb files from trajectory),\
                                                                            hb2 (use HBplus to get Hydrogen bonds list from each pdb), \
                                                                            ddbb (Scan the trajectory and build the Waters database )  all ( performe all those analysis) By default all")
    group.add_argument("-a","--analysis", action="store_true", help="Extract  information from the Water database")
    parser.add_argument("-d","--dump", action="store", help="Dump a pdb file from each frame where \
                                                                            only wil be the water molecules selected by \
                                                                            ZAHORI. You must indicate the final frame or write \'all\' if you want to dump all the trajectory")

    parser.add_argument("-w","--wetspots", action="store", help="Wetspots description file, default name wetspots.inp, if does not exit the program stops")
    parser.add_argument("-s","--setup", action="store_const", const="cfg.inp", help="Configuration file, if does not exist, it will use default values")
    group.add_argument("--temp", action="store_true", help="")
    parser.add_argument("--hb", action="store", help="")



    #outputs name hb2 and pdbs.. helper... copy and work in a directory

    #Ussues of the ddbb
    #Add a la confugracion el tema del nombre de las moleculas de agua... flexibilizar

    #parser.add_argument
    

    if len(sys.argv[1:]) == 0:

        print "For help execute ZAHORI with -h optional argument"
        sys.exit()

    args = parser.parse_args()

    #LOAD CFG
    print 'LOADING Configuration'
    if args.setup:

        cfg = Cfg(args.setup)
    else:
        #Load default configuration

        cfg = Cfg()

    #ACTIONS

    #Found potential wetspots from a pdb file
    if args.pdb_filename:

        print 'Building wetspots definitions...'
        found_wetspots(args.pdb_filename,cfg)
        sys.exit()

    #Load description of wetspots (default file wetspots.inp)
    print 'LOADING Wetspots file description'
    if args.wetspots:

        wetspots = parse_wetspot_cfg(args.wetspots)
    else:
        #Get working folder
        print 'ERROR'
        print 'Wetspot definition need it'
        sys.exit()

    if args.temp:
        for wetspot_key,wetspot_data in wetspots.iteritems():
            cluster_edf(wetspot_key,cfg)

    if args.hb:
        get_hb_selected_waters(args.hb,cfg)
        sys.exit()


        #FROM TRAJ FILE TO DDBB

    #Convert trajectory file to a pdb files 
    #Get hydrogen bound information from pdb files
    #Load all the useful information to the ddbb ()
    t0 = time.time()

    if args.process:

        print 'Start processing trajectory'
        
        if args.process == 'pdb':getpdb()
        if args.process == 'hb2':gethb2()

        if args.process == 'ddbb':

            #Write in ddbb info : waters with hbond on the wetspot and wetspot centroids
            check_wet_definitions(cfg,wetspots)
            z_ddbb.create_ddbb(cfg,wetspots)

            trj_to_ddbb(cfg,wetspots)
            time_lapse(t0,'Trajectory to DDBB')

        if args.process == 'all':
            
            getpdb()
            time_lapse(t0,'Trajectory to PDB files')
            gethb2()
            time_lapse(t0,'From PDB files Hydrogen bond info')
            #Write in ddbb info : waters with hbond on the wetspot and wetspot centroids
            check_wet_definitions(cfg,wetspots)
            z_ddbb.create_ddbb(cfg,wetspots)
            trj_to_ddbb(cfg,wetspots)
            time_lapse(t0,'Trajectory to DDBB')


    #CALC WATERS INFO FROM DDBB

        
    if args.analysis:
            
        #Calculate for each wetspot the occupancie time of the waters, track waters a long trajectory

        for wetspot_key,wetspot_data in wetspots.iteritems():
            
            #jobtime(wetspot_key)
            
            #Here is the point to introduce a variable for a wetspots cutoff
            
            get_water_intervals(wetspot_key,wetspot_data,cfg)
        
        time_lapse(t0,'Determination of the intervals')
        #Calculate the regular stats for each wetspot
        get_wetspots_stats(wetspots,cfg)
        get_interactions_stats(wetspots,cfg)
        time_lapse(t0,'Stadistics')

    #DUMP a PDB file for each Frame with only water molecules selected by zahori 

    if args.dump:
        if args.dump == 'all' or args.dump == 'All':
                dump_pdbs(cfg,wetspots,start_frame=1,wetspot_id='all')
        else:
                dump_pdbs(cfg,wetspots,start_frame=int(args.dump),wetspot_id='all')

        time_lapse(t0,'Dump pdb files')
        


    #output program
    #print out wetfile name, num of frames, etc.


        
if __name__ == "__main__":

    
    main()


