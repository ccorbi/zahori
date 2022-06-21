'''
Created on 18/03/2010

@author: ccorbi
'''
#import pysqlite2
#from pysqlite2 import dbapi2 as sqlite3
import sqlite3
import sys, os


def link_thisdb(ddbb):

	try:
		db = sqlite3.connect(ddbb)
	except sqlite3.Error, errmsg:
		print 'DB not available ' + str(errmsg)
		sys.exit()
    	
	db_cursor = db.cursor()
        return db_cursor, db

def link_db(cfg):
    '''initializes the database file'''
    
    
    db_path = cfg.ddbb
    
    try:
        db = sqlite3.connect(db_path)
    except sqlite3.Error, errmsg:
        print 'DB not available ' + str(errmsg)
        sys.exit()

    db_cursor = db.cursor()
    return db_cursor, db

def create_ddbb(cfg,wetspots):
    
    
    db_path = cfg.ddbb
    
    if os.path.isfile(db_path):
        print 'Error, DATABASE already exist'
        sys.exit()
        
    con = sqlite3.connect(db_path)
    
    ''' table where store interaction spot '''
     
    con.execute("""CREATE TABLE interaction_spots (
               id INTEGER PRIMARY KEY AUTOINCREMENT,
               wetspot_id VARCHAR(50),
               interaction_id VARCHAR(50),
               res_id INTEGER,
               res_name VARCHAR(3),
               res_atomtype VARCHAR(3))""")
    
    ''' table where store waters that interact with some interaction spot '''
    
    con.execute("""CREATE TABLE interaction_waters (
               id INTEGER PRIMARY KEY AUTOINCREMENT,
               wetspot_id VARCHAR(50),
               interaction_id VARCHAR(50),
               water_id INTEGER,
               distance REAL,
               angle REAL,
               frame INTEGER,
               coord_x REAL,
               coord_y REAL,
               coord_z REAL,
               FOREIGN KEY (interaction_id) REFERENCES interaction_spots(interaction_id))""")
    
    ''' table where store centroid coordinates for each frame '''
    
    con.execute("""CREATE TABLE centroids (
               id INTEGER PRIMARY KEY AUTOINCREMENT,
               wetspot_id VARCHAR(50),
               coord_x REAL,
               coord_y REAL,
               coord_z REAL)""")
    
    ''' table where store distances and interaction for water for each frame and interactio spot '''
    
    con.execute("""CREATE TABLE water_filter_param (
       id INTEGER PRIMARY KEY AUTOINCREMENT,
       wetspot_id VARCHAR(50),
       water_id INTEGER,
       frame INTEGER,
       distance REAL,
       num_interactions INTEGER,
       k_value REAL)""")
    
    
    con.execute("""CREATE TABLE water_intervals (
               id INTEGER PRIMARY KEY AUTOINCREMENT,
               wetspot_id VARCHAR(50),
               cutoff REAL,
               water_id INTEGER,
               start_interval INTEGER,
               end_interval INTEGER,
               occ_time INTEGER)""")

    
    con.execute("""CREATE TABLE stats_wetspots (
           id INTEGER PRIMARY KEY AUTOINCREMENT,
           wetspot_id VARCHAR(50),
           cutoff REAL,
           wet_time INTEGER,
           max_time INTEGER,
           avg_time REAL,
           std_time REAL,
           occupancy REAL,
           num_waters INTEGER)""")

    con.execute("""CREATE TABLE stats_wetspots_slow (
           id INTEGER PRIMARY KEY AUTOINCREMENT,
           wetspot_id VARCHAR(50),
           cutoff REAL,
           wet_time INTEGER,
           max_time INTEGER,
           avg_time REAL,
           std_time REAL,
           occupancy REAL,
           num_waters INTEGER)""")

    con.execute("""CREATE TABLE stats_interactions (
               id INTEGER PRIMARY KEY AUTOINCREMENT,
               wetspot_id VARCHAR(50),
               interaction_id VARCHAR(50),
               time REAL,
               mean_dist REAL,
               std_dist REAL,
               mean_angle REAL,
               std_angle REAL)""")

    
    con.execute("""CREATE TABLE wetspots_histogram (
       id INTEGER PRIMARY KEY AUTOINCREMENT,
       wetspot_id VARCHAR(50),
       time_1 INTEGER,prob_1 REAL,
       time_2 INTEGER,prob_2 REAL,
       time_3 INTEGER,prob_3 REAL,
       time_4 INTEGER,prob_4 REAL,
       time_5 INTEGER,prob_5 REAL,
       time_6 INTEGER,prob_6 REAL,
       time_7 INTEGER,prob_7 REAL,
       time_8 INTEGER,prob_8 REAL,
       time_9 INTEGER,prob_9 REAL,
       time_10 INTEGER,prob_10 REAL,
       value_25 REAL,
       value_50 REAL,
       value_75 REAL)""")


    con.close()
    
    print 'Creating the DataBase ', cfg.ddbb
    
    write_interaction_spots(wetspots,cfg)
    
    return



''' ################## 
    ##################
    ##################
    ##################
    ##################
    HANDELING DATABASE '''




def write_interaction_water(water_info,cfg):
    

    
    (c,db) = link_db(cfg)
    
    c.execute('INSERT INTO interaction_waters VALUES (null, ?,?,?,?,?,?,?,?,?)',(water_info['wetspot_id'],
                                                                 water_info['interaction_id'],
                                                                 water_info['res_id'],
                                                                 water_info['DAdistance'],
                                                                 water_info['DAangle'],
                                                                 water_info['frame'],
                                                                 water_info['coord_x'],
                                                                 water_info['coord_y'],
                                                                 water_info['coord_z']))
            
    db.commit()        
    db.close()
    
def write_interaction_spots(wetspots,cfg):
        
    
    (c,db) = link_db(cfg)
      
       
    
    for wetspot_key,wetspot_data in wetspots.iteritems():
        for interaction in wetspot_data.interactions:
            c.execute('INSERT INTO interaction_spots VALUES (null, ?,?,?,?,?)',(wetspot_key,
                                                                                interaction['interaction_id'],
                                                                                interaction['res_id'],
                                                                                interaction['res_name'],
                                                                                interaction['atom_type']))
    db.commit()
    db.close()        
