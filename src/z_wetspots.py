'''
1BBZ
'''

wetspots = {}

interaction = {}
interac_1 = {}
interac_2 = {}
interac_3 = {}
interac_4 = {}
interac_5 = {}
interac_6 = {} #Wetspot inventado para usar de control
interac_7 = {}
interac_8 = {}
interac_9 = {}
interac_10 = {}
interac_11 = {}
interac_12 = {}
interac_13 = {}

### WATER 4 (MIN)
interaction[1] = {'interaction_id':'W4-GLU98:O','res_id':35,'res_name':'GLU','res_atomtype':'O'}
interaction[2] = {'interaction_id':'W4-SER113:N','res_id':50,'res_name':'SER','res_atomtype':'N'}
interaction[3] = {'interaction_id':'W4-ASN114:N','res_id':51,'res_name':'ASN','res_atomtype':'N'}

### WATER 4 (EXT)
interac_7[1] = {'interaction_id':'W4-GLU98:O','res_id':35,'res_name':'GLU','res_atomtype':'O'}
interac_7[2] = {'interaction_id':'W4-SER113:N','res_id':50,'res_name':'SER','res_atomtype':'N'}
interac_7[3] = {'interaction_id':'W4-ASN114:N','res_id':51,'res_name':'ASN','res_atomtype':'N'}
interac_7[4] = {'interaction_id':'W4-ASN114:OD1','res_id':51,'res_name':'ASN','res_atomtype':'OD1'}


### WATER 5 (NEW)

interac_5[1] = {'interaction_id':'W5-ASN114:ND','res_id':51,'res_name':'ASN','res_atomtype':'ND2'}
interac_5[2] = {'interaction_id':'W5-PRO6:O','res_id':66,'res_name':'PRO','res_atomtype':'O'}

### WATER 1 (FYN)
interac_11[1] = {'interaction_id':'W1-HIS95:O','res_id':32,'res_name':'HIS','res_atomtype':'ND1'}
interac_11[2] = {'interaction_id':'W1-PRO2:O','res_id':61,'res_name':'PRO','res_atomtype':'O'}


#### WATER 33 
interac_1[1] = {'interaction_id':'W33-LEU17:O','res_id':17,'res_name':'LEU','res_atomtype':'O'}
interac_1[3] = {'interaction_id':'W33-LEU17:N','res_id':17,'res_name':'LEU','res_atomtype':'N'}
interac_1[4] = {'interaction_id':'W33-TRP47','res_id':47,'res_name':'TRP','res_atomtype':'O'}
interac_1[5] = {'interaction_id':'W33-VAL48','res_id':48,'res_name':'VAL','res_atomtype':'N'}

### WATER 2
interac_2[1] = {'interaction_id':'W2-ASN31:ND2','res_id':31,'res_name':'ASN','res_atomtype':'ND2'}
interac_2[2] = {'interaction_id':'W2-PRO2:O','res_id':61,'res_name':'PRO','res_atomtype':'O'}
interac_2[3] = {'interaction_id':'W2-TYR4:O','res_id':63,'res_name':'TYR','res_atomtype':'O'}
interac_2[4] = {'interaction_id':'W2-TRP36:NE1','res_id':36,'res_name':'TRP','res_atomtype':'NE1'}

### WATER 2 EXT
interac_4[1] = {'interaction_id':'W2-ASN31:ND2','res_id':31,'res_name':'ASN','res_atomtype':'ND2'}
interac_4[2] = {'interaction_id':'W2-PRO2:O','res_id':61,'res_name':'PRO','res_atomtype':'O'}
interac_4[3] = {'interaction_id':'W2-TYR4:O','res_id':63,'res_name':'TYR','res_atomtype':'O'}
interac_4[4] = {'interaction_id':'W2-TRP36:NE1','res_id':36,'res_name':'TRP','res_atomtype':'NE1'}
interac_4[5] = {'interaction_id':'W2-SER3:=','res_id':62,'res_name':'SER','res_atomtype':'O'}


### WATER548  (index 548)
interac_3[1] = {'interaction_id':'W548-ser3:OG','res_id':62,'res_name':'SER','res_atomtype':'OG'}
interac_3[2] = {'interaction_id':'W548-ASP14:O','res_id':14,'res_name':'ASP','res_atomtype':'OD1'}
interac_3[3] = {'interaction_id':'W548-ALA:O','res_id':60,'res_name':'ALA','res_atomtype':'O'}



### water 2041 (similar to WAT35 in FYN)

interac_6[1] = {'interaction_id':'w2041-tyr63','res_id':63,'res_name':'TYR','res_atomtype':'OH'}
interac_6[2] = {'interaction_id':'w2041-ser64','res_id':64,'res_name':'SER','res_atomtype':'O'}
interac_6[3] = {'interaction_id':'w2041-ser12','res_id':12,'res_name':'SER','res_atomtype':'OG'}
interac_6[4] = {'interaction_id':'w2041-val10','res_id':10,'res_name':'VAL','res_atomtype':'O'}

### water 2077 (Wetspot expuesta, control)

interac_13[1] = {'interaction_id':'s1','res_id':8,'res_name':'ASP','res_atomtype':'O'}
interac_13[2] = {'interaction_id':'s4','res_id':8,'res_name':'ASP','res_atomtype':'N'}
interac_13[3] = {'interaction_id':'s5','res_id':52,'res_name':'TYR','res_atomtype':'OH'}
interac_13[4] = {'interaction_id':'s6','res_id':67,'res_name':'PRO','res_atomtype':'N'}


### WATER 3 

interac_12[1] = {'interaction_id':'W3-GLU98:OE1','res_id':35,'res_name':'GLU','res_atomtype':'OE1'}
interac_12[2] = {'interaction_id':'W3-ASN96:ND2','res_id':33,'res_name':'ASN','res_atomtype':'OD1'}
interac_12[3] = {'interaction_id':'W3-TYR4:O','res_id':63,'res_name':'TYR','res_atomtype':'O'}

### WATER 7 

#interac_4[1] = {'interaction_id':'W7-SER98:O','res_id':18,'res_name':'SER','res_atomtype':'O'}
#interac_4[2] = {'interaction_id':'W7-ASN96:O','res_id':43,'res_name':'ASN','res_atomtype':'OD1'}


### WATER 6 

#interac_7[1] = {'interaction_id':'W6-THR16:O','res_id':16,'res_name':'THR','res_atomtype':'O'}
#interac_7[2] = {'interaction_id':'W6-THR16:OG1','res_id':43,'res_name':'THR','res_atomtype':'OG1'}
#interac_7[3] = {'interaction_id':'W6-GLY13:N','res_id':13,'res_name':'GLY','res_atomtype':'N'}
#interac_7[4] = {'interaction_id':'W6-ALA11:O','res_id':11,'res_name':'ALA','res_atomtype':'O'}




#wetspots = {'W4':interaction, 'W2':interac_2, 'W5':interac_5, 'W1':interac_11, 'W3':interac_12, }
#wetspots = {'W5':interac_5 }
wetspots = {'W33':interac_1, 'W2077':interac_13, 'W2041':interac_6, 'W5':interac_5, 'W1':interac_11,  'W3':interac_12,'WAT548':interac_3, 'W4EXT':interac_7, 'W2EXT':interac_4 }
wetspots_cutoff_mod = { 'W5':0, 'W4':0, 'W2':0, 'W1':0, 'W3':0,'W4-FYN':0,'W33':0, 'W2077':0, 'W2041':0,'WAT548':0,'W7':0,'W6':0,'W4EXT':0,'W2EXT':0 }

coordw5 = [-3.730,-1.558,8.467]
coordw2 = [3.260,-8.789,8.085]
coordw3 = [1.670,-9.253,10.297]
coordw4 = [-3.296,-3.891,5.697]
coordw1 = [5.005,-10.249,6.030]

user_centroids = {'W1':coordw1,'W2':coordw2,'W3':coordw3,'W4':coordw4,'W5':coordw5,'W4-FYN':coordw4}

INTERACTION_SPOT_INDEX = []
for wetspot_key,wetspot in wetspots.iteritems():
            for interaction in wetspot.itervalues():
                INTERACTION_SPOT_INDEX.append(interaction)
