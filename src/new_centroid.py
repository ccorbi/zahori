import math

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


def is_inclination(v1,v2):
 	if angle(v1,v2) <= math.pi:
 		return True
 	else:
 		return False

def is_azimuth(v1,v2):

 	if angle(v1,v2) <= 2*math.pi:
 		return True
 	else:
 		return False

def euclidian_distance(x,y,z,a,b,c):

        '''Calc euclidian distance in 3D space, and return float'''

        return math.sqrt((x-a)**2+(y-b)**2+(z-c)**2)

def conv_coord(radius,incltn,azith):


 	x = radius * math.sin(incltn) * math.cos(azith)
 	y = radius * math.sin(incltn) * math.sin(azith)
 	z = radius * math.cos(incltn)
 	coor = []
 	coor.append(x)
 	coor.append(y)
 	coor.append(z)

 	return coor

def trans_coor(origin,v2):
	coor = []
	for i in range(3):
		coor.append(v2[i] - origin[i])

	return coor

def retrans_coor(origin,v2):

	coor = []
	for i in range(3):
		coor.append(v2[i] + origin[i] )

	return coor


def radius(v1):
	value = 0.0
	for i in range(3):
		value += v1[i] **2

	return math.sqrt(value)

def get_shperical_coor(v1):

	shperical = []
	r = radius(v1)
	shperical.append(r)
	shperical.append(math.acos(v1[2]/r))
	
	shperical.append(math.atan(v1[1]/v1[0]))

	return shperical


def new_shperical(v1,r):

	shperical = []

	shperical.append(r)
	shperical.append(math.acos(v1[2]/r))
	
	shperical.append(math.atan(v1[1]/v1[0]))

	return shperical


def trans_spherical_ref(coord, radius,incltn,azith):
	a = conv_coord(radius,incltn,azith)
	return rebuild_ref(coord,a)


 #OG1 THR    90 
refpoint = [4.897, -3.240, -2.669]

 #CG2 THR    90       3.630  -1.191  -2.748
refpoint2 = [3.630,  -1.191 , -2.748]
water3 = [2.704, -3.934, -0.363]
refpoint3 = [4.630,  -1.191 , -2.748]

new_origin = trans_coor(refpoint,water3)
print new_origin
c = get_shperical_coor(new_origin)
print c
a = conv_coord(c[0],c[1],c[2])
print a
final = retrans_coor(refpoint,a)
print final

transformation = trans_coor(refpoint2,water3)
print 'trans', transformation
c = get_shperical_coor(transformation)
print 'geo_data', c
a = conv_coord(c[0],c[1],c[2])
print 'new_vector', a
final = retrans_coor(refpoint2,a)
print 'convert',final

d = euclidian_distance(refpoint2[0],refpoint2[1],refpoint2[2],water3[0],water3[1],water3[2])
vector = [water3[0]-refpoint2[0],water3[1]-refpoint2[1],water3[2]-refpoint2[2]]
she_vec = new_shperical(vector,d)
print 'distance ',d,' vector ',vector,'shperic values ', she_vec
print ' conversion to coordinates ', conv_coord(she_vec[0],she_vec[1],she_vec[2])
print retrans_coor(refpoint2,she_vec)





