# Convert a radial velocity to the Galactic Standard of Rest(GSR)
# 방사형 속도를 GSR(은하 표준 정지 좌표계)
# from : https://docs.astropy.org/en/stable/generated/examples/coordinates/rv-to-gsr.html#sphx-glr-generated-examples-coordinates-rv-to-gsr-py
import astropy.units as u
import astropy.coordinates as coord

# defaults can be changed at runtime by setting the parameter
# Sun-Galactic center distance value can be different because this default value is from 2014
coord.galactocentric_frame_defaults.set('latest')

# star - HD 155967
icrs = coord.SkyCoord(ra = 258.58356362*u.deg, dec = 14.55255619*u.deg, 
                      radial_velocity = -16.1*u.km/u.s, frame = 'icrs')

# coordinates.Galactocentric's document : https://docs.astropy.org/en/stable/api/astropy.coordinates.Galactocentric.html#astropy.coordinates.Galactocentric
# convert it to a 'CartesianRepresentation' object using the .to_cartesian() method of the 'CartesianDifferential object galcen_v_sun
# to_cartesian() : convert the representation to its Cartesian form
v_sun = coord.Galactocentric().galcen_v_sun.to_cartesian()

# making the unit vector
# transform_to(frame, merge_attributes = Ture) : transform this coordinate to a new frame
gal = icrs.transform_to(coord.Galactic)
# to_cartesian() : convert spherical polar coordinates to 3-dimension rectangular cartesian coordinates
cart_data = gal.data.to_cartesian()
unit_vector = cart_data / cart_data.norm()

v_proj = v_sun.dot(unit_vector)

rv_gsr = icrs.radial_velocity + v_proj
print(rv_gsr)

def rv_to_gsr(c, v_sun = None):
    # c : The radial velocity which is associated with a sky coordinates to be transformed
    # v_sun : optional parameter. 3-dimension velocity of the solar system barycenter in the GSR frame. Defaults is comparable with solar motion in the coordinates.Galactocentric frame.
    if v_sun is None:
        v_sun = coord.Galactocentric().galcen_v_sun.to_cartesian()
        
    gal = c.transform_to(coord.Galactic)
    cart_data = gal.data.to_cartesian()
    unit_vector = cart_data / cart_data.norm()
    
    v_proj = v_sun.dot(unit_vector)
    
    return c.radial_velocity + v_proj

rv_gsr = rv_to_gsr(icrs)
print(rv_gsr)

'''
mathematics behind the transformation from ICRS to Galactocentric coordinates
-> https://docs.astropy.org/en/stable/coordinates/galactocentric.html#coordinates-galactocentric
'''