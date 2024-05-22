def filter_cartopy_warnings():
    global APPLIED_FILTER_WARNINGS

    if not APPLIED_FILTER_WARNINGS:
        warnings.filterwarnings("ignore", message="Approximating coordinate system")
        APPLIED_FILTER_WARNINGS = True
        
def sECEFtoECIpatch(epc, x):


    """Transforms an Earth fixed state into an Inertial state

    The transformation is accomplished using the IAU 2006/2000A, CIO-based
    theory using classical angles. The method as described in section 5.5 of
    the SOFA C transformation cookbook.

    Args:
        epc (Epoch): Epoch of transformation
        x (np.ndarray): Earth-fixed state (position, velocity) [*m*; *m/s*]

    Returns:
        x_ecef (np.ndarray): Inertial state (position, velocity)
    """

    # Ensure input is array-like
    x = np.asarray(x)

    # Set state variable size
    dim_x = len(x)
    x_eci = np.zeros((dim_x,))

    # Extract State Components
    r_ecef = x[0:3]

    if dim_x >= 6:
        v_ecef = x[3:6]

    # Compute Sequential Transformation Matrices
    rot = bh.earth_rotation(epc)

    # Create Earth's Angular Rotation Vector
    omega_vec = np.array([0, 0, bh.constants.OMEGA_EARTH]) # Neglect LOD effect

    # Calculate ECEF State
    x_eci[0:3] = ( rot ).T @ r_ecef
    # x_eci[0:3] = (pm @ rot @ bpn).T @ r_ecef

    if dim_x >= 6:
        x_eci[3:6] = (rot ).T @ (v_ecef + bh.utils.fcross(omega_vec, r_ecef))

    return x_eci

def get_tle(epc0, alt, ecc, inc, raan, argp, M, ndt2=0.0, nddt6=0.0, bstar=0.0, norad_id=99999):
    '''Get a TLE object from the given orbital elements

    Args:
    - epc0 (Epoch): Epoch of the orbital elements / state
    - alt (float): Altitude of the orbit [km]
    - ecc (float): Eccentricity of the orbit
    - inc (float): Inclination of the orbit [deg]
    - raan (float): Right Ascension of the Ascending Node [deg]
    - argp (float): Argument of Perigee [deg]
    - M (float): Mean Anomaly [deg]

    Returns:
    - tle (TLE): TLE object for the given orbital elements
    '''

    alt *= 1e3 # Convert to meters

    # Get semi-major axis
    sma = bh.R_EARTH + alt

    # Get mean motion
    n = bh.mean_motion(sma)/(2*np.pi)*86400

    tle_string = bh.tle_string_from_elements(epc0, np.array([n, ecc, inc, raan, argp, M, ndt2, nddt6, bstar]), norad_id)
    tle = bh.TLE(*tle_string)
    return tle

def compute_earth_interior_angle(ele=0.0, alt=525):
    '''This function computes the earth interior angle for a given elevation angle and altitude.
    This is the angle between the satellite and the observer on the Earth's surface viewing the satellite
    at the given elevation angle. This is more useful for plotting than the look angle, since it
    can used alonside Earth's Radius to draw a circle around the subsatellite point to get
    the view cone of observers on the Earth's surface that would be able to see the satellite.


    Args:
    - ele (float): Elevation angle of the satellite [deg]
    - alt (float): Altitude of the satellite [km]

    '''
    ele = ele * math.pi / 180.0

    rho = math.asin(bh.R_EARTH/(bh.R_EARTH + alt * 1e3))

    eta = math.asin(math.cos(ele)*math.sin(rho))
    lam = math.pi/2.0 - eta - ele

    return lam

def proparea(height_incl):
    num = len(height_incl) / 2
    epc0 = bh.Epoch(2024, 5, 20, 0, 0, 0)
    ecc = 0.001
    raan = 45
    argp = 90
    M = 0
    norad_id = 99999
    
    tle_list_initial = []
    for i in num:
        tle_list_initial.append(get_tle(epc0, height_incl(i), ecc, height_incl(num + i), raan, argp, M, norad_id=norad_id))

    time_step = 1
    time = 1
    while time < 11:
        t = epc0 + time

        for i in num:
            x_eci = tle_list_initial[i].state_eci(t)
            x_ecef = tle_list_initial[i].state_ecef(t)
            x_geod = bh.sECEFtoGEOD(x_ecef[0:3], use_degrees=True)
            alt = x_geod[2]

    








