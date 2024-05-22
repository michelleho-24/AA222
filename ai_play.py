from datetime import datetime, timedelta
from typing import List, Tuple


class Epoch:
    def __init__(self, datetime_obj):
        self.datetime_obj = datetime_obj

    def __repr__(self):
        return f"Epoch({self.datetime_obj.isoformat()})"

def generate_initial_conditions(epc0: Epoch, alt: float, ecc: float, inc: float, raan: float, argp: float, M: float) -> List[Tuple[Epoch, float, float, float, float, float, float]]:
    """
    Generate initial conditions for 5 satellites that are evenly spaced.
    
    Args:
    - epc0 (Epoch): Epoch of the orbital elements / state
    - alt (float): Altitude of the orbit [km]
    - ecc (float): Eccentricity of the orbit
    - inc (float): Inclination of the orbit [deg]
    - raan (float): Right Ascension of the Ascending Node [deg]
    - argp (float): Argument of Perigee [deg]
    - M (float): Mean Anomaly [deg]
    
    Returns:
    - List of tuples containing initial conditions for each satellite
    """
    initial_conditions = []
    number_of_satellites = 5
    increment = 360 / number_of_satellites
    
    for i in range(number_of_satellites):
        inclination = inc + i * increment
        # Ensure mean anomaly is within 0-360 degrees
        inclination = inclination % 360
        initial_conditions.append((epc0, alt, ecc, inclination, raan, argp, M))
    
    return initial_conditions

# Example usage
epc0 = Epoch(datetime.utcnow())
alt = 500  # Altitude in km
ecc = 0.001  # Eccentricity
inc = 98.0  # Inclination in degrees
raan = 45.0  # Right Ascension of the Ascending Node in degrees
argp = 90.0  # Argument of Perigee in degrees
M = 0.0  # Mean Anomaly in degrees

initial_conditions = generate_initial_conditions(epc0, alt, ecc, inc, raan, argp, M)
for i, conditions in enumerate(initial_conditions):
    print(f"Satellite {i+1}: {conditions}")


