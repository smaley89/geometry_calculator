import argparse
import numpy as np


def calculate_distance(atom1: np.ndarray, atom2:np.ndarray, 
                       percision: int) -> float:
    
    '''
    Calculate the Euclidean distance of two atoms in Cartesian space

    Input: Numpy Array for atom 1 and atom 2

    Return: Distance between atoms in Angstroms
    '''

    squared_distance = np.sum((atom1-atom2)**2, axis=0)
    distance = np.sqrt(squared_distance)

    return round(distance, percision)


def calculate_angle(atom1: np.ndarray, atom2: np.ndarray, 
                    atom3: np.ndarray, percision: int) -> float:
    
    '''
    Calculate the angle between three atoms in Cartesian space

    Input: Numpy array for atom 1, atom 2, and atom 3

    Return: Angle between atoms in Degrees
    '''

    atom1_atom2 = atom1 - atom2
    atom3_atom2 = atom3 - atom2

    cosine_angle = np.dot(atom1_atom2, atom3_atom2) / (np.linalg.norm(atom1_atom2) * np.linalg.norm(atom3_atom2))
    angle = np.arccos(cosine_angle) *180.0/np.pi

    return round(angle, percision)


def calculate_dihedral(atom1: np.ndarray, atom2: np.ndarray, atom3: np.ndarray, 
                       atom4: np.ndarray, percision: int) -> float:
    '''
    Calculate the dihedral angle

    Input: Numpy array for atom 1, atom 2, atom 3, and atom 4

    Return: Dihedral angle of specified four atoms
    
    '''

    a = -1.0*(atom2 - atom1)
    b = atom3 - atom2
    c = atom4 - atom3

    b /= np.linalg.norm(b)
    v = a - np.dot(a, b)*b
    w = c - np.dot(c, b)*b

    x = np.dot(v, w)
    y = np.dot(np.cross(b, v), w)

    return round(np.degrees(np.arctan2(y, x)), percision)


def read_xyz(filename: str) -> np.array:
    '''
    Read the xyz file and return only the cartesian coordinates
    
    Input: Filename of the XYZ file

    Return: Cartesian coordinates as NumPy array
    '''
    with open(filename, 'r') as f:
        lines = [line.rstrip().split() for line in f]

    coordinates = [line[1:] for line in lines[2:]]

    return np.float_(coordinates)


def main(filename: str, geom_value: list, percision: int) -> float:
    
    geometry = read_xyz(filename)

    if len(geom_value) == 2:
        atom1, atom2 = geometry[geom_value[0:2]]
        distance = calculate_distance(atom1, atom2, percision)
        print(f'{geom_value[0]+1}-{geom_value[1]+1} distance is {distance}')

    elif len(geom_value) == 3:
        atom1, atom2, atom3 = geometry[geom_value[0:3]]
        angle = calculate_angle(atom1, atom2, atom3, percision)
        print(f'{geom_value[0]+1}-{geom_value[1]+1}-{geom_value[2]+1} angle is {angle}')

    elif len(geom_value) == 4:
        atom1, atom2, atom3, atom4 = geometry[geom_value[0:4]]
        dihedral = calculate_dihedral(atom1, atom2, atom3, atom4, percision)
        print(f'{geom_value[0]+1}-{geom_value[1]+1}-{geom_value[2]+1}-{geom_value[3]+1} angle is {dihedral}')

    else:
        print('Not recognized')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, required=True, help='Enter name of XYZ file')
    parser.add_argument('-g', type=int, nargs='+', help='Atom indices for geometry parameter of interest')
    parser.add_argument('-p', type=int, nargs='?', required=False, default=3, help='Set percision of value (default = 3)')
    args = parser.parse_args()
    geom_value = args.g
    geom_value = [i-1 for i in geom_value] # Subtract one from inputs to accomodate zero indexing
    filename = args.f
    percision = args.p
    main(filename, geom_value, percision)