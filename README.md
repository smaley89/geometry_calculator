# geometry_calculator
A python code for calculating bond lengths, angles and dihedrals from an XYZ file

## Usage
Call geometry_calculator.py on an xyz file, specifying the atom indices of interest.
2 atoms will calculate the distance (angstroms)
3 atoms will calculate the angle (degrees)
4 atoms will calculate the dihedral angle (degrees)

Example: python3 geometry_calculator.py -f your_xyz_file.xyz -g 1 2
This will calculate the distance between atom 1 and atom 2

The percision of the result can be specified using the optional '-p' flag. The default percision is 3 decimal places.

Example: python3 geometry_calculator.py -f your_xyz_file.xyz -g 1 2 3 4 -p 5
This will calculate the dihedral between atom 1, atom 2, atom 3, and atom 4 to 5 decimal places.
