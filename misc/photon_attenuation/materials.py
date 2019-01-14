import constants

class Element:

    name = ""
    density = 0
    atomic_number = 0
    molar_mass = 0
    electron_density = 0

    def __init__(self, name, density, atomic_number, molar_mass):

        self.name = name
        self.density = density
        self.atomic_number = atomic_number
        self.molar_mass = molar_mass

        atoms_in_kg = constants.N_a / molar_mass
        self.electron_density = atoms_in_kg * atomic_number * density

class Compound:

    name = ""
    elements = []
    density = 0
    molar_mass = 0
    electron_density = 0
    atomic_number = 0

    def __init__(self, name, elements, density, molar_mass):

        self.name = name
        self.density = density
        self.molar_mass = molar_mass
        self.elements = []

        self.atomic_number = sum([x.atomic_number for x in elements])

        atoms_in_kg = constants.N_a / molar_mass
        self.electron_density = atoms_in_kg * self.atomic_number * density

Si = Element("Si", 2320, 14, 0.02808550)
Cd = Element("Cd", 8650, 48, 0.112411)
Na = Element("Cd", 968, 11, 0.022989770)
I = Element("I", 4933, 53, 0.253808940)
Te = Element("Te", 5850, 52, 0.12760)
CdTe = Compound("CdTe", [Cd, Te], 5850, 0.2400110)
NaI = Compound("NaI", [Na, I], 3670, 0.149894239)
