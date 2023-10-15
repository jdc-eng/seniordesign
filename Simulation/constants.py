"""Inventory of numerical constants, their units, and what they are used for."""


G = 6.6743*10**(-11)        # Grav Constant

earthMass = 5.9722 * 10 ** 24 # mass of Earth
earthD = 12742 * 10 ** 3 
moonMass = .07346 * 10 ** 24
moonD = 3474.8 * 10 ** 3
mustar = moonMass/(moonMass + earthMass)
lstar = 3.844 * 10 ** 8
tstar = 382981


EMmu = G*(earthMass+moonMass)
