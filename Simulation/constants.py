"""Inventory of numerical constants, their units, and what they are used for."""


G = 6.6743*10**(-11)        # Grav Constant
earthMass = 5.9722 * 10**(24) # mass of Earth
moonMass = .07346 * 10 ** 24
mustar = moonMass/(moonMass + earthMass)
lstar = 3.844 * 10 ** 8
tstar = 382981


EMmu = G*(earthMass+moonMass)
