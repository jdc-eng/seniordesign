import constants as c

def StateDict(arg):
    states = {
        'Pointy Star': [0.5, 0, 0, 0, -0.1, 0],
        'Example 1':   [1-c.mustar, .0455, 0, -0.5, 0.5, 0]
    }

    return states[arg]
