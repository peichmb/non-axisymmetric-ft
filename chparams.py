from sys import argv

def replaceLine(fname, lineNum, newLine):

    if newLine[:2] != '\n':
        newLine = newLine+'\n'
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    
    lines[lineNum-1] = newLine

    f = open(fname, 'w')
    f.writelines(lines)
    f.close()



vv, dd, mm, tt = argv[1:]

vvNum = 42
vvStr = 'dRotModel = ' + vv + ' # Differential rotation model (see flows.py)'
replaceLine('params.py',vvNum,vvStr)

ddNum = 32
ddStr = 'eta0 = ' + dd + 'e10 # Diffusion coefficient [cm**2/s]'
replaceLine('params.py',ddNum,ddStr)

mmNum = 33
mmStr = 'mflow_return = ' + mm + ' # Magnitude of return meridional flow [cm/s]'
replaceLine('params.py',mmNum,mmStr)

if tt == 'xxxx':
    replaceLine('params.py',35,'limitMagneticMemory = False')
else:
    replaceLine('params.py',35,'limitMagneticMemory = True')
    ttNum = 36
    ttStr = 'decayTime = ' + tt + ' # Magnetic memory of the system'
    replaceLine('params.py',ttNum,ttStr)





