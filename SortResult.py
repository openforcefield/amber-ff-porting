import sys, os

def print_group(D):
  print('  Type     Amber      OFF      Diff.   Significant?')
  for key in D:
    print(key)
    flag = '      '
    if (abs(D[key][0] - D[key][4]) > 5.0e-3):
      flag = '   ***'
    print('  Bond:  %9.4f %9.4f %9.4f %s' % (D[key][0], D[key][4], D[key][0] - D[key][4], flag))
    flag = '   '
    if (abs(D[key][1] - D[key][5]) > 5.0e-3):
      flag = '   ***'
    print('  Angl:  %9.4f %9.4f %9.4f %s' % (D[key][1], D[key][5], D[key][1] - D[key][5], flag))
    flag = '   '
    if (abs(D[key][2] - D[key][6]) > 5.0e-3):
      flag = '   ***'
    print('  Dihe:  %9.4f %9.4f %9.4f %s' % (D[key][2], D[key][6], D[key][2] - D[key][6], flag))
    flag = '   '
    if (abs(D[key][3] - D[key][7]) > 5.0e-3):
      flag = '   ***'
    print('  Nonb:  %9.4f %9.4f %9.4f %s' % (D[key][3], D[key][7], D[key][3] - D[key][7], flag))

fi = open(sys.argv[1], 'r')
fmem = [ line for line in fi ]
print(fmem)
mainch = {}
nterm  = {}
cterm  = {}
for i, line in enumerate(fmem):
  path = os.path.normpath(line.split('\n')[0])
  tags = path.split(os.sep)
  if (len(tags) < 3):
    continue
  residue = tags[1]
  print(line, tags)
  ambBonds = 0.0
  ambAngls = 0.0
  ambDihes = 0.0
  ambNonbs = 0.0
  for j in range(1, 6):
    ttl = fmem[i+j].split()
    print(ttl)
    if (ttl[0] == 'HarmonicBondForce'):
      ambBonds = float(ttl[1])
    elif (ttl[0] == 'HarmonicAngleForce'):
      ambAngls = float(ttl[1])
    elif (ttl[0] == 'PeriodicTorsionForce'):
      ambDihes = float(ttl[1])
    elif (ttl[0] == 'NonbondedForce'):
      ambNonbs = float(ttl[1])
  offBonds = 0.0
  offAngls = 0.0
  offDihes = 0.0
  offNonbs = 0.0
  for j in range(6, 11):
    ttl = fmem[i+j].split()
    print(ttl)
    if (ttl[0] == 'HarmonicBondForce'):
      offBonds = float(ttl[1])
    elif (ttl[0] == 'HarmonicAngleForce'):
      offAngls = float(ttl[1])
    elif (ttl[0] == 'PeriodicTorsionForce'):
      offDihes = float(ttl[1])
    elif (ttl[0] == 'NonbondedForce'):
      offNonbs = float(ttl[1])
  if (tags[0] == 'MainChain'):
    mainch[tags[1]] = (ambBonds, ambAngls, ambDihes, ambNonbs,
                       offBonds, offAngls, offDihes, offNonbs)
  elif (tags[0] == 'NTerminal'):
    nterm[tags[1]] = (ambBonds, ambAngls, ambDihes, ambNonbs,
                      offBonds, offAngls, offDihes, offNonbs)
  elif (tags[0] == 'CTerminal'):
    cterm[tags[1]] = (ambBonds, ambAngls, ambDihes, ambNonbs,
                      offBonds, offAngls, offDihes, offNonbs)

print('Main Chain')
print_group(mainch)
print('\nN-Terminal')
print_group(nterm)
print('\nC-Terminal')
print_group(cterm)
