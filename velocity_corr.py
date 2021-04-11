import math
import matplotlib.pyplot as plt
import sys
f = open(sys.argv[1])
lines = f.readlines()
frames = []
O = []
vs = []
for i in range(1, len(lines)):
    tokens = lines[i].split()    
    if tokens[0] == 'END':
        frames.append(O)
        O = []
        continue
    x = float(tokens[5])
    y = float(tokens[6])
    z = float(tokens[7])
    O.append([x, y, z, tokens[1]])

vs = frames

n = len(frames) # Number of time steps
N = len(frames[0]) # Number of atoms.

print(n)
print(N)

res = []

for i in range(n - 1):
    s = 0
    cnt = 0
    for j in range(n - 1 - i):
        tmp = 0
        for k in range(N):
            tmp = tmp + (vs[j][k][0] * vs[j + i][k][0] + vs[j][k][1] * vs[j + i][k][1] + vs[j][k][2] * vs[j + i][k][2])
        tmp = tmp / N
        s = s + tmp
        cnt = cnt + 1
    if cnt == 0:
        break
    s = s / cnt
    res.append(s)
plt.title('Velocity correlation vs Time interval plot')
plt.xlabel('Time interval(t)')
plt.ylabel('Velocity correlation')
plt.plot(res)
plt.show()