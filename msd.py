import math
import matplotlib.pyplot as plt
import sys
if len(sys.argv) != 2:
    print('Incorrect usage. Refer README.')
    exit()

f = open(sys.argv[1])
lines = f.readlines()
frames = []
O = []

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

n = len(frames)
N = len(frames[0])

disp = 0.0
res = []
for i in range(1, n):
    cnt = 0
    s = 0
    for j in range(n - i):
        tmp = 0
        for k in range(N):
            assert(frames[i][k][3] == frames[j][k][3])
            tmp = tmp + (((frames[j + i][k][0] - frames[j][k][0])**2) + ((frames[j + i][k][1] - frames[j][k][1])**2) + ((frames[j + i][k][2] - frames[j][k][2])**2))
        tmp = tmp / N
        s = s + tmp
        cnt = cnt + 1
    if cnt == 0:
        break
    s = s / cnt
    res.append(s)

plt.title('Mean Square Displacement vs Time interval plot')
plt.xlabel('Time interval(t)')
plt.ylabel('Mean Square Displacement')
plt.plot(res)
plt.show()

print("Diffusion Coefficient = ", (res[1900] - res[880]) / (1900-880))