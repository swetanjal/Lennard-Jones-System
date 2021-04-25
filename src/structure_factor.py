import math
import matplotlib.pyplot as plt
import sys
import numpy as np
T = 2000
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

def dist(a, b):
    return math.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))

rs = []
for i in range(35):
    k = i
    for j in  range(10):
        k = k + 0.1
        rs.append(k)

def compute(t):
    ys = []
    for r in rs:
        res = 0.0
        freq = 0
        for o in range(0, 1):
            tau = o + t
            if tau >= n:
                break
            hit = 0
            for i in range(N):
                for j in range(N):
                    if abs(dist(frames[o][i],frames[tau][j]) - r) <= 0.01:
                        hit = hit + 1
            res = res + hit
            freq = freq + 1
        res = res * 1.0 / freq
        res = res * 1.0 / (N*N*N)
        
        ys.append(res)
    return ys

As = []
for i in range(0, T, 80):
    print(i)
    ys = compute(i)
    As.append(np.array(ys))

As = np.array(As)
res = np.fft.fft2(As)

print(res.shape)

for i in range(As.shape[0]):
    xs = []
    ys = []
    for j in range(len(rs)):
        ys.append(np.absolute(res[i, j]))
        xs.append(rs[j])
    plt.plot(xs, ys, label = "t = " + str(i))
# plt.legend(loc = 'upper right')    
plt.show()