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
        if freq == 0:
            break
        res = res * 1.0 / freq
        res = res * 1.0 / (N*N*N)
        ys.append(res)
    return ys


plt.title('Van Hove Correlation Function')
plt.xlabel('r')
plt.ylabel('Van Hove Correlation value')
# ys = compute(0)
# plt.plot(ys, label = "t = 0")
ys = compute(10)
plt.plot(ys, label = "t = 10")
ys = compute(25)
plt.plot(ys, label = "t = 25")
ys = compute(50)
plt.plot(ys, label = "t = 50")
ys = compute(100)
plt.plot(ys, label = "t = 100")
ys = compute(150)
plt.plot(ys, label = "t = 150")
ys = compute(200)
plt.plot(ys, label = "t = 200")
ys = compute(250)
plt.plot(ys, label = "t = 250")
plt.legend(loc = 'upper right')
plt.show()
# plt.savefig('vann_hove.png')
