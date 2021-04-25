import matplotlib.pyplot as plt
import sys
f = open(sys.argv[1])
lines = f.readlines()

xs = []
ys = []

for line in lines[3:]:
    toks = line.split()
    xs.append(float(toks[0]))
    ys.append(float(toks[1]))

plt.title('Dynamic Structure factor vs k')
plt.ylabel('S(k)')
plt.xlabel('k values')
plt.plot(xs, ys)
plt.show()