import matplotlib.pyplot as plt
import numpy as np

folder = './cah20_b/'

interest = 0
with open(folder+'egs5job.out', 'r') as file:
    lines = file.readlines()
    for i in range(len(lines)):
        line = lines[i].strip()
        if line == 'ENERGY DEPOSITION SUMMARY FOR PARTICLES WITH IQ=-1':
            interest = i-1
            #print(line)

for i in range(interest,interest+len(lines[interest:])):
    line = lines[i].strip()
    if(line!=''):
        print(line)

u, v, w = [], [], []
with open(folder+'trajectories.out', 'r') as file:
    lines = file.readlines()
    for i in range(len(lines)):
        line = lines[i].strip().split()
        u.append(float(line[1]))
        v.append(float(line[2]))
        w.append(float(line[3]))

#print(u,v,w)
        
# plot results
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)
ax.scatter(u, v, w, marker='o')
#ax.quiver(0, 0, 0, u, v, w,linewidth=1)
ax.set_xlabel('U')
ax.set_ylabel('V')
ax.set_zlabel('W')
ax.set_title('EGS5 uin, vin, win')
plt.show()