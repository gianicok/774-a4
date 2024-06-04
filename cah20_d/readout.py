import matplotlib.pyplot as plt
import numpy as np

interest = 0
file_path = 'egs5job.out'
with open(file_path, 'r') as file:
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