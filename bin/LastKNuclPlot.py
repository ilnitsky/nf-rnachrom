from sys import argv
import matplotlib.pyplot as plt

#alphavet desctipton
alpha = ['A', 'G', 'C', 'T', 'N']
n = argv[2]
fl = argv[3]

#generating a dictionary
d = {}
if n == '3':
    for i in alpha:
        for j in alpha:
            for k in alpha:
                d[i+j+k] = 0
elif n == '2':
    for i in alpha:
        for j in alpha:
            d[i+j] = 0

#for k in d.keys():
#       print(f'key {k} value {d[k]}')

#defining counting function
def step(s):
    if len(s) != int(n):
        return
    else:
        d[s] += 1

#counting
f = open(argv[1], 'r')
m = []
s = ""

if fl == "l":
    for line in f:
        if line[0] == '@':
            if (len(m) == 3 and m[2][0] != '+') or (len(m) == 4):
                s = m[1]
                s = s[-(int(n)+1):-1]
                step(s)
                m = []
        m.append(line)

if fl == "f":
    for line in f:
        if line[0] == '@':
            if (len(m) == 3 and m[2][0] != '+') or (len(m) == 4):
                s = m[1]
                if len(s) >= int(n):
                    s = s[:(int(n))]
                    if s in d.keys():
                        step(s)
            m = []
        m.append(line)

step(m)

#converting to precentage
c = 0
for key in d:
    c += d[key]
for key in d:
    d[key] = 100 * d[key] / c

#generating image
courses = list(d.keys())
values = list(d.values())
fig = plt.figure(figsize=(20,10))
plt.bar(courses, values, color='maroon', width=0.4)

#labeling and saving file
if fl == "l":
    plt.xlabel(f'Last {n} nucleotides')
    plt.xticks(rotation = 90)
    plt.ylabel('Percentage')
    plt.title(f'Last {n} nucleotides distribution of {argv[1]}')
    plt.savefig(f'{argv[1]}_last{n}.png')
elif fl == "f":
    plt.xlabel(f'First {n} nucleotides')
    plt.xticks(rotation = 90)
    plt.ylabel('Percentage')
    plt.title(f'First {n} nucleotides distribution of {argv[1]}')
    plt.savefig(f'{argv[1]}_first{n}.png')