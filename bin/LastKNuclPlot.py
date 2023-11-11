from sys import argv
import matplotlib.pyplot as plt

#alphavet desctipton
alpha = ['A', 'G', 'C', 'T', 'N']
n = argv[2]

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

#defining counting function
def step(m):
    if len(m) == 0:
        return
    else:
        d[m[1][-(int(n)+1):-1]] += 1

#counting
f = open(argv[1], 'r')
m = []

for line in f:
    if line[0] == '@':
        if (len(m) == 3 and m[2][0] != '+') or (len(m) == 4):
            step(m)
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
plt.xlabel(f'Last {n} nucleotides')
plt.xticks(rotation = 90)
plt.ylabel('Percentage')
plt.title(f'Last {n} nucleotides distribustion of {argv[1]}')
plt.savefig(f'{argv[1]}{n}.png')
