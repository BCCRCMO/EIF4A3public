#!/Users/alborzmazloomian/anaconda/bin/python

psi = {}

with open("all_events", 'r') as f:
    content = f.readlines()
    for i in range(len(content)):
        cols = content[i].rstrip().split(' ')
        psi[(cols[0], cols[1], cols[2], int(cols[5]))] = abs(float(cols[3]))


with open("diff_events_per_drug", 'w') as f:
    for key, val in psi.items():
        drugName = key[0][-5:]
        if key[3] == 1:
            f.write("{} {} {} {}\n".format(drugName, key[1], key[3], val))
        elif (key[0], key[1], key[2], key[3]-1) in psi:
            f.write("{} {} {} {}\n".format(drugName, key[1], key[3]-1, val-psi[(key[0], key[1], key[2], key[3]-1)]))
