f = open('queries.txt','w+')

numStates = 5;
numTimesteps = 5;

# single states
for t in xrange(1,numTimesteps+1):
  for i in xrange(1,numStates+1):
    f.write(str(t)+"="+str(i)+","+"\n")
# transitions
for t in xrange(2,numTimesteps+1):
  for i in xrange(1,numStates+1):
    for j in xrange(1,numStates+1):
      f.write(str(t-1)+"="+str(i)+","+str(t)+"="+str(j)+"\n")
