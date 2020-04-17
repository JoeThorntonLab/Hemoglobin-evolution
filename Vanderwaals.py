def distance (a,b):
  distance = ((b[0]-a[0])**2+(b[1]-a[1])**2+(b[2]-a[2])**2)**0.5
  return distance
infile = open("B14.pdb")  ##Enter file name here
chain = []
residues = []
coordinates = [] 
for line in infile:
  line1 = line.strip()
  if "ATOM" in line1 and not "REMARK" in line1:
    k = line.split()
  #  if k[2] == "N" or k[2] == "CA":
 #   print ((k[3]+k[5]),k[4],k[6]+k[7]+k[8])
    chain.append(k[4])
    coordinates.append((float(k[6]),float(k[7]),float(k[8])))
    residues.append((k[3]+k[5]))
    
c1 = "A"
c2 = "B"
pairs = []
a_residue_number = []
b_residue_number = []

for i in range (len(residues)):
  for j in range (len(residues)):
    if (chain[i] == c1 and chain[j] ==c2): #or (chain[i] == c2 and chain[j] ==c1):
      if distance(coordinates[i],coordinates[j])<3.5:
         pairs.append((chain[i]+residues[i],chain[j]+residues[j]))
         a_residue_number.append(residues[i]+chain[i])
         b_residue_number.append(residues[j]+chain[j])
         
       
                           
pairs = set(pairs)
a_residue_number = set(a_residue_number)
b_residue_number = set(b_residue_number)

#print (len(pairs))
print(pairs)
print (a_residue_number)
print (b_residue_number)
