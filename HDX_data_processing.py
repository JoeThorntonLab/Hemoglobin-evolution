import os
import random
import math
import numpy as np
import pylab
import copy
from itertools import combinations as comb
from scipy import stats

def squared_error(ys_orig,ys_line):
    a = np.array(ys_orig)
    b = np.array(ys_line)
    return sum((b - a) * (b - a))

def coefficient_of_determination(ys_orig,ys_line):
    y_mean_line = [np.mean(ys_orig) for y in ys_orig]
    squared_error_regr = squared_error(ys_orig, ys_line)
    squared_error_y_mean = squared_error(ys_orig, y_mean_line)
    return 1 - (squared_error_regr/squared_error_y_mean)

def burial_value(peptide,t):
    infile = open(t,"r")
    d = {}
    i1_sum = 0
    i2_sum = 0
    for line in infile:
    #  print (line)
      k = line.split()
      #print (k)
      if float(k[3]) == 0 and float(k[4]) == 0:
        d[int(k[2])] = 0
      else:
        d[int(k[2])] = float(k[4])

      
    infile.close()
    peptide1=peptide.split()
    positions = peptide1[0].split(",")
  #  print (d)
    length = int(float(positions[1]))+1-int(float(positions[0]))
    for i in range(int(float(positions[0])),int(float(positions[1]))+1):
       if i in d:
        i2_sum += float(d[i])
    i2_average = i2_sum#/float(len(positions))
    infile.close()
    
    return float(i2_sum)/float(length)


def overlap(a,b):
  if (a[0]<b[0] and a[1]<b[0]) or (a[0]>b[1] and a[1]>b[1]):
    return False
  else:
    return True
def whatrep(x):
  if "_1_" in x:
    return " 1"
  if "_2_" in x:
    return " 2"
  if "_3_" in x:  
    return " 3"

def whattime(x):
  times = ["1min","0.25min","30.000002min","5min","10min","60.000004min"]
  for elt in times:
    if elt in x:
      return elt.replace("min","")
def statistic(ab, a):
    
    sumab, suma = sum(ab), sum(a)
    return ( suma / len(a) -
             (sumab -suma) / (len(ab) - len(a)) ) 

def colors(x):
    
  color = {"0.25":(1.0/6.0,0,0),"1":(2.0/6.0,0,0),"5":(3.0/6.0,0,0),"10":(4.0/6.0,0,0),"30.000002":(5.0/6.0,0,0),"60.000004":(1.0,0,0)}
  if x in color:
     return color[x]

def permutationTest(a, b):
    
    ab = a + b
    #print (a,b)
    Tobs = statistic(ab, a)
    under = 0
    for count, perm in enumerate(comb(ab, len(a)), 1):
        if statistic(ab, perm) <= Tobs:
            under += 1
    return (1-(float(under) / float(count)))

def burial(start,stop):
    infile = open("A1B2_pisa.txt","r") #Text file containing A1B2 data
    d = {}
    i1_sum = 0
    i2_sum = 0
    for line in infile:
    #  print (line)
      k = line.split()
      if float(k[3]) == 0 and float(k[4]) == 0:
        d[k[2]] = 0
      else:
        d[k[2]] = float(k[4])

      
    infile.close()
    for i in d:
       if int(i)>=start and int(i)<=stop:
        i2_sum += float(d[i])
    i2_average = i2_sum/float(stop-start)
    infile.close()

    infile = open("A1B1_pisa.txt","r")  #Text file containing A1B1 data
    d = {}
    for line in infile:
      line = line.strip()
      k = line.split()
      d[k[2]] = float(k[4])
    infile.close()
    for i in d:
       if int(i)>=start and int(i)<=stop:
        i1_sum += float(d[i])
    i1_average = i1_sum/float(stop-start)
    infile.close()
    print (start,stop,i1_average,i2_average)

    if i1_average>8 and i2_average>8:
      return "both"
    
    if i1_average>8 and not i2_average>8:
        return "interface1"
    if not i1_average>8 and i2_average>8:
      return "interface2"
            
    else:
      return "none"

def sd_means_difference(a,b):
    sd1 = np.std(a)
    sd2 = np.std(b)
    na = float(len(a))
    nb = float(len(b))
    return (((sd1**2)/na) + (((sd2**2)/nb)))**0.5


color = {"0.250000":(1.0/6.0,0,0),"1.000000":(2.0/6.0,0,0),"5.000000":(3.0/6.0,0,0),"10.000000":(4.0/6.0,0,0),"30.000002":(5.0/6.0,0,0),"60.000004":(1.0,0,0)}

peptides_monomer_values = {}
peptides_dimer_values = {}
for filename in os.listdir("C:/Users/Arvind Pillai/Downloads/repeat uptake data_Arvind_020719/repeat uptake data_Arvind_020719/"):
        print (filename)
        infile = open("C:/Users/Arvind Pillai/Downloads/repeat uptake data_Arvind_020719/repeat uptake data_Arvind_020719/"+filename,"r")
        for line in infile:
          print (line)
   
          a = line.split()
          if not "start" in line.lower() and "67" in filename and len(a)>0:
            if not (a[0]+","+a[1]+" "+whattime(filename)) in peptides_monomer_values:
              peptides_monomer_values[a[0]+","+a[1]+" "+whattime(filename)] = [float(a[2])]            
            else:
              peptides_monomer_values[a[0]+","+a[1]+" "+whattime(filename)].append(float(a[2]))

          if not "start" in line.lower() and "75" in filename and len(a)>0:
            if not (a[0]+","+a[1]+" "+whattime(filename)) in peptides_dimer_values:
              peptides_dimer_values[a[0]+","+a[1]+" "+whattime(filename)] = [float(a[2])]            
            else:
              peptides_dimer_values[a[0]+","+a[1]+" "+whattime(filename)].append(float(a[2]))
        infile.close()
threshhold = 0.5
#print (peptides_monomer_values)             
#print (peptides_dimer_values)

new_peptides_monomer_values = {}
new_peptides_dimer_values = {}
for elt in peptides_monomer_values:

    t,p = stats.ttest_1samp(peptides_monomer_values[elt],0.0)#1samp(peptide_differences[elt],0.0))
  #  print (elt, p)
    if p<0.01:
      new_peptides_monomer_values[elt] = peptides_monomer_values[elt]
    else:
      print ("Hello",elt,p)
for elt in peptides_dimer_values:

    t,p = stats.ttest_1samp(peptides_dimer_values[elt],0.0)#1samp(peptide_differences[elt],0.0))
 #   print (elt, p)
    if p<0.1:
      new_peptides_dimer_values[elt] = peptides_dimer_values[elt]

#peptides_dimer_values = new_peptides_dimer_values
peptides_monomer_values = new_peptides_monomer_values
peptide_errors = {}
peptide_differences = {}
for elt in peptides_monomer_values:
  if elt in peptides_dimer_values: #and np.mean(peptides_monomer_values[elt])>threshhold:
    peptide_errors[elt] = sd_means_difference(peptides_monomer_values[elt],peptides_dimer_values[elt])
    peptide_differences[elt] = []

    for i in range(len(peptides_monomer_values[elt])):
      monomer = peptides_monomer_values[elt][i]
      dimer = peptides_dimer_values[elt][i]
      peptide_differences[elt].append((monomer-dimer)/dimer)


#print (peptide_differences)
peptide_interface = {}
for elt in peptide_differences:
  key = elt.split()
  positions = key[0].split(",")
  peptide_interface[elt] = burial(float(positions[0]),float(positions[1]))

#print (peptide_interface)

for elt in peptide_differences:
  x = 0   
  if peptide_interface[elt] =="none":
    x = 1
  if peptide_interface[elt] == "interface1":
    x = 2
  if peptide_interface[elt] == "interface2":
    x = 3
  if peptide_interface[elt] == "both":
    x = 4
  key = elt.split()
  time = float(key[1])
  if time ==0.25:
    x += 0.1
  if time == 1.0:
    x += 0.2
  if time == 5.0:
    x += 0.3
  if time == 10.0:
    x += 0.4
  if time == 30.000002:
    x += 0.5
  if time == 60.000004:
    x += 0.6
    
#  print (x,interfaces[i],normalized_differences_threshhold[i])

  pylab.scatter(x,np.mean(peptide_differences[elt]), color=colors(key[1]))
  pylab.errorbar(x,np.mean(peptide_differences[elt]),yerr=peptide_errors[elt],color=colors(key[1]), linestyle="None")#stats.sem(peptide_differences[elt]),color=colors(key[1]), linestyle="None")
 # if time == 60.000004:
  #  pylab.text(x,np.mean(peptide_differences[elt]),elt)

s = "normalized_diff_thresh.eps"
#pylab.ylim(-0.2,0.7)
pylab.savefig(s)
pylab.close()

######################################################


times = ["0.25","5","10","30.000002","60.000004"]

for time in times:
 
 peptides = []

 for elt in peptides_monomer_values:
  if elt in peptides_dimer_values and np.mean(peptides_monomer_values[elt])>0.5:
    k = elt.split(" ")
    if k[1] == time:
      positions = k[0].split(",")
      peptides.append((float(positions[0]),float(positions[1])))          

 print (peptides)
 
 
   
 p_values = []
 p_values2 = []
 for i in range(1000):
 #  print (i)
   poppable = copy.deepcopy(peptides) 
   p = random.choice(peptides)
   
   poppable.pop(poppable.index(p))
   curated_peptides = []
   curated_peptides.append(p)
   while (len(poppable)>0):
     p = random.choice(poppable)
 
     overlapping = 0
     for elt in curated_peptides:
       if overlap(elt,p):
         overlapping +=1
     if overlapping == 0:
       curated_peptides.append(p)
     poppable.pop(poppable.index(p))

   interface2 = []
   interface1 = []
#  print (curated_peptides)
   for elt in curated_peptides:
     monomer = random.choice(peptides_monomer_values[str(elt[0])+","+str(elt[1])+" "+str(time)])#random.choice(peptides_monomer_values[elt+" "+str(time)])
     dimer = random.choice(peptides_dimer_values[str(elt[0])+","+str(elt[1])+" "+str(time)])#random.choice(peptides_dimer_values[elt+" "+str(time)])
     if burial(elt[0],elt[1]) == "interface2":
       
       interface1.append(float(monomer-dimer)/float(monomer))
   #   report+=1
     if burial(elt[0],elt[1]) == "none" or burial(elt[0],elt[1]) == "interface1":
       interface2.append(float(monomer-dimer)/float(monomer))      
 # print (len(interface1),len(interface2))       
   if len(interface2)>0 and len(interface1)>0:                             
     p_values.append(permutationTest(interface1,interface2))

   interface2 = []
   interface1 = []
  
   for elt in curated_peptides:
     monomer = random.choice(peptides_monomer_values[str(elt[0])+","+str(elt[1])+" "+str(time)])#random.choice(peptides_monomer_values[elt+" "+str(time)])
     dimer = random.choice(peptides_dimer_values[str(elt[0])+","+str(elt[1])+" "+str(time)])#random.choice(peptides_dimer_values[elt+" "+str(time)])

     if burial(elt[0],elt[1]) == "interface1":
       interface1.append((monomer-dimer)/dimer)   #   report+=1
     if burial(elt[0],elt[1]) == "none" or burial(elt[0],elt[1]) == "interface2":
       interface2.append((monomer-dimer)/dimer) # print (len(interface1),len(interface2))       
   if len(interface2)>0 and len(interface1)>0:                             
     p_values2.append(permutationTest(interface1,interface2))




 print (time, np.mean(p_values),np.mean(p_values2))
 
 
 pylab.hist(p_values,bins=np.arange(0.0, 1.0, 0.02),color="red",normed=True)
 pylab.hist(p_values2,bins=np.arange(0.0, 1.0, 0.02),color="blue",normed=True),pylab.axvline(0.05, color='k', linestyle='dashed', linewidth=1)
 s = "p_values2 "+time+".eps"
 pylab.title("Histogram of p-values, Interface 1 vs all, Interface 2 vs all")
 pylab.xlabel("p-values")
 exchange = []
 buried = []
 pylab.savefig(s)
 pylab.close()

#print (peptides_monomer_values)             
#print (peptides_dimer_values
print (peptides)
buried_i1 = []
buried_i2 = []
  
for elt in peptide_differences:
 if "60.000004" in elt:
  exchange.append(np.mean(peptide_differences[elt]))
  buried_i2.append(burial_value(elt,"A1B2_PISA.txt"))
  buried_i1.append(burial_value(elt,"A1B1_HADDOCK.txt"))
 # peptide.append(elt)
  
pylab.scatter(buried_i1, exchange,color="black")
par = np.polyfit(buried_i1, exchange, 1, full=True)

slope=par[0][0]
intercept=par[0][1]
y1 = [slope*xx + intercept  for xx in buried_i1]
pylab.plot(buried_i1, y1, '-r')
slope, intercept, r_value, p_value, std_err = stats.linregress(buried_i1,exchange)

s = "i1_correlation.eps"
pylab.savefig(s)
pylab.close()
print("r_value",r_value)
pylab.scatter(buried_i2, exchange,color="black")
par = np.polyfit(buried_i2, exchange, 1, full=True)

slope=par[0][0]
intercept=par[0][1]
y1 = [slope*xx + intercept  for xx in buried_i2]
pylab.plot(buried_i2, y1, '-r')
slope, intercept, r_value, p_value, std_err = stats.linregress(buried_i2,exchange)

s = "i2_correlation.eps"
pylab.savefig(s)
pylab.close()
print("r_value",r_value)
