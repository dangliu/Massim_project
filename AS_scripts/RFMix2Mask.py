#!/usr/bin/python3

usage = """
This script is for producing masked genotype files (.geno format) from RFMix v1.54 outputs.
It was written by Dang Liu. Last updated: May 23 2022.
usage:
python3 RFMix2Mask.py Alleles.txt ForwardBackward.txt prob_cutoff

output:
anc1masked.geno anc2masked.geno

"""

# modules here
import sys, re, math
import itertools

# if input number is incorrect, print the usage
if len(sys.argv) < 4:
        print(usage)
        sys.exit()

# Prob_cutoff
prob_cutoff = float(sys.argv[3])

# Output files
out_f_anc1M = open(sys.argv[2].replace(".ForwardBackward.txt","") + ".anc1masked.geno", "w")
out_f_anc2M = open(sys.argv[2].replace(".ForwardBackward.txt","") + ".anc2masked.geno", "w")
#out_f_anc1E = open(sys.argv[2].replace(".ForwardBackward.txt","") + ".anc1extract.geno", "w")
#out_f_anc2E = open(sys.argv[2].replace(".ForwardBackward.txt","") + ".anc2extract.geno", "w")

# Read files
A_f = open(sys.argv[1], 'r') # alleles
#V_f = open(sys.argv[2], 'r')
F_f = open(sys.argv[2], 'r') # forward backward prob

# Remove headers
line_A_f = A_f.readline()
#line_V_f = V_f.readline()
line_F_f = F_f.readline()

print("There are " + str(len(list(line_A_f.replace('\n',"")))) + " haplotypes!")

# Use to count how many SNPs
n=1


while(line_A_f):
	line_A_f_s = list(line_A_f.replace('\n',""))
	#print(len(line_A_f_s))
	#print(line_A_f_s[-1])
	#line_V_f_s = re.split(r'\s+', line_V_f.replace(' \n',""))
	#print(len(line_V_f_s))
	#print(line_V_f_s[-1])
	line_F_f_s = re.split(r'\s+', line_F_f.replace(' \n',""))
	#print(len(line_F_f_s))
	#print(line_F_f_s[-1])
	#print("Reading " + str(len(line_A_f_s)) + " haplotypes!")
	# Read the two haplotypes per individual over the A, V, F inputs
	for A_a, A_b, F_a1, F_a2, F_b1, F_b2 in zip(line_A_f_s[0::2], line_A_f_s[1::2], line_F_f_s[0::4], line_F_f_s[1::4], line_F_f_s[2::4], line_F_f_s[3::4]):
		#print(A_a, A_b, F_a1, F_a2, F_b1, F_b2, sep="\t")
		# Mask a particular ancestry's genotype (by writting it to missing data in the eigenstrat geno format '9'), if:
		# 1. the forward backward probs of that ancestry passes the prob_cutoff for both alleles
		# 2. the forward backward prob of that ancestry passes the prob_cutoff for one of the alleles and the other ancestry prob passes for the other allele
		if (float(F_a1)>=prob_cutoff and float(F_b1)>=prob_cutoff):
			if (int(A_a)==0 and int(A_b)==0):
				print(9, end = '', file=out_f_anc1M)
				print(0, end = '', file=out_f_anc2M)
				#print(0, end = '', file=out_f_anc1E)
				#print(9, end = '', file=out_f_anc2E)				
			elif (int(A_a)!=int(A_b)):
				print(9, end = '', file=out_f_anc1M)
				print(1, end = '', file=out_f_anc2M)
				#print(1, end = '', file=out_f_anc1E)
				#print(9, end = '', file=out_f_anc2E)				
			else:
				print(9, end = '', file=out_f_anc1M)
				print(2, end = '', file=out_f_anc2M)
				#print(2, end = '', file=out_f_anc1E)
				#print(9, end = '', file=out_f_anc2E)				
		elif ((float(F_a1)>=prob_cutoff and float(F_b2)>=prob_cutoff) or (float(F_a2)>=prob_cutoff and float(F_b1)>=prob_cutoff)):
			if (int(A_a)==0 and int(A_b)==0):
				print(9, end = '', file=out_f_anc1M)
				print(9, end = '', file=out_f_anc2M)
				#print(0, end = '', file=out_f_anc1E)
				#print(0, end = '', file=out_f_anc2E)
			elif (int(A_a)!=int(A_b)):
				print(9, end = '', file=out_f_anc1M)
				print(9, end = '', file=out_f_anc2M)
				#print(1, end = '', file=out_f_anc1E)
				#print(1, end = '', file=out_f_anc2E)
			else:
				print(9, end = '', file=out_f_anc1M)
				print(9, end = '', file=out_f_anc2M)
				#print(2, end = '', file=out_f_anc1E)
				#print(2, end = '', file=out_f_anc2E)	
		elif (float(F_a2)>=prob_cutoff and float(F_b2)>=prob_cutoff):
			if (int(A_a)==0 and int(A_b)==0):
				print(0, end = '', file=out_f_anc1M)
				print(9, end = '', file=out_f_anc2M)
				#print(9, end = '', file=out_f_anc1E)
				#print(0, end = '', file=out_f_anc2E)
			elif (int(A_a)!=int(A_b)):
				print(1, end = '', file=out_f_anc1M)
				print(9, end = '', file=out_f_anc2M)
				#print(9, end = '', file=out_f_anc1E)
				#print(1, end = '', file=out_f_anc2E)
			else:
				print(2, end = '', file=out_f_anc1M)
				print(9, end = '', file=out_f_anc2M)
				#print(9, end = '', file=out_f_anc1E)
				#print(2, end = '', file=out_f_anc2E)
		else:
			if (int(A_a)==0 and int(A_b)==0):
				print(0, end = '', file=out_f_anc1M)
				print(0, end = '', file=out_f_anc2M)
				#print(9, end = '', file=out_f_anc1E)
				#print(9, end = '', file=out_f_anc2E)
			elif (int(A_a)!=int(A_b)):
				print(1, end = '', file=out_f_anc1M)
				print(1, end = '', file=out_f_anc2M)
				#print(9, end = '', file=out_f_anc1E)
				#print(9, end = '', file=out_f_anc2E)
			else:
				print(2, end = '', file=out_f_anc1M)
				print(2, end = '', file=out_f_anc2M)
				#print(9, end = '', file=out_f_anc1E)
				#print(9, end = '', file=out_f_anc2E)
	# Old version
	# for A_a, A_b, V_a, V_b, F_a1, F_a2, F_b1, F_b2 in zip(line_A_f_s[0::2], line_A_f_s[1::2], line_V_f_s[0::2], line_V_f_s[1::2], line_F_f_s[0::4], line_F_f_s[1::4], line_F_f_s[2::4], line_F_f_s[3::4]):
	# 	print(A_a, A_b, V_a, V_b, F_a1, F_a2, F_b1, F_b2, sep="\t")
	# 	if (int(V_a)==1 and int(V_b)==1):
	# 		if (float(F_a1)>=prob_cutoff and float(F_b1)>=prob_cutoff):
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(0, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(1, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(2, end = '', file=out_f_anc2M)
	# 		else:
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(0, end = '', file=out_f_anc1M)
	# 				print(0, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(1, end = '', file=out_f_anc1M)
	# 				print(1, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(2, end = '', file=out_f_anc1M)
	# 				print(2, end = '', file=out_f_anc2M)
	# 	elif (int(V_a)!=int(V_b)):
	# 		if (float(F_a1)>=prob_cutoff and float(F_b2)>=prob_cutoff):
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 		else:
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(0, end = '', file=out_f_anc1M)
	# 				print(0, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(1, end = '', file=out_f_anc1M)
	# 				print(1, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(2, end = '', file=out_f_anc1M)
	# 				print(2, end = '', file=out_f_anc2M)				
	# 	elif (int(V_a)==1 and int(V_b)==2):
	# 		if (float(F_a1)>=prob_cutoff and float(F_b2)>=prob_cutoff):
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 		else:
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(0, end = '', file=out_f_anc1M)
	# 				print(0, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(1, end = '', file=out_f_anc1M)
	# 				print(1, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(2, end = '', file=out_f_anc1M)
	# 				print(2, end = '', file=out_f_anc2M)			
	# 	elif (int(V_a)==2 and int(V_b)==1):
	# 		if (float(F_a2)>=prob_cutoff and float(F_b1)>=prob_cutoff):
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(9, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 		else:
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(0, end = '', file=out_f_anc1M)
	# 				print(0, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(1, end = '', file=out_f_anc1M)
	# 				print(1, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(2, end = '', file=out_f_anc1M)
	# 				print(2, end = '', file=out_f_anc2M)
	# 	else:
	# 		if (float(F_a2)>=prob_cutoff and float(F_b2)>=prob_cutoff):
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(0, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(1, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(2, end = '', file=out_f_anc1M)
	# 				print(9, end = '', file=out_f_anc2M)
	# 		else:
	# 			if (int(A_a)==0 and int(A_b)==0):
	# 				print(0, end = '', file=out_f_anc1M)
	# 				print(0, end = '', file=out_f_anc2M)
	# 			elif (int(A_a)!=int(A_b)):
	# 				print(1, end = '', file=out_f_anc1M)
	# 				print(1, end = '', file=out_f_anc2M)
	# 			else:
	# 				print(2, end = '', file=out_f_anc1M)
	# 				print(2, end = '', file=out_f_anc2M)
	
	# Move to a new line at the end
	print("\n", end='', file=out_f_anc1M)
	print("\n", end='', file=out_f_anc2M)
	#print("\n", end='', file=out_f_anc1E)
	#print("\n", end='', file=out_f_anc2E)
	n+=1			
	line_A_f = A_f.readline()
	#line_V_f = V_f.readline()
	line_F_f = F_f.readline()
A_f.close()
#V_f.close()
F_f.close()
out_f_anc1M.close()
out_f_anc2M.close()
#out_f_anc1E.close()
#out_f_anc2E.close()

print("Done! Total " + str(n) + " sites processed!")



#last_v220523