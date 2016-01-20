import re
import statistics as stats
import math

prim_len = 31

while True:
	mut_seq = raw_input("Enter the target DNA sequence for mutation: ").upper()
	mut_seq_len = len(mut_seq)
	if not re.match("^[A, T, C, G]*$", mut_seq):
		print "Input is restricted to nucleotide symbols (A, T, C or G)"
	elif mut_seq_len < prim_len:
		print ("The input sequence is too short, it must equal or exceed " 
			   + str(prim_len) + " bases")
	else:
		break

percent_comp = []
nucleotides = [mut_seq.count("A"), mut_seq.count("T"),
		   mut_seq.count("C"), mut_seq.count("G")]
for n in nucleotides:
	percent_comp.append("%.1f" %(float(n) / mut_seq_len * 100))

nuc_content = tuple(percent_comp)

print "The sequence to mutatate is", mut_seq_len, '''bp with the following 
nucleotide composition: %s%% A, %s%% T, %s%% C and %s%% G.''' %(nuc_content)

#mut_seq = "ATGCGATCGTAGCTGTTACGAGCCAGTTGAAATGTCGAATGCTATGGCCCATAAGCTTA"

while True:
	mut_target = raw_input("Select nucleotide type "
						   "to replace (A, T, C or G): ").upper()
	if not re.match("^[A, T, C, G]*$", mut_target):
		print "Base type is restricted to nucleotide symbols (A, T, C or G)"
	elif len(mut_target) > 1:
		print "Input a single base type, with no additional spaces/characters"
	else: 
		break

while True:
	mut_choice = raw_input("Select nucleotide type to insert "
						   "(A, T, C or G): ").upper()
	if not re.match("^[A, T, C, G]*$", mut_choice):
		print "Base type is restricted to nucleotide symbols (A, T, C or G)"
	elif mut_choice == mut_target:
		print "Select a nucleotide other than " + (mut_target)
	elif len(mut_choice) > 1:
		print "Input a single base type, with no additional spaces/characters"
	else: 
		break

prim_name = []
prim_bases = []
start_loc = 1
for i in range (0, len(mut_seq) - prim_len):
	end_loc = start_loc + prim_len
	prim_bases.append(mut_seq[start_loc:end_loc])
	label = (mut_target, str(end_loc - (prim_len / 2)), mut_choice, "_forward")
	prim_name.append("".join(label))
	start_loc += 1

while True:
	sel_prim_idx = []
	selected_primers = []
	mut = 0
	for idx, val in enumerate(prim_bases):
		if val[prim_len / 2] == mut_target:
			q = list(val)
			q[prim_len / 2] = mut_choice
			selected_primers.append("".join(q))
			sel_prim_idx.append(idx)
			mut += 1
	if mut < 1:
		print ("No nucleotides matching your mutagenesis selection " + 
			   "criteria were idenitifed in the specified sequence.")
		break

	primer_Tm = []
	for w in selected_primers:
		A = w.count("A")
		T = w.count("T")
		G = w.count("G")
		C = w.count("C")
		Tm = 64.9 + 41 * (G + C - 16.4) / (A + T + G + C)
		primer_Tm.append(float("%.1f" % Tm))
	# Makes the following assumptions:-
	#     annealing is at 50 nm primer, 50 mM Na+ and pH 7.0

	opt_primers = []
	multi_temp = []
	for idx, val in enumerate(primer_Tm):
		if (val - stats.mean(primer_Tm)) ** 2 <= 4:
			multi_temp.append(val)
			opt_primers.append(idx)

	print ("The Tm for application of this primer set is " 
		   + "%.1f" % (stats.mean(multi_temp) - 2) + u"\u2103" + ".")
	# Makes the following assumptions:-
	#     the mismatch will reduce the working Tm by ~2 deg.C
	#     which is incorporated into the displayed temperature

	for j in opt_primers:
		print prim_name[sel_prim_idx[j]]
		print selected_primers[j]
		rev_com = []
		v = list(reversed(selected_primers[j]))
		for i in v:
		 	if i == "A":
		 		i = "T"
			elif i == "T":
				i = "A"
			elif i == "G":
				i = "C"
			elif i == "C":
				i = "G"
			rev_com.append(i)
		print prim_name[sel_prim_idx[j]][0:-8] + "_reverse"
		print ("".join(rev_com))
	break







