import re
import statistics as stats
import math

prim_len = 31

while True:
	mut_seq = raw_input("Enter the target DNA sequence for mutation: ")
	mut_seq_len = len(mut_seq)
	bases = mut_seq.upper()
	if not re.match("^[A, T, C, G]*$", bases):
		print "Input is restricted to nucleotide symbols (A, T, C or G)"
	elif mut_seq_len < prim_len:
		print ("The input sequence is too short, it must equal or exceed " 
			   + str(prim_len) + " bases.")
	else:
		break

percent_comp = []
nucleotides = [bases.count("A"), bases.count("T"),
		   bases.count("C"), bases.count("G")]
for n in nucleotides:
	percent_comp.append("%.1f" %(float(n) / mut_seq_len * 100))

nuc_content = tuple(percent_comp)

print "The sequence to mutatate is", mut_seq_len, '''bp with the following 
nucleotide composition: %s%% A, %s%% T, %s%% C and %s%% G.''' %(nuc_content)

#bases = "ATGCGATCGTAGCTGTTACGAGCCAGTTGAAATGTCGAATGCTATGGCCCATAAGCTTA"

prim_name = []
prim_bases = []
start_loc = 1
for i in range (0, len(bases) - prim_len):
	end_loc = start_loc + prim_len
	prim_bases.append(bases[start_loc:end_loc])
	label = ("A", str(end_loc - (prim_len / 2)), "T_forward")
	prim_name.append("".join(label))
	start_loc += 1

sel_prim_idx = []
selected_primers = []
for idx, val in enumerate(prim_bases):
	if val[prim_len / 2] == "A":
		q = list(val)
		q[prim_len / 2] = "T"
		selected_primers.append("".join(q))
		sel_prim_idx.append(idx)

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








