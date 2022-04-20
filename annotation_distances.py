#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Imports
import sys
import urllib


from Bio import SeqIO

# Global variables


# Function definitions

def loadnamevariants():
	output = {}
	url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
	for line in urllib.request.urlopen(url):
		line = line.decode('utf-8').strip()
		name = line.split(";")[0]
		annotype = line.split(":")[0].split(";")[1]
		variants = line.split(":")[1].split(",")
		for v in variants:
			for g in ['', ' ']:
				v = v.replace(g, '')
				for s in ['',' GENE', ' '+annotype.upper()]:
					output[v+s] = name
	return(output)


def get_feat_name(feat):
	featname = "unknown"
	nametags = ['gene', 'product', 'label', 'standard_name']
	if(any(t in feat.qualifiers.keys() for t in nametags)):
		for t in nametags:
			if(t in feat.qualifiers.keys()):
				featname = feat.qualifiers[t][0].upper()
				break
	return(featname)


# Main

if __name__ == "__main__":
	
	namevariants = loadnamevariants()
	namevariants['unknown'] = 'unknown'
	
	sys.stdout.write("sequence\tname1\tpos1\tname2\tpos2\n")
	
	for seq_record in SeqIO.parse(sys.stdin, "genbank"):
		seqname = seq_record.name
		
		feats = seq_record.features
		
		feats = [feat for feat in feats if feat.type not in ['source', 'misc_feature', 'repeat_region', 'D-loop', 'rep_origin','gap']]
		
		for i, feat in enumerate(feats):
			#i = 5
			#feat = seq_record.features[i]
			
			if(feat.location == feats[i-1].location):
				continue
			prev_name, prev_end = [namevariants[get_feat_name(feats[i-1])], str(int(feats[i-1].location.end))] if i > 0 else ['NA', 'NA']
			
			line = "\t".join([seqname, prev_name, prev_end, namevariants[get_feat_name(feat)], str(int(feat.location.start))])
			
			sys.stdout.write(line + "\n")
			
	exit()