#!/zzh_gpfs/apps/python/bin/python
import os,sys,csv
from Bio import SeqIO
from glob import glob
from math import *

def create_cytoband_fole(vfile,dfile,jfile,outname): # germline fasta files
	handle = open("%s.CytoBand.txt"%outname,'wb')
	out = csv.writer(handle,delimiter='\t')
	out.writerow(["gene","Start","End","isotype","Stain"])
	myset = set()
	colors = ["gneg","acen","stalk","gpos100","gpos25","stalk","acen","gpos75","gneg","stalk"]
	n=1000000
	bb = ""
	for ind,file in enumerate([vfile,dfile,jfile]):
		for rec in SeqIO.parse(file,"fasta"):
			myid = rec.id.split("*")[0]
			if myid[2:4]!=bb:
				i = 0
				bb=myid[2:4]
			family = int(myid[4])
			if myid not in myset:
				out.writerow([myid[2:4],10*i*n,10*i*n+10*n,myid,colors[family]])
				i += 1
				myset.add(myid)
	handle.close()
def create_label_file(cytodict,outname):
	out1 = csv.writer(open('%s.Label.txt'%outname,'wb'),delimiter='\t')
	header = ["gene","Start","End","isotype"]
	for ind,file in enumerate([hvexpre,hdexpre,hjexpre,kvexpre,kjexpre,lvexpre,ljexpre]):
		vnum,dnum,jnum = 0,0,0
		for rec in csv.reader(open(file,'rU'),delimiter='\t'):
			if ind == 0 and rec[0].startswith("gene"):
				out1.writerow(header)
				continue
			if ind != 0 and rec[0].startswith("gene"):
				continue
			if "D" in cytodict[rec[0]][0]:
				if dnum < 20:
					out1.writerow(cytodict[rec[0]])
					dnum += 1
			elif "J" in cytodict[rec[0]][0]:
				if jnum < 3:
					out1.writerow(cytodict[rec[0]])
					jnum += 1
			elif "LV" in cytodict[rec[0]][0]:
				if vnum < 32:
					out1.writerow(cytodict[rec[0]])
					vnum += 1
			elif "KV" in cytodict[rec[0]][0]:
				if vnum < 47:
					out1.writerow(cytodict[rec[0]])
					vnum += 1
			elif "HV" in cytodict[rec[0]][0]:
				if vnum < 60:
					out1.writerow(cytodict[rec[0]])
					vnum += 1
def create_Heatmap_file(cytodict,outname):
	out2 = csv.writer(open('%s.Heatmap.txt'%outname,'wb'),delimiter='\t')
	header = ["gene","Start","End","isotype"]
	for ind,file in enumerate([hvexpre,hdexpre,hjexpre,kvexpre,kjexpre,lvexpre,ljexpre]):
		for rec in csv.reader(open(file,'rU'),delimiter='\t'):
			if ind == 0 and rec[0].startswith("gene"):
				out2.writerow(header+[rec[-2]])
				continue
			if ind != 0 and rec[0].startswith("gene"):
				continue
			out2.writerow(cytodict[rec[0]]+[log(1000*float(rec[-2])+1,2)])
def create_link_file(hrecom,krecom,lrecom,cytodict,outname):
	out3 = csv.writer(open('%s.Link.txt'%outname,'wb'),delimiter='\t')
	header = ["gene","Start","End","gene","Start","End"]
	out3.writerow(header)
	for rec in csv.reader(open(hrecom,'rU'),delimiter='\t'):
		if rec[0].startswith("gene"):
			continue
		if 'N/A' in rec[0]:
			continue
		if float(rec[-2])<0.001:
			continue
		genes = rec[0].split(":")
		out3.writerow(cytodict[genes[0]][:-1]+cytodict[genes[1]][:-1])
		out3.writerow(cytodict[genes[1]][:-1]+cytodict[genes[2]][:-1])
	for rec in csv.reader(open(krecom,'rU'),delimiter='\t'):
		if rec[0].startswith("gene"):
			continue
		if 'N/A' in rec[0]:
			continue
		if float(rec[-2])<0.001:
			continue
		genes = rec[0].split(":")
		out3.writerow(cytodict[genes[0]][:-1]+cytodict[genes[1]][:-1])
	for rec in csv.reader(open(lrecom,'rU'),delimiter='\t'):
		if rec[0].startswith("gene"):
			continue
		if 'N/A' in rec[0]:
			continue
		if float(rec[-2])<0.001:
			continue
		genes = rec[0].split(":")
		out3.writerow(cytodict[genes[0]][:-1]+cytodict[genes[1]][:-1])
	
def main():
	create_cytoband_fole(vfile,dfile,jfile,outname)
	cytodict = {}
	for rec in csv.reader(open("%s.CytoBand.txt"%outname,'rU'),delimiter='\t'):
		if rec[0].startswith('gene'):
			continue
		cytodict[rec[3]] = rec[:4]
	create_label_file(cytodict,outname)
	create_Heatmap_file(cytodict,outname)
	create_link_file(hrecom,krecom,lrecom,cytodict,outname)
	
if __name__=='__main__':
	vfile = "InputFiles/human_gl_v"
	dfile = "InputFiles/human_gl_d"
	jfile = "InputFiles/human_gl_j"
	hvexpre = "InputFiles/expression.IgH.V.unique.sort.txt"
	hdexpre = "InputFiles/expression.IgH.D.unique.sort.txt"
	hjexpre = "InputFiles/expression.IgH.J.unique.sort.txt"
	kvexpre = "InputFiles/expression.IgK.V.unique.sort.txt"
	kjexpre = "InputFiles/expression.IgK.J.unique.sort.txt"
	lvexpre = "InputFiles/expression.IgL.V.unique.sort.txt"
	ljexpre = "InputFiles/expression.IgL.J.unique.sort.txt"
	hrecom = "InputFiles/recombination.IgH.VDJ.unique.sort.txt"
	krecom = "InputFiles/recombination.IgK.VJ.unique.sort.txt"
	lrecom = "InputFiles/recombination.IgL.VJ.unique.sort.txt"
	outname = "OutputFiles/human.VDJ"
	main()
