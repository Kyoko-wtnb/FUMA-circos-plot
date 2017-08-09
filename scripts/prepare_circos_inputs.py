#!/usr/bin/python
import sys
import os
import subprocess
import re
import pandas as pd
import numpy as np
import math
import argparse
import time

##### Return index of a1 which exists in a2 #####
def ArrayIn(a1, a2):
	results = np.where(np.in1d(a1, a2))[0]
	return results

##### return unique element in a list #####
def unique(a):
	unique = []
	[unique.append(s) for s in a if s not in unique]
	return unique

##### get arguments #####
parser = argparse.ArgumentParser()
parser.add_argument('--configdir', default=None, type=str, help="Directory of base.conf file.")
parser.add_argument('--indir', default=None, type=str, help="Directory of text files downloded from FUMA. Use 'configdir' if not provided.")
parser.add_argument('--outdir', default=None, type=str, help="Directory of output files. Use 'configdir' if not provided.")
parser.add_argument('--sumstats', default=None, type=str, help="Path to GWAS summary statistics file. If not provided, only SNPs in the text file downloaded from FUMA will be used.")
parser.add_argument('--delim', default='\t', type=str, help="Delimiter of GWAS summary statistics. Default is tab.")
parser.add_argument('--chrcol', default=None, type=str, help="Header name of chromosome column in the sumstats file. Mandatory if --sumstats argument is given.")
parser.add_argument('--poscol', default=None, type=str, help="Header name of SNP position column in the sumstats file. Mandatory if --sumstats argument is given.")
parser.add_argument('--pcol', default=None, type=str, help="Header name of P-value column in the sumstats file. Mandatory if --sumstats argument is given.")
parser.add_argument('--chrom', default=None, type=str, help="Chromosome index to create circos plot. All chromosomes with at least one risk locus will be processed by default. Multiple chromsomes can be specified sepated by comma (e.g. --chrom 1,3,6).")
parser.add_argument('--max-N-snps', default=150000, type=int, help="The maximum number of SNPs per plot. Default is 150000.")
parser.add_argument('--max-snp-p', default=0.05, type=float, help="The maximum P-value of SNPs to plot. Default is 0.05.")
parser.add_argument('--max-N-links', default=100000, type=int, help="The maximum number of links (eQTLs/chromatin interactions) to plot. Default is 100000.")

##### create condig files #####
def createConfig(args, c, loci, ci, snps, allsnps, genes):
	regions = []
	breaks = ""
	loci = loci[loci[:,3].argsort()]

	loci = np.c_[loci, loci[:,4], loci[:,5]]
	for l in ci:
		if min(l[2], l[5]) < loci[loci[:,0]==l[0],6]:
			loci[loci[:,0]==l[0],6] = min(l[2], l[5])
		if max(l[3], l[6]) > loci[loci[:,0]==l[0],7]:
			loci[loci[:,0]==l[0],7] = max(l[3], l[6])
	for l in genes:
		if l[1] < loci[loci[:,0]==int(l[4]),6]:
			loci[loci[:,0]==int(l[4]),6] = l[1]
		if l[2] > loci[loci[:,0]==int(l[4]),7]:
			loci[loci[:,0]==int(l[4]),7] = l[2]
	cur_pos = 0
	tmp_start = []
	tmp_end = []
	for l in loci:
		if cur_pos == 0:
			if int((l[6]-1000)/1000000)<=0:
				tmp_start.append(0)
			else:
				breaks = "-hs"+str(c)+":0-"+str(int((l[6]-1000)/1000000)-1)
				tmp_start.append((int((l[6]-1000)/1000000)-1)*1000000)
			cur_pos = l[7]
		elif (int((l[6]-1000)/1000000)-1)-(int((cur_pos+1000)/1000000)+1) <= 1:
			cur_pos = max(cur_pos, l[7])
		else:
			if len(breaks) > 0:
				breaks += ";"
			breaks += "-hs"+str(c)+":"+str(int((cur_pos+1000)/1000000)+1)+"-"+str(int((l[6]-1000)/1000000)-1)
			tmp_end.append((int((cur_pos+1000)/1000000)+1)*1000000)
			tmp_start.append((int((l[6]-1000)/1000000)-1)*1000000)
			cur_pos = l[7]
	breaks += ";-hs"+str(c)+":"+str(int((cur_pos+1000)/1000000)+1)+"-)"
	tmp_end.append((int((cur_pos+1000)/1000000)+1)*1000000)
	regions = np.c_[tmp_start, tmp_end]

	tmp_snps = []
	if len(allsnps)>0:
		allsnps = allsnps[allsnps[:,0].astype(int)==c]
		for l in regions:
			tmp = allsnps[np.where((allsnps[:,1].astype(int)>=l[0]) & (allsnps[:,1].astype(int)<=l[1]))]
			if len(tmp_snps)==0:
				tmp_snps = tmp
			else:
				tmp_snps = np.r_[tmp_snps, tmp]
	if len(tmp_snps)>0:
		tmp_snps = np.c_[tmp_snps, [0]*len(tmp_snps)]
		snps = np.r_[snps, tmp_snps]
	snps[:,2] = [float(-1*x) for x in np.log10(snps[:,2].astype(float))]

	##### take top 150000 SNPs per chromosome #####
	if len(snps) > args.max_N_snps:
		snps = snps[snps[:,2].argsort()[::-1]]
		snps = snps[0:args.max_N_snps]
		snps = snps[snps[:,1].argsort()]

	maxlogP = int(max(snps[:,2]))+1
	minlogP = 0
	snps[:,0] = ["hs"+str(x) for x in snps[:,0]]
	snps = np.c_[snps[:,0:2], [x+1 for x in snps[:,1].astype(int)], snps[:,2:]]
	for l in snps:
		if float(l[4]) >= 0.8:
			l[4] = "id=1"
		elif float(l[4]) >= 0.6:
			l[4] = "id=2"
		elif float(l[4]) >= 0.4:
			l[4] = "id=3"
		elif float(l[4]) >= 0.2:
			l[4] = "id=4"
		else:
			l[4] = "id=5"
	with open(args.configdir+"/base.conf", 'r') as fin:
		cfg = fin.read();
	cfg = cfg.replace("[chr]", str(c))
	cfg = cfg.replace("[breaks]", breaks)
	cfg = cfg.replace("[maxlogP]", str(maxlogP))
	cfg = cfg.replace("[minlogP]", "0")

	with open(args.outdir+"/circos_chr"+str(c)+".conf", 'w') as o:
		o.write(cfg)

	regions = np.c_[[c]*len(regions), regions]
	return [snps, regions];

def main(args):
	start_time = time.time()
	##### check argument #####
	if args.configdir is None:
		parser.print_help()
		sys.exit('\nERROR: The --configdir flag is required.')

	if args.sumstats is not None:
		if args.chrcol is None or args.poscol is None or args.pcol is None:
			parser.print_help()
			sys.exit('\nERROR: The --chrcol, --poscol and --pcol flags are required when --sumstats is specified.')

	if not os.path.isabs(args.configdir):
		args.configdir = os.getcwd()+"/"+args.configdir

	if args.indir is None:
		args.indir = args.configdir
	if args.outdir is None:
		args.outdir = args.configdir

	if not os.path.isabs(args.indir):
		args.indir = os.getcwd()+"/"+args.indir
	if not os.path.isabs(args.outdir):
		args.outdir = os.getcwd()+"/"+args.outdir
	if args.sumstats is not None and not os.path.isabs(args.sumstats):
		args.sumstats = os.getcwd()+"/"+args.sumstats

	##### check files #####
	if not os.path.isfile(args.configdir+"/base.conf"):
		sys.exit('ERROR: '+args.configdir+'/base.conf does not exist.')
	if not os.path.isfile(args.indir+"/GenomicRiskLoci.txt"):
		sys.exit('ERROR: '+args.indir+'/GenomicRiskLoci.txt does not exist.')
	if not os.path.isfile(args.indir+"/snps.txt"):
		sys.exit('ERROR: '+args.indir+'/snps.txt does not exist.')
	if not os.path.isfile(args.indir+"/genes.txt"):
		sys.exit('ERROR: '+args.indir+'/genes.txt does not exist.')

	##### prepare directory #####
	if not os.path.isdir(args.outdir):
		os.makedirs(args.outdir)

	##### risk loci #####
	loci = pd.read_table(args.indir+"/GenomicRiskLoci.txt", delim_whitespace=True)
	loci = np.array(loci)
	loci = loci[:,[0,2,3,4,6,7]] #loci,rsID,chr,pos,start,end

	##### snps #####
	snps = pd.read_table(args.indir+"/snps.txt", delim_whitespace=True)
	snpshead = list(snps.columns.values)
	snps = np.array(snps)
	snps = snps[:,[2,3,7,snpshead.index("r2")]]
	snps = snps[np.where(np.isfinite(snps[:,2].astype(float)))]

	##### mapped genes #####
	genes = pd.read_table(args.indir+"/genes.txt", delim_whitespace=True)
	geneshead = list(genes.columns.values)
	genes = np.array(genes)

	##### 3D genome  #####
	ci = []
	if os.path.isfile(args.indir+"/ci.txt"):
		ci = pd.read_table(args.indir+"/ci.txt", delim_whitespace=True)
		ci = np.array(ci)
		ci = ci[ci[:,0].argsort()]
		ci = ci[ci[:,7]=="intra"]
		chr1 = [int(x.split(":")[0]) for x in ci[:,1]]
		chr2 = [int(x.split(":")[0]) for x in ci[:,2]]
		pos1min = [int(x.split(":")[1].split("-")[0]) for x in ci[:,1]]
		pos1max = [int(x.split(":")[1].split("-")[1]) for x in ci[:,1]]
		pos2min = [int(x.split(":")[1].split("-")[0]) for x in ci[:,2]]
		pos2max = [int(x.split(":")[1].split("-")[1]) for x in ci[:,2]]
		ci = np.c_[ci[:,0], chr1, pos1min, pos1max, chr2, pos2min, pos2max, ci[:,3:7]]
		### take top 100000 links per chromosome
		ci_chrom = unique(ci[:,1])
		ci_tmp = []
		for c in ci_chrom:
			tmp = ci[ci[:,1]==c]
			if len(tmp)>args.max_N_links:
				tmp = tmp[tmp[:,7].astype(float).argsort()]
				tmp = tmp[0:args.max_N_links]
			if len(ci_tmp)==0:
				ci_tmp = tmp
			else:
				ci_tmp = np.r_[ci_tmp, tmp]
		ci = ci_tmp
	else:
		print("WARNING: "+args.indir+"/ci.txt does not exit. "
		+"Circos plot will not include chromatin interactions. "
		+"Please download from FUMA to include them.")

	##### eqtl #####
	eqtl = []
	if os.path.isfile(args.indir+"/eqtl.txt"):
		eqtl = pd.read_table(args.indir+"/eqtl.txt", delim_whitespace=True)
		eqtl = np.array(eqtl)
		eqtl = eqtl[ArrayIn(eqtl[:,3], genes[:,0])]
		### take top 100000 links per chromosome
		if len(eqtl)>0:
			chrcol = np.array([int(x.split(":")[0]) for x in eqtl[:,0]])
			e_chrom = unique(chrcol)
			eqtl_tmp = []
			for c in e_chrom:
				tmp = eqtl[np.where(chrcol==c)]
				if len(tmp)>args.max_N_links:
					tmp = tmp[tmp[:,5].astype(float).argsort()]
					tmp = tmp[0:args.max_N_links]
				if len(eqtl_tmp)==0:
					eqtl_tmp = tmp
				else:
					eqtl_tmp = np.r_[eqtl_tmp, tmp]
			eqtl = eqtl_tmp
	else:
		print("WARNING: "+args.indir+"/eqtl.txt does not exit. "
		+"Circos plot will not include eQTLs. "
		+"Please download from FUMA to include them.")

	##### get chromosomes to process #####
	if args.chrom is None:
		chrom = unique(loci[:,2])
	else:
		chrom = [int(x) for x in args.chrom.split(",")]
		ci = ci[ArrayIn(ci[:,1].astype(int), chrom)]
		eqtl = eqtl[ArrayIn([int(x.split(":")[0]) for x in eqtl[:,0]], chrom)]
		genes = genes[ArrayIn(genes[:,2].astype(int), chrom)]
		loci = loci[ArrayIn(loci[:,2], chrom)]

	##### get all SNPs #####
	allsnps = []
	if args.sumstats is not None:
		if ".gz" in args.sumstats:
			header = subprocess.check_output("gzip -cd "+args.sumstats+" | head -1", shell=True)
			header = header.strip().split(args.delim)
		else:
			with open(args.sumstats, 'r') as fin:
				header = fin.readline()
				header = header.split(args.delim)
		col_idx = [None,None,None]
		for i in range(0, len(header)):
			if header[i]==args.chrcol:
				col_idx[0] = i
			elif header[i]==args.poscol:
				col_idx[1] = i
			elif header[i]==args.pcol:
				col_idx[2] = i
		if None in col_idx:
			sys.exit("ERROR: Either chromosome, position or P-value column was not found in the sumstats file. "
			+"Please check column names and delimiter are correct. "
			+"Note that column names are case sensitive.")
		for chunk in pd.read_table(args.sumstats, sep=args.delim, dtype=str, chunksize=100*10**6):
			chunk = np.array(chunk)
			if len(allsnps)==0:
				allsnps = chunk[:,col_idx]
				allsnps = allsnps[allsnps[:,2].astype(float)<args.max_snp_p]
				if len(allsnps)>0:
					allsnps = allsnps[ArrayIn(allsnps[:,0].astype(int), chrom)]
			else:
				tmp = chunk[:,col_idx]
				tmp = tmp[tmp[:,2].astype(float)<args.max_snp_p]
				if len(tmp)>0:
					tmp = tmp[ArrayIn(tmp[:,0].astype(int), chrom)]
				if len(tmp)>0:
					allsnps = np.r_[allsnps, tmp]

	##### process per chromosome #####
	snpsout = []
	regions = []
	for c in chrom:
		tmp_genes = genes[genes[:,2]==c]
		tmp_genes = tmp_genes[:,[2,3,4,1,geneshead.index("GenomicLocus")]]
		tmp_genes[:,4] = [int(x.split(":")[-1]) for x in tmp_genes[:,4].astype(str)]
		[tmp_snps, tmp_regions] = createConfig(args, c, loci[loci[:,2].astype(int)==c], ci[np.where((ci[:,1]==c) & (ci[:,4]==c))], snps[snps[:,0]==c], allsnps, tmp_genes)
		if len(snpsout)==0:
			snpsout = tmp_snps
			regions = tmp_regions
		else:
			snpsout = np.r_[snpsout, tmp_snps]
			regions = np.r_[regions, tmp_regions]

	##### write SNPs #####
	with open(args.outdir+"/circos_snps.txt", 'w') as o:
		np.savetxt(o, snpsout, delimiter=" ", fmt="%s")

	##### write regions #####
	regions[:,0] = regions[:,0].astype(str)
	c = ["hs"+str(x) for x in regions[:,0]]
	regions = np.c_[c, regions[:,[1,2]]]
	with open(args.outdir+"/circos_regions.txt", "w") as o:
		np.savetxt(o, regions, delimiter=" ", fmt="%s")

	##### write 3d genome links #####
	ci[:,1] = ["hs"+str(x) for x in ci[:,1]]
	ci[:,4] = ["hs"+str(x) for x in ci[:,4]]
	ci = ci[:,1:7]
	with open(args.outdir+"/ci_links.txt", "w") as o:
		np.savetxt(o, ci, delimiter=" ", fmt="%s")

	##### eqtl write out #####
	if len(eqtl) >0 :
		c = ["hs"+x.split(":")[0] for x in eqtl[:,0]]
		pos = [int(x.split(":")[1]) for x in eqtl[:,0]]
		gstart = list(map(lambda x: genes[genes[:,0]==x,3], eqtl[:,3]))
		gend = list(map(lambda x: genes[genes[:,0]==x,4], eqtl[:,3]))
		eqtl = np.c_[c, pos, [x+1 for x in pos], c, gstart, gend]
	with open(args.outdir+"/eqtl_links.txt", "w") as o:
		np.savetxt(o, eqtl, delimiter=" ", fmt="%s")

	### write genes
	if "eqtlMapSNPs" in geneshead:
		genes = genes[np.where((genes[:,geneshead.index("ciMap")]=="Yes") | (genes[:,geneshead.index("eqtlMapSNPs")]>0))]
	else:
		genes = genes[genes[:, geneshead.index("ciMap")]=="Yes"]
	gid = []
	if len(genes) > 0:
		gid = np.array(["id=0"]*len(genes))
		gid[np.where(genes[:, geneshead.index("ciMap")]=="Yes")] = "id=1"
		if "eqtlMapSNPs" in geneshead:
			gid[np.where((genes[:,geneshead.index("ciMap")]=="Yes") & (genes[:,geneshead.index("eqtlMapSNPs")]>0))] = "id=2"
		genes = np.c_[genes[:,[2,3,4,1]], gid]
		genes[:,0] = ["hs"+str(x) for x in genes[:,0]]
	with open(args.outdir+"/circos_genes.txt", "w") as o:
		np.savetxt(o, genes, delimiter=" ", fmt="%s")

	### write loci and rsID
	tmp = np.c_[["hs"+str(x) for x in loci[:,2]], loci[:,[4,5]]]
	with open(args.outdir+"/highlights.txt", "w") as o:
		np.savetxt(o, tmp, delimiter=" ", fmt="%s")
	tmp = np.c_[["hs"+str(x) for x in loci[:,2]], loci[:,3], [x+1 for x in loci[:,3]], loci[:,1]]
	with open(args.outdir+"/circos_rsID.txt", "w") as o:
		np.savetxt(o, tmp, delimiter=" ", fmt="%s")

	print "total time: "+str(time.time()-start_time)

if __name__=="__main__": main(parser.parse_args())
