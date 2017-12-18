# file format: date_chirp_virus_cells_techrep#_biorep#.xlsx
# csv file: All technical replicates should go on the same, unique line with comma-separation. 

from xlrd import open_workbook
import sys
import csv
import os
import argparse
from numpy import std
from xlwt import Workbook

def parseArguments():
	parser = argparse.ArgumentParser(description = "MS analysis")
	parser.add_argument('--csv', dest = "csv", type = str, nargs = '?', required = True, help = "Input csv file")
	parser.add_argument('--mc', action = "store_true", help = "Include if multiple chirp")
	parser.add_argument('--tr', dest = "tr", type = int, nargs = '?', required = True, help = "Number of tech reps")
	parser.add_argument('--br', dest = "br", type = int, nargs = '?', required = True, help = "Number of bio reps")
	parser.add_argument('--cr', dest = "cr", type = float, nargs = '?', default = 1.0, help = "Correction; number of peptides to be added to avoid 0; default is 1")
	parser.add_argument('--mp', dest = "mp", type = int, nargs = '?', default = 1, help =  "Minimum number of peptides per Bio rep")
	parser.add_argument('--mtr', dest = "mtr", type = int, nargs = '?', default = 1, help = "Minimum number of tech reps")
	parser.add_argument('--mbr', dest = "mbr", type = int, nargs = '?', default = -1, help = "Minimum number of bio reps")
	parser.add_argument('--s', dest = "s", type = float, nargs = '?', default = 100.0, help = "Minimum score")
	parser.add_argument('--lp', dest = "lp", type = float, nargs = '?', default = 5.0, help = "Minimum log prob")
	args = parser.parse_args()	
	return args.csv, args.mc, args.tr, args.br, args.cr, args.mp, args.mtr, args.mbr, args.s, args.lp

def extractName(protein_info, protein_info_list):
	try:
		if protein_info.find("GN") != -1:
			prot_id = protein_info_list[2].split(" ")[0]
		else:
			prot_id = protein_info_list[0][1: len(protein_info_list[0])]
	except:
		prot_id = protein_info_list[0][1: len(protein_info_list[0])]

	return prot_id

def combineTechFiles(tech_reps, num_tech_reps, num_bio_reps, num_comb_tech_exp, file_prefix, min_tech, min_score, min_log):
	# Function that combines data across technical replicates and creates a new spreadsheet for each set of technical replicates
	# comb_tech_exp[i][prot_id][j] entries:
	# j = 0: UniProt ID
	# j = 1: Spectra List
	# j = 2: Number of Amino Acids
	# j = 3: Normalized Spectra List
	# j = 4: Summed Spectra
	# j = 5: Summed Normalized Spectra
	# j = 6: Comma Separated Spectra List

	comb_tech_exp = [{} for i in range(num_comb_tech_exp)]
	comb_tech_files = []
	for i in range(num_comb_tech_exp):
		for j in range(num_tech_reps):
			prot_name_set = set()
			prot_UniID = {}
			prot_count = {}
			prot_num_aa = {}

			try:
				cur_wb = open_workbook(tech_reps[i][j])
			except:
				print "Can't open file " + tech_reps[i][j]
				return comb_tech_exp, False
			sheet_spectra = cur_wb.sheet_by_index(1)
			sheet_prot = cur_wb.sheet_by_index(2)

			prev_rank = -1
			for row_index in range(1, sheet_prot.nrows):
				cur_rank = int(sheet_prot.cell(row_index, 0).value)
				if cur_rank == prev_rank:
					continue
				else:
					prev_rank = cur_rank
				protein_info = str(sheet_prot.cell(row_index, 1).value)
				if protein_info.find("Reverse") != -1:
					continue
				protein_info_list = protein_info.split("=")
				prot_id = extractName(protein_info, protein_info_list)
				num_aa = float(sheet_prot.cell(row_index, 10).value)
				prot_num_aa[prot_id] = num_aa

			for row_index in range(1, sheet_spectra.nrows):
				if sheet_spectra.cell(row_index, 13).value >= min_score and sheet_spectra.cell(row_index, 16).value >= min_log:
					protein_info = str(sheet_spectra.cell(row_index, 18).value)
					if protein_info.find("Reverse") != -1:
							continue
					protein_info_list = protein_info.split("=")
					prot_id = extractName(protein_info, protein_info_list)
					try:
						prot_UniID[prot_id] = protein_info.split("|")[1]
					except: 
						prot_UniID[prot_id] = protein_info_list[0][1: len(protein_info_list[0])]
						
					if prot_id in prot_name_set:
						prot_count[prot_id] += 1.0
					else:
						prot_name_set.add(prot_id)
						prot_count[prot_id] = 1.0

			prot_count_keys = sorted(prot_count.keys(), key = lambda x: prot_count[x])
			for prot_id in prot_count:
				if prot_id not in prot_num_aa:
					prot_count.pop(prot_id)
				else:
					spectra = prot_count[prot_id]
					if prot_id in comb_tech_exp[i].keys():
						comb_tech_exp[i][prot_id][1].append(spectra)
					else:
						comb_tech_exp[i][prot_id] = []
						num_aa = prot_num_aa[prot_id]
						uniprot = prot_UniID[prot_id]
						comb_tech_exp[i][prot_id].append(uniprot)
						comb_tech_exp[i][prot_id].append([spectra])
						comb_tech_exp[i][prot_id].append(num_aa)
						comb_tech_exp[i][prot_id].append([])

			total_spectra = sum(prot_count.values())
			for prot_id in prot_count_keys:
				comb_tech_exp[i][prot_id][3].append(prot_count[prot_id] / total_spectra * 1000.0)

		for prot_id in comb_tech_exp[i]:
			if len(comb_tech_exp[i][prot_id][1]) >= min_tech:
				num_zeroes = num_tech_reps - len(comb_tech_exp[i][prot_id][1])
				while num_zeroes > 0:
					comb_tech_exp[i][prot_id][1].append(0.0)
					comb_tech_exp[i][prot_id][3].append(0.0)
					num_zeroes -= 1
			else:
				comb_tech_exp[i].pop(prot_id)

	for i in range(num_comb_tech_exp):
		name = csv_prefix + "_" + file_prefix[i] + "_%s" % (str(i % num_bio_reps + 1)) + "_tech" + ".xls"
		wb = Workbook()
		sheet1 = wb.add_sheet('Sheet 1')
		count = 1
		sheet1.write(0, 1, file_prefix[i] + ": ID")
		sheet1.write(0, 2, file_prefix[i] + ": UniProtID")
		sheet1.write(0, 3, file_prefix[i] + ": Spectras")
		sheet1.write(0, 4, file_prefix[i] + ": Sum_Spectra")
		sheet1.write(0, 5, file_prefix[i] + ": Normalized Spectras")
		sheet1.write(0, 6, file_prefix[i] + ": Sum Normalized Spectras")
		sheet1.write(0, 7, file_prefix[i] + ": num AAs")
		for prot_id in comb_tech_exp[i]:
			sum_spectra = sum(comb_tech_exp[i][prot_id][1])
			sum_normalized = sum(comb_tech_exp[i][prot_id][3])
			comb_tech_exp[i][prot_id].append(sum_spectra)
			comb_tech_exp[i][prot_id].append(sum_normalized)
			spectras = ", ".join(str(x) for x in comb_tech_exp[i][prot_id][1])
			norm_spectras = ", ".join(str(x) for x in comb_tech_exp[i][prot_id][3])
			comb_tech_exp[i][prot_id].append(spectras)
			sheet1.write(count, 0, count)
			sheet1.write(count, 1, prot_id)
			sheet1.write(count, 2, comb_tech_exp[i][prot_id][0])
			sheet1.write(count, 3, spectras)
			sheet1.write(count, 4, sum_spectra)
			sheet1.write(count, 5, norm_spectras)
			sheet1.write(count, 6, sum_normalized)
			sheet1.write(count, 7, comb_tech_exp[i][prot_id][2])
			count += 1
		wb.save(name)

	return comb_tech_exp, True 

def combineBioFiles(comb_tech_exp,num_bio_reps, num_comb_tech_exp, file_prefix, min_pep, min_bio):
	# Function that combines data across biological replicates and outputs Excel files for each set of biological replicates
	# comb_tech_exp[i][prot_id][j] entries:
	# j = 0: UniProt ID
	# j = 1: Comma Separated Spectra List
	# j = 2: Summed Normalized Spectra List (summed/normalized across technical)
	# j = 3: Average of Summed Normalized Spectra (summed/normalized across technical and then averaged across biological)
	# j = 4: Summed Spectra List (summed across technical)
	# j = 5: Number of Amino Acids
	# j = 6: Standard deviation of Summed Normalized Spectra
	# j = 7: Comma Separated Summed Normalized Spectra

	num_comb_bio_exp = num_comb_tech_exp / num_bio_reps
	comb_bio_files = []
	comb_bio_exp = [{} for i in range(num_comb_bio_exp)]

	start_index = 0
	end_index = num_bio_reps
	for i in range(num_comb_bio_exp):
		while start_index < end_index:
			for prot_id in comb_tech_exp[start_index]:
				if prot_id in comb_bio_exp[i].keys():
					comb_bio_exp[i][prot_id][1] += "," + comb_tech_exp[start_index][prot_id][6]
					comb_bio_exp[i][prot_id][2].append(comb_tech_exp[start_index][prot_id][5])
					comb_bio_exp[i][prot_id][3] += comb_tech_exp[start_index][prot_id][5]
					comb_bio_exp[i][prot_id][4].append(comb_tech_exp[start_index][prot_id][4])
				else:
					comb_bio_exp[i][prot_id] = []
					comb_bio_exp[i][prot_id].append(comb_tech_exp[start_index][prot_id][0])
					comb_bio_exp[i][prot_id].append(comb_tech_exp[start_index][prot_id][6])
					comb_bio_exp[i][prot_id].append([comb_tech_exp[start_index][prot_id][5]])
					comb_bio_exp[i][prot_id].append(comb_tech_exp[start_index][prot_id][5])
					comb_bio_exp[i][prot_id].append([comb_tech_exp[start_index][prot_id][4]])
					comb_bio_exp[i][prot_id].append(comb_tech_exp[start_index][prot_id][2])
			start_index += 1
		end_index += num_bio_reps

		for prot_id in list(comb_bio_exp[i]):
			num_above_cutoff = 0
			for j in range(len(comb_bio_exp[i][prot_id][4])):
				if comb_bio_exp[i][prot_id][4][j] >= min_pep:
					num_above_cutoff += 1
			if num_above_cutoff >= min_bio:
				num_zeroes = num_bio_reps - len(comb_bio_exp[i][prot_id][2])
				while num_zeroes > 0:
					comb_bio_exp[i][prot_id][2].append(0.0)
					num_zeroes -= 1
				comb_bio_exp[i][prot_id][3] = comb_bio_exp[i][prot_id][3] / num_bio_reps
				comb_bio_exp[i][prot_id].append(std(comb_bio_exp[i][prot_id][2]))
				norm_spectras = ", ".join(str(round(x, 5)) for x in comb_bio_exp[i][prot_id][2])
				comb_bio_exp[i][prot_id].append(norm_spectras)
			else:
				comb_bio_exp[i].pop(prot_id)

	cur_index = 0
	for i in range(num_comb_bio_exp):
		name = csv_prefix + "_" + file_prefix[i * num_bio_reps] + "_bio" + ".xls"
		comb_bio_files.append(name)
		cur_index += num_bio_reps
		if num_bio_reps <= 1:
			continue
		wb = Workbook()
		sheet1 = wb.add_sheet('Sheet 1')
		count = 1
		sheet1.write(0, 1, file_prefix[i * num_bio_reps] + ": ID")
		sheet1.write(0, 2, file_prefix[i * num_bio_reps] + ": UniProtID")
		sheet1.write(0, 3, file_prefix[i * num_bio_reps] + ": Average spectra")
		sheet1.write(0, 4, file_prefix[i * num_bio_reps] + ": Standard deviation")
		sheet1.write(0, 5, file_prefix[i * num_bio_reps] + ": Normalized spectras")
		sheet1.write(0, 6, file_prefix[i * num_bio_reps] + ": num AAs")
		for prot_id in comb_bio_exp[i]:
			sheet1.write(count, 0, count)
			sheet1.write(count, 1, prot_id)
			sheet1.write(count, 2, comb_bio_exp[i][prot_id][0])
			sheet1.write(count, 3, round(comb_bio_exp[i][prot_id][3], 5))
			sheet1.write(count, 4, round(comb_bio_exp[i][prot_id][6], 5))
			sheet1.write(count, 5, comb_bio_exp[i][prot_id][7])
			sheet1.write(count, 6, comb_bio_exp[i][prot_id][5])
			count += 1
		wb.save(name)

	return num_comb_bio_exp, comb_bio_files, comb_bio_exp

def multipleChirp(num_tech_reps, num_bio_reps, num_comb_bio_exp, comb_bio_files, comb_bio_exp, file_prefix, correction):
	# Function that creates a spreadsheet to compare data from different experiments
	
	multiple_chirp_comb = {}
	for i in range(num_comb_bio_exp):
		for prot_id in comb_bio_exp[i]:
			if prot_id not in multiple_chirp_comb:
				multiple_chirp_comb[prot_id] = [0 for x in range(3 * num_comb_bio_exp + 2)]
				multiple_chirp_comb[prot_id][3 * num_comb_bio_exp] = comb_bio_exp[i][prot_id][0]
				multiple_chirp_comb[prot_id][3 * num_comb_bio_exp + 1] = comb_bio_exp[i][prot_id][5]

	for prot_id in multiple_chirp_comb:
		for i in range(num_comb_bio_exp):
			if prot_id in comb_bio_exp[i]:
				multiple_chirp_comb[prot_id][i] = comb_bio_exp[i][prot_id][1]
				multiple_chirp_comb[prot_id][i + num_comb_bio_exp] = comb_bio_exp[i][prot_id][3] + correction
			else: 
				multiple_chirp_comb[prot_id][i] = ",".join(['0.0' for x in range(num_tech_reps * num_bio_reps)])
				multiple_chirp_comb[prot_id][i + num_comb_bio_exp] =  correction

			multiple_chirp_comb[prot_id][2 * num_comb_bio_exp + i] = (multiple_chirp_comb[prot_id][i + num_comb_bio_exp] / 
				multiple_chirp_comb[prot_id][3 * num_comb_bio_exp])

	wb = Workbook()
	sheet1 = wb.add_sheet('Sheet 1')
	sheet1.write(0, 0, "Protein ID")
	sheet1.write(0, 1, "UniProtID")
	for i in range(num_comb_bio_exp):
		sheet1.write(0, i + 2,  file_prefix[i * num_bio_reps] + ": Spectras")
		sheet1.write(0, i + 2 + num_comb_bio_exp,  file_prefix[i * num_bio_reps] + ": Average Spectra")
		sheet1.write(0, i + 3 + 2 * num_comb_bio_exp, file_prefix[i * num_bio_reps] + ": Average Spectra / numAAs")
	sheet1.write(0, 2 + 2 * num_comb_bio_exp, "numAAs")

	count = 1
	for prot_id in multiple_chirp_comb:
		sheet1.write(count, 0, prot_id)
		sheet1.write(count, 1, multiple_chirp_comb[prot_id][3 * num_comb_bio_exp])
		for i in range(num_comb_bio_exp):
			sheet1.write(count, i + 2, multiple_chirp_comb[prot_id][i])
			sheet1.write(count, i + 2 + num_comb_bio_exp, 
				round(multiple_chirp_comb[prot_id][i + num_comb_bio_exp], 5))
			sheet1.write(count, i + 2 + 2 * num_comb_bio_exp, 
				round(multiple_chirp_comb[prot_id][i + 2 * num_comb_bio_exp],5))
		sheet1.write(count, 2 + 3 * num_comb_bio_exp, multiple_chirp_comb[prot_id][3 * num_comb_bio_exp + 1])
		count += 1
	wb.save(csv_prefix + "_Combined_Experiments.xls")

def main():
	input_files, multiple_chirp, num_tech_reps, num_bio_reps, correction, min_pep, min_tech, min_bio, min_score, min_log = parseArguments()
	
	if min_bio == -1:
		min_bio = num_bio_reps

	if input_files.find("csv") == -1:
		print "First argument is not a csv file."
		return

	global csv_prefix
	csv_prefix = input_files.split(".")[0]

	try:
		sample_sheet = open(input_files, 'r')
	except: 
		print "Not a valid csv file/path"
		return
	else:
		tech_reps = list(csv.reader(sample_sheet))
		num_comb_tech_exp = len(tech_reps)
		for i in range(num_comb_tech_exp):
			for j in range(len(tech_reps[0])):
				if not os.path.exists(tech_reps[i][j]):
					print "%s is not a valid file/path" % tech_reps[i][j]
					return
		file_prefix = []
		for i in range(num_comb_tech_exp):
			prefix = tech_reps[i][0].split('_')
			file_prefix.append(prefix[1] + "_" + prefix[2])
		sample_sheet.close()

	comb_tech_exp, boolean = combineTechFiles(tech_reps, num_tech_reps, num_bio_reps, num_comb_tech_exp, file_prefix, min_tech, min_score, min_log)
	if not boolean:
		return 

	num_comb_bio_exp, comb_bio_files, comb_bio_exp = combineBioFiles(comb_tech_exp, 
			num_bio_reps, num_comb_tech_exp, file_prefix, min_pep, min_bio)

	if multiple_chirp:
		multipleChirp(num_tech_reps, num_bio_reps, num_comb_bio_exp, 
			comb_bio_files, comb_bio_exp, file_prefix, correction)

if __name__ == '__main__':
	main()
