#!/usr/bin/env python

import sys
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
import toolz
import time
import os
import requests
from glob import glob
from bs4 import BeautifulSoup
from gpcr_variables import COCOMAPS_SUBMIT, COCOMAPS_OUTPUT, \
							INTERPROSURF, DEFAULT_EMAIL, STRUCTURAL_FEATURES, \
							INTERPROSURF_START, RESULTS_FOLDER, COCOMAPS_START, \
							STRUCTURAL_PDBS_FOLDER, DEFAULT_FOLDER
requests.packages.urllib3.disable_warnings()

__author__ = "J.G. Almeida & A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

if 'darwin' in sys.platform:
	sys_sep = '/'
	CHROMEDRIVER_PATH = DEFAULT_FOLDER + sys_sep + "chromedriver"
elif 'win' in sys.platform:
	sys_sep = '\\'
	CHROMEDRIVER_PATH = DEFAULT_FOLDER + sys_sep + "chromedriver.exe"
else:
	sys_sep = '/'
	CHROMEDRIVER_PATH = DEFAULT_FOLDER + sys_sep + "chromedriver"

cocomaps_url = COCOMAPS_SUBMIT
cocomaps_url_finished = COCOMAPS_OUTPUT
interprosurf_url = INTERPROSURF

def open_lines(file):

	"""
	Open and read an input file
	"""
	o = open(file)
	lines = o.readlines()
	o.close()
	return lines

def open_write(file_name, data):

	"""
	Open an output file and write the data into it
	"""
	o = open(file_name,'w')
	o.write(data)
	o.close()

def get_atoms(lines):

	"""
	Isolate the atoms in the pdb file
	"""
	atom_lines = []
	for line in lines:
		if line[0:4] == 'ATOM':
			atom_lines.append(line)
	return atom_lines

def get_chains(atom_lines):

	"""
	Identify the individual chains in the pdb
	"""
	chains = set()
	for line in atom_lines:
		chain = line[21]
		chains.add(chain)
	chains = list(chains)
	return sorted(chains)

def cocomaps_digger(pdb_file,identifier,download_path, sleep_time = 10):

	"""
	Retrieve CoCoMaps data with selenium
	"""
	print("CoCoMaps data is being retrieved...")
	atoms = toolz.pipe(pdb_file,open_lines,get_atoms)
	chains = get_chains(atoms)
	
	params = {'pdb_file': pdb_file,'local_chain_molecule1': chains[0],\
					'local_chain_molecule2': chains[1], 'cutoff': 5}
	new_names = ["interaction_overview","minimum_distances","hbonds","asa_all","asa_mol1","asa_mol2"]
	options = Options()
	options.add_argument('--headless')
	options.add_argument('--disable-gpu')
	driver = webdriver.Chrome(CHROMEDRIVER_PATH, chrome_options = options)
	driver.get(cocomaps_url)
	assert "CoCoMaps" in driver.title
	driver.find_elements_by_name("type_pdb")[1].click()
	for param in params:
		pdb = driver.find_element_by_name(param)
		pdb.send_keys(params[param])
	driver.find_element_by_name("submit").click()
	time.sleep(sleep_time)
	download_list = driver.find_elements_by_class_name("downloadTable")[3:]

	for i in range(0,len(download_list)):
		curr_download = download_list[i]
		url_name = curr_download.get_attribute("href")
		new_name = download_path + sys_sep + COCOMAPS_START + '_' + '_'.join([identifier,new_names[i]]) + '.txt'
		r = requests.get(url_name,verify = False)
		with open(new_name, "wb") as table:
			table.write(r.content)

	dist_element = driver.find_element_by_xpath('//a[@target ="_table"]')
	dist_element.click()
	time.sleep(sleep_time)
	while len(driver.window_handles) == 1:
		time.sleep(time_sleep)
	driver.switch_to_window(driver.window_handles[1])
	time.sleep(sleep_time)
	download_dist = driver.find_element_by_class_name("downloadTable")
	table_url = download_dist.get_attribute("href")
	r = requests.get(table_url,verify = False)
	table_name = download_path + sys_sep + COCOMAPS_START + '_' + identifier + "_dist_table.txt"
	with open(table_name, "wb") as table:
		table.write(r.content)
	driver.close()
	driver.switch_to_window(driver.window_handles[0])
	

def interprosurf_digger(pdb_file,identifier,download_path):

	"""
	Retrieve pdb data
	"""
	print("InterProSurf data is being retrieved...")
	complex_col = STRUCTURAL_FEATURES
	all_data = []
	complex_data = []
	residue_data = []
	chains = set()
	options = Options()
	options.add_argument('--headless')
	options.add_argument('--disable-gpu')
	driver = webdriver.Chrome(CHROMEDRIVER_PATH, chrome_options=options)
	driver.get(interprosurf_url)
	assert "InterProSurf" in driver.title

	name = driver.find_element_by_name("name")
	email = driver.find_element_by_name("email")
	file = driver.find_element_by_name("PDBfile")
	
	name.send_keys(identifier)
	email.send_keys(email_address)
	file.send_keys(pdb_file)
	driver.find_element_by_xpath("//input[@type='submit']").click()
	time.sleep(5)
	driver.switch_to_window(driver.window_handles[1])
	all_html = driver.page_source
	soup = BeautifulSoup(all_html,"html.parser")
	for mess in soup.find_all('table'):
		bowl = mess.get_text()
		data = bowl.split('\n')
		data = [x.strip() for x in data if x.strip() != '']
		all_data.append(data)
	for data in all_data[4:]:
		if len(data[-2]) == 1: 
			chains.add(data[-2])
			residue_data.append(data)
	chains = sorted(list(chains))
	nchains = len(chains)
	for i in range(0,nchains):
		new_data = [chains[i]] + [x.split('=')[-1].strip() for x in all_data[i][1:]]
		complex_data.append(new_data)
	complex_full = [''.join(chains)] + [x.split('=')[-1].strip() for x in all_data[nchains][1:]]
	complex_data.append(complex_full)
	residue_col = all_data[nchains + 1]
	residue_name = download_path + sys_sep + INTERPROSURF_START + '_' + identifier + '_residues.csv'
	complex_name = download_path + sys_sep + INTERPROSURF_START + '_' + identifier + '_complex.csv'
	residue_write = ','.join(residue_col) + '\n' + '\n'.join([','.join(data) for data in residue_data])
	complex_write = ','.join(complex_col) + '\n' + '\n'.join([','.join(data) for data in complex_data])
	open_write(residue_name,residue_write)
	open_write(complex_name,complex_write)
	driver.close()
	driver.switch_to_window(driver.window_handles[0])

def parse_pdb_folder(folder):

	"""
	Identify all *pdb files in folder
	"""
	folder = folder.rstrip(sys_sep)
	reg = folder + sys_sep + '*.pdb'
	all_pdbs = glob(reg)
	return sorted(all_pdbs)

def wraper(in_folder,download_path):

	"""
	Initialize chromedriver instance, make sure the adequate chromedriver is in place
	"""
	all_pdbs = parse_pdb_folder(DEFAULT_FOLDER + sys_sep + STRUCTURAL_PDBS_FOLDER)
	options = Options()
	options.add_argument('--headless')
	options.add_argument('--disable-gpu')
	driver = webdriver.Chrome(CHROMEDRIVER_PATH, chrome_options = options)
	if os.path.isdir(download_path) == False:
		sys.exit("Input folder does not exist.")
	if os.path.isdir(download_path) == False:
		os.makedirs(download_path)
		print("Output folder will be created:",download_path)
	if len(all_pdbs) == 0:
		sys.exit("No PDBs in the supplied input folder.")

	for pdb in all_pdbs:
		pre_identifier = pdb.split(sys_sep)[-1]
		identifier = pre_identifier[:-4]
		print("Current PDB:",identifier)
		interprosurf_digger(pdb,identifier,download_path)
		success, current_sleep_time, output = False, 5, None
		output = cocomaps_digger(pdb,identifier, download_path, sleep_time = current_sleep_time)
		print("Data for",identifier,"has been retrieved!")
	driver.quit()

email_address = DEFAULT_EMAIL
in_folder = os.getcwd()
out_folder = in_folder + '/' + RESULTS_FOLDER
wraper(in_folder, out_folder)
