#!/usr/bin/python3
import csv
import traceback,os
from zipfile import ZipFile
import io
import argparse
from datetime import datetime
from time import time
import gzip
def higher(x,y):
    if x>y:
        return x
    elif x<y:
        return y
#METADATA
VERSION = "5.1"

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("inp",help="Path to the input file(vcf)")
parser.add_argument("outp",help="Path to the output file(csv)")
parser.add_argument("col1",help="Column names for set1")
parser.add_argument("col2",help="Column names for set2")
parser.add_argument("-l","--log_location",help = "Preferential log location/directory")
parser.add_argument("-s","--strict",action="store_true")
args = parser.parse_args()
args = vars(args)


readfile = args["inp"]
writefile = args["outp"]
contrastcolumns_1 = []
contrastcolumns_2 = []
optionallogpath= args["log_location"]
strict = args["strict"]
print(f"Starting vcftoallelecount version {VERSION}")

if readfile[-4:]==".zip": # whether we should first unzip it
    fi = ZipFile(readfile)
    print("Ooops, looks like your file is zipped, please enter the name of the vcf file inside:")
    fi.printdir()
    file_in_zip = input()
    main_vcf = io.TextIOWrapper(fi.open(file_in_zip,"r"))
elif readfile[-3:]==".gz":
    print("Unzippping gunzip file")
    main_vcf = gzip.open(readfile,mode="rt")
else:
    main_vcf = open(readfile,"r",encoding='latin-1') # main data vcf file
writeto = open(writefile,"w")
output_bool_table_file = csv.writer(writeto,delimiter="\t")
file_name = writefile.split("/")[-1] # name of the file
pth = f"log-{file_name}.txt"
if optionallogpath==None: # in case there is no optionallogpath
    vers = 0
    while os.path.exists(pth): # skip already written logs
        vers += 1
        pth = f"log-{file_name}({vers}).txt"
else:
    vers = 0
    while os.path.exists(optionallogpath + "/" + pth): # skip already written logs
        vers += 1
        pth = f"log-{file_name}({vers}).txt"
    pth = optionallogpath + "/" + pth

log = open(pth,"w") # metadata and run data
log.writelines(f"Starting program... version = {VERSION}\n")
log.writelines(f"input file = {readfile}, input_file_size = {round(os.path.getsize(readfile)/1000000,2)}MB\n")
log.writelines(f"output_file = {writefile}\n")
log.writelines(f"columns of set1 = {args['col1'].split(',')}\n")
log.writelines(f"columns of set2 = {args['col2'].split(',')}\n")
log.writelines(f"Date of the run: {datetime.now()}\n")
print(f"Starting program... input file = {readfile}, input_file_size = {round(os.path.getsize(readfile)/1000000,2)}MB, Strict mode={strict}")
print(f"output_file = {writefile}")
# statistics and logging
not_1_on_1 = 0  # A on AAC and sim.
total_diag = 0  # number of diagnostic reads
position_not_readable = 0 # number of N nucleotides
linecount = 0 # number of lines in total
secs = int(time()) # variable to derive run time later

ploid = None # is it haploid, triploid or more????

while True: # getting to the start of the file
    try:
        read = main_vcf.readline()
    except UnicodeDecodeError as e:
        raise ValueError("Seems like there was a decoding error, are you using the correct encoding?")
    if read[0:2]!="##":
        break
headers = read.split("\t")
headers[-1] = headers[-1][:-1] # removing the last character from the last header - it is a \n and would cause issues for column assignment
ploidity_table = {}

cols1 : list = args["col1"].split(",")
cols2 : list = args["col2"].split(",")

for n,col in enumerate(headers):
    if col in cols1:
        contrastcolumns_1.append(n)
        cols1.remove(col)
    elif col in cols2:
        contrastcolumns_2.append(n)
        cols2.remove(col)

if len(cols1)>0 or len(cols2)>0:
    raise ValueError(f"Please correct your column names, the following were not found: {[*cols1,*cols2]} ")

def nextline()->list: # returns next line of reads as a list
    return main_vcf.readline().split("\t")
def convert_to_map(lst:list)->dict: # converts read line to a map
    map = {}
    for i in range(len(headers)):
        map[headers[i]]=lst[i]
    return map
def single_certain(nuc1:list,nuc2:list)->int: # if no alleles are shared within the two species for given read, return 1 - contrast is undeniably different
    # nuc1 = [1,0,0,1] ex. value
    # is it diagnostic?
    bol = 1
    for i in range(len(nuc1)):
        if nuc1[i]*nuc2[i]!=0:
            bol = 0
            break
    if bol:
        global total_diag
        total_diag+=1
    return bol
def posread(set,col,line,alt_ref_list): # counting nucleotides in given columns for given animal
    found = False
    setslice = []
    for i in col:
        setslice.append((line[i],i)) # choosing the right columns alt_ref_list on previously selected ones  
    for sample_index in setslice: # all column values of the specific set
        sample = sample_index[0]
        col_index = sample_index[1]
        if ploidity_table[col_index]>2: # quick fix should not be the way to do this. TODO
            continue
        snps_list = sample[:ploidity_table[col_index]*2:2]# haplotypes of this sample
        for x in snps_list:
            if x==".":
                continue
            set[alt_ref_list[int(x)]]=1 #  adding to the set table
            found=True
    return found

output_bool_table_file.writerow(["CHROM","POS","set1_A","set1_C","set1_G","set1_T","set2_A","set2_C","set2_G","set2_T","diag"]) # headers in output
# the first time we are kind of doing a do while loop here
line = nextline()
ploid = len(line[contrastcolumns_1[0]].replace("|","/").split("/"))# deciding whether it's haploid or not
for x in contrastcolumns_1:
    ploidity_table[x]=len(line[x].replace("|","/").split("/")) # written like 4
for x in contrastcolumns_2:
    ploidity_table[x]=len(line[x].replace("|","/").split("/"))


while True:

    try:
        if len(line)==1:
            break
        linecount+=1
        current_map = convert_to_map(line)

        if len(current_map["REF"])>1 or (len(current_map["ALT"])>1 and (current_map["ALT"][1]!="," or len(current_map["ALT"])!=3)): # skipping those not 1 on 1 -> but not A G,C, just A TTA
            not_1_on_1+=1
        elif current_map["ALT"]=="N" or current_map["REF"]=="N" or "*" in current_map["ALT"].split(","): # skipping unreadable
            position_not_readable+=1
        else:
            alt_ref_list = [current_map["REF"],*current_map["ALT"].split(",")]
            set2 = {"A":0,"C":0,"G":0,"T":0} # tables of certainty - one for each set
            set1 = {"A":0,"C":0,"G":0,"T":0}
            posread(set1,contrastcolumns_1,line,alt_ref_list) # reading the columns and adding to set tables
            posread(set2,contrastcolumns_2,line,alt_ref_list)
            s1 = list(set1.values()) # the bool values of nucleotides for each set derived from tables of certainty
            s2 = list(set2.values())
            diag = single_certain(s1,s2) # determining if the nucleotide is diagnostic
            if diag==True or strict==False:
                output_bool_table_file.writerow([current_map["#CHROM"],current_map["POS"],*s1,*s2,diag,])

        line = nextline()
    except:
        log.writelines(f"ERROR: CRASHED on line {linecount} : {line}\n",)
        log.writelines(f"Exception: {traceback.format_exc()}\n")
        print("Program crashed, oh no... ", traceback.format_exc())
        break



writeto.close()
print("bool writing successful")
log.writelines(f"program was running for {int(time())-secs} seconds\n")
log.writelines(f"with parameters strict={strict}, optional_log_location={optionallogpath}\n")
log.writelines(f"It skipped {not_1_on_1} reads, because their contrast was not 1 on 1 => A on ACCT or similar\n")
log.writelines(f"{total_diag}: Number of diagnostic reads \n")
log.write(f"It converted {linecount-not_1_on_1-position_not_readable} reads in total\n")
log.write(f"{position_not_readable} lines had an unreadable nucleotide\n")

