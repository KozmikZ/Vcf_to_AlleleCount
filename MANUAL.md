Version 5 of the program.
This is a program designed to convert large vcf files to csv/tsv files comparing **two** species. It will also create a **log** file, in which it will leave valuable statistics and data about the run. The file gets created in the running directory by default. Check parameters for more.
Given indexes of the two sets columns, the program will create a table that counts alleles on each position with a format like this:

###	CHROM	POS	set1_A	set1_C	set1_G	set1_T	set2_A	set2_C	set2_G	set2_T	diag

A typical csv/tsv would look like this:
CHROM	POS	set1_A	set1_C	set1_G	set1_T	set2_A	set2_C	set2_G	set2_T	diag
cott-wg5-ch1	1501	0	0	0	0	0	0	0	2	1
cott-wg5-ch1	1543	0	0	0	0	2	0	0	0	1
cott-wg5-ch1	4226	0	0	0	0	0	0	0	0	1
cott-wg5-ch1	5305	0	0	0	0	0	0	0	0	1
cott-wg5-ch1	11866	0	2	0	0	0	0	0	2	1
cott-wg5-ch1	11870	0	0	0	0	0	0	0	2	1
cott-wg5-ch1	11880	0	0	0	0	0	0	2	0	1
cott-wg5-ch1	12172	2	0	0	0	0	0	0	0	1
cott-wg5-ch1	12191	0	2	0	0	0	0	0	0	1



# Columns:

**CHROM** (string) - The chromosome this snp was read on.

**POS** (integer) - the coordinate of this snp on that chromosome

**set1_N** (integer) - the specific base and how much of it was found in this set. This counts in **each haplotype** and adds it to the total sum of the alleles found. Meaning that if you have 0/0 in your vcf, it adds 2 to the total sum -> How many haplotypes of set1 had this allele.

**set2_N**(integer) - the specific base and how much of it was found in this set. This counts in **each haplotype** and adds it to the total sum of the alleles found.   Meaning that if you have 0/0 in your vcf, it adds 2 to the total sum -> How many haplotypes of set2 had this allele.

**diag** (boolean) - Was this line *diagnostic*? -> can we explicitly determine that according to this row, bases on this position are different in the two species. example -> a row like **0 0 3 0** , **0 2 0 0**  is diagnostic, because guanine is present only in set1 and cytosine only in set2, therefore the bases are always different.

This program accounts for **multiple alternatives**, but not for snps longer than 1 
So lines like
**REF : A , ALT: C,G** won't get skipped. But lines, such as **REF: TTC, ALT: GTTTTA**, **will** get skipped.


# Parameters:

### Needed:

**inp** - Input file path. The *.vcf* you want to convert

**outp** - Output file path. The *.csv* or *.tsv* you want to write to.

**col1** - Column names of the first species/set.

**col2** - Column names of the second species/set.


### Optional:

**-s(--strict)** - Strict mode -> During strict mode, only lines which are diagnostic get written to the submitted output file. All other lines are skipped. The data type is a boolean. False by default.

**-l(--log_location)** - Optional path to the log file. By default, the file gets created in the running directory.

An example command would be:
./vcftoallelecount.py Paces_test_SNP_Cob.g.vcf Paces_test_SNP_Cob.g.vcf.out.tsv 080118_TT1_sort,C280717_5_TT_sort 201218EE_sort,210617EE1M_sort
