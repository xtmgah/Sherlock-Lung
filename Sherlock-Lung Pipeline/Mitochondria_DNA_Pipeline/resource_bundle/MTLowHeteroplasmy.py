import sys

input=open(sys.argv[1],"r")
header=[]
variants=[]
for line in input.readlines():
    if line[0]=="#":
        header.append(line)
    else:
        variants.append(line)
threshold=float(sys.argv[2])
min_vars=int(sys.argv[3])
prefix=sys.argv[4]
outfile=open(prefix+"-low-het-filtered.vcf","w")
low_hets=0
tumor_name=""
tumor_column=-1

for line in header:
    outfile.write(line)
    if line.find("##tumor_sample") > -1:
        tumor_name=line.split("=")[1].strip()
    if line.find("#CHROM") > -1:
        if line.split("\t")[9].strip()==tumor_name:
            tumor_column=9
        elif line.split("\t")[10].strip()==tumor_name:
            tumor_column=10
        else:
            sys.exit("Tumor sample not in file")
    if line.find("##FILTER=<ID=PASS") > -1:
        outfile.write('##FILTER=<ID=many_low_het,Description="Low heteroplasmy variant in sample with too many low heteroplasmy variants (low heteroplasmy threshold=' + str(threshold) + ', required low het variants=' + str(min_vars) + ')">\n')

for line in variants:
    if line.split("\t")[6]=="PASS" or line.split("\t")[6]==".":
        tumor_vaf = line.split("\t")[tumor_column]
        AD = tumor_vaf.split(":")[1]
        ref_allele = int(AD.split(",")[0])
        alt_allele = int(AD.split(",")[1])
        if (float(alt_allele)/float(ref_allele+alt_allele))<threshold:
            low_hets+=1
    if low_hets >= min_vars:
        break

if low_hets>=min_vars:
    for line in variants:
        tumor_vaf = line.split("\t")[tumor_column]
        AD = tumor_vaf.split(":")[1]
        ref_allele = int(AD.split(",")[0])
        alt_allele = int(AD.split(",")[1])
        if (float(alt_allele)/float(ref_allele+alt_allele))>=threshold:
            outfile.write(line)
        elif line.split("\t")[6]=="PASS" or line.split("\t")[6]==".":
            splitLine = line.split("\t")
            splitLine[6] = "many_low_het"
            outfile.write("\t".join(splitLine))
        else:
            splitLine = line.split("\t")
            splitLine[6] = "many_low_het;" + splitLine[6]
            outfile.write("\t".join(splitLine))
else:
    for line in variants:
        outfile.write(line)

input.close()
outfile.close()
