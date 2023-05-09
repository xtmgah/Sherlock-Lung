#!/usr/bin/python
import sys

def main():
    annovar=open(sys.argv[1],"r")
    vcf=open(sys.argv[2],"r")
    prefix=sys.argv[3]
    annotated=[]
    vcf_variants=[]
    outfile=open(prefix+"-annovar-filtered.vcf","w")
    tumor_name=""
    tumor_column=-1

    for line in annovar.readlines():
        annotated.append(line)

    for line in vcf.readlines():
        if line[0]=="#":
            outfile.write(line)
            if line.find("##FILTER=<ID=PASS") > -1:
                outfile.write('##FILTER=<ID=annovar_filtered,Description="MT variant present in population above frequency 0.01">\n')
        else:
            vcf_variants.append(line)

    for line in annotated:
        annotation=line.split("\t")[1]
        AF_hom=float(annotation.split(";AF_hom=")[1].split(";")[0])
        AF_het=float(annotation.split(";AF_het=")[1].split(";")[0])
        if AF_hom+AF_het>=0.01:
            index_of=naive_search(line.split("\t"), vcf_variants)
            if index_of>-1:
                vcf_variants[index_of] = edit_filter(vcf_variants[index_of])
            else:
                print ('An annotated variant:" '+ " ".join(line.split("\t")[2:7]) + ' " is missing from the VCF file; double check file input')

    for line in vcf_variants:
        outfile.write(line)

    outfile.close()

def edit_filter(variant):
    if variant.split("\t")[6] == "PASS" or variant.split("\t")[6] == ".":
        splitLine = variant.split("\t")
        splitLine[6] = "annovar_filtered"
    else:
        splitLine = variant.split("\t")
        splitLine[6] = "annovar_filtered;" + splitLine[6]
    return "\t".join(splitLine)

# Use naive search until complexity becomes a problem or a robust solution to binary search problem is found
def naive_search(annotated_variant, variants):
    annotated_var_position = int(annotated_variant[3])
    for index in range(len(variants)):
        variant_position = int(variants[index].split("\t")[1])
        variant_ref = variants[index].split("\t")[3]
        variant_alt = variants[index].split("\t")[4]
        if annotated_var_position==variant_position and annotated_variant[5]==variant_ref and annotated_variant[6]==variant_alt:
            return index
    return -1

# Abandoned binary search for now; not sure how to handle cases where several variants have same position,
# but different snv/indel
#
#def binary_search(annotated_variant, variants):
#    first_entry=0
#    last_entry=len(variants)-1
#    annotated_var_position=annotated_variant[3]
#    while first_entry<=last_entry:
#        index=(first_entry+last_entry)/2
#        variant_position=int(variants[index].split("\t")[1])
#        variant_ref=variants[index].split("\t")[3]
#        variant_alt=variants[index].split("\t")[4]
#        if annotated_var_position==variant_position and annotated_variant[5]==variant_ref and annotated_variant[6]==variant_alt:
#            return index
#        elif annotated_var_position>variant_position:
#            first_entry=index+1
#        else:
#            last_entry=index-1

if __name__ == "__main__":
    main()
