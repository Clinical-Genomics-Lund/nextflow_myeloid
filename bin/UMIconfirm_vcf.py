#!/usr/bin/env python3
from pprint import pprint
import re
import sys
import subprocess
import gzip
from collections import Counter

from codecs import (open, getreader)

def parse_pileup(pile, read_pos, read_ids, read_counts, read_mismatches, duplex):

    bases = []
    i = 0

    while i < len(pile):
        entry = {}
        c = pile[i]

        next_c = ""
        if i < len(pile)-1:
            next_c = pile[i+1]

        if c in "ATGCN.":
            entry = {'base':c, 'dir':'fwd'}
            i += 1
        elif c in "atgcn,":
            entry = {'base':c.upper(), 'dir':'rev'}
            i += 1
        elif c == "$":
            i += 1
        elif c == "*":
            entry = {'base':c, 'dir':''}
            i += 1
        elif c == "^":
            entry = {'base':pile[i+2].upper(), 'dir':'rev', 'start':True}
            i += 3

        if next_c == "$":
            entry['end'] = True
        elif c != "^" and next_c != "" and next_c in "+-":
            # determine length of indel
            digit = re.findall('[\+-](\d+)[ACGTNacgtn*]+',pile[i:])[0]

            # get indel sequence
            start = i + len(digit) + 1
            indel_seq = pile[start:start+int(digit)]
            
            i += 1 + len(digit) + int(digit)
            
            dir = "rev"
            if indel_seq == indel_seq.upper():
               dir = "fwd"
               
            if next_c == "+":
                entry['ins'] = indel_seq.upper()
            elif next_c == "-":
                entry['del'] = indel_seq.upper()

        if entry:
            bases.append(entry)


    read_counts_arr = read_counts.split(',')
    read_ids_arr = read_ids.split(',')
    read_pos_arr = read_pos.split(',')
    read_mismatches_arr = read_mismatches.split(',')

    for i in range(0,len(read_pos_arr)):
        if duplex:
            bases[i]["rcnt"] = int(read_counts_arr[i*2])+int(read_counts_arr[i*2+1])
        else:
            bases[i]["rcnt"] = read_counts_arr[i]

        #bases[i]["read_id"] = read_ids_arr[i]
        bases[i]["read_pos"] = read_pos_arr[i]
        bases[i]["NM"] = read_mismatches_arr[i]
        
    return bases



def run(cmd, outfile = None):

    out = None
    if outfile:
        out = open(outfile, "w+")

    rc = subprocess.check_call(cmd, stdout=out)


    
def get_UMI_info_position(pileline, loc, duplex):

    # Get relevant data from pileup
    fields = pileline.rstrip().split("\t")
    pile =  fields[4]
    read_pos = fields[6]
    read_ids = fields[7]
    read_cnt = fields[8]
    read_mismatches = fields[9]

    # TODO: Check that loc and fields[0..1] agree
#    print("----- "+loc+ " ----- "+ fields[0]+" "+fields[1]+" "+fields[2])
    # Parse the data
    bases = parse_pileup(pile, read_pos, read_ids, read_cnt, read_mismatches, duplex)

    return bases

def main():
    
    bam = sys.argv[1]
    vcf = sys.argv[2]
    ref = sys.argv[3]
    sample_id = sys.argv[4] # FIXME: Get this from bam

    duplex = True
    
    if vcf.endswith(".gz"):
        vcf_handle = getreader('utf-8')(gzip.open(vcf), errors='replace')
    else:
        vcf_handle = open(vcf, mode='r', encoding='utf-8', errors='replace')

    # Create position file for all variants
    pos_handle = open("positions.tmp", "w")
    for vcf_line in vcf_handle:
        if vcf_line.startswith("#"):
            continue

        fields = vcf_line.split("\t")
        pos_handle.write(fields[0] + "\t" + fields[1] + "\n")
    pos_handle.close()

    vcf_handle.seek(0)
    # Run samtools mpileup
    run( ["samtools", "mpileup", bam,
          "-f", ref,
          "-l", "positions.tmp",
          "--output-QNAME",
          "--output-extra", "XZ,NM",
          "-O",
          '-o', 'piles.tmp',
          '-A',
          '-Q', '0'] )

    pile_handle = open("piles.tmp", mode='r', encoding='utf-8', errors='replace')

        
    header = []
    sample_idx = 0

    prev_chr = ""
    prev_pos = 0
    prev_umi_data = []
    for vcf_line in vcf_handle:
        vcf_line = vcf_line.rstrip()
        
        if vcf_line.startswith("##"):
            print(vcf_line)
            continue

        elif vcf_line.startswith("#"):
            print(vcf_line)
            header = vcf_line.split("\t")
            sample_idx = header.index(sample_id)

        else:
            variant = vcf_line.split("\t")

            ref_bases = variant[3]
            alt_bases = variant[4]
            #print(vcf_line)
            pos = variant[0]+":"+variant[1]+"-"+variant[1]

            if variant[0] != prev_chr or variant[1] != prev_pos:
                pileline = pile_handle.readline()
                umi_data = get_UMI_info_position(pileline, pos, duplex)

            family_size_counts = Counter()
            
            # SNVs
            if len(ref_bases) == 1 and len(alt_bases) == 1:
                for read in umi_data:
                    if read["base"] == alt_bases:
                        family_size_counts.update(str(read['rcnt']))
                        #pprint(read)

            # Insertions
            elif len(ref_bases) == 1 and len(alt_bases) > 1:
                ins_seq = alt_bases[1:]
                for read in umi_data:
                    if read.get("ins") == ins_seq:
                        family_size_counts.update(str(read['rcnt']))
                        #pprint(read)

            # Deletions
            elif len(alt_bases) == 1 and len(ref_bases) > 1:
                del_seq = ref_bases[1:]
                for read in umi_data:
                    if read.get("del") == del_seq:
                        family_size_counts.update(str(read['rcnt']))
                        #pprint(read)

            # Summarize UMI family sizes
            counts = []
            for el in sorted(family_size_counts.items()):
                counts.append(el[0]+"|"+str(el[1]))

            # Add to sample fields of vcf
            umi_counts_str = ",".join(counts)
            if not variant[8].endswith(":UMI"):
                variant[8] += ":UMI"

            if not umi_counts_str:
                umi_counts_str = "0"
            variant[sample_idx] += ":"+umi_counts_str

            # Print VCF line with new info
            print("\t".join(variant))

            
            prev_chr = variant[0]
            prev_pos = variant[1]
            prev_umi_data = umi_data
            
if __name__ == "__main__":
    main()



