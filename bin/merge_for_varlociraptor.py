##!/usr/bin/python3
import os 
import fnmatch
import gzip 
from optparse import OptionParser

def merge_vcfs(all_vcfs,analysis_dirPath,fout_name):
    """
    This  function (merge_vcfs) merges a list of vcfs for varlociraptor sw. 
    From header lines only lines with contig,INFO and ALT and from columns 
    CHROM, POS, ID, REF, ALT and INFO are written into outfile.
    """
    #print(merge_vcfs.__doc__)
    print('Merging...')
    header_lines=[]
    lines= []
    tags=[]
    format_line = None
    for  file in all_vcfs:
        file_path = analysis_dirPath + '/' + file
        #print(file_path)
        with gzip.open(file_path, 'rt') as fin, open(fout_name , 'w') as merged_file:
            for line in fin:
                    line = line.rstrip()
                    #try to extract INFO/contig tag lines from headers.
                    if  line.startswith('##'):
                        if line.startswith('##fileformat'):
                            if  format_line is None:
                                format_line=line
                                header_lines.append(format_line)
                        if line.startswith('##contig') or line.startswith('##INFO')  or  line.startswith('ALT') :
                            if line.split('<')[1].split(',')[0] not in tags:
                                tags.append(line.split('<')[1].split(',')[0])
                                header_lines.append(line)
                        
                    elif  line.startswith('#CHROM'):
                            continue
                    #try to get file contents
                    else:
                        cols=line.split('\t')
                        new_line = '\t'.join(cols[0:5])
                        #check if  'SVLEN' is present
                        if  'SVLEN' in cols[7]:
                            new_line = new_line + '\t' + '.' + '\t' + '.' + '\t'+ cols[7]
                        else:
                            new_line =  new_line + '\t' + '.' + '\t' + '.' + '\t' + '.'
                        if new_line not in lines:
                            lines.append(new_line)
            #collect all lines
            all_lines = header_lines + ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'] +  lines
            for line in all_lines:
                print(line, file= merged_file) #print lines into output file
    print("Merging vcf files has been done.\n#headers-->{}\n#lines-->{}".format(len(header_lines),len(lines)))        

    


def main():
    """    
    This script merges the vcf files from different callers which are
    given as options. Please notice that the name of each vcf file has to satisfy this pattern:
    '*nameofcaller.vcf.gz (e.g. *vardict.vcf.gz)'.
    The name of out file and the directory where all vcf are located in are required to be passed as other options.
    In the output file, All unique varinat lines and contig,INFO,ALT header lines are written.        
    """ 
    print(main.__doc__)
    Parser = OptionParser()
    Parser.add_option("--callers", action='store', dest='callers',help="A list of variant callers to be included.")
    Parser.add_option("--output", action='store', dest='output_name', help = "Name of outfile")
    Parser.add_option("--dir", action='store', dest='analysis_dirPath', help= "The PATH of the directory where vcf files locate.\nAccepted pattern for vcf files: *nameofcaller.vcf.gz (e.g. *vardict.vcf.gz)")
    

    #(options, args) = Parser.parse_args()
    options=  Parser.parse_args()[0]

    callers= options.callers
    callers = callers.split(',')
    analysis_dirPath = options.analysis_dirPath
    fout_name = options.output_name

    all_vcfs =[]

    for name in callers:
        pattern = '*' + name + '.vcf.gz'
        #globals()[name] = fnmatch.filter(os.listdir(analysis_dirPath), pattern)
        files = fnmatch.filter(os.listdir(analysis_dirPath), pattern)
        all_vcfs = all_vcfs + files
    print('Number of vcf files--> ', len(all_vcfs))
    merge_vcfs(all_vcfs=all_vcfs,analysis_dirPath=analysis_dirPath,fout_name=fout_name)

    #analysis_dirPath = '/fs1/results/myeloid/vcf/'
    #freebayes = fnmatch.filter(os.listdir(analysis_dirPath), '*freebayes.vcf.gz')
    #vardict = fnmatch.filter(os.listdir(analysis_dirPath), '*vardict.vcf.gz')
    #tnscope = fnmatch.filter(os.listdir(analysis_dirPath), '*tnscope.vcf.gz')
    #all_vcfs= freebayes +  vardict +  tnscope


if __name__ == '__main__':
    main()
