'''
Download complete record from GeneBank in .gb format. Also download fasta. Run
this script to get gtf, and with fasta, run cellranger mkref.

Neal G. Ravindra, 200402
'''


import os
from Bio import SeqIO

def gb2gtf(fname,included_types=['CDS'],outfile=None):
    """Download complete GeneBank record to get gb file and feed it here.

    NOTE: not sure about \n carriage in gtf line


    Args:
      fname (str): full filepath
      included_types (list): list of type names to include in annotation file (.gtf)

    Returns:
      saved gtf to same directory as fname

    References:
      https://useast.ensembl.org/info/website/upload/gff.html
    """
    print('Converting {} to .gtf'.format(os.path.split(fname)[1]))
    if outfile is None:
        outfile = fname.split('.gb')[0]+'.gtf'
    with open(outfile,'a') as out:
        for gb in SeqIO.parse(fname,'gb'):
            for f in gb.features:
                seqname = gb.id # seqqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly.

                # source isn't well defined, take acc and skip source
                if f.type=='source':
                    source = 'GeneBank'
                    continue
                elif f.type not in included_types:
                    continue # skip gene type because it seems to get duplicated
                else:
                    if (f.type=='CDS') or (f.type=='mat_peptide'):
                        feature = 'exon' # treat CDS as exons, though exons technically contain UTRs, CDS does not
                    else:
                        feature = f.type
                    start = f.location.start+1
                    end = f.location.end
                    score = '.' # a floating point value
                    if f.strand > 0:
                        strand = '+'
                    else:
                        strand = '-'

                    # get attributes
                    if 'codon_start' in f.qualifiers:
                        frame = f.qualifiers['codon_start'][0]
                    else:
                        frame = '.'
                    if ('gene' not in f.qualifiers) and ('product' in f.qualifiers):
                        if ('protein_id' in f.qualifiers):
                            # capture hrv protein for eff
                            gene_id = f.qualifiers['product'][0]
                            transcript_id = gene_id
                            gene_name = '{}_{}'.format(f.qualifiers['protein_id'][0],gene_id)
                            attribute = 'gene_id "%s";transcript_id "%s";gene_name "%s";' % (gene_id,transcript_id,gene_name)
                        else :
                            gene_id = f.qualifiers['product'][0]
                            transcript_id = f.qualifiers['product'][0]
                            gene_name = f.qualifiers['product'][0]
                            attribute = 'gene_id "%s";transcript_id "%s";gene_name "%s";' % (gene_id,transcript_id,gene_name)
                    elif ('protein_id' in f.qualifiers) and ('gene' in f.qualifiers) and ('product' in f.qualifiers):
                        gene_id = f.qualifiers['gene'][0]
                        transcript_id = f.qualifiers['protein_id'][0]
                        gene_name = f.qualifiers['product'][0]
                        attribute = 'gene_id "%s";transcript_id "%s";gene_name "%s";' % (gene_id,transcript_id,gene_name)
                    else :
                        attribute = ''
                        # collect all other attributes
                        for i in f.qualifiers:
                            if (i != 'translation') and (i != 'codon_start') :
                                if i=='gene':
                                    attribute+='gene_id "%s";' % f.qualifiers[i][0]
                                    attribute+='transcript_id "%s";' % f.qualifiers[i][0]
                                    attribute+='gene_name "%s";' % f.qualifiers[i][0]
                                else:
                                    attribute+='%s "%s";' % (i,f.qualifiers[i][0])
                    # add spaces to attribute
                    attribute = attribute[:-1].replace(';','; ') + ';' # don't replace last semicolon

                # add feature to file
                gtf = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (seqname,source,feature,start,end,score,strand,frame,attribute)
                out.write(gtf)
        out.close()
        print('... converted {}'.format(os.path.split(outfile)[1]))


if __name__ == '__main__' :
    # target for rhv
    ## KC894166.1      genebank        exon    592     7065    .       +       .       gene_id "HRV-1A_P1"; transcript_id "HRV-1A_P1"; gene_name "HRV-1A_P1";
    sc2 = '/home/ngr/ushare/sccovid/data/processed/MT020880.1.gb'
    rhv = '/home/ngr/ushare/sccovid/data/processed/HRV-1A_P1.gb'
    rhv_subreads = '/home/ngr/ushare/sccovid/data/processed/HRV-1A_P1_mat_peptide.gtf'
    gb2gtf(sc2)
    gb2gtf(rhv)
    gb2gtf(rhv,included_types=['mat_peptide'],outfile=rhv_subreads)
