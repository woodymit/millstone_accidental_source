import subprocess
import sys
import os

from django.conf import settings
from utils import convert_fasta_to_fastq
import millstone_de_novo_fns

BWA_BINARY = settings.TOOLS_DIR + '/bwa/bwa'
VELVETH_BINARY = settings.TOOLS_DIR + '/velvet/velveth'
VELVETG_BINARY = settings.TOOLS_DIR + '/velvet/velvetg'

"""
Workflow:

ref.fa reads.fq

    RUN make_bam

raw_alignment.bam

    RUN get_unmapped_reads
    RUN get_split_reads

split_unmapped.bam

    RUN add_mate_pairs

reads_to_assemble.bam

    RUN run_velvet (on bam)

contigs.fa

    RUN make_bam

contigs.alignment.sorted.bam

    Evaluate assembly with IGV

"""

def truncate_ext(path, count=1):
    count = count -1
    if count == 0:
        return path[:path.rindex(".")]
    else:
        return truncate_ext(path[:path.rindex(".")], count)


def bwa_align(reads, ref_fasta, data_dir):
    
#   1. python fasta_to_fastq.py reads.fa
    # reads_fq = []
    # if len(reads) == 1:
    #     if reads[0].endswith(".fq"):
    #         reads_fq = reads
    #     elif reads[0].endswith(".fa"):
    #         fastq = truncate_ext(reads[0]) + ".fq"
    #         reads_fq.append(fastq)
    #         if not os.path.exists(fastq):
    #             convert_fasta_to_fastq(reads[0], fastq)
    #     else:
    #         raise Exception("Reads must be in fastq or fasta format")

    # elif reads[0].endswith(".fq") and reads[1].endswith(".fq"):
    #     reads_fq = reads

    # elif reads[0].endswith(".fa") and reads[0].endswith(".fa"):
    #     for r in reads:
    #         fastq = truncate_ext(r) + ".fq"
    #         reads_fq.append(fastq)
    #         if not os.path.exists(fastq):
    #             convert_fasta_to_fastq(r, fastq)
    # else:
    #     raise Exception("Both reads must be in fastq or fasta format")

    # Cleaner

    if all([r.endswith(".fa") for r in reads]):
        reads_fq = [os.path.join(data_dir, truncate_ext(r) + ".fq") for r in reads]
        for i,r in enumerate(reads):
            convert_fasta_to_fastq(r, reads_fq[i])
    elif all([r.endswith(".fq") for r in reads]):
        reads_fq = reads
    else:
        raise(Exception("All reads must have file extension .fq or .fa"))

#   2. bwa index ref.fa #TODO: Check if already indexed
    cmd = "{bwa} index {ref_path}".format(
        bwa=BWA_BINARY,
        ref_path=ref_fasta)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

#   3. bwa mem ref.fa contigs.fq > alignment.sam
    alignment_unsorted_sam = os.path.join(data_dir, "bwa_align.alignment.unsorted.sam")

    cmd = "{bwa} mem {ref_fasta} {contigs_fastq} > {alignment_sam}".format(
        bwa=BWA_BINARY,
        contigs_fastq=" ".join(reads_fq),
        ref_fasta=ref_fasta,
        alignment_sam=alignment_unsorted_sam)
    
    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

#   4. samtools view -b alignment.sam > alignment.bam
    alignment_unsorted_bam = truncate_ext(alignment_unsorted_sam) + ".bam"
    
    cmd = "{samtools} view -b -S {alignment_sam} > {alignment_bam}".format(
        samtools=settings.SAMTOOLS_BINARY,
        alignment_sam=alignment_unsorted_sam,
        alignment_bam=alignment_unsorted_bam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

#   5. samtools sort alignment.bam alignment.sorted
    alignment_sorted_bam = truncate_ext(alignment_unsorted_bam, 2) + ".sorted.bam"
    
    cmd = "{samtools} sort {alignment_bam} {alignment_sorted}".format(
        samtools=settings.SAMTOOLS_BINARY,
        alignment_bam=alignment_unsorted_bam,
        alignment_sorted=truncate_ext(alignment_sorted_bam))

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


#   6. samtools index alignment.sorted.bam
    # TODO: Check if already existing index
    cmd = "{samtools} index {alignment_sorted_bam}".format(
        samtools=settings.SAMTOOLS_BINARY,
        alignment_sorted_bam=alignment_sorted_bam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

#   7. Index it
    index_bam(alignment_sorted_bam)

    return alignment_sorted_bam


def index_bam(bam):
    cmd = "{samtools} index {bam}".format(
        samtools=settings.SAMTOOLS_BINARY,
        bam=bam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

def sort_bam(input_bam, output_bam=None):
    if output_bam == None:
        output_bam = input_bam

    cmd = "{samtools} sort {alignment_bam} {alignment_sorted}".format(
        samtools=settings.SAMTOOLS_BINARY,
        alignment_bam=input_bam,
        alignment_sorted=output_bam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def make_sam(bam, sam_filename = None):
    if sam_filename == None:
        sam_filename = truncate_ext(bam) + ".sam"
    
    cmd = "{samtools} view -h {bam} > {sam}".format(
        samtools=settings.SAMTOOLS_BINARY,
        bam=bam,
        sam=sam_filename)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def make_bam(sam, bam_filename = None):
    if bam_filename == None:
        bam_filename = truncate_ext(sam) + ".bam"
    
    cmd = "{samtools} view -b -S {sam} > {bam}".format(
        samtools=settings.SAMTOOLS_BINARY,
        sam=sam,
        bam=bam_filename)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def concatenate_bams(bam_list, output):
    
    cmd = "{samtools} cat -o {output} {bam_files}".format(
        samtools=settings.SAMTOOLS_BINARY,
        bam_files=" ".join(bam_list),
        output=output)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


def rmdup(bam, output = None):
    if output == None:
        output = bam
    
    size_init = os.stat(bam).st_size

    cmd = "{samtools} rmdup {bam} {output}".format(
        samtools=settings.SAMTOOLS_BINARY,
        bam=bam,
        output=output)

    subprocess.call(cmd, shell = True, executable=settings.BASH_PATH)

    size_fin = os.stat(output).st_size

    print "Remove duplicates\n\tInitial bam size:\t " + str(size_init) + "\n\tFinal bam size:\t " + str(size_fin)


def run_velvet(reads, output_dir, opt_dict = None):

    default_opts = {
        'velveth': {
            'hash_length': 21
        },
        # 'velvetg': {
        #     'cov_cutoff': 20,
        #     'ins_length': 200,
        #     'ins_length_sd': 90
        # }
    }

    if not opt_dict:
        opt_dict = default_opts

    if not 'hash_length' in opt_dict['velveth']:
        raise Exception("If passing an option dictionary, an output_dir_name key \
        must exist in the velveth sub-dictionary")

    l = [
        VELVETH_BINARY,
        output_dir,
        str(opt_dict['velveth']['hash_length'])] + ["-" + k + " " + str(opt_dict['velveth'][k])
            for k in opt_dict['velveth']
            if k not in ['hash_length']] + ["-bam"] +   [reads]
    cmd = ' '.join(l)
    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    cmd = ' '.join([
    VELVETG_BINARY,
    output_dir] +
    ["-" + k + " " + str(opt_dict['velvetg'][k])
        for k in opt_dict['velvetg']])
    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


    # if not no_defaults:
    #     # First run velveth to build hash of reads.
    #     # TODO: What if we used shortPaired read type?
    #     cmd = '{velveth} {output_dir} {hash_length} -bam {reads}'.format(
    #             velveth=VELVETH_BINARY,
    #             output_dir=output_dir_name,
    #             hash_length=hash_length,
    #             reads=reads)



    #     subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    #     # Then run velvetg to build graph and find contigs.
    #     cmd = '{velvetg} {output_dir} -cov_cutoff {cov_cutoff} -ins_length {ins_length} -ins_length_sd {ins_length_sd}'.format(
    #             velvetg=VELVETG_BINARY,
    #             output_dir=output_dir_name,
    #             cov_cutoff=cov_cutoff,
    #             ins_length=ins_length,
    #             ins_length_sd=ins_length_sd)

    #     subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # else:
    #     # First run velveth to build hash of reads.
    #     # TODO: What if we used shortPaired read type?
    #     cmd = '{velveth} {output_dir} {hash_length} -bam -shortPaired {reads}'.format(
    #             velveth=VELVETH_BINARY,
    #             output_dir=output_dir_name,
    #             hash_length=hash_length,
    #             reads=reads)

    #     subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    #     # Then run velvetg to build graph and find contigs.
    #     cmd = '{velvetg} {output_dir} -max_branch_length 19'.format(
    #             velvetg=VELVETG_BINARY,
    #             output_dir=output_dir_name)

    #     subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)
