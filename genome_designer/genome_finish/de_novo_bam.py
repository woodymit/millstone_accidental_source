import subprocess
import sys
import os

from django.conf import settings
from utils import convert_fasta_to_fastq
import millstone_de_novo_fns

BWA_BINARY = settings.TOOLS_DIR + '/bwa/bwa'
VELVETH_BINARY = settings.TOOLS_DIR + '/velvet/velveth'
VELVETG_BINARY = settings.TOOLS_DIR + '/velvet/velvetg'

def truncate_ext(path, count=1):
    count = count -1
    if count == 0:
        return path[:path.rindex(".")]
    else:
        return truncate_ext(path[:path.rindex(".")], count)


def bwa_align(reads, ref_fasta, data_dir):
    
    # 0. Interpret reads file type
    if all([r.endswith(".fa") for r in reads]):
        reads_fq = [os.path.join(data_dir, truncate_ext(r) + ".fq") for r in reads]
        for i,r in enumerate(reads):
            convert_fasta_to_fastq(r, reads_fq[i])
    elif all([r.endswith(".fq") for r in reads]):
        reads_fq = reads
    else:
        raise(Exception("All reads must have file extension .fq or .fa"))

    # 1. bwa index ref.fa #TODO: Check if already indexed
    cmd = "{bwa} index {ref_path}".format(
        bwa=BWA_BINARY,
        ref_path=ref_fasta)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # 2. bwa mem ref.fa contigs.fq > alignment.sam
    alignment_unsorted_sam = os.path.join(data_dir, "bwa_align.alignment.unsorted.sam")

    cmd = "{bwa} mem {ref_fasta} {contigs_fastq} > {alignment_sam}".format(
        bwa=BWA_BINARY,
        contigs_fastq=" ".join(reads_fq),
        ref_fasta=ref_fasta,
        alignment_sam=alignment_unsorted_sam)
    
    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # 3. samtools view -b alignment.sam > alignment.bam
    alignment_unsorted_bam = truncate_ext(alignment_unsorted_sam) + ".bam"
    
    cmd = "{samtools} view -b -S {alignment_sam} > {alignment_bam}".format(
        samtools=settings.SAMTOOLS_BINARY,
        alignment_sam=alignment_unsorted_sam,
        alignment_bam=alignment_unsorted_bam)

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # 4. samtools sort alignment.bam alignment.sorted
    alignment_sorted_bam = truncate_ext(alignment_unsorted_bam, 2) + ".sorted.bam"
    
    cmd = "{samtools} sort {alignment_bam} {alignment_sorted}".format(
        samtools=settings.SAMTOOLS_BINARY,
        alignment_bam=alignment_unsorted_bam,
        alignment_sorted=truncate_ext(alignment_sorted_bam))

    subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)


    # # 5. samtools index alignment.sorted.bam
    # cmd = "{samtools} index {alignment_sorted_bam}".format(
    #     samtools=settings.SAMTOOLS_BINARY,
    #     alignment_sorted_bam=alignment_sorted_bam)

    # subprocess.call(cmd, shell=True, executable=settings.BASH_PATH)

    # 6. Index it
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


def rmdup(input_bam_file, output_bam_file):
    # Store input sam header
    output_sam = truncate_ext(output_bam_file) + ".sam"
    subprocess.check_call(
        ' '.join([settings.SAMTOOLS_BINARY, 'view', '-H', '-o', output_sam,  input_bam_file]),
        shell=True, executable=settings.BASH_PATH)

    # Convert input to sam, sort, remove duplicate adjacent lines, and append to header
    cmd = ' | '.join([
        settings.SAMTOOLS_BINARY + ' view ' + input_bam_file,
        'sort',
        'uniq'
        ]) + ' >> ' + output_sam
    subprocess.check_call(cmd, shell=True, executable=settings.BASH_PATH)

    make_bam(output_sam, output_bam_file)


def run_velvet(reads, output_dir, opt_dict = None):

    default_opts = {
        'velveth': {
            'hash_length': 21
        },
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
