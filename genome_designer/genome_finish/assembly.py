import sys
import os

from genome_finish import de_novo_bam
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.model_utils import get_dataset_with_type

VELVET_COVERAGE_CUTOFF = 3
VELVET_KMER_LIST = [21]

def truncate_ext(path, count=1):
    count = count -1
    if count == 0:
        return path[:path.rindex(".")]
    else:
        return truncate_ext(path[:path.rindex(".")], count)

def generate_contigs(experimentSampleToAlignment, contig_ref_genome):

    # Get fasta read files
    experimentSample = experimentSampleToAlignment.experiment_sample
    fastq1 = experimentSample.dataset_set.get(type=Dataset.TYPE.FASTQ1).get_absolute_location()
    fastq2 = experimentSample.dataset_set.get(type=Dataset.TYPE.FASTQ2).get_absolute_location()
    read_fastqs = [fastq1, fastq2]

    # Get fasta reference genome file
    referenceGenome = experimentSampleToAlignment.alignment_group.reference_genome
    ref_fasta = referenceGenome.dataset_set.get(type = Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    
    # Make data_dir directory to house genome_finishing files
    contig_dir = contig_ref_genome.get_model_data_dir()
    data_dir = os.path.join(contig_dir, "genome_finishing")

    # Make data_dir directory if it does not exist
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    # ---------- Begin Genome Finishing ---------- #

    # Retrieve bwa mem .bam alignment if exists otherwise generate it
    if experimentSampleToAlignment.dataset_set.filter(type=Dataset.TYPE.BWA_ALIGN).exists():
        alignment_sorted = get_dataset_with_type(experimentSampleToAlignment, Dataset.TYPE.BWA_ALIGN).get_absolute_location()
    else:
        alignment_sorted = de_novo_bam.bwa_align(read_fastqs, ref_fasta, data_dir)    

    alignment_prefix = os.path.join(data_dir, 'bwa_align')

    # Extract SV indicating reads
    unmapped_reads = alignment_prefix + ".unmapped.bam"
    de_novo_bam.millstone_de_novo_fns.get_unmapped_reads(alignment_sorted, unmapped_reads)

    split_reads = alignment_prefix + ".split.bam"
    de_novo_bam.millstone_de_novo_fns.get_split_reads(alignment_sorted, split_reads)

    clipped_reads = alignment_prefix + ".clipped.bam"
    de_novo_bam.millstone_de_novo_fns.get_clipped_reads(alignment_sorted, clipped_reads)

    # Aggregate SV indicants
    SV_indicants_bam = alignment_prefix + ".SV_indicants.bam"
    de_novo_bam.concatenate_bams([unmapped_reads, split_reads, clipped_reads], SV_indicants_bam)

    # Remove duplicates -- samtools rmdup has known bugs
    SV_indicants_no_dups_bam = alignment_prefix + ".SV_indicants_no_dups.bam"
    de_novo_bam.rmdup(SV_indicants_bam, SV_indicants_no_dups_bam) 

    # Convert SV indicants bam to sam
    SV_indicants_sam = alignment_prefix + ".SV_indicants.sam"
    de_novo_bam.make_sam(SV_indicants_no_dups_bam, SV_indicants_sam)

    # Add mate pairs to SV indicants sam
    SV_indicants_with_pairs_sam = alignment_prefix + ".SV_indicants_with_pairs.sam"
    de_novo_bam.millstone_de_novo_fns.add_paired_mates(SV_indicants_sam, alignment_sorted, SV_indicants_with_pairs_sam)

    # Make bam of SV indicants w/mate pairs
    SV_indicants_with_pairs_bam = alignment_prefix + ".SV_indicants_with_pairs.bam"
    de_novo_bam.make_bam(SV_indicants_with_pairs_sam, SV_indicants_with_pairs_bam)

    # Velvet assembly
    contig_files = []
    opt_dict = {
        'velveth': {
            'shortPaired': ''
        },
        'velvetg': {
            'cov_cutoff':VELVET_COVERAGE_CUTOFF
        }
    }
    kmer_list = VELVET_KMER_LIST
    for kmer_length in kmer_list:
        # Set hash length argument for velveth
        opt_dict['velveth']['hash_length'] = kmer_length

        # Run velvet assembly on SV indicants
        velvet_dir = os.path.join(data_dir, "velvet_k" + str(kmer_length))
        de_novo_bam.run_velvet(
                SV_indicants_with_pairs_bam,
                velvet_dir,
                opt_dict)

        # Collect resulting contigs fasta
        contigs_fasta = os.path.join(velvet_dir, "contigs.fa")
        contig_files.append(contigs_fasta)

    return contig_files
