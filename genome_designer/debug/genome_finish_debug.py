# Debug script for genome finishing
import os

from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.models import AlignmentGroup
from main.models import ExperimentSample
from main.models import ReferenceGenome
from utils.import_util import copy_and_add_dataset_source
from utils.jbrowse_util import add_bam_file_track
from genome_finish.de_novo_bam import bwa_align
from genome_finish.de_novo_bam import index_bam

def run_debug(sample_alignment_name = "contig_alignment_sample"):
    ref = ReferenceGenome.objects.get(label="gf_test_ref")
    ins = ReferenceGenome.objects.get(label="1kb_insertion")
    align_contigs_to_actual_genome(ref, ins, sample_alignment_name)

def run_debug_2(sample_alignment_name = "contig_alignment_sample"):
    ref = ReferenceGenome.objects.get(label="insertion_1_ref")
    ins = ReferenceGenome.objects.get(label="insertion_1_transformed")
    align_contigs_to_actual_genome(ref, ins, sample_alignment_name)

def align_contigs_to_actual_genome(reference_genome, actual_genome, sample_alignment_name = "contig_alignment_sample"):
    """
    Args:
        reference_genome: The ReferenceGenome that the contigs belong to
        actual_genome: The ReferenceGenome with structural variants
    """

    # Get reference_genome fasta
    ref_fasta = actual_genome.dataset_set.get(type = Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # Make experiment sample
    experiment_sample = ExperimentSample.objects.get_or_create(
        label=sample_alignment_name,
        project=reference_genome.project)[0]

    # Make alignment grouo
    alignment_group = AlignmentGroup.objects.get_or_create(
        label="contig_alignment",
        reference_genome=actual_genome)[0]

    # Make experiment sample to alignment
    experiment_sample_to_alignment = ExperimentSampleToAlignment.objects.get_or_create(
        alignment_group=alignment_group,
        experiment_sample=experiment_sample)[0]

    print "\n\tDEBUG: An experiment_sample_to_alignment was made\n"

    # Make data_dir directory to house genome_finishing files
    project_dir = reference_genome.project.get_model_data_dir()
    data_dir = os.path.join(project_dir, "genome_finishing", experiment_sample_to_alignment.uid)
 
    # Make data_dir directory if it does not exist
    if not os.path.exists(os.path.join(project_dir, "genome_finishing")):
        os.mkdir(os.path.join(project_dir, "genome_finishing"))
        os.mkdir(data_dir)
    elif not os.path.exists(data_dir):
        os.mkdir(data_dir)

    print "data_dir: " + str(data_dir)

    # Get contigs
    contigs = reference_genome.dataset_set.get(type = Dataset.TYPE.CONTIGS_FASTA).get_absolute_location()
    print "contigs: " + str(contigs)
    # Run alignment of contigs to actual_genome
    output_bam = bwa_align([contigs], ref_fasta, data_dir)
    print "output_bam: " + str(output_bam)

    if not experiment_sample_to_alignment.dataset_set.filter(type=Dataset.TYPE.BWA_ALIGN).exists():
        print "\n\n\n DEBUG: bam dataset being copied\n\n\n"
        copy_and_add_dataset_source(experiment_sample_to_alignment,
            "contig_aligned_bam",
            Dataset.TYPE.BWA_ALIGN,
            output_bam)

        # index bam file for jbrowse
        bam_path = experiment_sample_to_alignment.dataset_set.get(type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
        index_bam(bam_path)

    else:
        raise Exception("There already is an associated BWA_ALIGN dataset with the actual genome")

    add_bam_file_track(actual_genome, experiment_sample_to_alignment, Dataset.TYPE.BWA_ALIGN)
    print"\n\nDEBUG: bam file track added\n\n"
