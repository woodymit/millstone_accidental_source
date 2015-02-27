"""
Utility functions for working with bam files.
"""

from functools import partial
import os
import shutil
import subprocess

from django.conf import settings

from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.models import AlignmentGroup
from main.models import ExperimentSample
from main.models import ReferenceGenome
from utils.import_util import copy_and_add_dataset_source
from utils.jbrowse_util import add_bam_file_track
from genome_finish.de_novo_bam import index_bam


def filter_bam_file_by_row(input_bam_path, filter_fn, output_bam_path):
    """Filters rows out of a bam file that don't pass a given filter function.

    This function keeps all header lines.

    Args:
        input_bam_path: Absolute path to input bam file.
        filter_fn: Function applied to each row of the input bam and returns a
            Boolean. If True, keeps the row.
        output_bam_path: Absolute path to the output bam file.
    """
    output_root = os.path.splitext(output_bam_path)[0]
    initial_sam_intermediate = output_root + '.sam'
    filtered_sam_intermediate = output_root + '.filtered.sam'
    final_bam = output_root + '.filtered.bam'

    # Convert to SAM (preserve header with -h option).
    with open(initial_sam_intermediate, 'w') as output_fh:
        p_samtools_view = subprocess.call(
                [settings.SAMTOOLS_BINARY, 'view', '-h', input_bam_path],
                stdout=output_fh)

    # Filter.
    with open(filtered_sam_intermediate, 'w') as output_fh:
        with open(initial_sam_intermediate) as input_fh:
            for line in input_fh:
                # Always write header lines.
                if line[0] == '@':
                    output_fh.write(line)
                    continue

                if filter_fn(line):
                    output_fh.write(line)
                    continue

    # Write final bam.
    with open(final_bam, 'w') as fh:
        p_samtools_view = subprocess.call(
                [settings.SAMTOOLS_BINARY, 'view', '-bS',
                        filtered_sam_intermediate],
                stdout=fh)

    # Move temp file to the original file location.
    shutil.move(final_bam, output_bam_path)

    # Delete intermediate files.
    os.remove(initial_sam_intermediate)
    os.remove(filtered_sam_intermediate)

def get_qnames_dict(file_path):
    '''
    Input: sam or bam filepath
    Retrieves a list of dicts of the format:
    {'qname':'Frag2412', '1':False, '2':True}
    '''

    file_ext = os.path.splitext(file_path)[1]
    if file_ext not in ['.bam', '.sam']:
        raise Exception('File must have extension .bam or .sam')

    if file_ext == '.bam':
        sam_path = os.path.splitext(file_path)[0] + '.temp.sam'

        # Convert to SAM (preserve header with -h option).
        with open(sam_path, 'w') as output_fh:
            p_samtools_view = subprocess.call(
                    [settings.SAMTOOLS_BINARY, 'view', '-h', file_path],
                    stdout=output_fh)
    else:
        sam_path = file_path

    qname_dict = {}
    data = open(sam_path, 'r')
    for line in data:
        if line[0] == '@':
            continue
        samList = line.strip().split('\t')
        sam = SAM(samList)
        read_data = {'qname':sam.query, 1:bool(64&sam.flag), 2:bool(128&sam.flag)}
        if not sam.query in qname_dict:
            qname_dict[sam.query] = [(bool(64&sam.flag), bool(128&sam.flag))]
        else:
            qname_dict[sam.query].append((bool(64&sam.flag), bool(128&sam.flag)))

    if file_ext == '.bam':
        os.remove(sam_path)
    return qname_dict

def qname_filter_function(qname_dict, sam_line):
    samList = sam_line.strip().split('\t')
    sam = SAM(samList)
    if sam.query in qname_dict:
        order_tups = qname_dict[sam.query]
        if (bool(64&sam.flag), bool(128&sam.flag)) in order_tups:
            return True
    return False

def get_bam_intersection(bam1, bam2, output_bam):

    output_root_1 = os.path.splitext(bam1)[0]
    intermediate_sam_1 = output_root + ".intermediate.sam"

    output_root_2 = os.path.splitext(bam2)[0]
    intermediate_sam_2 = output_root + ".intermediate.sam"

    with open(intermediate_sam_1, 'w') as output_fh:
        p_samtools_view = subprocess.call(
                [settings.SAMTOOLS_BINARY, 'view', '-h', bam1],
                stdout=output_fh)
        qname_dict=get_qname_dict(output_fh)

    filter_bam_file_by_row(bam2, partial(qname_filter_function, qname_dict=qname_dict), output_bam)

def qnames_difference(qnames1, qnames2):
    outqnames = {}
    for qname in qnames1:
        if not qname in qnames2:
            outqnames[qname] = qnames1[qname]
        else:
            remainder = set(qnames1[qname])-set(qnames2[qname])
            if remainder:
                outqnames[qname] = list(remainder)
    return outqnames


def filter_bam_by_qname_dict(input_file_path, qname_dict, output_bam_path):
    filter_bam_file_by_row(input_file_path, partial(qname_filter_function, qname_dict=qname_dict), output_bam_path)

def import_bam(reference_genome, alignment_group_name = "contig_alignment", sample_alignment_name = "contig_alignment_sample", dataset_label = "contig_aligned_bam"):
    """
    Args:
        reference_genome: The ReferenceGenome that the contigs belong to
    """

    # Make experiment sample
    experiment_sample = ExperimentSample.objects.get_or_create(
        label=sample_alignment_name,
        project=reference_genome.project)[0]

    # Make alignment grouo
    alignment_group = AlignmentGroup.objects.get_or_create(
        label=alignment_group_name,
        reference_genome=reference_genome)[0]

    # Make experiment sample to alignment
    experiment_sample_to_alignment = ExperimentSampleToAlignment.objects.get_or_create(
        alignment_group=alignment_group,
        experiment_sample=experiment_sample)[0]

    print "\n\tDEBUG: An experiment_sample_to_alignment was made\n"

    if not experiment_sample_to_alignment.dataset_set.filter(type=Dataset.TYPE.BWA_ALIGN).exists():
        print "\n\n\n DEBUG: bam dataset being copied\n\n\n"
        copy_and_add_dataset_source(experiment_sample_to_alignment,
            dataset_label,
            Dataset.TYPE.BWA_ALIGN,
            output_bam)

        # index bam file for jbrowse
        bam_path = experiment_sample_to_alignment.dataset_set.get(type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
        index_bam(bam_path)

    else:
        raise Exception("There already is an associated BWA_ALIGN dataset with the actual genome")

    add_bam_file_track(reference_genome, experiment_sample_to_alignment, Dataset.TYPE.BWA_ALIGN)
    print"\n\nDEBUG: bam file track added\n\n"

class SAM (object):
    """
    __very__ basic class for SAM input.
    """
    def __init__(self, samList = []):
        if len(samList) > 0:
            self.query    = samList[0]
            self.flag     = int(samList[1])
            self.ref      = samList[2]
            self.pos      = int(samList[3])
            self.mapq     = int(samList[4])
            self.cigar    = samList[5]
            self.matRef   = samList[6]
            self.matePos  = int(samList[7])
            self.iSize    = int(samList[8])
            self.seq      = samList[9]
            self.qual     = samList[10]
            self.tags     = samList[11:]#tags is a list of each tag:vtype:value sets
            self.valid    = 1
        else:
            self.valid = 0
            self.query = 'null'

def test_script():

    all_reads_bam = 'temp_data/projects/a04e843e/alignment_groups/c4bd8c0e/sample_alignments/8cd2c703/bwa_align.sorted.grouped.withmd.bam'
    picked_reads_bam = 'temp_data/projects/a04e843e/alignment_groups/5029f514/sample_alignments/759f3c96/bwa_align.sorted.grouped.withmd.bam'

    qnames_all = get_qname_dict(all_reads_bam)
    qnames_picked = get_qnames_dict(picked_reads_bam)

    qnames_diff = qnames_difference(qnames_all, qnames_picked)

    output_root = os.path.splitext(all_reads_bam)[0]
    ouput_bam_path = output_root + '.diff.temp.bam'

    filter_bam_by_qnames_dict(all_reads_bam, qnames_diff, output_bam_path)



   # import_bam(