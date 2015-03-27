"""Script to combine ReferenceGenomes into a single one.
"""

import os
import re

from Bio import SeqIO
from django.core.exceptions import ObjectDoesNotExist

from main.models import Dataset
from main.models import ReferenceGenome
from utils.import_util import add_dataset_to_entity
from utils import generate_safe_filename_prefix_from_label


def combine_list(reference_genome_list, new_ref_genome_label, project):
    """Combine ReferenceGenomes into a new single ReferenceGenome
    composed of the component parts.

    Args:
        reference_genome_list: List of ReferenceGenome objects.
        new_ref_genome_label: Label for the new ReferenceGenome.
        project: Project to which the new ReferenceGenome will be added.

    Returns:
        Object with keys:
            * is_success
            * new_reference_genome (when is_success = True)
            * error_msg (when is_success = False)
    """
    ref_genome_dataset_type = Dataset.TYPE.REFERENCE_GENOME_GENBANK
    # Grab the Genbank dataset.
    genbank_dataset_list = []
    for ref_genome in reference_genome_list:
        try:
            genbank_dataset = ref_genome.dataset_set.get(
                    type=ref_genome_dataset_type)
            assert os.path.exists(genbank_dataset.get_absolute_location())
            genbank_dataset_list.append(genbank_dataset)
        except ObjectDoesNotExist:
            return {
                'is_success': False,
                'error_msg': 'Only support combining Genbank genomes.'
            }

    # Read the datasets into Biopython SeqRecord objects.
    seq_record_list = []
    for dataset in genbank_dataset_list:
        with open(dataset.get_absolute_location()) as input_fh:
            for record in SeqIO.parse(input_fh, 'genbank'):
                seq_record_list.append(record)

    # Create a new ReferenceGenome.
    new_ref_genome = ReferenceGenome.objects.create(
            project=project,
            label=new_ref_genome_label,
            num_chromosomes=len(seq_record_list),
            num_bases=sum([len(seq) for seq in seq_record_list]))

    # Generate a filename from the label with non-alphanumeric characters replaced
    # by underscores.
    filename_prefix = generate_safe_filename_prefix_from_label(
            new_ref_genome_label)
    filename = filename_prefix + '.gb'
    new_genbank_dest = os.path.join(new_ref_genome.get_model_data_dir(), filename)

    # Write the result.
    with open(new_genbank_dest, 'w') as output_fh:
        SeqIO.write(seq_record_list, output_fh, 'genbank')

    # Create a dataset which will point to the file. This step must happen after
    # writing the file because a signal will be triggered which requires the
    # Genbank to exist already.
    new_ref_genome_dataset = add_dataset_to_entity(new_ref_genome,
            ref_genome_dataset_type, ref_genome_dataset_type, new_genbank_dest)

    return {
        'is_success': True,
        'new_reference_genome': new_ref_genome
    }

def combine_list_allformats(reference_genome_list, new_ref_genome_label, project):
    """Combine ReferenceGenomes into a new single ReferenceGenome
    composed of the component parts.

    Args:
        reference_genome_list: List of ReferenceGenome objects.
        new_ref_genome_label: Label for the new ReferenceGenome.
        project: Project to which the new ReferenceGenome will be added.

    Returns:
        Object with keys:
            * is_success
            * new_reference_genome (when is_success = True)
            * error_msg (when is_success = False)
    """

    dataset_to_seqio_format = {
        Dataset.TYPE.REFERENCE_GENOME_GENBANK:'genbank',
        Dataset.TYPE.REFERENCE_GENOME_FASTA: 'fasta'
        }

    dataset_list = []
    for ref_genome in reference_genome_list:
        dataset = None
        for dataset_type in dataset_to_seqio_format.keys():
            filter_result = ref_genome.dataset_set.filter(type=dataset_type)
            if len(filter_result):
                dataset=filter_result[0]
                break
        if not dataset or not os.path.exists(dataset.get_absolute_location()):
            return {
                'is_success': False,
                'error_msg': 'All referenge genomes must have an associated Fasta or Genbank dataset'
            }
        else:
            dataset_list.append(dataset)

    # Read the datasets into Biopython SeqRecord objects.
    seq_record_list = []
    for dataset in dataset_list:
        with open(dataset.get_absolute_location()) as input_fh:
            for record in SeqIO.parse(input_fh, dataset_to_seqio_format[dataset.type]):
                seq_record_list.append(record)

    # Create a new ReferenceGenome.
    new_ref_genome = ReferenceGenome.objects.create(
            project=project,
            label=new_ref_genome_label,
            num_chromosomes=len(seq_record_list),
            num_bases=sum([len(seq) for seq in seq_record_list]))

    # Generate a filename from the label with non-alphanumeric characters replaced
    # by underscores.
    filename_prefix = generate_safe_filename_prefix_from_label(
            new_ref_genome_label)
    list_includes_genbank = Dataset.TYPE.REFERENCE_GENOME_GENBANK in [dataset.type for dataset in dataset_list]
    if list_includes_genbank:
        filename = filename_prefix + '.gb'
    else:
        filename = filename_prefix + '.fa'
    new_file_dest = os.path.join(new_ref_genome.get_model_data_dir(), filename)

    # Write the result.
    ref_genome_dataset_type = Dataset.TYPE.REFERENCE_GENOME_GENBANK if list_includes_genbank else Dataset.TYPE.REFERENCE_GENOME_FASTA
    output_file_format = dataset_to_seqio_format[ref_genome_dataset_type]
    with open(new_file_dest, 'w') as output_fh:
        SeqIO.write(seq_record_list, output_fh, output_file_format)

    # Create a dataset which will point to the file. This step must happen after
    # writing the file because a signal will be triggered which requires the
    # Genbank to exist already.

    print 'new_ref_genome:',new_ref_genome
    print 'ref_genome_dataset_type:',ref_genome_dataset_type
    print 'new_file_dest:', new_file_dest

    new_ref_genome_dataset = add_dataset_to_entity(new_ref_genome,
            ref_genome_dataset_type, ref_genome_dataset_type, new_file_dest)

    return {
        'is_success': True,
        'new_reference_genome': new_ref_genome
    }
