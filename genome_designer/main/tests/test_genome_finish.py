"""
Tests for genome finishing features
"""

import json
import os
import subprocess

from django.conf import settings
from django.contrib.auth import authenticate
from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.http.request import HttpRequest
from django.test import TestCase
from django.test import Client

from genome_finish.de_novo_bam import bwa_align
import genome_finish.assembly as assembly
import genome_finish.de_novo_bam as de_novo_bam
from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.models import ReferenceGenome
from main.model_utils import get_dataset_with_type
import main.xhr_handlers as xhr_handlers
from pipeline.pipeline_runner import run_pipeline
from utils import convert_fasta_to_fastq
from utils.import_util import import_reference_genome_from_local_file
from utils.import_util import add_dataset_to_entity

from genome_finish.millstone_de_novo_fns import get_match_counts

TEST_USERNAME = 'testuser'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@example.com'
TEST_PROJECT_NAME = 'testModels_project'
TEST_REF_GENOME_NAME = 'mg1655_partial'
TEST_REF_GENOME_PATH = os.path.join(settings.PWD,
    'test_data/full_vcf_test_set/mg1655_tolC_through_zupT.gb')

TEST_FASTA_1_PATH = os.path.join(settings.PWD,
    'test_data/genome_finish_test/random_fasta_1.fa')
TEST_FASTA_2_PATH = os.path.join(settings.PWD,
    'test_data/genome_finish_test/random_fasta_2.fa')
INS_1KB_REF_GENOME_PATH = os.path.join(settings.PWD,
    'test_data/genome_finish_test/ins_1kb.fa')
INS_1KB_FQ_1_PATH = os.path.join(settings.PWD,
    'test_data/genome_finish_test/ins_1kb.1.fq')
INS_1KB_FQ_2_PATH = os.path.join(settings.PWD,
    'test_data/genome_finish_test/ins_1kb.2.fq')
INS_1KB_INSERTION_SEQUENCE_PATH = os.path.join(settings.PWD,
    'test_data/genome_finish_test/ins_1kb_insertion.fa')
INSERTION_LENGTH = 1000

class TestGenomeConcatenation(TestCase):

    def setUp(self):
        # Useful models.
        self.user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.project = Project.objects.create(owner=self.user.get_profile(),
                title='Test Project')
        self.ref_genome = ReferenceGenome.objects.create(project=self.project,
                label='refgenome')

        # Fake web browser client used to make requests.
        self.client = Client()
        self.client.login(username=TEST_USERNAME, password=TEST_PASSWORD)

    def test_fasta_concatenation(self):

        fasta1_ref=import_reference_genome_from_local_file(self.project, 'fasta1', TEST_FASTA_1_PATH,
        'fasta', move=False)
        fasta2_ref=import_reference_genome_from_local_file(self.project, 'fasta2', TEST_FASTA_2_PATH,
        'fasta', move=False)

        request_data = {
            'newGenomeLabel': 'fasta12_concat',
            'refGenomeUidList': [fasta1_ref.uid, fasta2_ref.uid]
        }

        request = HttpRequest()
        request.POST = {'data':json.dumps(request_data)}
        request.method = 'POST'
        request.user = self.user

        authenticate(username=TEST_USERNAME, password=TEST_PASSWORD)
        self.assertTrue(request.user.is_authenticated())

        response = xhr_handlers.ref_genomes_concatenate(request)

        concat_ref=ReferenceGenome.objects.get(label='fasta12_concat')

        # Assert correct number of chromosomes
        assert(concat_ref.num_chromosomes==2)

        # Assert correct number of bases
        assert(concat_ref.num_bases == sum([rg.num_bases for rg in
                [fasta1_ref,fasta2_ref]]))

    def test_generate_contigs_xhr(self):
        """Tests generate_contigs xhr request handling
        """
        # Make ReferenceGenome
        ins_1kb_ref=import_reference_genome_from_local_file(self.project, 'ins_1kb', INS_1KB_REF_GENOME_PATH,
        'fasta', move=False)

        # Make ExperimentSample
        ins_1kb_reads=ExperimentSample.objects.create(
            project=self.project,
            label='ins_1kb_reads')
        add_dataset_to_entity(ins_1kb_reads, 'reads_fq_1', Dataset.TYPE.FASTQ1,
            filesystem_location=INS_1KB_FQ_1_PATH)
        add_dataset_to_entity(ins_1kb_reads, 'reads_fq_2', Dataset.TYPE.FASTQ2,
            filesystem_location=INS_1KB_FQ_2_PATH)

        # Make Alignment group
        alignment_group = AlignmentGroup.objects.create(
            reference_genome=ins_1kb_ref)

        # Make resulting ExperimentSampleToAlignment
        reads_align = ExperimentSampleToAlignment.objects.create(
            alignment_group=alignment_group,
            experiment_sample=ins_1kb_reads)

        contig_label='contigs_0'
        request_data = {
            'contig_label': contig_label,
            'experiment_sample_uid': reads_align.uid
        }

        request = HttpRequest()
        request.GET = request_data
        request.method = 'GET'
        request.user = self.user

        authenticate(username=TEST_USERNAME, password=TEST_PASSWORD)
        self.assertTrue(request.user.is_authenticated())

        response = xhr_handlers.generate_contigs(request)

        # Assert contigs were generated
        contig_ref_label = ' :: '.join([ins_1kb_ref.label, contig_label])
        assert len(ReferenceGenome.objects.filter(label=contig_ref_label))==1

        contig_ref = ReferenceGenome.objects.get(label=contig_ref_label)
        assert contig_ref.num_bases > 0

    def test_generate_contigs_with_existing_alignment(self):
        """Tests generate_contigs xhr request handling
        """
        # Make ReferenceGenome
        ins_1kb_ref=import_reference_genome_from_local_file(self.project, 'ins_1kb', INS_1KB_REF_GENOME_PATH,
        'fasta', move=False)

        # Make ExperimentSample
        ins_1kb_reads=ExperimentSample.objects.create(
            project=self.project,
            label='ins_1kb_reads')
        add_dataset_to_entity(ins_1kb_reads, 'reads_fq_1', Dataset.TYPE.FASTQ1,
            filesystem_location=INS_1KB_FQ_1_PATH)
        add_dataset_to_entity(ins_1kb_reads, 'reads_fq_2', Dataset.TYPE.FASTQ2,
            filesystem_location=INS_1KB_FQ_2_PATH)

        # Run alignment of reads to reference
        alignment_group_label = 'fq_reads_pipeline_alignment'
        ref_genome = ins_1kb_ref
        sample_list = [ins_1kb_reads]
        alignment_group, _ = run_pipeline(alignment_group_label, ref_genome, sample_list,
            perform_variant_calling=False, alignment_options={})

        # Get resulting ExperimentSampleToAlignment
        reads_align = ExperimentSampleToAlignment.objects.get(
            alignment_group=alignment_group,
            experiment_sample=ins_1kb_reads)

        contig_label='contigs_0'
        request_data = {
            'contig_label': contig_label,
            'experiment_sample_uid': reads_align.uid
        }

        request = HttpRequest()
        request.GET = request_data
        request.method = 'GET'
        request.user = self.user

        authenticate(username=TEST_USERNAME, password=TEST_PASSWORD)
        self.assertTrue(request.user.is_authenticated())

        response = xhr_handlers.generate_contigs(request)

        # Assert contigs were generated
        contig_ref_label = ' :: '.join([ins_1kb_ref.label, contig_label])
        assert len(ReferenceGenome.objects.filter(label=contig_ref_label))==1

        contig_ref = ReferenceGenome.objects.get(label=contig_ref_label)
        assert contig_ref.num_bases > 0


    def test_1kb_insertion_detection(self):
        # Make ReferenceGenome
        ins_1kb_ref=import_reference_genome_from_local_file(self.project, 'ins_1kb', INS_1KB_REF_GENOME_PATH,
        'fasta', move=False)

        # Make ExperimentSample
        ins_1kb_reads=ExperimentSample.objects.create(project=self.project,
            label='ins_1kb_reads')
        add_dataset_to_entity(ins_1kb_reads, 'reads_fq_1', Dataset.TYPE.FASTQ1,
            filesystem_location=INS_1KB_FQ_1_PATH)
        add_dataset_to_entity(ins_1kb_reads, 'reads_fq_2', Dataset.TYPE.FASTQ2,
            filesystem_location=INS_1KB_FQ_2_PATH)

        # Run alignment of reads to reference
        alignment_group_label = 'fq_reads_alignment'
        ref_genome = ins_1kb_ref
        sample_list = [ins_1kb_reads]
        alignment_group, _ = run_pipeline(alignment_group_label, ref_genome, sample_list,
            perform_variant_calling=False, alignment_options={})

        # Get resulting ExperimentSampleToAlignment
        reads_align = ExperimentSampleToAlignment.objects.get(
            alignment_group=alignment_group,
            experiment_sample=ins_1kb_reads)


        # #DEBUG:
        # alignment_bam = get_dataset_with_type(reads_align, Dataset.TYPE.BWA_ALIGN).get_absolute_location()
        # alignment_sam = alignment_bam[:-3] + 'sam'
        # de_novo_bam.make_sam(alignment_bam)
        # with open(alignment_sam,'r') as fh:
        #     for line in fh:
        #         print line

 
        # House the contigs in a ReferenceGenome
        contigs_ref_genome = ReferenceGenome.objects.create(project=self.project,
            label='ins_contigs')

        contig_files = assembly.generate_contigs(reads_align, contigs_ref_genome)
        contigs_0 = contig_files[0]

        # Add the generated contig fastas to the contig ReferenceGenome
        add_dataset_to_entity(contigs_ref_genome, 'ins_contigs_fasta', Dataset.TYPE.REFERENCE_GENOME_FASTA,
            filesystem_location=contigs_0)

        # Create Experiment Sample for the inserted sequence
        insertion_sequence_sample = ExperimentSample.objects.create(project=self.project,
            label='insertion_sequence')

        # Convert inserted sequence to fastq for alignment
        fastq_path = os.path.join(insertion_sequence_sample.get_model_data_dir(), 'insertion.fq')
        convert_fasta_to_fastq(INS_1KB_INSERTION_SEQUENCE_PATH, fastq_path)

        # Add fastq dataset to Experiment Sample
        add_dataset_to_entity(insertion_sequence_sample, 'insertion_fastq', Dataset.TYPE.FASTQ1,
            filesystem_location=fastq_path)

        # Previous allignment code using millstone's pipeline does not work for one input fastq
            # # Align insertion against contigs
            # alignment_group_label = 'insertion_to_contigs_alignment'
            # ref_genome = contigs_ref_genome
            # sample_list = [insertion_sequence_sample]
            # alignment_group, _ = run_pipeline(alignment_group_label, ref_genome, sample_list,
            #     perform_variant_calling=False, alignment_options={})

        # Select insertion sequence fastq and contigs fasta for bwa alignment 
        reads = [get_dataset_with_type(insertion_sequence_sample, Dataset.TYPE.FASTQ1).get_absolute_location()]
        ref_fasta = get_dataset_with_type(contigs_ref_genome, Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
        data_dir = os.path.join(contigs_ref_genome.get_model_data_dir(), 'alignment_data')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        alignment_bam = bwa_align(reads, ref_fasta, data_dir)

        # Find the maximum number of matches to the insertion captured by any contig
        match_counts = get_match_counts(alignment_bam)
        max_matches = 0
        for k in match_counts:
            if max(match_counts[k])>max_matches:
                max_matches = max(match_counts[k])

        # Ensure thorough coverage of the insertion
        INSERTION_COVERAGE_CUTOFF = 0.95
        max_cov_fraction = float(max_matches)/INSERTION_LENGTH
        print 'Contigs covered %.3f of the %d base inserted sequence\n' % (max_cov_fraction, INSERTION_LENGTH)
        assert max_cov_fraction >= INSERTION_COVERAGE_CUTOFF, (
                'The maximum fraction of the insertion captured by any contig was: %.3f, \
                less than the passing cutoff: %.3f ' % (max_cov_fraction, INSERTION_COVERAGE_CUTOFF))
