from django.conf.urls import include
from django.conf.urls import patterns
from django.conf.urls import url
from django.views.generic import RedirectView

import settings

urlpatterns = patterns('',
    url(r'^$', 'genome_designer.main.views.home_view'),

    # Project-specific views
    url(r'^projects$',
            'genome_designer.main.views.project_list_view'),
    url(r'^projects/create$',
            'genome_designer.main.views.project_create_view'),
    url(r'^projects/([\w-]+)$',
            'genome_designer.main.views.project_view'),
    url(r'^projects/([\w-]+)/delete$',
            'genome_designer.main.views.project_delete'),


    # Tab base views.
    url(r'^projects/([\w-]+)/data$',
            'genome_designer.main.views.project_view'),
    url(r'^projects/([\w-]+)/align$',
            'genome_designer.main.views.tab_root_align'),
    url(r'^projects/([\w-]+)/analyze$',
            'genome_designer.main.views.tab_root_analyze'),
    url(r'^projects/([\w-]+)/analyze/([\w-]+)$',
            'genome_designer.main.views.tab_root_analyze'),
    url(r'^projects/([\w-]+)/analyze/([\w-]+)/([\w-]+)$',
            'genome_designer.main.views.tab_root_analyze'),


    # Reference genomes
    url(r'^projects/([\w-]+)/refgenomes$',
            'genome_designer.main.views.reference_genome_list_view'),
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)$',
            'genome_designer.main.views.reference_genome_view'),

    # Alignments
    url(r'^projects/([\w-]+)/alignments$',
            'genome_designer.main.views.alignment_list_view'),
    url(r'^projects/([\w-]+)/alignments/create$',
            'genome_designer.main.views.alignment_create_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)$',
            'genome_designer.main.views.alignment_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)/samplealign/([\w-]+)/error$',
            'genome_designer.main.views.sample_alignment_error_view'),


    # Variant sets
    url(r'^projects/([\w-]+)/sets$',
            'genome_designer.main.views.variant_set_list_view'),
    url(r'^projects/([\w-]+)/sets/([\w-]+)$',
            'genome_designer.main.views.variant_set_view'),


    # Samples
    url(r'^projects/([\w-]+)/samples$',
            'genome_designer.main.views.sample_list_view'),

    # Variants
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)/variants/([\w-]+)$',
            'genome_designer.main.views.single_variant_view'),


    ############################################################################
    # Templates
    ############################################################################

    url(r'^templates/sample_list_targets_template.tsv$',
            'genome_designer.main.views.sample_list_targets_template'),

    url(r'^templates/variant_set_upload_template.vcf$',
            'genome_designer.main.views.variant_set_upload_template'),


    ############################################################################
    # Auth
    ############################################################################

    # django-registration defaults (further delgates to django.contrib.auth.url)
    (r'^accounts/', include('registration.backends.simple.urls')),

    # The default behavior of registration is redirect to 'users/<username>'.
    # For now let's catch this request here and just redirect to '/'.
    (r'^users/', RedirectView.as_view(url='/')),


    ############################################################################
    # XHR Actions
    ############################################################################

    url(r'^_/sets$',
            'genome_designer.main.xhr_handlers.get_variant_set_list'),
    url(r'^_/sets/exportcsv$',
            'genome_designer.main.xhr_handlers.export_variant_set_as_csv'),
    url(r'^_/variants$',
            'genome_designer.main.xhr_handlers.get_variant_list'),
    url(r'^_/variants/modify_set_membership$',
            'genome_designer.main.xhr_handlers.modify_variant_in_set_membership'),
    url(r'^_/variants/refresh_materialized_variant_table$',
            'genome_designer.main.xhr_handlers.refresh_materialized_variant_table'),


    ############################################################################
    # Template XHR's
    # TODO: Replace this with client-side templating.
    ############################################################################

    url(r'^_/templates/variant_filter_controls$',
            'genome_designer.main.template_xhrs.variant_filter_controls'),
    url(r'^_/templates/variant_set_controls$',
            'genome_designer.main.template_xhrs.variant_set_controls'),

)

if settings.S3_ENABLED:
    urlpatterns += patterns('',
        url(r'^_/genes$',
                'genome_designer.main.xhr_handlers.get_gene_list'),
        url(r'^_/projects/([\w-]+)/refgenomes/import_s3$',
                'genome_designer.main.xhr_handlers.import_reference_genome_s3',
                name="import_reference_genome_s3"),
        url(r'^_/projects/([\w-]+)/samples/parse_targets_file_s3$',
                'genome_designer.main.xhr_handlers.parse_targets_file_s3',
                name="parse_targets_file_s3"),
        url(r'^_/projects/([\w-]+)/samples/process_sample_files_s3$',
                'genome_designer.main.xhr_handlers.process_sample_files_s3',
                name="process_sample_files_s3"),
        url(r'^s3/signature', 'genome_designer.main.xhr_uploader.handle_s3', name="s3_signature"),
        url(r'^s3/delete', 'genome_designer.main.xhr_uploader.handle_s3', name='s3_delete'),
        url(r'^s3/success', 'genome_designer.main.xhr_uploader.success', name="s3_success")
    )

if settings.RUNNING_ON_EC2:
    urlpatterns += patterns('',
        url(r'^ec2/info$', 'genome_designer.main.views.ec2_info_view', name="ec2_info")
    )

if settings.DEBUG:
    from django.conf.urls.static import static
    urlpatterns += static('jbrowse', document_root=settings.JBROWSE_ROOT, show_indexes=True)
