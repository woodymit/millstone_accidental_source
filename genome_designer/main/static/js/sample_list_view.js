/**
 * @fileoverview View for samples, including upload.
 */


gd.SampleListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {
    $('#gd-sidenav-link-samples').addClass('active');

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-sample-list-view-datatable-hook'),
        serverTarget: '/_/samples',
        controlsTemplate: '/_/templates/sample_list_controls',
        requestData: {projectUid: this.model.get('uid')},
    });

    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  decorateControls: function() {
    this.samplesControlComponent = new gd.SamplesControlsComponent({
      el: '#gd-sample-list-view-datatable-hook-control'
    });
  }
});
