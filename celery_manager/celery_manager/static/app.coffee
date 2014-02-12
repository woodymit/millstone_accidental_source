@CeleryManager = @CeleryManager || {}
@CeleryManager.config = @CeleryManager.config || {}
@CeleryManager.classes = @CeleryManager.classes || {}

class CeleryManager.HomeView extends Backbone.View
  events:
    "click .test-button": "test"
  initialize: =>
    @template = _.template $("#home-view-template").text()
    @render()

  success: (fieldset, message) =>
    alert message
    fieldset.parents(".panel").removeClass("panel-danger").addClass("panel-success")

  error: (fieldset, message) =>
    alert message
    fieldset.parents(".panel").removeClass("panel-success").addClass("panel-danger")

  test: (e) =>
    fieldset = $(e.target).parents("fieldset")
    $.ajax
      url: fieldset.data("test-url")
      type: 'post'
      data: fieldset.serialize()
      success: (data, textStatus, jqXHR) =>
        if data.status == "error"
          @error fieldset, data.message
        else
          @success fieldset, "Passed test!"
      beforeSend: ->
        fieldset.block(message: null, css: {opacity: 0}).spin()
      complete: ->
        fieldset.unblock().spin false
    return false

  render: =>
    @$el.append @template()

    jQuery.validator.addMethod("domainOrIpv4", (value, element, param) ->
        ipv4regex = /^(25[0-5]|2[0-4]\d|[01]?\d\d?)\.(25[0-5]|2[0-4]\d|[01]?\d\d?)\.(25[0-5]|2[0-4]\d|[01]?\d\d?)\.(25[0-5]|2[0-4]\d|[01]?\d\d?)$/i
        ec2regex = /^ec2-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)-(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.compute-1\.amazonaws\.com$/i
        return this.optional(element) || ec2regex.test(value) || ipv4regex.test(value)
      , "Please enter a valid IPv4 or domain name of EC2 instance.")
    @$("form").validate rules:
      host: 
        required: true,
        domainOrIpv4: true
      user: 
        required: true,
      password: 
        required: false

class CeleryManager.Router extends Backbone.Router
  routes:
    ""  :   "home"

  initialize: =>
    jQuery.validator.setDefaults
      errorPlacement: (error, element) ->
        if(element.parent().hasClass('input-prepend') || element.parent().hasClass('input-append'))
          error.insertAfter(element.parent())
        else
          error.insertAfter(element)
      errorElement: "small"
      wrapper: "div"
      highlight: (element) ->
        $(element).closest('.form-group').addClass('has-error')
      success: (element) ->
        $(element).closest('.form-group').removeClass('has-error')

    CeleryManager.config.urlRoot = window.location.href.split("?")[0]

  home: =>
    @currentView = new CeleryManager.HomeView el: $("#startup-view")

$ ->
  CeleryManager.app = new CeleryManager.Router()
  Backbone.history.start()
