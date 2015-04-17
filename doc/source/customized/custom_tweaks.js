$(document).ready(function() {
  var MakeFirstParagraphLead = function(elem) {
    var lead = $(elem).children("p:first");
    lead.addClass("lead");
    if (lead.html()) {
      lead.html(lead.html().charAt(0).toUpperCase() + lead.html().slice(1));
    }
    return elem;
  };

  var RemoveEmptyParagraphs = function(parent) {
    $(parent).find("p").each(
      function() {
        if ($(this).html() == "") {
          $(this).detach();
        }
      }
    );
    return parent;
  };

  var GetMetablock = function(elem) {
    if ($(elem).hasClass("memdoc")) {
      return $(elem).parent().find(".meta-info-block");
    } else {
      return $(elem).find(".meta-info-block");
    }
  };

  var AddToMetablock = function(elem, css) {
    var parent = $(elem).parent();
    var metablock = GetMetablock(parent);

    if (metablock.length == 0) {
      if (parent.hasClass("memdoc")) {
        parent.after("<div class='panel-footer'><ul class='meta-info-block'></ul></div>");
      } else {
        // we are in a details textblock
        $(elem).before("<ul class='meta-info-block'></ul>");
      }
    }
    metablock = GetMetablock(parent);

    $(elem).detach();
    metablock.append("<li class='" + css + "'>" + $(elem).html() + "</li>");
  };

  var TweakMemberDocStyle = function(elem) {
    var definition = $(elem).find("p").filter(function() {
      return ($(this).text()).match('^Definition at line.*\.$');
    });
    var overloaded = $(elem).find("p").filter(function() {
      return ($(this).text()).match('^This is an overloaded.*\.$');
    });
    var referenced = $(elem).find("p").filter(function() {
      return ($(this).text()).match('^Referenced by.*\.$');
    });
    var references = $(elem).find("p").filter(function() {
      return ($(this).text()).match('^References.*\.$');
    });
    var reimpl     = $(elem).find("p").filter(function() {
      return ($(this).text()).match('^Reimplemented (in|from).*\.$');
    });
    var impl       = $(elem).find("p").filter(function() {
      return ($(this).text()).match('^Implement(ed in|s).*\.$');
    });

    AddToMetablock(definition, 'meta-definition');
    AddToMetablock(overloaded, 'meta-overload');
    AddToMetablock(referenced, 'meta-referenced');
    AddToMetablock(references, 'meta-references');
    AddToMetablock(reimpl, 'meta-reimplemented');
    AddToMetablock(impl, 'meta-implemented');

    var metablock = GetMetablock($(elem));

    metablock.find("li a.el").addClass("code").removeClass("el");

    // make the first paragraph stand out
    MakeFirstParagraphLead($(elem));

    // remove empty paragraphs
    RemoveEmptyParagraphs($(elem));
  };

  var MakePanel = function(container, title, content, footer) {
    $(container).addClass("panel panel-default");
    $(container).find(title).addClass("panel-heading");
    $(container).find(content).addClass("panel-body");
    if (typeof(footer) !== 'undefined') {
      $(container).find(footer).addClass("panel-footer");
      if ($(container).find(content).text() === "" 
          && $(container).find(footer).text() === "") {
        $(container).detach();
      }
    } else {
      if ($(container).find(content).text() == "") {
        $(container).detach();
      }
    }
  };

  var TweakMemberSections = function(dl) {
    MakePanel($(dl), 'dt', 'dd');
    dl.addClass("member-section").removeClass("panel-default");
    var dt = dl.find("dt");
    switch (dt.text()) {
      case "Template Parameters":
        dl.addClass("panel-primary");
        break;
      case "Parameters":
        dl.addClass("panel-primary");
        break;
      case "Returns":
        dl.addClass("panel-success");
        break;
      case "Precondition":
        dl.addClass("panel-info");
        break;
      case "Note":
        dl.addClass("panel-info");
        break;
      case "Exceptions":
        dl.addClass("panel-danger");
        break;
      case "Todo:":
        dl.addClass("panel-warning");
        break;
      default:  // 'See also', 'Since'
        dl.addClass("panel-default");
        break;
    }

    var param_table = dl.find("table.params, table.tparams, table.exception");
    if (param_table.length > 0) {
      param_table.addClass("table");
      param_table.unwrap();
    }

    var dds = dl.find("dd");
    if (dds.length > 1) {
      dt.after("<ul class='list-group'></ul>");
      var new_dds = dt.next();
      dds.each(function() {
        new_dds.append("<li class='list-group-item'>" + $(this).html() + "</li>");
        $(this).detach();
      });
    } else {
      dds.replaceWith(function() { return "<div class='" + $(this).attr('class') + "'>" + $(this).html() + "</div>"; });
    }
    dt.replaceWith("<div class='" + dt.attr('class') + "'>" + dt.html() + "</div>");
    dl.replaceWith("<div class='" + dl.attr('class') + "'>" + dl.html() + "</div>");
  };

  // detailed description text
  $("h2.groupheader:contains('Detailed Description')").each(
    function() {
      TweakMemberDocStyle($(this).next());
    }
  );

  // member details doc
  $(".memdoc").each(function() {
    TweakMemberDocStyle($(this));
    $(this).find("dl").each(function() {
      TweakMemberSections($(this));
    });
  });

  $(".memitem").each(function() {
    MakePanel($(this), '.memproto', '.memdoc');
  });

  $(".memitem").each(function() {
    var header = $(this).children(".panel-heading");
    var footer = $(this).children(".panel-footer").find(".meta-info-block");
    if (footer.find(".meta-overload").is(".meta-overload")) {
      header.find(".mlabels td.mlabels-right > .mlabels").prepend("<span class='mlabel label label-default'>overload</span>");
    }
  });

  $(".mlabel").each(function() {
    switch ($(this).text()) {
      case "delete":
        $(this).addClass("label label-danger");
        break;
      case "private":
        $(this).addClass("label label-warning");
        break;
      case "protected":
        $(this).addClass("label label-primary");
        break;
      case "static":
        $(this).addClass("label label-success");
        break;
      case "virtual":
        $(this).addClass("label label-info");
        break;
      default:  // 'override', 'inline'
        $(this).addClass("label label-default");
    }
  });

  // directory tables
  $(".directory table.directory").each(function() {
    $(this).addClass("table");
    $(this).find("span.icon").unwrap().addClass("badge").removeClass("icon");
  });

  $("table.doxtable").addClass("table");

  // reflists (i.e. ToDos)
  $("div.textblock dl.reflist").each(function() {
    $(this).find("dt").each(function() {
      $(this).wrap("<div class='panel panel-default'></div>");
      $(this).replaceWith("<div class='panel-heading'>" + $(this).html() + "</div>");
    });
    $(this).find("dd").each(function() {
      $(this).prev().append($(this));
      $(this).replaceWith("<div class='panel-body'>" + $(this).html() + "</div>");
    });
    $(this).find(".panel-body p").removeClass();
  });

  // panels/blocks in detailed descriptions
  $("div.textblock").each(function() {
    $(this).find("dl").each(function() {
      TweakMemberSections($(this));
    });
  });
});
