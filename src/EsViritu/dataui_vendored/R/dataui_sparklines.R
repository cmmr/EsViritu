# Vendored sparkline-relevant subset of the dataui R package.
# Source: https://github.com/timelyportfolio/dataui (master branch, commit circa 2020)
# License: MIT
# Vendored for EsViritu to remove the dependency on a GitHub-only package.
# Only the functions required by reactablefmtr::react_sparkline() are included.

# ---------------------------------------------------------------------------
# From R/utils.R
# ---------------------------------------------------------------------------

#' @keywords internal
wrap_components <- function(components) {
  if (!is.null(components) && length(components) > 0 && !is.list(components[[1]])) {
    components <- list(components)
  }
  components
}

# ---------------------------------------------------------------------------
# From R/components.R  (sparkline components only)
# ---------------------------------------------------------------------------

#' Line Series for 'Sparklines'
#'
#' @param ... see data-ui series
#' @return reactR component
#' @export
dui_sparklineseries <- function(...) {
  reactR::React$SparklineLineSeries(...)
}

#' Bar Series for 'Sparklines'
#'
#' @param ... see data-ui series
#' @return reactR component
#' @export
dui_sparkbarseries <- function(...) {
  reactR::React$SparklineBarSeries(...)
}

#' Point Series for 'Sparklines'
#'
#' @param ... see data-ui series
#' @return reactR component
#' @export
dui_sparkpointseries <- function(...) {
  reactR::React$SparklinePointSeries(...)
}

#' Horizontal Reference Line for 'Sparklines'
#'
#' @param ... see data-ui reference lines
#' @return reactR component
#' @export
dui_sparkrefline <- function(...) {
  reactR::React$HorizontalReferenceLine(...)
}

#' Vertical Reference Line for 'Sparklines'
#'
#' @param ... see data-ui reference lines
#' @return reactR component
#' @export
dui_sparkverticalrefline <- function(...) {
  reactR::React$VerticalReferenceLine(...)
}

#' Horizontal Reference Line for 'Sparklines' (alias)
#'
#' @param ... see data-ui reference lines
#' @return reactR component
#' @export
dui_sparkhorizontalrefline <- function(...) {
  reactR::React$HorizontalReferenceLine(...)
}

#' Band Lines for 'Sparklines'
#'
#' @param ... see data-ui reference bands
#' @return reactR component
#' @export
dui_sparkbandline <- function(...) {
  reactR::React$BandLine(...)
}

#' Pattern Fill for 'Sparklines'
#'
#' @param ... see data-ui patterns
#' @return reactR component
#' @export
dui_sparkpatternlines <- function(...) {
  reactR::React$PatternLines(...)
}

#' Linear Gradient for 'Sparklines'
#'
#' @param ... see data-ui gradients
#' @return reactR component
#' @export
dui_sparklineargradient <- function(...) {
  reactR::React$LinearGradient(...)
}

#' Labels for 'Sparklines'
#'
#' @param ... see data-ui label
#' @return reactR component
#' @export
dui_sparklabel <- function(...) {
  reactR::React$Label(...)
}

#' Tooltip Container for 'Sparklines'
#'
#' @param components list of children to include in the tooltip
#' @return reactR component
#' @export
dui_tooltip <- function(components) {
  component <- reactR::React$TooltipComponent()
  components <- wrap_components(components)
  component$children <- components
  component
}

# ---------------------------------------------------------------------------
# From R/duisparkline.R
# ---------------------------------------------------------------------------

#' 'data-ui Sparklines'
#'
#' @param data vector, list, or data.frame of data
#' @param className character css class name
#' @param margin list of the form list(top=, right=, bottom=, left=)
#' @param max,min numeric maximum/minimum for y axis
#' @param onMouseMove,onMouseLeave JS function for mouse events
#' @param renderTooltip JS function for tooltip rendering
#' @param preserveAspectRatio character svg preserveAspectRatio attribute
#' @param valueAccessor JS function for y value accessor
#' @param viewBox character svg viewBox attribute
#' @param ariaLabel character accessibility label
#' @param components list of child components
#' @param responsive logical, default TRUE
#' @param width,height numeric or css size
#' @param elementId character css identifier
#'
#' @import htmlwidgets
#' @return react htmlwidget
#' @export
dui_sparkline <- function(
  data = NULL,
  className = NULL,
  margin = NULL,
  max = NULL,
  min = NULL,
  onMouseMove = NULL,
  onMouseLeave = NULL,
  renderTooltip = NULL,
  preserveAspectRatio = NULL,
  valueAccessor = NULL,
  viewBox = NULL,
  ariaLabel = NULL,
  components = list(),
  responsive = TRUE,
  width = "100%", height = 100, elementId = NULL
) {
  if (responsive == TRUE || is.null(width)) {
    tagname <- "SparklineResponsive"
  } else {
    tagname <- "SparklineWithTooltip"
  }
  if (responsive == FALSE && !is.numeric(width)) {
    message(paste0(
      "Responsive is FALSE which often means width should be numeric to work properly.  ",
      paste0("In this case we see responsive = FALSE and width = ", width, ".  "),
      "Try a different width if the sparkline does not appear.",
      collapse = "\n"
    ))
  }
  component <- reactR::component(
    tagname,
    Filter(Negate(is.null), list(
      data = data,
      className = className,
      margin = margin,
      max = max,
      min = min,
      onMouseMove = onMouseMove,
      onMouseLeave = onMouseLeave,
      renderTooltip = renderTooltip,
      preserveAspectRatio = preserveAspectRatio,
      valueAccessor = valueAccessor,
      viewBox = viewBox,
      ariaLabel = ariaLabel,
      height = height,
      width = width
    ))
  )

  components <- wrap_components(components)

  hw <- htmlwidgets::createWidget(
    name = 'dataui',
    reactR::reactMarkup(component),
    width = width,
    height = height,
    package = 'dataui',
    elementId = elementId
  )

  hw$x$tag$children <- components
  hw
}

# ---------------------------------------------------------------------------
# From R/dependency.R
# ---------------------------------------------------------------------------

#' 'JavaScript' Dependencies for Standalone Usage
#'
#' @return htmltools::htmlDependency
#' @export
html_dependency_dataui <- function() {
  htmltools::htmlDependency(
    name = 'data-ui',
    version = '0.8.4',
    src = c(file = system.file("www", package = "dataui")),
    script = "dataui.standalone.js"
  )
}

# ---------------------------------------------------------------------------
# From R/reactable_helpers.R
# ---------------------------------------------------------------------------

#' @keywords internal
recurse <- function(l, func, ...) {
  l <- func(l, ...)
  if (is.list(l) && length(l) > 0) {
    lapply(
      l,
      function(ll) {
        recurse(ll, func, ...)
      }
    )
  } else {
    l
  }
}

#' @keywords internal
enclose_jseval <- function(l) {
  if (inherits(l, "JS_EVAL")) {
    l <- glue::glue("{{{{!!{code}!!}}}}", code = l)
  }
  return(l)
}

#' Convert 'dataui' to 'reactable' Custom Render Function
#'
#' @param dui dataui htmlwidget to convert
#' @param jsarg character argument name for the JS function
#' @return htmlwidgets::JS of class JS_EVAL
#' @export
dui_for_reactable <- function(dui, jsarg = "cellInfo") {
  dui <- recurse(dui, enclose_jseval)
  dui_js <- gsub(
    x = jsonlite::toJSON(
      dui$x$tag,
      force = TRUE,
      auto_unbox = TRUE
    ),
    pattern = "(\"\\{\\{!!)|(!!\\}\\}\")",
    replacement = "",
    perl = TRUE
  )

  message("
For this to work, please add the dataui dependency to your reactable instance.
dui_add_reactable_dep(...reactable widget...)
  ")

  htmlwidgets::JS(htmltools::HTML(
    gsub(
      x = glue::glue(
"
function({jsarg}) {{
  var spk = {dui_js}
  return window.reactR.hydrate(dataui, spk);
}}
",
        jsarg = jsarg,
        dui_js = dui_js
      ),
      pattern = "(\\\\n)",
      replacement = "\n"
    )
  ))
}

#' Add 'dataui' Dependency to 'reactable'
#'
#' @param rt reactable htmlwidget
#' @return reactable htmlwidget with dataui dependency attached
#' @export
dui_add_reactable_dep <- function(rt) {
  if (!inherits(rt, "reactable") && !inherits(rt, "htmlwidget")) {
    stop("Please supply a reactable or htmlwidget as the rt argument to this function.")
  }
  dep <- rt$dependencies
  if (!is.list(dep)) dep <- list()
  dep[[length(dep) + 1]] <- html_dependency_dataui()
  rt$dependencies <- dep
  rt
}
