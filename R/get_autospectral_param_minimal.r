# get_autospectral_param_minimal.r


#' Get Minimal Autospectral Parameters
#'
#' This function returns a minimal set of autospectral parameters.
#'
#' @title Get Minimal Autospectral Parameters
#' @description Returns a minimal set of autospectral parameters.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom parallelly availableCores
#' @return A list of minimal autospectral parameters.
#' @export


get.autospectral.param.minimal <- function()
{
    color.pal <- brewer.pal( 9, "Set1" )

    list(

      ### cytometer parameters
      # these may be updated by calling a specific cytometer
      cytometer = NULL,

      scatter.data.min.x = NULL,
      scatter.data.max.x = NULL,
      scatter.data.min.y = NULL,
      scatter.data.max.y = NULL,

      expr.data.min = NULL,
      expr.data.max = NULL,

      default.scatter.parameter = NULL,
      default.time.parameter = "Time",
      default.transformation.param = NULL,

      non.spectral.channel = NULL,

      af.channel = NULL,

      data.step = NULL,

      viability.gate.scaling = 0.67,

      large.gate.quantile = 0.25,
      large.gate.scaling.x = 3,
      large.gate.scaling.y = 6,

      ### general parameters
      verbose = TRUE,

      parallel = FALSE,

      worker.process.n = availableCores() - 1,

      max.memory.n = 2 * 1024^3,

      antigen.autof = "AF",

      marker.forbidden.char = " !\"#$%&'*,/:;?@[\\]^{|}~",
      marker.substitution.char = "-",

      database.dir = "./source/fluorophore_database/",

      similarity.warning.n = 0.95,

      ### autofluorescence and control cleaning parameters
      # peacoqc
      peacoqc.method = "MAD",

      # Autofluorescence gating for removal
      af.density.threshold = 0.75,

      af.gate.param = list(
        density.threshold = 0.001,
        region.auto = TRUE ),

      af.figure.gate.scale.expand = 0.01,

      af.gate.target.max = 1,
      af.gate.density.bw.factor = 3,
      af.gate.bound.density.neigh.size = 3,
      af.gate.bound.density.grid.n = 150,
      af.gate.bound.strict = TRUE,

      af.spline.x.bound.factor.low = 0.05,
      af.spline.x.bound.factor.high = 0.95,
      af.spline.y.bound.factor.low = 0.10,
      af.spline.y.bound.factor.high = 0.95,
      af.spline.maxit = 100,
      af.spline.sd.n = 2,

      af.remove.pop = 1,

      af.plot.bw.factor = 5,
      af.plot.density.grid.n = 100,
      af.plot.define.filename = "AF identification",
      af.plot.filename = "AF removal",

      # universal negative selection parameters
      positivity.threshold = 0.995,
      positivity.threshold.af = 0.95,
      scatter.match.threshold = 0.8,
      negative.n = 1000,
      positive.n = 500,
      scatter.match.plot.width = 12,
      scatter.match.plot.height = 6,
      scatter.match.plot.filename = "universal negative scatter plot.jpg",
      scatter.match.plot.text.size = 15,
      scatter.match.plot.text.face = "bold",

      ### settings to control what happens if you don't have enough data
      min.cell.warning.n = 500,
      min.cell.stop.n = 50,

      ### spectral ribbon plot parameters
      ribbon.plot.min = -1e3,
      ribbon.breaks = NULL,
      ribbon.scale.colors = c( NA, "#440154FF", "#238A8DFF", "#55C667FF",
                               "#B8DE29FF", "#FDE725FF" ),
      ribbon.scale.values = c( 0, 0.1, 0.2, 0.3, 0.4, 1 ),
      ribbon.bins = 200,
      ribbon.plot.width = 15,
      ribbon.plot.height = 10,
      ribbon.plot.filename = "spectral ribbon plot.jpg",
      ribbon.plot.axis.text.angle = 45,
      ribbon.plot.strip.text.size = 15,
      ribbon.plot.strip.text.face = "bold",

      ### gating parameters
      # gate parameters for cells
      default.gate.param.cells = list(
        density.threshold = 0.05,
        region.auto = TRUE,
        region.factor.x.low = 0.05,
        region.factor.x.high = 0.80,
        region.factor.y.low = 0.05,
        region.factor.y.high = 0.80
      ),

      gate.data.trim.factor.x.min.cells = 0.01,
      gate.data.trim.factor.x.max.cells = 0.99,
      gate.data.trim.factor.y.min.cells = 0.01,
      gate.data.trim.factor.y.max.cells = 0.99,

      gate.bound.density.bw.factor.cells = 3.0,
      gate.bound.density.grid.n.cells = 100,
      gate.bound.density.neigh.size.cells = 3,

      gate.bound.density.max.target.cells = 1,
      gate.bound.density.max.exclusion.x.cells = 0.1,
      gate.bound.density.max.exclusion.y.cells = 0.05,
      gate.bound.density.max.mad.factor.cells = 2.0,

      gate.region.density.bw.factor.cells = 2.0,
      gate.region.density.grid.n.cells = 100,
      gate.region.density.neigh.size.cells = 2,

      gate.region.max.density.bw.factor.cells = 1.0,
      gate.region.max.density.grid.n.cells = 100,
      gate.downsample.n.cells = 100000,

      # gate parameters for beads
      default.gate.param.beads = list(
        density.threshold = 0.05,
        region.auto = TRUE,
        region.factor.x.low = 0.02,
        region.factor.x.high = 0.50,
        region.factor.y.low = 0.02,
        region.factor.y.high = 0.50
      ),

      gate.data.trim.factor.x.min.beads = 0.01,
      gate.data.trim.factor.x.max.beads = 0.99,
      gate.data.trim.factor.y.min.beads = 0.01,
      gate.data.trim.factor.y.max.beads = 0.99,

      gate.bound.density.bw.factor.beads = 3.0,
      gate.bound.density.grid.n.beads = 100,
      gate.bound.density.neigh.size.beads = 3,

      gate.bound.density.max.target.beads = 1,
      gate.bound.density.max.exclusion.x.beads = 0.1,
      gate.bound.density.max.exclusion.y.beads = 0.05,
      gate.bound.density.max.mad.factor.beads = 3.0,

      gate.region.density.bw.factor.beads = 2.0,
      gate.region.density.grid.n.beads = 100,
      gate.region.density.neigh.size.beads = 2,

      gate.region.max.density.bw.factor.beads = 1.0,
      gate.region.max.density.grid.n.beads = 100,
      gate.downsample.n.beads = 10000,

      ### refine spillover (unmixing) parameters
      rlm.iter.max = 100,
      rlm.trim.factor = 0.003,
      rlm.downsample.n = 25000,

      rs.iter.max = 100,

      rs.lambda.coarse = 1.0,
      rs.lambda.fine = 0.1,

      rs.delta.history.n = 10,

      rs.delta.threshold.untr = 1e-2,
      rs.delta.threshold.tran = 1e-4,
      rs.delta.threshold.change = 1e-6,

      convergence.color.delta = color.pal[ 7 ],    # brown
      convergence.color.delta.max = color.pal[ 5 ],    # orange
      convergence.color.delta.change = color.pal[ 8 ],    # pink
      convergence.shape.linear = "triangle",
      convergence.shape.biexp = "circle",
      convergence.shape.posnegpop = "triangle open",

      ### fix my unmix parameters
      fix.iter.max = 20,
      fix.downsample.n = 30000,
      fix.positivity.threshold = 0.99,

      fix.spectra.filename = "fixed_spectra",
      fix.spillover.filename = "fixed_spillover",
      fix.compensation.filename = "fixed_compensation",

      ### directory parameters
      unmixed.fcs.dir = "AutoSpectral_unmixed",

      figure.scatter.dir.base = NULL,

      figure.gate.dir = NULL,
      figure.af.dir = NULL,
      figure.peacoqc.dir = NULL,
      figure.clean.control.dir = NULL,
      figure.spectral.ribbon.dir = NULL,
      figure.convergence.dir = NULL,
      figure.spectra.dir = NULL,
      figure.slope.error.dir = NULL,
      figure.skewness.dir = NULL,
      figure.similarity.heatmap.dir = NULL,

      table.convergence.dir = NULL,
      table.spectra.dir = NULL,
      table.slope.error.dir = NULL,
      table.skewness.dir = NULL,

      ### filename parameters
      # files you can use to load in information
      marker.file.name = NULL,
      gate.parameter.file.name = "fcs_gate_parameter.csv",
      scatter.parameter.file.name = "fcs_scatter_parameter.csv",
      transformation.parameter.file.name = "fcs_transformation_parameter.csv",

      # how the output files will be called
      convergence.file.name = "autospectral_convergence",
      af.file.name = "autospectral autofluorescence",
      spectra.file.name = "autospectral spectra",
      slope.error.file.name = "autospectral_slope_error",
      skewness.file.name = "autospectral_skewness",
      similarity.heatmap.file.name = "autospectral_similarity_matrix",
      ssm.heatmap.file.name = "autospectral_spread_matrix",

      ### plotting parameters
      # color palette for dot plots
      density.color.single = "blue3",
      density.color.initial = color.pal[ 3 ],    # green
      density.color.final = color.pal[ 2 ],    # blue
      density.color.posnegpop = color.pal[ 1 ],    # red
      density.palette.n = 1000,
      density.palette.base.n = 1000000,
      density.palette.base.color = c( "blue", "cyan", "green", "yellow", "red" ),

      # gating plot figure parameters
      gate.tesselation.color = "blue3",
      gate.downsample.seed = 5000,

      # parameters for main figures
      figure.width = 8.0,
      figure.height = 6.0,
      figure.margin = 4.0,

      figure.panel.line.size = 0.5,

      figure.axis.text.size = 12.0,
      figure.axis.title.size = 12.0,

      figure.convergence.point.size = 2.0,
      figure.convergence.line.size = 0.8,

      figure.density.line.size = 0.2,

      figure.gate.scale.expand = 0.01,
      figure.gate.point.size = 0.8,
      figure.gate.line.size = 0.5,
      figure.gate.bar.width = 1.0,
      figure.gate.bar.height = 25.0,
      figure.gate.bar.margin = 2.0,

      figure.matrix.point.size = 2.5,
      figure.matrix.line.size = 0.8,

      figure.scatter.alpha.gate.in = 0.8,
      figure.scatter.alpha.gate.out = 0.1,
      figure.scatter.point.size = 0.8,
      figure.scatter.line.size = 0.6,
      figure.scatter.error.label.size = 4.0,
      figure.scatter.error.label.pos.x = 0.90,
      figure.scatter.error.label.pos.y = 0.05,
      figure.scatter.axis.text.size = 12.0,
      figure.scatter.axis.title.size = 12.0,

      figure.spectra.line.size = 1,
      figure.spectra.point.size = 1

    )
}

