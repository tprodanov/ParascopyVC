library(zeallot)

WDIR <- 'ParascopyVC/reproducing/plots/'
source(file.path(WDIR, 'common.r'))
setwd('~/Data/proj/Benchmarks/')
plotdir <- file.path(getwd(), 'plots')
thresholds <- read_delim(file.path(plotdir, 'thresholds.csv'), '\t')

tools <- c(
    'gatk' = 'GATK',
    'freebayes' = 'FreeBayes',
    'deepvariant' = 'DeepVariant',
    'parascopy' = 'ParascopyVC')
eval_ty <- 'ref_pscn'
version <- 'v1.10.6'
samples <- sprintf('HG00%d', 1:7)

c(curves, metadata) %<-% load_curves(
    '%s/analysis/cmp/diploid-%s/%s/eval-%s',
    sample = samples, version = version, eval_ty = eval_ty, tool = tools)
curves <- add_f_metrics(curves, c(1/2, 1)) |>
    select(!c('version', 'eval_ty', 'qual', 'f1'))

tab <- left_join(curves, thresholds) |> filter(score >= thresh) |>
    group_by(sample, tool) |>
    slice_min(score, n = 1, with_ties = F)

latex_tab <- local({
    df <- inner_join(
        select(metadata, sample, baseline_vars) |> unique(),
        select(tab, sample, tool, precision, recall, F_0.5, F_1) |>
            pivot_wider(
                names_from = tool,
                values_from = c(precision, recall, F_0.5, F_1)),
        by = 'sample')
    colnames <- expand.grid(
        c('precision', 'recall', 'F_0.5', 'F_1'),
        c('GATK', 'FreeBayes', 'DeepVariant', 'ParascopyVC'))
    df[, c('sample', 'baseline_vars',
        sprintf('%s_%s', colnames$Var1, colnames$Var2))]
})
latex_tab <- rbind(latex_tab,
    apply(latex_tab, 2, function(x) mean(as.numeric(x))))
rownames(latex_tab) <- NULL

print(xtable::xtable(latex_tab, type = 'latex', digits = 3),
    file = file.path(plotdir, 'tab2.tex'))

filter(tab, sample == 'HG002')

####################

metadata
range(unique(metadata$sum_len) / 1e6)
range(unique(metadata$baseline_vars))
latex_tab0 <- latex_tab[1:7,]
with(latex_tab0, precision_ParascopyVC - precision_GATK) |> mean() |> round(3)
with(latex_tab0, precision_ParascopyVC - precision_FreeBayes) |> mean() |> round(3)
with(latex_tab0, precision_ParascopyVC - precision_DeepVariant) |> mean() |> round(3)

with(latex_tab0, recall_ParascopyVC - recall_GATK) |> mean() |> round(3)
with(latex_tab0, recall_ParascopyVC - recall_FreeBayes) |> mean() |> round(3)
with(latex_tab0, recall_ParascopyVC - recall_DeepVariant) |> mean() |> round(3)

(av <- tail(latex_tab, n = 1)) |> select(!sample) |> round(3)

filter(tab, sample == 'HG002') # |> as.data.frame()
filter(curves, sample == 'HG002' & tool == 'ParascopyVC'
    & precision >= 0.983 & recall >= 0.861)

####################

(sfig1 <- ggplot(
        filter(curves, sample != 'HG002'),
        aes(recall, precision, color = tool, fill = tool)) +
    geom_path(size = 1) +
    geom_point(data = filter(tab, sample != 'HG002'),
        size = 2, alpha = 0.8, shape = 21, color = 'gray20') +
    facet_wrap(~ sample, ncol = 2) +

    coord_cartesian(xlim = c(0.55, 1), ylim = c(0.8, 1)) +
    scale_x_continuous('Recall', expand = expansion(add = 0),
        breaks = scales::pretty_breaks(4),
        labels = function(x) sprintf('%.2g', x)) +
    scale_y_continuous('Precision',
        breaks = scales::pretty_breaks(3)) +
    scale_color_manual(NULL, values = tool_colors) +
    scale_fill_manual(NULL, values = tool_colors) +

    guides(
        color = guide_legend(override.aes = list(shape = 21), order = 1),
        fill = guide_legend(order = 1),
        shape = guide_legend(override.aes = list(color = 'black'), order = 2),
    ) +
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(margin = margin(t = 0)),
        axis.text.y = element_text(margin = margin(r = 0)),

        legend.position = 'top',
        legend.margin = margin(b = -7),
        legend.box = 'vertical',
        legend.box.margin = margin(l = 0, b = 0),
        legend.spacing.y = unit(4, 'pt'),
        legend.title = element_text(size = 9),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.text = element_text(size = 9, margin = margin(l = -3)),

        plot.margin = margin(l = 3, r = 8),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(margin = margin(t = 1, b = 1), face = 'bold'),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'gray50'),
        panel.spacing.y = unit(0.3, "lines")
    ))
ggsave(file.path(plotdir, '6wgs.pdf'),
    width = 10, height = 9, dpi = 400, scale = 0.45,
    title = '6 WGS samples')
rm(sfig1)

####################

(fig3a <- ggplot(
        filter(curves, sample == 'HG002'),
        aes(recall, precision, color = tool, fill = tool)) +
    geom_path(size = 1) +
    geom_point(data = filter(tab, sample == 'HG002'),
        size = 2, alpha = 0.8, shape = 21, color = 'gray20') +

    coord_cartesian(xlim = c(0.58, 1), ylim = c(NA, 1)) +
    scale_x_continuous('Recall', expand = expansion(add = 0),
        labels = function(x) sprintf('%.2g', x),
        breaks = scales::pretty_breaks(4)) +
    scale_y_continuous('Precision', breaks = scales::pretty_breaks(4)) +
    scale_color_manual(NULL, values = tool_colors) +
    scale_fill_manual(NULL, values = tool_colors) +

    guides(
        color = guide_legend(override.aes = list(shape = 21), order = 1),
        fill = guide_legend(order = 1),
        shape = guide_legend(override.aes = list(color = 'black'), order = 2),
    ) +
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),

        legend.position = 'top',
        legend.margin = margin(b = -7),
        legend.box = 'vertical',
        legend.box.margin = margin(l = 0, b = -5),
        legend.spacing.y = unit(4, 'pt'),
        legend.title = element_text(size = 9),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.text = element_text(size = 9, margin = margin(l = -3)),

        plot.margin = margin(l = 3, r = 8),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(margin = margin(t = 1, b = 1), face = 'bold'),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'gray50'),
    ))

#######################

similarity <- c(
    '0.97-0.98' = '0.97 - 0.98',
    '0.98-0.99' = '0.98 - 0.99',
    '0.99-0.995' = '0.99 - 0.995',
    '0.995-1.0' = '0.995 - 1.0',
    '0.99-1.0' = '0.99 - 1.0')
c(curves_ss, metadata) %<-% load_curves(
    '%s/analysis/cmp/diploid-v1.10.6-simil/%s/%s/eval-%s',
    sample = 'HG002', simil = similarity, eval_ty = eval_ty, tool = tools)
tab_ss <- left_join(curves_ss, thresholds) |> filter(score >= thresh) |>
    group_by(sample, tool, simil) |>
    slice_min(score, n = 1, with_ties = T)

(fig3b <- ggplot(filter(curves_ss, simil != '0.99 - 1.0'),
        aes(recall, precision, color = tool, fill = tool)) +
    # geom_point(aes(color = NA), # Dummy data (set y-limits to SIM-R)
    #     data = data.frame(recall = 0.8, precision = 0.75),
    #     alpha = 0) +
    geom_path(size = 1) +
    geom_point(data = filter(tab_ss, simil != '0.99 - 1.0'),
        shape = 21, size = 2, alpha = 0.8, color = 'gray20') +
    facet_wrap(~ sprintf('Similarity:  %s', simil)) +

    coord_cartesian(xlim = c(0.28, 1), ylim = c(0.78, 1.0)) +
    scale_x_continuous('Recall', expand = expansion(add = 0),
        labels = function(x) sprintf('%.2g', x)) +
    scale_y_continuous('Precision', expand = expansion(mult = 0.05),
        breaks = scales::pretty_breaks(3)) +
    scale_color_manual('Variant caller', values = tool_colors) +
    scale_fill_manual('Variant caller', values = tool_colors) +

    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(margin = margin()),
        axis.text.y = element_text(margin = margin(r = -1)),

        legend.position = 'none',
        plot.margin = margin(t = 7, l = 3, r = 8),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(margin = margin(t = 1, b = 1), face = 'bold'),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'gray50'),
    ))

cowplot::plot_grid(fig3a, fig3b,
    ncol = 1,
    labels = c('(a)', '(b)'),
    rel_heights = c(0.9, 1),
    label_size = 10, label_x = -0.01, label_y = c(1.01, 1.01)
    )
ggsave(file.path(plotdir, 'fig3.pdf'),
    width = 10, height = 11, dpi = 400, scale = 0.45,
    title = 'Figure 3')

metadata
filter(tab_ss, simil == '0.99 - 1.0') |>
    ungroup() |>
    select(!c(qual, sample, eval_ty, thresh, simil))

##############

var <- c('snp' = 'SNPs', 'non_snp' = 'Indels')
c(curves, metadata) %<-% load_curves(
    '%s/analysis/cmp/diploid-%s/%s/eval-%s/%s_roc.tsv.gz',
    sample = 'HG002', version = version, eval_ty = eval_ty, tool = tools,
    var = var)
curves <- add_f_metrics(curves, c(1/2, 1)) |>
    select(!c('version', 'eval_ty', 'qual', 'f1'))

tab <- left_join(curves, thresholds) |> filter(score >= thresh) |>
    group_by(sample, tool, var) |>
    slice_min(score, n = 1, with_ties = T) |>
    add_f_metrics(c(0.5, 1))

filter(tab, var == 'Indels') |> as.data.frame()
tab[, c('var', 'tool', 'fp', 'fn', 'precision', 'recall', 'F_0.5', 'F_1')] |>
    arrange(var, tool)
metadata

#######################

c(curves, metadata) %<-% load_curves(
    'chr15-17/HG002/analysis/cmp/diploid-%s/%s/eval-%s',
    version = version, eval_ty = eval_ty, tool = tools)
curves <- select(curves, !c('version', 'eval_ty', 'qual'))

tab <- left_join(curves, thresholds) |> filter(score >= thresh) |>
    group_by(tool) |>
    slice_min(score, n = 1, with_ties = F) |>
    ungroup()

tab[, c('tool', 'fp', 'fn', 'precision', 'recall', 'f1')]
metadata
