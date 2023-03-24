library(zeallot)

WDIR <- '~/Code/ParascopyVC-extra/plots/'
source(file.path(WDIR, 'common.r'))
setwd('~/Data/proj/Benchmarks/')
plotdir <- file.path(getwd(), 'plots')

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
curves <- select(curves, !c('version', 'eval_ty'))
rm(metadata)

curves1 <- local({
    df <- data.frame()
    for (i in 0:1) {
        filter(curves, score >= i) |>
            group_by(sample, tool) |>
            slice_min(score, n = 1, with_ties = F) |>
            ungroup()
    }
})

curves1 <- lapply(0:60, function(i) filter(curves, score >= i) |>
    group_by(sample, tool) |>
    slice_min(score, n = 1, with_ties = F) |>
    ungroup() |>
    mutate(thresh = i)
) %>% do.call(rbind, .) |>
    select(!c(score, qual))
curves1 <- add_f_metrics(curves1) |>
    mutate(
        YoudenJ = precision + recall - 1,
        GMean = sqrt(precision * recall)
    )

means <- group_by(curves1, tool, thresh) |>
    summarize(
        F_1 = mean(f1),
        F_2 = mean(F_2),
        #F_3 = mean(F_3),
        # F_4 = mean(F_4),
        # F_0.25 = mean(F_0.25),
        F_0.333 = mean(F_0.333),
        F_0.5 = mean(F_0.5),
        YoudenJ = mean(YoudenJ),
        GMean = mean(GMean),
        .groups = 'keep',
    ) |>
    pivot_longer(!c(tool, thresh), names_to = 'metric', values_to = 'value')

filter(means, metric == 'F_0.5') |>
    group_by(tool) |>
    slice_max(value, n = 1)

filter(means, metric == 'F_1') |>
    group_by(tool) |>
    slice_max(value, n = 1)

ggplot(means, aes(thresh, value, color = tool)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2,
        data = group_by(means, tool, metric) |> slice_max(value, with_ties = F)) +
    facet_wrap(~ metric, scales = 'free_x') +
    scale_x_continuous('Threshold', breaks = seq(0, 60, 10)) +
    scale_y_continuous('Value', breaks = seq(0, 1, 0.25)) +
    scale_color_manual('Variant caller', values = tool_colors) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = 'top',
    )
ggsave(file.path(plotdir, 'thresh_selection.png'), width = 8, height = 6, dpi = 400)

filter(means, metric == 'F_0.5') |>
    group_by(tool) |>
    slice_max(value, n = 1, with_ties = F) |>
    select(tool, thresh) |>
    write_delim(file.path(plotdir, 'thresholds.csv'), '\t')
