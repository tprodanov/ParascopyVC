library(zeallot)
library(scales)
library(cowplot)

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
database <- c('v1.10.6' = 'EUR', 'v1.10.6-EAS' = 'EAS')
version <- 'v1.10.6'
samples <- c('008' = 'SIM-Fixed', '108' = 'SIM-Poly')

c(curves, metadata) %<-% load_curves(
    'sim%s/analysis/cmp/diploid-%s/%s/eval-%s',
    sample = samples, version = version, eval_ty = eval_ty, tool = tools)
c(curves_eas, metadata_eas) %<-% load_curves(
    'sim%s/analysis/cmp/diploid-%s-EAS/%s/eval-%s',
    sample = samples[2], version = version, eval_ty = eval_ty, tool = tools[4])
curves <- rbind(
    mutate(curves, database = 'EUR'),
    mutate(curves_eas, database = 'EAS')
)
curves <- select(curves, !c('version', 'eval_ty', 'qual'))
tab <- left_join(curves, thresholds) |> filter(score >= thresh) |>
    group_by(sample, tool, database) |>
    slice_min(score, n = 1, with_ties = F) |>
    ungroup()

plots <- list()
for (s in samples) {
    plots[[s]] <- ggplot(
            filter(curves, sample == s),
            aes(recall, precision, color = tool, fill = tool)) +
        geom_point(aes(color = NA, fill = NA), # Dummy data (set y-limits)
            data = data.frame(recall = 0.8, precision = 0.97),
            alpha = 0) +
        geom_path(aes(linetype = database), linewidth = 1) +
        geom_point(data = filter(tab, sample == s),
            shape = 21, size = 2, alpha = 0.8, color = 'gray20') +

        ggtitle(s) +
        coord_cartesian(xlim = c(0.6, 1), ylim = c(NA, 1)) +
        scale_x_continuous('Recall', expand = expansion(add = 0),
            breaks = pretty_breaks(4),
            labels = function(x) sprintf('%.1g', x)) +
        scale_y_continuous('Precision', breaks = pretty_breaks(3),
            expand = expansion(mult = 0.07)) +
        scale_color_manual('Variant caller', values = tool_colors) +
        scale_fill_manual('Variant caller', values = tool_colors) +
        scale_shape_manual('Quality threshold', values = c('10' = 22, '20' = 21)) +
        scale_linetype_manual(NULL, values = c('EUR' = 'solid', 'EAS' = '22')) +

        guides(
            color = guide_legend(override.aes = list(shape = 21), order = 1),
            fill = guide_legend(order = 1),
            shape = guide_legend(override.aes = list(color = 'black'), order = 2),
            linetype = 'none',
        ) +
        theme_bw() +
        theme(
            axis.ticks = element_blank(),
            axis.title.x = element_text(size = 9),
            axis.title.y = element_text(size = 9),

            plot.title = element_text(face = 'bold', hjust = 0.5, size = 9,
                margin = margin()),
            plot.margin = margin(l = 3, r = 8),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = 'gray50'),
        )
}

plot_grid(get_legend(plots[[1]] +
        theme(
            legend.position = 'top',
            legend.margin = margin(b = -7),
            legend.box = 'vertical',
            legend.box.margin = margin(b = 15, l = 5),
            legend.spacing.y = unit(4, 'pt'),
            legend.title = element_text(size = 9),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.text = element_text(size = 9, margin = margin(l = -3)),
        )),
    plot_grid(
        plots[[1]] + theme(legend.position = 'none'),
        plots[[2]] + theme(legend.position = 'none'),
        labels = c('(a)', '(b)'), label_size = 9,
        label_x = c(-0.015, -0.035), label_y = 1.05
    ),
    ncol = 1,
    rel_heights = c(0.15, 1.0)
)
ggsave(file.path(plotdir, 'fig2.pdf'),
    width = 10, height = 5, dpi = 400, scale = 0.45,
    title = 'Figure 2')

metadata
filter(tab, sample == 'SIM-Fixed') |> as.data.frame()
filter(tab, sample == 'SIM-Poly' & database == 'EUR') |>
    select(!c(database, sample, thresh)) # |> as.data.frame()
filter(curves, tool == 'DeepVariant' & sample == 'SIM-Poly') |>
    arrange(score) |> head()
filter(tab, sample == 'SIM-Poly' & database == 'EAS') |>
    select(!c(database, sample))
