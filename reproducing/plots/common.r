library(tidyverse)

load_curve <- function(filename, thresholds = seq(0, 200, 10)) {
    curve <- read_delim(gzfile(filename),
            '\t', comment = '#', col_names = F, show_col_types = F)
    col_names <- c('score', 'tp_base', 'fp', 'tp_call', 'fn', 'precision', 'recall', 'f1')
    if (nrow(curve) == 0 || ncol(curve) == length(col_names)) {
        colnames(curve) <- col_names
    } else {
        comment_lines <- readLines(gzfile(filename))
        stopifnot('#score\ttrue_positives_baseline\tfalse_positives\ttrue_positives_call' %in% comment_lines)
        colnames(curve) = col_names[1:4]
        curve <- cbind(curve, data.frame(fn = 0, precision = NA, recall = NA, f1 = NA))
    }
    curve$qual <- NA
    for (t in thresholds) {
        curve$qual[tail(which(curve$score >= t), n = 1)] <- t
    }
    curve
}

load_metadata <- function(dir, curves) {
    regions <- read.csv(file.path(dirname(dir), 'comparison.bed'),
        sep = '\t', comment = '#', header = F)
    curves1 <- tail(curves, n=1)
    data.frame(
        sum_len = sum(regions$V3 - regions$V2),
        n_regions = nrow(regions),
        baseline_vars = curves1$tp_base + curves1$fn,
        call_vars = curves1$tp_call + curves1$fp
        )
}

load_empty_metadata <- function(dir) {
    if (!file.exists(dir)) {
        return (data.frame(sum_len = 0, n_regions = 0,baseline_vars = 0, call_vars = 0))
    }
    comment_lines <- readLines(gzfile(file.path(dir, 'weighted_roc.tsv.gz')))
    var_counts <- strsplit(comment_lines[grep('variants', comment_lines)], ': ')
    stopifnot(var_counts[[1]][1] == '#total baseline variants')
    stopifnot(var_counts[[2]][1] == '#total call variants')
    
    regions <- suppressMessages(read_delim(file.path(dirname(dir), 'comparison.bed'), '\t', comment = '#', col_names = F))
    data.frame(
        sum_len = sum(regions$X3 - regions$X2),
        n_regions = nrow(regions),
        baseline_vars = as.numeric(var_counts[[1]][2]),
        call_vars = as.numeric(var_counts[[2]][2]))
}

.names_or_itself <- function(x) {
    if (is.null(names(x))) { x } else { names(x) }
}

.factorize_columns <- function(df, args) {
    for (i in seq_along(args)) {
        col_name <- names(args)[i]
        values <- args[[col_name]]
        if (is.null(names(values))) {
            df[[col_name]] <- factor(df[[col_name]], levels = values)
        } else {
            df[[col_name]] <- factor(values[df[[col_name]]], levels = values)
        }
    }
    df
}

load_curves <- function(fmt_path, ...) {
    args <- list(...)
    args_names <- lapply(args, .names_or_itself)
    grid <- do.call(expand_grid, args_names)
    
    n = nrow(grid)
    curves <- data.frame()
    metadata <- data.frame()
    for (i in 1:n) {
        grid_i <- grid[i,]
        path <- do.call(sprintf, c(fmt_path, as.list(grid_i)))
        if (file_test('-f', path)) {
            path_file <- path
            path_dir <- dirname(path)
        } else {
            path_file <- file.path(path, 'weighted_roc.tsv.gz')
            path_dir <- path
        }
        
        if (!file.exists(file.path(path_dir, 'done'))) {
            cat(sprintf('WARN: Problem with %s\n', path))
            metadata <- rbind(metadata,
                cbind(
                    data.frame(sum_len = 0, n_regions = 0,baseline_vars = 0, call_vars = 0),
                    grid_i))
            next
        }
        cat(sprintf('[%3d / %3d] Loading %s\n', i, n, path))

        curve <- suppressMessages(load_curve(path_file))
        if (nrow(curve) == 0) {
            metadata <- rbind(metadata, cbind(load_empty_metadata(path_dir), grid_i))
            next
        }
        metadata <- rbind(metadata, load_metadata(path_dir, curve) |> cbind(grid_i))
        curves <- rbind(curves, cbind(curve, grid_i))
    }
    list(
        curves = .factorize_columns(curves, args),
        metadata = .factorize_columns(metadata, args))
}

fmt_len <- function(len, thin_space = T) {
    space <- ifelse(thin_space, ' ', ' ')
    case_when(
        len < 1000 ~ sprintf('%d%sbp',   len, space),
        len < 1e4  ~ sprintf('%.2f%skb', len / 1e3, space),
        len < 1e5  ~ sprintf('%.1f%skb', len / 1e3, space),
        len < 1e6  ~ sprintf('%.0f%skb', len / 1e3, space),
        len < 1e7  ~ sprintf('%.2f%sMb', len / 1e6, space),
        len < 1e8  ~ sprintf('%.1f%sMb', len / 1e6, space),
        T          ~ sprintf('%.0f%sMb', len / 1e6, space),
    )
}

tableau10 <- ggthemes::tableau_color_pal()(10)
hue_circle <- ggthemes::tableau_color_pal('Hue Circle')(19)
tool_fills <- c(
    'GATK'        = '#5687cc', # hue_circle[19] made a bit brighter
    'FreeBayes'   = hue_circle[13],
    'DeepVariant' = hue_circle[5],
    'ParascopyVC' = hue_circle[9]
)

tool_colors <- c(
    'GATK'        = '#5687cc', # hue_circle[19] made a bit brighter
    'FreeBayes'   = tableau10[3],
    'DeepVariant' = hue_circle[5],
    'ParascopyVC' = '#e6a817'
)

add_f_metrics <- function(df, betas = c(1/4, 1/3, 1/2, 1, 2, 3, 4)) {
    for (beta in betas) {
        beta2 <- beta * beta
        df[[sprintf('F_%.3g', beta)]] <- (1 + beta2) *
            with(df, precision * recall / (beta2 * precision + recall))
    }
    df
}
