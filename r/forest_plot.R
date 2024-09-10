#
# Purpose: forest_plot and model_fit functions
# These functions automate plotting of multi-state and Cox model hazards,
# as well as model fit statistics
#

library(broom)
library(dplyr)
library(forcats)
library(ggnewscale)
library(ggpp)
library(gt)
library(lubridate)
library(patchwork)
library(tidyr)

# forest_plot function
forest_plot <- function(fit,
                        cox_fit = NULL,
                        covars = NULL,
                        transition1 = "State 1 - State 2",
                        transition2 = NULL,
                        state_labels = c("COVID symptoms", "Non-COVID symptoms or asymptomatic"),
                        table = FALSE) {

    # extract the covs and covnames
    covnames <- fit$data$mf |> colnames()
    covs <- tibble(
        name = fit$data$mf[, 1] |> levels(),
        group = covnames[1],
        value_id = transition1,
        n = as.numeric(fit$data$mf[, 1] |> table())
    )
    i <- 2
    while (covnames[i] != "(subject)") {
        covs <- covs |> bind_rows(tibble(
            name = fit$data$mf[, i] |> levels(),
            group = covnames[i],
            value_id = transition1,
            n = as.numeric(fit$data$mf[, i] |> table())
        ))
        i <- i + 1
    }

    if (!is.null(transition2)) {
        covs <- covs |>
            bind_rows(
                covs |> mutate(value_id = transition2)
            )
    }

    covtable <- covs |>
        full_join(
            hazard.msm(fit) |> enframe() |> unnest_longer(col = value) |>
                mutate(
                    group = gsub(
                        paste0(
                            ".*(",
                            paste(covs$group |> unique(), collapse = "|"),
                            ").*"
                        ),
                        "\\1", name
                    ),
                    int_group = gsub(
                        paste0(
                            ".*(",
                            paste(covs$group |> unique(), collapse = "|"),
                            ").*(",
                            paste(covs$group |> unique(), collapse = "|"),
                            ").*"
                        ),
                        "\\1:\\2", name
                    ),
                    group = if_else(int_group != name, int_group, group),
                    name = gsub(paste(covs$group |> unique(), collapse = "|"), "", name),
                ),
            by = c("name", "group", "value_id")
        ) |>
        mutate(
            name = factor(name, levels = unique(name), labels = gsub("_", " ", unique(name))),
            group = str_to_sentence(gsub("_", " ", group)),
            est = value[, 1],
            lower = value[, 2],
            upper = value[, 3],
            est = if_else(is.na(est), 1, est)
        ) |>
        select(-value)

    if (!is.null(cox_fit)) {
        covtable <- covtable |>
            left_join(
                cox_fit |>
                    tidy(exponentiate = TRUE, conf.int = TRUE) |>
                    mutate(
                        name = term,
                        group = gsub(
                            paste0(
                                ".*(",
                                paste(covs$group |> unique(), collapse = "|"),
                                ").*"
                            ),
                            "\\1", name
                        ),
                        int_group = gsub(
                            paste0(
                                ".*(",
                                paste(covs$group |> unique(), collapse = "|"),
                                ").*(",
                                paste(covs$group |> unique(), collapse = "|"),
                                ").*"
                            ),
                            "\\1:\\2", name
                        ),
                        group = if_else(int_group != name, int_group, group),
                        name = gsub(paste(covs$group |> unique(), collapse = "|"), "", name)
                    ) |>
                    mutate(
                        name = factor(name, levels = unique(name), labels = gsub("_", " ", unique(name))),
                        group = str_to_sentence(gsub("_", " ", group)),
                        c_est = estimate,
                        c_lower = conf.low,
                        c_upper = conf.high
                    ) |>
                    select(name, c_est, c_lower, c_upper),
                by = c("name")
            )
    } else {
        covtable <- covtable |>
            cbind(tibble(c_est = NA, c_lower = NA, c_upper = NA))
    }

    if (!is.null(covars)) {
        covtable <- covtable |> filter(group %in% str_to_sentence(gsub("_", " ", covars)))
    }

    if (!is.null(cox_fit)) {
        p1 <- covtable |>
            filter(value_id == transition1) |>
            mutate(n = if_else(is.na(n), 1, n), group = if_else(group == "Vaccine short", "Vaccine", group)) |>
            ggplot() +
            aes(y = name, x = est, xmin = lower, xmax = upper) +
            labs(x = "Hazard ratio (95% CI, log scale)", y = "") +
            scale_x_continuous(trans = "log10") +
            scale_y_discrete(limits = rev) +
            geom_tile(aes(fill = factor(group, levels = group |> unique()), width = Inf), alpha = 0.3) +
            geom_point(shape = 15, aes(color = "Multi-state model", size = n, alpha = 0.2), position = position_nudge(y = 0.1)) +
            geom_point(aes(x = c_est, color = "Cox model"), shape = 17, position = position_nudge(y = -0.1)) +
            geom_errorbarh(height = 0.2, aes(color = "Multi-state model"), position = position_nudge(y = 0.1)) +
            geom_errorbarh(aes(xmin = c_lower, xmax = c_upper, color = "Cox model"), height = 0.2, position = position_nudge(y = -0.1)) +
            geom_vline(xintercept = 1, linetype = "longdash", colour = "black") +
            scale_fill_manual(values = rep(c("white", "grey"), 9)) +
            scale_color_manual(
                name = "",
                values = c(`Cox model` = "#d24c4e", `Multi-state model` = "#008ecd")
            ) +
            theme_classic(14) +
            theme(
                strip.background = element_blank(),
                strip.placement = "outside",
                panel.spacing = unit(0, "lines"),
                axis.title.x = element_text(),
                axis.title.y = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "top",
                panel.grid.major.x = element_line(colour = "lightgrey", linewidth = 0.5)
            ) +
            guides(size = "none", alpha = "none", fill = "none") +
            facet_grid(rows = vars(factor(group, levels = group |> unique())), scales = "free_y", space = "free_y", switch = "y")
    } else {
        p1 <- covtable |>
            filter(!(est == 1 & !is.na(lower))) |>
            mutate(n = if_else(is.na(n), 1, n), group = if_else(group == "Vaccine short", "Vaccine", group)) |>
            # filter(value_id == transition1) |>
            ggplot() +
            aes(
                y = name, x = est, xmin = lower, xmax = upper,
                label = if_else(!is.na(lower),
                    paste0(sprintf("%.2f", round(est, 2)), " (", sprintf("%.2f", round(lower, 2)), "-", sprintf("%.2f", round(upper, 2)), ")"),
                    "1 (baseline)"
                )
            ) +
            labs(x = "Hazard ratio (95% CI, log scale)", y = "") +
            scale_x_continuous(trans = "log10") +
            scale_y_discrete(limits = rev) +
            geom_tile(aes(fill = factor(group, levels = group |> unique()), width = Inf), alpha = 0.3) +
            scale_fill_manual(values = rep(c("white", "grey"), 9)) +
            new_scale_fill() +
            geom_point(shape = 22, position = position_dodge(width = 0.6), aes(size = n, alpha = 0.2, color = value_id, fill = value_id)) +
            geom_errorbarh(height = 0.2, position = position_dodge(width = 0.6), aes(color = value_id)) +
            geom_vline(xintercept = 1, linetype = "longdash", colour = "black") +
            coord_cartesian(clip = "off") +
            theme_classic(14) +
            theme(
                strip.background = element_blank(),
                strip.placement = "outside",
                panel.spacing = unit(0, "lines"),
                axis.title.x = element_text(),
                axis.title.y = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "none",
                panel.grid.major.x = element_line(colour = "lightgrey", linewidth = 0.5)
            ) +
            facet_grid(rows = vars(factor(group, levels = group |> unique())), scales = "free_y", space = "free_y", switch = "y")
    }
    if (!is.null(transition2)) {
        p1 <- p1 +
            scale_color_discrete(label = state_labels) +
            theme(legend.position = "bottom", legend.title = element_blank()) +
            guides(shape = "none", fill = "none", size = "none", alpha = "none", fill_new = "none")
    }
    if (table == TRUE) {
        if (!is.null(transition2)) {
            p1 <- covtable |>
                filter(!(est == 1 & !is.na(lower))) |>
                select(value_id, group, name, est, lower, upper)
        } else {
            p1 <- covtable |>
                filter(!(est == 1 & !is.na(lower))) |>
                select(group, name, est, lower, upper)
        }
    }
    p1
}

# model fit
model_fit <- function(model, prev, dataset = NA, n_states = 2, mintime = 0, timezero = 0, maxtime, ...) {
    # model assessment
    # cannot be run with the full dataset, but can make this quicker by down-sampling the number of subjects to use
    # to do this need to specify BOTH the subset and the initstates options (perhaps something to raise with Chris)
    # subset is a vector of study_ids, initstates is the number of individuals occupying each state at time zero

    out <- list()

    # to plot
    sample <- dataset |>
        slice_sample(prop = 0.2) |>
        distinct(study_id) # 10% sample to compare against(~100,000 obs)
    states <- dataset |>
        arrange(time) |>
        distinct(study_id, .keep_all = TRUE) |>
        count(state) |>
        select(n)

    # to get prev data
    prev <- prevalence.msm(model,
        timezero,
        times = seq(timezero, maxtime, by = 1),
        mintime, subset = sample$study_id, initstates = states$n,
        covariates = "mean",
        ci = "normal",
        ...
    )

    if (n_states == 2) {
        gg_tibble <- tibble(
            x = prev$`Observed percentages`[, 1] |> labels() |> as.numeric(),
            State1 = prev$`Observed percentages`[, 1],
            State2 = prev$`Observed percentages`[, 2],
            group = "Observed"
        ) |>
            bind_rows(
                tibble(
                    x = prev$`Expected percentages`$estimates[, 1] |> labels() |> as.numeric(),
                    State1 = prev$`Expected percentages`$estimates[, 1],
                    State2 = prev$`Expected percentages`$estimates[, 2],
                    group = "Expected"
                )
            ) |>
            pivot_longer(-c(x, group))


        gg_tibble_ci <- tibble(
            x = prev$`Expected percentages`$estimates[, 1] |> labels() |> as.numeric(),
            value = prev$`Expected percentages`$estimates[, 1],
            lcl = prev$`Expected percentages`$ci[, 1, "2.5%"],
            ucl = prev$`Expected percentages`$ci[, 1, "97.5%"],
            name = "State1",
            group = "Expected"
        ) |>
            bind_rows(
                tibble(
                    x = prev$`Expected percentages`$estimates[, 1] |> labels() |> as.numeric(),
                    value = prev$`Expected percentages`$estimates[, 2],
                    lcl = prev$`Expected percentages`$ci[, 2, "2.5%"],
                    ucl = prev$`Expected percentages`$ci[, 2, "97.5%"],
                    name = "State2",
                    group = "Expected"
                )
            )
    } else if (n_states == 3) {
        gg_tibble <- tibble(
            x = prev$`Observed percentages`[, 1] |> labels() |> as.numeric(),
            State1 = prev$`Observed percentages`[, 1],
            State2 = prev$`Observed percentages`[, 2],
            State3 = prev$`Observed percentages`[, 3],
            group = "Observed"
        ) |>
            bind_rows(
                tibble(
                    x = prev$`Expected percentages`$estimates[, 1] |> labels() |> as.numeric(),
                    State1 = prev$`Expected percentages`$estimates[, 1],
                    State2 = prev$`Expected percentages`$estimates[, 2],
                    State3 = prev$`Expected percentages`$estimates[, 3],
                    group = "Expected"
                )
            ) |>
            pivot_longer(-c(x, group))


        gg_tibble_ci <- tibble(
            x = prev$`Expected percentages`$estimates[, 1] |> labels() |> as.numeric(),
            value = prev$`Expected percentages`$estimates[, 1],
            lcl = prev$`Expected percentages`$ci[, 1, "2.5%"],
            ucl = prev$`Expected percentages`$ci[, 1, "97.5%"],
            name = "State1",
            group = "Expected"
        ) |>
            bind_rows(
                tibble(
                    x = prev$`Expected percentages`$estimates[, 1] |> labels() |> as.numeric(),
                    value = prev$`Expected percentages`$estimates[, 2],
                    lcl = prev$`Expected percentages`$ci[, 2, "2.5%"],
                    ucl = prev$`Expected percentages`$ci[, 2, "97.5%"],
                    name = "State2",
                    group = "Expected"
                )
            ) |>
            bind_rows(
                tibble(
                    x = prev$`Expected percentages`$estimates[, 1] |> labels() |> as.numeric(),
                    value = prev$`Expected percentages`$estimates[, 3],
                    lcl = prev$`Expected percentages`$ci[, 3, "2.5%"],
                    ucl = prev$`Expected percentages`$ci[, 3, "97.5%"],
                    name = "State3",
                    group = "Expected"
                )
            )
    }

    p1 <- gg_tibble |>
        filter(name == "State1") |>
        ggplot() +
        aes(x, value, color = group) +
        geom_line() +
        geom_ribbon(
            data = gg_tibble_ci |> filter(name == "State1"),
            aes(x = x, ymin = lcl, ymax = ucl, fill = group),
            alpha = 0.2,
            linetype = 0
        ) +
        scale_x_continuous(limits = c(0, maxtime), expand = expansion(c(0,0))) +
        scale_y_continuous(limits = c(NA, 100)) +
        labs(color = "", x = "Week of study", y = "Prevalence") +
        theme_bw() +
        theme(plot.subtitle = element_text(hjust = 0.5)) +
        guides(fill = "none")

    p2 <- gg_tibble |>
        filter(name == "State2") |>
        ggplot() +
        aes(x, value, color = group) +
        geom_line() +
        geom_ribbon(
            data = gg_tibble_ci |> filter(name == "State2"),
            aes(x = x, ymin = lcl, ymax = ucl, fill = group),
            alpha = 0.2,
            linetype = 0
        ) +
        scale_x_continuous(limits = c(0, maxtime), expand = expansion(c(0,0))) +
        scale_y_continuous(limits = c(0, NA)) +
        labs(color = "", x = "Week of study", y = "Prevalence") +
        theme_bw() +
        theme(plot.subtitle = element_text(hjust = 0.5)) +
        guides(fill = "none")

    # patchwork plot with shared legend at the bottom
    if (n_states == 2) {
        out$fig <- p1 + p2 + plot_layout(guides = "collect") + theme(legend.position = "top")
    } else if (n_states == 3) {
        p3 <- gg_tibble |>
            filter(name == "State3") |>
            ggplot() +
            aes(x, value, color = group) +
            geom_line() +
            geom_ribbon(
                data = gg_tibble_ci |> filter(name == "State2"),
                aes(x = x, ymin = lcl, ymax = ucl, fill = group),
                alpha = 0.2,
                linetype = 0
            ) +
            scale_x_continuous(limits = c(0, maxtime), expand = expansion(c(0,0))) +
            scale_y_continuous(limits = c(0, NA)) +
            labs(color = "", x = "Week of study", y = "Prevalence") +
            theme_bw() +
            theme(plot.subtitle = element_text(hjust = 0.5)) +
            guides(fill = "none")

        out$fig <- p1 + p2 + p3 + plot_layout(guides = "collect") + theme(legend.position = "top")
    }

    # just the second plot
    out$fig1 <- p2 + theme(legend.position = "top")

    return(out)
}
