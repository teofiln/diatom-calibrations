library(jsonlite)
library(CladeDate)
library(purrr)
library(dplyr)
library(rrapply)

# Method used to get calibration densities
method <- "Beta"

# load the nested-calibrations.json file
nested <- read_json("nested-calibrations.json")

# find all the clades. we'll recurse over them to find calibration points
# nested within

multiples <- rrapply::rrapply(
  nested,
  condition = \(x, .xsiblings, .xname) {
    "multiple" %in% .xsiblings && "name" %in% .xname
  },
  classes = "ANY",
  how = "melt"
)

# recurse over clades to find the calibrations nested within
get_calibs_for_clade <- function(nested, clade) {
  lower <- rrapply::rrapply(
    nested,
    condition = \(x, .xparents) clade %in% .xparents,
    classes = "ANY",
    how = "prune"
  ) |>
    rrapply::rrapply(
      condition = \(x, .xname) .xname == "lower",
      classes = "ANY",
      how = "melt"
    )

  upper <- rrapply::rrapply(
    nested,
    condition = \(x, .xparents) clade %in% .xparents,
    classes = "ANY",
    how = "prune"
  ) |>
    rrapply::rrapply(
      condition = \(x, .xname) .xname == "upper",
      classes = "ANY",
      how = "melt"
    )

  name <- rrapply::rrapply(
    nested,
    condition = \(x, .xparents) clade %in% .xparents,
    classes = "ANY",
    how = "prune"
  ) |>
    rrapply::rrapply(
      condition = \(x, .xname) .xname == "name",
      classes = "ANY",
      how = "melt"
    ) |>
    dplyr::filter(!grepl("clade", value))

  bracket1 <- rrapply::rrapply(
    nested,
    condition = \(x, .xparents) clade %in% .xparents,
    classes = "ANY",
    how = "prune"
  ) |>
    rrapply::rrapply(
      condition = \(x, .xname) .xname == "1",
      classes = "ANY",
      how = "melt"
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(value)

  bracket2 <- rrapply::rrapply(
    nested,
    condition = \(x, .xparents) clade %in% .xparents,
    classes = "ANY",
    how = "prune"
  ) |>
    rrapply::rrapply(
      condition = \(x, .xname) .xname == "2",
      classes = "ANY",
      how = "melt"
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(value)

  dplyr::bind_cols(
    calibration_name = name$value,
    upper = upper$value,
    lower = lower$value
  ) |>
    dplyr::mutate(clade = clade) |>
    dplyr::mutate(bracket_1 = bracket1) |>
    dplyr::mutate(bracket_2 = bracket2)
}

all_clades_list <- purrr::map(
  multiples$value,
  get_calibs_for_clade,
  nested = nested
) |>
  setNames(multiples$value)

all_clades_df <- dplyr::bind_rows(all_clades_list)

readr::write_csv(
  all_clades_df, file.path("output", "dates-by-calibration-point.csv")
)
jsonlite::write_json(
  all_clades_list, file.path("output", "dates-by-calibration-point.json")
)

pdf_name <- sprintf(
  "nested-calibrations-%s-distinct-dates-only.pdf",
  method
)
pdf(file.path("output", pdf_name))
all_clades_bounds <- purrr::map(
  all_clades_list, function(x) {
    # if (x$clade[[1]] == "diatom_clade") {
    #   return(NULL)
    # }
    message(x$clade[[1]])
    ages <- x[, 2:3] |>
      dplyr::distinct() |>
      dplyr::arrange(dplyr::desc(lower))
    str(ages)
    CD <- CladeDate::clade.date(
      ages = ages,
      method = method,
      plot = TRUE,
      KStest = TRUE
    )

    title(
      main = paste(
        x$clade[[1]], "\n",
        paste0(round(CD$Quantiles, 2), collapse = ", ")
      ),
      cex.main = 1,
      adj = 1,
      line = -5
    )
    # print(CD$Quantiles)
    c(
      Clade = x$clade[[1]],
      CD$Quantiles,
      method = method,
      fit_model = CD$PDFfit.model
    )
  }
)
dev.off()

all_clades_bounds_df <- dplyr::bind_rows(all_clades_bounds) |>
  dplyr::mutate_at(2:4, as.numeric) |>
  dplyr::mutate_if(is.numeric, round, 2)

clade_brackets <- all_clades_df |>
  dplyr::group_by(clade) |>
  dplyr::slice_head(n = 1) |>
  dplyr::select(clade, bracket_1, bracket_2) |>
  dplyr::distinct()

all_clades_bounds_df <- dplyr::left_join(
  all_clades_bounds_df,
  clade_brackets,
  by = c("Clade" = "clade")
)

csv_name <- sprintf(
  "nested-calibrations-%s-distinct-dates-only.csv",
  method
)
write.csv(
  all_clades_bounds_df,
  file.path("output", csv_name)
)
