library(dplyr)
names <- readxl::read_xlsx("extdata/names.xlsx", sheet = 1)
names$id <- paste0("ui_", names$id)

# ad translation for hc
names <- names %>%
  dplyr::filter(!id %in% c("ui_3hc", "ui_3hc2"))

names <- dplyr::bind_rows(
  names,
  dplyr::tibble(
    id = "ui_3conc",
    name = "Concentration"
  ),
  dplyr::tibble(
    id = "ui_thresh_type",
    name = "Get estimate by"
  ),
  dplyr::tibble(
    id = "ui_2dlrds",
    name = "Plot .rds"
  ),
  dplyr::tibble(
    id = "ui_checkHc",
    name = "Plot Threshold/Concentration"
  ),
  dplyr::tibble(
    id = "ui_3hc",
    name = "The model averaged estimate of the concentration that affects {percent} % of species is {conc}"
  ),
  dplyr::tibble(
    id = "ui_3hc2",
    name = "The model averaged estimate of the fraction affected by a concentration of {conc} is {percent} % of species"
  ),
  dplyr::tibble(
    id = "ui_1table1",
    name = "Table"
  ),
  dplyr::tibble(
    id = "ui_xmax",
    name = "X-axis maximum"
  ),
  dplyr::tibble(
    id = "ui_xmin",
    name = "X-axis minimum"
  ),
  dplyr::tibble(
    id = "ui_xlog",
    name = "Log x-axis"
  ),
  dplyr::tibble(
    id = "ui_sizeLabel",
    name = "Label size"
  ),
  dplyr::tibble(
    id = "ui_size",
    name = "Text size"
  ),
  dplyr::tibble(
    id = "ui_adjustLabel",
    name = "Shift label"
  ),
  dplyr::tibble(
    id = "ui_xbreaks",
    name = "X-axis ticks"
  )
)

chk::check_key(names, "id")

pal <- RColorBrewer::brewer.pal.info
pals <- pal[which(pal$category == "qual"), ] %>% row.names()