library(dplyr)
names <- readxl::read_xlsx("extdata/names.xlsx", sheet = 1)
names$id <- paste0("ui_", names$id)

names <- dplyr::bind_rows(
  names,
  dplyr::tibble(
    id = "ui_2dlrds",
    name = "Plot .rds"
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

salinity <- BNPdensity::salinity
acidity <- tibble(ind=1:length(BNPdensity::acidity), Conc=BNPdensity::acidity)

pal <- RColorBrewer::brewer.pal.info
pals <- pal[which(pal$category == "qual"), ] %>% row.names()
