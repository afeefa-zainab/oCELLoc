## code to prepare `DATASET` dataset goes here

# In a script inside data-raw/
human_ref <- read.delim("inst/references/human_ref.txt", sep="\t", stringsAsFactors = FALSE)
mouse_ref <- read.delim("inst/references/mouse_ref.txt", sep="\t", stringsAsFactors = FALSE)

# Save to the package's data/ directory
usethis::use_data(human_ref, overwrite = TRUE)
usethis::use_data(mouse_ref, overwrite = TRUE)

usethis::use_data(DATASET, overwrite = TRUE)