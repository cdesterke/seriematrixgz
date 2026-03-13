library(tidyverse)
library(janitor)

# -----------------------------
# 1) Lire directement le fichier .gz
# -----------------------------
file <- "GSE4412-GPL97_series_matrix.txt.gz"
raw  <- readLines(gzfile(file))

# -----------------------------
# 2) Extraire toutes les lignes !Sample_
# -----------------------------
sample_lines <- raw[grepl("^!Sample_", raw)]

parse_geo_line <- function(line) {
  parts  <- strsplit(line, "\t")[[1]]
  field  <- sub("^!Sample_", "", parts[1])
  values <- gsub('"', '', parts[-1])
  tibble(field = field, values = list(values))
}

parsed <- map_dfr(sample_lines, parse_geo_line)

# Nombre d'échantillons
n_samples <- length(parsed$values[[1]])

# -----------------------------
# 3) Séparer caractéristiques et autres champs
# -----------------------------
has_char <- "characteristics_ch1" %in% parsed$field

if (has_char) {
  char_rows <- parsed %>% filter(field == "characteristics_ch1")
  meta_rows <- parsed %>% filter(field != "characteristics_ch1")
} else {
  char_rows <- NULL
  meta_rows <- parsed
}

# -----------------------------
# 4) Construire les métadonnées simples
# -----------------------------
meta_rows_unique <- meta_rows %>%
  group_by(field) %>%
  slice(1) %>%
  ungroup()

meta_list <- lapply(seq_len(nrow(meta_rows_unique)), function(i) {
  v <- meta_rows_unique$values[[i]]
  if (length(v) == n_samples) v else rep(NA_character_, n_samples)
})

names(meta_list) <- meta_rows_unique$field
meta_df <- as_tibble(meta_list)

# -----------------------------
# 5) Construire les caractéristiques (si présentes)
# -----------------------------
if (has_char) {

  # Matrice N × M
  char_matrix <- do.call(rbind, char_rows$values)
  char_df     <- as_tibble(t(char_matrix))  # M × N

  # Fonction pour extraire le préfixe visible
  extract_prefix <- function(x) {
    prefix <- sub(":.*$", "", x[which(x != "")][1])
    make_clean_names(prefix)
  }

  # Fonction pour extraire la valeur
  extract_value <- function(x) {
    sub("^.*?: ", "", x)
  }

  # Construire les colonnes automatiquement
  char_clean <- map_dfc(seq_len(ncol(char_df)), function(i) {
    tibble(val = extract_value(char_df[[i]]))
  })

  # Noms de colonnes basés sur les préfixes visibles
  colnames(char_clean) <- map_chr(seq_len(ncol(char_df)), ~ extract_prefix(char_df[[.x]]))

  # Fusion
  pheno_full <- bind_cols(meta_df, char_clean)

} else {

  pheno_full <- meta_df
}

# -----------------------------
# 6) Résultat final
# -----------------------------
pheno <- as.data.frame(pheno_full)

pheno



write.csv(pheno,file="pheno.csv",row.names=F)

