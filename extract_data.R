# Hypothesis creation.
library(readr)
library(nzilbb.labbcat)
library(here)

labbcat.url <- "https://labbcat.canterbury.ac.nz/demo"

ids <- getTranscriptIdsInCorpus(labbcat.url, "QB")

qb_words <- getMatches(
  labbcat.url,
  pattern = ".+",
  anchor.confidence.min = 50,
  transcript.expression = expressionFromIds(ids)
)

to_collect <- c(
  "word frequency", "syllables per minute", "participant_gender",
  "participant_age_category", "keyness", "orthography",
  "pos", "syllables", "syllable count"
)

qb_words <- qb_words |>
  appendLabels(
    layer.ids = to_collect
  ) |>
  # get prev word
  appendLabels(
    layer.ids = c(
      "orthography"
    ),
    target.offset = -1
  )

# save
write_rds(qb_words, 'data/raw_data.rds')

