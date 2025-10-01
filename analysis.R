library(tidyverse)
library(here)
library(broom.mixed)
library(modelbased)
library(lme4)
library(lmerTest)
library(patchwork)

theme_set(theme_bw())

qb_words <- read_rds('data/raw_data.rds')

# renaming variables
qb_words <- qb_words |>
  rename(
    transcript = Transcript,
    participant = Participant,
    corpus = Corpus,
    line = Line,
    line_end = LineEnd,
    match_id = MatchId,
    gender = participant_gender,
    age = participant_age_category,
    speech_rate = syllables.per.minute,
    frequency = word.frequency,
    start = Target.word.start,
    end = Target.word.end,
    prev_word = Token.minus.1.orthography,
    syll_count = syllable.count
  )

# arrange by participant, line, target onset.
qb_words <- qb_words |>
  arrange(participant, line, start)

# do corpus counts and prep for joining (want count of 'y', i.e. prev word)
prev_counts <- qb_words |>
  count(orthography) |>
  rename(
    prev_word = orthography,
    prev_count = n
  )

corpus_counts <- qb_words |>
  count(orthography) |>
  rename(
    word_count = n
  )

# do paired counts
context_counts <- qb_words |>
  count(orthography, prev_word) |>
  rename(context_count = n)

# count word within utterance
qb_words <- qb_words |>
  group_by(participant, line) |>
  mutate(
    a_word = 1,
    word_no = cumsum(a_word),
    final = word_no == max(word_no)
  ) |>
  select(-a_word) |>
  ungroup()

qb_words <- qb_words |>
  left_join(corpus_counts) |>
  left_join(context_counts) |>
  left_join(prev_counts)

# create new variables.
qb_words <- qb_words |>
  mutate(
    word_duration = end - start,
    prev_pred = context_count / prev_count,
    back_pred = context_count / word_count
  )

# informativity
informativity <- qb_words |>
  group_by(orthography, prev_word) |>
  summarise(
    prev_pred = first(prev_pred, na_rm = TRUE),
    back_pred = first(back_pred, na_rm = TRUE)
  ) |>
  group_by(orthography) |>
  summarise(
    prev_info = - sum(back_pred * log(prev_pred), na.rm = TRUE)
  )

qb_words <- qb_words |>
  left_join(informativity)

# get repetitions. Possibly there's a less computationally intensive way to
# do this.
qb_words <- qb_words %>%
  mutate(
    rep_30 = map2_lgl(
      start,
      orthography,
      ~ .y %in% (qb_words |>
        filter(
          between(start, max(c(0, .x - 30)), .x-0.01)
        ) |>
        pull(orthography))
    )
  )

# Create word duration column
qb_words <- qb_words |>
  mutate(
    word_duration = end - start
  )

## Filtering from pre-reg
# (1) We will remove unfinished words (marked with '~' in the corpus).
# (2) We will remove function words (i.e. we include nouns, adjectives, adverbs and verbs)
# (3) We remove tokens ending in "'s".
# (4) We will exclude tokens of words whose duration are more than 2.5SD from the mean duration given their number of syllables.

qb_words_filtered <- qb_words |>
  filter(
    !str_detect(orthography, "~"),
    str_detect(pos, "(^JJ)|(^NN)|(^RB)|(^V)"),
    !str_detect(orthography, "'s"),
    between(
      word_duration,
      mean(word_duration) - 2.5 * sd(word_duration),
      mean(word_duration) + 2.5 * sd(word_duration)
    )
  )

# 23737 -> 10850

# 1. Plot

hist(qb_words$word_duration)
hist(qb_words_filtered$word_duration)

qb_words_filtered |>
  ggplot(
    aes(
      x = final,
      y = word_duration
    )
  ) +
  geom_boxplot()

qb_words_filtered |>
  ggplot(
    aes(
      x = prev_info,
      y = word_duration
    )
  ) +
  geom_point(alpha = 0.5)

qb_words_filtered |>
  ggplot(
    aes(
      x = word_count,
      y = word_duration
    )
  ) +
  geom_point(alpha = 0.5)

# 2. Fit

qb_words_filtered <- qb_words_filtered |>
  mutate(
    # Scale quant vars
    freq_s = scale(word_count)[,1],
    info_s = scale(prev_info)[,1],
    speech_rate_s = scale(speech_rate)[,1],
    syll_count_s = scale(syll_count)[,1],
    # make logical vectors factors
    final = factor(final),
    rep_30 = factor(rep_30)
  )

prereg_fit <- lmer(
  word_duration ~ freq_s + info_s + final + pos + speech_rate_s + rep_30 +
    syll_count_s +
    (1 + freq_s + info_s + final|participant) + (1|orthography),
  data = qb_words_filtered
)

# Uh oh, doesn't converge.

# From prereg - we will first increase the iterations.

prereg_fit_2 <- lmer(
  word_duration ~ freq_s + info_s + final + pos + speech_rate_s + rep_30 +
    syll_count_s +
    (1 + freq_s + info_s + final|participant) + (1|orthography),
  data = qb_words_filtered,
  control=lmerControl(optCtrl=list(maxeval=20000))
)

# Still no! Simplify random effects

prereg_fit_3 <- lmer(
  word_duration ~ freq_s + info_s + final + pos + speech_rate_s + rep_30 +
    syll_count_s +
    (1 + freq_s + info_s + final||participant) + (1|orthography),
  data = qb_words_filtered,
  control=lmerControl(optCtrl=list(maxeval=1000))
)

# Still no! Simplify random effects

prereg_fit_4 <- lmer(
  word_duration ~ freq_s + info_s + final + pos + speech_rate_s + rep_30 +
    syll_count_s +
    (1 + info_s + final||participant) + (1|orthography),
  data = qb_words_filtered,
  control=lmerControl(optCtrl=list(maxeval=1000))
)

# Still no! Simplify random effects
prereg_fit_5 <- lmer(
  word_duration ~ freq_s + info_s + final + pos + speech_rate_s + rep_30 +
    syll_count_s +
    (1 + final||participant) + (1|orthography),
  data = qb_words_filtered,
  control=lmerControl(optCtrl=list(maxeval=1000))
)


# Still no! Simplify random effects
prereg_fit_6 <- lmer(
  word_duration ~ freq_s + info_s + final + pos + speech_rate_s + rep_30 +
    syll_count_s +
    (1|participant) + (1|orthography),
  data = qb_words_filtered,
  control=lmerControl(optCtrl=list(maxeval=1000))
)

# save model
write_rds(prereg_fit_6, here('models', 'prereg_fit.rds'))


# 3. Check

## Some residual checks
qqnorm(resid(prereg_fit_6))
qqline(resid(prereg_fit_6))

## Looks bad at high end, but not unexpected given distribution of word
## durations. We have much more room for error at the upper end of the
## distribution.

plot(prereg_fit_6)

# 4. Interpret

# Look at summary
summary(prereg_fit_6)

# Look at tidy table of fixed effects
tidy(prereg_fit_6, effects = 'fixed')

# Look at fixed effects with associated hypotheses.
tidy(prereg_fit_6) |>
  filter(
    term %in% c("freq_s", "info_s", "finalTRUE")
  ) |>
  select(term, estimate, statistic, std.error, p.value)

# Plot

final_predictions <- predict_response(
  prereg_fit_6,
  terms = "final",
  margin = "empirical"
)

final_plot <- plot(final_predictions) + labs(title=NULL)
final_plot

freq_predictions <- predict_response(
  prereg_fit_6,
  terms = "freq_s",
  margin = "empirical"
)

freq_plot <- plot(freq_predictions) + labs(title=NULL)
freq_plot

info_predictions <- predict_response(
  prereg_fit_6,
  terms = "info_s",
  margin = "empirical"
)

info_plot <- plot(info_predictions) + labs(title=NULL)
info_plot

combined_plot <- freq_plot + info_plot + final_plot +
  plot_annotation(tag_levels = "A")
combined_plot

ggsave(
  here('plots', 'prediction_plot.png'),
  combined_plot,
  units = "mm",
  width = 200,
  height = 100
)

# Bonus: What if we want to go back to original scale? Just reverse the scaling applied
# to the explanatory variable. Scaling replaces a variable with
# scaled_x = (x - mean(x)) / sd(x), so reverse with x = scaled_x * sd(x) + mean(x)
# e.g.:
info_predictions <- info_predictions |>
  mutate(
    x = x * sd(qb_words_filtered$prev_info) + mean(qb_words_filtered$prev_info)
  )
plot(info_predictions)
