
ReturnDifferences <- . %>%
  #mutate(freq = scale(freq, center = FALSE)) %>%
  group_by(status) %>%
  summarise(freq_mean = mean(freq)) %>%
  ungroup() %>%
  summarise(value = max(freq_mean) - min(freq_mean),
            value2 = status[which.max(freq_mean)]) %>%
  mutate(value3 = if_else(value2 == 1, value, -1 * value)) %>%
  pull(value3)
