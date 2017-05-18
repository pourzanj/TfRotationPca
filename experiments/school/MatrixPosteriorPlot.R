library(tidyr)

Ws0G <- Ws0 %>%
  gather(Element, Sample) %>%
  mutate_at(vars(Element), as.factor) %>%
  group_by(Element) %>%
  summarize(Q2 = quantile(Sample, 0.025),
            Q50 = quantile(Sample, 0.5),
            Q97 = quantile(Sample, 0.975)) %>%
  mutate(Column = "1")

Ws1G <- Ws1 %>%
  gather(Element, Sample) %>%
  mutate_at(vars(Element), as.factor) %>%
  group_by(Element) %>%
  summarize(Q2 = quantile(Sample, 0.025),
            Q50 = quantile(Sample, 0.5),
            Q97 = quantile(Sample, 0.975)) %>%
  mutate(Column = "2")

Ws2G <- Ws2 %>%
  gather(Element, Sample) %>%
  mutate_at(vars(Element), as.factor) %>%
  group_by(Element) %>%
  summarize(Q2 = quantile(Sample, 0.025),
            Q50 = quantile(Sample, 0.5),
            Q97 = quantile(Sample, 0.975)) %>%
  mutate(Column = "3")

Ws <- rbind(Ws0G, Ws1G) %>%
  mutate(Element = str_extract(Element,"[0-9]+")) %>%
  #mutate_at(vars(Column), as.factor) %>%
  mutate(Element = factor(Element, levels = c("12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"))) %>%
  mutate_at(vars(Q2, Q50,Q97), exp)

Ws %>%
  ggplot(aes(x = Q50, y = Element)) +
  geom_point() +
  geom_errorbarh(aes(xmin = Q2, xmax = Q97), width = 0.3) +
  facet_grid(.~Column) +
  ylab("Matrix Element") +
  xlab("Interaction Rate") +
  theme(text = element_text(size=16))
