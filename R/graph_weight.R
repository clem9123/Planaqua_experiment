ggplot(BDD_a %>% filter(Year == "2019" & Tag_year == "2019"), aes ( x = Weight))+
  geom_bar(aes(x = Weight, color = Weight == 58))

ggplot(BDD_a %>% filter(Year == "2019"), aes ( x = Weight))+
  geom_bar(aes(x = Weight))+
  xlim(60,200)

ggplot(BDD_a %>% filter( Tag_year == "2016"))+
  geom_density(aes(x  = Log_size, color = Year))

ggplot(BDD_a %>% filter( Tag_year == "2016"))+
  geom_boxplot(aes(x = Year, y = Log_weight, color = Treatment))

ggplot(BDD_a %>% filter(Tag_year == "2019"))+
  geom_boxplot(aes(x = Year, y = Log_weight, color = Treatment))

ggplot(BDD_a %>% filter( Tag_year == "2016"), aes(x = Year, y  = Weight, color = Treatment))+
  #geom_point()+
  geom_line(aes(group = Tag_id))+
  facet_wrap(~Treatment)
