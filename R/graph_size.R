ggplot(BDD_a, aes (x = Year, y = Weight, color = Tag_year))+
  geom_line(aes(group = Tag_id))+
  geom_point()+
  facet_wrap(~ Tag_year)

ggplot(BDD_a %>% filter (Tag_year == "2016" & Year != "2022"), aes (x = Year, y = Size))+
  geom_line(aes(group = Tag_id))+
  facet_wrap(~ Lake)

ggplot(BDD_a %>% filter (Tag_year == Year), aes (x = Weight, color = Tag_year))+
  geom_density()

ggplot(BDD_a %>% filter(Tag_year == "2016"), aes(x =Size , color = Year))+
         geom_density()
