size_break <- c(0,120,180,300)
BDD_a <- BDD_a %>% mutate (s_group = discretize(Size,method = "fixed",
                                                   breaks = size_break, 
                                                   labels= c(1:(length(size_break)-1))))
a <- BDD_g %>% filter (Year == "2021") %>% select(Tag_id,s_group,Size,Weight)
colnames(a) <- c("Tag_id","sg_2021","size_2021", "Weight")
b <- BDD_g %>% filter (Year == "2022") %>% select(Tag_id,s_group,Size,Weight)
colnames(b) <- c("Tag_id","sg_2022","size_2022", "Weight")

c <- merge(a,b,all=T)
nrow(c)

b %>% filter(Tag_id == "no_tag" & sg_2022 == 1) %>% mutate(sg_2021 = NA, size_2021 = NA) %>% 
  select("Tag_id","sg_2021","size_2021","sg_2022","size_2022")

c <- rbind(c,b %>% filter(Tag_id == "no_tag" & sg_2022 == 1) %>% 
             mutate(sg_2021 = NA, size_2021 = NA) %>% 
             select("Tag_id","sg_2021","size_2021","sg_2022","size_2022"))

ggplot(c ,aes(axis1 = sg_2021, axis2 = sg_2022)) +
  geom_alluvium(aes(fill = sg_2022),width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("2021", "2022"), expand = c(.05, .05))
