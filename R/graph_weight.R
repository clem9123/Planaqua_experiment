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

ggplot(BDD_a %>% filter( Tag_year == "2016" & Year == "2021"))+
  geom_density(aes(x = Weight, color = Year))+
  geom_density(aes(x= rnorm(length(Weight),15,8)))

L = 71.830
K = 0.666
t0 = 2015.5

x = c(2016:2022)
y = L*(1-exp(-K*(x-t0)))
et = c(8,10,13,16,15,16,16)
#ggplot(data = data.frame(x,y), aes(x = factor(x), y = y))+
#  geom_point()

ggplot(BDD_a %>% filter(Tag_year == "2016"))+ #Year !="2022" & 
  #geom_point()+
  geom_line(aes(x = Year, y  = Weight, color = Treatment, group = Tag_id))+
  geom_point (data = data.frame(x,y), aes(x = factor(x), y = y))#+
  #facet_wrap(~Lake)

ggplot(BDD_a %>% filter(Tag_id %in% tag))+ #Year !="2022" & 
  #geom_point()+
  geom_line(aes(x = Year, y  = Weight, color = Treatment, group = Tag_id))+
  geom_point (data = data.frame(x,y), aes(x = factor(x), y = y))+
  geom_errorbar(data = data.frame(x,y,var), aes(x = factor(x), y = y, ymax =y+et,ymin = y-et))

ggplot(BDD_a %>% filter(Tag_year == "2016"))+ #Year !="2022" & 
  #geom_point()+
  geom_line(aes(x = Year, y  = Weight, group = Tag_id, color = Treatment))+
  geom_point (data = data.frame(x,y), aes(x = factor(x), y = y))+
  geom_errorbar(data = data.frame(x,y,var), aes(x = factor(x), y = y, ymax =y+et,ymin = y-et))+
  facet_wrap(~ Tag_id %in% tag)

etg = c(8,8,13,9.5,11,8,8)
yg = c(19)
for (i in 1:6){
  yg <- append(yg, yg[length(yg)] + 0.357*(78-yg[length(yg)]))
}

ggplot(BDD_a %>% filter(Tag_id %in% tag))+ #Year !="2022" & 
  #geom_point()+
  geom_line(aes(x = Year, y  = Weight, color = Treatment, group = Tag_id))+
  geom_point (data = data.frame(x,yg), aes(x = factor(x), y = yg))+
  geom_errorbar(data = data.frame(x,yg,var), aes(x = factor(x), y = yg, ymax =yg+etg,ymin = yg-etg))



yg2 = data.frame(subset(t, !is.na(t[,1]))[,1])
for (i in 1: 6){
  yg2 <- cbind(yg2,0)
  yg2 <- apply(yg2, rnorm(1,yg2[,ncol(yg2)] + 0.274*(94-yg2[,ncol(yg2)]), 15))
}
colnames(yg2) <- c("2016","2017","2018","2019","2020","2021","2022")
#yg2 <- cbind(yg2,ind = c(1:nrow(yg2)))
yg2 <- pivot_longer(yg2, cols = colnames(yg2))

yg2 <- yg2 %>% arrange(name)
yg2 <- yg2 %>% mutate(name = factor(name, levels = c("2016","2017","2018","2019","2020","2021","2022")))
yg2 <- cbind(yg2, Tag_id = c(1:1926))

ggplot(BDD_a %>% filter(Tag_id %in% tag))+ #Year !="2022" & 
  #geom_point()+
  geom_line(aes(x = Year, y  = Weight, color = Treatment, group = Tag_id))+
  geom_point (data = data.frame(x,yg2), aes(x = factor(x), y = yg2))+
  geom_errorbar(data = data.frame(x,yg2,var), aes(x = factor(x), y = yg2, ymax =yg2+15,ymin = yg2-15))


ggplot(BDD_a %>% filter(Tag_year == "2016" ))+ #Year !="2022" & 
  #geom_point()+
  geom_line(aes(x = Year, y  = Weight, color = Treatment, group = Tag_id))+
  geom_line (data = yg2, aes(x = factor(name), y = value ,group = Tag_id))

ggplot(yg2)+ #Year !="2022" & 
  #geom_point()+
  #geom_line(aes(x = Year, y  = Weight, color = Treatment, group = Tag_id))+
  geom_line (data = yg2, aes(x = factor(name), y = value ,group = Tag_id))

ggplot(yg2 %>% filter(name == "2016"))+
  geom_density(aes( x = value))













