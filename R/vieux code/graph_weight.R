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

ggplot(BDD_a %>% filter( Tag_year == "2016"), aes(x = Year, y  = Size, color = Treatment))+
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


ggplot(BDD_a, aes(x = Year, y  = Size, color = Treatment))+
  geom_line(aes(group = Tag_id))

yg2 = data.frame(subset(t, !is.na(t[,1]))[,1])
yg2 = cbind(yg2, "2017" = 1,"2018"= 1,"2019"= 1,"2020"= 1,"2021"= 1,"2022"= 1)
for (i in 2:7){
  for (j in 1: nrow(yg2)){
    mean = yg2[j,i-1] + 0.274*(94-yg2[j,i-1])
    yg2[j,i] <- rtruncnorm(1, yg2[j,i-1],Inf ,mean, 15)
  }
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





L = 160
K = 1
t0 = 2012

x = c(2016:2021)
y = c(125)
for (t in x[-1]){
  y = append(y, y[length(y)]+K*(L-y[length(y)]))
}
y = L*(1-exp(-K*(x-t0)))


ggplot(BDD_a %>% filter( Tag_year =="2016"))+ #%>% filter(Year != 2022))+
  #geom_point()+
  geom_line(aes(group = Tag_id, x = Year, y  = Size, color = Treatment))+
  geom_point (data = data.frame(x,y), aes(x = factor(x), y = y))+
  facet_wrap(~Treatment)


yg2 = data.frame(subset(s, !is.na(s[,1]))[,1])
yg2 = cbind(yg2, "2017" = 1,"2018"= 1,"2019"= 1,"2020"= 1,"2021"= 1,"2022"= 1)
for (i in 2:7){
  for (j in 1: nrow(yg2)){
    mean = yg2[j,i-1] + 0.161*(161 - yg2[j,i-1])
    yg2[j,i] <- rtruncnorm(1, yg2[j,i-1],Inf ,mean, 15)
  }
}
colnames(yg2) <- c("2016","2017","2018","2019","2020","2021","2022")
#yg2 <- cbind(yg2,ind = c(1:nrow(yg2)))
yg2 <- pivot_longer(yg2, cols = colnames(yg2))

yg2 <- yg2 %>% arrange(name)
yg2 <- yg2 %>% mutate(name = factor(name, levels = c("2016","2017","2018","2019","2020","2021","2022")))
yg2 <- cbind(yg2, Tag_id = c(1:1926))



L = 200
K = 0.261
t0 = 2014.5


#y = c(125)
#for (t in x[-1]){
#  y = append(y, y[length(y)]+K*(L-y[length(y)]))
#}

x = c(2016:2021)
vB <- function(x,K,L,t0) {return (L*(1-exp(-K*(x-t0))))}


ggplot(BDD_a %>% filter(Tag_year =="2016" & Year != "2022"))+
  #geom_point()+
  #facet_wrap(~Treatment)+
  #geom_line(aes(group = Tag_id, x = Year, y  = Size, color = Treatment))+
  geom_point (data = data.frame(x,y=vB(x,0.617,189,2014.5), Treatment = 1), aes(x = factor(x), y = y), color = "red")+
  geom_point (data = data.frame(x,y=vB(x,0.756,180,2014.5), Treatment = 2), aes(x = factor(x), y = y), color = "green")+
  geom_point (data = data.frame(x,y=vB(x,0.85,163.4,2014.5), Treatment = 3), aes(x = factor(x), y = y), color = "blue")+
  geom_point (data = data.frame(x,y=vB(x,0.743,173.4,2014.5), Treatment = 4), aes(x = factor(x), y = y), color = "violet")
  


ggplot(BDD_a %>% filter(Year != "2022"))+ #Tag_year =="2016" & 
  #geom_point()+
  facet_wrap(~Treatment)+
  geom_line(aes(group = Tag_id, x = Year, y  = Size, color = Treatment))+
  geom_point (data = data.frame(x,y=vB(x,1.839,145,2014.5), Treatment = 1), aes(x = factor(x), y = y), color = "red")+
  geom_point (data = data.frame(x,y=vB(x,1.839,152,2014.5), Treatment = 2), aes(x = factor(x), y = y), color = "green")+
  geom_point (data = data.frame(x,y=vB(x,1.839,139.4,2014.5), Treatment = 3), aes(x = factor(x), y = y), color = "blue")+
  geom_point (data = data.frame(x,y=vB(x,1.839,151.4,2014.5), Treatment = 4), aes(x = factor(x), y = y), color = "violet")+
  geom_point (data = data.frame(x,y=vB(x,0.964,172,2014.5), Treatment = 1), aes(x = factor(x), y = y), color = "red")+
  geom_point (data = data.frame(x,y=vB(x,0.964,171,2014.5), Treatment = 2), aes(x = factor(x), y = y), color = "green")+
  geom_point (data = data.frame(x,y=vB(x,0.964,162.4,2014.5), Treatment = 3), aes(x = factor(x), y = y), color = "blue")+
  geom_point (data = data.frame(x,y=vB(x,0.964,166.4,2014.5), Treatment = 4), aes(x = factor(x), y = y), color = "violet")


ggplot(BDD_a %>% filter(Tag_year =="2016" & Year != "2022"))+
  #geom_point()+
  #facet_wrap(~Treatment)+
  #geom_line(aes(group = Tag_id, x = Year, y  = Size, color = Treatment))+
  geom_point (data = data.frame(x,y=vB(x,0.8,176,2014.5), Treatment = 1), aes(x = factor(x), y = y), color = "red")+
  geom_point (data = data.frame(x,y=vB(x,1,172,2014.5), Treatment = 2), aes(x = factor(x), y = y), color = "green")+
  geom_point (data = data.frame(x,y=vB(x,1.1,160,2014.5), Treatment = 3), aes(x = factor(x), y = y), color = "blue")+
  geom_point (data = data.frame(x,y=vB(x,2,200,2014.5), Treatment = 4), aes(x = factor(x), y = y), color = "violet")

ggplot(BDD_a %>% filter(Tag_year =="2016"  & Treatment %in% c(1,3)))+ #Tag_year =="2016" & 
  #geom_point()+
  #facet_wrap(~Lake)+
  geom_line(aes(group = Tag_id, x = Year, y  = log(Size), color = Treatment))

Aux <- BDD_a %>% filter(Tag_year =="2016") %>% group_by(Year, Treatment) %>% summarise (mean = mean(Size), max = max(Size), min = min(Size))

ggplot(Aux)+ #Tag_year =="2016" & 
  #geom_point()+
  #facet_wrap(~Lake)+
  geom_line(aes(group = Treatment, x = Year, y  = mean, color = Treatment))+
  geom_line(aes(group = Treatment, x = Year, y  = min, color = Treatment))+
  geom_line(aes(group = Treatment, x = Year, y  = max, color = Treatment))


ggplot(BDD_a %>% filter(Tag_year =="2016" & Year != "2022"))+ #Tag_year =="2016" & 
  #geom_point()+
  #facet_wrap(~Treatment)+
  geom_boxplot(aes(x = Year, y  = Size, color = Treatment))














