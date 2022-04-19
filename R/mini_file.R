t <- BDD_a %>% group_by(Year,Lake) %>%summarize(n()) %>% arrange(Lake) %>% filter(Year %in% c("2016","2017",2021,"2022"))

j = 1
chose = FALSE
for (i in 1 : (nrow(t)/2)){
  print(chose)
  if(chose){
    print(t[j,3]/t[j+1,3])
  }
  else {print(t[j+1,3]/t[j,3])}
  chose = !(chose)
  j = j+2
}
