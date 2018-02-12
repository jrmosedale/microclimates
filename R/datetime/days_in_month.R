# Calculates days in month (from 1980-2015) from day, month, year
# Used to write monthly files of imputed values 
days.in.month<-function(day,month,year){
  y<-year-1979 # so 1=1980
  feb.d<-c(29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28,29,28,28,28,
           29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28) # days of Feb from 1980 to 2015
  monthdays<-c(31,feb.d[y],31,30,31,30,31,31,30,31,30,31)
  return(monthdays[month])
}
