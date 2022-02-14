# Current Use of Contraception Among Married Women
# by Age, Education and Desire for More Children
# Fiji Fertility Survey, 1975
# Reference: Little, R. J. A. (1978). 
# Generalized Linear Models for Cross-Classified Data from the WFS. 
# World Fertility Survey Technical Bulletins, Number 5.
# see also: 
# Rodríguez, G. (2007). Lecture Notes on Generalized Linear Models.
# URL: http://data.princeton.edu/wws509/notes/

FujiFerilitySurvey<-function(){
  age=c("<25","<25","<25","<25",
        "25-29","25-29","25-29","25-29",
        "30-39","30-39","30-39","30-39",
        "40-49","40-49","40-49","40-49")
  age_num=c(20,20,20,20,
            27.5,27.5,27.5,27.5,
            35,35,35,35,
            45,45,45,45)
  education=c("low","low","high","high",
              "low","low","high","high",
              "low","low","high","high",
              "low","low","high","high")
  wantsMore=c("yes","no","yes","no",
              "yes","no","yes","no",
              "yes","no","yes","no",
              "yes","no","yes","no")
  notUsing=c(53,10,212,50,
             60,19,155,65,
             112,77,118,68,
             35,46,8,12)
  using=c(6,4,52,10,
          14,10,54,27,
          33,80,46,78,
          6,48,8,31)
  df=data.frame(age,age_num,education,wantsMore,notUsing,using)
  return(df)
}

