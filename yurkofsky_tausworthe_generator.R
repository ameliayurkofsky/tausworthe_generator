#ISYE 6644
#Summer 2021
#Mini-Project #2

#A Yurkofsky
#July 19 2021

####INPUTS YOU CAN SET####
#based on primitive polynomial x^p + x^r + 1
p = 9
r = 4
l = 8

#Derived variables
q = p-r
B = l * 10000 #size of simulation

####CREATE PRNs####

#Initialize list of 1s for B1, B2, B3, ..., Bq
initial_list <- rep(1, p) 

#Find Bi and append to the list for B reps
for (i in 1:B){
  nbit_i <- length(initial_list) + 1
  nbit <- ifelse(initial_list[nbit_i - p] == initial_list[nbit_i - p+q], 0, 1)
  initial_list <- append(initial_list, nbit)
}
initial_list <- initial_list[(p+1):length(initial_list)] #remove initial list of p bits

#First make a list of lists to store bits
bits_list <- split(initial_list, ceiling(seq_along(initial_list) / l))

#Then turn bits into base10 using 'unbinary' package in compositions library
library(compositions)
PRN_list <- NULL #initial list to hold PRNs
for (i in 1:length(bits_list)){
  PRN <- (unbinary(paste(bits_list[[i]], collapse=""))) / 2^l
  PRN_list <- append(PRN_list, PRN)
}
###PLOTS###
par(mfrow=c(1,2))
#hist
hist(PRN_list, freq = FALSE, 
     main = paste0("Histogram of PRNs (B=",B,", r=",r,", q=",q,", l=",l,")"),
     xlab = "PRNs")

#plot adjacent PRNs on unit square 
#turn PRNs into two vectors with coordinates (Ui, Ui+1), (Ui+1, Ui+2), etc.
x <- PRN_list[1:(length(PRN_list) -1)]
y <- PRN_list[-1]
plot(x, y, main = paste0("Adjacent PRNs in Unit Square (B=",B,", r=",r,", q=",q,", l=",l,")"))

####STATISTICAL TESTS ON PRNs#####

#goodness of fit test
set.seed(1124)
k = sample(3:51,1) #randomly select k

#create thresholds for cuts based on value of k
cuts <- NULL #initialize list to hold cut values
for (i in 1:(k-1)){
  cuts <- c(cuts, i/k)
}

#Create list of n per k interval
n_list <- length(subset(PRN_list, PRN_list<cuts[1]))
stop <- length(cuts) -1
for (i in 1:stop){
  n1x <- length(subset(PRN_list, PRN_list>=cuts[[i]] & PRN_list<cuts[[i+1]]))
  n_list <- append(n_list,  n1x)
}
n_list <- append(n_list, length(subset(PRN_list, PRN_list>=cuts[stop+1])))

#Obtain expected value of n
E = length(PRN_list)/k

#calculate test-statistic
test_stat_list <- NULL #initialize list
for (i in 1:length(n_list)){
  test_stat_list <- append(test_stat_list, ((n_list[[i]] - E)^2) / E)
}

#reject null hypothesis that PRNs are uniformly distributed? 
#"TRUE" = reject
#"FALSE" = fail to reject
test_statistic <- sum(test_stat_list)
critical_value <- qchisq(p=0.05, df=(k-1), lower.tail=FALSE)
test_statistic > critical_value

#tests of independence
#run test 1: above and below the mean
mean_runlist <- NULL #initialize
for (i in 1:length(PRN_list)){
  mean_runlist <- append(mean_runlist, ifelse(PRN_list[i] >=0.5, "+", "-"))
}
#count runs
mean_runs <- rle(mean_runlist) 
Btest <- length(mean_runs$lengths) 

#obtain expected value and variation of B
n = length(PRN_list)
n1 = length(subset(PRN_list, PRN_list>=0.5))
n2 = n - n1

E_B <- ((2*n1*n2)/n) + 0.5
E_Bv <- ((2*n1*n2)*((2*n1*n2)-n))/((n^2)*(n-1))

#calculate test statistic
test_statistic2 <- (Btest - E_B)/(sqrt(E_Bv))

#set critical value give alpha = 0.05
critical_value2 <- 1.96

#reject null hypothesis that PRNs are independent? 
#"TRUE" = reject
#"FALSE" = fail to reject
test_statistic2 > abs(critical_value2)

#run test 2: up and down
updown_runlist <- NULL #initialize
for (i in 1:length(PRN_list)){
  updown_runlist <- append(updown_runlist, ifelse(PRN_list[i+1] < PRN_list[i], "-", "+"))
}
#count runs
updown_runs <- rle(updown_runlist) 
A_updown <- length(updown_runs$lengths) 

#calculate expected value and variance of A
E_A <- ((2*(length(PRN_list)))-1)/3
E_Av <- ((16*(length(PRN_list)))-29)/90

#calculate test statistic
test_statistic3 <- (A_updown - abs(E_A)) / sqrt(E_Av)

#reject null hypothesis that PRNs are independent? 
#"TRUE" = reject
#"FALSE" = fail to reject
abs(test_statistic3) > 1.96

### Nor(0,1) deviates
#generate a few Nor(0,1) deviates (anyway you want) using Unif(0,1)'s from your Tausworthe generator
norm_list <- NULL
for (i in 1:length(PRN_list){
  z1 <- (((-2) * log(PRN_list[i])))^(.5) * (cosine(2 * 3.1415 * PRN_list[i+1]))
  norm_list <- append(norm_list, z1)
}
  