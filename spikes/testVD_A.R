
N = 100000

treatment = rnorm(N,12,2)
control = rnorm(N,13,3)
d = c(treatment,control)
f = factor(rep(c("Treatment","Control"),c(length(treatment),length(control))),
           levels=c("Treatment","Control"),ordered=T)


res=list()
for( i in 1:30){
res$fact = c(res$fact,system.time({  
  mask = f==levels(f)[1]
  r = sum(rank(d)[mask])
  #r = aggregate(rank(d),by=list(f),FUN=sum)[1,2]
  m = sum(mask)
  n = sum(!mask)
  A = (r1/m - (m+1)/2)/n
  #n = table(f)
  #A = (r/n[1] - (n[1]+1)/2)/n[2]
})["elapsed"])

res$rng = c(res$rng,system.time({
  tc = split(d,f)
  treatment = tc[[1]]
  control = tc[[2]]
  # Compute the rank sum
m = length(treatment)
n = length(control)
r = rank(c(treatment,control))
r1 = sum(r[seq_len(m)])

# Compute the measure
A = (r1/m - (m+1)/2)/n
})["elapsed"])

res$blog=c(res$blog,system.time({
  # Compute the rank sum (Eqn 13)
r = rank(c(treatment,control))
r1 = sum(r[seq_along(treatment)])

# Compute the measure (Eqn 14) 
m = length(treatment)
n = length(control)
A = (r1/m - (m+1)/2)/n
})["elapsed"])
}

par(mar=c(3,3,.5,.5))
boxplot(res,ylim=c(0,0.08))

kruskal.test(res)
wilcox.test(res$fact,res$rng)
wilcox.test(res$rng,res$blog)
