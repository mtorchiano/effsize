#Setup some dummy data
set.seed(1234)
numSubjects <- 20
subject <- 1:numSubjects
before <- runif(numSubjects)
after <- before - 0.35*runif(length(before))
df <- data.frame(
  subject=rep(subject,2),
  group=c(
    rep("before",numSubjects),
    rep("after",numSubjects)
  ),
  value=c(before, after)
)

##################################

valueBefore <- df[df$group == "before",]$value
valueAfter <- df[df$group == "after",]$value

# d1 == d2
d1 <- effsize::cohen.d(valueAfter, valueBefore, paired=T)
d2 <- effsize::cohen.d(value~group, data=df, paired=T)

##################################

# Sorted by subjectId
df <- df[order(df$subject),]
valueBefore <- df[df$group == "before",]$value
valueAfter <- df[df$group == "after",]$value

# d1==d2==d3 is the same: d estimate: -0.5307307 (medium)
# but d4 is not: -0.145407 (negligible)

d3 <- effsize::cohen.d(valueAfter, valueBefore, paired=T)
d4 <- effsize::cohen.d(value~group, data=df, paired=T)

##################################

# Sorted by group
df <- df[order(df$group),]
valueBefore <- df[df$group == "before",]$value
valueAfter <- df[df$group == "after",]$value

# d1 == d2 == d3 == d5 == d6 are the same (although still false due to #49)
d5 <- effsize::cohen.d(valueAfter, valueBefore, paired=T)
d6 <- effsize::cohen.d(value~group, data=df, paired=T)

##################################
