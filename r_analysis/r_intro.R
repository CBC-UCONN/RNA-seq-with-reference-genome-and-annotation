#############################################
# a really incredibly brief introduction to R. 
#############################################

# this is an R script. 

# lines beginning with "#" are comment lines. R will ignore the contents. 
# you can use these lines to document your code. 

# you can execute lines from this script interactively by entering them in the R console (bottom left pane in RStudio)
  # you can copy and paste them, or press command-enter


#############################################
# objects
#############################################

# everything in R is an object. there are several kinds of objects. 
# we can assign values to objects using the assignment operator "<-" like this:

a <- 5

# you can see the contents of many objects by simply entering them on the console

a

# you can do math in the console and assign results to objects

5 + 5

a <- 5 + 5

a

# you can update objects

a <- a + 5

a

# there are many kinds of objects:

# vectors:
a <- c(1,2,3,4,5)
b <- c("bear","dog","fox","cat","rat")
a
b

  # access subsets of vectors like this:
  a[2]
  b[2:4]
  
# lists:
g <- list(x=a,y=b)
g

  # access parts of lists like this
  g[[2]]

# arrays (or matrices)
d <- matrix(nrow=5,ncol=2,data=1:10)
d
  
  # access parts of matrices like this:
  d[1,]
  d[,1]
  d[1,2]
  
# data frames
e <- data.frame(x=a,y=b)
e

  # access parts of data frames like matrices, but also:
  e$x

# R objects can contain several basic data types:

# numeric: 1,2,3...
# character: "bear","dog","fox"...
# logical: TRUE/FALSE
# factor (categorical): "bear","dog","fox"...


# R does *vector* arithmetic!

a + a
a^2
a + 100

# TRUE/FALSE vectors are super useful

f <- c(TRUE,FALSE,TRUE,FALSE,TRUE)

a[f]
e[f,]
d[!f,]

# you can create logical vectors like this:
f <- d[,1] > 2

# you can do arithmetic on logical vectors

TRUE + 1
FALSE + 1

#############################################
# functions
#############################################

# functions are objects that contain code
# you use them to manipulate other objects

# the function is()
# you can ask what an object is:

is(a)

is(d)

is(is)

# you can often see the code in a function by typing its name:

is

# you can write your own functions
x10 <- function(x){
  z <- x * 10
  return(z)
  }

x10(10)

# most functions have documentation!

?sum

sum(a)

#############################################
# analysis and plotting
#############################################

# an absurdly simple example. 
# simulate some data

x <- rnorm(n=300,mean=10,sd=10)
y <- x * 5 + rnorm(n=300,mean=25,sd=10)

z <- data.frame(x,y)

# plot the data
hist(z$x)
hist(z$y)
plot(z$x,z$y)

# fit a linear model
model <- lm(y ~ x, z)

# summarize linear model
summary(model)

# add trendline to plot
abline(model,col="red")


#############################################
# tidyverse and "pipes"
#############################################

# above we have used 'base' R. a very commonly used set of packages
# fall under the umbrella of the 'tidyverse' project. 
# in the DE analysis we'll use some tidyverse functions,
# but one major piece of syntax to be aware of is the "pipe"

library(tidyverse)

sum(a)

a %>% sum()

a %>% sum() %>% sqrt()

w <- a %>% 
    sum() %>% 
    sqrt() %>% 
    (function(x){x^5})

w
