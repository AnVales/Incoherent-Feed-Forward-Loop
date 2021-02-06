#######################################
# incoherent type 1 feed forward loop #
#######################################

x<-0 # first TF
y<-0 # second TF
w<-0 # target gene
d1<-0.1 # degradation rate for x
d2<-0.01 # degradation rate for y
d3<-0.1 # degradation rate for w
ON<-10 # time when signal starts
OFF<-60 # time when signal ends

data<-matrix(0,150,3) # container

# signal
for (t in 1:150) {
  if (t>ON & t<OFF){
    swtch<-1
  }else{swtch<-0}
  
  # first transcription factor
  r<-(1-(x^2)/1000)*swtch 
  x<-x+r-d1*x
  data[t,1]<-x
  
  # second transcription factor
  r<-(0.3-(x)/100)*swtch 
  y<-y+r-d2*y
  data[t,2]<-y
  
  # target gene, regulated by the TFs
  r<-0.1 + x/(1+y^2)
  
  w<-w+r-d3*w
  data[t,3]<-w
  
}

# plot the network
plot(1:t, data[,1], col='blue', type='b', ylim=c(0,max(data,na.rm='T')), main="Incoherent Forward Loop Type I", xlab="Time", ylab="Protein level" )
lines(1:t, data[,2], col='red', type='b')
lines(1:t, data[,3], col='black', type='b')
legend(x = "topright", legend = c("TF1", "TF2","Target"), fill = c("blue", "red","black"), 
       title = "Network")
abline(v=ON, lty=2)
abline(v=OFF, lty=2)


###############################
# Simple regulation vs. IIFFL #
###############################


x<-0 # first TF
y<-0 # second TF
w<-0 # target gene 1
o<-0 # target gene 2
q<-0 # target gene 3
d1<-0.1 # degradation rate for x
d2<-0.01 # degradation rate for y
d3<-0.1 # degradation rate for w
ON<-10 # time when signal starts
OFF<-60 # time when signal ends

data<-matrix(0,150,5) # container

# signal
for (t in 1:150) {
  if (t>ON & t<OFF){
    swtch<-1
  }else{swtch<-0} 
  
  # first transcription factor
  r<-(1-(x^2)/1000)*swtch 
  x<-x+r-d1*x 
  data[t,1]<-x 
  
  # second transcription factor
  r<-(0.3-(x)/100)*swtch
  y<-y+r-d2*y
  data[t,2]<-y
  
  # target gene 1, regulated by TF1
  r<-0.1+(1-(x)/30)*swtch
    
    # with a denominator smaller than 17 the hump is smaller and the decay is slower but it decay earlier 
    # with a denominator greater than 17, the hump is bigger with the biggest denominators
    # the decay is very similar, because with a smaller hump, the target ends at the same time that the transcription factor
    # which have the samen decay curve
  
  w<-w+r-d3*w
  data[t,3]<-w
  
  # target gene 2, regulated by TF1
  r<-0.1+(1-(x)/500)*swtch

  o<-o+r-d3*o
  data[t,4]<-o
  
  # target gene 3, regulated by TF1
  r<-0.1+(1-(x)/10)*swtch

  q<-q+r-d3*q
  data[t,5]<-q
}

# plot the simple regulation network
plot(1:t, data[,1], col='blue', type='b', ylim=c(0,max(data,na.rm='T')), main="Simple Regulation", xlab="Time", ylab="Protein level" )
lines(1:t, data[,2], col='red', type='b')
lines(1:t, data[,3], col='mediumseagreen', type='b')
legend(x = "topright", legend = c("TF1", "TF2","Target"), fill = c("blue", "red","mediumseagreen"), 
       title = "Network")
abline(v=ON, lty=2)
abline(v=OFF, lty=2)

# plot 3 different simple regulation genes
plot(1:t, data[,3], col='mediumvioletred', type='b', ylim=c(0,max(data,na.rm='T')), main="Simple Regulation", xlab="Time", ylab="Protein level" )
lines(1:t, data[,4], col='mediumspringgreen', type='b')
lines(1:t, data[,5], col='mediumslateblue', type='b')
legend(x = "topright", legend = c("Target 1", "Target 2","Target 3"), fill = c("mediumvioletred", "mediumspringgreen","mediumslateblue"), 
       title = "Network")
abline(v=ON, lty=2)
abline(v=OFF, lty=2)


######################################################################
# Noise in Simple Regulation and incoherent type 1 feed forward loop #
######################################################################

# number of replicas
replicas <- 150

data<-matrix(0,150,4) # container to save the data

# define replicas loop
for (s in 1:replicas){
  
  x<-0 # first TF
  y<-0 # second TF
  w<-0 # target gene 1
  q<-0 # target gene 2
  d1<-0.1 # degradation rate for x
  d2<-0.01 # degradation rate for y
  d3<-0.1 # degradation rate for w
  ON<-10 # time when signal starts
  OFF<-60 # time when signal ends
  
  # define the times and the noise
  for (t in 1:150){
    noise<-rnorm(1,1,0.3)

  # signal
  if (t>ON & t<OFF){
    swtch<-1
  }else{swtch<-0}
    
    # first transcription factor
    r<-(1-(x^2)/1000)*swtch
    x<-x+r*noise-d1*x
    data[t,1]<-x
    
    # second transcription factor
    r<-(0.3-(x)/100)*swtch
    y<-y+r*noise-d2*y
    data[t,2]<-y
    
    # target gene 1, regulated by TF1 and TF2
    r<-0.1 + x/(1+y^2)
    
    w<-w+r*noise-d3*w
    data[t,3]<-w
    
    # target gene 1, regulated by TF1
    r<-0.1+(1-(x)/30)*swtch

    q<-q+r*noise-d3*q
    data[t,4]<-q
    
  }
  
  
}

# plot the noise effect in the network
plot(1:t, data[,1], col='blue', type='b', ylim=c(0,max(data,na.rm='T')), main="Incoherent Forward Loop Type I and Simple Regulation with noise", xlab="Time", ylab="Protein level" )
lines(1:t, data[,2], col='red', type='b')
lines(1:t, data[,3], col='black', type='b')
lines(1:t, data[,4], col='mediumseagreen', type='b')
legend(x = "topright", legend = c("TF1", "TF2","Target IIFFL","Target SR"), fill = c("blue", "red","black","mediumseagreen"), 
       title = "Network")
abline(v=ON, lty=2)
abline(v=OFF, lty=2)


#######################################################################
# Noise in Simple Regulation and incoherent type 1 feed forward loop: #
# BOXPLOT                                                             #
#######################################################################

# number of replicas
replicas <- 150

data<-matrix(0,replicas,4) # container

# define replicas loop
for (s in 1:replicas){
  
  x<-0 # first TF
  y<-0 # second TF
  w<-0 # target gene 1
  q<-0 # target gene 2
  d1<-0.1 # degradation rate for x
  d2<-0.01 # degradation rate for y
  d3<-0.1 # degradation rate for w
  ON<-10 # time when signal starts
  OFF<-60 # time when signal ends
  
  # define the times and the noise
  for (t in 1:62){ # 12, 22, 32, 42, 52, 62
    noise<-rnorm(1,1,0.3)
    
    if (t>ON & t<OFF){
      swtch<-1
    }else{swtch<-0}
    
    
    # first transcription factor
    r<-(1-(x^2)/1000)*swtch
    x<-x+r*noise-d1*x
    
    # second transcription factor
    r<-(0.3-(x)/100)*swtch
    y<-y+r*noise-d2*y
    
    # target gene 1, regulated by TF1 and TF2
    r<-0.1 + x/(1+y^2)
    w<-w+r*noise-d3*w
    
    # target gene 1, regulated by TF1
    r<-0.1+(1-(x)/30)*swtch
    q<-q+r*noise-d3*q
    
  }
  # save the data 
  data[s,1]<-x
  data[s,2]<-y
  data[s,3]<-w
  data[s,4]<-q
  
}

# boxplot of the dispersion in IIFFL and simple regulation
boxplot(data,  ylim = c(0, 15), main="Incoherent Forward Loop Type I and Simple Regulation with noise at time 62",col=c("blue", "red","black","mediumseagreen"),xlab="Transcription factors and genes",ylab="Dispersion", xaxt="n") #, xaxt="n"
legend(x = "topright", legend = c("TF1", "TF2","Target IIFFL","Target SR"), fill = c("blue", "red","black","mediumseagreen"),title = "Network")

