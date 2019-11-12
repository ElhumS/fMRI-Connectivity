library(FIACH)
library(RNiftyReg)
library(pracma)
library(Matrix)
library(oro.nifti)
library(AnalyzeFMRI)

########################
######## Inputs ########
########################
sigthresh <- 0.001 
clusterthresh <- 10 
ROI_file <- 'E:/Patients/Sub_01/fMRI_Preprocessing/T1/t1.nii'
TR <- 1.99 # in seconds
noise_basis_file <- 'E:/Patients/Sub_01/fMRI_Preprocessing/rp/noise_basis6.txt'
folder_with_data <- 'E:/Patients/Sub_01/fMRI_Preprocessing/safilt_ufmri'
output_direc <- 'E:/Patients/Sub_01/fMRI_Analysis/'
template_output_file <- 'E:/Patients/Sub_01/fMRI_Preprocessing/safilt_ufmri/safilt_ufmri.nii'
wm_timecourse <- 'E:/Patients/Sub_01/fMRI_Preprocessing/WM/wm_timecourse.txt'
csf_timecourse <- 'E:/Patients/Sub_01/fMRI_Preprocessing/CSF/csf_timecourse.txt'

##########################
######## Get Data ########
##########################
confound_filename <- read.table(file = noise_basis_file,header = FALSE,sep = ' ')
files <- list.files(path = folder_with_data,pattern = '*.nii', full.names = TRUE)
func_files <- readNii(input = files, fourD = TRUE)
func_filename <- func_files[,,,1:dim(func_files)[4]]

###########################
###### Detrending #########
###########################
reshape <- matrix(func_filename, dim(func_filename)[1]*dim(func_filename)[2]*dim(func_filename)[3], byrow = F) # time x space
reshape2 <- t(reshape) ## time x space ##

d_data <- detrend(x = reshape2, tt = 'linear')
d_timecourse <- array(t(d_data), dim=c(dim(func_filename)[1], dim(func_filename)[2], dim(func_filename)[3], dim(func_filename)[4]))

#########################################
######### Regressing out Noise ##########
#########################################
xx <- matrix(nrow = dim(d_timecourse)[4],ncol = 3) #constructing a linear quatratic
xx[,1] <- array(1,dim(d_timecourse)[4])
xx[,2] <- t(1:dim(d_timecourse)[4])/dim(d_timecourse)[4]
xx[,3] <- t((1:dim(d_timecourse)[4])^2)/dim(d_timecourse)[4]^2

wm_filename <- read.table(file = wm_timecourse,header = FALSE,sep = ' ')
noise_wm_avg <- wm_filename
noise_wm_demean <- noise_wm_avg-colMeans(noise_wm_avg)

csf_filename <- read.table(file = csf_timecourse,header = FALSE,sep = ' ')
noise_csf_avg <- csf_filename
noise_csf_demean <- noise_csf_avg-colMeans(noise_csf_avg)

noise_avg <- confound_filename[1:dim(func_filename)[4],]
#noise_demean <- noise_avg-rowMeans(noise_avg)

X <- cbind(xx,noise_wm_demean,noise_csf_demean,noise_avg)
y <- d_timecourse

y_regress_ss <- function(y,X){   #given the same name as related function in matlab
  n <- dim(X)[1]
  ncolX <- dim(X)[2]
  signature(qr = "sparseQR") # makes it 'economy-size decomposition'
  z <- qr(X,LAPACK=T)
  Q <- qr.Q(z)
  R <- qr.R(z)
  perm <- z$pivot
  p <- sum(abs(diag(R)) > max(n,ncolX)*eps(R[1]))
  
  if (p < ncolX) { #previously used isTRUE()
    R <- R[1:p,1:p]
    Q <- Q[,1:p]
    perm <- perm[1:p]
  }
  
  b <- zeros(ncolX,1)
  b[perm] <- mldivide(A = R,B = (t(Q)%*%y))
  yhat <- as.matrix(X)%*%b
  r <- y-yhat
  return(r)
}

V0Cov <- array(0,dim(y))
for (i in 1:dim(y)[1]){
  for (j in 1:dim(y)[2]){
    for (k in 1:dim(y)[3]){
      res <- y_regress_ss(drop(y[i,j,k,]),X)
      V0Cov[i,j,k,] <- res
    }}}
V0Cov[is.nan(V0Cov)]=0

#########################################
#### Bandpass filter data 0.01-0.1Hz ####
#########################################

## high pass ##
hp <- array(0,c(dim(func_filename)[1],dim(func_filename)[2],dim(func_filename)[3],dim(func_filename)[4]))
for (i in 1:dim(V0Cov)[1]){
  for (j in 1:dim(V0Cov)[2]){
    for (k in 1:dim(V0Cov)[3]){
      hp[i,j,k,] <- highPass(V0Cov[i,j,k,],freq = 100,tr = TR)
    }
  }
}

## low pass ##
kern <- kaiserWin(fl=0.1,tw=.025,sf=1/TR,d.sa=70,d.pbr=.1,type="low")
hplp_d_timecourse <- array(0,c(dim(hp)[1],dim(hp)[2],dim(hp)[3],dim(hp)[4]))
for (i in 1:dim(hp)[1]){
  for (j in 1:dim(hp)[2]){
    for (k in 1:dim(hp)[3]){
      hplp_d_timecourse[i,j,k,] <- convolve1d(x=hp[i,j,k,],fir=kern,subtractMed=TRUE)
    }
  }
}

#########################################
######## Find timecourse of ROI #########
#########################################

## Extract Timecourse ##
ROI_transform<- readNii(input = ROI_file, fourD = F)

## Compute average timecourse of ROI ##
nTimePoints <- dim(hplp_d_timecourse)[4]
tc <- array(0,c(1,nTimePoints))
dx <- dim(hplp_d_timecourse)[1]
dy <- dim(hplp_d_timecourse)[2]
dz <- dim(hplp_d_timecourse)[3]

for (i in 1:dx) {
  for (j in 1:dy) {
    for (k in 1:dz) {
      if (ROI_transform[i,j,k]>0) {
        tc <- tc + t(drop(hplp_d_timecourse[i,j,k,]))
      }
    }}}

tc_focus <- tc/nnzero(ROI_transform)

################################################ 
######## Perform Seed Based Correlation ######## 
################################################ 

## compute correlation
### positive
final_cor_pos <- array(0,c(dim(hplp_d_timecourse)[1],dim(hplp_d_timecourse)[2],dim(hplp_d_timecourse)[3],dim(hplp_d_timecourse)[4]))
for (i in 1:dx){
  for (j in 1:dy){
    for (k in 1:dz){
      x1 <- cor.test(x = t(tc_focus),y = hplp_d_timecourse[i,j,k,],alternative = 'greater',na.action = 'na.omit')
      final_cor_pos[i,j,k,] <-x1$p.value
    }}}
### negative
final_cor_neg <- array(0,c(dim(hplp_d_timecourse)[1],dim(hplp_d_timecourse)[2],dim(hplp_d_timecourse)[3],dim(hplp_d_timecourse)[4]))
for (i in 1:dx){
  for (j in 1:dy){
    for (k in 1:dz){
      x2 <- cor.test(x = t(tc_focus),y = hplp_d_timecourse[i,j,k,],alternative = 'less',na.action = 'na.omit')
      final_cor_neg[i,j,k,] <-x2$p.value 
    }}}

###################################### 
######## Significance Testing  #######
######################################
## Find out where pvalues are <sigthresh ##
sig_fz_cor_pos <- final_cor_pos[,,,1]#not looking over time/dynamic, so (x,y,z,1)
for (i in 1:dim(sig_fz_cor_pos)[1]) {
  for (j in 1:dim(sig_fz_cor_pos)[2]) {
    for (k in 1:dim(sig_fz_cor_pos)[3]) {
      if (is.na(final_cor_pos[i,j,k,1])) {
        sig_fz_cor_pos[i,j,k] <- 0} 
    }}}

for (i in 1:dx) {
  for (j in 1:dy) {
    for (k in 1:dz) {
      if (isTRUE(final_cor_pos[i,j,k,1]<=sigthresh)) {
        sig_fz_cor_pos[i,j,k] <- 1 
        print("Yup")} 
      else{
        sig_fz_cor_pos[i,j,k] <- 0
        print("Nope")}
    }}}

sig_fz_cor_pos_cluster_thresh <- cluster.threshold(x = sig_fz_cor_pos,size.thr = clusterthresh)


sig_fz_cor_neg <- final_cor_neg[,,,1]#not looking over time/dynamic, so (x,y,z,1)
for (i in 1:dim(sig_fz_cor_neg)[1]) {
  for (j in 1:dim(sig_fz_cor_neg)[2]) {
    for (k in 1:dim(sig_fz_cor_neg)[3]) {
      if (is.na(final_cor_neg[i,j,k,1])) {
        sig_fz_cor_neg[i,j,k] <- 0} 
    }}}

for (i in 1:dx) {
  for (j in 1:dy) {
    for (k in 1:dz) {
      if (isTRUE(final_cor_neg[i,j,k,1]<=sigthresh)) {
        sig_fz_cor_neg[i,j,k] <- 1
        print("Yup")} 
      else{
        sig_fz_cor_neg[i,j,k] <- 0
        print("Nope")}
    }}}

sig_fz_cor_neg_cluster_thresh <- cluster.threshold(x = sig_fz_cor_neg,size.thr = clusterthresh)


#############################
##### Writing out Files #####
#############################
template<-template_output_file
output <- paste(output_direc,"fc_pos_cluster_thresh",sep="")                                             
hmma<-updateNifti(sig_fz_cor_pos_cluster_thresh,template=template)
writeNifti(image = hmma,file = output,template = template,datatype = 'float')

template<-template_output_file
output <- paste(output_direc,"fc_neg_cluster_thresh",sep="")                                             
hmma<-updateNifti(sig_fz_cor_neg_cluster_thresh,template=template)
writeNifti(image = hmma,file = output,template = template,datatype = 'float')
