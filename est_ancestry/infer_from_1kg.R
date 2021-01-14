suppressMessages(library(caret))
suppressMessages(library(argparser))
p = arg_parser('Read in the plink .eigenvec & .eigenval files produced by merge_1kg_geno_PCA.sh and use the known ancestry associated with the 1kg samples to assign ancestry to the merged in data we supplied.')
p = add_argument(p, 'pcadir', help='Base path containing the PCA results')
p = add_argument(p, 'pcafile', help='Base file name of PCA results')
p = add_argument(p, 'kgpopfile', help='1kg phase 3 annotation file containing the following 3 columns: FID	IID	ancestry. The ancestry values in column 3 are the labels to estimate in our data. For files prepared by me, the name of that column may vary, so this script reassigns the column 3 name to \'ancestry\'.')
p = add_argument(p, '--centerscale', type='flag', flag=T, help='Flag to preprocess principal components by centering and scaling before model training. FALSE if not provided, which defaults to no preprocessing of PCA input.')
args = parse_args(p)

pcadir = args$pcadir
pcafile = args$pcafile
kgpopfile = args$kgpopfile
cs = args$centerscale

setwd(pcadir)
options(scipen=100, digits=3)

##############################################################

# read in the eigenvectors, produced in PLINK
eigenvec = data.frame(read.table(paste0(pcadir,'/',pcafile,'.eigenvec'), header=FALSE, skip=0, sep=' '))
eigenval = t(data.frame(read.table(paste0(pcadir,'/',pcafile,'.eigenval'), header=FALSE, skip=0, sep=' ')))
rownames(eigenvec) = eigenvec[,2]
eigenvec = eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) = colnames(eigenval) = paste0('PC', c(1:20))

# read in the 1kg metadata & read population from 3rd column (col name varies between metadata files I prepared)
ped = data.frame(read.table(kgpopfile, header=TRUE, skip=0, sep='\t'))
colnames(ped)[3] = 'ancestry'
ped = ped[which(ped$IID %in% rownames(eigenvec)), ]
rownames(ped) = ped$IID

# Read in metadata
# Add ancestry column to eigenvalues dataframe
eigenvec = cbind(eigenvec, ped[, 'ancestry'][match(rownames(eigenvec), rownames(ped))])
colnames(eigenvec)[21] = 'ancestry'

# Split eigenvectors file into 1KG (training) and our data (testing) datasets
test = eigenvec[is.na(eigenvec['ancestry']),]
train = eigenvec[!is.na(eigenvec['ancestry']),]

# Train the KNN model
if(cs){
  print('Centering and scaling PCA input...')
  preprocess_method=c('center', 'scale')
} else {
  preprocess_method=c()
}
trctrl = trainControl(method = 'repeatedcv', number = 10, repeats = 4)
knn_fit = train(ancestry ~., data = train, method = 'knn',
 trControl = trctrl,
 preProcess = preprocess_method,
 tuneLength = 10)

# Get predictions for the new data
est = predict(knn_fit, newdata = test)
# Add estimated ancestry to test data and write out results
test$ancestry = est
write.table(test,paste0(pcadir,'/',pcafile,'.knn.ancestry_est'),sep='\t',quote=F)
