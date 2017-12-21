require(e1071)
require(neuralnet)
require(genalg)

############
# 
# Search for the following UPPERCASE subjects to see its contents:
# 
# GET_DATASET - functions to read and organize the datset
# CLASSIF_TOOLS - classifiers and cross-validation
# GRID_SEARCH - Implementation of grid search
# GA_FEATURE_SELECTION - Genetic algorithm for feature selection
# PLOTTING - Plotting subspaces and the SVM hyperplane
# 
############

## GET_DATASET
dataDir = "./testFolder"
datasetFile = "merged_dataset_moments.csv"

getData <- function(dir,exc=c()){
	dataset = list()
	data = read.csv(paste(c(dir,datasetFile),collapse="/"))
	
	#without patient identification
	#dataset$x = data[,1:ncol(data)-1]
	#dataset$y = data[,ncol(data)]

	dataset$x = data[,1:(ncol(data)-2)]
	dataset$y = data[,ncol(data)-1]
	dataset$patient = data[,ncol(data)]
	dataset$pat = substr(dataset$patient,1,10)

	if(length(exc)==0){
		exc = c(
			"meanpy",
			
			#"meanvx","meanvy",
			
			"meanfy",
			"sdevpx","sdevpy","sdevvx",
			"skewpy","skewvx",
			"kurtpx","kurtpy","kurtvx","kurtvy","kurtfx","kurtfy"
		)
	}

	#if(length(exc)==0){
	#	exc = c("meanpy","sdevpy","kurtpy","meanpx","meanvy","skewpy","meanfy","kurtfx","kurtpx","sdevpx","skewvx")
	#}
	#exc = c()

	colsToExc = colnames(dataset$x) %in% exc
	dataset$x = dataset$x[,!colsToExc]

	return(dataset)
}

## End of GET_DATASET


## CLASSIF_TOOLS

getModel <- function(dataset, seed=0){
	set.seed(seed)
	model = svm(x=dataset$x, y=dataset$y,
		kernel="linear",
		#kernel="poly",
		#degree=2,
		#coef0=1.0,
		#gamma=0.096,

		#cost=0.4,
		cost=0.866667,

		cross=10
	)

	return(model)
}

getModelForPCV_SVM <- function(x, y, params=c()){
	if(class(params) == "list"){
		model = svm(
			x=x,
			y=y,

			kernel=params$kernel,
			degree=params$degree,
			coef0=params$coef0,
			gamma=params$gamma,

			cost=params$cost,

			cross=0
		)
		return(model)
	} else {
		model = svm(
			x=x,
			y=y,
			
			kernel="linear",
			##kernel="poly",
			##degree=2,
			##coef0=1.0,
			##gamma=0.096,

			cost=0.4,
			#cost=0.86667, #grid search with all features

			cross=0
		)
		return(model)
	}
}

getModelForPCV_MLP <- function(x, y, params=c()){
	fo = as.formula(paste("aside ~",paste(colnames(x),collapse=" + ")))
	theData = cbind(x,2-as.numeric(y))
	colnames(theData) <- c(colnames(x),"aside")

	#resilient bacpropagation rprop
	if(class(params)=="list"){
		model = neuralnet(formula=fo, data=theData, hidden=1, linear.output=FALSE,
			threshold=params$threshold,
			stepmax=1e8, #max epochs
			lifesign="full", #verbose
			err.fct=params$err.fct
		)

		return(model)
	} else {
		model = neuralnet(formula=fo, data=theData, hidden=1, linear.output=FALSE,
			threshold=1.3e-3,
			#threshold=0.1,
			stepmax=1e8, #max epochs
			lifesign="full", #verbose
			err.fct="ce" #cross-entropy
		)

		return(model)
	}
}

getModelInfoPatientCV <- function(dataset, seed=0, getModelFunc, params=c()){
	set.seed(seed)
	modelInfo = list()
	modelInfo$folds = unique(dataset$pat)
	modelInfo$accs = c()
	modelInfo$percs = c()

	for(p in 1:length(modelInfo$folds)){
		remainingInstances = which(dataset$pat != modelInfo$folds[p])
		model = getModelFunc(
			x=dataset$x[remainingInstances,],
			y=dataset$y[remainingInstances],
			params=params
		)
		#model = svm(
		#	x=dataset$x[remainingInstances,],
		#	y=dataset$y[remainingInstances],

		#	kernel="linear",
		#	#kernel="poly",
		#	#degree=2,
		#	#coef0=1.0,
		#	#gamma=0.096,

		#	cost=0.4,

		#	cross=0
		#)

		xTest = dataset$x[-remainingInstances,]
		yTest = dataset$y[-remainingInstances]

		correct = length(which(doPredict(model,xTest) == yTest))
		modelInfo$accs = c(modelInfo$accs, correct/nrow(xTest))
		modelInfo$percs = c(modelInfo$percs, nrow(xTest))
		modelInfo$model = getModelFunc(x=dataset$x,y=dataset$y,params=params)

	}

	#modelInfo$acc = mean(modelInfo$accs)
	modelInfo$acc = sum(modelInfo$accs*modelInfo$percs)/nrow(dataset$x) #mean weighting the fold size
	modelInfo$accDisp = sum(modelInfo$percs*(modelInfo$acc - modelInfo$accs)^2)*(sum(modelInfo$percs)/(sum(modelInfo$percs)^2-sum(modelInfo$percs^2))) #variance weighting the fold size (dispersion measure)

	return(modelInfo)
}

doPredict <- function(model, input){
	if(class(model)=="nn"){
		preds = c()
		for(pr in 1:nrow(input)){
			pred = round(compute(model,input[pr,])$net.result[[1]])
			if(pred=="0"){
				preds = c(preds,"esquerdo")
			} else {
				preds = c(preds,"direito")
			}
		}
		#cat("PRED: ",preds,"\n")
		return(preds)
	} else {
		#cat("PRED: ",predict(model,input),"\n")
		return(predict(model, input))
	}
}

# CALL THIS FUNC!
getAcc <- function(){
	data = getData(dataDir)
	cat(ncol(data$x)," features: ",colnames(data$x),"\n")
	cat("Acc (cross-validation): ",getModel(data)$tot.accuracy,"\n\n")
}

# OR CALL THIS FUNC!
getMi <- function(){
	#data = getData(dataDir)
	data = getData(dataDir,exc=rep(0,24))

	cat(ncol(data$x)," features: ",colnames(data$x),"\n")

	#theMi = getModelInfoPatientCV(data,getModelFunc=getModelForPCV_SVM)
	theMi = getModelInfoPatientCV(data,getModelFunc=getModelForPCV_MLP)

	cat("Acc (patient-fold cross-validation): ",theMi$acc,"\n\n")

	#cat("Acc (patient-fold cross-validation): ",getModelInfoPatientCV(data, getModelFunc=getModelForPCV_SVM)$acc,"\n\n")
	#cat("Acc (patient-fold cross-validation): ",getModelInfoPatientCV(data, getModelFunc=getModelForPCV_MLP)$acc,"\n\n")

	return(theMi)
}


## End of CLASSIF_TOOLS


## GA_FEATURE_SELECTION

getAttrNames <- function(x,included=TRUE){
	attrList = c(
		"meanpx","meanpy","meanvx","meanvy","meanfx","meanfy",
		"sdevpx","sdevpy","sdevvx","sdevvy","sdevfx","sdevfy",
		"skewpx","skewpy","skewvx","skewvy","skewfx","skewfy",
		"kurtpx","kurtpy","kurtvx","kurtvy","kurtfx","kurtfy"
	)
	return(attrList[which((x==0)==included)])
}

# fitness function for attribute selection
dim=-1
suggs = NULL
attrFitness <- function(x){#attributes with 1 are to be disregarded
	if(dim > -1 && sum(1-x) >= dim){
		return(0)
	}
	attrList = c(
		"meanpx","meanpy","meanvx","meanvy","meanfx","meanfy",
		"sdevpx","sdevpy","sdevvx","sdevvy","sdevfx","sdevfy",
		"skewpx","skewpy","skewvx","skewvy","skewfx","skewfy",
		"kurtpx","kurtpy","kurtvx","kurtvy","kurtfx","kurtfy"
	)

	exc = attrList[which(x==1)]

	data = getData(dir=dataDir,exc=exc)

	#with patient bias (10-fold cross-validation)
	return(-getModel(data)$tot.accuracy)

	#without patient bias (patient-wise cross-validation)
	#return(-getModelInfoPatientCV(data,getModelFunc=getModelForPCV_SVM)$acc)

}

doGA <- function(){
	gaModel = rbga.bin(
		size=24, #size of a chromosome
		popSize=150, #100, #size of a population
		iters=100, #300, #generations #STOPPING CRITERION
		mutationChance = 0.1,
		elitism=TRUE,
		evalFunc=attrFitness,
		verbose=TRUE,
		suggestions=suggs
	)
	cat(summary(gaModel))
	return(gaModel)
}

selectFeaturesWithGA2 <- function(tol,x=NULL,highest=0,prevgm=NULL){
	suggs <<- NULL
	dim <<- -1
	nextgm=NULL
	current = 0

	if(length(x)==24){
		#highest = -attrFitness(x)
		#current = highest
		current = -attrFitness(x)
		if(highest==0){
			highest=current
		}
		dim <<- sum(1-x)
		suggs <<- rbind(suggs,x)
		for(i in which(x==0)){
			thisX = x
			thisX[i]=1
			suggs <<- rbind(suggs,thisX)
		}
	} else {
		nextgm = doGA();
		highest = -attrFitness(nextgm$population[1,])
		current = highest
		dim <<- sum(1-nextgm$population[1,])
	}

	cat("Attempting to fetch less than ",dim,"\n")
	nextgm = doGA()
	current = -attrFitness(nextgm$population[1,])
	#dim <<- sum(1-nextgm$population[1,]) #nao precisa.. testar!

	plot(nextgm)
	if(!(highest-current < tol)){
		cat("FALED: ",current," (highest: ",highest,")\n")
		return(prevgm)
	} else {
		cat("SUCCEEDED: ",current," (highest: ",highest,")\n")
		nextgm=selectFeaturesWithGA(tol,nextgm$population[1,],highest,nextgm)
	}

	
	return(nextgm)
}

selectFeaturesWithGA <- function(tol,x=NULL,highest=0){
	suggs <<- NULL
	dim <<- -1

	if(length(x)==0){
		x = doGA()$population[1,];

		#simulating output:
		#x = c(0,1,0,1,0,1,1,0,1,0,0,0,0,1,1,0,0,0,1,0,0,1,0,1)
	}
	for(i in which(x==0)){
		thisX = x
		thisX[i]=1
		suggs <<- rbind(suggs,thisX)
	}

	current = -attrFitness(x)
	if(highest==0){
		highest=current
		cat("HIGHEST: ",highest,"\n")
	}
	dim <<- sum(1-x)

	cat("Attempting to fetch less than ",dim,"\n")
	nextgm = doGA()$population[1,]
	current = -attrFitness(nextgm)

	if(current > highest){#TESTAR ISTO!
		highest=current
	}

	if(!(highest-current < tol)){
		cat("FALED: ",current," (highest: ",highest,")\n")
		return(x)
	} else {
		cat("SUCCEEDED: ",current," (highest: ",highest,")\n")
		nextgm=selectFeaturesWithGA(tol,nextgm,highest)
	}

	
	return(nextgm)
}

## End of GA_FEATURE_SELECTION




## GRID_SEARCH

doRawGridSearch <- function(){ #raw feature
	data = getData(dataDir,exc=rep(0,24))
	return(gridSearchSVM(data))
}

doGridSearch <- function(){ #reduced features, result of FeatureSelection algorithm
	data = getData(dataDir,exc=getAttrNames(c(0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,0,1,1,0,1,0,1),included=FALSE)) 
	return(gridSearchSVM(data))
}

gridSearchSVM <- function(data){ #configuration of Table 4
	kernel = c("linear","poly")
	coef0 = c(0,1)
	gamma = 2:15/200+0.05
	degree = c(2,3)
	cost = 5:24/15

	bestAcc = 0
	bestParams = list()

	progress = 0
	maxProgress = length(kernel)*length(coef0)*length(gamma)*length(degree)*length(cost)


	for(k in kernel){

		if(k=="poly"){
			for(c in coef0){
				for(g in gamma){
					for(d in degree){
						for(c in cost){
							cat(progress,"th iteration out of ",maxProgress,"\n")
							progress = progress+1


							params = list()
							params$kernel = k
							params$coef0 = c
							params$gamma = g
							params$degree = d
							params$cost = c

							theMi = getModelInfoPatientCV(data,getModelFunc=getModelForPCV_SVM, params=params)
							if(theMi$acc > bestAcc){
								bestAcc = theMi$acc
								params$mi = theMi
								bestParams = params
							}

						}
					}
				}
			}
		} else {
			for(c in cost){
				cat(progress,"th iteration out of ",maxProgress,"\n")
				progress = progress+1

				params = list()
				params$kernel = "linear"
				params$coef0 = 0
				params$gamma = 0
				params$degree = 1
				params$cost = c
				
				theMi = getModelInfoPatientCV(data,getModelFunc=getModelForPCV_SVM, params=params)
				if(theMi$acc > bestAcc){
					bestAcc = theMi$acc
					params$mi = theMi
					bestParams = params
				}
			}
		}

	}
	return(bestParams)
}

gridSearchMLP <- function(data){
	errorFunc = c("ce","sse") # cross-entropy and sum of squared errors
	threshold = 1.3*10^(-4:-1)

	bestAcc = 0
	bestParams = list()

	progress = 0
	maxProgress = length(errorFunc)*length(threshold)

	for(ef in errorFunc){
		for(t in threshold){
			cat(progress,"th iteration out of ",maxProgress,"\n")
			progress = progress+1

			params = list()
			params$err.fct = ef
			params$threshold = t

			theMi = getModelInfoPatientCV(data,getModelFunc=getModelForPCV_MLP, params=params)
			if(theMi$acc > bestAcc){
				bestAcc = theMi$acc
				params$mi = theMi
				bestParams = params
			}
		}
	}
	return(bestParams)
}

## End of GRID_SEARCH


## PLOTTING

plotSVM <- function(dataset){

	aside = dataset$y
	allData = cbind(dataset$x,aside)
	model = svm(aside~.,data=allData,kernel="linear",cost=0.4)
	
	values = rnorm(ncol(dataset$x),0,2000)/100000

	#plot(x=model, data=allData, formula=skewfy ~ skewfx, fill=T, grid=50, slice=list(
	#meanpx=values[1],
	#meanvx=values[2],
	#meanvy=values[3],
	#meanfx=values[4],
	#sdevvy=values[5],
	#sdevfx=values[6],
	#sdevfy=values[7],
	#skewpx=values[8],
	#skewvy=values[9],
	#skewfx=values[10],
	#skewfy=values[11]
	#))
	
	cat("auheoauhe\n")
	plot(model,mode="pca")

}



## End of PLOTTING
