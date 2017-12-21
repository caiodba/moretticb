source("classif.r")

tsDataset = "hemiparesis.csv"

getFinalData <- function(){
	tsCsv = read.csv(paste(c(dataDir,tsDataset),collapse="/"))
	data=list()
	data$x=tsCsv[,getAttrNames(c(0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,0,1,1,0,1,0,1))]
	data$y=tsCsv$aside
	data$timeoffset=tsCsv$timeoffset
	data$pat=tsCsv$patient
	data$device=tsCsv$device

	return(data)

}

getFinalModel <- function(data=getFinalData()){ 
	model = svm(x=data$x,y=data$y,kernel="linear",cost=0.866667,cross=0)
	return(model)
}

getCyclesByPatient <- function(data=getFinalData(),patient){
	indexes = which(data$pat==patient)
	patData=list()
	patData$x = data$x[indexes,]
	patData$timeoffset = data$timeoffset[indexes]
	
	orderIndexes = order(patData$timeoffset)
	patData$x = patData$x[orderIndexes,]
	patData$timeoffset = patData$timeoffset[orderIndexes]
	return(patData)
}

getTsByPatient <- function(patData,model){
	ts=attr(predict(model,patData$x,decision.values=T),"decision.values")
	return(abs(ts))
}

getLinearCoefFromTs <- function(ts){
	x = 1:length(ts)
	lin = lm(ts~x)
	return(lin$coefficients[2])
}

getAllCoeffs <- function(){
	data = getFinalData()
	patients = unique(data$pat)
	model = getFinalModel(data)
	
	coeffData = list()
	coeffData$coeffs = c()
	coeffData$patients = c()
	coeffData$ncycles = c()

	for(p in patients){
		cbp = getCyclesByPatient(data,p)
		if(nrow(cbp$x)>10){
			coeffData$coeffs = c(coeffData$coeffs,getLinearCoefFromTs(getTsByPatient(cbp,model)))
			coeffData$patients = c(coeffData$patients,p)
			coeffData$ncycles = c(coeffData$ncycles,nrow(cbp$x))
		}
	}

	return(coeffData)
}

plotPatientTs <- function(patient){
	data = getFinalData()
	model = getFinalModel(data)
	cyc=getCyclesByPatient(data,patient)
	ts=getTsByPatient(cyc,model)
	
	plot(ts)
	lines(ts)
}
