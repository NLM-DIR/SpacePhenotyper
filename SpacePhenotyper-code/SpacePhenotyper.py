import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sksurv.metrics import concordance_index_censored
import os
import numba
import time



def readTableAsDF(path):
	return pd.read_csv(path, index_col=0, header=0, sep='\t' if (path.endswith(".tsv") or path.endswith(".tsv.gz")) else ',')



def predictFromExpressionNumpy(npExpr, npEP, npResults):
	np.dot(npExpr, npEP, out=npResults)



def findFullFilenameTSV(pathFolder, prefix):
	for filename in os.listdir(pathFolder):
		if filename.startswith(prefix) and filename.endswith(".tsv"):
			return filename
	return None



def remove_extensions(filePath):
	extensions = [".zip", ".gz", ".tar", ".csv", ".tsv"]
	for ext in extensions:
		if (filePath.endswith(ext)):
			filePath = filePath[: -len(ext)]
	return filePath



def rot2d(x, y, theta, unit='degree'):
	a =np.vstack([x,y]).T
	if unit == 'degree':
		theta = theta*np.pi/180 # convert to radian
	R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
	ar=a.dot(R.T)
	return ar[:,0], ar[:,1]



def mirrorHorizontally(x, y):
	a = np.vstack([x,y]).T
	R = np.array([[-1, 0], [0, 1]])
	ar = a.dot(R.T)
	return ar[:,0], ar[:,1]



def logAndPlot_BulkExpressionPhenotypeAssociation(dfCorr, path):
	dfCorr.to_csv(path + ".tsv", sep='\t')
	mpl.rcdefaults()
	plt.figure(figsize=(8,6))
	plt.hist(dfCorr, bins=50, edgecolor="black")
	plt.savefig(path + ".pdf", dpi = 75, bbox_inches = 'tight')
	plt.close()



def PlotHeatmap(xloc, yloc, predictions, fig, ax, pointSize=9, minQuant=0.02, maxQuant=0.98):
	mpl.rcdefaults()
	vmin = predictions.quantile(minQuant).item()
	vmax = predictions.quantile(maxQuant).item()
	# vmin = predictions.min()
	# vmax = predictions.max()
	
	colorMap = mcolors.LinearSegmentedColormap.from_list(
		"white_yellow_red",
		[(0, 0, 0), (1, 1, 125/255), (220/255, 0, 0)],	# black, olive, light-yellow, orange, dark red
	)
	s = sns.scatterplot(
		x = xloc,
		y = yloc,
		hue = predictions,
		palette = colorMap,
		hue_norm = mcolors.Normalize(
			vmin = vmin,
			vmax = vmax,
		),
		legend=False,
		# kwargs for matplotlib.axes.Axes.scatter() below
		s = pointSize,	# size of all points
		marker = 'h',	# hexagon-shaped point marker
		edgecolors = None,
		ax = ax,
	)
	ax.collections[0].set_rasterized(True)
	ax._cbar = fig.colorbar(
		mpl.cm.ScalarMappable(
			norm = mcolors.Normalize(vmin, vmax),
			cmap = colorMap,
		),
		ax = ax,
		orientation = 'vertical',
		shrink=0.75,
		pad = 0.02,
	)
	ax._cbar.set_ticks([vmin, vmax])
	ax._cbar.set_ticklabels(["min", "max"])
	ax._cbar.ax.tick_params(labelsize=14)
	ax.set_facecolor((0.92, 0.92, 0.92))	# light gray
	
	s.set(xticklabels = []); s.set(yticklabels = [])
	s.set(xlabel = None); s.set(ylabel = None)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.spines['top'].set_linewidth(0.95)
	ax.spines['bottom'].set_linewidth(0.95)
	ax.spines['left'].set_linewidth(0.95)
	ax.spines['right'].set_linewidth(0.95)
	ax.set_box_aspect(1)
	sns.set(font_scale=1.5, style='white', font='Arial') # Change the font size depending on figure size



def PlotPrediction(xloc, yloc, predictions, pathFigure, pointSize=9):
	mpl.rcdefaults()
	fig, ax = plt.subplots(figsize=(6.5, 5))
	PlotHeatmap(xloc=xloc, yloc=yloc, predictions=predictions, fig=fig, ax=ax, pointSize=pointSize)
	plt.savefig(pathFigure, dpi = 75, bbox_inches = 'tight')
	plt.close()



def makeNames(prefix, n):
	names = []
	for i in range(n):
		names.append(prefix + str(i+1))
	return names



# Minimizes the negative log likelihood by finding the root of the derivative with Newton's algorithm
@numba.njit
def coxUnivariateNewton_assumingSorted(x, time, event, max_iter=25, tol=1e-9):
	n = len(x)
	coeff = 0.0

	for _ in range(max_iter):
		riskScores = coeff * x
		maxRiskSoFar = -np.inf
		sumExp = 0.0
		sumXjExp = 0.0
		sumXjXjExp = 0.0
		derivative = 0.0
		derivative2 = 0.0
		# Assuming that [x, time, event] are sorted (descending) with time as the key
		i = 0
		while i < n:
			j = i
			eventXj = 0.0
			d = 0
			while ((j < n) and (time[j] == time[i])):
				# Changing the exponents for the sake of numerical stabilization (to avoid overflow). The subtracted number in the exponents cancels out when computing the ratios of sums in the first and second derivatives, so the end result is the same.
				# By dividing by the greatest exponent seen so far, we avoid having huge numbers that can cause numberical instability via overflow and/or loss of precision
				if (riskScores[j] > maxRiskSoFar):
					if (maxRiskSoFar != -np.inf):
						diffFactor = np.exp(maxRiskSoFar - riskScores[j])
						# rescale the cumulative sums
						sumExp *= diffFactor
						sumXjExp *= diffFactor
						sumXjXjExp *= diffFactor
					maxRiskSoFar = riskScores[j]
					expRiskScore = 1.0	# exp(riskScores[j] - riskScores[j]) = 1
				else:
					expRiskScore = np.exp(riskScores[j] - maxRiskSoFar)
				
				sumExp += expRiskScore
				sumXjExp += expRiskScore * x[j]
				sumXjXjExp += expRiskScore * x[j] * x[j]
				if event[j] == 1:
					eventXj += x[j]
					d += 1
				j += 1
			# All events with time up to and including the j-th (remember, sorted descending) are in the risk set
			if d > 0:
				# Breslow stuff
				mean = sumXjExp / sumExp
				variance = sumXjXjExp / sumExp - mean * mean
				derivative += eventXj - d * mean
				derivative2 += d * variance
			i = j
		if (derivative2 < tol):
			break
		step = derivative / derivative2
		coeff += step
		if (abs(step) < tol):
			break
	return coeff



@numba.njit
def coxMultivariateNewton_assumingSorted(X, time, event, max_iter=250, tol=1e-9):
	n, m = X.shape
	beta = np.zeros(m)
	grad = np.empty(m)
	hess = np.empty((m, m))
	eventXj = np.empty(m)
	sumXjExp = np.empty(m)
	sumXjXjExp = np.empty((m, m))
	for _ in range(max_iter):
		riskScores = np.dot(X, beta)
		maxRiskSoFar = -np.inf
		sumExp = 0.0
		sumXjExp.fill(0.0)
		sumXjXjExp.fill(0.0)
		grad.fill(0.0)
		hess.fill(0.0)
		# Assuming that [x, time, event] are sorted (descending) with time as the key
		i = 0
		while i < n:
			j = i
			eventXj.fill(0.0)
			d = 0
			while ((j < n) and (time[j] == time[i])):
				# Changing the exponents for the sake of numerical stabilization (to avoid overflow). The subtracted number in the exponents cancels out when computing the ratios of sums in the gradient and hessian, so the end result is the same.
				# By dividing by the greatest exponent seen so far, we avoid having huge numbers that can cause numberical instability via overflow and/or loss of precision
				if (riskScores[j] > maxRiskSoFar):
					if (maxRiskSoFar != -np.inf):
						diffFactor = np.exp(maxRiskSoFar - riskScores[j])
						# rescale the cumulative sums
						sumExp *= diffFactor
						sumXjExp *= diffFactor
						sumXjXjExp *= diffFactor
					maxRiskSoFar = riskScores[j]
					expRiskScore = 1.0	# exp(riskScores[j] - riskScores[j]) = 1
				else:
					expRiskScore = np.exp(riskScores[j] - maxRiskSoFar)
				
				sumExp += expRiskScore
				sumXjExp += expRiskScore * X[j]
				sumXjXjExp += expRiskScore * np.outer(X[j], X[j])
				if (event[j] == 1):
					eventXj += X[j]
					d += 1
				j += 1
			# All events with time up to and including the j-th (remember, sorted descending) are in the risk set
			if (d > 0):
				# Breslow
				means = sumXjExp / sumExp
				covariances = np.outer(means, means) - sumXjXjExp / sumExp
				grad += eventXj - d * means
				hess += d * covariances
			i = j
		if (np.linalg.norm(hess) < tol):
			break
		if (np.linalg.det(hess) < 1e-16):
			# Adding 1.0 to diagonal entries of the hessian to prevent it being a singular matrix. It's fine to be slightly off the exact solution since it's assuming a 2nd degree Taylor approximation of the function anyway.
			step = np.linalg.solve(hess + 1.00*np.eye(m), -grad)
		else:
			step = np.linalg.solve(hess, -grad)
		beta += step
		if (np.linalg.norm(step) < tol):
			break
	return beta



@numba.njit
def getMaxNumPCs(npBulk, npGenePhiCorr, minTotalVariance):
	idxGenesPos = np.flatnonzero(npGenePhiCorr > 0)
	idxGenesNeg = np.flatnonzero(npGenePhiCorr < 0)
	if ((len(idxGenesPos) > 0) and (len(idxGenesNeg) > 0)):
		Up, Sp, VTp = np.linalg.svd(npBulk[:, idxGenesPos], full_matrices = False)
		Un, Sn, VTn = np.linalg.svd(npBulk[:, idxGenesNeg], full_matrices = False)
		cumSum = np.cumsum(np.square(Sp))
		lastIdxP = 0
		while ((lastIdxP+1 < len(Sp)) and (cumSum[lastIdxP+1] / cumSum[-1] <= minTotalVariance)):
			lastIdxP += 1
		cumSum = np.cumsum(np.square(Sn))
		lastIdxN = 0
		while ((lastIdxN+1 < len(Sn)) and (cumSum[lastIdxN+1] / cumSum[-1] <= minTotalVariance)):
			lastIdxN += 1
		return 1 + min(lastIdxP, lastIdxN)
	else:
		U, S, VT = np.linalg.svd(npBulk, full_matrices = False)
		cumSum = np.cumsum(np.square(S))
		lastIdx = 0
		while ((lastIdx+1 < len(S)) and (cumSum[lastIdx+1] / cumSum[-1] <= minTotalVariance)):
			lastIdx += 1
		return 1 + lastIdx



@numba.njit
def getEigenPatientContinuousSplit(npBulk, npPhi, npGenePhiCorr, numTopPCs=1):
	idxGenesPos = np.flatnonzero(npGenePhiCorr > 0)
	idxGenesNeg = np.flatnonzero(npGenePhiCorr < 0)
	vecOnes = np.ones((npBulk.shape[0], 1))
	if ((len(idxGenesPos) > 0) and (len(idxGenesNeg) > 0)):
		Up, Sp, VTp = np.linalg.svd(npBulk[:, idxGenesPos], full_matrices = False)
		Un, Sn, VTn = np.linalg.svd(npBulk[:, idxGenesNeg], full_matrices = False)
		nPCp = min(numTopPCs, Up.shape[1])
		nPCn = min(numTopPCs, Un.shape[1])
		Up = np.ascontiguousarray(Up[:, 0:nPCp])
		Un = np.ascontiguousarray(Un[:, 0:nPCn])
		VpinvSp = np.dot(np.ascontiguousarray(VTp.T[:, 0:nPCp]), np.diag(1.0/Sp[0:nPCp]))
		VninvSn = np.dot(np.ascontiguousarray(VTn.T[:, 0:nPCn]), np.diag(1.0/Sn[0:nPCn]))
		X = np.column_stack((Up, Un))
		regCoeffs = np.linalg.lstsq(np.column_stack((X, vecOnes)), npPhi)[0]
		EG = np.dot(X, regCoeffs[:-1])
		EP = np.empty(npBulk.shape[1])
		EP[idxGenesPos] = np.dot(VpinvSp, regCoeffs[:nPCp])
		EP[idxGenesNeg] = np.dot(VninvSn, regCoeffs[nPCp:-1])
	else:
		U, S, VT = np.linalg.svd(npBulk, full_matrices = False)
		nPC = min(numTopPCs, U.shape[1])
		U = np.ascontiguousarray(U[:, 0:nPC])
		VinvS = np.dot(np.ascontiguousarray(VT.T[:, 0:nPC]), np.diag(1.0/S[0:nPC]))
		X = U
		regCoeffs = np.linalg.lstsq(np.column_stack((X, vecOnes)), npPhi)[0]
		EG = np.dot(X, regCoeffs[:-1])
		EP = np.dot(VinvS, regCoeffs[:-1])
	return EP, EG



@numba.njit
def getEigenPatientCoxSplit_assumingSorted(npBulk, npEvent, npTime, npGenePhiCorr, numTopPCs=1):
	idxGenesPos = np.flatnonzero(npGenePhiCorr > 0)
	idxGenesNeg = np.flatnonzero(npGenePhiCorr < 0)
	if ((len(idxGenesPos) > 0) and (len(idxGenesNeg) > 0)):
		Up, Sp, VTp = np.linalg.svd(npBulk[:, idxGenesPos], full_matrices = False)
		Un, Sn, VTn = np.linalg.svd(npBulk[:, idxGenesNeg], full_matrices = False)
		nPCp = min(numTopPCs, Up.shape[1])
		nPCn = min(numTopPCs, Un.shape[1])
		Up = np.ascontiguousarray(Up[:, 0:nPCp])
		Un = np.ascontiguousarray(Un[:, 0:nPCn])
		VpinvSp = np.dot(np.ascontiguousarray(VTp.T[:, 0:nPCp]), np.diag(1.0/Sp[0:nPCp]))
		VninvSn = np.dot(np.ascontiguousarray(VTn.T[:, 0:nPCn]), np.diag(1.0/Sn[0:nPCn]))
		X = np.column_stack((Up, Un))
		regCoeffs = coxMultivariateNewton_assumingSorted(X, time=npTime, event=npEvent)
		EG = np.dot(X, regCoeffs)
		EP = np.empty(npBulk.shape[1])
		EP[idxGenesPos] = np.dot(VpinvSp, regCoeffs[:nPCp])
		EP[idxGenesNeg] = np.dot(VninvSn, regCoeffs[nPCp:])
	else:
		U, S, VT = np.linalg.svd(npBulk, full_matrices = False)
		nPC = min(numTopPCs, U.shape[1])
		U = np.ascontiguousarray(U[:, 0:nPC])
		VinvS = np.dot(np.ascontiguousarray(VT.T[:, 0:nPC]), np.diag(1.0/S[0:nPC]))
		regCoeffs = coxMultivariateNewton_assumingSorted(U, time=npTime, event=npEvent)
		EG = np.dot(U, regCoeffs)
		EP = np.dot(VinvS, regCoeffs)
	return EP, EG



@numba.njit
def getEigenPatientCoxSplit(npBulk, npEvent, npTime, npGenePhiCorr, numTopPCs=1):
	time = npTime.copy()
	event = npEvent.copy()
	bulk = npBulk.copy()
	order = np.argsort(-time)
	event = event[order]
	time = time[order]
	bulk = bulk[order]
	EP, EG = getEigenPatientCoxSplit_assumingSorted(npBulk=bulk, npEvent=event, npTime=time, npGenePhiCorr=npGenePhiCorr, numTopPCs=numTopPCs)
	originalOrder = np.argsort(order)
	return EP, EG[originalOrder]



@numba.njit
def getPhenotypeAssociationContinuous(npBulk, npPhi):
	meanPhi = npPhi.sum() / len(npPhi)
	stdPhi = np.sqrt(np.sum((npPhi - meanPhi)**2))
	meanGenes = np.sum(npBulk, axis=0) / npBulk.shape[0]
	stdGenes = np.sqrt(np.sum((npBulk - meanGenes)**2, axis=0))
	return np.dot(npPhi - meanPhi, npBulk - meanGenes) / (stdPhi * stdGenes)



@numba.njit
def getPhenotypeAssociationHazard_assumingSorted(npBulk, npTime, npEvent):
	results = np.empty(npBulk.shape[1], dtype=np.float64)
	for gIdx in range(npBulk.shape[1]):
		beta = coxUnivariateNewton_assumingSorted(
			x=npBulk[:, gIdx] / np.linalg.norm(npBulk[:, gIdx]),
			time=npTime,
			event=npEvent
		)
		results[gIdx] = beta
	return results



@numba.njit(parallel=True)
def getPhenotypeAssociationHazard(npBulk, npTime, npEvent):
	bulk = npBulk.copy()
	time = npTime.copy()
	event = npEvent.copy()
	order = np.argsort(-time)
	bulk = bulk[order]
	time = time[order]
	event = event[order]
	return getPhenotypeAssociationHazard_assumingSorted(npBulk=bulk, npTime=time, npEvent=event)



def scorePredictionsHazard(dfPredictions, dfPhi):
	dfPredictions = dfPredictions.loc[dfPhi.index]
	
	result = concordance_index_censored(event_indicator=dfPhi['vital_status'].astype(bool), event_time=dfPhi['survival_time'], estimate=dfPredictions['SpacePhenotyper'])
	return result[0]



@numba.njit
def scorePredictionsContinuous(npPredicted, npGroundTruth):
	# Pearson correlation
	return np.corrcoef(npPredicted, npGroundTruth)[0, 1].item()



@numba.njit(parallel=True)
def foldCV_numba_continuous(npBulk, npPhi, percentile=0.05, maxNumPCs=50, CVNumFolds=10, CVNumShuffles=1, seed=1337):
	# determine the maximum number of PCs
	genePhiCorr = getPhenotypeAssociationContinuous(npBulk, npPhi)
	idx = np.argsort(np.abs(genePhiCorr))
	numGenesUsed = min(npBulk.shape[1], int(np.ceil(npBulk.shape[1]*percentile)))
	usedGeneIndices = idx[-numGenesUsed:]
	maxNumPCs = min(maxNumPCs, getMaxNumPCs(npBulk=npBulk[:, usedGeneIndices], npGenePhiCorr=genePhiCorr[usedGeneIndices], minTotalVariance=0.99))
	
	CVNumShuffles = int(CVNumShuffles)
	maxNumPCs = int(maxNumPCs)
	CVNumFolds = int(CVNumFolds)
	numJobs = maxNumPCs * CVNumFolds * CVNumShuffles
	CVResults = np.full((maxNumPCs, CVNumFolds*CVNumShuffles), np.nan, dtype=np.float64)
	
	numSamples = npBulk.shape[0]
	bucketSize = numSamples // CVNumFolds
	dataBulk = [npBulk] * CVNumShuffles
	dataPhi = [npPhi] * CVNumShuffles
	order = np.arange(numSamples)
	np.random.seed(seed)
	for i in range(CVNumShuffles):
		np.random.shuffle(order)
		dataBulk[i] = npBulk[order]
		dataPhi[i] = npPhi[order]
	
	for jIdx in numba.prange(numJobs):
		numPCs = 1 + jIdx // (CVNumFolds*CVNumShuffles)
		rest = jIdx % (CVNumFolds*CVNumShuffles)
		foldIdx = rest // CVNumShuffles
		shuffleIdx = rest % CVNumShuffles
		leftEndpoint = foldIdx * bucketSize
		rightEndpoint = (foldIdx+1)*bucketSize if (foldIdx < CVNumFolds-1) else numSamples
		maskTrain = np.ones(numSamples, dtype=np.bool_)
		maskTrain[leftEndpoint:rightEndpoint] = False
		maskTest = ~maskTrain
		bulkTrain = (dataBulk[shuffleIdx])[maskTrain]
		phiTrain = (dataPhi[shuffleIdx])[maskTrain]
		bulkTest = (dataBulk[shuffleIdx])[maskTest]
		phiTest = (dataPhi[shuffleIdx])[maskTest]
		# Remove genes with no expression left after cutting out test samples
		colSums = np.sum(bulkTrain, axis=0)
		geneMask = (colSums > 0)
		bulkTrain = bulkTrain[:, geneMask]
		### calculate gene-to-phenotype association
		genePhiCorr = getPhenotypeAssociationContinuous(bulkTrain, phiTrain)
		### select top percentile% genes by absolute value of association
		idx = np.argsort(np.abs(genePhiCorr))
		numGenesUsed = min(bulkTrain.shape[1], int(np.ceil(bulkTrain.shape[1]*percentile)))
		usedGeneIndices = idx[-numGenesUsed:]
		bulkTrain = bulkTrain[:, usedGeneIndices]
		bulkTest = bulkTest[:, usedGeneIndices]
		genePhiCorr = genePhiCorr[usedGeneIndices]
		EP, _ = getEigenPatientContinuousSplit(npBulk=bulkTrain, npPhi=phiTrain, npGenePhiCorr=genePhiCorr, numTopPCs=numPCs)
		### predict on the test samples
		CVResults[numPCs-1, foldIdx*CVNumShuffles + shuffleIdx] = scorePredictionsContinuous(npPredicted=np.dot(bulkTest, EP), npGroundTruth=phiTest)
	return CVResults



@numba.njit(parallel=True)
def foldCV_numba_hazard(npBulk, npEvent, npTime, percentile=0.05, maxNumPCs=50, CVNumFolds=10, CVNumShuffles=1, seed=1337):
	# determine the maximum number of PCs
	genePhiCorr = getPhenotypeAssociationHazard(npBulk=npBulk, npTime=npTime, npEvent=npEvent)
	idx = np.argsort(np.abs(genePhiCorr))
	numGenesUsed = min(npBulk.shape[1], int(np.ceil(npBulk.shape[1]*percentile)))
	usedGeneIndices = idx[-numGenesUsed:]
	maxNumPCs = min(maxNumPCs, getMaxNumPCs(npBulk=npBulk[:, usedGeneIndices], npGenePhiCorr=genePhiCorr[usedGeneIndices], minTotalVariance=0.99))
	
	CVNumShuffles = int(CVNumShuffles)
	maxNumPCs = int(maxNumPCs)
	CVNumFolds = int(CVNumFolds)
	sizeFoldsShuffles = int(CVNumFolds*CVNumShuffles)
	numJobs = maxNumPCs * CVNumFolds * CVNumShuffles
	CVResults = np.full((maxNumPCs, CVNumFolds*CVNumShuffles), np.nan, dtype=np.float64)
	
	numSamples = npBulk.shape[0]
	idxDead = np.flatnonzero(npEvent == 1)
	idxAlive = np.flatnonzero(npEvent == 0)
	bucketSizeDead = idxDead.shape[0] // CVNumFolds
	bucketSizeAlive = idxAlive.shape[0] // CVNumFolds
	
	dataBulk = [npBulk] * CVNumShuffles
	dataEvent = [npEvent] * CVNumShuffles
	dataTime = [npTime] * CVNumShuffles
	dataIdxDead = [idxDead] * CVNumShuffles
	dataIdxAlive = [idxAlive] * CVNumShuffles
	order = np.arange(numSamples)
	np.random.seed(seed)
	for i in range(CVNumShuffles):
		np.random.shuffle(order)
		dataBulk[i] = npBulk[order]
		dataEvent[i] = npEvent[order]
		dataTime[i] = npTime[order]
		dataIdxDead[i] = np.flatnonzero(dataEvent[i] == 1)
		dataIdxAlive[i] = np.flatnonzero(dataEvent[i] == 0)
	
	for jIdx in numba.prange(numJobs):
		numPCs = 1 + (jIdx // sizeFoldsShuffles)
		rest = int(jIdx % sizeFoldsShuffles)
		foldIdx = int(rest // CVNumShuffles)
		shuffleIdx = int(rest % CVNumShuffles)
		
		leftEndpointDead = foldIdx * bucketSizeDead
		rightEndpointDead = (foldIdx+1)*bucketSizeDead if (foldIdx < CVNumFolds-1) else idxDead.shape[0]
		leftEndpointAlive = foldIdx * bucketSizeAlive
		rightEndpointAlive = (foldIdx+1)*bucketSizeAlive if (foldIdx < CVNumFolds-1) else idxAlive.shape[0]
		idxTest = np.concatenate(((dataIdxDead[shuffleIdx])[leftEndpointDead:rightEndpointDead], (dataIdxAlive[shuffleIdx])[leftEndpointAlive:rightEndpointAlive]))
		maskTrain = np.ones(numSamples, dtype=np.bool_)
		maskTrain[idxTest] = False
		maskTest = ~maskTrain
		bulkTrain = (dataBulk[shuffleIdx])[maskTrain]
		eventTrain = (dataEvent[shuffleIdx])[maskTrain]
		timeTrain = (dataTime[shuffleIdx])[maskTrain]
		bulkTest = (dataBulk[shuffleIdx])[maskTest]
		eventTest = (dataEvent[shuffleIdx])[maskTest]
		timeTest = (dataTime[shuffleIdx])[maskTest]
		# sort by reverse time to make the cox regression below correct
		order = np.argsort(-timeTrain)
		bulkTrain = bulkTrain[order]
		eventTrain = eventTrain[order]
		timeTrain = timeTrain[order]
		# Remove genes with no expression left after cutting out test samples
		colSums = np.sum(bulkTrain, axis=0)
		geneMask = (colSums > 0)
		bulkTrain = bulkTrain[:, geneMask]
		### calculate gene-to-phenotype association
		genePhiCorr = getPhenotypeAssociationHazard_assumingSorted(npBulk=bulkTrain, npTime=timeTrain, npEvent=eventTrain)
		### select top percentile% genes by absolute value of association
		idx = np.argsort(np.abs(genePhiCorr))
		numGenesUsed = min(bulkTrain.shape[1], int(np.ceil(bulkTrain.shape[1]*percentile)))
		usedGeneIndices = idx[-numGenesUsed:]
		bulkTrain = bulkTrain[:, usedGeneIndices]
		bulkTest = bulkTest[:, usedGeneIndices]
		genePhiCorr = genePhiCorr[usedGeneIndices]
		# ### Calculate the Eigen-Patient
		EP, _ = getEigenPatientCoxSplit_assumingSorted(npBulk=bulkTrain, npEvent=eventTrain, npTime=timeTrain, npGenePhiCorr=genePhiCorr, numTopPCs=numPCs)
		### predict on the test samples
		# First sort the test samples by reverse time to make the concordance index calculation below correct
		order = np.argsort(-timeTest)
		bulkTest = bulkTest[order]
		eventTest = eventTest[order]
		predTest = np.dot(bulkTest, EP)
		## full concordance index (ignoring ties)
		numComparable = 0
		numCorrectOrder = 0
		for i in np.flatnonzero(eventTest == 1):
			numComparable += i
			numCorrectOrder += np.sum(predTest[:i] < predTest[i])
		if (numComparable > 0):
			CVResults[numPCs-1, foldIdx*CVNumShuffles + shuffleIdx] = numCorrectOrder / numComparable
	return CVResults



















def SpacePhenotyper_Train(dfBulk, dfPhi, phenotype, percentile=0.05, numPCs=5):
	# STEP 1: Find gene associations with the phenotype
	genePhiCorr = (
		getPhenotypeAssociationHazard(npBulk=dfBulk.to_numpy(), npTime=dfPhi['survival_time'].to_numpy(), npEvent=dfPhi['vital_status'].to_numpy())
		if (phenotype == "hazard")
		else getPhenotypeAssociationContinuous(dfBulk.to_numpy(), dfPhi.to_numpy().ravel())
	)
	# print(1)
	dfCorr = pd.DataFrame(genePhiCorr.T, columns=['assocWithPhenotype'], index=dfBulk.columns)
	dfCorr.sort_values(by=dfCorr.columns[0], key=lambda x: x.abs(), ascending=False, inplace=True)
	
	numGenesUsed = min(dfBulk.shape[1], int(np.ceil(dfBulk.shape[1]*percentile)))
	idx = np.argsort(np.abs(genePhiCorr))
	usedGeneIndices = idx[-numGenesUsed:]
	npBulkFiltered = (dfBulk.iloc[:, usedGeneIndices]).to_numpy()
	
	# print(2)
	# STEP 3: Construct the Eigen-Patient based on the chosen gene subset
	EP, EG = (
		getEigenPatientCoxSplit(npBulk=npBulkFiltered, npEvent=dfPhi['vital_status'].to_numpy(), npTime=dfPhi['survival_time'].to_numpy(), npGenePhiCorr=genePhiCorr[usedGeneIndices], numTopPCs=numPCs)
		if (phenotype == "hazard")
		else getEigenPatientContinuousSplit(npBulk=npBulkFiltered, npPhi=dfPhi.iloc[:, 0].to_numpy(), npGenePhiCorr=genePhiCorr[usedGeneIndices], numTopPCs=numPCs)
	)
	# print(3)
	return EP, EG, usedGeneIndices, dfCorr



def SpacePhenotyper_TrainAndPredict(dfBulk, dfPhi, dfSRT, dfSRTCoords, phenotype="continuous", CV=None, CVNumShuffles=1, percentile=0.05, numPCs=10, maxNumPCs=0):
	assert (dfPhi.index.is_unique), "ERROR: Phenotype data cannot have duplicate row (bulk sample) names"
	assert (dfBulk.index.is_unique), "ERROR: Bulk expression data cannot have duplicate row (bulk sample) names"
	assert (dfBulk.columns.is_unique), "ERROR: Bulk expression data cannot have duplicate column (gene) names"
	assert (dfSRT.index.is_unique), "ERROR: SRT expression data cannot have duplicate row (spot) names"
	assert (dfSRT.columns.is_unique), "ERROR: SRT expression data cannot have duplicate column (gene) names"
	assert (dfSRTCoords.index.is_unique), "ERROR: SRT coordinate data cannot have duplicate row (spot) names"
	
	# Remove samples whose expression is ==0 across all genes
	dfBulk = dfBulk[dfBulk.sum(axis=1) > 0]
	dfSRT = dfSRT[dfSRT.sum(axis=1) > 0]
	# Remove genes whose expression is ==0 across all samples
	dfBulk = dfBulk.loc[:, dfBulk.sum(axis=0) > 0]
	dfSRT = dfSRT.loc[:, dfSRT.sum(axis=0) > 0]
	# Remove genes that aren't expressed in at least 5% SRT spots
	dfSRT = dfSRT.loc[:, (dfSRT == 0).mean() <= 0.95]
	# Intersect genes between bulk and SRT
	commonGenes = dfSRT.columns.intersection(dfBulk.columns)
	dfSRT = dfSRT[commonGenes]
	dfBulk = dfBulk[commonGenes]
	# Intersect samples between bulk expression and bulk phenotype
	commonBulkSamples = dfPhi.index.intersection(dfBulk.index)
	dfPhi = dfPhi.loc[commonBulkSamples]
	dfBulk = dfBulk.loc[commonBulkSamples]
	print("\nbulk shape: ", dfBulk.shape)
	print("\nphi shape: ", dfPhi.shape)
	
	dfCVResults = None
	if (CV > 1):	# n-fold cross validation
		npBulk = np.ascontiguousarray(dfBulk.to_numpy())
		startTime = time.perf_counter()
		if (phenotype == "hazard"):
			npEvent = np.ascontiguousarray(dfPhi['vital_status'].to_numpy())
			npTime = np.ascontiguousarray(dfPhi['survival_time'].to_numpy())
			CVResults = foldCV_numba_hazard(npBulk, npTime=npTime, npEvent=npEvent, percentile=percentile, maxNumPCs=maxNumPCs, CVNumFolds=CV, CVNumShuffles=CVNumShuffles)
		else:
			npPhi = np.ascontiguousarray(dfPhi.to_numpy().ravel())
			CVResults = foldCV_numba_continuous(npBulk, npPhi, percentile=percentile, maxNumPCs=maxNumPCs, CVNumFolds=CV, CVNumShuffles=CVNumShuffles)
		maxNumPCs = CVResults.shape[0]
		endTime = time.perf_counter()
		print(f"{CV}-fold CV for all cases of up to {CVResults.shape[0]} PCs completed in {endTime - startTime:.3f} seconds")
		dfCVResults = pd.DataFrame(CVResults, index=makeNames("PCs_", CVResults.shape[0]))
		npMeans = np.mean(CVResults, axis=1)
		npStds = np.std(CVResults, axis=1)
		dfCVResults['means'] = npMeans
		dfCVResults['stds'] = npStds
		numPCs = 1 + dfCVResults['means'].to_numpy().argmax()
		print(f"{numPCs} selected as the optimal number of PCs per gene set by CV")
		dfCVResults.sort_values(by="means", ascending=False, inplace=True)
	# else (CV <= 1) no cross validation
	
	npEP, npEG, usedGenesIndices, dfCorr = SpacePhenotyper_Train(dfBulk, dfPhi, phenotype=phenotype, percentile=percentile, numPCs=numPCs)
	
	# Filter out unused genes
	dfSRT = dfSRT.iloc[:, usedGenesIndices]
	dfBulk = dfBulk.iloc[:, usedGenesIndices]
	# Make sure that after filtering out some genes, only SRT samples with >0 remaining gene expression are used
	dfSRT = dfSRT[dfSRT.sum(axis=1) > 0]
	# Intersect samples between SRT expression and coords
	commonSRTSamples = dfSRT.index.intersection(dfSRTCoords.index)
	dfSRT = dfSRT.loc[commonSRTSamples]
	dfSRTCoords = dfSRTCoords.loc[commonSRTSamples]
	
	# Predict phenotype on Bulk
	predictionsOnBulk = np.zeros(dfBulk.shape[0])
	predictFromExpressionNumpy(dfBulk.to_numpy(copy=False), npEP, predictionsOnBulk)
	dfPredictionsOnBulk = pd.DataFrame(predictionsOnBulk.T, index=dfBulk.index, columns = ['SpacePhenotyper'])
	assert (np.linalg.norm(predictionsOnBulk - npEG) < 1e-4), "ERROR: The Eigen-Gene and predictions computed via the Eigen-Patient have a difference norm >= 1e-4"
	
	# Predict phenotype on SRT
	predictions = np.zeros(dfSRT.shape[0])
	predictFromExpressionNumpy(dfSRT.to_numpy(copy=False), npEP, predictions)
	dfPredictions = pd.DataFrame(predictions.T, columns=['SpacePhenotyper'], index=dfSRT.index)
	dfPredictions = pd.concat([dfSRTCoords, dfPredictions], axis=1)
	
	dfVariablePCPredictions = []
	for pcNum in range(1, maxNumPCs+1):
		npGenePhiCorr = dfCorr.loc[dfBulk.columns].T.to_numpy()
		EPtemp, EGtemp = (
			getEigenPatientCoxSplit(dfBulk.to_numpy(), npEvent=dfPhi['vital_status'].to_numpy(), npTime=dfPhi['survival_time'].to_numpy(), npGenePhiCorr=npGenePhiCorr, numTopPCs=pcNum)
			if (phenotype=="hazard")
			else getEigenPatientContinuousSplit(npBulk=dfBulk.to_numpy(), npPhi=dfPhi.iloc[:, 0].to_numpy(), npGenePhiCorr=npGenePhiCorr, numTopPCs=pcNum)
		)
		predictionsSRT = np.zeros(dfSRT.shape[0])
		predictFromExpressionNumpy(npExpr=dfSRT.to_numpy(copy=False), npEP=EPtemp, npResults=predictionsSRT)
		dfPredictionsSRTTemp = pd.DataFrame(predictionsSRT.T, columns=['SpacePhenotyper'], index=dfSRT.index)
		bulkScore = -np.inf
		if (phenotype == "hazard"):
			dfPredictionsBulkTemp = pd.DataFrame(EGtemp.T, columns=['SpacePhenotyper'], index=dfBulk.index)
			bulkScore = scorePredictionsHazard(dfPredictions=dfPredictionsBulkTemp, dfPhi=dfPhi)
		else:
			bulkScore = scorePredictionsContinuous(npPredicted=EGtemp, 
			npGroundTruth=dfPhi.T.to_numpy().ravel())
		dfVariablePCPredictions.append([pd.concat([dfSRTCoords, dfPredictionsSRTTemp], axis=1), bulkScore])
	
	return dfPredictions, npEP, usedGenesIndices, dfCorr, dfPredictionsOnBulk, dfCVResults, dfVariablePCPredictions, dfBulk.index



def SpacePhenotyper_IOWrapper(inputData, outputRoot, phenotype="continuous", phenoNameTP=None, geneFilter=None, pointSize=14, CVNumFolds=0, CVNumShuffles=1, percentile=0.05, numPCs=10, maxNumPCs=0):
	for pathBulkExp, pathBulkPheno, pathSRTExp, pathSRTCoords in inputData:
		fileBulkExp = remove_extensions(os.path.basename(pathBulkExp))
		fileBulkPheno = remove_extensions(os.path.basename(pathBulkPheno))
		fileSRTExp = remove_extensions(os.path.basename(pathSRTExp))
		print("")
		print(fileBulkExp)
		print(fileBulkPheno)
		print(fileSRTExp, flush=True)
		
		outputFolder = outputRoot + fileBulkPheno + "/"
		if not os.path.exists(outputFolder):
			os.makedirs(outputFolder)
		
		pathDataset = fileSRTExp + "_" + fileBulkExp
		pathOutGenePhenotypeAssociationNoExt = outputFolder + "GenePhenotypeAssociation_" + pathDataset
		pathOutEigenPatientNoExt = outputFolder + "EigenPatient_" + pathDataset
		pathOutBulkPredictionsNoExt = outputFolder + "BulkPredictions_" + pathDataset
		pathOutPredictionsNoExt = outputFolder + "Predictions_" + pathDataset
		pathOutCVNoExt = outputFolder + "CV_" + pathDataset
		
		# Rows are samples, columns are genes or metadata
		dfBulk = readTableAsDF(pathBulkExp).astype(np.float64)
		dfPhi = readTableAsDF(pathBulkPheno)
		dfSRT = readTableAsDF(pathSRTExp)
		dfSRTCoords = readTableAsDF(pathSRTCoords)
		dfSRTCoords = dfSRTCoords[['array_row', 'array_col']]
		if (fileSRTExp.startswith("BRCA-")):
			if (fileSRTExp.startswith("BRCA-2.0.0")):
				dfSRTCoords['array_row'], dfSRTCoords['array_col'] = mirrorHorizontally(dfSRTCoords['array_row'], dfSRTCoords['array_col'])
			else:
				dfSRTCoords['array_row'], dfSRTCoords['array_col'] = rot2d(dfSRTCoords['array_row'], dfSRTCoords['array_col'], -90)
		elif (fileSRTExp.startswith("COAD-")):
			dfSRTCoords['array_row'], dfSRTCoords['array_col'] = rot2d(dfSRTCoords['array_row'], dfSRTCoords['array_col'], 90)
			dfSRTCoords['array_row'], dfSRTCoords['array_col'] = mirrorHorizontally(dfSRTCoords['array_row'], dfSRTCoords['array_col'])
		
		print(f"{dfSRT.shape[1]} genes in SRT data")
		print(f"{dfBulk.shape[1]} genes in bulk data")
		commonGenes = (dfSRT.columns).intersection(dfBulk.columns)
		print(f"{len(commonGenes)} genes in common")
		if (geneFilter is not None):
			commonGenes = commonGenes.intersection(geneFilter.index)
			print(f"{len(commonGenes)} genes in common with the given gene filter list")
		dfBulk = dfBulk[commonGenes]
		dfSRT = dfSRT[commonGenes]
		
		if (phenotype == "hazard"):
			# relevantColumns=['vital_status','days_to_death','days_to_last_follow_up', 'age_at_diagnosis', 'treatment_or_therapy']
			relevantColumns=['vital_status','days_to_death','days_to_last_follow_up']
			CoxData = dfPhi[relevantColumns].copy()
			CoxData['days_to_death'] = CoxData['days_to_death'].replace(np.nan, 0)
			CoxData['days_to_last_follow_up'] = CoxData['days_to_last_follow_up'].replace(np.nan, 0)
			CoxData['vital_status'] = CoxData['vital_status'].apply(lambda x: 1 if x == 'Dead' else 0 )
			CoxData['survival_time'] = CoxData['days_to_death'] + CoxData['days_to_last_follow_up']
			CoxData = CoxData[['vital_status','survival_time']].dropna()
			dfPhi = CoxData.copy()
			print("Dead patients:", dfPhi[dfPhi['vital_status']==1].shape[0])
			print("Censored patients:", dfPhi[dfPhi['vital_status']==0].shape[0])
		else:
			dfPhi = dfPhi.astype(np.float64)
		
		if (("tumorpurity" not in fileBulkPheno.lower()) and (phenoNameTP is not None)):
			# pathTumorPurityPredictions = outputRoot + phenoNameTP + "/Predictions_" + pathDataset + ".tsv"
			pathTumorPurityPredictions = outputRoot + phenoNameTP + "/" + findFullFilenameTSV(pathFolder = outputRoot + phenoNameTP, prefix = "Predictions_" + fileSRTExp)
			# SRTTumorPurity = pd.read_csv(pathTumorPurityPredictions, sep='\t', index_col=0, header=0)
			SRTTumorPurity = readTableAsDF(pathTumorPurityPredictions)
			minTP = SRTTumorPurity['SpacePhenotyper'].min()
			maxTP = SRTTumorPurity['SpacePhenotyper'].max()
			cutoffPoint = minTP + (maxTP - minTP)*0.6
			index = SRTTumorPurity[SRTTumorPurity['SpacePhenotyper'] >= cutoffPoint].index
			dfSRT = dfSRT.loc[index]
			dfSRTCoords = dfSRTCoords.loc[index]
		
		dfPredictions, npEP, usedGenesIndices, dfCorr, dfPredictionsOnBulk, dfCVResults, dfVariablePCPredictions, newSampleIndex = SpacePhenotyper_TrainAndPredict(dfBulk, dfPhi, dfSRT, dfSRTCoords, phenotype=phenotype, CV=CVNumFolds, CVNumShuffles=CVNumShuffles, percentile=percentile, numPCs=numPCs, maxNumPCs=maxNumPCs)
		
		dfBulk = dfBulk.loc[newSampleIndex]
		dfPhi = dfPhi.loc[newSampleIndex]
		
		# bugged; these indices do not correspond to indices of genes actually used to construct the eigen-patient
		dfEP = pd.DataFrame(npEP.T, index=dfBulk.columns[usedGenesIndices], columns=['Eigen-Patient'])
		
		print(f"\n", dfCorr)
		print(f"\nBulk shape: ", dfBulk.shape)
		print(f"\nPhi shape: ", dfPhi.shape)
		print(f"\n", dfEP, flush=True)
		
		logAndPlot_BulkExpressionPhenotypeAssociation(dfCorr, pathOutGenePhenotypeAssociationNoExt)
		# logAndPlot_gridResults(dfGeneSubsetScores, pathOutGridSearchNoExt, phenotype=phenotype)
		dfEP.to_csv(pathOutEigenPatientNoExt + '.tsv', sep='\t')
		dfPredictionsOnBulk.to_csv(pathOutBulkPredictionsNoExt + '.tsv', sep='\t')
		if (dfCVResults is not None):
			# print(f"\n", dfCVResults, flush=True)
			dfCVResults.to_csv(pathOutCVNoExt + '.tsv', sep='\t')
		dfPredictions.to_csv(pathOutPredictionsNoExt + '.tsv', sep='\t')
		
		
		
		if (len(dfVariablePCPredictions) > 0):
			fig, axes = plt.subplots(
				nrows=len(dfVariablePCPredictions),
				ncols=1,
				figsize=(6, len(dfVariablePCPredictions)*4),
				sharex=False,
				sharey=False,
				# constrained_layout=True,
				squeeze=False,
			)
			plt.subplots_adjust(wspace=0.02, hspace=0.05)
			for i in range(len(dfVariablePCPredictions)):
				pred, bulkScore = dfVariablePCPredictions[i]
				# print(pred)
				PlotHeatmap(
					xloc = pred['array_row'],
					yloc = pred['array_col'],
					predictions = pred['SpacePhenotyper'],
					fig = fig,
					ax = axes[i, 0],
					pointSize = pointSize,
					minQuant = 0.02,
					maxQuant = 0.98,
				)
				infoTxt = str(i+1) + " PCs per gene set" + "\n\n"
				infoTxt += "Full bulk score: " + str(round(bulkScore, 3)) + "\n"
				if ((dfCVResults is not None) and (("PCs_" + str(i+1)) in dfCVResults.index)):
					mean = dfCVResults['means'].loc["PCs_" + str(i+1)]
					sd = dfCVResults['stds'].loc["PCs_" + str(i+1)]
					infoTxt += "Bulk CV mean: " + str(round(mean, 3)) + "\n"
					infoTxt += "Bulk CV SD: " + str(round(sd, 3))
				axes[i, 0].text(
					-0.1, 0.5,
					infoTxt,
					transform=axes[i, 0].transAxes,
					ha="right",
					va="center",
					fontsize=14,
				)
			if (dfCVResults is not None):
				maxMean = dfCVResults['means'].max()
				# print(maxMean)
				for i in range(len(dfVariablePCPredictions)):
					if ((("PCs_" + str(i+1)) in dfCVResults.index) and (abs(dfCVResults['means'].loc["PCs_" + str(i+1)].item() - maxMean) < 0.0001)):
						# print(i, dfCVResults['means'].loc["PCs_" + str(i+1)].item())
						axes[i, 0].plot(
							[-0.8, -0.05],
							[0.3, 0.3],
							transform=axes[i, 0].transAxes,
							color="red",
							linestyle="-",
							linewidth=2.0,
							clip_on=False,
						)
						axes[i, 0].plot(
							[-0.8, -0.05],
							[0.7, 0.7],
							transform=axes[i, 0].transAxes,
							color="red",
							linestyle="-",
							linewidth=2.0,
							clip_on=False,
						)
						axes[i, 0].plot(
							[-0.8, -0.8],
							[0.3, 0.7],
							transform=axes[i, 0].transAxes,
							color="red",
							linestyle="-",
							linewidth=2.0,
							clip_on=False,
						)
						axes[i, 0].plot(
							[-0.05, -0.05],
							[0.3, 0.7],
							transform=axes[i, 0].transAxes,
							color="red",
							linestyle="-",
							linewidth=2.0,
							clip_on=False,
						)
			axes[0, 0].text(
				0.5, 1.2,
				fileBulkPheno + "\n" + fileSRTExp,
				transform=axes[0, 0].transAxes,
				ha="center",
				va="center",
				fontsize=20,
			)
			plt.savefig(pathOutPredictionsNoExt + "_increasingPCs.pdf", dpi = 75, bbox_inches = 'tight')
			plt.close()
		
		print(pathOutPredictionsNoExt, flush=True)
		PlotPrediction(dfPredictions['array_row'], dfPredictions['array_col'], dfPredictions['SpacePhenotyper'], pathOutPredictionsNoExt + ".pdf", pointSize)
		
