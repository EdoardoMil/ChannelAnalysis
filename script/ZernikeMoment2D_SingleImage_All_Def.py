import os
import numpy as np
import pandas as pd
import math
import mahotas
import cv2
from PIL import Image
import glob
import os

order = 20

ListJpeg = glob.glob('../results_temp_singleIm/*png')
print(ListJpeg)

conta = 1

ZernAll = []
for i in range(0,121):
    ZernAll.append(0)
    
for filepath in glob.iglob('../results_temp_singleIm/*.png'):

	print(filepath)
	name = int(''.join(filter(str.isdigit, filepath)))
	print(name)	
	#filename_output = "/Users/edo/Documents/Progetti/ImageGPCR/Zernike_Best_C4/" + "Zernike_GPCR_" +  str(name) + ".txt"
	filename_output = "../Zernike_Best_C4/" + "Zernike_GPCR.txt"

	bs_Ab = cv2.imread(filepath)
	bs_Ab = cv2.cvtColor(bs_Ab, cv2.COLOR_BGR2GRAY)

	Cm_ab = mahotas.center_of_mass(bs_Ab, labels=None)
	R_Ab_x = np.int(np.shape(bs_Ab)[0])
	R_Ab_y = np.int(np.shape(bs_Ab)[1])
	R_Ab = np.int(np.max([R_Ab_x,R_Ab_y])*0.45)
	
	dist_max = 0
	for i in range(0,np.shape(bs_Ab)[0]-1):
		for j in range(0,np.shape(bs_Ab)[1]-1):
			el_aus = bs_Ab[i,j]
			if el_aus > 0:
				dist_aus = np.linalg.norm(np.array([i,j])-Cm_ab)
				if dist_aus > dist_max:
					dist_max = dist_aus
			
	R_Ab = 	dist_max + 0.05*dist_max
		
	print(R_Ab)
	Z_ab = mahotas.features.zernike_moments(bs_Ab, R_Ab, degree=order, cm=Cm_ab)

	#np.savetxt(filename_output, Z_ab, delimiter=",")
	#conta = conta + 1
	
	ZernAll = np.row_stack([ZernAll,Z_ab])

ZernAll = np.delete(ZernAll, (0), axis=0)
#print(ZernAll)
#print(np.shape(ZernAll))

# in R:
#ZeroInvarients = [2,4,5,6,8,10,11,12,14,16,17,18,19,20,22,24,26,27,28,29,30,32,34,36,37,38,39,40,41,42,44,46,48,50,51,52,53,54,55,56,58,60,62,64,65,66,67,68,69,70,71,72,74,76,78,80,82,83,84,85,86,87,88,89,90,92,94,96,98,100,101,102,103,104,105,106,107,108,109,110,112,114,116,118,120,121]

# in python:
ZeroInvarients = [1,3,4,5,7,9,10,11,13,15,16,17,18,19,21,23,25,26,27,28,29,31,33,35,36,37,38,39,40,41,43,45,47,49,50,51,52,53,54,55,57,59,61,63,64,65,66,67,68,69,70,71,73,75,77,79,81,82,83,84,85,86,87,88,89,91,93,95,97,99,100,101,102,103,104,105,106,107,108,109,111,113,115,117,119,120]

#print(ZeroInvarients)

ZernAll_red = ZernAll[:,ZeroInvarients]
#print(ZernAll_red)
#print(np.shape(ZernAll_red))

ave_vet = np.mean(ZernAll_red, axis = 1)
#print(ave_vet)
index_min = np.argmin(ave_vet)
print(index_min)

ZernBestC4 = ZernAll[index_min,:]
#print(ZernBestC4)

np.savetxt(filename_output, ZernBestC4, delimiter=",")

bestFigFile = ListJpeg[index_min]
print(bestFigFile)
os.system('mv ' + str(bestFigFile) + ' ../Figure_Best_C4/FIG_GPCR.png')
