import numpy
from matplotlib import pyplot as plt

#########
# define rewrap matrix
##############
def wrap(t, gridwidth):
    twoDt=[[0]*gridwidth for i in range (gridwidth)]
    i = 0
    j= 0
    while i < gridwidth:
        while j < gridwidth:
	    twoDt[i][j]= complex(t[i*gridwidth+j])
            j = j + 1
        i = i + 1
        j = 0
    return twoDt
mapcount=0
map_number=5
gridwidth=128


# Create file names
datafiles=[]
CPUfiles=[]
DFEfiles=[]
dataname = "5maps/data"
DFEname = "DFEresultsreal"
CPUname = "CPUresultsreal"
while mapcount < map_number:
	datafiles.append( str( dataname + str(mapcount+1)  ) )
	CPUfiles.append( str( CPUname + str(mapcount+1)+ ".txt"  ) )
	DFEfiles.append( str( DFEname + str(mapcount+1)+ ".txt"  ) )
	mapcount=mapcount+1

#read files in
data=[]
cpu=[]
dfe=[]
mapcount=0;
while mapcount < map_number:
	readfile1 = open(datafiles[mapcount], 'r')
	read1 = readfile1.readlines()
	readfile2 = open(CPUfiles[mapcount], 'r')
	read2 = readfile2.readlines()
	readfile3 = open(DFEfiles[mapcount], 'r')
	read3 = readfile3.readlines()
	readfile1.close()
	readfile2.close()
	readfile3.close()

	i=0
	while i < gridwidth*gridwidth:
		data.append(float(read1[i]))
		cpu.append(float(read2[i]))
		dfe.append(float(read3[i]))
		i=i+1
	print((mapcount+1), " map imported")
	mapcount=mapcount+1

mapcount=0;
while mapcount < map_number:
	data2=[]
	dfe2=[]
	cpu2=[]
	i=0
	while i < gridwidth*gridwidth:
		cpu2.append(cpu[i+gridwidth*gridwidth*mapcount])
		dfe2.append(dfe[i+gridwidth*gridwidth*mapcount])
		data2.append(data[i+gridwidth*gridwidth*mapcount])
		i=i+1
	

	data3=numpy.real(numpy.array(wrap(data2, gridwidth)))
	dfe3=numpy.real(numpy.array(wrap(dfe2, gridwidth)))
	cpu3=numpy.real(numpy.array(wrap(cpu2, gridwidth)))

	textsize=20
	plt.figure(figsize=(27.0,9.0))
	plt.subplot(1,3,1),plt.imshow((data3), cmap = 'jet', interpolation='nearest')
	plt.ylabel('Data', fontsize=textsize)
	plt.colorbar()
	plt.subplot(1,3,2),plt.imshow(cpu3, cmap = 'jet', interpolation='nearest')
	plt.ylabel('CPU filtered', fontsize=textsize)
	plt.colorbar()
	plt.subplot(1,3,3),plt.imshow(dfe3, cmap = 'jet', interpolation='nearest')
	plt.ylabel('DFE filtered', fontsize=textsize)
	plt.colorbar()
	name = str( "python_result_maps/Results_map_" + str(mapcount+1) + ".png" )
	plt.savefig(name, dpi = 600)
	#plt.show()
	print((mapcount+1), " map saved")
	mapcount=mapcount+1
	


