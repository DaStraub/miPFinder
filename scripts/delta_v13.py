#returns % blast hits for: (1) <=miPmaxlength, (2) miPmaxlength<x<ancestorminlength, (3) >=ancestorminlength
def percentZones(hits,fastadb,miPmaxlength,ancestorminlength):
	length = []
	for agi in hits:
		length.append(len(fastadb[1][fastadb[0].index(agi)]))
	temp = []
	if length == []:
		temp = ['failed','failed','failed']
	else:
		mipZone = sorted(i for i in length if i <= miPmaxlength)
		greyZone = sorted(i for i in length if i < ancestorminlength and i > miPmaxlength)
		ancZone = sorted(i for i in length if i >= ancestorminlength)
		temp.append(float(len(mipZone))/len(length)*100)
		temp.append(float(len(greyZone))/len(length)*100)
		temp.append(float(len(ancZone))/len(length)*100)
	return temp


#does NOT ignore splice variants
def percentsmall(hits,fastadb,ancestorminlength):
	length = []
	for agi in hits:
		length.append(len(fastadb[1][fastadb[0].index(agi)]))
	if length == []:
		nonancestorPerCent = 'failed'
	else:
		nonancestorhits = sorted(i for i in length if i < ancestorminlength)
		nonancestorPerCent = float(len(nonancestorhits))/len(length)*100
	return nonancestorPerCent
	

