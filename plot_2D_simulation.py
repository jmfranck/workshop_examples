from pyspecdata import *
from mayavi import mlab
import re

fp = open('TMPDMP17.fit','r')
datalist = []
thisshape = [0,0]
thisdimnames = ["",""]
thisunits = ["",""]
for j,line in enumerate(fp):
    if j>3:
        datalist.append(array([float(line[k*10:k*10+10]) for k in range(len(line)/10)]))
    elif j>2:
        start1,step1,start2,step2,datamin,datamax = [float(line[k*10:k*10+10]) for k in range(len(line)/10)]
    elif j>0:
        thisshape[j-1],_,_,thisdimnames[j-1],thisunits[j-1] = [int32(line[k*10:k*10+10]) if k<3 else line[30+5:30+7] if k<4 else line[30+8:] for k in range(5)]
    else:
        print lsafen([line[k*10:k*10+10] for k in range(len(line)/10)])
fp.close()
unit_regexp = re.compile('\((.*)\)')
for j in range(len(thisunits)):
    m = unit_regexp.search(thisunits[j])
    if m:
        print m.groups()
        thisunits[j] = m.groups()[0]
    else:
        raise ValueError("I can't parse the units in "+repr(thisunits[j]))
print lsafen("shape:",thisshape,"dimnames",thisdimnames,"units",thisunits)
if thisunits != ['MHz','MHz']:
    raise ValueError('The units are not MHz')
datalist = r_[tuple(datalist)]#.reshape(thisshape) #likely the header gives the dimension, though
data = nddata(datalist,thisshape,thisdimnames)
for j in [0,1]:
    data.setaxis(thisdimnames[j],1e6*(start1 + step1*r_[0:ndshape(data)[thisdimnames[j]]])) # 1e6 for MHz
    data.set_units(thisdimnames[j],'Hz')
#fl = figlist_var(mlab = mlab)
fl = figlist_var()
fl.next('2D simulation')
#data = data['f1':(-60e6,70e6)]['f2':(-70e6,60e6)]
fl.mesh(data)
fl.show('test_mlab_140717.png')
