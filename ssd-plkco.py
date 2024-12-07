#ssd-plkco -- sypsim_struct_discrete, plk only comes on at CO
#020210226pmc

import numpy as np
from matplotlib import pyplot as pl
from PIL import Image as im

nruns=5000
nsyps=3000
nplks=300
sbsmax=8
sypmax=3
phosinit=0.1   #initial syp-1-phos rate
ph_off=0.025     # phos-SYP off-rate
np_off=0.05     # non-phos SYP off-rate
syp_on=0.9      # SYP on-rate
fr_dep=1      # chance of getting dephosphorylated if free
bo_dep=0.0    # chance of getting dephosphorylated if bound
bo_plkloss=0.0      # chance of losing plk-2 if bound
co_pho=1      # chance of getting phosphorylated near CO
syp_dc=1     # chance of SYP moving from its current position
sbs_mult=0.1    # modifier of SYP attractiveness for phosphorylated SBS
sbs_deph=0.1   # rate of SBS dephosphorylation

class syp1:
    gen=0
    xs=-1
    pos=-1
    phos=False
    plk=False
    bindcount=0
    phoscount=0
    def __init__(self):
        self.gen=syp1.gen
        syp1.gen+=1
        self.phos=(np.random.random()<phosinit)

def getconf(fname):
    xsl=[];copos=[];
    with open(fname) as conf:
        for line in list(conf):
            if(line[0] != '#'):
                words=line.strip()
                words=words.split(' ')
                xsl.append(int(words[0]))
                copos.append([int(x) for x in words[1:]])
    return(xsl,copos)

def deph_sbs():
    for i in sbs:
        i-=((np.random.random(i.shape)<sbs_deph)*(i>0))

def get_random_loc():
    #pick a random location among all possible chromosomes
    tloc=int(np.random.random()*txs)
    thexs=np.min(np.where(cxs>tloc))-1
    thepos=tloc%xsl[thexs]
    return(thexs,thepos)

def get_random_syps(n):
    #pick n random SYP instances without repeating
    pos=np.random.permutation(nsyps)
    return([syplist[x] for x in pos[0:n]])

def get_random_syp():
    #pick a random SYP instance
    pos=int(np.random.random()*len(syplist))
    return(syplist[pos])

def set_syp(xsl):
    x=[np.zeros(i) for i in xsl]
    xp=[np.zeros(i) for i in xsl]
    for i in syplist:
        if(i.xs>=0):
            x[i.xs][i.pos]+=1
            if(i.phos):
                xp[i.xs][i.pos]+=1
    return x,xp

def set_sbs():
    for i in syplist:
        r=sbs[i.xs][i.pos]
        if(r<sbsmax and i.plk):
            sbs[i.xs][i.pos]+=1


def syp_bind():
    for i in syplist:
        if(i.xs < 0):                                       # is it a free SYP?
            if(np.random.random()<syp_on):                  # it has a chance to bind...
                thexs,thepos=get_random_loc()
                if(syp[thexs][thepos]<sypmax):
                    i.xs=thexs
                    i.pos=thepos
                    i.bindcount+=1
            else:                                           # otherwise it can get dephosphorylated
                if(np.random.random()<fr_dep):
                    i.phos=False

def syp_step():
    global nplks
    for i in syplist:
        if(i.xs >= 0):                                       # is it a bound SYP?
            cos=copos[i.xs]
            if((i.pos in cos) or (i.pos-1 in cos)):             # it's on either side of a CO, maybe phosphorylate
                if((i.phos == False) and (np.random.random()<co_pho)):
                    i.phos = True
                    i.phoscount += 1
                    if((i.plk == False) and (nplks>0)):
                        i.plk = True ## experimenting with adding PLK always at CO. 
                        i.bindcount+=1
                        #beware! syp_step goes in array order, so there may be order effects here that make things wack. But since they bind randomly, it might be ok.
                        nplks-=1

            sbsfactor=sbs_mult ** sbs[i.xs][i.pos]
            if(i.phos):
                takeoff=ph_off*sbsfactor
            else:
                #takeoff=np_off
                takeoff=np_off*sbsfactor #20210302pmc--to make SYP-nonphos insensitive to sbs

            if(np.random.random()<takeoff):         #SYP comes off the chromosome, PLK comes off by necessity
                i.xs=-1
                i.pos=-1
                if(i.plk):
                    i.plk=False
                    nplks+=1

            else:                                   #lateral diffusion    
                if(np.random.random() < syp_dc):    #it's going to move...
                    if(i.pos==0 or (i.pos-1) in cos): #move right only
#                        if(syp[i.xs][i.pos]>=syp[i.xs][i.pos+1]):
                        if(np.random.random()<0.5):
                            i.pos+=1
                    elif(i.pos==(xsl[i.xs]-1)  or (i.pos in cos)): #move left only
#                        if(syp[i.xs][i.pos]>=syp[i.xs][i.pos-1]):
                        if(np.random.random()<0.5):
                            i.pos-=1
                    else: #not at end or CO, need to look right and left
                        if(np.random.random()<.5):
                            i.pos-=1
                        else:
                            i.pos+=1
            if(np.random.random() < bo_dep):     #experimenting with bound-SYP dephos
                if(i.phos):
                    i.phos=False
                    if(i.plk):
                        i.plk=False
                        nplks+=1
            if(np.random.random() < bo_plkloss): #experimenting with PLK removal from SYP-1phos on SC
                if(i.plk):
                    i.plk=False
                    nplks+=1

xsspacer=np.array(([10,10,10],[0,0,0],[0,0,0]))

xsl,copos=getconf('ssd.conf')
txs=sum(xsl) #total number of positions
cxs=np.cumsum(np.c_[0,[xsl]])

syp=[np.zeros(l) for l in xsl]   # SYP-1 locations
sypp=[np.zeros(l) for l in xsl]   # SYP-1phos locations
sbs=[np.zeros(l) for l in xsl]  # PLK-2 substrate

syplist=[syp1() for x in range(nsyps)]

xs2=np.zeros((nruns,txs+(3*(len(xsl)-1)),3))
sb2=np.zeros((nruns,txs+(3*(len(xsl)-1)),3))
for i in range(nruns):
    syplist=get_random_syps(nsyps)
    deph_sbs()
    set_sbs()
    syp,sypp=set_syp(xsl)
    syp_bind()
    syp_step()

    xs2[i,:,0]+=np.concatenate([np.concatenate((x,xsspacer[0])) for x in syp])[:-3]
    xs2[i,:,1]+=np.concatenate([np.concatenate((x,xsspacer[1])) for x in sypp])[:-3]
    xs2[i,:,2]+=np.concatenate([np.concatenate((x,xsspacer[2])) for x in syp])[:-3]
    for color in range(3):
        sb2[i,:,color]=np.concatenate([np.concatenate((x,xsspacer[color])) for x in sbs])[:-3]

pl.imsave("xs.png",xs2/np.max(xs2))
pl.imsave("sbs.png",sb2/np.max(sb2))

print("Bound+nonphos",sum([i.xs>=0 and i.phos==False for i in syplist]))    
print("Bound+phos",sum([i.xs>=0 and i.phos for i in syplist]))    
print("Bound+phos+plk",sum([i.xs>=0 and i.phos and i.plk for i in syplist]))    
print("Free+nonphos",sum([i.xs<0 and i.phos==False for i in syplist]))
print("Free+phos",sum([i.xs<0 and i.phos for i in syplist]))
print("gen ",i, "PLKs: ",nplks)
