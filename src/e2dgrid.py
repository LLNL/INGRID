def ReadTranData(path,restore=0,save=0):
    """ Function to read all data in subfolders of path
    
    Subfolders are expected to contain ASCII files produced by the fortran routine readtran_holm.f
    
    Parameters:
    path    -   Absolute or relative path to data to be read
    restore -   Switch whether to restore already-read pickled data or re-read it
                =1 (default): the data is restored from a pickle file "data.p"
                =0 :    the data is read from the files and returned in a list of 
                        dicts according to the number of subfolders
    save    -   Switch telling the function to write the data 
                =1 : the read data is saved to pickle "data.p" if restore=0
                =0 : no writing to pickle

    Output
    Returns dictionary with all read data
    """
    from os import walk
    from os.path import isdir
    from io import open
    from csv import reader
    from numpy import asarray
    from csv import field_size_limit
    from sys import maxsize
    #from cPickle import load,dump


    def findparname(filename):
        """Returns parameter name if file to be read, otherwise returns False"""
        if "tran_" in filename:                     # Assert right file type
            for nm in ["NIMB","readtran","indices"]:    # Exclusion list
                if nm in filename:
                    return False                        # Dont read data not conforming to E2D tran files
            else:
                return filename.split("tran_",1)[1]     # Return variable name
        else: 
            return False                            # If not, return False

    field_size_limit(maxsize)
    if restore:
        dat=load(open(path+"/data.p","rb"))
        print("Loaded variables:")
        if isinstance(dat,list):
            a=dat[0].keys()
            a.sort()
            print(a)
        else:
            a=dat.keys()
            a.sort()
            print(a)
        return dat
    else:
        # Read subfolder names into variable
        casefolders=[x[0] for x in walk(path)]
        casefiles=[x[2] for x in walk(path)]

        # Loop through all folders, saving the data from tran files to dict data
        if len(casefolders)>1:
            dat=[]
            for i in range(1,len(casefolders)):
                casebuff=dict()
                for j in casefiles[i]:
                    if findparname(j):
                        with open(casefolders[i]+"/"+j,newline="") as f:
                            reader=csv.reader(f)
                            listbuff=next(reader)
                            listbuff=[float(n) for n in listbuff[0].split()]
                        casebuff[findparname(j)]=asarray(listbuff)
                dat.append(casebuff)
            print("Read variables:")
            a=dat[0].keys()
            a.sort()
            print(a)
        else:
            dat=dict()
            for param in casefiles[0]:
                if findparname(param):
                    with open(casefolders[0]+"/"+param,newline="") as f:
                        rdr=reader(f)
                        listbuff=next(rdr)
                        if param=="tran_POT":
                            stop=1
                        listbuff=[float(n) for n in listbuff[0].split()]
                    dat[findparname(param)]=asarray(listbuff)
            rdr=reader(open(casefolders[0]+"/tran_indices",newline=""))
            params=next(rdr)[0].split()
            values=next(rdr)[0].split()
            for i in range(len(params)):
                dat[params[i]]=int(values[i])
            
            print("Read variables:")
            a=dat.keys()
            sorted(a)
            print(a)
                        

        if save:
            dump(dat,open(path+"/data.p","wb"))
            print("Saved pickle")
        return dat




def calcUEindices(data,cut=False):
    """ Function parsing all UEDGE indices from the data dictionary """
    from numpy import reshape,transpose,roll,zeros, fliplr,append

    # TODO not verified to work for all cases

    # Quirks:
    #   - E2D core is rotated (core cut does not start above PFR cut)
    #   - Sometimes there is overlap for the core cells: hence the "cut" oarameter to alleviate this.
    #   - E2D only includes nodes for RVERTP and ZVERTP - all other values CC values
    #   - 

    ret={}
    names=["GRD",""]
    # NP - total number of cells, incl guard cells
    # IOPEN - First open FT
    # NXW - wall ring index: core boundary to wall
    # NC - number of rings in total
    # NROW - number of rows/poloidal indices
    # JPRGT - right normal to cut core (ixpt2?)
    # JPLFT - left normal to cut core (ixpt1?)

    for guard in [0,1]:
        solrow=data["NROW"]-2*(guard==False)                  # Poloidal cells, considering guard cells
        odivrow=data["JPRGT"]-1-(guard==False)                # Number of rows in outer leg (GC cons.)
        idivrow=data["NROW"]-data["JPLFT"]-(guard==False)   # Number of rows in inner leg (GC cons.)
        pfrrow=odivrow+idivrow                          # N.o. PFR rows
        corerow=data["JPLFT"]-data["JPRGT"]+1-cut+(guard==True)   # N.o. core rows, considering cut 
        solring=data["NXW"]-data["IOPEN"]+(guard==True)           # Number of rings in SOL
        corering=data["IOPEN"]-1-(guard==False)               # Number of rings in core
        pfrring=data["NC"]-data["NXW"]-(guard==False)       # Number of rings in PFR

        datlen=int((guard==False)*(len(data["RVERTP"])/5)+guard*(data["NP"]))    # Number of cells depending on GC or not
        core=transpose(reshape(range(corerow*corering),(corering,corerow)))   # Core indices
        # If grid index, roll core region
        # TODO: does this actually hold only for non.GC or in general? Check whether to roll all values!
        if guard==False: core=roll(core,1,axis=0)
        sol=transpose(reshape(range(corerow*corering,solring*solrow+corerow*corering),(solring,solrow))) # SOL indices
        pfr=transpose(reshape(range(solring*solrow+corerow*corering,datlen),(pfrring,pfrrow))) # PFR indices
        ret[names[guard]+"ind"]=zeros((solrow,solring+pfrring))     # Index array
        ret[names[guard]+"ind"][:,pfrring:]=sol                        # Add SOL as is
        ret[names[guard]+"ind"][:odivrow,:pfrring]=pfr[:odivrow,:]     # Add odiv
        ret[names[guard]+"ind"][odivrow:-idivrow,:pfrring]=core[:len(core)-guard,-pfrring:] # Add core w/o extra inner most cells
        ret[names[guard]+"ind"][-idivrow:,:pfrring]=pfr[-idivrow:,:]   # Add outer linner div leg
        ret[names[guard]+"ind"]=transpose(fliplr(transpose(ret[names[guard]+"ind"].astype(int))))   # Flip to get UE direcctions
    
        # Store most important geom. parameters
        ret[names[guard]+"nx"]=solrow
        ret[names[guard]+"ny"]=solring+pfrring
        ret[names[guard]+"ixpt1"]=idivrow-1
        ret[names[guard]+"ixpt2"]=solrow-odivrow-1
        ret[names[guard]+"iysptrx"]=pfrring-1

        # Store vessel grid points
    ret["rves"]=append(data["RVESM1"],data["RVESM1"][0])
    ret["zves"]=append(-data["ZVESM1"],-data["ZVESM1"][0])

    return ret


def addGC(geo,delta):
    """ Adds delta thick guard cells around the "rm" and "zm" of geo """
    from numpy import zeros,sqrt,mod


    # Pad in 2D w/ zeros
    rm,zm=geo["rm"],geo["zm"]  # Store arrays before extending
    shp=tuple(map(sum,zip(rm.shape,(2,2,0))))
    geo["rm"],geo["zm"]=zeros(shp),zeros(shp)
    geo["rm"][1:-1,1:-1],geo["zm"][1:-1,1:-1]=rm,zm


    # Outer wall GC
    bound,prevring=-1,-2
    for i in range(1,shp[0]-1):               # Loop poloidally, omit corners
        for var in ["rm","zm"]:
            for node in [1,2]:  # SW and SE corners are NW and NE of the cell below, respectively
                geo[var][i,bound,node]=geo[var][i,prevring,node+2]
            for node in [3,4]:  # Norther corners are extrapolated d in radial direction at the W/E cell faces
                # Calculate the length of each cell face
                l=sqrt((geo["rm"][i,prevring,node]-geo["rm"][i,prevring,node-2])**2+(geo["zm"][i,prevring,node]-geo["zm"][i,prevring,node-2])**2)
                # Extrapolate d in direction radially out towards boundary using inner ring's sides
                geo[var][i,bound,node]=geo[var][i,bound,node-2]+(delta/l)*(geo[var][i,prevring,node]-geo[var][i,prevring,node-2])
            geo[var][i,bound,0]=0.25*sum(geo[var][i,bound,:])    # Center is weighed mean of nodes

    # Inner wall GC
    bound,prevring,=0,1
    for i in range(1,shp[0]-1):               # Loop poloidally, omit corners
        for var in ["rm","zm"]:
            for node in [3,4]:  # SW and SE corners are NW and NE of the cell below, respectively
                geo[var][i,bound,node]=geo[var][i,prevring,node-2]
            for node in [1,2]:  # Norther corners are extrapolated d in radial direction at the W/E cell faces
                # Calculate the length of each cell face
                l=sqrt((geo["rm"][i,prevring,node]-geo["rm"][i,prevring,node+2])**2+(geo["zm"][i,prevring,node]-geo["zm"][i,prevring,node+2])**2)
                # Extrapolate d in direction radially out towards boundary using inner ring's sides
                geo[var][i,bound,node]=geo[var][i,bound,node+2]+(delta/l)*(geo[var][i,prevring,node]-geo[var][i,prevring,node+2])
            geo[var][i,bound,0]=0.25*sum(geo[var][i,bound,:])    # Center is weighed mean of nodes

    # Inner plate GC
    bound,prevring=0,1
    for i in range(1,shp[1]-1):               # Loop poloidally, omit corners
        for var in ["rm","zm"]:
            for node in [2,4]:  # SW and SE corners are NW and NE of the cell below, respectively
                geo[var][bound,i,node]=geo[var][prevring,i,node-1]
            for node in [1,3]:  # Norther corners are extrapolated d in radial direction at the W/E cell faces
                # Calculate the length of each cell face
                l=sqrt((geo["rm"][prevring,i,node]-geo["rm"][prevring,i,node+1])**2+(geo["zm"][prevring,i,node]-geo["zm"][prevring,i,node+1])**2)
                # Extrapolate d in direction radially out towards boundary using inner ring's sides
                geo[var][bound,i,node]=geo[var][bound,i,node+1]+(delta/l)*(geo[var][prevring,i,node]-geo[var][prevring,i,node+1])
            geo[var][bound,i,0]=0.25*sum(geo[var][bound,i,:])    # Center is weighed mean of nodes


    # Inner plate GC
    bound,prevring=-1,-2
    for i in range(1,shp[1]-1):               # Loop poloidally, omit corners
        for var in ["rm","zm"]:
            for node in [1,3]:  # SW and SE corners are NW and NE of the cell below, respectively
                geo[var][bound,i,node]=geo[var][prevring,i,node+1]
            for node in [4,2]:  # Norther corners are extrapolated d in radial direction at the W/E cell faces
                # Calculate the length of each cell face
                l=sqrt((geo["rm"][prevring,i,node]-geo["rm"][prevring,i,node-1])**2+(geo["zm"][prevring,i,node]-geo["zm"][prevring,i,node-1])**2)
                # Extrapolate d in direction radially out towards boundary using inner ring's sides
                geo[var][bound,i,node]=geo[var][bound,i,node-1]+(delta/l)*(geo[var][prevring,i,node]-geo[var][prevring,i,node-1])
            geo[var][bound,i,0]=0.25*sum(geo[var][bound,i,:])    # Center is weighed mean of nodes


    # Set corner's known boundaries
    corners=[[0,0],[0,-1],[-1,0],[-1,-1]]           # Corner indices
    setpoints=[[2,3,4],[1,2,4],[1,3,4],[1,2,3]]     # Each corner's corresponding already known nodes
    polcell=[[1,0,0],[0,0,1],[-2,-1,-1],[-1,-1,-2]] # Poloidal neighbor of cell for node
    radcell=[[0,1,1],[-2,-2,-1],[0,1,1],[-2,-2,-1]] # Radial neighbor cell for node
    sharednode=[[1,1,2],[3,4,3],[2,1,2],[3,4,4]]    # Shared node n of [corners[0],corners[1],node[n]] in cell [polcell[c],radcell[c]]
    for c in range(4):  # Loop over corners
        for n in range(3):  # Loop over shared nodes
            for var in ["rm","zm"]: # Loop over coords
                # Sets vector of the three setpoints to the values defined by neighboring cells
                geo[var][corners[c][0],corners[c][1],setpoints[c][n]]=geo[var][polcell[c][n],radcell[c][n],sharednode[c][n]]

    # Do corners by mirroring in diagonal
    mirrorpoints=[[4,1],[2,3],[3,2],[1,4]]  # Mirror points: [point mirrorred, point to be placed]
    vecs=[[2,3],[1,4],[4,1],[3,2]]          # Indices for vectors pointing to the two set vectors
    for c in range(4):  # Loop over each corner
        for var in ["rm","zm"]: # Loop over each coordinate
            # Place last point symmetrically by vactor operation: p_placed=p_adjacent1+p_adjacent2-p_opposite
            geo[var][corners[c][0],corners[c][1],mirrorpoints[c][1]]=geo[var][corners[c][0],corners[c][1],vecs[c][0]]+geo[var][corners[c][0],corners[c][1],vecs[c][1]]-geo[var][corners[c][0],corners[c][1],mirrorpoints[c][0]]
            # Place center point as weighed arithmetic mena
            geo[var][corners[c][0],corners[c][1],0]=0.25*sum(geo[var][corners[c][0],corners[c][1],:])
    
    # Modify the indices of ixpt1, ixpt2, iysptrx, nx and ny if needed as cells are added!
    shp=geo["rm"].shape
    if mod(geo["nx"]-shp[0],2)==1 or mod(geo["ny"]-shp[1],2)==1:
        print("Uneven number of cells added in either pol or rad dir comp. no nx/ny. Terminating!")
        exit(0)
    else:
        geo["nx"]+=int(shp[0]-geo["nx"])
        geo["ny"]+=int(shp[1]-geo["ny"])
        geo["ixpt1"]+=int((shp[0]-geo["nx"])/2)
        geo["ixpt2"]+=int((shp[0]-geo["nx"])/2)
        geo["iysptrx"]+=int((shp[1]-geo["ny"])/2)

def createNeigh(geo): 
    """ Creates stencils "neighx" and "neighy" for neighboring cell indices.
    
    The 9-point stencil is seen below:
    ------------------
    |  8  |  1  |  5  |
    ------------------
    |  4  |  0  |  2  |
    ------------------
    |  7  |  3  |  6  |
    ------------------
    """
    from numpy import zeros,asarray 

    if geo["rm"].shape[0]!=geo["nx"] or geo["rm"].shape[1]!=geo["ny"]:
        print("'rm' grid size does not coincide with 'nx' and/or 'ny' sizes! Terminating!")
        exit(0)


    stencil=[   ["ix","iy"],                                                        # Center
                ["ix","iyp1"],["ixp1","iy"],["ix","iym1"],["ixm1","iy"],            # Major directions
                ["ixp1","iyp1"],["ixp1","iym1"],["ixm1","iym1"],["ixm1","iyp1"]]    # Minor directions

   # CREATE NEIGHBORING INDICES
    geo["shp"]=geo["rm"].shape # Get correct shape

    # Allocate space
    tmp={}
    for var in ["ixm1","iym1","ixp1","iyp1","ix","iy"]:
        tmp[var]=zeros(geo["ind"].shape)

    # Poloidal indices
    for i in range(geo["shp"][0]):    # Loop over cells poloidally
        tmp["ix"][i,:]=i                        # ix returns current cell
        tmp["ixm1"][i,:]=tmp["ix"][i,:]-1     # ixm1 returns previous pol cell: BOUNDS STILL BORK
        tmp["ixp1"][i,:]=tmp["ix"][i,:]+1     # ixp1 returns next pol cell: BOUNDS STILL BORK
    tmp["ixm1"][0,:]=tmp["ix"][0,:]       # Fix LB
    tmp["ixp1"][-1,:]=tmp["ix"][-1,:]     # Fix RB

    # Fix poloidal core cut
    for j in range(geo["iysptrx"]+1): # Loop over radial cells up to separatrix
        tmp["ixm1"][geo["ixpt1"]+1,j]=geo["ixpt2"]    # Close core on itself from the left
        tmp["ixm1"][geo["ixpt2"]+1,j]=geo["ixpt1"]    # Close core on itself from the right
        tmp["ixp1"][geo["ixpt1"],j]=geo["ixpt2"]+1    # Connect PFR legs from the left
        tmp["ixp1"][geo["ixpt2"],j]=geo["ixpt1"]+1    # Connect PFR legs from the right

    # Radial indices
    for j in range(geo["shp"][1]): # Loop over cells radially
        tmp["iy"][:,j]=j                      # iy returns current index
        tmp["iym1"][:,j]=tmp["iy"][:,j]-1   # iym1 returns previous radial cell: BOUNDS STILL BORK
        tmp["iyp1"][:,j]=tmp["iy"][:,j]+1   # iyp1 returns next radial cell: BOUNDS STILL BORK
    tmp["iym1"][:,0]=tmp["iy"][:,0]         # Fix LB
    tmp["iyp1"][:,-1]=tmp["iy"][:,-1]       # Fix RB
    # Make indices to ints
    for var in ["ix","iy","ixm1","ixp1","iym1","iyp1"]:
        tmp[var]=tmp[var].astype(int)

    # Allocate arrays
    geo["neighx"]=zeros(geo["ind"].shape+(9,)) # Exted 2D grid to four dimensions
    geo["neighy"]=zeros(geo["ind"].shape+(9,)) # Exted 2D grid to four dimensions

    # Populate the arrays
    for i in range(geo["shp"][0]):    # Loop poloidally
        for j in range(geo["shp"][1]):    # Loop radially
            var=["neighx","neighy"]         # Variables to be stored to
            for coord in range(2):          # Loop over coords
                buff=[]                         # Placeholder for each coord
                for neig in stencil:                # Loop over stencil indices
                    buff.append(tmp[neig[coord]][i,j])    # Append indices to buff
                geo[var[coord]][i,j,:]=asarray(buff)       # Add coordinates to array
    geo["neighx"]=geo["neighx"].astype(int)
    geo["neighy"]=geo["neighy"].astype(int)


def writeGridue(path,geo):
    """ Writes gridue file from geo """
    f=open(path+"/gridue","w+")       # Open gridue for writing
    # Write header line containing: nx ny ixpt1 ixpt2 iysptrx (BASIS INDEXING!)
    f.write("  "+str(geo["GRDnx"])+"  "+str(geo["GRDny"])+"  "+str(geo["ixpt1"])+"  "+str(geo["ixpt2"])+"   "+str(geo["iysptrx"]))
    f.write('\n\n') # Blank line required
    for var in ["rm","zm","psi","br","bz","bpol","bphi","b"]: 
                                                # Loop through variables to be written
        output=[]                               # Output array
        for n in range(5):                      # Loop through nodes
            for j in range(geo["ny"]):            # Loop through radial indices (write FT-wise)
                for i in range(geo["nx"]):            # Loop through poloidal indices
                    output.append(geo[var][i,j,n])     # Append data to
        for i in range(len(output)/3):          # Loop through each line in gridue, whic is 3 entries long
            string=""                               # String to be written
            for j in range(3):                      # Loop through each value for the row
                k=i*3+j                             # Current index in out array
                if output[k]<0:                         # If negative value, one less whitespace to accommodate for sign
                    string=string+" %1.15E" % (output[k])
                else:                                   # Else, include the additional whitespace
                    string=string+"  %1.15E" % (output[k])
            string=string+'\n'                          # End with new row
            f.write(string.replace("E","D"))            # Write to gridue, replacing "E" w/ "D" (BASIS standard)
        if var!="b":                            # If this is not the last parameter, add empty line for read
            f.write('\n')   
    f.write("   EFITD   E2D/HOLM   #160299   2240") # Write "runid"-line
    f.close()                               # Close file


def writeUEdata(data,geo,keylist):
    from numpy import zeros
    ret={}
    for var in keylist:         # Loop over coordinates 
        ret[var[1]]=zeros((geo["nx"],geo["ny"])) 
        for i in range(geo["nx"]):  # Loop poloidally
            for j in range(geo["ny"]):  # Loop radially
                ret[var[1]][i,j]=data[var[0]][geo["ind"][i,j]]
    return ret

    

def grid(path,plot=False):
    """ Reads a E2D grid file and returns [rm, zm, nx,ny, ixpt1, ixpt2, iysptrx]
    
    Parameters:
    path -          Absolute path to the E2D file to be converted. 
                    The data must have been extracted from thetran 
                    file using the script "readtran_holm.f".
    interp -        Choose model for interpolation to cell corner nodes.
                    =0 : Arithmetic mean
                    =1 (default) : Weighted L1 arithmetic mean
                    =2 : Harmonic mean (not implemented yet)

    Output:
    gridue -        Writes a UEDGE grid file under the name gridue
    E2D_2_UE.pckl - Pickle file containing the following dicts:
                    UEgeo: grid and geometric properties in UE standard
                    UEdata: any written plasma parameters in UE standard
    
    """
    from matplotlib.pyplot import figure
    from numpy import zeros

    # Read tran data to dict
    data=ReadTranData(path)



    """ RESTORE INDEX DATA """
    # Core cells: 32*14 = (NXW+2)*(JPRGT-2)
    # SOL cells: 72*14 = (NROW-2)*IOPEN(-2)
    # PFR cells: 40*5 = ((NROW-2)-(COREROWS))*
    # ODIV=14 cells
    # IDIV=26 cells
    # GRID INDICES
    UEgeo=calcUEindices(data)
    # Create geo dict with node positional data
    mpts=int(len(data["RVERTP"])/5.0)        # Number of mesh points
    geo={}
    geo["rm"]=zeros((mpts,5))        # Allocate space for meshes
    geo["zm"]=zeros((mpts,5))        

    for i in range(5):  # Loop through each node
        geo["rm"][:,i]=data["RVERTP"][i::5] # Store node data to geo
        geo["zm"][:,i]=-data["ZVERTP"][i::5]

    # Indices now correspond to:
    # 0=ul, 1=ll, 2=lr, 3=ur, 4=cc

    # UE index order_
    # cc, ll, lr, ul, ur
    # Conversion from E2D order to UE order:
    nodes=(4,3,0,2,1)
    # Convert to UE indexing
    geo["rm"]=geo["rm"][:,nodes]
    geo["zm"]=geo["zm"][:,nodes]


    # ARRANGE THE DATA INTO UE-STANDARD ARRAYS
    for var in ["rm","zm"]:         # Loop over coordinates
        UEgeo[var]=zeros((UEgeo["GRDnx"],UEgeo["GRDny"],5))
        for node in range(5):       # Loop over UE nodes
            for i in range(UEgeo["GRDnx"]):
                for j in range(UEgeo["GRDny"]):
                    UEgeo[var][i,j,node]=geo[var][UEgeo["GRDind"][i,j],node]

    # Add guard cells with width 2e-9
    addGC(UEgeo,2e-9)

    ret=[]
    for var in ['rm','zm','GRDnx','GRDny','ixpt1','ixpt2','iysptrx']:
        ret.append(UEgeo[var])

    if plot is True:
        f=figure()
        ax=f.add_subplot(111)
        for x in range((UEgeo['GRDnx'])):
            for y in range((UEgeo['GRDny'])):
                bx,by=[],[]
                for i in [1,2,4,3,1]:
                    bx.append(UEgeo['rm'][x,y,i])
                    by.append(UEgeo['zm'][x,y,i])
                ax.plot(bx,by,'k-')
                ax.plot(UEgeo['rm'][x,y,0],UEgeo['zm'][x,y,0],'k.')
            

    return ret




def recreateGrid(path,psi0=1,interp=1,restore=1,save=0):
    """ Reads a E2D grid file and writes a gridue file and picke with variables.
    
    Parameters:
    path -          Absolute path to the E2D file to be converted. 
                    The data must have been extracted from thetran 
                    file using the script "readtran_holm.f".
    interp -        Choose model for interpolation to cell corner nodes.
                    =0 : Arithmetic mean
                    =1 (default) : Weighted L1 arithmetic mean
                    =2 : Harmonic mean (not implemented yet)

    Output:
    gridue -        Writes a UEDGE grid file under the name gridue
    E2D_2_UE.pckl - Pickle file containing the following dicts:
                    UEgeo: grid and geometric properties in UE standard
                    UEdata: any written plasma parameters in UE standard
    
    """

    from numpy import zeros,sin,cos,pi,sqrt,savez
    from cPickle import dump
    from calc_holm import calcang, interpolate

    # Read tran data to dict
    if "data" not in vars():
        data=ReadTranData(path,restore=restore,save=save)



    """ RESTORE INDEX DATA """
    # Core cells: 32*14 = (NXW+2)*(JPRGT-2)
    # SOL cells: 72*14 = (NROW-2)*IOPEN(-2)
    # PFR cells: 40*5 = ((NROW-2)-(COREROWS))*
    # ODIV=14 cells
    # IDIV=26 cells
    # GRID INDICES
    UEgeo=calcUEindices(data)
    # Create geo dict with node positional data
    mpts=int(len(data["RVERTP"])/5.0)        # Number of mesh points
    geo={}
    geo["rm"]=zeros((mpts,5))        # Allocate space for meshes
    geo["zm"]=zeros((mpts,5))        

    for i in range(5):  # Loop through each node
        geo["rm"][:,i]=data["RVERTP"][i::5] # Store node data to geo
        geo["zm"][:,i]=-data["ZVERTP"][i::5]

    # Indices now correspond to:
    # 0=ul, 1=ll, 2=lr, 3=ur, 4=cc

    # UE index order_
    # cc, ll, lr, ul, ur
    # Conversion from E2D order to UE order:
    nodes=(4,3,0,2,1)
    # Convert to UE indexing
    geo["rm"]=geo["rm"][:,nodes]
    geo["zm"]=geo["zm"][:,nodes]


    # ARRANGE THE DATA INTO UE-STANDARD ARRAYS
    for var in ["rm","zm"]:         # Loop over coordinates
        UEgeo[var]=zeros((UEgeo["GRDnx"],UEgeo["GRDny"],5))
        for node in range(5):       # Loop over UE nodes
            for i in range(UEgeo["GRDnx"]):
                for j in range(UEgeo["GRDny"]):
                    UEgeo[var][i,j,node]=geo[var][UEgeo["GRDind"][i,j],node]

    # Add guard cells with width 2e-9
    addGC(UEgeo,2e-9)

    # Create neighboring cell indices arrays "neighx" and "neighy" in UEgeo
    createNeigh(UEgeo)

    # Calculate angles and store them to the geo file
    calcang(UEgeo)


    # Arrange grid data into arrays
    data["sxcc"]=2*pi*data["RMESH"]*data["HRHO"]
    # Loop through translations E2D-UE
    for var in [["PSI","psi"],["BFI","bphi"],["SH","SH"],["sxcc","sxcc"]]:  # Loop through existing cc grid data
        UEgeo[var[1]]=zeros(UEgeo["shp"])
        for i in range(UEgeo["shp"][0]):     # Loop over pol coords
            for j in range(UEgeo["shp"][1]): # Loop over rad coords
                UEgeo[var[1]][i,j,0]=data[var[0]][UEgeo["ind"][i,j]] # Store cc data in correct order

    # Allocate arrays
    for var in ["bpol","b","br","bz"]: # Loop for creating arrays
        UEgeo[var]=zeros(UEgeo["shp"])   # Create empty arrays

    # Calculate CC values
    UEgeo["bpol"][:,:,0]=(UEgeo["SH"][:,:,0]*UEgeo["bphi"][:,:,0])/sqrt(1-UEgeo["SH"][:,:,0]**2) # Calc bpol from bphi and bpol/btot-ratio
    UEgeo["b"][:,:,0]=UEgeo["bpol"][:,:,0]/sqrt(1-UEgeo["SH"][:,:,0]**2)  # Calc btot from bpol and bpol/btot-ratio
    # TODO
    """
    In grid2d source file (/home/sim/cmg/jams/v121218/depot/grid2d/source/fort/grid2d.f:L7876):
    PSI=PSI0*PSI, i.e. PSI0 is the inverse of the normalization!
    -> psi=PSI/PSI0
    For grid /home/mgroth/cmg/catalog/grid2d/d3d/160299/aug1816/seq#1: PSI0=2.283
    Units still unclear: PSI likely taken straight from EFIT equilibrium. Here PSI is given as weber/rad...
    """
    UEgeo["psi"][:,:,0]=UEgeo["psi"][:,:,0]/(2*pi*psi0) # Calculate cc-PSI from normalized data and normalization factor PSI0
    UEgeo["bz"][:,:,0]=UEgeo["bpol"][:,:,0]*sin(UEgeo["angpol"])   # Calc poloidal B-field using angles
    UEgeo["br"][:,:,0]=UEgeo["bpol"][:,:,0]*cos(UEgeo["angpol"])   # Calc radial B-field using angles 




    # Interpolate to corner nodes
    interpolate(UEgeo,interp=interp) 
    # TODO what is the E2D polodal direction? Clockwise of counterclockwise?
    # TODO what happens when the grid is flipped (ZVERTP needs to be flipped) -> ccw becomes cw! -> could explain the signs
    # TODO what is the poloidal direction in the PFR

    # Output gridue file
    writeGridue(path,UEgeo)
            
    UErawdata=writeUEdata(data,UEgeo,
        [   ["DENEL","ne"],["DEN","ni0"],["DA","ni1"],
            ["DENZ01","ni2"],["DENZ02","ni3"],["DENZ03","ni4"],["DENZ04","ni5"],["DENZ05","ni6"],["DENZ06","ni7"],
            ["VPI","up0"],["VA0","up1"],
            ["VPZ01","up2"],["VPZ02","up3"],["VPZ03","up4"],["VPZ04","up5"],["VPZ05","up6"],["VPZ06","up7"],
            ["DA","ng0"],["DM","ng1"],["DZ","ng2"],
            ["TEV","ti"],["TEVE","te"],
            ["ENEUTA","tg0"],["ENEUTM","tg1"],["ENEUTZ","tg2"],
            ["POT","phi"],
            ["HRHO","sx"]   ]   )

    # TODO : What is VA0: polidal/parallel direction? Calculate from VA0R/VA0Z/VA0T?


    UEtemp=writeUEdata(data,UEgeo, [["HRHO","sx"],["RMESH","rmesh"],["ZMESH","zmesh"]])
    UEgeo["rxpt"]=data["RPX"]
    UEgeo["zxpt"]=-data["ZPX"]
    for k in UEtemp.keys():
        UEgeo[k]=UEtemp[k]
    


    eV=1.602e-19
    # COMPILE DATA TO UEDGE ARRAYS
    UEdata={}
    # Electron temps and dens
    UEdata["ne"]=UErawdata["ne"]
    UEdata["te"]=UErawdata["te"]*eV
    # Plasma potential
    # TODO: Reverse sign due to poloidal directionality convention in E2D v UE??
    UEdata["phi"]=UErawdata["phi"]
    # Ion temps and dens
    UEdata["ti"]=UErawdata["ti"]*eV
    UEdata["ni"]=zeros(UEdata["ne"].shape+(8,))
    UEdata["ni"][:,:,0]=UErawdata["ni0"]
    UEdata["ni"][:,:,1]=UErawdata["ni1"]
    UEdata["ni"][:,:,2]=UErawdata["ni2"]
    UEdata["ni"][:,:,3]=UErawdata["ni3"]
    UEdata["ni"][:,:,4]=UErawdata["ni4"]
    UEdata["ni"][:,:,5]=UErawdata["ni5"]
    UEdata["ni"][:,:,6]=UErawdata["ni6"]
    UEdata["ni"][:,:,7]=UErawdata["ni7"]
    # Neutral densities
    UEdata["ng"]=zeros(UEdata["ne"].shape+(3,))
    UEdata["ng"][:,:,0]=UErawdata["ng0"]
    UEdata["ng"][:,:,1]=UErawdata["ng1"]
    UEdata["ng"][:,:,2]=UErawdata["ng2"]
    # Parallel velocities
    # E2D uses out-to-in directionality: reverse sign for UEDGE compatibility!!
    UEdata["up"]=zeros(UEdata["ne"].shape+(8,))
    UEdata["up"][:,:,0]=-UErawdata["up0"]
    UEdata["up"][:,:,1]=-UErawdata["up1"]
    UEdata["up"][:,:,2]=-UErawdata["up2"]
    UEdata["up"][:,:,3]=-UErawdata["up3"]
    UEdata["up"][:,:,4]=-UErawdata["up4"]
    UEdata["up"][:,:,5]=-UErawdata["up5"]
    UEdata["up"][:,:,6]=-UErawdata["up6"]
    UEdata["up"][:,:,7]=-UErawdata["up7"]
    # Neutral temperatures
    UEdata["tg"]=zeros(UEdata["ne"].shape+(3,))
    UEdata["tg"][:,:,0]=UErawdata["ti"]*eV
    UEdata["tg"][:,:,1]=UErawdata["tg1"]*eV
    UEdata["tg"][:,:,2]=UErawdata["ti"]*eV
    # Neutral densities, E2D
    UEdata["tg_all"]=zeros(UEdata["ne"].shape+(3,))
    UEdata["tg_all"][:,:,0]=UErawdata["tg0"]*eV
    UEdata["tg_all"][:,:,1]=UErawdata["tg1"]*eV
    UEdata["tg_all"][:,:,2]=UErawdata["tg2"]*eV
    # Poloidal area
    UEdata["sx"]=UErawdata["sx"]
    # Vessel data
    UEdata["xlim"]=UEgeo["rves"]
    UEdata["ylim"]=UEgeo["zves"]




    # Pretend to create a UEDGE case on the grid for writing data
    from uedge import bbb,com,grd,svr,aph,api
    from uedge.hdf5 import hdf5_save    # Import functions for handling hdf5 files
    from os import getcwd,chdir

    aph.aphdir="/fusion/projects/codes/uedge/aph"   # Hydrogen rates
    api.apidir="/fusion/projects/codes/uedge/api"   # Impurity rates 


    bbb.ngrid=1             # Number of grids to allocate space for? odesetup.m L6488
    bbb.mhdgeo=1            # Grid geometry:
    bbb.gengrid=0           #1= generates grid, 0=restores grid from gridue 
    com.nhsp=2          # N.o. hydrogenic species
    com.ngsp=3          # N.o. hydrogenic gaseous specie
    com.nhgsp=2     # Number of hydrogenic gas species
    bbb.ishymol=1       # Includes molecules as 2nd gaseous species (index 1)
    bbb.isimpon=6       # Switch for activating impurities
    species,imps=0,6            # Impurity species and impurity charge states
    com.nzsp[species]=imps          # Number of impurity species for gas species species+1
    bbb.ismctab=2       # Define data to be used for multi-charge-state rates
    com.mcfilename="C_rates.adas"   # Rate data to be used
    bbb.minu[2:8]=12      # Atomic mass unit species mass
    bbb.znuclin[2:8]=6  # Nuclear charge of impurities
    bbb.ziin[:8]=[1,0,1,2,3,4,5,6]

    curdir=getcwd()
    chdir(path)
    
    bbb.issfon=0;bbb.ftol=1e100;bbb.exmain()
    chdir(curdir)

    bbb.tis=UEdata["ti"]
    bbb.tes=UEdata["te"]
    bbb.tgs=UEdata["tg"]
    bbb.nis=UEdata["ni"]
    bbb.ngs=UEdata["ng"]
    bbb.ups=UEdata["up"]
    bbb.phis=UEdata["phi"]
    
        

    from uedge.hdf5 import hdf5_save
    hdf5_save(path+"/E2Dexport.h5")

 #   outfile=open(path+"/E2Dexport.npz","wb")
 #   savez(outfile,tis=UEdata["ti"], tes=UEdata["te"], tgs=UEdata["tg"], nis=UEdata["ni"], ngs=UEdata["ng"], ups=UEdata["up"], phis=UEdata["phi"]) 
 #   outfile.close()

    with open(path+"/E2D_2_UE.pckl", "wb") as f:
        dump(UEgeo, f)
        dump(UEdata, f)


  
