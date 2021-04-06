def getDensity(filename,filename_out):
    datanum=0
    density_l=[]
    density_v=[]
    density_h=[]
    molV_l=[]
    molV_v=[]
    molV_h=[]
    cp_l=[]
    cp_v=[]
    cp_h=[]
    comp_l=[]
    comp_v=[]
    comp_h=[]
    alldata=linecache.getlines(filename)
    totals=len(alldata)
    # fpout=open(filename_out,'w')
    # pbar.printProgressBar(0, totals, prefix = 'GetDens:', suffix = 'Complete', length = 50)
    for i in range(0,len(alldata)):
        if(alldata[i][0:19]=='Please enter T in C'):
            i=i+1
            flag_str=alldata[i][0:49]
            # *  Fluid is in the single phase state           *
            # *  Fluid is in the V + L phase state            *
            # *  Fluid is in the V + H phase state            *
            # *  Fluid is in the L + H phase state            *
            if(flag_str=='*  Fluid is in the single phase state           *'):
                # print('single phase')
                i=i+3
                str_line=alldata[i]
                density_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                molV_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                cp_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                comp_l.append(str_line.split()[2])

                density_v.append(0)
                density_h.append(0)
                
                molV_v.append(0)
                molV_h.append(0)
                
                cp_v.append(0)
                cp_h.append(0)
                
                comp_v.append(-100)
                comp_h.append(-100)
                # datanum=datanum+1
            elif(flag_str=='*  Fluid is in the V + L phase state            *'):
                # print('V+L')
                #vapor
                i=i+3
                str_line=alldata[i]
                density_v.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                molV_v.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                cp_v.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                comp_v.append(str_line.split()[2])
                #liquid
                i=i+3
                str_line=alldata[i]
                density_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                molV_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                cp_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                comp_l.append(str_line.split()[2])
                #halit
                density_h.append(0)
                molV_h.append(0)
                cp_h.append(0)
                comp_h.append(-100)
            elif(flag_str=='*  Fluid is in the V + H phase state            *'):
                # print('V+H')
                #vapor
                i=i+3
                str_line=alldata[i]
                density_v.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                molV_v.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                cp_v.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                comp_v.append(str_line.split()[2])
                #salt
                i=i+3
                str_line=alldata[i]
                density_h.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                molV_h.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                cp_h.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                comp_h.append(str_line.split()[2])
                #liquid
                density_l.append(0)
                molV_l.append(0)
                cp_l.append(0)
                comp_l.append(-100)
            elif(flag_str=='*  Fluid is in the L + H phase state            *'):
                # print('L+H')
                # datanum=datanum+1
                #Liquid
                i=i+3
                str_line=alldata[i]
                density_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                molV_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                cp_l.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                comp_l.append(str_line.split()[2])
                #salt
                i=i+3
                str_line=alldata[i]
                density_h.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                molV_h.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                cp_h.append(str_line.split()[2])
                i=i+1
                str_line=alldata[i]
                comp_h.append(str_line.split()[2])
                #vapor
                density_v.append(0)
                molV_v.append(0)
                cp_v.append(0)
                comp_v.append(-100)
            else:
                print(flag_str)
        # pbar.printProgressBar(i+1, totals, prefix = 'GetDens:', suffix = 'Complete', length = 50)
    # print (len(density_h))
    linecache.clearcache()
