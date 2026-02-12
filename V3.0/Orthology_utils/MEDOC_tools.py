import numpy as np
import math
global refs


refs=[['GLX','LYX','ASX','TYX','SXP','TXP','YXP','HDX','HEX','HIX'],
     [['E','e'],['K','k'],['D','d'],['Y','y'],['U','u'],['X','x'],['Z','z'],['H','9'],['H','8'],['H','h']],
     ['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN'],['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN'],
     ['','','','','','','','',''],[1,-1,1,1,1,1,1,-1,-1,-1],[4.34,10.34,3.86,9.76,5.96,6.3,5.96,7.15,6.55,6.45],
     [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN'],
     ['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN'],['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN'],
     ['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN']]

def convert_AA_1_to_3(seq,mode=1) :
    dic={   "A" : "ALA",
            "C" : "CYS",
            "D" : "ASP",
            "d" : "ASH",
            "E" : "GLU",
            "e" : "GLH",
            "F" : "PHE",
            "G" : "GLY",
            "H" : "HIP",
            "h" : "HID",
            "8" : "HEX",
            "9" : "HDX",
            "I" : "ILE",
            "K" : "LYS",
            "k" : "LYD",
            "L" : "LEU",
            "M" : "MET",
            "N" : "ASN",
            "P" : "PRO",
            "Q" : "GLN",
            "R" : "ARG",
            "S" : "SER",
            "T" : "THR",
            "V" : "VAL",
            "W" : "TRP",
            "y" : "TYR",
            "Y" : "TYO",
            "z" : "NME",
            "_" : ""}
    if mode==1:
        W=''
        for i in range(len(seq)):
            try :
                W+=dic[seq[i]]+"\n"
            except :
                print("Could not find "+seq[i]+", exiting.")
                exit(1)
    if mode==2:
        W=[]
        for i in range(len(seq)):
            try :
                W+=[dic[seq[i]]]
            except :
                print("Could not find "+seq[i]+", exiting.")
                exit(1)

    return W

def convert_AA_3_to_1_letter(input_arr,silent=False):
    dico = {"ALA" : "A",
            "CYS" : "C",
            "ASP" : "D",
            "ASH" : "d",
            "ASX" : "D",
            "GLU" : "E",
            "GLH" : "e",
            "GLX" : "E",
            "PHE" : "F",
            "GLY" : "G",
            "HIS" : "H",
            "HID" : "h",
            "HIE" : "h",
            "HIP" : "H",
            "HIX" : "H",
            "HEX" : "8",
            "HDX" : "9",
            "ILE" : "I",
            "LYS" : "K",
            "LYD" : "k",
            "LYX" : "K",
            "LEU" : "L",
            "MET" : "M",
            "ASN" : "N",
            "PRO" : "P",
            "GLN" : "Q",
            "SER" : "S",
            "TYR" : "y",
            "TYO" : "Y",
            "TYX" : "Y",
            "THR" : "T",
            "VAL" : "V",
            "TRP" : "W",
            "ACE" : "z",
            "ARG" : "R",
            "NME" : "z"}
    output_arr=[]
    for i in range(len(input_arr)):
        try :
            output_arr+=[dico[input_arr[i]]]
        except :
            if not silent:
                print("Did not convert "+str(input_arr[i]))
    return output_arr

def get_base_contexts(map,neigh,seq_1,seq_data_q):
    # So the first part : get all the context. Read only once.
    contexts=[]
    states=[]
    #Now going to get actual pattern
    for i in range(len(map)):
        tmp=''
        tmp2=[]
        for j in range(-neigh,neigh+1):
            #if not outside of the sequence
            if map[i]+j>=1 and map[i]+j<len(seq_1)-1:
                tmp+=seq_1[map[i]+j].upper()
                tmp2+=[seq_data_q[map[i]+j]]
            # Just padding
            else :
                tmp+='_'
                tmp2+=[float('nan')]
        contexts+=[tmp]
        states+=[tmp2]
    return contexts,states

def get_G_sum(Garr,T):
    #New now gets rid of infinites
    Garr=Garr[np.isfinite(Garr)]
    if len(Garr)==0:
        return np.longdouble('+inf')
    # This offset is to guaranty there is no float overflow :
    # Since the energy being compared are always within mesostates, they are always close to each other.
    # Using the first element is faster than mean
    # The offset is then reapplied, which check out mathematically (eq. 7 and  in the manuscript)
    offset=Garr[0]
    Garr=Garr-offset
    out=-R*T*np.log(np.sum(np.exp(-(Garr)/(R*T)),axis=0))+offset
    return out

def get_all_contexts(contexts,states,neigh,refs,unshifted,reverse,penta,additive,T,base_rep,reduced):
    import itertools
    missing_pat=False
    p=0
    pat_seq=[[] for i in range(len(contexts))]
    pat_ste=[[] for i in range(len(contexts))]
    pat_E=[[] for i in range(len(contexts))]

    # Now to get the contexts combinatorics
    # 0 is proton bound
    # 1 is proton unbound
    for i in range(len(contexts)):
        count=0
        new_seq=''
        titratable_pre=[]
        for j in range(neigh):
            if contexts[i][j]=='E' or contexts[i][j]=='K' or contexts[i][j]=='H' or contexts[i][j]=='D' or contexts[i][j].upper()=='Y':
                titratable_pre+=[1]
            else :
                titratable_pre+=[0]
        titratable_post=[]
        for j in range(neigh+1,neigh*2+1):
            if contexts[i][j]=='E' or contexts[i][j]=='K' or contexts[i][j]=='H' or contexts[i][j]=='D' or contexts[i][j].upper()=='Y':
                titratable_post+=[1]
            else :
                titratable_post+=[0]
        for j in range(neigh):
            if states[i][j]==-1 or states[i][j]==1 :
                count+=1
        #inverted 0 and 1
        for curr in itertools.product([1,0],repeat=count):
            tmp=''
            new_seq=''
            offset=0
            txt_ste=''
            h=0

            for k in range(len(titratable_pre)):
                if titratable_pre[k]==1 :
                    if contexts[i][k]=='E' or contexts[i][k]=='D' or contexts[i][k].upper()=='Y':
                        if curr[h]==0:
                            tmp+=contexts[i][k].upper()
                        elif curr[h]==1:
                            tmp+=contexts[i][k].lower()
                    elif contexts[i][k]=='K' or contexts[i][k]=='H' :
                        if curr[h]==1:
                            tmp+=contexts[i][k].upper()
                        elif curr[h]==0:
                            tmp+=contexts[i][k].lower()
                    txt_ste+=str(curr[h])
                    h+=1

                else :
                    txt_ste+='1'
                    tmp+=contexts[i][k]

            new_seq+=tmp
            #central residue is always written charged
            new_seq+=contexts[i][neigh].upper()

            for k in range(neigh+1,neigh*2+1):
                #inverted lower upper
                if contexts[i][k]=='E' or contexts[i][k]=='D' or contexts[i][k].upper()=='Y':
                    new_seq+=contexts[i][k].lower()
                elif contexts[i][k]=='K' or contexts[i][k]=='H' :
                    new_seq+=contexts[i][k].upper()
                else :
                    new_seq+=contexts[i][k]

            pat_ste[i]+=[txt_ste]
            pat_seq[i]+=[new_seq]

            #This is just to have a good unshifted graph that does not skip over any of the combinatorics
            if unshifted==True :
                for j in range(len(refs[1])):
                    if refs[1][j][1].upper()==new_seq[neigh].upper() :
                        break

                F_mdcp_expt=R*T*np.log(10**(refs[6][j]))
                if reverse==False:
                    pat_E[i]+=[-(F_mdcp_expt)]
                else :
                    pat_E[i]+=[F_mdcp_expt]
                continue
            else :
                #Now we also need to correct this values so that it represent X instead of GXG
                if new_seq[neigh]=='H':
                    He_seq=''
                    Hd_seq=''
                    for l in range(2*neigh+1):
                        if l==neigh:
                            He_seq+='8'
                            Hd_seq+='9'
                        else :
                            He_seq+=new_seq[l]
                            Hd_seq+=new_seq[l]
                    out_e=get_pattern_energy_from_database(He_seq,refs,reverse,penta,additive,T,base_rep,neigh,silent=True,reduced=reduced)
                    out_d=get_pattern_energy_from_database(Hd_seq,refs,reverse,penta,additive,T,base_rep,neigh,silent=True,reduced=reduced)

                    if out_e==None or out_d==None :
                        missing_pat=True
                        continue
                    else :
                        F_pat_e,F_pat_err,F_mdcp_e,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat_e,mdcp_offset_e,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt_e,T_pat=out_e
                        F_pat_d,F_pat_err,F_mdcp_d,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat_d,mdcp_offset_d,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt_d,T_pat=out_d

                else :
                    out=get_pattern_energy_from_database(new_seq,refs,reverse,penta,additive,T,base_rep,neigh,silent=True,reduced=reduced)

                    if out==None :
                        missing_pat=True
                        continue
                    else :
                        F_pat,F_pat_err,F_mdcp,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat,mdcp_offset,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt,T_pat=out

            if reverse==False:
                sign2=-1
            else :
                sign2=1
            if new_seq[neigh]=='H':
                pat_E[i]+=[get_G_sum(np.array([sign2*(F_mdcp_expt_e-sign*((F_pat_e-offset_pat_e)-(F_mdcp_e-mdcp_offset_e))),sign2*(F_mdcp_expt_d-sign*((F_pat_d-offset_pat_d)-(F_mdcp_e-mdcp_offset_e)))]),T)]
            else :
                pat_E[i]+=[sign2*(F_mdcp_expt-sign*((F_pat-offset_pat)-(F_mdcp-mdcp_offset)))]


        pat_ste[i]=np.array(pat_ste[i])

    if missing_pat==False:
        for i in range(len(pat_E)):
            if math.isnan(np.sum(pat_E[i])):
                missing_pat=True
                break

    if (missing_pat):
        print("Some pattern free energy are missing from the database.")
        quit()
    return pat_E,pat_ste,pat_seq


def get_pattern_additive_F(pattern,T,neigh,ent_corr=False):
    global refs
    global AA_type_QC
    global groups
    global R
    global DF_V,DF_Q,DU_V,DU_Q

    signs=refs[5]
    F_pat=0.

    for a0 in range(len(refs[0])):
        if pattern[neigh]==AA_type_QC[a0]:
            break

    #ARG is not part of this because no uncharged state, hence the DV already has the sign in it
    ions='EDH89KY'

    seq_F=[]
    seq_F_q=[]

    if ent_corr:
        DX_V=np.copy(DU_V)
        DX_Q=np.copy(DU_Q)
    else:
        DX_V=np.copy(DF_V)
        DX_Q=np.copy(DF_Q)

    #This is a correction to account for the context of the model compound
    for a1 in range(len(AA_type_QC)):
        if 'G'==AA_type_QC[a1]:
            break

    F_pat+=-DX_V[a0,a1,0]*2

    for i in range(len(pattern)):
        for a1 in range(len(AA_type_QC)):
            if AA_type_QC[a1]==pattern[i].upper() and pattern[i]!='_' and i!=neigh:
                F_pat+=DX_V[a0,a1,abs(i-neigh)-1]
                seq_F+=[DX_V[a0,a1,abs(i-neigh)-1]]
                #Here the five corresponds to the 6 ionizable amino acid
                if a1<=5 and pattern[i].isupper():
                    F_pat+=DX_Q[a0,a1,abs(i-neigh)-1]
                    seq_F_q+=[DX_Q[a0,a1,abs(i-neigh)-1]]

    if math.isnan(F_pat):
        print(DX_V)
        print(pattern," has a problem")
        input(F_pat)
    return F_pat

def check_sign(AA,refs,reverse):
    for j in range(len(refs[1])):
        if refs[1][j][1].upper()==AA.upper()  :
            break

    if j==len(refs[2]):
        return 0,0
    if reverse==True :
        sign2=-refs[5][j]
    else :
        sign2=refs[5][j]
    sign=refs[5][j]

    return sign,sign2
def get_pattern_energy_from_database(pattern,refs,reverse,penta,additive,T,data_dir,neigh,silent=True,reduced=True,
                                     unshifted=True,ent_corr=False):
    found=False
    for j in range(len(refs[1])):
        if refs[1][j][1].upper()==pattern[neigh].upper():
            found=True
            break

    if found==False:
        print(pattern[neigh])
        input("This should not happen")

    sign,sign2=check_sign(pattern[neigh],refs,reverse)
    #I will have to introduce temperature renormalaization
    F_mdcp_expt=R*T*np.log(10**(refs[6][j]))

    if additive:
        T_pat=T
        F_pat=get_pattern_additive_F(pattern,T,neigh,ent_corr=ent_corr)
        F_pat_err=0.
        F_mdcp=0.
        S_pat=0.
        S_mdcp=0.
        U_pat=F_pat
        U_mdcp=0.
        offset_pat=0.
        mdcp_offset=0.
        U_mdcp_expt=0.
        S_mdcp_expt=0.

    else:
        #This is if penta and additive is false (i.e. model compound values)
        T_pat=T
        F_pat=0.
        F_pat_err=0.
        F_mdcp=0.
        S_pat=0.
        S_mdcp=0.
        U_pat=0.
        U_mdcp=0.
        offset_pat=0.
        mdcp_offset=0.
        U_mdcp_expt=0.
        S_mdcp_expt=0.

    if F_pat_err>0.1:
        print("pentpeptide "+pattern+" has a large error ("+str(F_pat_err)+")")
    return F_pat,F_pat_err,F_mdcp,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat,mdcp_offset,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt,T_pat


def main_prediction(T,neigh,contexts,pat_E,pat_ste,pat_seq,map,prun_during,per_res,base_E,test_time,max_diff,silent=False,state_prun=False):
    ## Revision of 23/10/2025 summary
    ## Previous code had a mistake in the computation of non per_res energy,
    ## where the kept states would be included in the discarded energy, but not for the two end states
    ## still trying to address the issue, but it seems I will have to alwasy put all of the states inthe discarded
    ## states array, which will circunvent the issue.
    ## We should still be able to decompose into Kept vs Discarded , but by substracting instead of adding
    import itertools
    import sys
    if test_time :
        import time
        t0=time.time()
    blocks=[]
    tmp=[]
    for i in range(len(contexts)):
        tmp+=[i]
        if i%2==1:
            blocks+=[np.array(tmp)]
            tmp=[]

    #This is the total energy associated with all discarded states up until this point
    dim=(len(contexts)+1,)
    # dim base added 23/10/2025 : is the coordinates of the 11 state, which is the base state
    dim_base=(0,)
    for i in range(neigh):
        dim=dim+(2,)
        dim_base+=(0,)

    if not per_res:
        # This is the energy associated with all the states.
        lvl_tot_E=np.zeros(dim,dtype=np.longdouble)
        lvl_tot_E[:]=np.longdouble('+inf')
        # New 23/10/2025, setting the first energy to the base energy
        lvl_tot_E[dim_base]=base_E

    # This is a binary, because if no discarded states are in this array position, it should not be taken into account
    # These 2 arrays are for site specific information to be kept.
    # Ok so the strategy is to keep adding
    dim=(len(contexts)+1,len(contexts))
    for i in range(neigh):
        dim=dim+(2,)

    # This is new : so that we only compute the combinatorics once
    all_contexts=list(itertools.product([0,1],repeat=neigh))

    if per_res:
        lvl_context0=np.zeros(dim,dtype=np.longdouble)
        lvl_context0[:]=float('+inf')
        lvl_context1=np.zeros(dim,dtype=np.longdouble)
        lvl_context1[:]=float('+inf')
    elif state_prun :
        # These are the energy and identities of the non discraded states
        lvl_ste=[np.array(['']) for i in range(len(contexts)+1)]
        lvl_E=[np.array([],dtype=np.longdouble) for i in range(len(contexts)+1)]
        lvl_E[0]=np.append(lvl_E[0],[base_E])
        for q in range(1,len(lvl_E)):
            lvl_E[q]=np.append(lvl_E[q],[np.longdouble("+inf")])
    ############################################################################
    # Looping thought position top find the correct patterns
    # This is the calcualtion main loop (i.e. the step or t in the publication)
    for t in range(len(pat_E)):
        if not silent :
            tmp=t+1
            sys.stdout.write("\r Main prediction : step %4.f" %  tmp+" out of "+str(len(pat_E))+'\n')

        ############################################################################
        # Looping through already populated mesostates
        # Now we are going to find the relevant pattern, and save it in temp_pat_E
        dim=()
        for l in range(neigh):
            dim=dim+(2,)

        temp_pat_E=np.zeros(dim,dtype=np.longdouble)
        temp_pat_E[:]=float('+inf')
        ############################################################################
        # Now_finding all the relevant actual pattern in the simplified sequence space,
        # That is instead of a linear array, we have a ndimensional array, where n is the number of neighbor neigh
        for l in range(len(pat_ste[t])):
            all_ind=[()]
            ind=[]
            offset=0
            for o in range(1,neigh+1):
                if t-o>=0:
                    while neigh-o-offset>=0:
                        if map[t]-o-offset==map[t-o]:
                            ind=[neigh-o-offset]+ind
                            break
                        else:
                            offset+=1
            true_ind=[]
            for curr in  itertools.product([0,1],repeat=neigh-len(ind)):
                ind_temp=()
                for h in range(len(ind)):
                    ind_temp+=(int(pat_ste[t][l][ind[h]]),)
                true_ind+=[curr+ind_temp]
            for h in range(len(true_ind)):
                temp_pat_E[true_ind[h]]=pat_E[t][l]

        ############################################################################
        # Now populating the temp arrays for the "half steps" (G'  in the publication)
        if not per_res :
            tot_temp_E=np.zeros(np.shape(lvl_tot_E),dtype=np.longdouble)
            tot_temp_E[:,:]=float('+inf')

        if per_res:
            lvl_temp_context0=np.zeros(np.shape(lvl_context0),dtype=np.longdouble)
            lvl_temp_context0[:,:]=float('+inf')
            lvl_temp_context1=np.zeros(np.shape(lvl_context1),dtype=np.longdouble)
            lvl_temp_context1[:,:]=float('+inf')
        ############################################################################
        #Looping throught mesostates
        for q in range(t+1):
            # In this bit I add the energy to the previous energies while keeping context information.
            if per_res==True:
                if q==0:
                    curr1=()
                    for k in range(neigh):
                        curr1+=(1,)
                    curr0=curr1[:-1]+(0,)
                    # If this is the first time the residue is taken into account, prepopulate the free energy of the two first mesostates
                    # This part should be recoded I think, because we should be able to set the first free energy once and for all when the array is declared
                    lvl_context0[(q+1,t,)+curr0]=base_E+temp_pat_E[curr1]
                    lvl_context1[(q,t,)+curr1]=base_E

                # This is a loop through the residue taken into account so far
                for k in range(t):
                    #For every residue before i add the free energy of the context associated with i
                    base1=(q,k)
                    base0=(q+1,k)
                    for curr in all_contexts:
                        cont=base1+curr
                        pat_cont=curr

                        new_cont0=base0+curr[1:]+(0,)
                        new_cont1=base1+curr[1:]+(1,)

                        lvl_temp_context0[new_cont0]=get_G_sum(np.array([lvl_temp_context0[new_cont0],lvl_context0[cont]+temp_pat_E[pat_cont]]),T)
                        lvl_temp_context1[new_cont0]=get_G_sum(np.array([lvl_temp_context1[new_cont0],lvl_context1[cont]+temp_pat_E[pat_cont]]),T)
                        lvl_temp_context0[new_cont1]=get_G_sum(np.array([lvl_temp_context0[new_cont1],lvl_context0[cont]]),T)
                        lvl_temp_context1[new_cont1]=get_G_sum(np.array([lvl_temp_context1[new_cont1],lvl_context1[cont]]),T)

                k=t
                if q==0 :
                    k1=k
                    k2=k

                if q!=0 :
                    k1=k-1
                    k2=k

                base1k1=(q,k1)
                base1k2=(q,k2)
                #base0k1=(q+1,k1)
                base0k2=(q+1,k2)

                # Now for the levels that are not the first (q!=0) and the last residue taken into account,(t==k), the
                # free energy is inherited from the previous residue (eq 12abcd in the publication)
                for curr in all_contexts:
                    cont=base1k1+curr
                    pat_cont=curr
                    new_cont0=base0k2+curr[1:]+(0,)
                    new_cont1=base1k2+curr[1:]+(1,)
                    #This the free energy of residue
                    lvl_temp_context0[new_cont0]=get_G_sum(np.array([lvl_temp_context0[new_cont0],lvl_context0[cont]+temp_pat_E[pat_cont],lvl_context1[cont]+temp_pat_E[pat_cont]]),T)
                    lvl_temp_context1[new_cont1]=get_G_sum(np.array([lvl_temp_context1[new_cont1],lvl_context1[cont],lvl_context0[cont]]),T)

            if not per_res:
                base0=(q+1,)
                base1=(q,)
                for curr in all_contexts:
                    cont=base1+curr
                    pat_cont=curr
                    new_cont0=base0+curr[1:]+(0,)
                    new_cont1=base1+curr[1:]+(1,)
                    tot_temp_E[new_cont0]=get_G_sum(np.array([tot_temp_E[new_cont0],lvl_tot_E[cont]+temp_pat_E[pat_cont]]),T)
                    tot_temp_E[new_cont1]=get_G_sum(np.array([tot_temp_E[new_cont1],lvl_tot_E[cont]]),T)

        if per_res :
            lvl_context0=lvl_temp_context0
            lvl_context1=lvl_temp_context1
            continue
        else :
            lvl_tot_E=tot_temp_E

        if state_prun:
            lvl_E_tmp0=[[] for q in range(len(lvl_E))]
            #lvl_E_tmp1=[[] for q in range(len(lvl_E))]
            lvl_ste_tmp1=[[] for q in range(len(lvl_E))]
            lvl_ste_tmp0=[[] for q in range(len(lvl_E))]

            for q in range(t+1):
                # Loop in microstates
                # This part is only to find the relevant pattern and index it based the reduced sequence
                # to save it in to the kept states only !!!!
                # pat_E, pat_ste and pat_seq is the pattern in the non-reduced sequence, it will only contain one value if none its
                # neighbors are ionizable
                for m in range(len(lvl_ste[q])):
                    dims=()
                    for n in range(neigh):
                        if len(lvl_ste[q][m])-neigh+n<0:
                            dims+=(1,)
                        else :
                            dims+=(int(lvl_ste[q][m][len(lvl_ste[q][m])-neigh+n]),)

                    lvl_E_tmp0[q+1]+=[lvl_E[q][m]+temp_pat_E[dims]]
                    lvl_ste_tmp0[q+1]+=[lvl_ste[q][m]+'0']

                    #lvl_E_tmp1[q]+=[lvl_E[q][m]]
                    lvl_ste_tmp1[q]+=[lvl_ste[q][m]+'1']
            # Now recombine
            for q in range(t+1):
                lvl_E[q]=np.append(lvl_E[q],lvl_E_tmp0[q])
                lvl_ste[q]=np.append(lvl_ste_tmp1[q],lvl_ste_tmp0[q])
            lvl_E[q+1]=np.array(lvl_E_tmp0[q+1])
            lvl_ste[q+1]=np.array(lvl_ste_tmp0[q+1])

            # If you are going to keep everything anyway, skip the prunning
            if prun_during  and not max_diff==0:
                for q in range(t+1):
                    # If you are going to keep everything anyway, skip the prunning
                    if len(lvl_ste[q])<2 :
                        continue

                    # Skip the yet unpopulated levels, and do not do anything until at least two elements(avoids coding for exceptions)
                    if lvl_ste[q][0]=='':
                        break
                    ind=np.argsort(lvl_E[q])
                    if not (max_diff<=0):

                        # New code (28/10/2025)
                        G_tot=get_G_sum(lvl_tot_E[q].flatten(),T)
                        # print(lvl_E[q][ind])

                        DG_kept_vs_tot=np.zeros((len(lvl_E[q])))
                        for m in range(len(lvl_E[q])):
                            DG_kept_vs_tot[m]=G_tot-get_G_sum(lvl_E[q][ind[:m+1]],T)
                        temp=np.where(max_diff+DG_kept_vs_tot>0)
                        # print(DG_kept_vs_tot)
                        # print(temp)
                        if len(temp[0])==0:  #If no states get discarded because no combination of states is enough
                            temp=[[len(ind)]]
                        elif len(temp[0])==len(ind):#If all states get dicarded, keep the most stable one
                            temp=[[1]]

                        elif temp[0][0]<len(ind):# Finally, you can add one states (due to the where function, we want to include the one that makes the threshold be crossed)
                            temp=[[temp[0][0]+1]]

                    elif (math.isinf(max_diff) and max_diff>0):# No state gets discarded because the max diff is 0
                        temp=[[len(ind)]]
                    lvl_ste[q]=lvl_ste[q][ind[:temp[0][0]]]
                    lvl_E[q]=lvl_E[q][ind[:temp[0][0]]]

    if test_time :
        t1=time.time()
        print(t1-t0)
    if per_res:
        return lvl_context0,lvl_context1
    elif not state_prun:
        return None,None,lvl_tot_E
    else :
        return lvl_ste,lvl_E,lvl_tot_E

def find_ind(ind,arr) :
    for i in range(len(arr)):
        if arr[i]==ind :
            return True
    return False
def HSQ_internal_sequence_create(seq,list_res,no_write=False):
    #New version : passes the sequence instead of the name of the sequence file
    data=seq
    new_W=''
    pos_res=0
    neg_res=0
    base_charge=0  # This is to take into account the phosphoresidues -1 base state
    arg_res=0

    seq_data_q=[]
    seq_data_id=[]
    raw_seq=[]
    raw_seq_2=[]
    map=[]

    for i in range(len(data)):

        line=data[i]
        do_it=find_ind(i,list_res)
        raw_seq_2+=[line]

        if (line=="GLU" or line=="GLH" or line=="GLX") and do_it==True:
            new_W+='GLX\n'
            raw_seq+=['GLX']
            neg_res+=1
            seq_data_id+=[0]
            seq_data_q+=[-1]
        elif (line=="ASP" or line=="ASH" or line=="ASX") and do_it==True:
            new_W+='ASX\n'
            raw_seq+=['ASX']
            neg_res+=1
            seq_data_id+=[2]
            seq_data_q+=[-1]
        elif (line=="LYD" or line=="LYS" or line=="LYX") and do_it==True:
            new_W+='LYX\n'
            raw_seq+=['LYX']
            pos_res+=1
            seq_data_q+=[1]
            seq_data_id+=[1]
        elif line=="ARG":
            new_W+='ARG\n'
            raw_seq+=['ARG']
            pos_res+=1
            arg_res+=1
            seq_data_q+=[0]  # Why 0 ?
            seq_data_id+=[float('NaN')]
        elif (line=="TYX" or line=="TYR" or line=="TYO") and do_it==True:
            new_W+='TYX\n'
            raw_seq+=['TYX']
            neg_res+=1
            seq_data_q+=[-1]
            seq_data_id+=[3]
        elif (line=="SXP" or line=="S1P" or line=="S2P") and do_it==True:
            new_W+='SXP\n'
            raw_seq+=['SXP']
            neg_res+=1
            base_charge+=-1
            seq_data_q+=[-1]
            seq_data_id+=[4]
        elif (line=="TXP" or line=="T1P" or line=="T2P") and do_it==True:
            new_W+='TXP\n'
            raw_seq+=['TXP']
            neg_res+=1
            base_charge+=-1
            seq_data_q+=[-1]
            seq_data_id+=[5]
        elif (line=="YXP" or line=="Y1P" or line=="Y2P") and do_it==True:
            new_W+='YXP\n'
            raw_seq+=['YXP']
            neg_res+=1
            base_charge+=-1
            seq_data_q+=[-1]
            seq_data_id+=[6]
        elif (
                line=="HIE" or line=="HID" or line=="HIS" or line=="HEX" or line=="HDX" or line=="HIP" or line=="HIX") and do_it==True:
            new_W+='HIX\n'
            raw_seq+=['HIX']
            pos_res+=1
            seq_data_q+=[1]
            seq_data_id+=[9]
        elif do_it==False and (line=="HIS" or line=="LYS"):

            new_W+=line+'\n'
            seq_data_q+=[0]
            seq_data_id+=[float('NaN')]
            raw_seq+=[line]
            base_charge+=1
        elif do_it==False and (line=="ASP" or line=="GLU" or line=="TYO"):
            new_W+=line+'\n'
            seq_data_q+=[0]
            seq_data_id+=[float('NaN')]
            base_charge+=-1
            raw_seq+=[line]
        elif line!="END":
            if line!="NA+" and line!="CL-" and line!="CLX" and line!="NAX":
                raw_seq+=[line]
            new_W+=line+'\n'
            seq_data_q+=[0]
            seq_data_id+=[float('NaN')]  #, for now we will ignore the context, so that is ok
        if seq_data_q[-1]!=0 and line!="END":
            map+=[i]

    seq_id_reduced=[]
    for i in range(len(seq_data_id)):
        if math.isnan(seq_data_id[i])==False:
            seq_id_reduced+=[seq_data_id[i]]

    return pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W


def read_sequence(seq_3,list_res):
    #This has to be caled twice : this is the first call
    pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W=(HSQ_internal_sequence_create(seq_3,list_res))
    titrable_residue_indexes=[]
    if titrable_residue_indexes==[]:
        list_res=np.array([i for i in range(len(seq_3))])
    else:
        list_res=np.array([titrable_residue_indexes[i] for i in range(len(titrable_residue_indexes))])

    if titrable_residue_indexes==[]:
        ind=0
        sites_num=neg_res+pos_res-arg_res
        titrable_residue_indexes=np.zeros((sites_num),dtype=int)
        for i in range(len(seq_data_q)):
            if seq_data_q[i]!=0:
                titrable_residue_indexes[ind]=i
                ind+=1
    else:
        sites_num=len(titrable_residue_indexes)
        titrable_residue_indexes=np.array(titrable_residue_indexes)

    pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W=HSQ_internal_sequence_create(seq_3,list_res)

    return sites_num,titrable_residue_indexes,pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W


def MEDOC_wrapper(seq_1,T=293,neigh=4,unshifted=False,additive=True,base_rep='./.'):
    # This is a wrapper to the MEDOC algorithm, for large data set deployement within python
    # It does not support per res, nor does it print any figure
    # if (not unshifted and not additive) or (unshifted and additive) :
    reduced=False
    reverse=True
    penta=False
    test_time=False
    prun_during=False
    state_prun=False
    per_res=False
    base_E=-709.0*(R*T)
    max_diff=0.
    if seq_1[0]!='z':
        seq_1='z'+seq_1
    if seq_1[-1]!='z':
        seq_1=seq_1+'z'
    seq_3=convert_AA_1_to_3(seq_1,mode=2)
    list_res=np.array([i for i in range(len(seq_3))])
    seq_1=convert_AA_3_to_1_letter(seq_3)

    sites_num,titrable_residue_indexes,pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W=read_sequence(seq_3,list_res)
    Nmes=pos_res+base_charge+1-(-neg_res+arg_res+base_charge)
    #layers_q is inverted to match the order of more to less protons used in the free energy (new in 1.4.5 where the dicrepancy was corrected.)
    layers_q=np.array([i for i in range(-neg_res+arg_res+base_charge,pos_res+base_charge+1)])[::-1]
    contexts,states=get_base_contexts(map,neigh,seq_1,seq_data_q)
    pat_E,pat_ste,pat_seq=get_all_contexts(contexts,states,neigh,refs,unshifted,reverse,penta,additive,T,base_rep,
                                              reduced)
    lvl_ste,lvl_E,lvl_tot_E=main_prediction(T,neigh,contexts,pat_E,pat_ste,pat_seq,map,prun_during,per_res,base_E,test_time,max_diff,state_prun=state_prun,silent=True)
    Meso_G_all=np.zeros((Nmes),dtype=np.longdouble)
    for i in range(Nmes):
        Meso_G_all[i]=get_G_sum(lvl_tot_E[i],T)

    return layers_q,Meso_G_all

def get_q_profile_from_F(pH,Fs,q,T):
    p=get_p_profile_from_F(pH,Fs,q,T)
    q_out=np.zeros((len(pH)))
    for i in range(len(Fs)):
        q_out[:]+=p[i,:]*q[i]
    return q_out

def get_p_profile_from_F(pH,Fs,q,T):
    W=np.ndarray((len(Fs),len(pH)),dtype=np.longdouble)  #dtype=float)#dtype=np.longdouble)
    # OK : For each pH values I will introduced a offset
    # that will correspond to
    exponent=np.zeros((len(Fs),len(pH)))
    for i in range(len(Fs)):
        # Think forthe quenched sequences I will have to get this working :
        c=max(q)-q[i]
        exponent[i,:]=-(Fs[i]-np.log(10.)*R*T*pH[:]*c)/(R*T)

    for p in range(len(pH)):
        W[:,p]=np.exp(exponent[:,p]-np.amax(exponent[:,p]))

    p=np.ndarray((len(Fs),len(pH)),dtype=float)
    p[:,:]=W[:,:]/np.sum(W,0)[:]
    return p

def  initialize_additive_DF_array_MEDOC_public():
    import numpy as np
    global DF_V,DF_Q,R,AA_type_QC,AA_sign_QC

    AA_type_QC=['E','D','H','K','Y','8','9','R','C','A','F','G','L','I','M','N','P','Q','T','V','S','W']
    AA_sign_QC=[-1,-1,1,1,-1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]


    R=np.longdouble(1.98720425864083*10**(-3))

    DF_V=np.array([[[-0.09814503006414027,-0.00829083674554501,0.02069052993698899,0.022764800034716755,0.026376024825803358,0.04079061538041999]\
    ,[-0.09304533546944738,-0.01168163063259197,0.029148162211159886,0.034573952585421716,0.036962516294832264,0.04389309879451119]\
    ,[-0.10917635227208095,-0.01508296471050332,0.021012423403039278,0.033690195118149475,0.04356625287618849,0.029829597140747288]\
    ,[-0.08448688459551854,0.0010452734081126264,0.013127658435950915,0.012953705538618056,-0.029404194130851484,0.055030180399804576]\
    ,[-0.21690303214175155,-0.0048112297775041585,0.05285639369305626,0.03221913291471838,0.05672514238466199,0.107241996117071]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[0.0960264723482048,0.14624034591723908,0.18347899068707088,0.18975875712655169,0.14936689949676674,0.09615404956002349]\
    ,[-0.06467454370648426,0.013380488692223039,0.04335885350741095,0.02316539326397672,0.049903954408949965,0.05498538753463285]\
    ,[-0.05512108276442692,0.01663051907476934,0.032942793205950234,0.029777445424963897,0.042157426628694025,0.042587166017544566]\
    ,[-0.11832732948563851,-0.00594652189278002,0.05161987726427455,-0.0032819449596649397,0.02802986568478801,-0.05330800036191256]\
    ,[-0.017883702304757478,0.02456300301888289,0.03806376083390284,0.031143225951056598,0.04316622527388991,0.03861323108798072]\
    ,[-0.0874480230112806,0.0005622223046289621,0.01942881006255066,-0.014285434717318975,-0.0312269154603095,0.06646508312597733]\
    ,[-0.09989478442823643,-0.01808123925823244,0.017072070544487825,0.0348368909924136,-0.03310474356655035,0.09507777627859972]\
    ,[-0.08158971639604323,0.003750311972060806,0.021579582095461788,0.02995580604310409,0.00500911632668722,-0.01839556768577185]\
    ,[-0.08951352186020228,-0.0006698252487210704,0.024195553599692513,0.03159319892243239,0.030212308735902412,0.05633540193579425]\
    ,[-0.236836729200542,-0.017997997084072716,-0.01553659428180449,0.030141611227043254,0.042040512786275526,0.056201476460303386]\
    ,[-0.09242179483084709,-0.013124673800463844,0.01576136836717485,0.02174268842319781,0.03447444655875749,0.009438608909709427]\
    ,[-0.07935310042767338,-0.0010719297463103114,0.03344422486666974,0.030610203233417295,0.03890124193586497,0.04473450961219479]\
    ,[-0.08695017047325584,-0.005373952979996946,0.029127857638777982,0.0312603918689327,0.04934931757260941,0.0368315168559229]\
    ,[-0.05628676396346869,0.018962415116782136,0.030200837655631817,0.04260530926805834,0.04881322480499954,0.04544388691431143]\
    ,[-0.1579309119001219,-0.02341443756830642,0.037225513676135386,0.0362792660383576,0.004737830912580683,0.04319702479472116]]\
    ,[[-0.1736875245498383,-0.03212340028874131,0.0018993119947976443,0.0059849245089558765,0.02124825095580544,0.01542135737912272]\
    ,[-0.1710557565478377,-0.03364870882799501,0.006962473577143866,0.014208279143679049,0.019176229117946218,0.018666569601631332]\
    ,[-0.18704258417458502,-0.03620398638246763,0.0021258605695279216,0.00958528489738269,0.015950595123349547,0.015888493456329162]\
    ,[-0.15978519863964097,-0.015193111655213011,-0.004543537228943732,-0.005810504902595361,0.001211831711295016,-0.029637568260411744]\
    ,[-0.265794333809359,-0.01098275661415156,0.023310980280308924,0.04604403647604573,-0.004383328252121028,0.15281322208498269]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[0.028327524005790344,0.12678580626802685,0.15053243868145888,0.1396583921394307,0.10704990032049617,0.08041056791800583]\
    ,[-0.14040388012430757,-0.004634556521311273,0.021221304206803588,0.025415163643996482,0.01300114155781146,0.035960326594166314]\
    ,[-0.11180005121071279,-0.00790661654901701,0.008065269220748556,0.01915074666059384,0.013501678469567737,0.02865841453001855]\
    ,[-0.18403285574313183,-0.015268232708448258,0.01906798353026109,0.01277357309186815,-0.005488314153653567,0.06417328536697009]\
    ,[-0.05591611979451071,0.002490637524058679,0.017408220939823003,0.017456553476757463,0.019898598245154074,0.016195291585171256]\
    ,[-0.17062346145813187,-0.018482405659511855,-0.006386155406357567,-0.0059814323124224265,0.030944887598309455,-0.08044967502170272]\
    ,[-0.1965033680490634,-0.044031667092552115,0.007872091844576625,0.012051754836853036,0.018100660878122947,0.014100063508717125]\
    ,[-0.15864127590529692,-0.012387549958144464,0.0074724624621604445,0.009305376582755747,0.011997383040821622,0.019479288206484018]\
    ,[-0.15805865104866434,-0.025678890170240803,0.0056716629622748475,0.011313294883195971,0.0231028999052599,0.025783803622022983]\
    ,[-0.3365687377038947,-0.03774988730777259,-0.025653389166144812,0.0422404998983396,-0.024724320160489067,-0.042838831623705706]\
    ,[-0.16928157649557038,-0.03166840373673335,-0.0012311576526688114,0.009746434166865894,0.006601513489848894,0.02649554221267829]\
    ,[-0.1595613855219334,-0.022950366888830474,0.010269451758506987,0.017551839786707345,0.027947727012862557,0.015999523029874063]\
    ,[-0.17522506486364278,-0.030856192130948194,0.014790293555891787,0.007538631121956855,0.015682699720808956,0.031195806945032087]\
    ,[-0.11783122007918154,-0.007190406917247816,0.013019100625245444,0.020529828574650338,0.022204857375337155,0.02236469649994579]\
    ,[-0.21939902018557844,-0.0236924257889203,0.018484603045736094,-0.005867413603460225,0.014840957243556803,0.0020237262457752613]]\
    ,[[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]\
    ,[[-0.0006560363800428431,0.017576246527726744,0.03008164293748003,0.03415531927321283,0.030141568877442886,0.044207304828451685]\
    ,[0.006982542206302393,0.025014332773102395,0.032886085087876785,0.03572446532684106,0.048716296570434164,0.033745325739614765]\
    ,[-0.007518042027014102,0.023770881049359188,0.03101304534592246,0.036380156633212375,0.043646569226534766,0.04024154617270405]\
    ,[-0.01007046276208682,0.012742062309855771,0.023814123942612335,-0.0021809941525407495,0.010467561390157747,-0.06573477481798018]\
    ,[0.055477009911087984,0.00627945480423403,0.0092874033522772,-0.038247383584260476,0.007466647906540322,0.04268048266274342]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[-0.1321694626787929,-0.0967316389078674,-0.06873179885546174,-0.016529858238667795,0.019671997748505726,0.029419251146583184]\
    ,[-0.010364462267885927,0.013790474140330963,0.026345199526713606,0.022818234394044866,0.025409162636675052,0.012760754702843317]\
    ,[0.005442484530153016,0.025787234454961675,0.03376592092627775,0.032877785177585404,0.031348968336479405,0.043241683575552595]\
    ,[-0.02251451717471522,0.007290575813393679,0.029400593000571247,0.029340306325182742,0.03981605687929736,0.01750648773562922]\
    ,[0.0166600334696883,0.03005643512105581,0.03641745659081505,0.03817866564826463,0.042348875288438774,0.041078310001161036]\
    ,[-0.00597323606546267,0.015142925079749282,0.022055801572964103,0.007399495395348531,-0.021541172512766056,0.07552783558902708]\
    ,[-0.008063015206341818,0.015665392013667137,0.022606845297376915,0.017825963721647624,0.06940144304651008,0.009221847125200112]\
    ,[-0.013793695757994474,0.009711396978727332,0.020258871518649742,-0.012595591191361784,0.024467224765768833,0.06333164148878584]\
    ,[0.006810928658946724,0.02783802048600452,0.033148777963438444,0.0351957271708956,0.03950335340840296,0.03897919554488183]\
    ,[-0.017184653241082218,0.05871433408116973,0.030341030778760117,0.05681139880223231,0.06634117688305469,0.09662211105241225]\
    ,[0.003242459770048896,0.022390925906296884,0.03127235396341395,0.03149620373316191,0.04106117199247397,0.030986738230998078]\
    ,[-0.0029020200448697936,0.01895285746039096,0.032008369025381,0.02894746133170481,0.03739628818550957,0.042540476297507034]\
    ,[-0.0033934732964372646,0.018971356384024366,0.0292200144718752,0.027318158555992473,0.03527353223490641,0.047494930505991816]\
    ,[-0.0020452791650115223,0.022734026778768082,0.03245626597561265,0.03171437996861728,0.037548846694976,0.04066331047112597]\
    ,[-0.029279333847121773,0.0075547513069324015,0.03964118638194607,0.031266639383598735,0.011548266870763664,0.08870967327908091]]\
    ,[[-0.17084210594904514,-0.022099918877378576,0.014120014237955923,0.01255769428691967,0.026583953212703074,0.01508086914077416]\
    ,[-0.15069516354151227,-0.017271791267881753,0.013818742554240177,0.025776460947306273,0.03734918248509099,0.03749059989028692]\
    ,[-0.16953746812172468,-0.02493478943963027,0.008277129893803796,0.022277924856971006,0.030544694025584632,0.02482303514003073]\
    ,[-0.171039209692462,-0.01946687683707459,-0.0074749602680714775,-0.023348622356317343,-0.05257503495583072,0.037288310594207286]\
    ,[-0.21650285796168306,-0.011581570239445858,0.033906637707397144,-0.005483527606754593,-0.050016790733475903,0.10194235048455046]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[-0.07662516120041804,0.13940191874848715,0.21611212509088626,0.30629152638941,0.23390459740311087,0.12167267060185707]\
    ,[-0.14018106245921416,0.0018208238718888312,0.01530542767703933,0.025531693374791897,0.02918652076400237,0.033180993443134514]\
    ,[-0.1030398571705276,0.00826837665424289,0.017135425858650803,0.02848090435760641,0.03280738721302972,0.02767118822745839]\
    ,[-0.18277744826286102,-0.022016125001255972,0.019320434845680724,0.022704868380946538,0.00352130636308514,-0.04920120272901054]\
    ,[-0.05074551177888,0.02043138928037861,0.02999078335203812,0.031082277613801064,0.037303725043066664,0.03535185363560578]\
    ,[-0.17133841843828002,-0.020458348940630693,-0.016589408014715888,-0.056413389743416066,-0.03451895611115501,-0.09671283790963942]\
    ,[-0.17791201755244088,-0.03257102792905086,-0.007922328089363398,0.021382457041577985,0.006235694645523211,0.07206103259101587]\
    ,[-0.17202345863873184,-0.017508965353317633,0.004866383509992529,-0.04401083967617992,0.004158337780145391,-0.12476469400782567]\
    ,[-0.14596461828626855,-0.014392467967268535,0.011380757575454856,0.0259447173353388,0.040060425921753226,0.02289681007288067]\
    ,[-0.2491078747670432,-0.03351869936540874,0.031363147841980424,0.02747909500284543,0.010956243002188969,0.05665228202486094]\
    ,[-0.16640801557942994,-0.027874538747053403,0.00043379734266717207,0.00399293616568637,0.01652422827784284,0.016853112800304597]\
    ,[-0.14642355668248966,-0.006516623643763118,0.017823607455771964,0.02637505361495418,0.031168036345400215,0.05807761918516998]\
    ,[-0.15553111788032722,-0.014889632976041199,0.018827175251755857,0.017106832380835508,-0.0003291527991291718,0.04658364918622336]\
    ,[-0.11400032152104517,0.006110822465742141,0.02408009581539319,0.030118137588815312,0.0384340040672249,0.042394350383246564]\
    ,[-0.2218388958182611,-0.030693963794283748,0.027490772302327198,-0.00823444092883456,0.003123786489651184,0.006554724627917757]]\
    ,[[-0.06954768366510868,0.030177504977789076,0.054770271364577294,0.04261030636559745,0.03775383628891438,0.05575076563937954]\
    ,[-0.055708066148441454,0.03725641633259451,0.05277781767756678,0.050113796544419206,0.04306540454092449,0.05956084128910389]\
    ,[-0.0787986979727334,0.03619307215733659,0.05236189562250809,0.05550988814525403,0.04176739806355216,0.0677442994364343]\
    ,[-0.0853426228335451,0.0280769108239566,0.048490052459715355,0.03894908998485083,0.052477352528545365,0.05047300245967379]\
    ,[0.04700648739448552,0.03935459306590028,0.047433501968440456,-0.02263925513439806,0.1191425610604626,0.22011421957512442]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[-0.3306146077056944,-0.18463624664290362,-0.12288530646194068,-0.04584689421456875,0.01855413513242482,0.030648706657669222]\
    ,[-0.09607090294139137,0.011935001221971518,0.03420555166470182,0.030355550886832124,0.042193093069604816,0.027961177914468144]\
    ,[-0.05019612964583711,0.03714754412064866,0.047726615121928224,0.041694545552436485,0.04913910317498142,0.0562924610341742]\
    ,[-0.09480370450444114,0.01879951232631937,0.04740579291844597,0.03390329629944078,0.10970039412260016,0.018237465680700865]\
    ,[-0.020860277017295384,0.042094489883172985,0.04364132984350217,0.0514679791048815,0.04804117451776813,0.05012549631127791]\
    ,[-0.07881329188054159,0.03058398532969841,0.052413373611755236,0.029208470689335093,0.013611902688318592,0.12047661263327764]\
    ,[-0.08839035468092035,0.028389191034801167,0.0513679896103843,0.01955461662590141,0.07450388335634055,0.04846708006735181]\
    ,[-0.09468332464926194,0.014729385661499006,0.04318280523427651,0.024583592299987182,0.07959910619937274,-0.03495810941190802]\
    ,[-0.04778421234029004,0.043440230076270285,0.04848710698126801,0.052152826234281725,0.05307107483319181,0.052211872090265524]\
    ,[-0.15953590511358468,0.09241466070454163,0.06451155147331174,0.10173037154073514,0.05604151973154802,0.076883009229945]\
    ,[-0.05857973226001911,0.04035072207494691,0.056966234697544714,0.043926088288897595,0.05460213058708099,0.05397138818586648]\
    ,[-0.08357920342077421,0.03092900602301457,0.04506631759431765,0.038634990717413116,0.047629581221684994,0.05477173150265985]\
    ,[-0.07947823192886266,0.030155277561817552,0.049350445192060034,0.03843069643496533,0.0523094672125382,0.025845428311250678]\
    ,[-0.07072944159839083,0.027201651739164907,0.04335906829644441,0.04064806080895663,0.04592628726211476,0.050792661897275916]\
    ,[-0.08329137203581599,0.026398779598624196,0.06347964745866509,0.048704222857452545,0.06842390576347482,0.10480425478276015]]\
    ,[[-0.06241764196708248,0.030258851929258535,0.058497951759854436,0.03543142885149094,0.055694904685946595,0.04973755035173888]\
    ,[-0.049871899316551074,0.03868122766364475,0.051551739673094725,0.04719869900017683,0.046431426368617086,0.061658814055687225]\
    ,[-0.06767352542436204,0.034069287121829196,0.05312766600416529,0.052570228616864556,0.04009483359560788,0.06964873858943778]\
    ,[-0.0759514469708133,0.027821978509968766,0.0594104641746174,0.03170396614939606,0.05504977666724973,0.10995166859727987]\
    ,[0.06306411997928503,0.056981090378274035,0.041489788159023225,0.259430698425717,0.09158564296497097,-0.016283221939314474]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
    ,[-0.3235061069678418,-0.1842539289857012,-0.12280615414625906,-0.044429593884948895,0.020936026179415687,0.03735161906539631]\
    ,[-0.08816306817279224,0.013827610542929851,0.03579976556657624,0.03628921986673164,0.03679358183660416,0.036376937040425036]\
    ,[-0.044018773965817395,0.03777359820865288,0.05337695426745953,0.042406261914143,0.0512259910133126,0.046167730710859736]\
    ,[-0.08717816826930877,0.02159139980237254,0.04266751465389509,0.05244245660732449,0.05523080239506714,0.05650872536106053]\
    ,[-0.015754064383376398,0.040656052258744614,0.04399057020102738,0.04714714651079062,0.04869226120207992,0.0548035218448261]\
    ,[-0.07032816789656965,0.03364919208313946,0.05938864713476791,0.04761186182260165,0.08641261778709591,0.1269246178151925]\
    ,[-0.07713146002135796,0.028857066131721832,0.055891628601149,0.033906344880834924,0.03614197037450594,0.11425205461784368]\
    ,[-0.08534918855170714,0.018139916957455128,0.0374870332875474,0.034780312129211184,0.07805103236546099,-0.01144277270919658]\
    ,[-0.04048455964478685,0.04434063707057383,0.050034275300582506,0.04603917028678425,0.0480787709732764,0.06796938365975187]\
    ,[-0.1668632793946712,0.07451305551296952,0.06455361020041653,0.09803899503255731,0.0875010043260986,0.05297868326911327]\
    ,[-0.048740498375005825,0.04209079981959804,0.06011331485568622,0.0508357994987634,0.04907782823721556,0.061750730256215865]\
    ,[-0.07068978756286647,0.02493631564298981,0.04914099994522264,0.03490048880350949,0.04771410986344599,0.04733895305420864]\
    ,[-0.07167569778366181,0.03181337229462032,0.051764565467234604,0.05255319173179999,0.050930784654692265,0.049469690472103346]\
    ,[-0.06220262655911074,0.026375444158202355,0.045660392521118004,0.03675975472204902,0.04715801231830443,0.049098366874776186]\
    ,[-0.07447997933624267,0.02916747879908975,0.06389228591248804,0.02170970257759925,0.1308768760056414,-0.022224698687782096]]],dtype=np.longdouble)

    DF_Q=np.array([[[-0.26912017526287124,-0.19790042588058968,-0.1771325119001334,-0.1390228005347286,-0.1157313592250951,-0.10623193041288091]\
         ,[-0.3560263520607516,-0.24627933025003113,-0.20934366100910642,-0.1830222389036811,-0.14392042913809658,-0.138717075048026]\
         ,[0.4256475738644047,0.372167923185914,0.252328295694633,0.2013536613799452,0.16072755318784787,0.15249780400932217]\
         ,[0.22843715557208621,0.2269313528526908,0.17840996635791126,0.1394174304403735,0.10357655822167937,0.09321639483857544]\
         ,[-0.09475537748184883,-0.19557242312043074,-0.22431237619319067,-0.17131551458749222,-0.16931225672131353,-0.1947467014079522]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]\
        ,[[-0.13248987068156245,-0.10049503323340414,-0.08685052124104743,-0.0701251457833387,-0.05702023922070306,-0.050942501573915534]\
         ,[-0.17017971064212944,-0.1253152356448843,-0.10228366950929935,-0.08397203855373263,-0.07077990144669889,-0.058161752283569616]\
         ,[0.20775097292477857,0.1654398880690938,0.1109290835441011,0.08640010633109754,0.07202897414994697,0.06512582577626504]\
         ,[0.10843464777552216,0.09739262094715187,0.07512217833872348,0.06024057868295137,0.054960462708227624,
           0.049159700886994465]\
         ,[-0.0540502395024895,-0.1084284349040225,-0.10078772557805978,-0.08528817049496785,-0.06503185730851604,
           -0.07861169364041086]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]\
        ,[[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]\
        ,[[0.17121877608402247,0.17594084065247703,0.12904059170967974,0.10599129686565438,0.0943658131688274,0.08009514546801189]\
         ,[0.22150694409377722,0.21535255584026441,0.159419594251915,0.12981785696298773,0.1091688918983498,0.09918833660533939]\
         ,[-0.18949274939405808,-0.14328125208896705,-0.13308323983500342,-0.10907394655331742,-0.09240760246281003,-0.08130748156420634]\
         ,[-0.12189804277217274,-0.08684060560503566,-0.08147177646216318,-0.061704775073613716,-0.05818944513703821,-0.03908642889654734]\
         ,[0.0901457700763897,0.18835106928036163,0.14681904386655734,0.14354934675687644,0.1362350484771623,0.17465564751635806]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]\
        ,[[-0.07645991225335493,-0.12870699114633158,-0.134525810533965,-0.10662665665725732,-0.08830409119430149,-0.07816545475323736]\
         ,[-0.10317973252714768,-0.16529233106043975,-0.1458598903772544,-0.12974338865199825,-0.11340577364567742,-0.10199699387084803]\
         ,[0.17147809499398675,0.28614938533479006,0.2024521716185936,0.1628948550877939,0.1379793587864573,0.12891520715142293]\
         ,[0.09516414231334615,0.1802548822952621,0.13379804760658265,0.13001112605475074,0.1278112177766325,0.10176741009077647]\
         ,[-0.0063066501821508305,-0.11717161716288614,-0.15125766498204232,-0.11602402373209846,-0.10881237599290189,-0.08702408889399628]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]\
        ,[[0.324381869187757,0.2896563269057624,0.21037396198031996,0.16785162077403384,0.15061781758824555,0.12587627126945658]\
         ,[0.4436283912859322,0.37054832958760164,0.2602175109665942,0.2076770991050208,0.17872218702103043,0.15183733805515598]\
         ,[-0.39755021304844146,-0.2791305921075633,-0.24608355614259247,-0.19942038873986354,-0.16606973931225189,-0.14613821155689652]\
         ,[-0.23051190847571815,-0.1665060086739995,-0.15035870954979552,-0.12817908241645828,-0.11366588399840559,-0.08828327497847531]\
         ,[0.1617488660702995,0.29542855228415854,0.23544614440777636,0.22903045500858785,0.23155974983657546,0.4235103664889902]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]\
        ,[[0.3293775540187506,0.2913838324984375,0.21159625457985348,0.17397836516386642,0.14959887418322163,0.13085586119282838]\
         ,[0.4546388456822364,0.3717559299691354,0.2684538320174004,0.21535473636227676,0.17999577597819516,0.15832215012571763]\
         ,[-0.4007088436796262,-0.28065196046897967,-0.25069862782592517,-0.20252049317423879,-0.1622791413840075,-0.14624337458185724]\
         ,[-0.23337398224807474,-0.16983815270401478,-0.15930132103271072,-0.12233479968121273,-0.11269969693784955,-0.10376378536059733]\
         ,[0.15723558477489155,0.29073213834403167,0.23113927194876213,0.19725037930799214,0.22571642137322856,0.2932908418138935]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]\
         ,[float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]]],dtype=np.longdouble)

