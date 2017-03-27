import numpy as np
    
def GetResults_sub(sub, ind, H, H1, H0, M, S, NT):
    '''Get QALY and events results for a particular subset of the population'''
    IQALY = H1.QALY[sub] - H.QALY[sub]
    ILY = H1.LY[sub] - H.LY[sub]
   
    days = 365.25 
    M[ind,0] = (H.QALY[sub]).mean() * days
    M[ind,1] = (H1.QALY[sub]).mean() * days 
    M[ind,2] = IQALY.mean() * days 
        
    M[ind,3] = (H.LY[sub]).mean() * days 
    M[ind,4] = (H1.LY[sub]).mean() * days 
    M[ind,5] = ILY.mean() * days 

    S[ind,0] = (H.QALY[sub]).std() * days 
    S[ind,1] = (H1.QALY[sub]).std() * days 
    S[ind,2] = IQALY.std() * days 
    S[ind,3] = (H.LY[sub]).std() * days 
    S[ind,4] = (H1.LY[sub]).std() * days 
    S[ind,5] = ILY.std() * days 
    
    ## QALY and LY gains per health check
    nhcbase = H.Attending.sum()
    nhc = H1.Attending.sum()
    IQALYscen_no = H1.QALY[sub] - H0.QALY[sub] # scenario vs no HC 
    IQALYbase_no = H.QALY[sub] - H0.QALY[sub]  # base case vs no HC
    miq_scenno = IQALYscen_no.sum() / nhc
    miq_baseno = 0 if (nhcbase==0) else IQALYbase_no.sum() / nhcbase
    M[ind,6] = (miq_scenno - miq_baseno)*days 
    ILYscen_no = H1.LY[sub] - H0.LY[sub] 
    ILYbase_no = H.LY[sub] - H0.LY[sub]  
    mil_scenno = ILYscen_no.sum() / nhc
    mil_baseno = 0 if (nhcbase==0) else ILYbase_no.sum() / nhcbase
    M[ind,7] = (mil_scenno - mil_baseno)*days
    vscen = ((IQALYscen_no - miq_scenno)**2).sum() / nhc
    vbase = 0 if (nhcbase==0) else ((IQALYbase_no - miq_baseno)**2).sum() / nhcbase
    varM6 = vscen + vbase
    S[ind,6] = np.sqrt(varM6)*days
    vscen = ((ILYscen_no - miq_scenno)**2).sum() / nhc
    vbase = 0 if (nhcbase==0) else ((ILYbase_no - miq_baseno)**2).sum() / nhcbase
    varM7 = vscen + vbase
    S[ind,7] = np.sqrt(varM7)*days
         
#    IQALY.sum() / nhc * days 
#    ILY.sum() / nhc * days
#    S[ind,6] = np.sqrt( ((IQALY - M[ind,6])*(IQALY - M[ind,6])).sum() / nhc ) * days 
#    S[ind,7] = np.sqrt( ((ILY - M[ind,7])*(ILY - M[ind,7])).sum() / nhc ) * days 

    H.dead = np.logical_not(H.alive)
    H1.dead = np.logical_not(H1.alive)    

    k = 8
    ## Counts of events by specific ages 
    outcomes = ['IHD', 'Stroke', 'Dementia', 'LungCancer']
    ages = [80]
    for i in ages:
        for j in outcomes:
            exec('D = (H.%s * (H.age == %s)).sum(axis=1)' % (j, i))
            exec('D1 = (H1.%s * (H1.age == %s)).sum(axis=1)' % (j, i))
            M[ind,k] = (D[sub]).mean()
            M[ind,k+1] = (D1[sub]).mean()
            M[ind,k+2] = M[ind,k] - M[ind,k+1]
            S[ind,k] = (D[sub]).std()
            S[ind,k+1] = (D1[sub]).std()
            S[ind,k+2] = (D[sub] - D1[sub]).std()
            k = k+3

    ## event counts by end of the model, i.e. age 99. Everyone dies by age 100.
    prevage = np.column_stack((H.age[:,0]-1, H.age[:,:-1]))
    first99 = (H.age == 99) * (prevage == 98) # get around ages after 99 being labelled 99, don't understand why this was done
    for j in outcomes:
        exec('D = (H.%s * first99).sum(axis=1)' % (j))
        exec('D1 = (H1.%s * first99).sum(axis=1)' % (j))
        M[ind,k] = (D[sub]).mean()
        M[ind,k+1] = (D1[sub]).mean()
        M[ind,k+2] = M[ind,k] - M[ind,k+1]
        S[ind,k] = (D[sub]).std()
        S[ind,k+1] = (D1[sub]).std()
        S[ind,k+2] = (D[sub] - D1[sub]).std()
        k = k+3

    outcomes = ['dead']
    ages = [75, 80]
    for i in ages:
        for j in outcomes:
            exec('D = (H.%s * (H.age == %s)).sum(axis=1)' % (j, i))
            exec('D1 = (H1.%s * (H1.age == %s)).sum(axis=1)' % (j, i))
            M[ind,k] = (D[sub]).mean()
            M[ind,k+1] = (D1[sub]).mean()
            M[ind,k+2] = M[ind,k] - M[ind,k+1]
            S[ind,k] = (D[sub]).std()
            S[ind,k+1] = (D1[sub]).std()
            S[ind,k+2] = (D[sub] - D1[sub]).std()
            k = k+3
                    
    if (ind>=4):
        NT[ind-4] = sub.sum()
    return M,S,NT
    

def GetResults_longterm(H, H1, H0, M, S, N):
    '''Define various subpopulations of interest and get their expected sizes, QALY and events results'''
    allpop = np.ones(H1.population_size, dtype=bool)
    ## eligible, offered, or attend at least once 
    el = (H1.eligible.sum(axis=1)>0) 
    hc_offered = (H1.OfferedHC.sum(axis=1)>0)
    att = (H1.Attending.sum(axis=1)>0)
    t_offered = (H1.OfferedTreatment.sum(axis=1)>0)

    ## Eligible for HC and for each kind of treatment at any time, using either Q>=20 or Q>=10
    el_trtq20 = ((H1.eligible * (H1.QRisk >= 20)).sum(axis=1) > 0)
    el_trtq10 = ((H1.eligible * (H1.QRisk >= 10)).sum(axis=1) > 0)
    el_bp = ((H1.eligible * (H1.q_sbp > 140)).sum(axis=1) > 0)
    el_bmi = ((H1.eligible * (H1.bmi >= 30)).sum(axis=1) > 0) 
    el_sm = ((H1.eligible * (H1.q_smoke_cat>1)).sum(axis=1) > 0)
    el_trt20 = el_trtq20 + el_bp + el_bmi + el_sm
    el_trt10 = el_trtq10 + el_bp + el_bmi + el_sm

    ## Attending and eligible for each kind of treatment
    att_trtq20 = ((H1.Attending * (H1.QRisk >= 20)).sum(axis=1) > 0)
    att_trtq10 = ((H1.Attending * (H1.QRisk >= 10)).sum(axis=1) > 0)
    att_bp = ((H1.Attending * (H1.q_sbp > 140)).sum(axis=1) > 0)
    att_bmi = ((H1.Attending * (H1.bmi >= 30)).sum(axis=1) > 0) 
    att_sm = ((H1.Attending * (H1.q_smoke_cat>1)).sum(axis=1) > 0)
    att_trt20 = att_trtq20 + att_bp + att_bmi + att_sm
    att_trt10 = att_trtq10 + att_bp + att_bmi + att_sm
        
    ## Most and least deprived quintile 
    dep = (H1.SES == 5)
    ndep = (H1.SES == 1)
        
    stat_off = (H1.Statins_Offered.sum(axis=1)>0)
    aht_off = (H1.Hypertensives_Offered.sum(axis=1)>0)
    sc_off = (H1.SmokingCessation_Offered.sum(axis=1)>0)
    wr_off = (H1.WeightReduction_Offered.sum(axis=1)>0)
    stat = (H1.Statins.sum(axis=1)>0)
    aht = (H1.Hypertensives.sum(axis=1)>0)
    sc = (H1.SmokingCessation.sum(axis=1)>0)
    wr = (H1.WeightReduction.sum(axis=1)>0)
    nhc = H1.Attending.sum()

    NT = np.zeros(5, dtype=int) # unused in this function
    M, S, NT = GetResults_sub(allpop,     0, H, H1, H0, M, S, NT)
    M, S, NT = GetResults_sub(el,     1, H, H1, H0, M, S, NT)
    M, S, NT = GetResults_sub(att, 2, H, H1, H0, M, S, NT)
    M, S, NT = GetResults_sub(t_offered,  3, H, H1, H0, M, S, NT)
    M, S, NT = GetResults_sub(dep,        4, H, H1, H0, M, S, NT)
    M, S, NT = GetResults_sub(ndep,       5, H, H1, H0, M, S, NT)

    N[1] = el.sum()
    N[2] = att.sum()
    N[3] = t_offered.sum()
    N[4] = dep.sum()
    N[5] = ndep.sum()
    N[6] = el_trt20.sum()
    N[7] = el_trt10.sum()
    N[8] = att_trt20.sum()
    N[9] = att_trt10.sum()
    N[10] = stat_off.sum()
    N[11] = aht_off.sum()
    N[12] = sc_off.sum()
    N[13] = wr_off.sum()
    N[14] = stat.sum()
    N[15] = aht.sum()
    N[16] = sc.sum()
    N[17] = wr.sum()
    N[18] = nhc

    return M,S,N



def FirstHC(H1):
    ''' returns time when first HC occurs for each person'''
    FHC = np.zeros((H1.population_size), dtype=np.int) - 1
    for i in range(H1.population_size):
        try:
            FHC[i] = np.where(H1.Attending[i]==1)[0][0]
        except:
            pass # no attendance for this individual
    return FHC


def PostHCOuts(H, H1):
    ''' 10-year improvement in QRisk and other factors, given by HC '''
    att = H1.Attending.sum(axis=1)>0
    HC1 = FirstHC(H1)[att]
    qcon =  H.QRisk[att, HC1 + 10] - H.QRisk[att, HC1]   # change in qrisk over 10 years, neg good, will be pos as increases with age 
    qhc =   H1.QRisk[att, HC1 + 10] - H1.QRisk[att, HC1]
    qdiff = qhc - qcon                # lower change under hc, negative is good 
    scon =  H.q_sbp[att, HC1 + 10] - H.q_sbp[att, HC1]
    shc =   H1.q_sbp[att, HC1 + 10] - H1.q_sbp[att, HC1]
    sdiff = shc - scon 
    dcon =  H.dia[att, HC1 + 10] - H.dia[att, HC1]
    dhc =   H1.dia[att, HC1 + 10] - H1.dia[att, HC1]
    ddiff = dhc - dcon 
    ccon =  H.chol[att, HC1 + 10] - H.chol[att, HC1]
    chc =   H1.chol[att, HC1 + 10] - H1.chol[att, HC1]
    cdiff = chc - ccon 
    
    posthcmeans = [ qdiff.mean(), sdiff.mean(), ddiff.mean(), cdiff.mean() ] 
    posthcsds = [ qdiff.std(), sdiff.std(), ddiff.std(), cdiff.std() ] 
    
    return np.column_stack((posthcmeans, posthcsds))


def GetResults_paper(H, H1, H0, M, S, N, P):
    M, S, N = GetResults_longterm(H, H1, H0, M, S, N)
    P = PostHCOuts(H, H1)
    return M, S, N, P


def BaselineCharsSub(H, sub):
    '''Demographic characteristics of simulated population under base case'''

    return [
        np.bincount(H.gender[sub]).tolist(),  # 0 female 1 male 
        np.bincount(H.eth[sub]).tolist(),
        np.bincount(H.SES[sub]).tolist(),
        np.bincount(H.educ[sub]).tolist(),
        [H.age[sub,0].mean(), np.percentile(H.age[sub,0], 25), np.percentile(H.age[sub,0], 75)],
        [H.QRisk[sub,0].mean(), np.percentile(H.QRisk[sub,0], 25), np.percentile(H.QRisk[sub,0], 75)],
        [(H.QRisk[sub,0]>=20).sum()],
        [np.logical_and(H.QRisk[sub,0]>=10, H.QRisk[sub,0]<20).sum()],
        [H.q_sbp[sub,0].mean(), np.percentile(H.q_sbp[sub,0], 25), np.percentile(H.q_sbp[sub,0], 75)],
        [H.dia[sub,0].mean(), np.percentile(H.dia[sub,0], 25), np.percentile(H.dia[sub,0], 75)],
        [np.logical_or(H.q_sbp[sub,0]>=140, H.dia[sub,0]>=90).sum()],
        [H.q_b_treatedhyp[sub,0].sum()],
        [H.chol[sub,0].mean(), np.percentile(H.chol[sub,0], 25), np.percentile(H.chol[sub,0], 75)],
        [H.q_rati[sub,0].mean(), np.percentile(H.q_rati[sub,0], 25), np.percentile(H.q_rati[sub,0], 75)],
        [H.bmi[sub,0].mean(), np.percentile(H.bmi[sub,0], 25), np.percentile(H.bmi[sub,0], 75)],
        [(H.bmi[sub,0]>=30).sum()],
        [H.glyhb[sub,0].mean(), np.percentile(H.glyhb[sub,0], 25), np.percentile(H.glyhb[sub,0], 75)],
        [(H.glyhb[:,0]>=6.5).sum()],
        [H.diabetes[sub,0].sum()],
        np.bincount(H.q_smoke_cat[sub,0]).tolist()
    ]


def BaselineChars(H1):
    '''Demographic characteristics of simulated population at baseline, all population and specific subsets'''
    allpop = np.ones(H1.population_size, dtype=bool)
    base_all = BaselineCharsSub(H1, allpop)

    ## People who are eligible for HC at any time
    el = H1.eligible.sum(axis=1)>0
    base_el = BaselineCharsSub(H1, el)
    
    ## People who go on to attend HC, still at baseline
    att = H1.Attending.sum(axis=1)>0
    base_att = BaselineCharsSub(H1, att)
    
    flatten = lambda l: [item for sublist in l for item in sublist]
    return np.column_stack((flatten(base_all), flatten(base_el), flatten(base_att)))
