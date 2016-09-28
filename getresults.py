import numpy as np

## TODO do we want different LY and QALY by gender?

def GetResults_longterm(H, H1, M, S, N):
    el_tot = (H1.eligible.sum(axis=1)>0)
    hc_offered = (H1.OfferedHC.sum(axis=1)>0)
    att = H1.Attending * H1.eligible  # H1.Attending should give same answer
    n_att_once = (att.sum(axis=1)>0)
    t_offered = (H1.OfferedTreatment.sum(axis=1)>0)
    stat = (H1.Statins.sum(axis=1)>0)
    aht = (H1.Hypertensives.sum(axis=1)>0)
    sc = (H1.SmokingCessation.sum(axis=1)>0)
    wr = (H1.WeightReduction.sum(axis=1)>0)
    total_hc = att.sum()

    M[0,0] = (H.QALY).mean()
    M[1,0] = (H.QALY[el_tot]).mean()
    M[2,0] = (H.QALY[n_att_once]).mean()
    M[3,0] = (H.QALY[t_offered]).mean()

    M[0,1] = (H1.QALY).mean()
    M[1,1] = (H1.QALY[el_tot]).mean()
    M[2,1] = (H1.QALY[n_att_once]).mean()
    M[3,1] = (H1.QALY[t_offered]).mean()

    M[0,2] = (H1.QALY - H.QALY).mean()
    M[1,2] = (H1.QALY[el_tot] - H.QALY[el_tot]).mean()
    M[2,2] = (H1.QALY[n_att_once] - H.QALY[n_att_once]).mean()
    M[3,2] = (H1.QALY[t_offered] - H.QALY[t_offered]).mean()

    M[0,3] = (H.LY).mean()
    M[1,3] = (H.LY[el_tot]).mean()
    M[2,3] = (H.LY[n_att_once]).mean()
    M[3,3] = (H.LY[t_offered]).mean()

    M[0,4] = (H1.LY).mean()
    M[1,4] = (H1.LY[el_tot]).mean()
    M[2,4] = (H1.LY[n_att_once]).mean()
    M[3,4] = (H1.LY[t_offered]).mean()

    M[0,5] = (H1.LY - H.LY).mean()
    M[1,5] = (H1.LY[el_tot] - H.LY[el_tot]).mean()
    M[2,5] = (H1.LY[n_att_once] - H.LY[n_att_once]).mean()
    M[3,5] = (H1.LY[t_offered] - H.LY[t_offered]).mean()

    S[0,0] = (H.QALY).std()
    S[1,0] = (H.QALY[el_tot]).std()
    S[2,0] = (H.QALY[n_att_once]).std()
    S[3,0] = (H.QALY[t_offered]).std()

    S[0,1] = (H1.QALY).std()
    S[1,1] = (H1.QALY[el_tot]).std()
    S[2,1] = (H1.QALY[n_att_once]).std()
    S[3,1] = (H1.QALY[t_offered]).std()

    S[0,2] = (H1.QALY - H.QALY).std()
    S[1,2] = (H1.QALY[el_tot] - H.QALY[el_tot]).std()
    S[2,2] = (H1.QALY[n_att_once] - H.QALY[n_att_once]).std()
    S[3,2] = (H1.QALY[t_offered] - H.QALY[t_offered]).std()

    S[0,3] = (H.LY).std()
    S[1,3] = (H.LY[el_tot]).std()
    S[2,3] = (H.LY[n_att_once]).std()
    S[3,3] = (H.LY[t_offered]).std()

    S[0,4] = (H1.LY).std()
    S[1,4] = (H1.LY[el_tot]).std()
    S[2,4] = (H1.LY[n_att_once]).std()
    S[3,4] = (H1.LY[t_offered]).std()

    S[0,5] = (H1.LY - H.LY).std()
    S[1,5] = (H1.LY[el_tot] - H.LY[el_tot]).std()
    S[2,5] = (H1.LY[n_att_once] - H.LY[n_att_once]).std()
    S[3,5] = (H1.LY[t_offered] - H.LY[t_offered]).std()

    N[1] = el_tot.sum()
    N[2] = n_att_once.sum()
    N[3] = t_offered.sum()
    N[4] = stat.sum()
    N[5] = aht.sum()
    N[6] = sc.sum()
    N[7] = wr.sum()
    N[8] = total_hc

    return M,S,N

def GetResults_stat(H, H1, M, S, NT):
    stat = (H1.Statins.sum(axis=1)>0)
    M[4,0] = (H.QALY[stat]).mean()
    M[4,1] = (H1.QALY[stat]).mean()
    M[4,2] = (H1.QALY[stat] - H.QALY[stat]).mean()
    M[4,3] = (H.LY[stat]).mean()
    M[4,4] = (H1.LY[stat]).mean()
    M[4,5] = (H1.LY[stat] - H.LY[stat]).mean()
    S[4,0] = (H.QALY[stat]).std()
    S[4,1] = (H1.QALY[stat]).std()
    S[4,2] = (H1.QALY[stat] - H.QALY[stat]).std()
    S[4,3] = (H.LY[stat]).std()
    S[4,4] = (H1.LY[stat]).std()
    S[4,5] = (H1.LY[stat] - H.LY[stat]).std()
    NT[0] = stat.sum()
    return M,S,NT

def GetResults_aht(H, H1, M, S, NT):
    aht = (H1.Hypertensives.sum(axis=1)>0)
    M[5,0] = (H.QALY[aht]).mean()
    M[5,1] = (H1.QALY[aht]).mean()
    M[5,2] = (H1.QALY[aht] - H.QALY[aht]).mean()
    M[5,3] = (H.LY[aht]).mean()
    M[5,4] = (H1.LY[aht]).mean()
    M[5,5] = (H1.LY[aht] - H.LY[aht]).mean()
    S[5,0] = (H.QALY[aht]).std()
    S[5,1] = (H1.QALY[aht]).std()
    S[5,2] = (H1.QALY[aht] - H.QALY[aht]).std()
    S[5,3] = (H.LY[aht]).std()
    S[5,4] = (H1.LY[aht]).std()
    S[5,5] = (H1.LY[aht] - H.LY[aht]).std()
    NT[1] = aht.sum()
    return M,S,NT

def GetResults_sc(H, H1, M, S, NT):
    sc = (H1.SmokingCessation.sum(axis=1)>0)
    M[6,0] = (H.QALY[sc]).mean()
    M[6,1] = (H1.QALY[sc]).mean()
    M[6,2] = (H1.QALY[sc] - H.QALY[sc]).mean()
    M[6,3] = (H.LY[sc]).mean()
    M[6,4] = (H1.LY[sc]).mean()
    M[6,5] = (H1.LY[sc] - H.LY[sc]).mean()
    S[6,0] = (H.QALY[sc]).std()
    S[6,1] = (H1.QALY[sc]).std()
    S[6,2] = (H1.QALY[sc] - H.QALY[sc]).std()
    S[6,3] = (H.LY[sc]).std()
    S[6,4] = (H1.LY[sc]).std()
    S[6,5] = (H1.LY[sc] - H.LY[sc]).std()
    NT[2] = sc.sum()
    return M,S,NT

def GetResults_wr(H, H1, M, S, NT):
    wr = (H1.WeightReduction.sum(axis=1)>0)
    M[7,0] = (H.QALY[wr]).mean()
    M[7,1] = (H1.QALY[wr]).mean()
    M[7,2] = (H1.QALY[wr] - H.QALY[wr]).mean()
    M[7,3] = (H.LY[wr]).mean()
    M[7,4] = (H1.LY[wr]).mean()
    M[7,5] = (H1.LY[wr] - H.LY[wr]).mean()
    S[7,0] = (H.QALY[wr]).std()
    S[7,1] = (H1.QALY[wr]).std()
    S[7,2] = (H1.QALY[wr] - H.QALY[wr]).std()
    S[7,3] = (H.LY[wr]).std()
    S[7,4] = (H1.LY[wr]).std()
    S[7,5] = (H1.LY[wr] - H.LY[wr]).std()
    NT[3] = wr.sum()
    return M,S,NT


def FirstHC(H1):
    ### returns array of size populationsize x simulation time of when first HC occurs
    FHC = np.zeros((H1.population_size,H1.simulation_time),dtype=bool)
    for i in range(H1.population_size):
        try:
            idx = np.where(H1.Attending[i]==1)[0][0]
            FHC[i,idx] = 1
        except:
            pass # no attendance for this individual
    return FHC


def Process(H, H1):

    HC1 = FirstHC(H1) # when was first HC
    popsize = H1.population_size
    el_tot = H1.eligible.sum(axis=1)>0
    eligible = el_tot.sum()
    att_tot = (H1.Attending * H1.eligible)
    attending = (att_tot.sum(axis=1)>0).sum()
    att_q = H1.QRisk[att_tot].mean()
    att_chol = H1.chol[att_tot].mean()
    att_hdl = H1.hdl[att_tot].mean()
    q20plus = (H1.QRisk>20) * HC1
    qhi = (q20plus).sum()
    qhi_presc_stat = (q20plus*H1.Statins_Offered).sum()
    qhi_stat = (q20plus*H1.Statins).sum()
    q20minus = (H1.QRisk<20) * HC1
    qlo = (q20minus).sum()
    qlo_presc_stat = (q20minus*H1.Statins_Offered).sum()
    qlo_stat = (q20minus*H1.Statins).sum()
    nonatt = H1.population_size - HC1.sum()
    natt_all = H1.Attending.sum(axis=1) == 0
    nonatt_q = H1.QRisk[natt_all].mean()
    nonatt_chol = H1.chol[natt_all].mean()
    nonatt_h = H1.hdl[natt_all].mean()
    
    att_sbp = H1.q_sbp[att_tot].mean()
    att_dbp = H1.dia[att_tot].mean()
    qhi_presc_aht = (q20plus*H1.Hypertensives_Offered).sum()
    qhi_aht = (q20plus*H1.Hypertensives).sum()
    qlo_presc_aht = (q20minus*H1.Hypertensives_Offered).sum()
    qlo_aht = (q20minus*H1.Hypertensives).sum()
    nonatt_sbp = H1.q_sbp[natt_all].mean()
    nonatt_dbp = H1.dia[natt_all].mean()
    hbp = H1.q_sbp>140
    att_hbp = hbp * HC1  # of attending individuals, how many have high bp at first HC?
    att_nhbp = att_hbp.sum()

    att_bmi = H1.bmi[att_tot].mean()
    bmi30plus = (H1.bmi>30) * HC1
    hbmi = (bmi30plus).sum()
    hbmi_presc_wr = (bmi30plus*H1.WeightReduction_Offered).sum()
    hbmi_wr = (bmi30plus*H1.WeightReduction).sum()
    bmi30minus = (H1.bmi<30) * HC1
    lobmi = (bmi30minus).sum()
    lobmi_presc_wr = (bmi30minus*H1.WeightReduction_Offered).sum()
    lobmi_wr = (bmi30minus*H1.WeightReduction).sum()
    nonatt_bmi = H1.bmi[natt_all].mean()

    smokers = (H1.q_smoke_cat>1) * HC1
    att_smk = (smokers).sum()
    att_sc = (smokers*H1.SmokingCessation_Offered).sum()
    att_quit1 = (smokers*H1.SmokingCessation).sum()
    nonsmokers = (H1.q_smoke_cat<=1) * HC1
    att_nsmk = (nonsmokers).sum()

    P = [popsize, eligible, attending,
        att_q, att_chol, att_hdl,
               qhi, qhi_presc_stat, qhi_stat, qlo, qlo_presc_stat, qlo_stat,
               nonatt, nonatt_q, nonatt_chol, nonatt_h,
               att_sbp, att_dbp, att_nhbp,
               qhi_presc_aht, qhi_aht, qlo_presc_aht, qlo_aht,
               nonatt_sbp, nonatt_dbp, 
               att_bmi, hbmi, hbmi_presc_wr, hbmi_wr,
               lobmi, lobmi_presc_wr, lobmi_wr,
               nonatt_bmi,
               att_smk, att_sc, att_quit1, att_nsmk
               ]
        
    return P

def meanRiskTrajectories(H, H1, treatment='statins',variable='QRisk',firstHC=True):
    ''' returns mean Risk Factor trajectories for control and HC groups after treatment:
    Output:
    * R_array: Risk factor array for Control group
    * R1_array: Risk factor array for HC group
    * age_array: array containing age at time of risk factor recording
    * takeup_array: array containing time point of first takeup of treatment

    treatments can be: 'statin', 'aht', 'wr', 'sc'

    if firstHC = True, then treatment is only counted if it occurs at first attended HC

    variable can be:
        'QRisk' - Qrisk
        'CVD_risks' - theoretical CVD Risks based upon annualised QRisks
        'chol' - total cholesterol
        'hdl' - HDL cholesterol
        'q_sbp' - Systolic Blood pressure
        'dia' - Diastolic Blood pressure
        'bmi' - BMI
        'q_smoke_cat' - smoking category

    '''

    ctr = 0

    HC1 = FirstHC(H1)

    if treatment == 'statins':
        if firstHC == False:
            subgroup = H1.Statins.sum(axis=1)>0
        else:
            subgroup = (H1.Statins * HC1).sum(axis=1)>0
        treatment_array = H1.Statins

    elif treatment == 'aht':
        if firstHC == False:
            subgroup = H1.Hypertensives.sum(axis=1)>0
        else:
            subgroup = (H1.Hypertensives * HC1).sum(axis=1)>0
            treatment_array = H1.Hypertensives

    elif treatment == 'wr':
        if firstHC == False:
            subgroup = H1.WeightReduction.sum(axis=1)>0
        else:
            subgroup = (H1.WeightReduction * HC1).sum(axis=1)>0
        treatment_array = H1.WeightReduction

    elif treatment == 'sc':
        if firstHC == False:
            subgroup = H1.SmokingCessation_Offered.sum(axis=1)>0
        else:
            subgroup = (H1.SmokingCessation_Offered * HC1).sum(axis=1)>0
        treatment_array = H1.SmokingCessation_Offered

    else:
        print('treatment should be one of the following: \'statins\', \'aht\', \'wr\', \'sc\'')

    # define output variable

    varnames =['QRisk','CVD_risks', 'chol','hdl','q_sbp','dia','bmi','q_smoke_cat']

    if variable in varnames:

        C_string = 'C_array = H.%s.copy()' % variable
        HC_string = 'HC_array = H1.%s.copy()' % variable

        # load these strings
        exec(C_string)
        exec(HC_string)

    else:
        print('unknown variable name \'%s' % variable)

    # now extract variable trajectories

    ctr = 0

    R_array = np.zeros((treatment_array.sum(),H.simulation_time)) # QRisk for statin takers
    R1_array = np.zeros((treatment_array.sum(),H.simulation_time)) # QRisk for control

    age_array = np.zeros(treatment_array.sum())
    takeup_array = np.zeros((H.population_size,H.simulation_time),dtype=bool)

    # loop through populations
    for i in range(H.population_size):
        if subgroup[i] == True:

            # first treatment takeup:

            t_idx = np.where(treatment_array[i]==1)[0][0]

            # enter when treatment was taken up
            takeup_array[i,t_idx] = 1
            if variable == 'q_smoke_cat':
                r = C_array[i,t_idx:] > 1
                R_array[ctr,:r.size] = np.array(r,dtype=int)
                r1 = HC_array[i,t_idx:] > 1
                R1_array[ctr,:r1.size] = np.array(r1,dtype=int)

            else:
                r = C_array[i,t_idx:]
                R_array[ctr,:r.size] = r
                r1 = HC_array[i,t_idx:]
                R1_array[ctr,:r1.size] = r1




            age_array[ctr] = H.age[i,t_idx]

            # increase counter
            ctr += 1

    return R_array,R1_array,age_array,takeup_array







def EventTrajectories(H, H1, treatment='statins',firstHC=True):
    ''' returns mean event trajectories for control and HC groups after treatment:
    * NE = number of CVD/Stroke/IHD events per year post treatment
    * NA = number of alive individuals, correspondingly
    * event rate = NE/NA

    if firstHC = True, then treatment is only counted if it occurs at first attended HC

    treatments can be: 'statins', 'aht', 'wr', 'sc'
    '''

    ctr = 0

    HC1 = FirstHC(H1)

    if treatment == 'statins':
        if firstHC == False:
            subgroup = H1.Statins.sum(axis=1)>0
        else:
            subgroup = (H1.Statins * HC1).sum(axis=1)>0
        treatment_array = H1.Statins

    elif treatment == 'aht':
        if firstHC == False:
            subgroup = H1.Hypertensives.sum(axis=1)>0
        else:
            subgroup = (H1.Hypertensives * HC1).sum(axis=1)>0
        treatment_array = H1.Hypertensives

    elif treatment == 'wr':
        if firstHC == False:
            subgroup = H1.WeightReduction.sum(axis=1)>0
        else:
            subgroup = (H1.WeightReduction * HC1).sum(axis=1)>0
        treatment_array = H1.WeightReduction

    elif treatment == 'sc':
        if firstHC == False:
            subgroup = H1.SmokingCessation_Offered.sum(axis=1)>0
        else:
            subgroup = (H1.SmokingCessation_Offered * HC1).sum(axis=1)>0
        treatment_array = H1.SmokingCessation_Offered

    else:
        print('treatment should be one of the following: \'statin\', \'aht\', \'wr\', \'sc\'')

    CVD_HC = np.zeros((subgroup.sum(),H.simulation_time)) # CVD events for HC simulation
    CVD_C = np.zeros((subgroup.sum(),H.simulation_time)) #  CVD events for control

    Stroke_HC = np.zeros((subgroup.sum(),H.simulation_time)) # Stroke events for HC simulation
    Stroke_C = np.zeros((subgroup.sum(),H.simulation_time)) #  Stroke events for control

    IHD_HC = np.zeros((subgroup.sum(),H.simulation_time)) # IHD events for HC simulation
    IHD_C = np.zeros((subgroup.sum(),H.simulation_time)) #  IHD events for control

    LC_HC = np.zeros((subgroup.sum(),H.simulation_time)) # lung cancer events for HC simulation
    LC_C = np.zeros((subgroup.sum(),H.simulation_time)) #  lung cancer events for control

    alive_HC =  np.zeros((subgroup.sum(),H.simulation_time),dtype =bool)
    alive_C = np.zeros((subgroup.sum(),H.simulation_time),dtype =bool)

    # loop through populations
    for i in range(H.population_size):
        if subgroup[i] == True:
            # first statins prescription:

            t_idx = np.where(treatment_array[i]==1)[0][0]

            # for an event, everyone counts as 'alive' who either is alive in a given year or dies in that year (as events precede death)
            aHC = (H1.alive[i,t_idx:]==1) + (H1.Death[i,t_idx:]==1)
            alive_HC[ctr,:aHC.size] = aHC

            aC = (H.alive[i,t_idx:]==1) + (H.Death[i,t_idx:]==1)
            alive_C[ctr,:aC.size] = aC

            # now count CVD events
            cv_hc = H1.CVD_events[i,t_idx:]
            CVD_HC[ctr,:cv_hc.size] = cv_hc

            cv_c = H.CVD_events[i,t_idx:]
            CVD_C[ctr,:cv_c.size] = cv_c


            # now count stroke events
            s_hc = H1.StrokeEvents[i,t_idx:]
            Stroke_HC[ctr,:s_hc.size] = s_hc

            s_c = H.StrokeEvents[i,t_idx:]
            Stroke_C[ctr,:s_c.size] = s_c


            # now count IHD events
            i_hc = H1.IHDEvents[i,t_idx:]
            IHD_HC[ctr,:i_hc.size] = i_hc

            i_c = H.IHDEvents[i,t_idx:]
            IHD_C[ctr,:i_c.size] = i_c

            # now count lung cancer events
            l_hc = H1.LungCancerEvents[i,t_idx:]
            LC_HC[ctr,:l_hc.size] = l_hc

            l_c = H.LungCancerEvents[i,t_idx:]
            LC_C[ctr,:l_c.size] = l_c


            ctr += 1

    A = alive_C.sum(axis=0)
    A1 = alive_HC.sum(axis=0)
    C = CVD_C.sum(axis=0)
    C1 = CVD_HC.sum(axis=0)
    S = Stroke_C.sum(axis=0)
    S1 = Stroke_HC.sum(axis=0)
    I = IHD_C.sum(axis=0)
    I1 = IHD_HC.sum(axis=0)
    L = LC_C.sum(axis=0)
    L1 = LC_HC.sum(axis=0)

    return A,A1,C,C1,S,S1,I,I1,L,L1



def DeathTrajectories(H, H1, treatment='statins',firstHC=True,CauseOfDeath='all'):
    ''' returns mean Death trajectories for control and HC groups after treatment:
    * NE = number of Deaths per year post treatment
    * NA = number of alive/just dying individuals, correspondingly
    * death rate = NE/NA

    if firstHC = True, then treatment is only counted if it occurs at first attended HC
    if CauseOfDeath == 'CVD', only CVD deaths are counted
    if CauseOfDeath == 'Lung Cancer',  only lung cancer deaths are counted


    treatments can be: 'statins', 'aht', 'wr', 'sc'
    '''

    ctr = 0

    HC1 = FirstHC(H1)

    if treatment == 'statins':
        if firstHC == False:
            subgroup = H1.Statins.sum(axis=1)>0
        else:
            subgroup = (H1.Statins * HC1).sum(axis=1)>0
        treatment_array = H1.Statins

    elif treatment == 'aht':
        if firstHC == False:
            subgroup = H1.Hypertensives.sum(axis=1)>0
        else:
            subgroup = (H1.Hypertensives * HC1).sum(axis=1)>0
        treatment_array = H1.Hypertensives

    elif treatment == 'wr':
        if firstHC == False:
            subgroup = H1.WeightReduction.sum(axis=1)>0
        else:
            subgroup = (H1.WeightReduction * HC1).sum(axis=1)>0
        treatment_array = H1.WeightReduction

    elif treatment == 'sc':
        if firstHC == False:
            subgroup = H1.SmokingCessation_Offered.sum(axis=1)>0
        else:
            subgroup = (H1.SmokingCessation_Offered * HC1).sum(axis=1)>0
        treatment_array = H1.SmokingCessation_Offered
    else:
        print('treatment should be one of the following: \'statin\', \'aht\', \'wr\', \'sc\'')

    D_HC = np.zeros((subgroup.sum(),H.simulation_time)) # Death events for HC simulation
    D_C = np.zeros((subgroup.sum(),H.simulation_time)) #  Deaths for control

    alive_HC =  np.zeros((subgroup.sum(),H.simulation_time),dtype =bool)
    alive_C = np.zeros((subgroup.sum(),H.simulation_time),dtype =bool)

    if CauseOfDeath == 'CVD':
        DiseaseDeath = (H.CauseOfDeath == 'Stroke') + (H.CauseOfDeath == 'IHD')
    if CauseOfDeath == 'Lung Cancer':
        DiseaseDeath = H.CauseOfDeath == 'Lung Cancer'
    if CauseOfDeath == 'Dementia':
        DiseaseDeath = H.CauseOfDeath == 'Dementia'
    if CauseOfDeath == 'Other':
        DiseaseDeath = H.CauseOfDeath == 'Other'
    if CauseOfDeath == 'None':
        DiseaseDeath = H.CauseOfDeath == 'None'

    # loop through populations
    for i in range(H.population_size):
        if subgroup[i] == True:

            if CauseOfDeath == 'all':
            # first statins prescription:

                t_idx = np.where(treatment_array[i]==1)[0][0]

                # for an event, everyone counts as 'alive' who either is alive in a given year or dies in that year (as events precede death)
                aHC = (H1.alive[i,t_idx:]==1) + (H1.Death[i,t_idx:]==1)
                alive_HC[ctr,:aHC.size] = aHC

                aC = (H.alive[i,t_idx:]==1) + (H.Death[i,t_idx:]==1)
                alive_C[ctr,:aC.size] = aC

                # now count Deaths
                d_hc = H1.Death[i,t_idx:]
                D_HC[ctr,:d_hc.size] = d_hc

                d_c = H.Death[i,t_idx:]
                D_C[ctr,:d_c.size] = d_c

                ctr += 1

            else:

                if DiseaseDeath[i] == True:

                    t_idx = np.where(treatment_array[i]==1)[0][0]

                    # for an event, everyone counts as 'alive' who either is alive in a given year or dies in that year (as events precede death)
                    aHC = (H1.alive[i,t_idx:]==1) + (H1.Death[i,t_idx:]==1)
                    alive_HC[ctr,:aHC.size] = aHC

                    aC = (H.alive[i,t_idx:]==1) + (H.Death[i,t_idx:]==1)
                    alive_C[ctr,:aC.size] = aC

                    # now count Deaths
                    d_hc = H1.Death[i,t_idx:]
                    D_HC[ctr,:d_hc.size] = d_hc

                    d_c = H.Death[i,t_idx:]
                    D_C[ctr,:d_c.size] = d_c

                    ctr += 1

    A = alive_C.sum(axis=0)
    A1 = alive_HC.sum(axis=0)
    D = D_C.sum(axis=0)
    D1 = D_HC.sum(axis=0)

    return A,A1,D,D1


def compute_mean_ignoring_zeros(array):
    a_mean = np.zeros(array.shape[1])
    a_std = np.zeros(array.shape[1]) 
    for i in range(array.shape[1]):
        a_pop = array[array[:,i]!=0,i]
        a_mean[i] = a_pop.mean()
        a_std[i] = a_pop.std()

    return a_mean, a_std


def ShortTermOutcomesTrt(H, H1, trt):
    Q,Q1,age_array,takeup_array = meanRiskTrajectories(H, H1, treatment=trt, variable='QRisk')
    CH,CH1,age_array,takeup_array = meanRiskTrajectories(H, H1, treatment=trt, variable='chol')
    HDL,HDL1,age_array,takeup_array = meanRiskTrajectories(H, H1, treatment=trt, variable='hdl')
    SBP,SBP1,age_array,takeup_array = meanRiskTrajectories(H, H1, treatment=trt, variable='q_sbp')
    DBP,DBP1,age_array,takeup_array = meanRiskTrajectories(H, H1, treatment=trt, variable='dia')
    BMI,BMI1,age_array,takeup_array = meanRiskTrajectories(H, H1, treatment=trt, variable='bmi')
    SM,SM1,age_array,takeup_array = meanRiskTrajectories(H, H1, treatment=trt, variable='q_smoke_cat')

    # compute mean trajectories of these
    Q,QSD = compute_mean_ignoring_zeros(Q)
    Q1,Q1SD = compute_mean_ignoring_zeros(Q1)
    CH,CHSD = compute_mean_ignoring_zeros(CH)
    CH1,CH1SD = compute_mean_ignoring_zeros(CH1)
    HDL,HDLSD = compute_mean_ignoring_zeros(HDL)
    HDL1,HDL1SD = compute_mean_ignoring_zeros(HDL1)
    SBP,SBPSD = compute_mean_ignoring_zeros(SBP)
    SBP1,SBP1SD = compute_mean_ignoring_zeros(SBP1)
    DBP,DBPSD = compute_mean_ignoring_zeros(DBP)
    DBP1,DBP1SD = compute_mean_ignoring_zeros(DBP1)
    BMI,BMISD = compute_mean_ignoring_zeros(BMI)
    BMI1,BMI1SD = compute_mean_ignoring_zeros(BMI1)
    SM = SM.sum(axis=0) # counting number of smokers here
    SM1 = SM1.sum(axis=0)

    A,A1,C,C1,S,S1,I,I1,L,L1  = EventTrajectories(H, H1, treatment=trt)
    AD,AD1,D,D1 = DeathTrajectories(H, H1, treatment=trt)
    ADcvd,ADcvd1,Dcvd,Dcvd1 = DeathTrajectories(H, H1, treatment=trt, CauseOfDeath='CVD')
    ADlc,ADlc1,Dlc,Dlc1 = DeathTrajectories(H, H1, treatment=trt, CauseOfDeath='Lung Cancer')
    ADdem,ADdem1,Ddem,Ddem1 = DeathTrajectories(H, H1, treatment=trt, CauseOfDeath='Dementia')
    ADoth,ADoth1,Doth,Doth1 = DeathTrajectories(H, H1, treatment=trt, CauseOfDeath='Other')

    # events
    cvd_events_t = np.cumsum(C1)
    cvd_events_c = np.cumsum(C)
    lc_events_t = np.cumsum(L1)
    lc_events_c = np.cumsum(L)
    deaths_t = np.cumsum(D1)
    deaths_c = np.cumsum(D)
    deaths_cvd_t = np.cumsum(Dcvd1)
    deaths_cvd_c = np.cumsum(Dcvd)
    deaths_lc_t = np.cumsum(Dlc1)
    deaths_lc_c = np.cumsum(Dlc)
    deaths_dem_t = np.cumsum(Ddem1)
    deaths_dem_c = np.cumsum(Ddem)
    deaths_oth_t = np.cumsum(Doth1)
    deaths_oth_c = np.cumsum(Doth)
            
    resn = np.column_stack((Q1, Q, CH1, CH, HDL1, HDL,
                           SBP1, SBP, DBP1, DBP, BMI1, BMI, SM1, SM,
                           Q1SD, QSD, CH1SD, CHSD, HDL1SD, HDLSD,
                           SBP1SD, SBPSD, DBP1SD, DBPSD, BMI1SD, BMISD, 
                           cvd_events_t, cvd_events_c,
                           lc_events_t, lc_events_c,
                           deaths_t, deaths_c,
                           deaths_cvd_t, deaths_cvd_c,
                           deaths_lc_t, deaths_lc_c,
                           deaths_dem_t, deaths_dem_c,
                           deaths_oth_t, deaths_oth_c))
    
    res = np.empty((resn.shape[0], resn.shape[1] + 1), dtype=object)
    res[:,0] = np.array([trt]*Q1.size) # assumes immutable
    res[:,1:res.shape[1]] = resn 
               
    return res[[0,1,4,9],]


def ShortTermOutcomes(H, H1):
    ST = ShortTermOutcomesTrt(H, H1, 'statins')
    AH = ShortTermOutcomesTrt(H, H1, 'aht')
    WM = ShortTermOutcomesTrt(H, H1, 'wr')
    SC = ShortTermOutcomesTrt(H, H1, 'sc')
    return np.vstack((ST, AH, WM, SC))


def GetResults_all(H, H1, M, S, N, P, ST):
    M, S, N = GetResults_longterm(H, H1, M, S, N)
    P = Process(H, H1)
    ST = ShortTermOutcomes(H, H1)
    return M, S, N, P, ST


## with artificially boosted treated subsets

def ShortTermOutcomesSub(H, HS, HA, HW, HC):
    ST = ShortTermOutcomesTrt(H, HS, 'statins')
    AH = ShortTermOutcomesTrt(H, HA, 'aht')
    WM = ShortTermOutcomesTrt(H, HW, 'wr')
    SC = ShortTermOutcomesTrt(H, HC, 'sc')
    return np.vstack((ST, AH, WM, SC))

def GetResults_allSub(H, H1, HS, HA, HW, HC, M, S, N, NT, P, ST):
    M, S, N = GetResults_longterm(H, H1, M, S, N)
    M, S, NT = GetResults_stat(H, HS, M, S, NT)
    M, S, NT = GetResults_aht(H, HA, M, S, NT)
    M, S, NT = GetResults_wr(H, HW, M, S, NT)
    M, S, NT = GetResults_sc(H, HC, M, S, NT)
    P = Process(H, H1)
    ST = ShortTermOutcomesSub(H, HS, HA, HW, HC)
    return M, S, N, NT, P, ST
