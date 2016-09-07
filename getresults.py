## TODO do we want different LY and QALY by gender?

def GetResults_all(H, H1, M, S, N):
    el_tot = (H1.eligible.sum(axis=1)>0)
    hc_offered = (H1.OfferedHC.sum(axis=1)>0)
    att_tot = H1.Attending * H1.eligible  # H1.Attending should give same answer
    n_att_once = (att_tot.sum(axis=1)>0)
    t_offered = (H1.OfferedTreatment.sum(axis=1)>0)
    stat = (H1.Statins.sum(axis=1)>0)
    aht = (H1.Hypertensives.sum(axis=1)>0)
    sc = (H1.SmokingCessation.sum(axis=1)>0)
    wr = (H1.WeightReduction.sum(axis=1)>0)
    total_hc = att_tot.sum()
    
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
