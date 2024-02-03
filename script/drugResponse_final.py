#source activate
#conda activate pytorch_gpu
#Author: Shi Yulong
#Date: 2023-04-18
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Drug sensitivity prediction for lung cancer patients with EGFR mutation.')
    parser.add_argument('-m', '--mutation', required=True, type=str,
                        help='EGFR mutation: L858R,L858R+T790M,etc.')
    args = parser.parse_args()

mutation = args.mutation

import pandas as pd
import numpy as np
from DeepPurpose import utils, dataset
from DeepPurpose import DTI as models
import warnings
import pickle
import os
warnings.filterwarnings("ignore")

EGFR_fasta = '''MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEV\
VLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALA\
VLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDF\
QNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGC\
TGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYV\
VTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFK\
NCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAF\
ENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKL\
FGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCN\
LLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVM\
GENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVV\
ALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGS\
GAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGI\
CLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAA\
RNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSY\
GVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPK\
FRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQ\
QGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTED\
SIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLN\
TVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRV\
APQSSEFIGA
'''

#'Olmutinib':'CN1CCN(C2=CC=C(C=C2)NC3=NC(OC4=CC=CC(NC(C=C)=O)=C4)=C5SC=CC5=N3)CC1',
EGFR_TKI_set = {'Gefitinib':'COC1=C(C=C2C(NC3=CC(Cl)=C(C=C3)F)=NC=NC2=C1)OCCCN4CCOCC4',
             'Erlotinib':'COCCOC1=CC2=C(C(NC3=CC(C#C)=CC=C3)=NC=N2)C=C1OCCOC',
             'Icotinib':'C#CC1=CC=CC(NC2=NC=NC3=CC4=C(OCCOCCOCCO4)C=C23)=C1',
             'Afatinib':'CN(C/C=C/C(NC1=C(C=C2N=CN=C(C2=C1)NC3=CC(Cl)=C(C=C3)F)O[C@H]4CCOC4)=O)C',
             'Dacomitinib':'COC1=C(C=C2C(NC3=CC(Cl)=C(C=C3)F)=NC=NC2=C1)NC(/C=C/CN4CCCCC4)=O',
             'Osimertinib':'COC1=C(C=C(C(N(CCN(C)C)C)=C1)NC(C=C)=O)NC2=NC=CC(C3=CN(C4=C3C=CC=C4)C)=N2',
             'Almonertinib':'COC1=CC(N(CCN(C)C)C)=C(C=C1NC2=NC(C3=CN(C4=C3C=CC=C4)C5CC5)=CC=N2)NC(C=C)=O',
             'Furmonertinib':'CN(C1=C2C=CC=C1)C=C2C3=NC(NC(C=C(C(N(CCN(C)C)C)=N4)NC(C=C)=O)=C4OCC(F)(F)F)=NC=C3'
            }

def point_mutation(fasta, old_amino, num, new_amino):
    if fasta[num] == old_amino:
        fasta[num] = new_amino
    else:
        print("Error: wrong mutation site.")
        os._exit(1)

def deletion_mutation_type1(fasta, amino1, num1, amino2, num2):
    if fasta[num1] == amino1 and fasta[num2] == amino2:
        for id in range(num1,num2+1):
            fasta.pop(id)
    else:
        print("Error: wrong mutation site.")
        os._exit(1)

def deletion_mutation_type2(fasta, amino1, num1):
    if fasta[num1] == amino1:
        fasta.pop(num1)
    else:
        print("Error: wrong mutation site.")
        os._exit(1)

def insert_mutation(fasta, amino, num, insert_amino):
    if fasta[num] == amino:
        fasta[num] = amino + insert_amino
    else:
        print("Error: wrong mutation site.")
        os._exit(1)

def duplicate_mutation(fasta, amino, num, dup_amino):
    #A767dupASV
    if fasta[num] == amino:
        amino2 = fasta[num+len(dup_amino)-1]
        fasta[num+len(dup_amino)-1] = amino2 + dup_amino
    else:
        print("Error: wrong mutation site.")
        os._exit(1)

def delins_mutation_type1(fasta, amino1, num1, amino2, num2, delins_amino):
    if fasta[num1] == amino1 and fasta[num2] == amino2:
        for id in range(num1+1,num2+1):
            fasta.pop(id)
        fasta[num1] = delins_amino
    else:
        print("Error: wrong mutation site.")
        os._exit(1)

def delins_mutation_type2(fasta, amino1, num1, delins_amino):
    if fasta[num1] == amino1:
        fasta[num1] = delins_amino
    else:
        print("Error: wrong mutation site.")
        os._exit(1)

def Ex19del_mutation(fasta):
    for id in range(729,761+1):
        fasta.pop(id)

def generate_mutated_fasta(wild_fasta, mutation_type):
    amino_id = dict()
    num = 1
    for amino in list(wild_fasta):
        amino_id[num] = amino
        num+=1    
    fasta = amino_id.copy()
    
    for mutation in mutation_type.split('+'):
        print('Mutation:', mutation)
        if mutation == 'WT':
            pass
        elif 'delins' in mutation:
            if mutation[9:15] == 'delins':
                delins_mutation_type1(fasta, mutation[0], int(mutation[1:4]), mutation[5], int(mutation[6:9]), mutation[15:])
            elif mutation[4:10] == 'delins':
                delins_mutation_type2(fasta, mutation[0], int(mutation[1:4]), mutation[10:])
        elif 'del' in mutation and 'Ex19' not in mutation:
            if mutation[9:] == 'del':
                deletion_mutation_type1(fasta, mutation[0], int(mutation[1:4]), mutation[5], int(mutation[6:9]))
            elif mutation[4:] == 'del':
                deletion_mutation_type2(fasta, mutation[0], int(mutation[1:4]))
        elif 'Ex19del' in mutation:
            Ex19del_mutation(fasta)
        elif 'ins' in mutation:
            insert_mutation(fasta, mutation[0], int(mutation[1:4]), mutation[7:])
        elif 'dup' in mutation:
            duplicate_mutation(fasta, mutation[0], int(mutation[1:4]), mutation[7:])
        elif len(mutation) == 5:
            point_mutation(fasta, mutation[0], int(mutation[1:4]), mutation[4])
        else:
            print("Error: unrecognized mutation type.")
            os._exit(1)
    out_fasta = ''.join(fasta.values()).strip()
    return out_fasta

#point_mutation(aa,'T', 790, 'M')                      #T790M
#deletion_mutation_type1(aa,'E',746,'A',750)           #E746_A750del
#deletion_mutation_type2(aa, 'V', '834')               #V834del
#insert_mutation(aa, 'D', 770, 'SVD')                  #D770insSVD
#duplicate_mutation(aa, 'A', 767, 'ASV')               #A767dupASV
#delins_mutation_type1(aa, 'L', 747, 'P', 753, 'S')    #L747_P753delinsS
#delins_mutation_type2(aa,'D', '770', 'GY')            #D770delinsGY

def predict_drug_score_response_probability(fasta, mutation, drug, DL_model,LR_model):
    mutation_fasta = generate_mutated_fasta(fasta, mutation)
    drug_smiles =[EGFR_TKI_set[drug]]
    y_pred = round(models.virtual_screening(drug_smiles, mutation_fasta, DL_model)[0],2)
    response_type = list(LR_model.predict(np.array([[y_pred]])))[0]
    if response_type == -1:
        response = 'CR/PR'
    elif response_type == 0:
        response = 'SD'
    elif response_type == 1:
        response = 'PD'
    probability = LR_model.predict_proba(np.array([[y_pred]]))
    print(probability)
    #CR_PR=round(probability[0,0],3)
    CR_PR = "%.1f%%" % (probability[0,0] * 100)
    SD="%.1f%%" % (probability[0,1] * 100)
    PD="%.1f%%" % (probability[0,2] * 100)
    print(y_pred,response,CR_PR,SD,PD)
    return y_pred,response,CR_PR,SD,PD

# load deep learning model
DL_model = models.model_pretrained(path_dir = 'model/Morgan_CNN')
# load Logistic regression model 
LR_model = pickle.load(open("model/Logistic_regression/D3EGFR_Logistic_Regression.dat","rb"))

mutation_fasta = generate_mutated_fasta(EGFR_fasta, mutation)
print('Mutant Fasta: ', mutation_fasta)
with open('Prediction/' + mutation + '.txt', 'w', encoding='utf-8') as f:
    f.write('Mutation Drug Predicted_response P(CR/PR) P(SD) P(PD) Score'+'\n')
    for drug in EGFR_TKI_set.keys():
        pred_score,pred_response,CR_PR,SD,PD = predict_drug_score_response_probability(EGFR_fasta,mutation,drug,DL_model,LR_model)
        f.write(mutation + ' ' + drug + ' ' + pred_response + ' ' + str(CR_PR) + ' ' + str(SD) + ' ' + str(PD) + ' ' + str(pred_score) + '\n')
        print(mutation + ' ' + drug + ' ' + pred_response + ' ' + str(CR_PR) + ' ' + str(SD) + ' ' + str(PD) + ' ' + str(pred_score) + '\n')
