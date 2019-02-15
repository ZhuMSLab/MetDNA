// export class Paras {
//     id = 0;
//     polarity: string[] = [];
//     lcType: string[] = [];
//     ce: string[] = [];
//     controlGroup: string[] = [];
//     caseGroup: string[] = [];
//     statType: string = '';
//     adjustP: string = '';
//     pCutoff: number = 0.05;
//     species: string[] = [];
// }
const POL: string[] = ['Positive', 'Negative', 'Both'];
export const polarity: Object = {'both': POL, 'pos': POL.slice(0, 1), 'neg': POL.slice(1, 2)} ;
export const lc: string[] = ['HILIC', 'RP'];
export const ce: string[] = ['10', '20', '30', '35±15', '40', '50', 'NCE15-30-45'];
export const statMethod: string[] = ['Student t-test', 'Wilcox test'];
export const adjustP: string[] = ['No', 'Yes'];
export const species: string[] = ['Homo sapiens (human)', 'Mus musculus (mouse)', 
                                  'Rattus norvegicus (rat)', 'Bos taurus (cow)', 
                                  'Drosophila melanogaster (fruit fly)', 'Gallus gallus (chicken)', 
                                  'Danio rerio (zebrafish)', 'Caenorhabditis elegans (nematode)', 
                                  'Saccharomyces cerevisiae (yeast)', 'Arabidopsis thaliana (thale cress)', 
                                  'Schistosoma mansoni', 'Plasmodum falciparum 3D7 (Malaria)', 
                                  'Trypanosoma brucei', 'Escherichia coli K-12 MG1655',
                                  'Pseudomonas putida KT2440', 'Synechococcus elongatus'];

export const instrument = ['Sciex TripleTOF', 'Agilent QTOF', 'Other QTOF', 'Thermo Orbitrap (HCD)']
export const ceInstt: Object = {'STT': ce.slice(0, 6), 'AQT': ce.slice(0, 3).concat(ce.slice(4, 6)),
                                'OQ': ce.slice(0, 3).concat(ce.slice(4, 6)), 'TOQ': ce.slice(0, 3).concat(ce.slice(4, 7))}