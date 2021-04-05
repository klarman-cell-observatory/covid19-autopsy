import sys
import pandas as pd

def classify(cmd):
    assert cmd in ['cta', 'wta', 'dsp']

    if cmd == 'cta':
        df = pd.read_csv("data/Broad-COVID_CTA_BioProbeCountMatrix.txt", sep='\t')
        df0 = df.loc[list(map(lambda s: s.startswith('RTS00'), df['ProbeDisplayName'].values))]
        df0 = df0.loc[df0['TargetName']!='Negative Probe']
        pd.DataFrame(df0['TargetName'].astype('category').cat.categories).to_csv("cta_target.csv", index=False, header=False)

        df1 = df.loc[list(map(lambda s: s.startswith('RTS01'), df['ProbeDisplayName'].values))]
        df1 = df1.loc[df1['TargetName']!='SARS-CoV-2 Neg']
        pd.DataFrame(df1['TargetName'].astype('category').cat.categories).to_csv("cta_sars_cov_2.csv", index=False, header=False)

    elif cmd == 'wta':
        df = pd.read_csv("data/Broad-COVID_WTA_BioProbeCountMatrix.txt", sep='\t')
        df0 = df.loc[list(map(lambda s: s.startswith('RTS00'), df['ProbeDisplayName'].values))]
        df0 = df0.loc[df0['TargetName']!='Neg Probe']
        pd.DataFrame(df0['TargetName'].astype('category').cat.categories).to_csv("wta_target.csv", index=False, header=False)

        df1 = df.loc[list(map(lambda s: s.startswith('RTS01'), df['ProbeDisplayName'].values))]
        df1 = df1.loc[df1['TargetName']!='SARS-CoV-2 Neg']
        pd.DataFrame(df1['TargetName'].astype('category').cat.categories).to_csv("wta_sars_cov_2.csv", index=False, header=False)

    else:
        df = pd.read_csv("data/NonNormProteins.csv", header=None)
        df = df.loc[(df[1]=='Protein') & (df[2]=='Endogenous')]
        df[[3]].to_csv("dsp_target.csv", index=False, header=False)

if __name__ == '__main__':
    cmd = sys.argv[1]
    classify(cmd)
