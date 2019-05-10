import pandas as pd
def convert_xlsx(xlsx_name,tsv_name):
    df = pd.read_excel(xlsx_name)
    df.to_csv(tsv_name,sep='\t',index=False)
if __name__ == "__main__":
    convert_xlsx('metadata.xlsx','metadata.tsv')