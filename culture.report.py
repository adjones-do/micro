import collections
import pandas as pd
import re

def culture_report():
    
    #Microbiology culture activity report
    source = "path/to/file"
    
    f = open(source, 'r')
    
    subject = f.read()
    
    match = re.findall(r'''(?<=icrobiology\sCulture\sActivity\sReport)(.*?)(?=\d\d:\d\d:\d\d\s*M)''', subject, flags=re.S)
    
    culture_accession = []
    pt_name = []
    final_bug = []
    culture_date = []
    source = []
    dept = []
    
    for i in match:
        a = re.search(r"(?<=Accession\sNumber:\s{8})(.*?)(?=\s*Body)", i)
        if a:
            culture_accession.append(a.group())
        if not a:
            culture_accession.append('None')
        
        b = re.search(r"(?<=Name:\s{3})(.*?)(?=\s*Provider|\s*\d)", i)
        if b:
            pt_name.append(b.group())
        if not b:
            pt_name.append('None')
            
        c = re.search(r"(?<=Final\s{36}Status:\s{3}Verified\s{74})(.*?)(?=\s*\d\d/\d\d/\d\d|=)", i, flags=re.S)
        c_amend = re.search(r"(?<=Amend.{36}Status:\s{3}Verified\s{74})(.*?)(?=\s*\d\d/\d\d/\d\d|=)", i, flags=re.S)
        c_performed = re.search(r"(?<=Final\s{36}Status:\s{3}Performed\s{73})(.*?)(?=\s*\d\d/\d\d/\d\d|=)", i, flags=re.S)
        c_corrected = re.search(r"(?<=Corr.{37}Status:\s{3}Verified\s{74})(.*?)(?=\s*\d\d/\d\d/\d\d|=)", i, flags=re.S)
        if c_amend:
            final_bug.append(c_amend.group())
        elif c_corrected:
            final_bug.append(c_corrected.group())
        elif c:
            final_bug.append(c.group())
        elif c_performed:
            final_bug.append(c_performed.group())
        elif not (c or c_amend or c_performed or c_corrected):
            final_bug.append('None')
            
        d = re.search(r"(?=\d*/\d*/\d*\s\d*:\d*\s[A-Z\.]*\s*Final)(\d*/\d*/\d*)", i)
        if d:
            culture_date.append(d.group())
        if not d:
            culture_date.append('None')
        
        e = re.search(r"(?<=:\s{18})([a-zA-Z\s]+?)(?=\s*Accession)", i)
        if e:
            source.append(e.group())
        if not e:
            source.append('None')
            
        f = re.search(r"(?<=Location:\s{6}0125\D/0125\D-)(.*?)(?=/)", i)
        if f:
            dept.append(f.group())
        if not f:
            dept.append('None')
            
    urine_culture_list = [list(x) for x in zip(pt_name, culture_accession, culture_date, final_bug, source, dept)]
    
    uc_df = pd.DataFrame(urine_culture_list, columns = ['name', 'culture_acc', 'culture_date', 'final_bug', 'culture_source', 'culture_ordering_dept'])

    uc_df = uc_df[uc_df.culture_acc != 'None']
    uc_df['culture_date'] = uc_df['culture_date'].astype('str')
    uc_df['culture_date'] = pd.to_datetime(uc_df['culture_date'], format='%m/%d/%y') 
    uc_df['final_bug'] = uc_df['final_bug'].replace(r'[\r\n]', ' ', regex=True) 
    uc_df['final_bug'] = uc_df['final_bug'].replace(r'\s+', ' ', regex=True)
    
    split_cdx_df = uc_df['final_bug'].apply(cdx_split)
    
    df = split_cdx_df.apply(pd.Series)
    
    df = df.rename(columns = lambda x : 'bug_' + str(x))
    df = df.astype(str)
    df['bug_0'] = df.bug_0.str.replace("'", "")
    df['bug_1'] = df.bug_1.str.replace("'", "")
    df['bug_2'] = df.bug_2.str.replace("'", "")
    df_1 = df.bug_0.str.strip("[]").str.split(pat=",", expand=True)
    df_1 = df_1.rename(columns={0 : 'cfu_1', 1 : 'bug_1'})
    df_2 = df.bug_1.str.strip("[]").str.split(pat=",", expand=True)
    df_2 = df_2.rename(columns={0: 'cfu_2', 1: 'bug_2'})
    df_3 = df.bug_2.str.strip("[]").str.split(pat=",", expand=True)
    df_3 = df_3.rename(columns={0: 'cfu_3', 1: 'bug_3'})
    df = pd.concat([df_1, df_2, df_3], axis=1, sort=False)
       
    uc_df = pd.concat([uc_df, df], axis=1, sort=False)
    
    uc_df['cfu_1'] = uc_df['cfu_1'].str.replace('nan', '0')
    uc_df['cfu_1'] = uc_df['cfu_1'].astype(int)
    uc_df[['bug_1', 'bug_2', 'bug_3']] = uc_df[['bug_1', 'bug_2', 'bug_3']].applymap(lambda x: x.strip() if isinstance(x, str) else x)

    #uc_df.to_csv('path/to/file')
    
    return uc_df


def cdx_split(bug):

    bugex = re.sub(r'(?<=\d),(?=\d+)','',bug)
    
    neg = re.search(r'^no\sgrowth.*', bugex, flags=re.IGNORECASE)
    pos = re.search(r'^>.*|^\d.*', bugex)
    other = re.search(r'^\w.*|^<.*', bugex)

    dict = collections.defaultdict(list)
        
    if neg:
        dict[0].append('0')
        
    elif pos:
        m = pos.group()
        n = re.findall(r'(\d+)(?=\scfu/ml)', m, flags=re.IGNORECASE)
        o = re.findall(r'(?<=cfu/ml\s)([a-zA-Z\s/]*?)(?=[0-9:;.,>(#-]|with|and\s|Routine|Ceph|Was|unable|Suscept|$|Specimen|PLEASE|worked)', m, flags=re.IGNORECASE)
        
        for i in range(0, len(n)):
            dict[i].append(n[i])
        
        for i in range(0, len(o)):
            dict[i].append(o[i])
        
    elif other:
        p = other.group()
        
        neg_other = re.search(r'no\sgrowth,\s<1000\scfu/ml', p, flags=re.IGNORECASE)
        
        if neg_other:
            dict[0].append('0')
        
        else:
            q = re.findall(r'(\d+)(?=\scfu/ml)', p, flags=re.IGNORECASE)
            r = re.findall(r'(?<=cfu/ml\s)([a-zA-Z\s]*?)(?=[0-9:;.,>(#]|with|and|Routine|Ceph|Was|unable|suscept|$|Specimen|PLEASE|worked)', p, flags=re.IGNORECASE)

            for i in range(0, len(q)):
                dict[i].append(q[i])

            for i in range(0, len(r)):
                dict[i].append(r[i])
                
    else:
        dict[0].append(bug)
        
    dict_list = [list(x) for x in dict.values()]
    
    return dict_list