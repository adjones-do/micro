import re
import pandas as pd

def urine_report():

    source = 'C:/Users/austin.d.jones43/Desktop/urine.pi.data/arkray.accession.au4050-1&2.14jan-29feb.01aug-30sep.txt'
    
    f = open(source,'r')

    subject = f.read()
    
    #Find each entry so that we can loop through them
    match = re.findall(r'''(?<=DEMOGRAPHICS)(.*?)(?=UA\spH)''', subject, flags=re.S|re.M)
    
    #Create the lists which we'll append to later
    name = []
    accession = []
    ua_date = [] #this is the collected date/time
    age = []
    mrn = []
    sex = []

    for i in match:
        
        a = re.search(r'''(.*,\s.*?)(?=\s*\d*\syears|\s*\d*\smonths)''', i)
        if a:
            name.append(a.group())
        if not a:
            name.append('None')
            
        b = re.search(r'''(\d{4}-20-\d{3}-\d{4})''', i)
        if b:
            accession.append(b.group())
        if not b:
            accession.append('None')
            
        c = re.search(r'''(\d{2}\w{3}20\s\d{4})(?=\s{21}UA\sAppear)''', i)
        if c:
            ua_date.append(c.group())
        if not c:
            ua_date.append('None')
            
        d = re.search(r'''(\d+\s\w+)(?=\s*Female|\s*Male)''', i)
        if d:
            age.append(d.group())
        if not d:
            age.append('None')
        
        e = re.search(r'''(\d+)(?=\s*\d{4}[A-Z])''', i)
        if e:
            mrn.append(e.group())
        if not e:
            mrn.append('None')
            
        f = re.search(r'''(Female|Male)''', i)
        if f:
            sex.append(f.group())
        if not f:
            sex.append('Not found')
    
    #Create list of lists
    name_acc_list = [list(x) for x in zip(name, accession, ua_date, age, sex, mrn)]
    
    na_df  = pd.DataFrame(name_acc_list, columns = ['name', 'ua_accession_no', 'ua_date', 'age', 'sex', 'ua_mrn'])
    
    na_df['ua_date'] = na_df['ua_date'].str.split(expand=True)
    na_df['ua_date'] = na_df['ua_date'].astype('str')
    na_df = na_df[na_df.ua_date != 'None']
    na_df['ua_date'] = pd.to_datetime(na_df['ua_date'], format='%d%b%y')
    na_df = na_df.drop_duplicates(subset=['ua_accession_no'])
    na_df[['age_yr', 'age_t']] = na_df['age'].str.split(expand=True)
    na_df['age_yr'] = na_df.apply(
                                lambda x: round(float(x['age_yr']) / float(12),2) 
                                if x['age_t'] == 'months' 
                                else x['age_yr'], axis=1
                                )
    na_df['age_yr'] = na_df.apply(
                                lambda x: round(float(x['age_yr']) / float(52),2) 
                                if x['age_t'] == 'weeks' 
                                else x['age_yr'], axis=1
                                )

    na_df = na_df.drop(columns=['age_t'])
    
    na_df.age_yr = na_df.age_yr.astype(float)
    return na_df
