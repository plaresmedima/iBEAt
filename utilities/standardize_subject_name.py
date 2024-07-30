import re

def subject_hifen (subject_name):

    pattern = r'\d{4}-\d{3}'
    matches = re.findall(pattern, subject_name)

    if matches:
        matches = matches[0]
        matches = matches.replace('-', '_')
        #print(matches)
        return matches
    else:
        return 0
    
def subject_underscore (subject_name):

    pattern = r'\d{4}_\d{3}'
    matches = re.findall(pattern, subject_name)

    if matches:
        return matches[0]
    else:
        return 0
    

def subject_seven_digitd (subject_name):

    pattern = r'\d{7}'
    matches = re.findall(pattern, subject_name)

    if matches:
        matches = matches[0]
        matches = matches[0:4]+"_"+matches[4:]
        return matches
    else:
        return 0