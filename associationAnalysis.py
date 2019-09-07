import re
import itertools
import pandas as pd

def read(file_path):
    file_object = open(file_path, "r")
    database = list()
    for line in file_object:
        line_list = re.split(r'\t+', line.rstrip('\n'))
        for col in range(len(line_list)-1):
            line_list[col] = "gene"+str(col+1)+"_"+line_list[col]
        line_list[len(line_list)-1] = "outcome_"+line_list[len(line_list)-1]
        database.append(set(line_list))
    return database


def frequent(database, min_support, candidate_Ck, k):
    frequent_Lk = list()
    cnt = 0
    for candidate in candidate_Ck:
        for trans in database:
            if trans.issuperset(candidate):
                cnt = cnt + 1
        if cnt/100 >= min_support:
            frequent_Lk.append(candidate)
        cnt = 0

    print("number of length-" + str(k) + " frequent itemsets:"+str(len(frequent_Lk)))
    return frequent_Lk



def selfjoinCandidate(frequent_Lk, k):
    m = k +1
    candidate_Cm = list()
    # first pick two sets in frequent_Lk nonrepeatedly to compare
    for i in range(len(frequent_Lk)-1):
        frequent_i = frequent_Lk[i]
        for j in range(i+1, len(frequent_Lk)):
            frequent_j = frequent_Lk[j]
            # compare k-1 elements to self join length k+1 candidates
            if sorted(frequent_i)[:k-1] == sorted(frequent_j)[:k-1]:
                candidate = frequent_i.union(frequent_j)
                candidate_Cm.append(candidate)

    print("*found " + str(len(candidate_Cm)) + " candidate itemsets with length m = " + str(m))
    return candidate_Cm

def prune(frequent_Lk, candidate_Cm, m):

    candidate_CmPrune = list()
    indexSet = set()

    for l in range(len(candidate_Cm)):
        for candidate_subk in set(itertools.combinations(candidate_Cm[l], m-1)):

            # check if this subset present as an element in list frequent_Lk
            i = 0
            while i < len(frequent_Lk):

                if frequent_Lk[i] == set(candidate_subk):
                    # subset being frequent, check next subset
                    break
                i = i + 1
            if i == len(frequent_Lk) and frequent_Lk[i-1]!= set(candidate_subk):
                # find infrequent subset
                indexSet.add(l)
                break
    for j in range(len(candidate_Cm)):
        if indexSet.issuperset({j}) == False:
            candidate_CmPrune.append(candidate_Cm[j])

    print("*after pruning found " + str(len(candidate_CmPrune)) + " candidates")
    return candidate_CmPrune


def apriori(database, min_support):
    print("support is set to be "+str(min_support))
    # generate all 204 item candidate set
    candidate_C1 = list()
    for gene in range(100):
        candidate_C1.append({"gene" + str(gene + 1) + "_Up"})
        candidate_C1.append({"gene" + str(gene + 1) + "_Down"})
    candidate_C1.append({"outcome_ALL"})
    candidate_C1.append({"outcome_Breast Cancer"})
    candidate_C1.append({"outcome_AML"})
    candidate_C1.append({"outcome_Colon Cancer"})

    frequent_L1 = frequent(database, min_support, candidate_C1, 1)

    frequent_Lk = frequent_L1
    frequent_L = list()

    k = 1
    # repeat until no new frequent itemsets at length k generated or k reaches the max
    while len(frequent_Lk) >= 1 and k <= 204:
        frequent_L = frequent_L + frequent_Lk
        # m = k + 1
        candidate_Cm = selfjoinCandidate(frequent_Lk, k)
        # add prune function here
        candidate_CmPrune = prune(frequent_Lk, candidate_Cm, k+1)

        k = k + 1
        frequent_Lk = frequent(database, min_support, candidate_CmPrune, k)

    print("number of all lengths frequent itemsets: " + str(len(frequent_L)))
    return frequent_L



def gen_first_rules(itemset):
    # do the first generation of a frequent itemset
    # ex: for freq_set{a,b,c,d}
    #     generate{{a,b,c}->{d},{a,b,d}->{c},{a,c,d}->{b},{b,c,d}->{a}}
    rules = list()
    for item in itemset:
        body = {item}
        head = itemset-body
        rule = (head,body)
        rules.append(rule)
    return rules

def selfjoin_rules(rules):
    new_rules = list()
    for i in range(len(rules)):
        rule_i = rules[i]
        if (len(rule_i[0]) > 1):
            for j in range(i+1, len(rules)):
                rule_j = rules[j]
                # join({head1},{body1}),({head2},{body2}) -> ({head1&head2},{body1|body2})
                if(len(rule_j[0]) > 1 and len(rule_i[0] & rule_j[0]) > 0):
                    new_rule = (rule_i[0] & rule_j[0], rule_i[1] | (rule_j[1]))
                    if (new_rule not in new_rules):
                        new_rules.append(new_rule)
    print("*generated " + str(len(new_rules)) + " candidate rules with HEAD length  = " + str(len(new_rules[0][0])))
    return new_rules

def prune_rules(rules, prev_rules):
    for rule in rules:
        (head, body) = rule
        for item in body:
            parent_rule = (head | {item}, body - {item})
            if (parent_rule not in prev_rules):
                rules.remove(rule)
                break
    print("*after pruning found " + str(len(rules)) + " candidate rules" )
    return rules

def select_rules(database, candidate_rules, min_confidence):
    high_conf_rules = list()
    # calculate the confidence of each rule
    # if rules have high confidence, put them in the list
    for cand_rule in candidate_rules:
        head = cand_rule[0]
        body = cand_rule[1]
        unionset = union(head,body)
        numerator = 0
        denominator = 0
        for trans in database:
            if trans.issuperset(unionset):
                numerator += 1
            if trans.issuperset(head):
                denominator += 1
        if denominator>0:
            if numerator/denominator >= min_confidence:
                high_conf_rules.append(cand_rule)
    return high_conf_rules
    

def rule_generation(database, freq_sets, min_confidence):
    asso_rules = list()
    for itemset in freq_sets:
        n = len(itemset)
        k = n
        #generate new rules until last selected rules empty or head size = 1
        if k>1:
            cand_rules = gen_first_rules(itemset)
            new_rules = select_rules(database, cand_rules, min_confidence)
            asso_rules+=new_rules
            k=k-1
            while len(new_rules)>=1 and k>1:
                cand_rules = selfjoin_rules(new_rules)
                cand_rules = prune_rules(cand_rules, new_rules)
                new_rules = select_rules(database, cand_rules, min_confidence)
                asso_rules+=new_rules
                k=k-1
            
        for rule in asso_rules:
            print(rule)
            
    return asso_rules
                

def print_rules_result(asso_rules):
    # to do
    for rule in asso_rules:
        print(rule)


def intersection(lst1, lst2): 
    return set(lst1) & set(lst2)

def union(lst1, lst2): 
    return set(lst1) | set(lst2)

def template1(df, rule, num, item_list):   
    temp_df = pd.DataFrame(columns = item_list)
    temp_df['result'] = pd.Series(True, index=df.index)

    for item in item_list:
        if rule == 'RULE':
            temp_df[item]= df.HEAD.str.contains(item) | df.BODY.str.contains(item)
        else:
            temp_df[item]= df[rule].str.contains(item)
        temp_df['result'] = temp_df[item] & temp_df.result
    
    if num == 'NONE':
        cnt = temp_df.result.value_counts()[0]  
        result_list = temp_df.index[temp_df['result'] == False].tolist()
    else:
        cnt = temp_df.result.value_counts()[1]
        result_list = temp_df.index[temp_df['result'] == True].tolist()
        
    return result_list, cnt


def template2(df, rule, num):
    temp_df = pd.Series(0, index=df.index)
    
    for i in range(len(temp_df.index)):
        if rule == 'RULE':
            temp_df[i] = df.iloc[i].HEAD.count(', ') + df.iloc[i].BODY.count(', ') + 2
        else:
            temp_df[i] = df.iloc[i][rule].count(', ') + 1
  
    result = temp_df.loc[lambda x : x >= num ]
    return result.index.tolist(), result.size


def template3(df, query_type, rule1, num1, rule2, num2, item_list1=[], item_list2=[]):
    combined_result = []
    
    if query_type.count('1') == 2:
        result1, cnt1 = template1(df, rule1, num1, item_list1)
        result2, cnt2 = template1(df, rule2, num2, item_list2)
    elif query_type.count('2') == 2:
        result1, cnt1 = template2(df, rule1, num1)
        result2, cnt2 = template2(df, rule2, num2)
    else:
        result1, cnt1 = template1(df, rule1, num1, item_list1)
        result2, cnt2 = template2(df, rule2, num2)

    combined_result = union(result1, result2) if query_type.count('or') == 1 else intersection(result1, result2)
    return len(combined_result)


def main():

    database = read("./associationruletestdata.txt")
    # for element in database:
    #     print(element)

    frequent_L = apriori(database, 0.5)
    conf_rules = rule_generation(database, frequent_L, 0.7)
    #print_rules_result(conf_rules)
    
    df = pd.DataFrame(conf_rules, columns = ['HEAD','BODY'])
    df[['HEAD','BODY']] = df[['HEAD','BODY']].astype(str)
    #print(template1(df,'RULE',1,['gene82_Down']))
    #print(template2(df,'RULE',3))
    #print(template3(df,'1and2','RULE',1,'RULE',3,['gene82_Down']))

    

if __name__ == '__main__':
    main()
