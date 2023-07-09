alpha_p = {"E":1.53, "A":1.45, "L":1.34, "H":1.24, "M":1.20, "Q":1.17, "W":1.14, "V":1.14, "F":1.12, "K":1.07, "I":1.00, "D":0.98, "T":0.82, "S":0.79, "R":0.79, "C":0.77, "N":0.73, "Y":0.61, "P":0.59, "G":0.53}   
beta_p = {"M":1.67, "V":1.65, "I":1.60, "C":1.30, "Y":1.29, "F":1.28, "Q":1.23, "L":1.22, "T":1.20, "W":1.19, "A":0.97, "R":0.90, "G":0.81, "D":0.80, "K":0.74, "S":0.72, "H":0.71, "N":0.65, "P":0.62, "E":0.26}    

seq = input("Enter sequence: ")

#example sequence: SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF

l = len(seq)

alpha_seq_mapping = []
beta_seq_mapping = []
helix_mapping = []
strand_mapping = []

#Initialising these lists
for i in range(l):
    alpha_seq_mapping.append(alpha_p[seq[i]])
    beta_seq_mapping.append(beta_p[seq[i]])
    helix_mapping.append('-')
    strand_mapping.append('-')
    
def extensionR (init_score, start, end, arr):
    score = init_score
    while(end<l):
        score -= arr[start]
        score += arr[end]
        if(score<4):
            return end
        start+=1
        end+=1
    return end

def extensionL (init_score, start, end, arr):
    score = init_score
    while(start>=0):
        score -= arr[end-1]
        score += arr[start-1]
        if(score<4):
            return start
        start-=1
        end-=1
    return start

#Since strings are immutable
def list_to_string(l):
    s=''
    for i in l:
        s+=i
    return s


#Identifying helices
for i in range(l-5):
    start = i
    end = i+6
    window = seq[start:end]

    #is valid? check
    c=0
    for j in range(start,end):
        if(alpha_seq_mapping[j]>=1):
            c+=1
    if(c<4):
        continue

    #found nucleation site
    Rscore = sum(alpha_seq_mapping[end-4:end])
    Lscore = sum(alpha_seq_mapping[start:start+4])

    new_start = extensionL(Lscore,start,start+4,alpha_seq_mapping)
    new_end = extensionR(Rscore,end-4,end,alpha_seq_mapping)

    for j in range(new_start,new_end):
        helix_mapping[j] = "H"

#Identifying strands
for i in range(l-4):
    start = i
    end = i+5
    window = seq[start:end]

    #is valid? check
    c=0
    for j in range(start,end):
        if(beta_seq_mapping[j]>=1):
            c+=1
    if(c<3):
        continue

    #found nucleation site
    Rscore = sum(beta_seq_mapping[end-4:end])
    Lscore = sum(beta_seq_mapping[start:start+4])

    new_start = extensionL(Lscore,start,start+4,beta_seq_mapping)
    new_end = extensionR(Rscore,end-4,end,beta_seq_mapping)

    for j in range(new_start,new_end):
        strand_mapping[j] = "S"

#secondary structure will be stored here
ans = ['-' for i in range(l)]

#helper method for conflict resolution
def count_sequence_length(i,arr,s):
    for j in range(i,l):
        if arr[j]!=s:
            return j
    return l

def count_dig(n):
    return len(str(n))

#Filling ans list
i=0
while i<l:
    if helix_mapping[i]=='-' and strand_mapping[i]=='S':
        ans[i] = 'S';
    elif helix_mapping[i]=='H' and strand_mapping[i]=='-':
        ans[i] = 'H';
    elif helix_mapping[i]=='-' and strand_mapping[i]=='-':
        ans[i] = '-';
    else:
        #overlap condition. We count the length of Helix and Strand from this point. Minimum
        #of the two is the overlaping part. We calculate the propensity scores for this length
        #and select the one with higher score to be displayed in ans. i moves to the end of overlap
        stop_H=count_sequence_length(i,helix_mapping,'H')
        stop_S=count_sequence_length(i,strand_mapping,'S')

        stop_overlap = min(stop_H,stop_S)

        #calculating respective scores from i to i+length_overlap
        s_H = sum(alpha_seq_mapping[i:stop_overlap])
        s_S = sum(beta_seq_mapping[i:stop_overlap])
        
        if(s_H>s_S):
            ans[i:stop_overlap] = 'H'*(stop_overlap-i)
        else:
            ans[i:stop_overlap] = 'S'*(stop_overlap-i)
        
        i = stop_overlap-1
    i+=1

seq = list_to_string(seq)
ans = list_to_string(ans)

#displaying output
print("\nThe predicted secondary structure for the given sequence is:\n")  

c=5                                      #number of parts to display the output in
for i in range(c):
    s = (i) * (l//c)
    g = count_dig(s)
    e = (i+1)*(l//c)
    print(s, (5-g)*" ", seq[s:e], "      ", e)
    print(6*" ", ans[s:e], "      ", "  ")
    print()
