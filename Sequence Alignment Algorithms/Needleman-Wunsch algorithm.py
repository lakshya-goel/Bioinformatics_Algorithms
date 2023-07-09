#example sequences: (1) G A T G C G C A G, (2) G G C A G T A

a = input("First sequence: ").split()
b = input("Second sequence: ").split()

#example values: match = 2, mismatch = -1, gap = -3

match_score = int(input("Match score: "))
mismatch_score = int(input("Mismatch score: "))
gap_penalty = int(input("Gap score: "))


adj_matrix = []
matching = []
sequences = []

m = len(a)
n = len(b)


#used to calculate the values for filling up the matrix
def optimum_score(i,j):
    matching_score = adj_matrix[i-1][j-1] + mismatch_score
    if(a[i-1] == b[j-1]):
        matching_score = adj_matrix[i-1][j-1] + match_score

    gap_score = adj_matrix[i-1][j] + gap_penalty
    if(adj_matrix[i][j-1]>adj_matrix[i-1][j]):
        gap_score = adj_matrix[i][j-1] + gap_penalty

    return max(matching_score,gap_score)

#used to generate all possible sequence pairings storing one possible sequence in 'matching' list at a time
def get_sequences(i,j):

    if(i<1 or j<1):                  #base-case
        list1 =[]
        list2 =[]
        for k in matching:
            list1.append(k[0])
            list2.append(k[1])

        list1.reverse()
        list2.reverse()
      
        if(i!=0):
            for k in range(i):
                list1.insert(k,a[k])
                list2.insert(k,"-")
        elif(j!=0):
            for k in range(j):
                list2.insert(k,b[k])
                list1.insert(k,"-")

        sequences.append(list1)
        sequences.append(list2)
        return

    current = adj_matrix[i][j]
    diagonal = adj_matrix[i-1][j-1]
    left = adj_matrix[i][j-1]
    up = adj_matrix[i-1][j]

    if( (a[i-1]==b[j-1]) and (current - diagonal == match_score) ) or ( (a[i-1]!=b[j-1]) and (current - diagonal == mismatch_score) ):
        matching.append([a[i-1],b[j-1]])

        get_sequences(i-1,j-1)

        matching.pop()

    if( up - current == -gap_penalty ):
        matching.append([a[i-1],'-'])

        get_sequences(i-1,j)

        matching.pop()

    if ( left - current == -gap_penalty):
        matching.append(['-',b[j-1]])

        get_sequences(i,j-1)

        matching.pop()

    return


#initializing adj_matrix with 0s
for i in range(m+1):
    row = []
    for j in range(n+1):
        row.append(0)
    adj_matrix.append(row)

#initializing adj_matrix's first row and first column
for i in range(m+1):
    adj_matrix[i][0] = gap_penalty*i
for i in range(n+1):
    adj_matrix[0][i] = gap_penalty*i

#filling up rest of the matrix
for i in range(1,m+1):
    for j in range(1,n+1):
        adj_matrix[i][j] = optimum_score(i,j)

print()

#printing the adjacency matrix
for i in adj_matrix:
    print(i)

print()

#generating sequences
get_sequences(m,n)

count = 0
p=0

#printing sequences
while(p<len(sequences)):
    str1=""
    str2=""
    for i in sequences[p]:
        str1+=i + " "
    for i in sequences[p+1]:
        str2+=i + " "
    print(str1)
    for j in range(len(sequences[p])):
        print("|",end=" ")

    print("")
    print(str2)
    # print(sequences[p+1])
    p+=2
    count+=1
    print()

print("Total sequences possible with score: ",adj_matrix[m][n]," are : ",count)
print()