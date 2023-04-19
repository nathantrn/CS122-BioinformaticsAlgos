#!/usr/bin/env python
# coding: utf-8

# Given a Sequence and a Pattern, count how many times the pattern occurs

# In[1]:


def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
                   if Text[i:i+len(Pattern)] == Pattern:
                       count = count + 1 
    return count


# Given a Sequence, return most frequent K-mer of length K

# In[2]:


def better_frequent_words(text, k):
    frequent_patterns = []
    freq_map = frequency_table(text, k)
    max_count = max(freq_map.values())
    for pattern, count in freq_map.items():
        if count == max_count:
            frequent_patterns.append(pattern)
    return frequent_patterns


# In[3]:


def frequency_table(text, k):
    freq_map = {}
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i+k]
        if pattern not in freq_map:
            freq_map[pattern] = 1
        else:
            freq_map[pattern] += 1
    return freq_map


# Given a sequence, return reverse complement

# In[20]:


def ReverseComplement(text):
    string = ""
    for i in range(len(text)):
        if text[i] == 'A':
            string = string + 'T'
        if text[i] == 'C':
            string = string + 'G'
        if text[i] == 'G':
             string = string + 'C'
        if text[i] == 'T':
             string = string + 'A'
    txt = string[::-1]
    return(txt)


# Given a Sequence and a Pattern, return all the locations where the Pattern Occurs

# In[6]:


def PatternLocationOccurences(Text, Pattern):
    string = ""
    for i in range(len(Text)-len(Pattern)):
                   if Text[i:i+len(Pattern)] == Pattern:
                       string = string + str(i) + " "
    return string


# Given a Sequence, return k-mers forming (L,t)-clumps

# In[7]:


def FindClumps(Text, k, L, t):
    Patterns = []
    n = len(Text)
    for i in range(n-L+1):
        Window = Text[i:i+L]
        freqMap = frequency_table(Window, k)
        for s in freqMap:
            if freqMap[s] >= t:
                Patterns.append(s)
    Patterns = list(set(Patterns))
    return Patterns


# Given a Sequence, return the Skew value (difference between C and G occurences)

# In[8]:


def Skew(Text):
    diff = 0
    diff_list = []
    diff_list.append(0)
    for i in range(len(Text)):
        if Text[i] == 'C':
            diff = diff - 1
        if Text[i] == 'G':
            diff = diff + 1
        diff_list.append(diff)
    return diff_list


# Given a Sequence, determine the points of inflection in the Skew (helps determine position of ORI)

# In[10]:


def PosMinSkew(Text):
    diff = 0
    diff_list = []
    diff_list.append(0)
    for i in range(len(Text)):
        if Text[i] == 'C':
            diff = diff - 1
        if Text[i] == 'G':
            diff = diff + 1
        diff_list.append(diff)
    x = min(diff_list)
    min_pos = []
    for i in range(len(diff_list)):
        if x == diff_list[i]:
            min_pos.append(i)
    return min_pos


# Given two sequences, determine Hamming Distance

# In[12]:


def HammingDist(Text1, Text2):
    dist = 0
    for i in range(len(Text1)):
        if Text1[i] != Text2[i]:
            dist += 1
    return dist
        


# Find the most frequent k-mers (with mismatches and reverse complements) in a string.

# In[13]:


def PatternHamming(Text, Pattern, d):
    count = 0
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count + 1 
            positions.append(i)
        elif HammingDist(Text[i:i+len(Pattern)], Pattern) <= d:
            count = count + 1
            positions.append(i)
    string = ""
    for i in range(len(positions)):
        string = string + str(positions[i]) + " "
    return string


# In[14]:


def PatternHammingCount(Text, Pattern, d):
    count = 0
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count + 1 
            positions.append(i)
        elif HammingDist(Text[i:i+len(Pattern)], Pattern) <= d:
            count = count + 1
            positions.append(i)
    string = ""
    for i in range(len(positions)):
        string = string + str(positions[i]) + " "
    return count


# In[16]:


def Suffix(Text, i):
    return Text[i:]


# In[17]:


def Neighbors(Pattern, d):
    if d == 0:
        return [Pattern]
    if len(Pattern) == 1:
        return ['A', 'C', 'G', 'T']
    Neighborhood = []
    nucleotides = ['A', 'T', 'G', 'C']
    SuffixNeighbors = Neighbors(Suffix(Pattern,1), d)
    for Text in SuffixNeighbors:
        if HammingDist(Suffix(Pattern,1), Text) < d:
            for nucleotide in ['A', 'C', 'G', 'T']:
                Neighborhood.append(nucleotide + Text)
        else:
            Neighborhood.append(Pattern[0] + Text)
    return Neighborhood


# In[21]:


def FrequentWordsWithMismatches(Text, k, d):
    Patterns = set()
    freqMap = {}
    n = len(Text)
    for i in range(n-k):
        Pattern = Text[i:i+k]
        ReverseComp = ReverseComplement(Pattern)
        neighborhood = Neighbors(Pattern, d)
        for neighbor in neighborhood:
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1
        neighborhood = Neighbors(ReverseComp,d)
        for neighbor in neighborhood:
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1
    m = max(freqMap.values())
    for Pattern, count in freqMap.items():
        if count == m:
            Patterns.add(Pattern)
    return Patterns


# Given a sequence of reads, return all (k,d) motifs in the DNA

# In[24]:


def MotifEnumeration(Dna, k, d):
    Patterns = set()
    reads = []
    s = ""
    for i in range(len(Dna)+1):
        if i == len(Dna):
            reads.append(s)
        elif Dna[i] == " ":
            reads.append(s)
            s = ""
        else:
            s = s + Dna[i]
    neighborhood = []
    for j in range(len(reads[0])-k+1):
        Pattern = reads[0][j:j+k]
        neighborhood.append(Pattern)
        for m in range(len(Neighbors(Pattern, d))): 
            neighborhood.append(Neighbors(Pattern, d)[m])
    for neighbor in neighborhood:
        count = 1
        for read in reads:
            for r in range(len(read)-k+1):
                if HammingDist(neighbor,read[r:r+k]) <= d:
                    count = count + 1
                    break
        if count == len(reads)+1:
            Patterns.add(neighbor)
    string = ""
    for items in Patterns:
        string = string + str(items) + " "
    return string


# Construct a Trie

# In[25]:


class TrieNode:
    def __init__(self,num):
        self.num = num
        self.children = {}
        self.is_end_of_word = False

class Trie:
    def __init__(self):
        self.root = TrieNode(0)

    def insert(self, word):
        current_node = self.root
        for char in word:
            if char not in current_node.children:
                current_node.children[char] = TrieNode()
            current_node = current_node.children[char]
        current_node.is_end_of_word = True

def TrieConstruction(patterns):
    trie = Trie()
    counter = 0
    for pattern in patterns:
        current_node = trie.root
        for char in pattern:
            if char in current_node.children:
                current_node = current_node.children[char]
            else:
                counter += 1
                new_node = TrieNode(counter)
                current_node.children[char] = new_node
                current_node = new_node
        current_node.is_end_of_word = True
    stack = []
    stack.append(trie.root)
    output = ""
    while stack: 
        current_node = stack.pop()
        for child in current_node.children.keys():
            stack.append(current_node.children[child])
            output += f"{current_node.num} {current_node.children[child].num} {child}\n"
    pyperclip.copy(output)
    return trie


# Construct a Burrow Wheeler Tranform and its Inverse

# In[26]:


def BWT(Text):
    rotations = []
    for i in range(len(Text)):
        string = Text[len(Text)-i:] + Text[:len(Text)-i]
        rotations.append(string)
    rotations = sorted(rotations)
    BWT = ""
    for items in rotations:
        BWT = BWT + items[len(items)-1]
    return BWT
    


# In[27]:


def BWT_inverse(s):
    string = []
    for c in s:
        string.append(c)
    tmp = string.copy()
    tmp.sort()
    for i in range(len(s)-1):
        for i in range(len(string)):
            string[i] += tmp[i][-1]
        tmp = string.copy()
        tmp.sort()
    for item in string:
        if item[len(item)-1] == "$":
            return item


# In[ ]:




