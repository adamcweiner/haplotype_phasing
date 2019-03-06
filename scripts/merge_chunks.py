import numpy as np
import random
from time import time

def generate_random():
	overlap = 5
	start = []
	end = []
	position = 0
	test_data = []
	for k in range(5):
		start.append(position)
		end.append(position+19)
		chunk = []
		for j in range(10):
			row = []
			for i in range(10):
				row.append(random.randint(0,1))
			chunk.append(row)
			position += 1
		position -= overlap
		test_data.append(chunk)

	for i in range(len(test_data)):
		test_data[i] = np.array(test_data[i])

	for i in range(1,len(test_data)):
		test_data[i][0:overlap] = test_data[i-1][len(test_data[i-1])-overlap:len(test_data[i-1])]		
	return test_data, start, end
		

def count_match(lst1, lst2):
	count = 0
	length = len(lst1)

	for i in range(length):
		if(lst1[i]==lst2[i]):
			count+=1

	return count


#function merge_chunks(haplotype_chunks)
#input:
#	haplotype_chunks: list of the nx2m matrices, start positions, end positions
#output:
#	haplotypes: one 2mxn matrix where n is the number of total snps
def merge_chunks(data, start, end):
	people = []
	for i in range(len(data[0][0])):
		person = [item[i] for item in data[0]]
		people.append(person)

	for i in range(1,len(data)):
		overlap = end[i-1] - start[i] +1
		for j in range(0,len(data[i-1][0]),2):
			lst1 = [item[j] for item in data[i-1][len(data[i-1])-overlap:len(data[i-1])]] #ending overlapping sequence of first haplotype in a group
			lst2 = [item[j] for item in data[i][0:overlap]] #beginning overlapping sequence of first haplotype
			lst3 = [item[j+1] for item in data[i][0:overlap]] #beginning overlapping sequence of second haplotype
			lst1v2 = count_match(lst1, lst2)
			lst1v3 = count_match(lst1,lst3)
			
			if(lst1v2 > lst1v3): #puts the first two together
				people[j] += [item[j] for item in data[i][overlap:len(data[i])]]
				people[j+1] += [item[j] for item in data[i][overlap:len(data[i])]]
			else: #puts the first haplotype with second
			#note that if they are matching it will default to this but shouldn't matter either way because it would be 50/50
				people[j+1] += [item[j] for item in data[i][overlap:len(data[i])]]
				people[j] += [item[j] for item in data[i][overlap:len(data[i])]]
				
	return people


#data, start, end = generate_random()
#merge_chunks(data,start,end)

