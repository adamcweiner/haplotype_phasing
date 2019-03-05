import numpy as np
import random

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
		for j in range(20):
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
		
#function merge_chunks(haplotype_chunks)
#input:
#	haplotype_chunks: list of the nx2m matrices
#output:
#	haplotypes: one nx2m matrix where n is the number of total snps
def merge_chunks(data, start, end):
	print(data[0])
	print(data[1])
	people = []
	for i in range(len(data[0][0])):
		person = [item[i] for item in data[0]]
		people.append(person)

	for i in range(1,len(data)):
		overlap = end[i-1] - start[i] +1
		for j in range(0,len(data[i-1][0]),2):
			lst1 = [item[0] for item in data[i-1][len(data[i-1])-5:len(data[i-1])]]
			lst2 = [item[0] for item in data[i][0:5]]
			if(lst1 == lst2):
				people[j] += [item[j] for item in data[i][overlap:len(data[i])]]
				people[j+1] += [item[j] for item in data[i][overlap:len(data[i])]]
			else:
				people[j+1] += [item[j] for item in data[i][overlap:len(data[i])]]
				people[j] += [item[j] for item in data[i][overlap:len(data[i])]]
				
	return people


#data, start, end = generate_random()
#merge_chunks(data,start,end)